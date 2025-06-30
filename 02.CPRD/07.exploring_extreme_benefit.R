
# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)
library(rms)
library(tableone)
library(writexl)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}


matched_cohort_overall_benefit <- function(data, 
                                                  drug_var, 
                                                  outcome_var, 
                                                  pred_cols = NULL, 
                                                  conc_tolerance = NULL,
                                                  matching_var = NULL, 
                                                  match.exact = NULL, 
                                                  match.antiexact = NULL) {
  
  # Load required libraries ----
  require(tidyverse)  # For data manipulation and piping
  require(MatchIt)    # For matching procedures
  
  # Input validation ----
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!all(pred_cols %in% colnames(data))) stop("Some pred_cols not found in data")
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("conc_tolerance must be numeric")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  
  # Prepare data ----
  pre_data <- data %>%
    rename_with(~ "dataset_drug_var", all_of(drug_var)) %>%    # Consistent treatment variable name
    rename_with(~ "dataset_outcome_var", all_of(outcome_var))  # Consistent outcome variable name
  
  # Extract common prefix from predicted outcome columns ----
  common_prefix <- function(strings) {
    min_len <- min(nchar(strings))
    prefix <- ""
    for (i in seq_len(min_len)) {
      current_chars <- substr(strings, i, i)
      if (length(unique(current_chars)) == 1) {
        prefix <- paste0(prefix, current_chars[1])
      } else {
        break
      }
    }
    return(prefix)
  }
  
  prefix <- common_prefix(pred_cols)            # Extract common prefix
  drug_names <- gsub(prefix, "", pred_cols)     # Extract drug names by removing prefix
  
  # Assign best predicted drug per patient using helper function ----
  interim_dataset <- get_best_drugs(
    data = pre_data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    rename("function_rank1_drug_name" := paste0(prefix, "rank1_drug_name"))
  
  # Define concordance labels ----
  if (is.null(conc_tolerance)) {
    # Concordant if received exactly the top predicted drug
    interim_dataset <- interim_dataset %>%
      mutate(conc_disc_label = ifelse(dataset_drug_var == function_rank1_drug_name, 1, 0))
  } else {
    # Concordant if received any drug within tolerance of best predicted outcome
    n_drugs_required <- length(pred_cols)
    tolerance_vars <- paste0("tolerance_drug_", 1:n_drugs_required)
    
    interim_dataset <- get_best_drugs(
      data = interim_dataset,
      tolerance = conc_tolerance,
      column_names = pred_cols,
      final_var_name = prefix
    ) %>%
      rename("function_tolerance_drug_name" := paste0(prefix, "within_", conc_tolerance, "_of_best_drug_name")) %>%
      mutate(
        conc_disc_label = ifelse(str_detect(function_tolerance_drug_name, paste0("\\b", dataset_drug_var, "\\b")), 1, 0)
      ) %>%
      mutate(drug_list = str_split(function_tolerance_drug_name, "\\s*[;,\\s]\\s*")) %>%
      mutate(drug_list = map(drug_list, ~ rep(.x, length.out = n_drugs_required))) %>%
      mutate(drug_list = map(drug_list, ~ set_names(.x, paste0("tolerance_drug_", seq_along(.x))))) %>%
      unnest_wider(drug_list) %>%
      mutate(across(
        all_of(tolerance_vars),
        ~ if_else(conc_disc_label == 0, dataset_drug_var, .x)
      ))
  }
  
  # Perform matching to control confounding ----
  interim_dataset <- interim_dataset %>% drop_na(all_of(matching_var))
  
  categorical_vars <- matching_var[sapply(interim_dataset[matching_var], \(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
  
  # Dynamically build matching formula ----
  matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))
  for (v in categorical_vars) {
    if (length(unique(interim_dataset[[v]])) > 1) {
      matching_formula <- paste(matching_formula, "+", v)
    }
  }
  
  # Antiexact matching variables ----
  match_model_antiexact_vars <- if (!is.null(conc_tolerance)) {
    unique(c(tolerance_vars, "dataset_drug_var", match.antiexact))
  } else {
    unique(c("dataset_drug_var", match.antiexact))
  }
  
  # Run MatchIt nearest neighbor matching ----
  match_model <- MatchIt::matchit(
    formula = as.formula(matching_formula),
    data = interim_dataset,
    method = "nearest",
    distance = "mahalanobis",
    replace = FALSE,
    exact = unique(c("function_rank1_drug_name", match.exact)),
    antiexact = match_model_antiexact_vars
  )
  
  # Extract matched pairs ----
  matched_data <- MatchIt::get_matches(match_model, data = interim_dataset)
  
  # Calculate observed and predicted differences within matched pairs ----
  processed_data <- matched_data %>%
    group_by(subclass) %>%
    mutate(
      concordant_drugclass  = dataset_drug_var[1],
      discordant_drugclass  = dataset_drug_var[2],
      calibration_obs       = diff(dataset_outcome_var)  # discordant - concordant outcome difference
    ) %>%
    ungroup() %>%
    distinct(subclass, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(
      calibration_pred = get(paste0(prefix, discordant_drugclass)) - get(paste0(prefix, concordant_drugclass))  # predicted difference
    ) %>%
    ungroup()
  
  # Return final results ----
  return(processed_data)
  
}




sapply(
  paste0("01.Functions/", list.files("01.Functions", pattern = "\\.R$")),
  source
)
source("02.CPRD/03.impute_missingness.R")
source("02.CPRD/04.model_predictions.R")

# load model
load("fivedrugmodel_5knot_share_20230823.Rdata")

# set up dataset
analysis_post_2020_raw <- analysis_post_2020_raw %>%
  analysis$cached("analysis_post_2020") %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

analysis_pre_2020_raw <- analysis_pre_2020_raw %>%
  analysis$cached("analysis_pre_2020") %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)


# Pre-processing datasets ########################################

## Post 2020-10-14 ----
analysis_post_2020 <- analysis_post_2020_raw %>%
  mutate(
    pated = paste(patid, dstartdate, drug_class, sep = "."),
    sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    agetx = dstartdate_age,
    ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
    ethnicity = factor(ethnicity, levels = c(0:5), labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),
    
    smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),
    smoke = factor(smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),
    imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
    imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5), labels = c("1 (least)", "2", "3", "4", "5 (most)")),
    
    ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
    ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
    ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),
    drugline = ifelse(drugline_all > 5, 5, drugline_all),
    drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),
    drugclass = drug_class
  ) %>%
  group_by(pated) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  select(-row) %>%
  select(all_of(c(
    "pated", "agetx", "sex", "t2dmduration", "ethnicity",
    "drug_substance", "drug_class",
    "imd5", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  )
  )) %>%
  rename("drugclass" = "drug_class") %>%
  as.data.frame()


## Pre 2020-10-14 ----
analysis_pre_2020 <- analysis_pre_2020_raw %>%
  mutate(
    pated = paste(patid, dstartdate, drug_class, sep = "."),
    sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    agetx = dstartdate_age,
    ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
    ethnicity = factor(ethnicity, levels = c(0:5), labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),
    
    smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),
    smoke = factor(smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),
    imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
    imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5), labels = c("1 (least)", "2", "3", "4", "5 (most)")),
    
    ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
    ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
    ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),
    drugline = ifelse(drugline_all > 5, 5, drugline_all),
    drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),
    drugclass = drug_class
  ) %>%
  group_by(pated) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  select(-row) %>%
  select(all_of(c(
    "pated", "agetx", "sex", "t2dmduration", "ethnicity",
    "drug_substance", "drug_class",
    "imd5", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  )
  )) %>%
  rename("drugclass" = "drug_class") %>%
  as.data.frame()



# Missing data imputation ######################################

## Post 2020-10-14 ----
### By subgroups
analysis_post_2020_group_imputation <- imputation_methods(data = analysis_post_2020,
                                                          method = "group")

### By mice
analysis_post_2020_mice_imputation <- imputation_methods(data = analysis_post_2020,
                                                         method = "mice",
                                                         mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))


# Predictions from datasets #########################################

## Post 2020-10-14 ----
### Original variables
analysis_post_2020_prediction_orig <- predict_5drugmodel(analysis_post_2020,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"),
                                                         pred_col = "pred.orig.")

### Imputation by groups
analysis_post_2020_prediction_group <- predict_5drugmodel(analysis_post_2020_group_imputation,
                                                          model = m1.5.final,
                                                          drug_var = "drugclass",
                                                          drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))

### Imputation by mice
analysis_post_2020_prediction_mice <- predict_5drugmodel(analysis_post_2020_mice_imputation,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))


### merge impute columns into main dataset
analysis_post_2020 <- analysis_post_2020 %>%
  cbind(
    analysis_post_2020_prediction_orig %>%
      select(contains("pred.orig")) %>%
      rename_with(~ str_replace(., "pred\\.orig\\.", "pred.orig.preclosed.")),
    analysis_post_2020_prediction_group %>%
      select(contains("group_impute")),
    analysis_post_2020_prediction_group %>%
      select(contains("pred.group")) %>%
      rename_with(~ str_replace(., "pred\\.group\\.", "pred.group.preclosed.")),
    analysis_post_2020_prediction_mice %>%
      select(contains("mice_impute")),
    analysis_post_2020_prediction_mice %>%
      select(contains("pred.mice")) %>%
      rename_with(~ str_replace(., "pred\\.mice\\.", "pred.mice.preclosed."))
  )



## Pre 2020-10-14 ----
### Original variables
analysis_pre_2020_prediction_orig <- predict_5drugmodel(analysis_pre_2020,
                                                        model = m1.5.final,
                                                        drug_var = "drugclass",
                                                        drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"),
                                                        pred_col = "pred.orig.")



### merge impute columns into main dataset
analysis_pre_2020 <- analysis_pre_2020 %>%
  cbind(
    analysis_pre_2020_prediction_orig %>%
      select(contains("pred.orig")) %>%
      rename_with(~ str_replace(., "pred\\.orig\\.", "pred.orig.preclosed."))
  )


# Closed loop test #################################################

## Post 2020-10-14 ----

### Original variables
#### SGLT2
closed_loop_test_results_SGLT2_post_2020_orig <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "SGLT2") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_post_2020_orig <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "GLP1") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_post_2020_orig <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "DPP4") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_post_2020_orig <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "TZD") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_post_2020_orig <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "SU") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)


### group variables
#### SGLT2
closed_loop_test_results_SGLT2_post_2020_group <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "SGLT2") %>%
    mutate(
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_post_2020_group <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "GLP1") %>%
    mutate(
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_post_2020_group <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "DPP4") %>%
    mutate(
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_post_2020_group <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "TZD") %>%
    mutate(
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_post_2020_group <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "SU") %>%
    mutate(
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)



### mice variables
#### SGLT2
closed_loop_test_results_SGLT2_post_2020_mice <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "SGLT2") %>%
    mutate(
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_post_2020_mice <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "GLP1") %>%
    mutate(
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_post_2020_mice <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "DPP4") %>%
    mutate(
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_post_2020_mice <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "TZD") %>%
    mutate(
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_post_2020_mice <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_post_2020 %>% filter(drugclass == "SU") %>%
    mutate(
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)



#### Make predictions ----
analysis_post_2020 <- analysis_post_2020 %>%
  mutate(
    # original variables
    pred.orig.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "SGLT2")),
    pred.orig.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "GLP1")),
    pred.orig.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "DPP4")),
    pred.orig.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "TZD")),
    pred.orig.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_orig, analysis_post_2020 %>% mutate(drugclass == "SU")),
    # group variables
    pred.group.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_group, analysis_post_2020 %>% mutate(
      drugclass = "SGLT2",
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    )),
    pred.group.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_group, analysis_post_2020 %>% mutate(
      drugclass = "GLP1",
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    )),
    pred.group.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_group, analysis_post_2020 %>% mutate(
      drugclass = "DPP4",
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    )),
    pred.group.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_group, analysis_post_2020 %>% mutate(
      drugclass = "TZD",
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    )),
    pred.group.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_group, analysis_post_2020 %>% mutate(
      drugclass = "SU",
      prebmi = prebmi_group_impute,
      preegfr = preegfr_group_impute,
      pretotalcholesterol = pretotalcholesterol_group_impute,
      prehdl = prehdl_group_impute,
      prealt = prealt_group_impute
    )),
    # mice variables
    pred.mice.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_mice, analysis_post_2020 %>% mutate(
      drugclass = "SGLT2",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.mice.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_mice, analysis_post_2020 %>% mutate(
      drugclass = "GLP1",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.mice.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_mice, analysis_post_2020 %>% mutate(
      drugclass = "DPP4",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.mice.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_mice, analysis_post_2020 %>% mutate(
      drugclass = "TZD",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.mice.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_mice, analysis_post_2020 %>% mutate(
      drugclass = "SU",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ))
  )



## Pre 2020-10-14 ----

### Original variables
#### SGLT2
closed_loop_test_results_SGLT2_pre_2020_orig <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_pre_2020 %>% filter(drugclass == "SGLT2") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_pre_2020_orig <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_pre_2020 %>% filter(drugclass == "GLP1") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_pre_2020_orig <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_pre_2020 %>% filter(drugclass == "DPP4") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_pre_2020_orig <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_pre_2020 %>% filter(drugclass == "TZD") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_pre_2020_orig <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_pre_2020 %>% filter(drugclass == "SU") %>% drop_na(),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)


#### Make predictions ----
analysis_pre_2020 <- analysis_pre_2020 %>%
  mutate(
    # original variables
    pred.orig.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_pre_2020_orig, analysis_pre_2020 %>% mutate(drugclass = "SGLT2")),
    pred.orig.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_pre_2020_orig, analysis_pre_2020 %>% mutate(drugclass = "GLP1")),
    pred.orig.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_pre_2020_orig, analysis_pre_2020 %>% mutate(drugclass = "DPP4")),
    pred.orig.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_pre_2020_orig, analysis_pre_2020 %>% mutate(drugclass = "TZD")),
    pred.orig.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_pre_2020_orig, analysis_pre_2020 %>% mutate(drugclass = "SU"))
  )




# Overall calibration (summary) #########################################

## Tolerance 3 mmol/mol ----

### Post 2020-10-14 ----

### Exact matching only on best drug / sex / hba1c_10
matched_cohort_tolerance3_post_2020 <- matched_cohort_overall_benefit(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
) %>%
  mutate(grouping = ntile(calibration_pred, 10))


### Pre 2020-10-14 ----

### Exact matching only on best drug / sex / hba1c_10
matched_cohort_tolerance3_pre_2020 <- matched_cohort_overall_benefit(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
) %>%
  mutate(grouping = ntile(calibration_pred, 10))



## Rank 1 ----

### Post 2020-10-14 ----

### Exact matching only on best drug / sex / hba1c_10
matched_cohort_rank1_post_2020 <- matched_cohort_overall_benefit(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
) %>%
  mutate(grouping = ntile(calibration_pred, 10))


### Pre 2020-10-14 ----

### Exact matching only on best drug / sex / hba1c_10
matched_cohort_rank1_pre_2020 <- matched_cohort_overall_benefit(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
) %>%
  mutate(grouping = ntile(calibration_pred, 10))



# Grouping: Tables of characteristics ----
vars <- c(
  "agetx", "sex", "t2dmduration", "ethnicity", 
  "concordant_drugclass", "discordant_drugclass",
  "imd5", "smoke",
  "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
  "drugline", "ncurrtx"
)

cat_vars <- c(
  "sex", "ethnicity", "drug_substance", "imd5", "smoke",
  "drugline", "ncurrtx"
)


## Tolerance 3 mmol/mol ----

### Post 2020-10-14 ----
top_group_value <- round(matched_cohort_tolerance3_post_2020 %>% filter(grouping>8) %>% select(calibration_pred) %>% min(), 2)
bottom_group_value <- round(matched_cohort_tolerance3_post_2020 %>% filter(grouping<3) %>% select(calibration_pred) %>% max(), 2)

interim_cohort <- matched_cohort_tolerance3_post_2020 %>%
  mutate(table_var = "Whole") %>%
  rbind(
    matched_cohort_tolerance3_post_2020 %>%
      filter(grouping > 8) %>%
      mutate(table_var = paste0("Benefit > ", top_group_value)),
    matched_cohort_tolerance3_post_2020 %>%
      filter(grouping < 3) %>%
      mutate(table_var = paste0("Benefit < ", bottom_group_value))
  )

table_characteristics_interim <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  strata = c("table_var"),
  data = interim_cohort,
  test = FALSE
)

table_characteristics_interim <- print(table_characteristics_interim, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_interim, "Outputs/CPRD/07.overall_benefit_groups_tolerance3_post2020.csv")


### Pre 2020-10-14 ----
top_group_value <- round(matched_cohort_tolerance3_pre_2020 %>% filter(grouping>8) %>% select(calibration_pred) %>% min(), 2)
bottom_group_value <- round(matched_cohort_tolerance3_pre_2020 %>% filter(grouping<3) %>% select(calibration_pred) %>% max(), 2)

interim_cohort <- matched_cohort_tolerance3_pre_2020 %>%
  mutate(table_var = "Whole") %>%
  rbind(
    matched_cohort_tolerance3_pre_2020 %>%
      filter(grouping > 8) %>%
      mutate(table_var = paste0("Benefit > ", top_group_value)),
    matched_cohort_tolerance3_pre_2020 %>%
      filter(grouping < 3) %>%
      mutate(table_var = paste0("Benefit < ", bottom_group_value))
  )

table_characteristics_interim <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  strata = c("table_var"),
  data = interim_cohort,
  test = FALSE
)

table_characteristics_interim <- print(table_characteristics_interim, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_interim, "Outputs/CPRD/07.overall_benefit_groups_tolerance3_pre2020.csv")



## Rank 1 ----

### Post 2020-10-14 ----
top_group_value <- round(matched_cohort_rank1_post_2020 %>% filter(grouping>8) %>% select(calibration_pred) %>% min(), 2)
bottom_group_value <- round(matched_cohort_rank1_post_2020 %>% filter(grouping<3) %>% select(calibration_pred) %>% max(), 2)

interim_cohort <- matched_cohort_rank1_post_2020 %>%
  mutate(table_var = "Whole") %>%
  rbind(
    matched_cohort_rank1_post_2020 %>%
      filter(grouping > 8) %>%
      mutate(table_var = paste0("Benefit > ", top_group_value)),
    matched_cohort_rank1_post_2020 %>%
      filter(grouping < 3) %>%
      mutate(table_var = paste0("Benefit < ", bottom_group_value))
  )

table_characteristics_interim <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  strata = c("table_var"),
  data = interim_cohort,
  test = FALSE
)

table_characteristics_interim <- print(table_characteristics_interim, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_interim, "Outputs/CPRD/07.overall_benefit_groups_rank1_post2020.csv")



### Pre 2020-10-14 ----
top_group_value <- round(matched_cohort_rank1_pre_2020 %>% filter(grouping>8) %>% select(calibration_pred) %>% min(), 2)
bottom_group_value <- round(matched_cohort_rank1_pre_2020 %>% filter(grouping<3) %>% select(calibration_pred) %>% max(), 2)

interim_cohort <- matched_cohort_rank1_pre_2020 %>%
  mutate(table_var = "Whole") %>%
  rbind(
    matched_cohort_rank1_pre_2020 %>%
      filter(grouping > 8) %>%
      mutate(table_var = paste0("Benefit > ", top_group_value)),
    matched_cohort_rank1_pre_2020 %>%
      filter(grouping < 3) %>%
      mutate(table_var = paste0("Benefit < ", bottom_group_value))
  )

table_characteristics_interim <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  strata = c("table_var"),
  data = interim_cohort,
  test = FALSE
)

table_characteristics_interim <- print(table_characteristics_interim, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_interim, "Outputs/CPRD/07.overall_benefit_groups_rank1_pre2020.csv")







# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)
library(rms)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
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




# Calibration of Predictions ##########################################

## Post 2020-10-14 ----
plot_pred_response_analysis_post_2020_preclosed <- analysis_post_2020 %>%
  pivot_longer(cols = contains("pred.orig.preclosed")) %>%
  select(name, value, drugclass, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.preclosed\\.", "", name)) %>%
  filter(drugclass == name) %>%
  ggplot(aes(x = value, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Post 2020 (Pre closed loop)")
plot_pred_response_analysis_post_2020_preclosed_glp1 <- analysis_post_2020 %>%
  pivot_longer(cols = contains("pred.orig.preclosed")) %>%
  select(name, value, drugclass, drug_substance, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.preclosed\\.", "", name)) %>%
  filter(drugclass == name) %>%
  filter(drugclass == "GLP1") %>%
  ggplot(aes(x = value, y = posthba1cfinal, colour = drug_substance)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Post 2020 (Pre closed loop)")
plot_pred_response_analysis_post_2020 <- analysis_post_2020 %>%
  pivot_longer(cols = contains("pred.orig")) %>%
  select(name, value, drugclass, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
  filter(drugclass == name) %>%
  ggplot(aes(x = value, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Post 2020 (Post closed loop)")
plot_pred_response_analysis_post_2020_glp1 <- analysis_post_2020 %>%
  pivot_longer(cols = contains("pred.orig")) %>%
  select(name, value, drugclass, drug_substance, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
  filter(drugclass == name) %>%
  filter(drugclass == "GLP1") %>%
  ggplot(aes(x = value, y = posthba1cfinal, colour = drug_substance)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Post 2020 (Post closed loop)")

plot_pred_response_analysis_post_2020_comparison <- analysis_post_2020 %>%
  pivot_longer(cols = contains("pred.orig.preclosed")) %>%
  select(name, value, drugclass, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.preclosed\\.", "", name)) %>%
  filter(drugclass == name & drugclass %in% c("GLP1", "SGLT2")) %>%
  mutate(method = "Pre recalibration") %>%
  rbind(
    analysis_post_2020 %>%
      pivot_longer(cols = contains("pred.orig")) %>%
      select(name, value, drugclass, posthba1cfinal) %>%
      mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
      filter(drugclass == name & drugclass %in% c("GLP1", "SGLT2")) %>%
      mutate(method = "Post recalibration"),
    analysis_post_2020 %>%
      pivot_longer(cols = contains("pred.orig")) %>%
      select(name, value, drugclass, posthba1cfinal) %>%
      mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
      filter(drugclass == name & !(drugclass %in% c("GLP1", "SGLT2"))) %>%
      mutate(method = "No recalibration")
  ) %>%
  ggplot(aes(x = value, y = posthba1cfinal, colour = method)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  theme_classic() +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", colour = "Recalibration method:") +
  theme(legend.position = "bottom")

## Pre 2020-10-14 ----
plot_pred_response_analysis_pre_2020 <- analysis_pre_2020 %>%
  pivot_longer(cols = contains("pred.orig")) %>%
  select(name, value, drugclass, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
  filter(drugclass == name) %>%
  ggplot(aes(x = value, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Pre 2020")
plot_pred_response_analysis_pre_2020_glp1 <- analysis_pre_2020 %>%
  pivot_longer(cols = contains("pred.orig")) %>%
  select(name, value, drugclass, drug_substance, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
  filter(drugclass == name) %>%
  filter(drugclass == "GLP1") %>%
  ggplot(aes(x = value, y = posthba1cfinal, colour = drug_substance)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Pre 2020")


# Find Optimal drug ########################################
## Rank 1 ----
### Orig
analysis_post_2020 <- get_best_drugs(
  data = analysis_post_2020,
  rank = 1,
  column_names = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.orig."
)

### Grouped imputation model predictions
analysis_post_2020 <- get_best_drugs(
  data = analysis_post_2020,
  rank = 1,
  column_names = paste0("pred.group.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.group."
)

### MICE imputation model predictions
analysis_post_2020 <- get_best_drugs(
  data = analysis_post_2020,
  rank = 1,
  column_names = paste0("pred.mice.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.mice."
)

## Tolerance 3 mmol/mol ----
### Orig
analysis_post_2020 <- get_best_drugs(
  data = analysis_post_2020,
  tolerance = 3,
  column_names = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.orig."
)

### Grouped imputation model predictions
analysis_post_2020 <- get_best_drugs(
  data = analysis_post_2020,
  tolerance = 3,
  column_names = paste0("pred.group.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.group."
)

### MICE imputation model predictions
analysis_post_2020 <- get_best_drugs(
  data = analysis_post_2020,
  tolerance = 3,
  column_names = paste0("pred.mice.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.mice."
)

### Orig pre 2020
analysis_pre_2020 <- get_best_drugs(
  data = analysis_pre_2020,
  tolerance = 3,
  column_names = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  final_var_name = "pred.orig."
)




# Overall calibration (summary) #########################################

## Tolerance 3 mmol/mol ----

### Post 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_calibration_tolerance3_post_2020_orig <- overall_predicted_benefit_performance(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_calibration_tolerance3_post_2020_orig_match_sex <- overall_predicted_benefit_performance(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_calibration_tolerance3_post_2020_orig_match_sex_hba1c <- overall_predicted_benefit_performance(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

### Pre 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_calibration_tolerance3_pre_2020_orig <- overall_predicted_benefit_performance(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_calibration_tolerance3_pre_2020_orig_match_sex <- overall_predicted_benefit_performance(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_calibration_tolerance3_pre_2020_orig_match_sex_hba1c <- overall_predicted_benefit_performance(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


## Rank 1 ----

### Post 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_calibration_rank1_post_2020_orig <- overall_predicted_benefit_performance(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_calibration_rank1_post_2020_orig_match_sex <- overall_predicted_benefit_performance(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_calibration_rank1_post_2020_orig_match_sex_hba1c <- overall_predicted_benefit_performance(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


### Pre 2020-10-14 ----

### Exact matching only on best drug
overall_benefit_calibration_rank1_pre_2020_orig <- overall_predicted_benefit_performance(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_calibration_rank1_pre_2020_orig_match_sex <- overall_predicted_benefit_performance(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_calibration_rank1_pre_2020_orig_match_sex_hba1c <- overall_predicted_benefit_performance(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)



## Summary ----

calibration_summary <- do.call(rbind, overall_benefit_calibration_tolerance3_post_2020_orig) %>%
  mutate(Matching = "Best drug") %>%
  rbind(
    do.call(rbind, overall_benefit_calibration_tolerance3_post_2020_orig_match_sex) %>%
      mutate(Matching = "Best drug, sex"),
    do.call(rbind, overall_benefit_calibration_tolerance3_post_2020_orig_match_sex_hba1c) %>%
      mutate(Matching = "Best drug, sex, hba1c")
  ) %>%
  mutate(Dataset = "Post 2020") %>%
  rbind(
    do.call(rbind, overall_benefit_calibration_tolerance3_pre_2020_orig) %>%
      mutate(Matching = "Best drug") %>%
      rbind(
        do.call(rbind, overall_benefit_calibration_tolerance3_pre_2020_orig_match_sex) %>%
          mutate(Matching = "Best drug, sex"),
        do.call(rbind, overall_benefit_calibration_tolerance3_pre_2020_orig_match_sex_hba1c) %>%
          mutate(Matching = "Best drug, sex, hba1c")
      ) %>%
      mutate(Dataset = "Pre 2020")
  ) %>%
  mutate(Method = "3mmol") %>%
  rbind(
    do.call(rbind, overall_benefit_calibration_rank1_post_2020_orig) %>%
      mutate(Matching = "Best drug") %>%
      rbind(
        do.call(rbind, overall_benefit_calibration_rank1_post_2020_orig_match_sex) %>%
          mutate(Matching = "Best drug, sex"),
        do.call(rbind, overall_benefit_calibration_rank1_post_2020_orig_match_sex_hba1c) %>%
          mutate(Matching = "Best drug, sex, hba1c")
      ) %>%
      mutate(Dataset = "Post 2020") %>%
      rbind(
        do.call(rbind, overall_benefit_calibration_rank1_pre_2020_orig) %>%
          mutate(Matching = "Best drug") %>%
          rbind(
            do.call(rbind, overall_benefit_calibration_rank1_pre_2020_orig_match_sex) %>%
              mutate(Matching = "Best drug, sex"),
            do.call(rbind, overall_benefit_calibration_rank1_pre_2020_orig_match_sex_hba1c) %>%
              mutate(Matching = "Best drug, sex, hba1c")
          ) %>%
          mutate(Dataset = "Pre 2020")
      ) %>%
      mutate(Method = "rank1")
  ) %>%
  rownames_to_column()

saveRDS(calibration_summary, "Outputs/CPRD/05.overvall_benefit_calibration_summary.rds")

# Overall calibration (plot benefit) #########################################

## Tolerance 3 mmol/mol ----

### Post 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_tolerance3_post_2020_orig <- overall_predicted_benefit(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_tolerance3_post_2020_orig_match_sex <- overall_predicted_benefit(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_tolerance3_post_2020_orig_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

### Pre 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_tolerance3_pre_2020_orig <- overall_predicted_benefit(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_tolerance3_pre_2020_orig_match_sex <- overall_predicted_benefit(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_tolerance3_pre_2020_orig_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

## Rank 1 ----

### Post 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_rank1_post_2020_orig <- overall_predicted_benefit(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_rank1_post_2020_orig_match_sex <- overall_predicted_benefit(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_rank1_post_2020_orig_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


### Pre 2020-10-14 ----

### Exact matching only on best drug
overall_benefit_rank1_pre_2020_orig <- overall_predicted_benefit(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_rank1_pre_2020_orig_match_sex <- overall_predicted_benefit(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_rank1_pre_2020_orig_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)




# Plots ####################################################

## Best drug combinations ----

optimal_drug_comparison_plot <- function(data, groups) {
  
  summary_number <- data.frame(drugs = data) %>%
    mutate(drug_count = if_else(drugs == "", 0L, str_count(drugs, ",") + 1)) %>%
    filter(drug_count != 0) %>%
    count(drug = drug_count, name = "Count") %>%
    mutate(
      Percentage = round(Count / length(data), 4)
    )
  
  summary_combinations <- purrr::map_dfr(
    names(groups),
    function(label) {
      tibble(drugs = data) %>%
        mutate(drug_count = if_else(drugs == "", 0L, str_count(data, ",") + 1)) %>%
        filter(drug_count %in% groups[[label]]) %>%
        count(drug = drugs, name = "Count") %>%
        mutate(
          Percentage = round(Count / length(data), 4),
          Combinations = label
        )
    }
  ) %>%
    group_by(Combinations) %>%
    arrange(Combinations, desc(Percentage), .by_group = TRUE) %>%
    mutate(Group_Percentage = round(sum(Count) / length(data), 4)) %>%
    ungroup() %>%
    filter(Percentage >= 0.001)
  
  plot_a <- summary_number %>%
    ggplot(aes(x = reorder(drug, rev(drug)), y = Percentage)) +
    geom_text(aes(label = drug), 
              hjust = -0.5,  # place label just outside bars to the right
              size = 8) +
    geom_col(fill = "#076fa2") +
    geom_hline(aes(yintercept = 0)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, max(summary_number$Percentage) + 0.05, by = 0.05)) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Number of predicted optimal drug", subtitle = "Proportion of individuals optimal drugs (defined by >3 mmol/mol predicted benefit)") + 
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank()
    )
  
  # Set the order of Combinations by descending Group_Percentage
  combi_order <- names(groups)
  
  # Prepare data with factor levels for ordering
  df2 <- summary_combinations %>%
    mutate(
      Combinations = factor(Combinations, levels = combi_order),
    ) %>%
    group_by(Combinations) %>%
    arrange(desc(Percentage), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(
      drug = factor(drug, levels = rev(unique(drug)))  # will reorder after reordering y axis below
    )
  
  # Now reorder drug factor within groups by Percentage, preserving group order
  # This is the trickier part:
  
  df2 <- df2 %>%
    arrange(Combinations, desc(Percentage)) %>%
    mutate(drug = factor(drug, levels = rev(drug)))
  
  # Plot
  plot_b <- ggplot(df2, aes(x = Percentage, y = drug)) +
    geom_text(aes(label = scales::percent(Percentage, accuracy = 0.1)), 
              hjust = -0.1,  # place label just outside bars to the right
              size = 3) +
    geom_col(fill = "#076fa2") +
    labs(
      title = "Predicted HbA1c-optimal drug classes"
    ) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, max(df2$Percentage) + 0.05, by = 0.05), expand = expansion(mult = c(0, 0.1))) +
    theme_minimal() +
    theme(
      axis.title = element_blank()
    )
  
  
  final_plot <- patchwork::wrap_plots(plot_a, plot_b, ncol = 2, nrow = 1)
  
  
  return(final_plot)
}


# Define groupings
groups <- list(
  "1-drug combination" = 1,
  "2-drug combinations" = 2,
  "3-drug combinations" = 3,
  "3/4/5-drug combinations" = 4:5
)


plot_tolerance3_analysis_post_2020 <- optimal_drug_comparison_plot(analysis_post_2020$pred.orig.within_3_of_best_drug_name, groups)

plot_tolerance3_analysis_pre_2020 <- optimal_drug_comparison_plot(analysis_pre_2020$pred.orig.within_3_of_best_drug_name, groups)



## Tolerance 3 mmol/mol ----

### Post 2020-10-14 ----
### Overall calibration
plot_overall_benefit_tolerance3_post_2020_orig <- overall_benefit_tolerance3_post_2020_orig %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (tolerance 3mmol) (matching best drug)")

plot_overall_benefit_tolerance3_post_2020_orig_match_sex <- overall_benefit_tolerance3_post_2020_orig_match_sex %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (tolerance 3mmol) (matching best drug / sex)")

plot_overall_benefit_tolerance3_post_2020_orig_match_sex_hba1c <- overall_benefit_tolerance3_post_2020_orig_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (tolerance 3mmol) (matching best drug / sex / hba1c_10)")


### Pre 2020-10-14 ----
### Overall calibration
plot_overall_benefit_tolerance3_pre_2020_orig <- overall_benefit_tolerance3_pre_2020_orig %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (tolerance 3mmol) (matching best drug)")

plot_overall_benefit_tolerance3_pre_2020_orig_match_sex <- overall_benefit_tolerance3_pre_2020_orig_match_sex %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (tolerance 3mmol) (matching best drug / sex)")

plot_overall_benefit_tolerance3_pre_2020_orig_match_sex_hba1c <- overall_benefit_tolerance3_pre_2020_orig_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (tolerance 3mmol) (matching best drug / sex / hba1c_10)")



## Rank 1 ----

### Post 2020-10-14 ----
### Overall calibration
plot_overall_benefit_rank1_post_2020_orig <- overall_benefit_rank1_post_2020_orig %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (matching best drug)")

plot_overall_benefit_rank1_post_2020_orig_match_sex <- overall_benefit_rank1_post_2020_orig_match_sex %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (matching best drug / sex)")

plot_overall_benefit_rank1_post_2020_orig_match_sex_hba1c <- overall_benefit_rank1_post_2020_orig_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (matching best drug / sex / hba1c_10)")


### Pre 2020-10-14 ----
### Overall calibration
plot_overall_benefit_rank1_pre_2020_orig <- overall_benefit_rank1_pre_2020_orig %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (matching best drug)")

plot_overall_benefit_rank1_pre_2020_orig_match_sex <- overall_benefit_rank1_pre_2020_orig_match_sex %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (matching best drug / sex)")

plot_overall_benefit_rank1_pre_2020_orig_match_sex_hba1c <- overall_benefit_rank1_pre_2020_orig_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (matching best drug / sex / hba1c_10)")


## PDF ----
pdf("Outputs/CPRD/05.calibration_predictions.pdf", width = 8, height = 5)
plot_pred_response_analysis_post_2020_preclosed
plot_pred_response_analysis_post_2020
plot_pred_response_analysis_pre_2020
patchwork::wrap_plots(
  plot_pred_response_analysis_post_2020_preclosed_glp1 +
    theme(legend.position = "bottom"), 
  plot_pred_response_analysis_post_2020_glp1 +
    theme(legend.position = "bottom"), 
  plot_pred_response_analysis_pre_2020_glp1 +
    theme(legend.position = "bottom")
) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
dev.off()


pdf("Outputs/CPRD/05.pre_drug_combinations.pdf", width = 13, height = 5)
plot_tolerance3_analysis_pre_2020
dev.off()

pdf("Outputs/CPRD/05.post_drug_combinations.pdf", width = 13, height = 5)
plot_tolerance3_analysis_post_2020
dev.off()

pdf("Outputs/CPRD/05.post_overall_calibration.pdf", width = 7, height = 5)
plot_overall_benefit_rank1_post_2020_orig
plot_overall_benefit_rank1_post_2020_orig_match_sex
plot_overall_benefit_rank1_post_2020_orig_match_sex_hba1c
plot_overall_benefit_tolerance3_post_2020_orig
plot_overall_benefit_tolerance3_post_2020_orig_match_sex
plot_overall_benefit_tolerance3_post_2020_orig_match_sex_hba1c
dev.off()

pdf("Outputs/CPRD/05.pre_overall_calibration.pdf", width = 7, height = 5)
plot_overall_benefit_rank1_pre_2020_orig
plot_overall_benefit_rank1_pre_2020_orig_match_sex
plot_overall_benefit_rank1_pre_2020_orig_match_sex_hba1c
plot_overall_benefit_tolerance3_pre_2020_orig
plot_overall_benefit_tolerance3_pre_2020_orig_match_sex
plot_overall_benefit_tolerance3_pre_2020_orig_match_sex_hba1c
dev.off()

pdf("Outputs/CPRD/05.pred_vs_obs_post_2020.pdf", width = 8, height = 5)
plot_pred_response_analysis_post_2020_comparison
dev.off()





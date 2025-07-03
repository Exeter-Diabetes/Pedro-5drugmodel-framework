
# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)
library(rms)
library(qpdf)

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
analysis_semaglutide_raw <- analysis_semaglutide_raw %>%
  analysis$cached("analysis_injectable_semaglutide") %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

analysis_oral_semaglutide_raw <- analysis_oral_semaglutide_raw %>%
  analysis$cached("analysis_oral_semaglutide") %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

analysis_post_2020_raw <- analysis_post_2020_raw %>%
  analysis$cached("analysis_post_2020") %>%
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

## Semaglutide ----
analysis_semaglutide <- analysis_semaglutide_raw %>%
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


## Oral semaglutide ----
analysis_oral_semaglutide <- analysis_oral_semaglutide_raw %>%
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

## Semaglutide ----
### By mice
analysis_semaglutide_mice_imputation <- imputation_methods(data = analysis_semaglutide,
                                                         method = "mice",
                                                         mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))


## Oral emaglutide ----
### By mice
analysis_oral_semaglutide_mice_imputation <- imputation_methods(data = analysis_oral_semaglutide,
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

### merge impute columns into main dataset
analysis_post_2020 <- analysis_post_2020 %>%
  cbind(
    analysis_post_2020_prediction_orig %>%
      select(contains("pred.orig")) %>%
      rename_with(~ str_replace(., "pred\\.orig\\.", "pred.orig.preclosed."))
  )


## Semaglutide ----
### Imputation by mice
analysis_semaglutide_prediction_mice <- predict_5drugmodel(analysis_semaglutide_mice_imputation,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))


### merge impute columns into main dataset
analysis_semaglutide <- analysis_semaglutide %>%
  cbind(
    analysis_semaglutide_prediction_mice %>%
      select(contains("mice_impute")),
    analysis_semaglutide_prediction_mice %>%
      select(contains("pred.mice")) %>%
      rename_with(~ str_replace(., "pred\\.mice\\.", "pred.mice.preclosed."))
  )


## Oral semaglutide ----
### Imputation by mice
analysis_oral_semaglutide_prediction_mice <- predict_5drugmodel(analysis_oral_semaglutide_mice_imputation,
                                                           model = m1.5.final,
                                                           drug_var = "drugclass",
                                                           drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))


### merge impute columns into main dataset
analysis_oral_semaglutide <- analysis_oral_semaglutide %>%
  cbind(
    analysis_oral_semaglutide_prediction_mice %>%
      select(contains("mice_impute")),
    analysis_oral_semaglutide_prediction_mice %>%
      select(contains("pred.mice")) %>%
      rename_with(~ str_replace(., "pred\\.mice\\.", "pred.mice.preclosed."))
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


## Semaglutide ----

### mice variables
#### GLP1
closed_loop_test_results_GLP1_semaglutide_mice <- closedtest_continuous_function(
  cohort = "Semaglutide subcohort",
  dataset = analysis_semaglutide %>%
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

## Adjustment used: https://doi.org/10.1016/j.diabres.2021.108904
## Difference between Fig5 Dulaglutide vs Semaglutide 0.48% (0.17,0.8)
injectable_semaglutide_hba1c_percent <- 0.48
injectable_semaglutide_hba1c_mmol <- injectable_semaglutide_hba1c_percent * 10.929

## Oral semaglutide ----

### mice variables
#### GLP1
closed_loop_test_results_GLP1_oral_semaglutide_mice <- closedtest_continuous_function(
  cohort = "Semaglutide subcohort",
  dataset = analysis_oral_semaglutide %>%
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



## Make predictions ----
analysis_post_2020 <- analysis_post_2020 %>%
  mutate(
    # original variables
    pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "SGLT2")),
    pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "GLP1")),
    pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "DPP4")),
    pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "TZD")),
    pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "SU")),
    pred.Sema = pred.orig.preclosed.GLP1 - injectable_semaglutide_hba1c_mmol,
    pred.Oral = predict_with_modelchoice_function(closed_loop_test_results_GLP1_oral_semaglutide_mice, analysis_post_2020 %>% mutate(drugclass = "GLP1"))
  )

analysis_semaglutide <- analysis_semaglutide %>%
  mutate(
    pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_orig, analysis_semaglutide %>% mutate(
      drugclass = "SGLT2",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_orig, analysis_semaglutide %>% mutate(
      drugclass = "GLP1",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_orig, analysis_semaglutide %>% mutate(
      drugclass = "DPP4",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_orig, analysis_semaglutide %>% mutate(
      drugclass = "TZD",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_orig, analysis_semaglutide %>% mutate(
      drugclass = "SU",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.Sema = pred.mice.preclosed.GLP1 - injectable_semaglutide_hba1c_mmol,
    pred.Oral = predict_with_modelchoice_function(closed_loop_test_results_GLP1_oral_semaglutide_mice, analysis_semaglutide %>% mutate(
      drugclass = "GLP1",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ))
  )


analysis_oral_semaglutide <- analysis_oral_semaglutide %>%
  mutate(
    pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_orig, analysis_oral_semaglutide %>% mutate(
      drugclass = "SGLT2",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_orig, analysis_oral_semaglutide %>% mutate(
      drugclass = "GLP1",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_orig, analysis_oral_semaglutide %>% mutate(
      drugclass = "DPP4",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_orig, analysis_oral_semaglutide %>% mutate(
      drugclass = "TZD",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_orig, analysis_oral_semaglutide %>% mutate(
      drugclass = "SU",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    )),
    pred.Sema = pred.mice.preclosed.GLP1 - injectable_semaglutide_hba1c_mmol,
    pred.Oral = predict_with_modelchoice_function(closed_loop_test_results_GLP1_oral_semaglutide_mice, analysis_oral_semaglutide %>% mutate(
      drugclass = "GLP1",
      prebmi = prebmi_mice_impute,
      preegfr = preegfr_mice_impute,
      pretotalcholesterol = pretotalcholesterol_mice_impute,
      prehdl = prehdl_mice_impute,
      prealt = prealt_mice_impute
    ))
  )


# Calibration plots ########################################

plot_calibration_semaglutide <- analysis_semaglutide %>%
  select(posthba1cfinal, pred.Sema, pred.GLP1, pred.mice.preclosed.GLP1) %>%
  rename("GLP1" = pred.GLP1, "Semaglutide" = pred.Sema, "Original" = pred.mice.preclosed.GLP1) %>%
  gather(key = "Predictions", value = "pred", -posthba1cfinal) %>%
  rename("obs" = posthba1cfinal) %>%
  ggplot(aes(x = pred, y = obs, colour = Predictions)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth() +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c (mmol/mol)") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

plot_calibration_oral_semaglutide <- analysis_oral_semaglutide %>%
  select(posthba1cfinal, pred.Oral, pred.GLP1, pred.mice.preclosed.GLP1) %>%
  rename("GLP1" = pred.GLP1, "Oral semaglutide" = pred.Oral, "Original" = pred.mice.preclosed.GLP1) %>%
  gather(key = "Predictions", value = "pred", -posthba1cfinal) %>%
  rename("obs" = posthba1cfinal) %>%
  ggplot(aes(x = pred, y = obs, colour = Predictions)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth() +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c (mmol/mol)") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Unified validation ########################################
interim_dataset <- analysis_post_2020 %>%
  rename("drugclass_old" = "drugclass") %>%
  mutate(drugclass = drugclass_old) %>%
  select(-matches("preclosed")) %>%
  rbind(
    analysis_semaglutide %>%
      rename("drugclass_old" = "drugclass") %>%
      mutate(drugclass = "Semaglutide") %>%
      mutate(
        prebmi = prebmi_mice_impute,
        preegfr = preegfr_mice_impute,
        pretotalcholesterol = pretotalcholesterol_mice_impute,
        prehdl = prehdl_mice_impute,
        prealt = prealt_mice_impute
      ) %>%
      select(-matches("mice")),
    analysis_oral_semaglutide %>%
      rename("drugclass_old" = "drugclass") %>%
      mutate(drugclass = "Oral semaglutide") %>%
      mutate(
        prebmi = prebmi_mice_impute,
        preegfr = preegfr_mice_impute,
        pretotalcholesterol = pretotalcholesterol_mice_impute,
        prehdl = prehdl_mice_impute,
        prealt = prealt_mice_impute
      ) %>%
      select(-matches("mice"))
  )

analysis_calibration <- unified_validation(
  data = interim_dataset, 
  drug_var = "drugclass",
  drugs = c("Semaglutide", "SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.Sema", "pred.SGLT2", "pred.GLP1", "pred.TZD", "pred.SU", "pred.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

analysis_semaglutide_calibration <- unified_validation(
  data = interim_dataset, 
  drug_var = "drugclass",
  drugs = c("Oral semaglutide", "Semaglutide", "SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.Oral", "pred.Sema", "pred.SGLT2", "pred.GLP1", "pred.TZD", "pred.SU", "pred.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

final_comparison <- analysis_calibration %>%
  filter(!(drug1 == "Semaglutide" & drug2 == "GLP1")) %>%
  filter(drug1 == "Semaglutide") %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) 


plot_unified_calibration <- final_comparison %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 3) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Linear adjustment")


final_semaglutide_comparison <- analysis_semaglutide_calibration %>%
  filter(!(drug1 == "Semaglutide" & drug2 == "GLP1") &
           !(drug1 == "Oral semaglutide" & drug2 == "Semaglutide") &
           !(drug1 == "Oral semaglutide" & drug2 == "GLP1")) %>%
  filter(str_detect(drug1, "Oral semaglutide")) %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) 


plot_unified_semaglutide_calibration <- final_semaglutide_comparison %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 3) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Linear adjustment")



# PDFs ########################################

pdf("plot_calibration_semaglutide.pdf", width = 6, height = 4)
plot_calibration_semaglutide
dev.off()

pdf("plot_unified_calibration_semaglutide.pdf", width = 6, height = 6)
plot_unified_calibration
dev.off()

# combine plots
pdf_files <- c("plot_calibration_semaglutide.pdf", "plot_unified_calibration_semaglutide.pdf")

# Combine PDFs
pdf_combine(input = pdf_files, output = "Outputs/CPRD/08.adding_injectable_semaglutide.pdf")
file.remove(pdf_files)


pdf("plot_calibration_semaglutide.pdf", width = 6, height = 4)
plot_calibration_oral_semaglutide
dev.off()

pdf("plot_unified_calibration_semaglutide.pdf", width = 6, height = 6)
plot_unified_semaglutide_calibration
dev.off()

# combine plots
pdf_files <- c("plot_calibration_semaglutide.pdf", "plot_unified_calibration_semaglutide.pdf")

# Combine PDFs
pdf_combine(input = pdf_files, output = "Outputs/CPRD/08.adding_oral_semaglutide.pdf")
file.remove(pdf_files)




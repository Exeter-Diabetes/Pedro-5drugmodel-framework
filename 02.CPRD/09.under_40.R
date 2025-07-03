
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
analysis_pre_2020_raw <- analysis_pre_2020_raw %>%
  analysis$cached("analysis_pre_2020") %>%
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

### By mice
analysis_post_2020_imputation <- imputation_methods(data = analysis_post_2020,
                                                           method = "mice",
                                                           mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))


# Predictions from datasets #########################################

### Original variables
analysis_post_2020_prediction <- predict_5drugmodel(analysis_post_2020_imputation,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))

### merge impute columns into main dataset
analysis_post_2020 <- analysis_post_2020 %>%
  select(-matches(
    analysis_post_2020_prediction %>%
      select(contains("mice_impute")) %>%
      rename_with(~ str_replace(., "_mice_impute", "")) %>%
      colnames()
  )) %>%
  cbind(
    analysis_post_2020_prediction %>%
      select(contains("mice_impute")) %>%
      rename_with(~ str_replace(., "_mice_impute", "")),
    analysis_post_2020_prediction %>%
      select(contains("pred.mice")) %>%
      rename_with(~ str_replace(., "pred\\.mice\\.", "pred.preclosed."))
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



## Make predictions ----
analysis_post_2020 <- analysis_post_2020 %>%
  mutate(
    # original variables
    pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "SGLT2")),
    pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "GLP1")),
    pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "DPP4")),
    pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "TZD")),
    pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "SU"))
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




analysis_cohort <- analysis_post_2020 %>%
  rename_with(~ str_replace(., ".orig", "")) %>%
  select(-contains("preclosed")) %>%
  rbind(
    analysis_pre_2020 %>%
      rename_with(~ str_replace(., ".orig", "")) %>%
      select(-contains("preclosed"))
  ) %>%
  filter(agetx <= 40)

# Overall calibration (plot benefit) #########################################

## Tolerance 3 mmol/mol ----

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_tolerance3_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_cohort %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 5,
  pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  conc_tolerance = 3,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


## Rank 1 ----

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_rank1_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_cohort %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 5,
  pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


# Unified validation ########################################

analysis_cohort_calibration_adj <- unified_validation(
  data = analysis_cohort, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)


# PDFs ----
plot_overall_benefit_tolerance3_match_sex_hba1c <- overall_benefit_tolerance3_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Best drugs (combination) defined by values closest to best within 3 mmol/mol")

plot_overall_benefit_rank1_match_sex_hba1c <- overall_benefit_rank1_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Best drug defined by best predicted response")


plot_unified_validation <- analysis_cohort_calibration_adj %>%
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
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")



pdf("plot_overall_benefits.pdf", width = 7, height = 5)
plot_overall_benefit_tolerance3_match_sex_hba1c
plot_overall_benefit_rank1_match_sex_hba1c
dev.off()

pdf("plot_unified_validation.pdf", width = 12, height = 5)
plot_unified_validation
dev.off()

# combine plots
pdf_files <- c("plot_overall_benefits.pdf", "plot_unified_validation.pdf")

# Combine PDFs
pdf_combine(input = pdf_files, output = "Outputs/CPRD/09.under_40.pdf")
file.remove(pdf_files)




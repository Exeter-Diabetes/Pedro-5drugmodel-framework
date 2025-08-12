
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

analysis_dataset <- analysis_pre_2020_raw %>%
  rbind(analysis_post_2020_raw) %>%
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
analysis_mice_imputation <- imputation_methods(data = analysis_dataset,
                                                    method = "mice",
                                                    mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))



# Predictions from datasets #########################################

### Imputation by mice
analysis_prediction_mice <- predict_5drugmodel(analysis_mice_imputation,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"))


### merge impute columns into main dataset
analysis_dataset <- analysis_dataset %>%
  select(-c("prebmi", "preegfr", "pretotalcholesterol", "prehdl", "prealt")) %>%
  cbind(
    analysis_prediction_mice %>%
      select(contains("mice_impute")) %>%
      rename_with(~ str_replace(., "_mice_impute", "")),
    analysis_prediction_mice %>%
      select(contains("pred.mice")) %>%
      rename_with(~ str_replace(., "pred\\.mice\\.", "pred."))
  )



# Split cohort into different ethnicities #########################################

analysis_white <- analysis_dataset %>%
  filter(ethnicity == "White")

analysis_asian <- analysis_dataset %>%
  filter(ethnicity == "South Asian")

analysis_black <- analysis_dataset %>%
  filter(ethnicity == "Black")

analysis_other <- analysis_dataset %>%
  filter(ethnicity %in% c("Other", "Mixed", "Missing"))


# Closed loop test #################################################

## White ----

#### SGLT2
closed_loop_test_results_SGLT2_white <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_white %>% filter(drugclass == "SGLT2"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_white <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_white %>% filter(drugclass == "GLP1"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_white <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_white %>% filter(drugclass == "DPP4"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_white <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_white %>% filter(drugclass == "SU"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_white <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_white %>% filter(drugclass == "TZD"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)


## Asian ----

#### SGLT2
closed_loop_test_results_SGLT2_asian <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_asian %>% filter(drugclass == "SGLT2"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_asian <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_asian %>% filter(drugclass == "GLP1"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_asian <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_asian %>% filter(drugclass == "DPP4"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_asian <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_asian %>% filter(drugclass == "SU"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_asian <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_asian %>% filter(drugclass == "TZD"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)


## Black ----

#### SGLT2
closed_loop_test_results_SGLT2_black <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_black %>% filter(drugclass == "SGLT2"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_black <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_black %>% filter(drugclass == "GLP1"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_black <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_black %>% filter(drugclass == "DPP4"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_black <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_black %>% filter(drugclass == "SU"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_black <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_black %>% filter(drugclass == "TZD"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)


## Other ----

#### SGLT2
closed_loop_test_results_SGLT2_other <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = analysis_other %>% filter(drugclass == "SGLT2"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### GLP1
closed_loop_test_results_GLP1_other <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = analysis_other %>% filter(drugclass == "GLP1"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### DPP4
closed_loop_test_results_DPP4_other <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = analysis_other %>% filter(drugclass == "DPP4"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### SU
closed_loop_test_results_SU_other <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = analysis_other %>% filter(drugclass == "SU"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

#### TZD
closed_loop_test_results_TZD_other <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = analysis_other %>% filter(drugclass == "TZD"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)


## Summary table ----

summary_table <- closed_loop_test_results_SGLT2_white$testing_results %>%
  mutate(Treatment = "SGLT2", Ethnicity = "White") %>%
  rbind(
    closed_loop_test_results_GLP1_white$testing_results %>%
      mutate(Treatment = "GLP1", Ethnicity = "White"),
    closed_loop_test_results_DPP4_white$testing_results %>%
      mutate(Treatment = "DPP4", Ethnicity = "White"),
    closed_loop_test_results_TZD_white$testing_results %>%
      mutate(Treatment = "TZD", Ethnicity = "White"),
    closed_loop_test_results_SU_white$testing_results %>%
      mutate(Treatment = "SU", Ethnicity = "White"),
    closed_loop_test_results_SGLT2_black$testing_results %>%
      mutate(Treatment = "SGLT2", Ethnicity = "Black"),
    closed_loop_test_results_GLP1_black$testing_results %>%
      mutate(Treatment = "GLP1", Ethnicity = "Black"),
    closed_loop_test_results_DPP4_black$testing_results %>%
      mutate(Treatment = "DPP4", Ethnicity = "Black"),
    closed_loop_test_results_TZD_black$testing_results %>%
      mutate(Treatment = "TZD", Ethnicity = "Black"),
    closed_loop_test_results_SU_black$testing_results %>%
      mutate(Treatment = "SU", Ethnicity = "Black"),
    closed_loop_test_results_SGLT2_asian$testing_results %>%
      mutate(Treatment = "SGLT2", Ethnicity = "South Asian"),
    closed_loop_test_results_GLP1_asian$testing_results %>%
      mutate(Treatment = "GLP1", Ethnicity = "South Asian"),
    closed_loop_test_results_DPP4_asian$testing_results %>%
      mutate(Treatment = "DPP4", Ethnicity = "South Asian"),
    closed_loop_test_results_TZD_asian$testing_results %>%
      mutate(Treatment = "TZD", Ethnicity = "South Asian"),
    closed_loop_test_results_SU_asian$testing_results %>%
      mutate(Treatment = "SU", Ethnicity = "South Asian"),
    closed_loop_test_results_SGLT2_other$testing_results %>%
      mutate(Treatment = "SGLT2", Ethnicity = "Other"),
    closed_loop_test_results_GLP1_other$testing_results %>%
      mutate(Treatment = "GLP1", Ethnicity = "Other"),
    closed_loop_test_results_DPP4_other$testing_results %>%
      mutate(Treatment = "DPP4", Ethnicity = "Other"),
    closed_loop_test_results_TZD_other$testing_results %>%
      mutate(Treatment = "TZD", Ethnicity = "Other"),
    closed_loop_test_results_SU_other$testing_results %>%
      mutate(Treatment = "SU", Ethnicity = "Other")
  ) %>%
  mutate(
    Ethnicity = factor(Ethnicity, levels = c("White" , "South Asian", "Black", "Other")),
    Treatment = factor(Treatment, levels = c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
    model = factor(model, levels = c("Original", "Recalibrated intercept", "Recalibrated"))
  ) %>%
  rename("Intercept" = "intercept_with_CI", "Model" = "model", "N" = "n_population", "Slope" = "slope_with_CI") %>%
  select(Treatment, Ethnicity, N, Model, Intercept, Slope) %>%
  arrange(Treatment, Ethnicity, Model)

plot_intercepts <- summary_table %>%
  extract(Intercept, into = c("coef", "coef_low", "coef_high"),
          regex = "(-?[\\d.]+) \\((-?[\\d.]+), (-?[\\d.]+)\\)", convert = TRUE) %>%
  filter(Model == "Recalibrated intercept") %>%
  mutate(
    cross_zero = if_else(coef_low < 0 & coef_high > 0, NA_character_, Treatment)
  ) %>%
  ggplot(aes(x = coef, xmin = coef_low, xmax = coef_high, y = Ethnicity, colour = cross_zero)) +
  geom_vline(aes(xintercept = 0), colour = "black") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~Treatment, nrow = 1) +
  labs(x = "Intercept adjustment") +
  theme_bw() +
  theme(
    legend.title = element_blank()
  )

plot_recalibration <- summary_table %>%
  extract(Intercept, into = c("coef", "coef_low", "coef_high"),
          regex = "(-?[\\d.]+) \\((-?[\\d.]+), (-?[\\d.]+)\\)", convert = TRUE) %>%
  filter(Model == "Recalibrated") %>%
  mutate(
    cross_zero = if_else(coef_low < 0 & coef_high > 0, NA_character_, Treatment),
    Type = "Intercept"
  ) %>%
  select(-Slope) %>%
  rbind(
    summary_table %>%
      extract(Slope, into = c("coef", "coef_low", "coef_high"),
              regex = "(-?[\\d.]+) \\((-?[\\d.]+), (-?[\\d.]+)\\)", convert = TRUE) %>%
      filter(Model == "Recalibrated") %>%
      mutate(
        cross_zero = if_else(coef_low < 1 & coef_high > 1, NA_character_, Treatment),
        Type = "Slope"
      ) %>%
      select(-Intercept)
  ) %>%
  ggplot(aes(x = coef, xmin = coef_low, xmax = coef_high, y = Ethnicity, colour = cross_zero)) +
  geom_vline(data = data.frame(intercept = c(0, 1), Type = c("Intercept", "Slope")), aes(xintercept = intercept), colour = "black") +
  geom_point() +
  geom_errorbar() +
  facet_grid(Type~Treatment, scales = "free") +
  labs(x = "Intercept adjustment") +
  theme_bw() +
  theme(
    legend.title = element_blank()
  )


# Overall calibration #########################################

## White ----
overall_benefit_white <- cateval::compute_overall_benefit(
  data = analysis_white %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


## South Asian ----
overall_benefit_asian <- cateval::compute_overall_benefit(
  data = analysis_asian %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

## Black ----
overall_benefit_black <- cateval::compute_overall_benefit(
  data = analysis_black %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 5,
  pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

## Other/Mixed/Missing ----
overall_benefit_other <- cateval::compute_overall_benefit(
  data = analysis_other %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 5,
  pred_cols = paste0("pred.", c("DPP4", "GLP1", "SGLT2", "SU", "TZD")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


# Unified validation ########################################

## White ----
unified_calibration_white <- cateval::unified_validation(
  data = analysis_white, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)


## South Asian ----
unified_calibration_asian <- cateval::unified_validation(
  data = analysis_asian, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

## Black ----
unified_calibration_black <- cateval::unified_validation(
  data = analysis_black, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)


## Other/Mixed/Missing ----
unified_calibration_other <- cateval::unified_validation(
  data = analysis_other, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)



# PDFs ----

## Overall benefit ----
plot_overall_benefit <- NULL
plot_overall_benefit[["white"]] <- overall_benefit_white %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: White")

plot_overall_benefit[["asian"]] <- overall_benefit_asian %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: South Asian")

plot_overall_benefit[["black"]] <- overall_benefit_black %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: Black")

plot_overall_benefit[["other"]] <- overall_benefit_other %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  # geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: Other / Mixed / Missing")


pdf("Outputs/CPRD/11.plot_overall_benefits.pdf", width = 7, height = 5)
plot_overall_benefit
dev.off()



## Unified validation ----
plot_unified_validation <- NULL

plot_unified_validation[["white"]] <- unified_calibration_white %>%
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
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: White (at least 100 individuals on either drug)")

plot_unified_validation[["asian"]] <- unified_calibration_asian %>%
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
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: South Asian (at least 100 individuals on either drug)")

plot_unified_validation[["black"]] <- unified_calibration_black %>%
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
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: Black (at least 100 individuals on either drug)")

plot_unified_validation[["other"]] <- unified_calibration_other %>%
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
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Ethnicity: Other / Mixed / Missing (at least 100 individuals on either drug)")


pdf("Outputs/CPRD/11.plot_unified_validation.pdf", width = 12, height = 5)
plot_unified_validation
dev.off()


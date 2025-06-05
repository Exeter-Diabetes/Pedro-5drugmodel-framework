
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

sapply(paste0("01.Functions/", list.files("01.Functions")), source)
source("02.CPRD/02.impute_missingness.R")
source("02.CPRD/03.model_predictions.R")

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
      select(contains("pred.orig")),
    analysis_post_2020_prediction_group %>%
      select(contains("group_impute")),
    analysis_post_2020_prediction_group %>%
      select(contains("pred.group")),
    analysis_post_2020_prediction_mice %>%
      select(contains("mice_impute")),
    analysis_post_2020_prediction_mice %>%
      select(contains("pred.mice"))
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
      select(contains("pred.orig"))
  )



# Per drug pair calibration ########################################

## Worked example: SGLT2vsDPP4 ----

# Testing out a single drug pair: SGLT2 vs DPP4
test_dataset <- analysis_post_2020 %>%
  mutate(benefit = pred.orig.SGLT2 - pred.orig.DPP4) # benefit variable

## Linear adjustment
analysis_SGLT2_DPP4_orig_calibration_adj <- heterogenous_effect_calibration(
  data = test_dataset, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "DPP4"),
  benefit_var = "benefit",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

## Matching
analysis_SGLT2_DPP4_orig_calibration_matching <- heterogenous_effect_calibration(
  data = test_dataset %>% mutate(hba1c_group = ntile(prehba1c, 10)), 
  drug_var = "drugclass",
  drugs = c("SGLT2", "DPP4"),
  benefit_var = "benefit",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  matching = TRUE,
  adjustment_var = NULL,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

## Matching + Linear adjustment
analysis_SGLT2_DPP4_orig_calibration_matching_adj <- heterogenous_effect_calibration(
  data = test_dataset %>% mutate(hba1c_group = ntile(prehba1c, 10)), 
  drug_var = "drugclass",
  drugs = c("SGLT2", "DPP4"),
  benefit_var = "benefit",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  matching = TRUE,
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)



# Unified validation ########################################

## Post 2020-10-14 ----

### Testing imputation methods ----
# Original
analysis_post_2020_orig_calibration <- unified_validation(
  data = analysis_post_2020, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

# Group imputation
analysis_post_2020_group_calibration <- unified_validation(
  data = analysis_post_2020, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.group.SGLT2", "pred.group.GLP1", "pred.group.TZD", "pred.group.SU", "pred.group.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  adjustment_var = c("t2dmduration", "prebmi_group_impute", "prehba1c", "agetx", "prealt_group_impute", "preegfr_group_impute", "pretotalcholesterol_group_impute", "prehdl_group_impute", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

# Mice imputation
analysis_post_2020_mice_calibration <- unified_validation(
  data = analysis_post_2020, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.mice.SGLT2", "pred.mice.GLP1", "pred.mice.TZD", "pred.mice.SU", "pred.mice.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  adjustment_var = c("t2dmduration", "prebmi_mice_impute", "prehba1c", "agetx", "prealt_mice_impute", "preegfr_mice_impute", "pretotalcholesterol_mice_impute", "prehdl_mice_impute", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)


### All drug classes ----
analysis_post_2020_orig_calibration_adj <- unified_validation(
  data = analysis_post_2020, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

analysis_post_2020_orig_calibration_matching <- unified_validation(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)), 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  matching = TRUE,
  adjustment_var = NULL,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

analysis_post_2020_orig_calibration_matching_adj <- unified_validation(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)), 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  matching = TRUE,
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


## Pre 2020-10-14 ----

### All drug classes ----
analysis_pre_2020_orig_calibration_adj <- unified_validation(
  data = analysis_pre_2020, 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

analysis_pre_2020_orig_calibration_matching <- unified_validation(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)), 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  matching = TRUE,
  adjustment_var = NULL,
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)

analysis_pre_2020_orig_calibration_matching_adj <- unified_validation(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)), 
  drug_var = "drugclass",
  drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
  prediction_vars = c("pred.orig.SGLT2", "pred.orig.GLP1", "pred.orig.TZD", "pred.orig.SU", "pred.orig.DPP4"),
  outcome_var = "posthba1cfinal",
  cal_groups = c(3, 5, 10),
  matching = TRUE,
  adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)




# Plots ####################################################

## Testing imputation methods ----
pdf("Outputs/CPRD/06.drug_pair_calibration_imputation.pdf", width = 12, height = 5)
analysis_post_2020_orig_calibration %>%
  mutate(Imputation = "Complete") %>%
  rbind(
    analysis_post_2020_group_calibration %>%
      mutate(Imputation = "Group"),
    analysis_post_2020_mice_calibration %>%
      mutate(Imputation = "MICE")
  ) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Imputation)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Linear adjustment")
dev.off()



## Post 2020-10-14 ----
pdf("Outputs/CPRD/06.drug_pair_calibration_matching_post_2020.pdf", width = 12, height = 5)
analysis_post_2020_orig_calibration_adj %>%
  mutate(Method = "Adjustment") %>%
  rbind(
    analysis_post_2020_orig_calibration_matching %>%
      mutate(Method = "Matching"),
    analysis_post_2020_orig_calibration_matching_adj %>%
      mutate(Method = "Matching + Adj")
  ) %>%
  filter(n_groups == 10) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")
analysis_post_2020_orig_calibration_adj %>%
  mutate(Method = "Adjustment") %>%
  rbind(
    analysis_post_2020_orig_calibration_matching %>%
      mutate(Method = "Matching"),
    analysis_post_2020_orig_calibration_matching_adj %>%
      mutate(Method = "Matching + Adj")
  ) %>%
  filter(n_groups == 5) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")
analysis_post_2020_orig_calibration_adj %>%
  mutate(Method = "Adjustment") %>%
  rbind(
    analysis_post_2020_orig_calibration_matching %>%
      mutate(Method = "Matching"),
    analysis_post_2020_orig_calibration_matching_adj %>%
      mutate(Method = "Matching + Adj")
  ) %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(Method, drugcombo, n_groups) %>%
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
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")
dev.off()


## Pre 2020-10-14 ----
pdf("Outputs/CPRD/06.drug_pair_calibration_matching_pre_2020.pdf", width = 12, height = 5)
analysis_pre_2020_orig_calibration_adj %>%
  mutate(Method = "Adjustment") %>%
  rbind(
    analysis_pre_2020_orig_calibration_matching %>%
      mutate(Method = "Matching"),
    analysis_pre_2020_orig_calibration_matching_adj %>%
      mutate(Method = "Matching + Adj")
  ) %>%
  filter(n_groups == 10) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")
analysis_pre_2020_orig_calibration_adj %>%
  mutate(Method = "Adjustment") %>%
  rbind(
    analysis_pre_2020_orig_calibration_matching %>%
      mutate(Method = "Matching"),
    analysis_pre_2020_orig_calibration_matching_adj %>%
      mutate(Method = "Matching + Adj")
  ) %>%
  filter(n_groups == 5) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")
analysis_pre_2020_orig_calibration_adj %>%
  mutate(Method = "Adjustment") %>%
  rbind(
    analysis_pre_2020_orig_calibration_matching %>%
      mutate(Method = "Matching"),
    analysis_pre_2020_orig_calibration_matching_adj %>%
      mutate(Method = "Matching + Adj")
  ) %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(Method, drugcombo, n_groups) %>%
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
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~title, nrow = 2) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")
dev.off()












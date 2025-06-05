
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

# Calibration of Predictions ##########################################

## Post 2020-10-14 ----
plot_pred_response_analysis_post_2020 <- analysis_post_2020 %>%
  pivot_longer(cols = contains("pred.orig")) %>%
  select(name, value, drugclass, posthba1cfinal) %>%
  mutate(name = gsub("pred\\.orig\\.", "", name)) %>%
  filter(drugclass == name) %>%
  ggplot(aes(x = value, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  stat_smooth(method='lm', formula = y~poly(x,2)) +
  facet_wrap(~name) +
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Post 2020")
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
  labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c", title = "Post 2020")
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




# Overall calibration (benefit) #########################################

## Post 2020-10-14 ----
### Exact matching only on best drug
overall_benefit_post_2020_orig <- overall_predicted_benefit(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_post_2020_orig_match_sex <- overall_predicted_benefit(
  data = analysis_post_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_post_2020_orig_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_post_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)


## Pre 2020-10-14 ----

### Exact matching only on best drug
overall_benefit_pre_2020_orig <- overall_predicted_benefit(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

### Exact matching only on best drug / sex
overall_benefit_pre_2020_orig_match_sex <- overall_predicted_benefit(
  data = analysis_pre_2020,
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex")
)

### Exact matching only on best drug / sex / hba1c_10
overall_benefit_pre_2020_orig_match_sex_hba1c <- overall_predicted_benefit(
  data = analysis_pre_2020 %>% mutate(hba1c_group = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  outcome_var = "posthba1cfinal",
  cal_groups = 10,
  pred_cols = paste0("pred.orig.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
  matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_group")
)




# Plots ####################################################

## Post 2020-10-14 ----
### Overall calibration
plot_overall_benefit_post_2020_orig <- overall_benefit_post_2020_orig %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (matching best drug)")

plot_overall_benefit_post_2020_orig_match_sex <- overall_benefit_post_2020_orig_match_sex %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (matching best drug / sex)")

plot_overall_benefit_post_2020_orig_match_sex_hba1c <- overall_benefit_post_2020_orig_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Post 2020 (matching best drug / sex / hba1c_10)")


## Pre 2020-10-14 ----
### Overall calibration
plot_overall_benefit_pre_2020_orig <- overall_benefit_pre_2020_orig %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (matching best drug)")

plot_overall_benefit_pre_2020_orig_match_sex <- overall_benefit_pre_2020_orig_match_sex %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Pre 2020 (matching best drug / sex)")

plot_overall_benefit_pre_2020_orig_match_sex_hba1c <- overall_benefit_pre_2020_orig_match_sex_hba1c %>%
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
pdf("Outputs/CPRD/04.calibration_predictions.pdf", width = 8, height = 5)
plot_pred_response_analysis_post_2020
plot_pred_response_analysis_pre_2020
patchwork::wrap_plots(
  plot_pred_response_analysis_post_2020_glp1 +
    theme(legend.position = "bottom"), 
  plot_pred_response_analysis_pre_2020_glp1 +
    theme(legend.position = "bottom")
) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
dev.off()

pdf("Outputs/CPRD/04.overall_calibration.pdf", width = 7, height = 5)
plot_overall_benefit_post_2020_orig
plot_overall_benefit_post_2020_orig_match_sex
plot_overall_benefit_post_2020_orig_match_sex_hba1c
plot_overall_benefit_pre_2020_orig
plot_overall_benefit_pre_2020_orig_match_sex
plot_overall_benefit_pre_2020_orig_match_sex_hba1c
dev.off()








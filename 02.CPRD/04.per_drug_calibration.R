# load libraries
library(tidyverse)
library(aurum)


# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

# set up functions
is.integer64 <- function(x){
  class(x)=="integer64"
}

sapply(paste0("01.Functions/", list.files("01.Functions")), source)
source("02.CPRD/02.impute_missingness.R")
source("02.CPRD/03.model_predictions.R")

# set up dataset
analysis_post_2020_raw <- analysis_post_2020_raw %>%
  analysis$cached("analysis_post_2020") %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

# load model
load("fivedrugmodel_5knot_share_20230823.Rdata")

#:------------------------------------------------------
# Pre-processing datasets
#:-- Post 2020-10-14
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


#:------------------------------------------------------
# Missing data imputation
analysis_post_2020_group_imputation <- imputation_methods(data = analysis_post_2020,
                                                          method = "group")

analysis_post_2020_mice_imputation <- imputation_methods(data = analysis_post_2020,
                                                         method = "mice",
                                                         mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))

#:------------------------------------------------------
# Predictions from datasets
library(rms)
analysis_post_2020_prediction_orig <- predict_5drugmodel(analysis_post_2020,
                                                    model = m1.5.final,
                                                    drug_var = "drugclass",
                                                    drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"),
                                                    pred_col = "pred.orig.")


analysis_post_2020_prediction_group <- predict_5drugmodel(analysis_post_2020_group_imputation,
                                                         model = m1.5.final,
                                                         drug_var = "drugclass",
                                                         drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))

analysis_post_2020_prediction_mice <- predict_5drugmodel(analysis_post_2020_mice_imputation,
                                                    model = m1.5.final,
                                                    drug_var = "drugclass",
                                                    drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))



# merge impute columns into main dataset
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




#:------------------------------------------------------
# Per drug calibration

#:----
# Getting the best drug
## Rank 1
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

## Tolerance 3 mmol/mol
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


#:----
# Calibration by drug
## Simple method
analysis_orig_per_drug_cal_SGLT2 <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "SGLT2",
  pred_var = "pred.orig.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.orig.rank1_drug_name",
  cal_groups = 10,
  match.vars = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_orig_per_drug_cal_GLP1 <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "GLP1",
  pred_var = "pred.orig.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.orig.rank1_drug_name",
  cal_groups = 5,
  match.vars = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_orig_per_drug_cal_TZD <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "TZD",
  pred_var = "pred.orig.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.orig.rank1_drug_name",
  cal_groups = 1,
  match.vars = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_orig_per_drug_cal_SU <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "SU",
  pred_var = "pred.orig.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.orig.rank1_drug_name",
  cal_groups = 5,
  match.vars = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

## Grouped imputation model predictions
analysis_group_per_drug_cal_SGLT2 <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "SGLT2",
  pred_var = "pred.group.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.group.rank1_drug_name",
  cal_groups = 10,
  match.vars = c("t2dmduration", "prebmi_group_impute", "agetx", "prealt_group_impute", "preegfr_group_impute", "pretotalcholesterol_group_impute", "prehdl_group_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_group_per_drug_cal_GLP1 <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "GLP1",
  pred_var = "pred.group.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.group.rank1_drug_name",
  cal_groups = 5,
  match.vars = c("t2dmduration", "prebmi_group_impute", "agetx", "prealt_group_impute", "preegfr_group_impute", "pretotalcholesterol_group_impute", "prehdl_group_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_group_per_drug_cal_TZD <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "TZD",
  pred_var = "pred.group.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.group.rank1_drug_name",
  cal_groups = 1,
  match.vars = c("t2dmduration", "prebmi_group_impute", "agetx", "prealt_group_impute", "preegfr_group_impute", "pretotalcholesterol_group_impute", "prehdl_group_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_group_per_drug_cal_SU <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "SU",
  pred_var = "pred.group.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.group.rank1_drug_name",
  cal_groups = 5,
  match.vars = c("t2dmduration", "prebmi_group_impute", "agetx", "prealt_group_impute", "preegfr_group_impute", "pretotalcholesterol_group_impute", "prehdl_group_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

## MICE imputation model predictions
analysis_mice_per_drug_cal_SGLT2 <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "SGLT2",
  pred_var = "pred.mice.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.mice.rank1_drug_name",
  cal_groups = 10,
  match.vars = c("t2dmduration", "prebmi_mice_impute", "agetx", "prealt_mice_impute", "preegfr_mice_impute", "pretotalcholesterol_mice_impute", "prehdl_mice_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_mice_per_drug_cal_GLP1 <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "GLP1",
  pred_var = "pred.mice.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.mice.rank1_drug_name",
  cal_groups = 5,
  match.vars = c("t2dmduration", "prebmi_mice_impute", "agetx", "prealt_mice_impute", "preegfr_mice_impute", "pretotalcholesterol_mice_impute", "prehdl_mice_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_mice_per_drug_cal_TZD <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "TZD",
  pred_var = "pred.mice.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.mice.rank1_drug_name",
  cal_groups = 1,
  match.vars = c("t2dmduration", "prebmi_mice_impute", "agetx", "prealt_mice_impute", "preegfr_mice_impute", "pretotalcholesterol_mice_impute", "prehdl_mice_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)

analysis_mice_per_drug_cal_SU <- calibration_per_drug(
  data = analysis_post_2020 %>%
    mutate(hba1c_groups = ntile(prehba1c, 10)),
  drug_var = "drugclass",
  drug = "SU",
  pred_var = "pred.mice.",
  outcome_var = "posthba1cfinal",
  best_drug_var = "pred.mice.rank1_drug_name",
  cal_groups = 5,
  match.vars = c("t2dmduration", "prebmi_mice_impute", "agetx", "prealt_mice_impute", "preegfr_mice_impute", "pretotalcholesterol_mice_impute", "prehdl_mice_impute", "smoke", "imd5", "ncurrtx", "drugline"),
  match.exact = c("sex", "hba1c_groups")
)







#:------------------------
# Plotting these results

pdf(file = "Outputs/CPRD/04.per_drug_calibration.pdf")
analysis_orig_per_drug_cal_SGLT2 %>%
  mutate(Method = "Discard") %>%
  rbind(
    analysis_group_per_drug_cal_SGLT2 %>%
      mutate(Method = "Simple"),
    analysis_mice_per_drug_cal_SGLT2 %>%
      mutate(Method = "MICE")
  ) %>%
  mutate(`Drug class` = "SGLT2") %>%
  rbind(
    analysis_orig_per_drug_cal_GLP1 %>%
      mutate(Method = "Discard") %>%
      rbind(
        analysis_group_per_drug_cal_GLP1 %>%
          mutate(Method = "Simple"),
        analysis_mice_per_drug_cal_GLP1 %>%
          mutate(Method = "MICE")
      ) %>%
      mutate(`Drug class` = "GLP1"),
    analysis_orig_per_drug_cal_TZD %>%
      mutate(Method = "Discard") %>%
      rbind(
        analysis_group_per_drug_cal_TZD %>%
          mutate(Method = "Simple"),
        analysis_mice_per_drug_cal_TZD %>%
          mutate(Method = "MICE")
      ) %>%
      mutate(`Drug class` = "TZD"),
    analysis_orig_per_drug_cal_SU %>%
      mutate(Method = "Discard") %>%
      rbind(
        analysis_group_per_drug_cal_SU %>%
          mutate(Method = "Simple"),
        analysis_mice_per_drug_cal_SU %>%
          mutate(Method = "MICE")
      ) %>%
      mutate(`Drug class` = "SU")
  ) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
  geom_vline(aes(xintercept = 0), colour = "grey", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black") +
  geom_point(alpha = 0.7) +
  geom_errorbar(alpha = 0.7) +
  stat_smooth(method = "lm", alpha = 0.16) +
  theme_minimal() +
  facet_wrap(~`Drug class`) +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)")
dev.off()









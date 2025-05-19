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




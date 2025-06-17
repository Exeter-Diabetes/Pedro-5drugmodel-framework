
# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)
library(rms)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

# set up dataset
t2d_1stinstance <- t2d_1stinstance %>%
  analysis$cached("20250327_t2d_1stinstance")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

# Flow diagram: inclusion/exclusions ----

analysis = cprd$analysis("pedro_mm")
## Datasets
#:-- Post 2020-10-14
#:-- Pre 2020-10-14 (development cohort)


## Post 2020-10-14 ----
analysis_post_2020 <- t2d_1stinstance %>%
  ## Exclusions based on time - initiations after 14th October 2020
  filter(dstartdate > as.Date("2020-10-14")) %>%
  ## Exclusions based on type of drugsubstance: 
  ### SU - non-gliclazide
  ### TZD - rosiglitazone
  ### GLP1-RA - lixisenatide, exenatide slow-release, semaglutide
  ### MFN
  filter(
    drug_class != "MFN" & drug_class != "INS" & 
      drug_class != "Glinide"& drug_class != "Acarbose"
  ) %>%
  filter(
    drug_substance != "Glimepiride" & drug_substance != "Lixisenatide" &
      drug_substance != "Glipizide" & drug_substance != "Ertugliflozin" &
      drug_substance != "Glibenclamide" & drug_substance != "Tolbutamide" &
      drug_substance != "Low-dose semaglutide" & drug_substance != "Oral semaglutide" &
      drug_substance != "Semaglutide, dose unclear"
  ) %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <80
  filter(is.na(INS) | INS != 1) %>%
  filter(drugline_all != 1) %>%
  filter(!(preckdstage %in% c("stage 5"))) %>%
  filter(dstartdate_age > 18 & dstartdate_age < 80) %>%
  ## Exclusions
  ### Multiple treatments on same day
  ### Start within 61 days since last therapy
  ### missing HbA1c
  ### baseline HbA1c<53
  ### baseline HbA1c>120
  filter(multi_drug_start_class == 0) %>%
  filter(timeprevcombo_class >= 61) %>%
  filter(!is.na(prehba1c)) %>%
  filter(prehba1c >= 53) %>%
  filter(prehba1c <= 110) %>%
  ## Exclusions
  ### missing baseline (BMI, eGFR, ALT, HDL, total cholesterol)
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  # filter(!is.na(prebmi)) %>%
  # filter(!is.na(preegfr)) %>%
  # filter(!is.na(prealt)) %>%
  # filter(!is.na(prehdl)) %>%
  # filter(!is.na(pretotalcholesterol)) %>%
  # ## Exclusions
  # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  filter(!is.na(posthba1cfinal)) %>%
  analysis$cached("analysis_post_2020", indexes=c("patid", "dstartdate", "drug_substance"))


## Pre 2020-10-14 ----
analysis_pre_2020 <- t2d_1stinstance %>%
  ## Exclusions based on time - initiations after 14th October 2020
  filter(dstartdate <= as.Date("2020-10-14")) %>%
  ## Exclusions based on type of drugsubstance:
  ### SU - non-gliclazide
  ### TZD - rosiglitazone
  ### GLP1-RA - lixisenatide, exenatide slow-release, semaglutide
  ### MFN
  filter(
    drug_class != "MFN" & drug_class != "INS" &
      drug_class != "Glinide"& drug_class != "Acarbose"
  ) %>%
  filter(
    drug_substance != "Glimepiride" & drug_substance != "Lixisenatide" &
      drug_substance != "Glipizide" & drug_substance != "Ertugliflozin" &
      drug_substance != "Glibenclamide" & drug_substance != "Tolbutamide" &
      drug_substance != "Low-dose semaglutide" & drug_substance != "Oral semaglutide" &
      drug_substance != "Semaglutide, dose unclear"
  ) %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <80
  filter(is.na(INS) | INS != 1) %>%
  filter(drugline_all != 1) %>%
  filter(!(preckdstage %in% c("stage 5"))) %>%
  mutate(dstartdate_age = (datediff(dstartdate, dob))/365.25) %>% # not created?
  filter(dstartdate_age > 18 & dstartdate_age < 80) %>%
  ## Exclusions
  ### Multiple treatments on same day
  ### Start within 61 days since last therapy
  ### missing HbA1c
  ### baseline HbA1c<53
  ### baseline HbA1c>120
  filter(multi_drug_start_class == 0) %>%
  filter(timeprevcombo_class >= 61) %>%
  filter(!is.na(prehba1c)) %>%
  filter(prehba1c >= 53) %>%
  filter(prehba1c <= 110) %>%
  ## Exclusions
  ### missing baseline (BMI, eGFR, ALT, HDL, total cholesterol)
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  filter(!is.na(prebmi)) %>%
  filter(!is.na(preegfr)) %>%
  filter(!is.na(prealt)) %>%
  filter(!is.na(prehdl)) %>%
  filter(!is.na(pretotalcholesterol)) %>%
  # ## Exclusions
  # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  filter(!is.na(posthba1cfinal)) %>%
  analysis$cached("analysis_pre_2020", indexes=c("patid", "dstartdate", "drug_substance"))




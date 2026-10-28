
# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("mm")

# set up dataset
t2d_1stinstance <- t2d_1stinstance %>%
  analysis$cached("20260205_t2d_1stinstance")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

# Flow diagram: inclusion/exclusions ----
analysis = cprd$analysis("pedro_mm")
## Datasets

#:-- All initiations with T2D treatment selection criteria
analysis_ethnicity <- t2d_1stinstance %>%
  ## Exclusions based on type of drugsubstance: 
  ### SU - non-gliclazide
  ### TZD - rosiglitazone
  ### GLP1-RA - lixisenatide, exenatide slow-release
  ### MFN
  filter(
    drug_class != "MFN" & drug_class != "INS" &
      drug_class != "Glinide" & drug_class != "Acarbose" & drug_class != "GIPGLP1"
  ) %>%
  filter(
    drug_substance != "Glimepiride" & drug_substance != "Lixisenatide" &
      drug_substance != "Glipizide" & drug_substance != "Ertugliflozin" &
      drug_substance != "Glibenclamide" & drug_substance != "Tolbutamide"
  ) %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <18 or >80
  filter(is.na(INS) | INS != 1) %>%
  filter(drugline_all != 1) %>%
  filter(is.na(preckdstage) | preckdstage != "stage_5") %>%
  filter(dstartdate_age > 18 & dstartdate_age < 80) %>%
  ## Exclusions
  ### Multiple treatments on same day
  ### Start within 61 days since last therapy
  ### missing HbA1c
  ### baseline HbA1c<53
  ### baseline HbA1c>110
  filter(multi_drug_start_class == 0) %>%
  filter(timeprevcombo_class >= 61) %>%
  filter(!is.na(prehba1c)) %>%
  filter(prehba1c >= 53) %>%
  filter(prehba1c <= 110) %>%
  ## Exclusions
  ### missing baseline (duration of diabetes)
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  ## Outcome derivation (no exclusion applied here)
  ### Initiations with a missing outcome HbA1c are kept in this cohort so that
  ### missingness can be described, and are only dropped in the analysis scripts
  ### once the cohort has been split by ethnicity.
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  analysis$cached("ethnicity", indexes=c("patid", "dstartdate", "drug_substance"))





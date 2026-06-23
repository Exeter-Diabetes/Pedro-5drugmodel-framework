
# load libraries
library(tidyverse)
library(aurum)
library(rms)
library(qpdf)
library(tableone)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("mm")

# set up dataset
t2d_1stinstance <- t2d_1stinstance %>%
  analysis$cached("20250327_t2d_1stinstance")


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
# analysis_semaglutide_raw <- analysis_semaglutide_raw %>%
#   analysis$cached("analysis_injectable_semaglutide") %>%
#   collect() %>%
#   mutate(patid=as.character(patid)) %>%
#   mutate_if(is.integer64, as.integer)


# Pre-processing datasets ########################################

## Semaglutide ----
analysis_semaglutide <- t2d_1stinstance %>%
  ## Inclusion based on type of drugsubstance:
  ### GLP1-RA - semaglutide
  filter(
    drug_substance == "Low-dose semaglutide" | drug_substance == "High-dose semaglutide" |
      drug_substance == "Semaglutide, dose unclear"
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
  # filter(!is.na(prebmi)) %>%
  # filter(!is.na(preegfr)) %>%
  # filter(!is.na(prealt)) %>%
  # filter(!is.na(prehdl)) %>%
  # filter(!is.na(pretotalcholesterol)) %>%
  collect() %>%
  # # ## Exclusions
  # # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), ifelse(is.na(posthba1c6m), NA, posthba1c6m), posthba1c12m)) %>%  # generating posthba1c
  # # generating HbA1c month
  # mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  # mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  # mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  # filter(!is.na(posthba1cfinal))  %>%
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
    "drugline", "ncurrtx", "posthba1cfinal",
    "stopdrug_12m_6mFU"
  )
  )) %>%
  rename("drugclass" = "drug_class") %>%
  as.data.frame()


# Table of characteristics ----

## Missing vars
analysis_semaglutide <- analysis_semaglutide %>% 
  mutate(
    prebmi_missing = ifelse(is.na(prebmi), "Yes", "No"),
    preegfr_missing = ifelse(is.na(preegfr), "Yes", "No"),
    pretotalcholesterol_missing = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
    prehdl_missing = ifelse(is.na(prehdl), "Yes", "No"),
    prealt_missing = ifelse(is.na(prealt), "Yes", "No"),
    posthba1cfinal_missing = ifelse(is.na(posthba1cfinal), "Yes", "No"),
    stopdrug_12m_6mFU = ifelse(is.na(stopdrug_12m_6mFU), 1, stopdrug_12m_6mFU),
    stopdrug_12m_6mFU = ifelse(stopdrug_12m_6mFU == 1, "Yes", "No"),
    outcome_9mplus_and_above_58 = ifelse(stopdrug_12m_6mFU == "No" & !is.na(posthba1cfinal) & posthba1cfinal > 58, "Yes", "No"),
    resphba1c = posthba1cfinal - prehba1c
  )


vars <- c(
  "agetx", "sex", "t2dmduration", "ethnicity", 
  "drug_substance",
  "imd5", "smoke",
  "prebmi", "prebmi_missing", "prehba1c", "preegfr", "preegfr_missing", 
  "pretotalcholesterol", "pretotalcholesterol_missing", "prehdl", "prehdl_missing", 
  "prealt", "prealt_missing",
  "drugline", "ncurrtx", "hba1cmonth",
  "posthba1cfinal",
  "posthba1cfinal_missing",
  "resphba1c",
  "stopdrug_12m_6mFU",
  "outcome_9mplus_and_above_58"
)

cat_vars <- c(
  "sex", "ethnicity", "drug_substance", "imd5", "smoke",
  "drugline", "ncurrtx", "prebmi_missing", "preegfr_missing", "pretotalcholesterol_missing", "prehdl_missing", "prealt_missing",
  "stopdrug_12m_6mFU", "outcome_9mplus_and_above_58"
)


table_characteristics <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  # strata = c("drugclass"),
  data = analysis_semaglutide,
  test = FALSE
)

table_characteristics_print <- print(table_characteristics, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)


write.csv(table_characteristics_print, file = "Outputs/CPRD/10.inject_semaglutide_outcomes.csv")





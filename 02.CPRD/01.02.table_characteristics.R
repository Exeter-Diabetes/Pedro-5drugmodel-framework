# load libraries
library(tidyverse)
library(aurum)
library(tableone)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

# set up functions
is.integer64 <- function(x){
  class(x)=="integer64"
}
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


#:------------------------------------------------------
# Pre-processing datasets
#:-- Pre 2020-10-14
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
  select(-row)

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
  select(-row)



#:-----------------------------------------------------------------------------
# Table of characteristics for each dataset
#:-----------------------------------------------------------------------------

vars <- c(
  "agetx", "sex", "t2dmduration", "ethnicity", 
  "drug_substance",
  "imd5", "smoke",
  "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
  "drugline", "ncurrtx", "hba1cmonth",
  "posthba1cfinal"
)

cat_vars <- c(
  "sex", "ethnicity", "drug_substance", "imd5", "smoke",
  "drugline", "ncurrtx"
)


#:-- Pre 2020-10-14
table_characteristics_pre2020 <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  strata = c("drug_class"),
  data = analysis_pre_2020,
  test = FALSE
)

table_characteristics_pre2020_print <- print(table_characteristics_pre2020, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)


write.csv(table_characteristics_pre2020_print, file = "Outputs/CPRD/01.02.characteristics_pre2020.csv")


#:-- Post 2020-10-14
table_characteristics_post2020 <- CreateTableOne(
  vars = vars,
  factorVars = cat_vars,
  includeNA = TRUE,
  strata = c("drug_class"),
  data = analysis_post_2020,
  test = FALSE
)

table_characteristics_post2020_print <- print(table_characteristics_post2020, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)


write.csv(table_characteristics_post2020_print, file = "Outputs/CPRD/01.02.characteristics_post2020.csv")




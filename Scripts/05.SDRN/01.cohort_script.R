## SDRN version of Scripts/01.cohort_scripts.R: builds the SDRN cohort used by the
##  validation scripts in this folder (02.validation_v1.R and 03.validation_v2.R).
##  The inclusions and exclusions are those of the CPRD cohort, applied to the SDRN
##  treatment response dataset (t2d_1stinstance, made by SDRN-Cohort-scripts,
##  03-Treatment-response-(MASTERMIND)/11_final_merge.R), all initiations of any date.
##  Where SDRN codes a variable differently from CPRD (ethnicity), it is recoded to
##  the CPRD coding, so that the validation scripts read it as the CPRD test runs do.
##
## The cohort is saved in Scripts/05.SDRN/Data/ (not in git). Patient-level data
##  never leaves the SDRN server: this script only prints counts (the flow diagram),
##  with counts below 5 shown as "<5".
##
## Run on the SDRN server before the validation scripts.


# Initial set-up ###########################################

# load libraries
library(tidyverse)

# the repository root on the SDRN server (every path below is relative to it)
setwd("/home/pcardoso/workspace/Pedro-5drugmodel-framework")

# SDRN treatment response dataset: change the date to the newest version
sdrn_dataset <- "mm_20250506_t2d_1stinstance"

# where the cohort is saved (not in git)
cohort_file <- "Scripts/05.SDRN/Data/sdrn_cohort.RData"

# load dataset
load(paste0("/home/pcardoso/workspace/SDRN-Cohort-scripts/Final_Datasets/", sdrn_dataset, ".RData")) # t2d_1stinstance

## functions ----
# counts below 5 shown as "<5" (PHS disclosure rule)
n5 <- function(n) ifelse(n > 0 & n < 5, "<5", format(n, big.mark = ",", trim = TRUE))

# number of initiations after each step of the flow diagram
flow <- tibble(step = character(), n = integer())
add_step <- function(data, step) {
  n <- nrow(data) # first, so that the earlier steps of the pipe are counted first
  flow <<- bind_rows(flow, tibble(step = step, n = n))
  data
}


# SDRN codings changed to the CPRD ones ----
## Ethnicity: SDRN has the 5 categories as names ("Not stated/Unknown" when not
## recorded), CPRD codes them 0-5 (0 White, 1 South Asian, 2 Black, 3 Other,
## 4 Mixed, 5 Missing). A dataset already coded 0-5 is kept as it is.
ethnicity_codes <- c("White" = 0, "South Asian" = 1, "Black" = 2, "Other" = 3,
                     "Mixed" = 4, "Not stated/Unknown" = 5)

unknown_ethnicity <- setdiff(unique(as.character(t2d_1stinstance$ethnicity_5cat)),
                             c(names(ethnicity_codes), as.character(0:5), NA))
if (length(unknown_ethnicity) > 0) {
  stop("ethnicity_5cat values with no CPRD code: ", paste(unknown_ethnicity, collapse = ", "))
}

t2d_1stinstance <- t2d_1stinstance %>%
  mutate(
    ethnicity_5cat = as.character(ethnicity_5cat),
    ethnicity_5cat = ifelse(ethnicity_5cat %in% names(ethnicity_codes),
                            unname(ethnicity_codes[ethnicity_5cat]),
                            suppressWarnings(as.numeric(ethnicity_5cat)))
  )


# Flow diagram: inclusion/exclusions ----

#:-- All initiations with T2D treatment selection criteria
analysis_cohort <- t2d_1stinstance %>%
  add_step("initiations in t2d_1stinstance") %>%
  ## Exclusions based on type of drugsubstance (SDRN also has some older substances,
  ## excluded with the others of their drug class):
  ### SU - non-gliclazide
  ### TZD - troglitazone (rosiglitazone is kept, as in Scripts/01.cohort_scripts.R)
  ### SGLT2 - ertugliflozin
  ### GLP1-RA - lixisenatide, albiglutide (exenatide, including prolonged-release, is
  ###   kept: the analysis scripts group it with liraglutide)
  ### MFN
  filter(
    drug_class != "MFN" & drug_class != "INS" &
      drug_class != "Glinide" & drug_class != "Acarbose" & drug_class != "GIPGLP1"
  ) %>%
  filter(
    !(drug_substance %in% c(
      "Glimepiride", "Glipizide", "Glibenclamide", "Tolbutamide", "Tolazamide",
      "Gliquidone", "Chlorpropamide", "Acetohexamide", "Glymidine",
      "Troglitazone", "Ertugliflozin", "Lixisenatide", "Albiglutide"
    ))
  ) %>%
  add_step("drug class / substance analysed") %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <18 or >80
  filter(is.na(INS) | INS != 1) %>%
  add_step("not on insulin") %>%
  filter(drugline_all != 1) %>%
  add_step("not first-line") %>%
  filter(is.na(preckdstage) | preckdstage != "stage_5") %>%
  add_step("no CKD stage 5") %>%
  filter(dstartdate_age > 18 & dstartdate_age < 80) %>%
  add_step("age 18-80") %>%
  ## Exclusions
  ### Multiple treatments on same day
  ### Start within 61 days since last therapy
  ### missing HbA1c
  ### baseline HbA1c<53
  ### baseline HbA1c>110
  filter(multi_drug_start_class == 0) %>%
  add_step("one drug class started that day") %>%
  filter(timeprevcombo_class >= 61) %>%
  add_step(">=61 days since last treatment change") %>%
  filter(!is.na(prehba1c)) %>%
  filter(prehba1c >= 53) %>%
  filter(prehba1c <= 110) %>%
  add_step("baseline HbA1c 53-110") %>%
  ## Exclusions
  ### missing baseline (duration of diabetes)
  mutate(t2dmduration = as.numeric(difftime(dstartdate, dm_diag_date_all, units = "days"))/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  add_step("diabetes duration known") %>%
  ## Outcome derivation (no exclusion applied here)
  ### Initiations with a missing outcome HbA1c are kept in this cohort so that
  ### missingness can be described, and are only dropped in the validation scripts
  ### (after imputation).
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = as.numeric(difftime(posthba1c12mdate, dstartdate, units = "days")) / 30) %>%
  mutate(hba1cmonth_6 = as.numeric(difftime(posthba1c6mdate, dstartdate, units = "days")) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12))


# Counts (rows of t2d_1stinstance) ----
for (i in seq_len(nrow(flow))) {
  cat(paste0("S", i, ": ", flow$step[i], " = ", n5(flow$n[i]), "\n"))
}

by_class <- analysis_cohort %>% count(drug_class)
by_glp1 <- analysis_cohort %>% filter(drug_class == "GLP1") %>% count(drug_substance)
cat("C1: by drug class:", paste(by_class$drug_class, n5(by_class$n), collapse = " | "), "\n")
cat("C2: GLP1 by substance:", paste(by_glp1$drug_substance, n5(by_glp1$n), collapse = " | "), "\n")
cat("C3: outcome HbA1c missing (dropped in the validation scripts) =",
    n5(sum(is.na(analysis_cohort$posthba1cfinal))), "\n")
cat("C4: people =", n5(n_distinct(analysis_cohort$serialno)), "| initiations from",
    format(min(analysis_cohort$dstartdate)), "to", format(max(analysis_cohort$dstartdate)), "\n")
cat("C5: ethnicity (0 White, 1 South Asian, 2 Black, 3 Other, 4 Mixed, 5/NA missing):",
    paste(names(table(analysis_cohort$ethnicity_5cat, useNA = "ifany")),
          n5(as.numeric(table(analysis_cohort$ethnicity_5cat, useNA = "ifany"))), collapse = " | "), "\n")


# Save the cohort ----
sdrn_cohort <- analysis_cohort
dir.create(dirname(cohort_file), showWarnings = FALSE, recursive = TRUE)
save(sdrn_cohort, file = cohort_file)

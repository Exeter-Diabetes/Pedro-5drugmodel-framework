## SDRN run of example_validation_script_v1.R: all SDRN initiations (any date) of
##  the published model's drug classes, with GLP1 as one class and semaglutide left
##  out (it is not part of the published model's GLP1 class), as in the CPRD test
##  run (Scripts/00.V1_Validation). The script is v1 unchanged apart from lines
##  marked "SDRN run": the working directory is the repository on the SDRN server,
##  the SDRN cohort (made by Scripts/05.SDRN/01.cohort_script.R) is loaded in
##  Part 2a, predictions come from the model file in Model_files/, each patient's
##  own deprivation (imd5, the SIMD quintile) is used instead of the middle quintile
##  (in the predictions, and adjusted and matched on), counts below 5 in the
##  characteristics table are shown as "<5", and outputs go to
##  Scripts/05.SDRN/Outputs/V1. When v1 changes, copy it again and redo these edits.
##
## This script shows how to validate the 5-drug model in a new cohort and how to
##  visualise the results.
##
## Model version: v1, the published model (Dennis et al. Lancet 2025), for these
##  drug classes: SGLT2i, GLP-1RA, DPP4i, SU, TZD
##  The validation can be done for a subset of these drug classes (drugs_selected,
##  Part 1). To validate the updated model (v2), use example_validation_script_v2.R
##  instead.
##
## Model predictions can be made using the online calculator or the model object
##  (both available on request; use_calculator, Part 1).
##
## The script has five parts:
## 1) Settings
## 2) Prepare your data
##    a) Format the variables
##    b) Characteristics table
##    c) Impute missing data and choose the drug classes analysed
##    d) File for the online calculator (calculator only)
## 3) Model predictions for each drug class, in four versions:
##    a) original predictions
##    b) minimal recalibration-in-the-large (a drug-specific intercept, only where
##       it is beyond 1.5 mmol/mol)
##    c) recalibration-in-the-large (a drug-specific intercept)
##    d) recalibration (a drug-specific intercept and slope)
## 4) Figure code (needs to be run; no need to read it)
## 5) Validation, done for each version of the predictions in turn, with every
##    figure saved as it is made. The files sort in the same order:
##      1_original_..., 2_minimal_recalibration_in_the_large_...,
##      3_recalibration_in_the_large_..., 4_recalibration_...,
##    and within each, by figure: 01_pairwise, 02_overall, 03_overall_by_drug,
##    04_loveplot, 05_pred_obs, 06_rank1, 07_tolerance3.
##    Then the summaries (5_summary_...).


#############################################################################

# load required libraries
library(tidyverse)
library(tableone)
# library(remotes)
# remotes::install_github("PM-Cardoso/cateval")
library(cateval) # Custom R package developed by Pedro for calibration evaluation


#############################################################################

# Part 1: Settings ------------------------------------------------------------

# How the model predictions are made (Part 3):
#   TRUE:  with the online calculator. Run the script twice: the first run
#          writes the file to upload (calculator_upload_file) and stops. Upload it in
#          the calculator's Validation Framework tab, save the file it returns as
#          calculator_results_file, and run the script again.
#   FALSE: from the model object, which needs to be in the working directory.
use_calculator <- FALSE # SDRN run: model files (Model_files/), patient-level data stays on the server

# SDRN run: the repository on the SDRN server (every path below is relative to it)
setwd("/home/pcardoso/workspace/Pedro-5drugmodel-framework")

# where the outputs are saved
output_path <- "Scripts/05.SDRN/Outputs/V1" # SDRN run
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# files for the online calculator (only used when use_calculator is TRUE)
calculator_upload_file  <- paste0(output_path, "/calculator_upload.csv") # written in Part 2d
calculator_results_file <- "calculator_predictions.csv"                   # returned by the calculator

# The drug classes the model predicts (do not change these)
model_drugs <- c("SGLT2", "GLP1", "DPP4", "SU", "TZD")

# The drug classes to analyse. You do not need initiations of every drug class: the
# drug classes analysed are those of drugs_selected with initiations in your data
# (drugs, set in Part 2c). To leave out a drug class, e.g. one with too few
# initiations, remove it here, e.g.
#   drugs_selected <- c("SGLT2", "GLP1", "DPP4", "SU")
# Its initiations are dropped after imputation (Part 2c), so they still help impute
# missing values. The model still predicts it (original predictions, not
# recalibrated), so it can be the best drug class in the "all" figures (Part 5).
drugs_selected <- model_drugs

if (!all(drugs_selected %in% model_drugs)) {
  stop("drugs_selected has drug classes the model does not predict: ",
       paste(setdiff(drugs_selected, model_drugs), collapse = ", "), call. = FALSE)
}

# drug class names as shown in the plots
drug_labels <- c("SGLT2" = "SGLT2i", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i",
                 "SU" = "SU", "TZD" = "TZD")

# drug class colours (Okabe-Ito), as in the online calculator
drug_colours <- c("SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "DPP4" = "#0072B2",
                  "SU" = "#CC79A7", "TZD" = "#D55E00")


#############################################################################

# Part 2: Prepare your data ------------------------------------------------------

# Part 2a: Format the variables ----
# Note: This portion of code can be adapted based on your data source and structure.
# The cohort (who is included, and how each variable is derived) is described in
# the cohort description document (available on request). In short: new users of
# one of the drug classes below, not on insulin, initiating only one
# glucose-lowering drug class that day.

# REQUIRED VARIABLES
# Your dataset (analysis_dataset) must include the following variables for validation:
#
# OUTCOME VARIABLES:
#   - posthba1cfinal: outcome HbA1c (numeric, mmol/mol): the HbA1c closest to 12
#     months after initiation (between 9 and 15 months) or, if there is none, the
#     HbA1c closest to 6 months (between 3 and 9 months)
#   - hba1cmonth: months from initiation to that outcome HbA1c (numeric, 3-15 months)
#
# TREATMENT VARIABLES:
#   - drugclass: Drug class being initiated (categorical). Code the drug classes in
#     your data with these names (you do not need all of them):
#       "SGLT2", "GLP1", "DPP4", "SU", or "TZD"
#     Initiations with any other value are dropped.
#
# BASELINE CHARACTERISTICS:
#   - agetx: age at initiation (numeric, 18-79 years)
#   - sex: "Male" or "Female"
#   - t2dmduration: duration of type 2 diabetes, from diagnosis to initiation
#     (numeric, 0-71 years)
#   - prehba1c: baseline HbA1c, closest to initiation within 6 months before to 7
#     days after (numeric, 53-110 mmol/mol)
#   The next five are the value closest to initiation within 2 years before to 7
#   days after:
#   - prebmi: BMI (numeric, 15-100 kg/m2)
#   - preegfr: eGFR (numeric, 45-300 mL/min per 1.73 m2)
#   - prealt: alanine aminotransferase (numeric, 0-200 U/L)
#   - pretotalcholesterol: total cholesterol (numeric, 0.5-20 mmol/L)
#   - prehdl: HDL cholesterol (numeric, 0.2-10 mmol/L)
#   - ethnicity: "White", "South Asian", "Black", "Mixed", "Other" or "Missing"
#   - smoke: smoking status, "Non-smoker", "Active smoker", "Ex-smoker" or
#     "Not recorded"
#   - ncurrtx: number of glucose-lowering drug classes being taken, including the
#     one initiated ("1", "2", "3" or "4+")
#   - drugline: number of glucose-lowering drug classes ever taken, including the
#     one initiated ("2", "3", "4" or "5+")

# The ranges above are the ones the online calculator accepts. Initiations outside
# them are dropped after imputation (Part 2c), whichever way the predictions are made.
valid_ranges <- list(
  agetx = c(18, 79), t2dmduration = c(0, 71), prehba1c = c(53, 110),
  prebmi = c(15, 100), preegfr = c(45, 300), prealt = c(0, 200),
  pretotalcholesterol = c(0.5, 20), prehdl = c(0.2, 10), hba1cmonth = c(3, 15)
)

## SDRN run: load the SDRN cohort (start) ---------------------------------------
## Made by Scripts/05.SDRN/01.cohort_script.R (the inclusions and exclusions of the
## CPRD cohort). All initiations, any date, formatted as in the CPRD test run
## (Scripts/00.V1_Validation): GLP1 as one class, semaglutide left out.

# Semaglutide substances (GLP1 drug_class in the cohort), left out: semaglutide is
# not part of the published model's GLP1 class. SDRN writes "Semaglutide dose
# unclear" without the comma, so both spellings are listed.
sema_substances <- c("Low-dose semaglutide", "High-dose semaglutide", "Oral semaglutide",
                     "Semaglutide, dose unclear", "Semaglutide dose unclear")

# load the cohort (sdrn_cohort) and select the variables needed
load("Scripts/05.SDRN/Data/sdrn_cohort.RData")
analysis_raw <- sdrn_cohort %>%
  select(
    serialno, dstartdate,
    gender, dstartdate_age,
    t2dmduration, ethnicity_5cat,
    drug_substance, drug_class,
    imd_decile, smoking_cat,
    prebmi, prehba1c, preegfr, pretotalcholesterol, prehdl, prealt,
    drugline_all, MFN, SGLT2, GLP1, DPP4, TZD, SU, hba1cmonth,
    posthba1cfinal
  )
rm(sdrn_cohort)

analysis_dataset <- analysis_raw %>%
  mutate(
    pated = paste(serialno, dstartdate, drug_class, sep = "."),
    agetx = dstartdate_age,

    # Sex
    sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),

    # Ethnicity (coded 0-5 in 01.cohort_script.R, as in CPRD)
    ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
    ethnicity = factor(ethnicity, levels = c(0:5),
                       labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),

    # Smoking
    smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),

    # Deprivation (SIMD quintile; imd_decile is the SIMD decile flipped, so 1 is the
    # least deprived, as in CPRD), used in the model predictions (Part 3)
    imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
    imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5),
                  labels = c("1 (least)", "2", "3", "4", "5 (most)")),

    # Number of current treatments
    ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
    ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
    ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),

    # Drug line
    drugline = ifelse(drugline_all > 5, 5, drugline_all),
    drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),

    # Drug class (GLP1 as one class)
    drugclass = drug_class
  ) %>%
  group_by(pated) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  select(-row) %>%
  filter(!(drug_substance %in% sema_substances)) %>%
  select(all_of(c(
    "dstartdate", "agetx", "sex", "t2dmduration", "ethnicity",
    "drugclass", "imd5", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  ))) %>%
  as.data.frame()

cat("cohort: initiations from",
    format(min(analysis_dataset$dstartdate, na.rm = TRUE)), "to",
    format(max(analysis_dataset$dstartdate, na.rm = TRUE)), "\n")
table(analysis_dataset$drugclass, useNA = "ifany")

# v1 imputes with every column of analysis_dataset (Part 2c), so drop the date
analysis_dataset <- analysis_dataset %>% select(-dstartdate)

## SDRN run: load the SDRN cohort (end) -----------------------------------------


# Initiations whose drugclass is not one of the drug class names above are not
# analysed: this prints how many there are of each value, to catch a drug class
# coded with another name (e.g. "SGLT2i" or "GLP-1RA") before it is dropped below.
analysis_dataset %>%
  filter(!(as.character(drugclass) %in% model_drugs)) %>%
  count(drugclass, name = "n_dropped") %>%
  print()


# Categorical variables must be properly formatted as factors with the specified levels.
# Example formatting code:
analysis_dataset <- analysis_dataset %>%
  mutate(
    # Sex
    sex = factor(sex,
                 levels = c("Male", "Female")),

    # Ethnicity
    ethnicity = factor(ethnicity,
                       levels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),

    # Smoking
    smoke = factor(smoke,
                   levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),

    # Number of current treatments
    ncurrtx = factor(ncurrtx,
                     levels = c("1", "2", "3", "4+")),

    # Drug line
    drugline = factor(drugline,
                      levels = c("2", "3", "4", "5+")),

    # Drug class
    drugclass = factor(drugclass,
                       levels = model_drugs)
  ) %>%
  # keep only the drug classes of this model version
  filter(!is.na(drugclass))

# the validation compares drug classes, so it needs initiations of at least two
if (n_distinct(analysis_dataset$drugclass) < 2) {
  stop("At least two drug classes with initiations are needed for the validation; found: ",
       paste(unique(analysis_dataset$drugclass), collapse = ", "), ".", call. = FALSE)
}


# Part 2b: Characteristics table ----

# This table describes the validation cohort at treatment initiation, stratified by
# initiated drug class.
# It includes baseline clinical characteristics, missingness indicators, the timing
# of the outcome, the outcome HbA1c (at 12 months, or at 6 months when the
# 12-month value is missing; see Part 2a) and the observed response
# (posthba1cfinal - prehba1c).


# Add observed response and missingness indicators
analysis_dataset_drugclass <- analysis_dataset %>%
  mutate(
    # only the drug classes in your data
    drugclass = droplevels(drugclass),
    observed_response = posthba1cfinal - prehba1c,
    # HbA1c and response also expressed as DCCT % for the characteristics table
    prehba1c_pct = (prehba1c / 10.929) + 2.15,
    posthba1cfinal_pct = (posthba1cfinal / 10.929) + 2.15,
    observed_response_pct = observed_response / 10.929,
    missing_agetx = as.numeric(is.na(agetx)),
    missing_t2dmduration = as.numeric(is.na(t2dmduration)),
    missing_prehba1c = as.numeric(is.na(prehba1c)),
    missing_prebmi = as.numeric(is.na(prebmi)),
    missing_preegfr = as.numeric(is.na(preegfr)),
    missing_pretotalcholesterol = as.numeric(is.na(pretotalcholesterol)),
    missing_prehdl = as.numeric(is.na(prehdl)),
    missing_prealt = as.numeric(is.na(prealt)),
    missing_posthba1cfinal = as.numeric(is.na(posthba1cfinal))
  )

## Table of characteristics stratified by drug class
drugclass_vars <- c(
  "agetx", "missing_agetx",
  "t2dmduration", "missing_t2dmduration",
  "sex", "ethnicity", "smoke",
  "drugline", "ncurrtx",

  "prebmi", "missing_prebmi",

  "prehba1c", "prehba1c_pct", "missing_prehba1c",

  "preegfr", "missing_preegfr",
  "pretotalcholesterol", "missing_pretotalcholesterol",
  "prehdl", "missing_prehdl",
  "prealt", "missing_prealt",

  "posthba1cfinal", "posthba1cfinal_pct", "missing_posthba1cfinal",
  "observed_response", "observed_response_pct",
  "hba1cmonth"
)

factor_vars <- c(
  "sex", "ethnicity", "smoke", "drugline", "ncurrtx",
  "missing_agetx", "missing_t2dmduration", "missing_prehba1c",
  "missing_prebmi", "missing_preegfr", "missing_pretotalcholesterol",
  "missing_prehdl", "missing_prealt", "missing_posthba1cfinal"
)

characteristics_table <- CreateTableOne(
  vars = drugclass_vars,
  strata = "drugclass",
  data = analysis_dataset_drugclass,
  factorVars = factor_vars,
  includeNA = TRUE,
  test = FALSE,
  addOverall = TRUE # adds a column for the whole cohort
)

## Export table
characteristics_table_df <- print(
  characteristics_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  includeNA = TRUE
)

## SDRN run: counts below 5 shown as "<5" (PHS disclosure rule) (start) --------
## The counts in the table are the column sizes (row n) and the "n (%)" cells of the
## categorical variables. Each count from 1 to 4 is replaced by "<5" (with its
## percentage, which would give it away), and every cell of a column with fewer
## than 5 initiations is replaced by "<5". The other summaries (mean (SD)) are kept.
suppress_small_counts <- function(tab) {
  # the variable of each row (the rows of a variable's other levels have no name)
  row_var <- rownames(tab)
  for (i in seq_along(row_var)[-1]) {
    if (row_var[i] == "") row_var[i] <- row_var[i - 1]
  }
  count_rows <- row_var == "n" | grepl("\\(%\\)$", row_var)
  leading_count <- function(x) suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", x)))
  for (j in which(colnames(tab) != "level")) {
    if (isTRUE(leading_count(tab["n", j]) < 5)) {
      tab[, j] <- "<5"
    } else {
      counts <- leading_count(tab[count_rows, j])
      tab[count_rows, j][!is.na(counts) & counts >= 1 & counts < 5] <- "<5"
    }
  }
  tab
}
characteristics_table_df <- suppress_small_counts(characteristics_table_df)
## SDRN run: counts below 5 (end) -----------------------------------------------

write.csv(
  characteristics_table_df,
  paste0(output_path, "/characteristics_by_drugclass.csv"),
  row.names = TRUE
)


# Part 2c: Impute missing data and choose the drug classes analysed ----

# Missing baseline values are filled in by multiple imputation (mice):
#   1. the incomplete dataset is imputed mice_m = 10 times, giving 10 complete
#      datasets
#   2. the model predictions are made in each of them and averaged (Part 3)
#   3. the analyses (Part 5) also use the clinical features themselves: the
#      pairwise calibration adjusts for them and the overall benefit calibration
#      matches on them. For these, each missing value is replaced by its average
#      over the 10 imputed datasets, giving one complete dataset (analysis_dataset).
# The outcome (posthba1cfinal) is not imputed.

# Example imputation code (using mice package):

# load require libraries
library(mice)

# Define variables to exclude from imputation
# The outcome (posthba1cfinal) and its timing (hba1cmonth) are excluded so they are
# neither imputed nor used to impute other variables
# Add any other variables in your dataset that should not be used (e.g., IDs)
vars_to_exclude <- c("hba1cmonth", "posthba1cfinal")

# Filter to only those that actually exist in your dataset
vars_to_exclude <- vars_to_exclude[vars_to_exclude %in% names(analysis_dataset)]

# Create dataset for imputation (excluding specified variables)
vars_to_impute <- setdiff(names(analysis_dataset), vars_to_exclude)
data_for_imputation <- analysis_dataset[, vars_to_impute]

# Identify which values are missing in the dataset to be imputed
missing_indicator <- is.na(data_for_imputation)

# Perform multiple imputation (m=10 imputed datasets)
mice_m = 10; mice_maxit = 20
imputed_data <- mice(data_for_imputation, m = mice_m, method = 'pmm', maxit = mice_maxit, seed = 500, printFlag = TRUE)

# Extract all imputed datasets and calculate mean of imputed values
# Start with a copy of the imputed data
analysis_dataset_complete <- data_for_imputation

# For each variable with missing data, replace with mean of imputed values
for (var in names(data_for_imputation)) {
  if (any(missing_indicator[, var])) {
    # Extract imputed values for this variable across all m imputations
    imputed_values <- sapply(1:mice_m, function(i) {
      complete(imputed_data, action = i)[[var]]
    })
    # Calculate row-wise mean of imputed values
    mean_imputed <- rowMeans(imputed_values)
    # Replace missing values in original data with mean imputed values
    analysis_dataset_complete[missing_indicator[, var], var] <- mean_imputed[missing_indicator[, var]]
  }
}

# Add back the excluded variables
for (var in vars_to_exclude) {
  analysis_dataset_complete[[var]] <- analysis_dataset[[var]]
}

analysis_dataset <- analysis_dataset_complete

# Each imputed dataset, with the excluded variables added back, for the model
# predictions (Part 3). These have the same rows, in the same order, as analysis_dataset.
imputed_list <- lapply(1:mice_m, function(i) {
  imputed <- complete(imputed_data, action = i)
  for (var in vars_to_exclude) {
    imputed[[var]] <- analysis_dataset[[var]]
  }
  imputed
})

# Initiations with a missing outcome HbA1c are dropped only now, so that
# missingness is described in the characteristics table and imputation uses
# the full cohort
keep <- !is.na(analysis_dataset$posthba1cfinal)
analysis_dataset <- analysis_dataset[keep, ]
imputed_list <- lapply(imputed_list, function(imputed) imputed[keep, ])

# Initiations outside the accepted ranges (valid_ranges, Part 2a) in any of the
# imputed datasets are dropped, so the same initiations are analysed whichever
# way the predictions are made
in_range <- Reduce(`&`, lapply(imputed_list, function(imputed) {
  Reduce(`&`, lapply(names(valid_ranges), function(var) {
    x <- imputed[[var]]
    !is.na(x) & x >= valid_ranges[[var]][1] & x <= valid_ranges[[var]][2]
  }))
}))
cat(sum(!in_range), "initiations outside the accepted ranges in at least one imputed dataset are dropped\n")
analysis_dataset <- analysis_dataset[in_range, ]
imputed_list <- lapply(imputed_list, function(imputed) imputed[in_range, ])

# Initiations of drug classes not in drugs_selected (Part 1) are dropped
selected <- analysis_dataset$drugclass %in% drugs_selected
cat(sum(!selected), "initiations of drug classes not in drugs_selected are dropped\n")
analysis_dataset <- analysis_dataset[selected, ]
imputed_list <- lapply(imputed_list, function(imputed) imputed[selected, ])

# The drug classes analysed: those of drugs_selected with initiations in your data.
drugs <- model_drugs[model_drugs %in% drugs_selected & model_drugs %in% analysis_dataset$drugclass]

analysis_dataset <- analysis_dataset %>%
  mutate(drugclass = factor(drugclass, levels = drugs))

drug_classes_analysed <- analysis_dataset %>%
  count(drugclass, name = "n_initiations")
print(drug_classes_analysed)
write.csv(drug_classes_analysed, paste0(output_path, "/drug_classes_analysed.csv"),
          row.names = FALSE)

if (length(drugs) < 2) {
  stop("At least two drug classes with initiations are needed for the validation; ",
       "found: ", paste(drugs, collapse = ", "), ".")
}


# Part 2d: File for the online calculator ------------------------------------

# Only when the predictions come from the online calculator (use_calculator, Part 1).
# The file holds every imputed dataset, one after the other (imputation says which),
# so it has mice_m rows per initiation, with numeric values rounded to 3 decimals. The calculator returns every column of this
# file in the same order, with the predictions added, so row_id links them back to
# these rows in Part 3.
if (use_calculator) {
  calculator_upload <- bind_rows(lapply(1:mice_m, function(i) {
    imputed_list[[i]] %>% mutate(imputation = i)
  })) %>%
    # numeric values rounded to 3 decimals, to keep the file small
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  calculator_upload$row_id <- seq_len(nrow(calculator_upload))

  calculator_columns <- c(
    "row_id", "imputation", "agetx", "sex", "ethnicity", "smoke", "t2dmduration",
    "prebmi", "prehba1c", "preegfr", "prealt", "pretotalcholesterol", "prehdl",
    "ncurrtx", "drugline", "drugclass", "hba1cmonth", "posthba1cfinal"
  )
  write.csv(calculator_upload[, calculator_columns], calculator_upload_file, row.names = FALSE)

  if (!file.exists(calculator_results_file)) {
    stop("Upload ", calculator_upload_file, " in the Validation Framework tab of the ",
         "online calculator, save the file it returns as ", calculator_results_file,
         " and run this script again.", call. = FALSE)
  }
}


#############################################################################

# Part 3: Model predictions ------------------------------------------------

# The model predictions are one column for each drug class in model_drugs (absolute
# values, not the change from baseline):
#   - pred.SGLT2: Predicted 12-month HbA1c on SGLT2 inhibitor (numeric, mmol/mol)
#   - pred.GLP1: Predicted 12-month HbA1c on GLP-1 receptor agonist (numeric, mmol/mol)
#   - pred.DPP4: Predicted 12-month HbA1c on DPP4 inhibitor (numeric, mmol/mol)
#   - pred.SU: Predicted 12-month HbA1c on Sulfonylurea (numeric, mmol/mol)
#   - pred.TZD: Predicted 12-month HbA1c on Thiazolidinedione (numeric, mmol/mol)

pred_cols <- paste0("pred.", model_drugs)

if (use_calculator) {

  ## Predictions from the online model calculator (Validation Framework tab)
  calculator_results <- read.csv(calculator_results_file, check.names = FALSE,
                                 stringsAsFactors = FALSE)

  # check the file is the one returned for this dataset: same rows, in the same order
  same_rows <- nrow(calculator_results) == nrow(calculator_upload) &&
    all(calculator_results$row_id == calculator_upload$row_id) &&
    all(as.character(calculator_results$drugclass) == as.character(calculator_upload$drugclass)) &&
    isTRUE(all.equal(calculator_results$prehba1c, calculator_upload$prehba1c, tolerance = 1e-6))
  if (!same_rows) {
    stop(calculator_results_file, " does not match this dataset. Upload the new ",
         calculator_upload_file, " to the calculator and replace ", calculator_results_file,
         " with the file it returns.", call. = FALSE)
  }

  missing_pred_cols <- setdiff(pred_cols, names(calculator_results))
  if (length(missing_pred_cols) > 0) {
    stop(calculator_results_file, " has no ", paste(missing_pred_cols, collapse = ", "),
         " column(s).", call. = FALSE)
  }

  # the initiation (row of analysis_dataset) of each row of the file
  initiation <- rep(seq_len(nrow(analysis_dataset)), times = mice_m)

  # average the predictions of each initiation across the imputed datasets
  for (col in pred_cols) {
    analysis_dataset[[col]] <- as.numeric(tapply(calculator_results[[col]], initiation, mean))
  }

  # rows the calculator could not predict have the reason in error_message. An
  # initiation is dropped if it could not be predicted in any of the imputed datasets.
  calculator_failed <- !is.na(calculator_results$error_message) &
    calculator_results$error_message != ""
  initiation_failed <- seq_len(nrow(analysis_dataset)) %in% initiation[calculator_failed]
  cat(sum(initiation_failed), "initiations could not be predicted by the calculator and are dropped\n")
  if (any(calculator_failed)) print(table(calculator_results$error_message[calculator_failed]))
  analysis_dataset <- analysis_dataset[!initiation_failed, ]

} else {

  ## Predictions from the model object, which needs to be in the working directory:
  ##   - fivedrugmodel_5knot_share_20230823.Rdata: the 5-drug model (m1.5.final)
  ## These give the same predictions as the online calculator.

  # load required libraries
  library(rms)

  # load model
  load("Model_files/fivedrugmodel_5knot_share_20230823.Rdata") # m1.5.final (SDRN run: path)

  # The model also uses deprivation (imd5), which is set to the middle quintile as
  # in the online calculator. hba1cmonth is each initiation's own.
  # (the model needs a data.frame, not a tibble)
  prediction_list <- lapply(imputed_list, function(imputed) {
    # SDRN run: each patient's own deprivation (imd5, the SIMD quintile, from the SDRN cohort)
    as.data.frame(imputed)
  })

  # Published model, for each of its 5 drug classes: predicted in each imputed
  # dataset and averaged across them
  for (drug in c("SGLT2", "GLP1", "DPP4", "SU", "TZD")) {
    analysis_dataset[[paste0("pred.", drug)]] <- rowMeans(sapply(prediction_list, function(prediction_data) {
      as.numeric(predict(m1.5.final, newdata = prediction_data %>% mutate(drugclass = drug)))
    }))
  }
}


## We then generate versions of the predictions that are recalibrated using intercept-only
## (only those beyond minimal_threshold, or all) and intercept + slope adjustments.

### Recalibration
# drugs (Part 2c) and model_drugs (Part 1)
outcome_var <- "posthba1cfinal"
pred_prefix <- "pred."

# Optional prognostic factors (leave empty if none). These are variables that
# affect the outcome irrespective of treatment, so they enter as a single overall
# (treatment-independent) effect: a shared main effect in the recalibration models
# below, plus covariate adjustment in the pairwise analysis and matching in the
# overall-benefit analysis (via matching_vars_prog). The original predictions are
# not adjusted for them.
prognostic_vars <- c()
prognostic_terms <- if (length(prognostic_vars) > 0) {
  paste0(" + ", paste(prognostic_vars, collapse = " + "))
} else {
  ""
}

# Minimal recalibration-in-the-large uses a drug class's intercept only when it is
# beyond minimal_threshold (mmol/mol): half the clinically relevant difference (3
# mmol/mol), so two drug classes whose intercepts are not used have a predicted
# difference off by at most 3 mmol/mol.
minimal_threshold <- 1.5

# summarising function: tidy a fitted lm into a coefficient table
summarize_lm <- function(fit, model_type) {
  smry <- summary(fit)
  coef_mat <- as.data.frame(smry$coefficients)
  coef_mat$term <- rownames(coef_mat)
  rownames(coef_mat) <- NULL
  names(coef_mat) <- c("estimate", "std_error", "t_value", "p_value", "term")

  ci_mat <- as.data.frame(confint(fit))
  ci_mat$term <- rownames(ci_mat)
  rownames(ci_mat) <- NULL
  names(ci_mat) <- c("ci_low", "ci_high", "term")

  # n per coefficient: for a drug-specific term (drugclass<DRUG>, or its slope
  # interaction) use the number of patients on that drug; for any shared term
  # (e.g. a prognostic factor) use the total number fitted.
  mf <- model.frame(fit)
  drug_counts <- if ("drugclass" %in% names(mf)) table(mf$drugclass) else integer(0)
  term_drug <- ifelse(
    startsWith(coef_mat$term, "drugclass"),
    sub(":.*$", "", sub("^drugclass", "", coef_mat$term)),
    NA_character_
  )
  coef_mat$n <- ifelse(
    !is.na(term_drug) & term_drug %in% names(drug_counts),
    as.integer(drug_counts[term_drug]),
    nobs(fit)
  )

  coef_mat %>%
    left_join(ci_mat, by = "term") %>%
    mutate(recalibration = model_type) %>%
    select(recalibration, term, estimate, std_error, ci_low, ci_high, n)
}

# Recalibration is done with a single pooled model per method rather than one
# model per drug. Each drug keeps its own intercept (and slope), but any
# prognostic_vars enter as a single shared main effect, so their influence on the
# outcome is the same regardless of which treatment a patient received. With
# prognostic_vars empty these pooled fits reproduce the per-drug fits exactly.

# Anchor each patient on the prediction for the drug they actually received; this
# is the offset (intercept model) / covariate (intercept + slope model).
recal_data <- analysis_dataset
recal_data$drugclass <- factor(recal_data$drugclass, levels = drugs)
pred_matrix <- as.matrix(analysis_dataset[, paste0(pred_prefix, drugs)])
recal_data$pred_assigned <- pred_matrix[
  cbind(seq_len(nrow(pred_matrix)), match(recal_data$drugclass, drugs))
]

# helper: predict the recalibrated value for every patient under a given drug
predict_for_drug <- function(fit, drug) {
  nd <- analysis_dataset
  nd$drugclass <- factor(drug, levels = drugs)
  nd$pred_assigned <- nd[[paste0(pred_prefix, drug)]]
  as.numeric(predict(fit, newdata = nd))
}

# helper: each patient's prediction (columns <prefix><drug>) for the drug they
# actually received
pred_received <- function(data, prefix) {
  pm <- as.matrix(data[, paste0(prefix, drugs)])
  pm[cbind(seq_len(nrow(pm)), match(data$drugclass, drugs))]
}

#### Recalibration-in-the-large (drug-specific intercept, shared prognostic effect)
formula_recal_large <- as.formula(
  paste0(outcome_var, " ~ 0 + drugclass", prognostic_terms, " + offset(pred_assigned)")
)
fit_recal_large <- lm(formula_recal_large, data = recal_data)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal_large.", drug)]] <- predict_for_drug(fit_recal_large, drug)
}
recalib_recal_large_params <- summarize_lm(fit_recal_large, "recal_large")

#### Minimal recalibration-in-the-large (drug-specific intercept beyond minimal_threshold, shared prognostic effect)
# The recalibration-in-the-large predictions without the intercepts within
# minimal_threshold. accepted marks the terms used (in recalibration_params_full.csv).
intercepts_recal_large <- coef(fit_recal_large)[paste0("drugclass", drugs)]
accepted <- setNames(abs(intercepts_recal_large) > minimal_threshold, drugs)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal_minimal.", drug)]] <- predict_for_drug(fit_recal_large, drug) -
    if (accepted[[drug]]) 0 else intercepts_recal_large[[paste0("drugclass", drug)]]
}
cat("Minimal recalibration-in-the-large: intercepts used for",
    if (any(accepted)) paste(drugs[accepted], collapse = ", ") else "no drug class", "\n")
recalib_recal_large_params <- recalib_recal_large_params %>%
  mutate(accepted = !startsWith(term, "drugclass") | abs(estimate) > minimal_threshold)

#### Recalibration (drug-specific intercept & slope, shared prognostic effect)
formula_recal <- as.formula(
  paste0(outcome_var, " ~ 0 + drugclass + drugclass:pred_assigned", prognostic_terms)
)
fit_recal <- lm(formula_recal, data = recal_data)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal.", drug)]] <- predict_for_drug(fit_recal, drug)
}
recalib_recal_params <- summarize_lm(fit_recal, "recal")

# drug classes not analysed (not in drugs) have no data to be recalibrated, so they
# keep their original predictions
for (drug in setdiff(model_drugs, drugs)) {
  analysis_dataset[[paste0("pred.recal_minimal.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
  analysis_dataset[[paste0("pred.recal_large.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
  analysis_dataset[[paste0("pred.recal.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
}

# Combine and save recalibration parameters
recalib_params_full <- bind_rows(recalib_recal_large_params, recalib_recal_params)

write_csv(recalib_params_full, paste0(output_path, "/recalibration_params_full.csv"))


#############################################################################

# Part 4: Figure code ----------------------------------------------------------
# This part only defines the functions that draw and save the figures. It needs
# to be run, but there is no need to read it: nothing in it needs changing.

# facet columns: one row for up to 4 panels, two rows above that
facet_ncol <- function(n_panels) if (n_panels <= 4) n_panels else ceiling(n_panels / 2)

# number of categories on a discrete axis ("x" or "y") of a plot
n_axis_levels <- function(plot, axis = "x") {
  length(ggplot2::ggplot_build(plot)$layout$panel_params[[1]][[axis]]$get_labels())
}

# how much the text grows in a figure whose panels are width inches wide in all
width_scale <- function(width) max(1, 1 + 0.2 * (width / 3.4 - 1))

# the largest font size (pt) at which every label fits on one line in width inches
fit_font_size <- function(labels, width) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  widest <- max(vapply(as.character(labels), function(label) {
    grid::convertWidth(grid::stringWidth(label), "in", valueOnly = TRUE)
  }, numeric(1)))
  # stringWidth is measured at the default 12 pt; width grows in proportion to size
  12 * width / widest
}

# width (inches) of each label at size pt. ggplot writes negative numbers with
# a minus sign, as wide as "+", not a hyphen.
text_width <- function(labels, size) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  vapply(gsub("-", "+", as.character(labels)), function(label) {
    grid::convertWidth(grid::stringWidth(label), "in", valueOnly = TRUE) * size / 12
  }, numeric(1), USE.NAMES = FALSE)
}

# the largest axis title size (pt) at which label fits on one line along
# length inches, between the 1-column size (16 pt) and 16 pt times text_scale
axis_title_size <- function(label, text_scale, length) {
  max(16, min(16 * text_scale, fit_font_size(label, length)))
}

# text wrapped onto as few lines as fit in width inches at size pt, never
# inside a drug name; text already on several lines is wrapped line by line
wrap_to_width <- function(text, size, width) {
  if (grepl("\n", text)) {
    return(paste(vapply(strsplit(text, "\n")[[1]], wrap_to_width, "", size = size, width = width),
                 collapse = "\n"))
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fits <- function(wrapped) {
    all(vapply(strsplit(wrapped, "\n")[[1]], function(line) {
      grid::convertWidth(grid::stringWidth(line), "in", valueOnly = TRUE) * size / 12 <= width
    }, logical(1)))
  }
  # drug names that must not be split across lines are joined while wrapping
  protected <- text
  for (k in c("S.c. ", "Oral Semaglutide")) protected <- gsub(k, gsub(" ", "§", k), protected, fixed = TRUE)
  restore <- function(x) gsub("§", " ", x, fixed = TRUE)
  for (n_chars in seq(nchar(protected), 10)) {
    wrapped <- restore(str_wrap(protected, n_chars))
    if (fits(wrapped)) return(wrapped)
  }
  restore(str_wrap(protected, 10))
}

# a title ending in "(step)" split into the title and the step (capitalised),
# so the step can go in the subtitle; step is NULL when there is none
split_title_step <- function(title) {
  step <- regmatches(title, regexec(" [(]([^()]*)[)]$", title))[[1]]
  if (length(step) != 2) return(list(title = title, step = NULL))
  list(title = sub(" [(][^()]*[)]$", "", title),
       step  = paste0(toupper(substr(step[2], 1, 1)), substring(step[2], 2)))
}

# how far (inches) the outermost x axis numbers stick out past the left and
# right edges of a panel panel_width inches wide, at size pt, for an axis
# trained on the values x (expanded by ggplot's default 5% each side)
x_overhang <- function(x, size, panel_width) {
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  lo_e <- lo - 0.05 * (hi - lo); hi_e <- hi + 0.05 * (hi - lo)
  breaks <- scales::extended_breaks()(c(lo, hi))
  breaks <- breaks[breaks >= lo_e & breaks <= hi_e]
  labels <- format(breaks, trim = TRUE)
  w <- text_width(labels[c(1, length(labels))], size)
  per_unit <- panel_width / (hi_e - lo_e)
  c(left  = max(0, w[1] / 2 - (min(breaks) - lo_e) * per_unit),
    right = max(0, w[2] / 2 - (hi_e - max(breaks)) * per_unit))
}

# the drug class names of one figure, given the drug classes it shows: with only
# one semaglutide form in the figure, it is called "Semaglutide" (the figure
# says which form); with both, they keep "Oral" / "S.c."
figure_labels <- function(figure_drugs) {
  labels <- drug_labels
  sema   <- intersect(c("Oral_Semaglutide", "Inje_Semaglutide"), names(labels))
  if (sum(sema %in% figure_drugs) == 1) labels[sema] <- "Semaglutide"
  labels
}

# save a ggplot or patchwork with panels of panel_width x panel_height inches
# (pdf, or png at 300 dpi for slides). Either size can be a function of the
# plot, e.g. to give each bar or each covariate the same room.
save_plot <- function(plot, file, panel_width, panel_height) {
  if (is.function(panel_width))  panel_width  <- panel_width(plot)
  if (is.function(panel_height)) panel_height <- panel_height(plot)
  grDevices::pdf(NULL)
  g <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)

  # panels: the relative ("null") rows and columns get their size in inches
  null_w <- grid::unitType(g$widths) == "null"
  null_h <- grid::unitType(g$heights) == "null"
  g$widths[null_w]  <- grid::unit(panel_width  * as.numeric(g$widths[null_w]), "in")
  g$heights[null_h] <- grid::unit(panel_height * as.numeric(g$heights[null_h]), "in")

  inches <- function(u) grid::convertWidth(u, "in", valueOnly = TRUE)
  col_edges <- c(0, cumsum(vapply(seq_along(g$widths), function(i) inches(g$widths[i]), numeric(1))))
  width <- max(col_edges)

  # titles, captions and legends wider than the plot widen the page instead of
  # being cut off (titles start at their left edge, captions end at their
  # right edge, legends are centred)
  pad_left <- 0; pad_right <- 0
  for (i in seq_len(nrow(g$layout))) {
    name <- g$layout$name[i]
    if (!grepl("^(title|subtitle|caption|guide-box)", name)) next
    grob_width <- inches(grid::grobWidth(g$grobs[[i]]))
    left  <- col_edges[g$layout$l[i]]
    right <- col_edges[g$layout$r[i] + 1]
    if (grepl("^caption", name)) {
      pad_left <- max(pad_left, grob_width - right)
    } else if (grepl("^guide-box", name)) {
      centre <- (left + right) / 2
      pad_left  <- max(pad_left, grob_width / 2 - centre)
      pad_right <- max(pad_right, centre + grob_width / 2 - width)
    } else {
      pad_right <- max(pad_right, left + grob_width - width)
    }
  }
  # a little room on the right for the last x axis label, which overhangs its panel
  pad_right <- max(pad_right, 0.15)
  if (pad_left > 0)  g <- gtable::gtable_add_cols(g, grid::unit(pad_left, "in"), pos = 0)
  if (pad_right > 0) g <- gtable::gtable_add_cols(g, grid::unit(pad_right, "in"), pos = -1)

  width  <- width + pad_left + pad_right
  height <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  grDevices::dev.off()

  if (grepl("[.]png$", file)) {
    ragg::agg_png(file, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  }
  grid::grid.draw(g)
  invisible(grDevices::dev.off())
}


## Pairwise drug calibration figure ----

# One square panel (3.4 inches) per pair of drug classes, labelled by `title`,
# all on one symmetric range: observed against predicted benefit, one point per
# group of predicted benefit, or a single point for a pooled pair. A trailing
# "(step)" in the title goes in the subtitle.
pairwise_plot <- function(plot_data, plot_title, plot_subtitle = NULL) {
  parts <- split_title_step(plot_title)
  if (is.null(plot_subtitle)) plot_subtitle <- parts$step
  plot_title <- parts$title
  panel_size <- 3.4
  n_panels   <- n_distinct(plot_data$title)
  n_cols     <- facet_ncol(n_panels)
  n_rows     <- ceiling(n_panels / n_cols)
  text_scale <- 1 + 0.2 * (n_cols - 1)

  # Facet titles: as large as the longest drug pair of this figure allows on one
  # line within a panel, but no larger than the axis titles; on two lines
  # ("Oral Semaglutide" / "vs SGLT2i") when one line would make them smaller
  # than 75% of that
  # panel titles in the order they first appear (a factor keeps its own order)
  if (!is.factor(plot_data$title)) plot_data$title <- factor(plot_data$title, levels = unique(plot_data$title))
  strip_size <- min(16 * text_scale, fit_font_size(levels(plot_data$title), panel_size - 0.25))
  if (strip_size < 0.75 * 16 * text_scale) {
    levels(plot_data$title) <- sub(" vs ", "\nvs ", levels(plot_data$title))
    strip_size <- min(16 * text_scale,
                      fit_font_size(unlist(strsplit(levels(plot_data$title), "\n")), panel_size - 0.25))
  }
  title_width <- n_cols * panel_size + 0.8

  # gap between panel columns: the outermost x axis numbers overhang the panel
  # edges (the axis runs 10% beyond the data); the gap fits the overhang of the
  # right number of one panel and the left number of the next, plus 0.12 inches
  axis_text_size <- 15 * text_scale
  lim       <- max(abs(c(plot_data$coef_low, plot_data$coef_high, plot_data$mean)), na.rm = TRUE)
  axis_end  <- 1.1 * lim
  breaks    <- scales::extended_breaks()(c(-lim, lim))
  breaks    <- breaks[abs(breaks) <= axis_end]
  edge_room <- (axis_end - max(abs(breaks))) / (2 * axis_end) * panel_size
  overhang  <- pmax(0, text_width(format(range(breaks), trim = TRUE), axis_text_size) / 2 - edge_room)
  panel_gap <- sum(overhang) + 0.12

  ggplot(plot_data, aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
    geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
    # error bars without caps; points and error bars grow with the text
    geom_errorbar(width = 0, colour = "#800020", linewidth = 0.5 * text_scale) +
    geom_point(colour = "#800020", size = 1.8 * text_scale) +
    facet_wrap(~title, ncol = n_cols) +
    # inject the predicted (mean) and observed-CI (coef_low/coef_high) values onto
    # both axes so x and y share one symmetric range; with coord_fixed this gives a
    # square panel of identical shape in every plot
    geom_blank(aes(x = coef_low, y = coef_low)) +
    geom_blank(aes(x = coef_high, y = coef_high)) +
    geom_blank(aes(x = mean, y = mean)) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
    scale_y_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
    theme_classic() +
    theme(
      # axis titles on one line along the height of the panels
      axis.title = element_text(size = axis_title_size("Observed HbA1c benefit (mmol/mol)", text_scale,
                                                       n_rows * (panel_size + 0.35))),
      axis.text = element_text(size = axis_text_size),
      # facet titles without a box
      strip.text = element_text(size = strip_size),
      strip.background = element_blank(),
      plot.title = element_text(size = 20 * text_scale),
      plot.subtitle = element_text(size = 17 * text_scale),
      plot.title.position = "plot",
      panel.spacing.x = grid::unit(panel_gap, "in"), panel.spacing.y = grid::unit(0.6, "lines")
    ) +
    labs(
      x = expression("Predicted HbA"[1*c]*" benefit (mmol/mol)"),
      y = expression("Observed HbA"[1*c]*" benefit (mmol/mol)"),
      # wrapped only where they are wider than the figure (panels plus the y axis)
      title = wrap_to_width(plot_title, 20 * text_scale, title_width),
      subtitle = if (!is.null(plot_subtitle)) wrap_to_width(plot_subtitle, 17 * text_scale, title_width)
    )
}

# helper: save a pairwise figure (nothing is saved when it was not made)
save_pairwise <- function(fig, file) {
  if (!is.null(fig)) save_plot(fig, file, panel_width = 3.4, panel_height = 3.4)
}


## Overall benefit calibration figures ----

# Points are the mean observed benefit per group (from benefit_by_group() /
# benefit_by_drug()), with a confidence interval, against the mean predicted
# benefit. The dashed blue curve is a LOESS fit through these group point
# estimates themselves (not the individual matched pairs), so the curve
# summarises exactly the points shown, not the noisier patient-level data.
# The `n` column from benefit_by_group() is the number of matched pairs in each
# group, so its sum is the total number of pairs reported in the subtitle.
# When several cal_groups settings are present (an n_groups column), the figure
# shows the largest number of groups whose groups all contain at least `min_n`
# matched pairs (or else the smallest number of groups); the number of groups is
# not written, as the points show it.
#   text_scale: how much the text, points and lines grow
#   title_size, subtitle_size, axis_title_size: font sizes (pt)
#   panel_width: if given, the plot gets room on its right for its outermost x
#     axis number, so it does not run into the plot next to it
plot_calibration <- function(deciles, title = NULL, min_n = 100, text_scale = 1,
                             title_size = 20 * text_scale, subtitle_size = 17 * text_scale,
                             axis_title_size = 16 * text_scale, panel_width = NULL) {
  if ("n_groups" %in% names(deciles)) {
    ok <- deciles %>%
      group_by(n_groups) %>%
      summarise(.min_group_n = min(n, na.rm = TRUE), .groups = "drop")
    chosen <- if (any(ok$.min_group_n >= min_n)) {
      max(ok$n_groups[ok$.min_group_n >= min_n])
    } else {
      min(ok$n_groups)
    }
    deciles <- deciles %>% filter(n_groups == chosen)
  }

  axis_text_size <- 15 * text_scale
  # the x axis also covers 0 (the dashed vertical line)
  right_margin <- if (!is.null(panel_width)) {
    x_overhang(c(deciles$pred.benefit, 0), axis_text_size, panel_width)[["right"]] + 0.12
  } else {
    18 / 72
  }

  n_pairs <- sum(deciles$n, na.rm = TRUE)
  ggplot(deciles, aes(x = pred.benefit, y = obs.benefit)) +
    geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    # LOESS only where there are more than 5 group points (in each drug class);
    # through 5 or fewer its confidence band is unstable and can swamp the axis
    geom_smooth(data = function(d) d %>% group_by(across(any_of("index_drug"))) %>% filter(n() > 5) %>% ungroup(),
                method = "loess", formula = y ~ x, span = 1, se = TRUE, colour = "blue",
                linetype = "dashed", linewidth = 0.75 * text_scale) +
    # error bars without caps; points and error bars grow with the text
    geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0, colour = "black",
                  linewidth = 0.5 * text_scale) +
    geom_point(size = 1.8 * text_scale) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      plot.title = element_text(size = title_size),
      plot.subtitle = element_text(size = subtitle_size),
      plot.margin = margin(5.5, right_margin * 72, 5.5, 5.5, unit = "pt")
    ) +
    labs(
      x = expression("Predicted HbA"[1*c]*" benefit (mmol/mol)"),
      y = expression("Observed HbA"[1*c]*" benefit (mmol/mol)"),
      title = title,
      subtitle = paste0("n = ", format(n_pairs, big.mark = ","), " matched pairs")
    )
}

# the concordant, discordant and combined matched pairs side by side (panels of
# 4.2 inches), with one x and one y axis title for the three; a trailing
# "(step)" in the title goes in the subtitle
overall_figure <- function(deciles_conc, deciles_disc, deciles_both, title) {
  panel_size <- 4.2
  text_scale <- 1 + 0.2 * (3 - 1)

  # panel titles: as large as the longest allows within a panel, but no larger
  # than the axis titles; the numbers of pairs (subtitles) no larger than them
  panel_titles <- c("Concordant matched to discordant", "Discordant matched to concordant", "Combined")
  title_size   <- min(16 * text_scale, fit_font_size(panel_titles, panel_size - 0.25))
  sizes <- list(text_scale = text_scale, title_size = title_size,
                subtitle_size = min(15 * text_scale, title_size),
                axis_title_size = axis_title_size("Observed HbA1c benefit (mmol/mol)", text_scale, panel_size),
                panel_width = panel_size)
  panel <- function(deciles, panel_title) do.call(plot_calibration, c(list(deciles, title = panel_title), sizes))

  parts <- split_title_step(title)
  panel(deciles_conc, panel_titles[1]) +
    panel(deciles_disc, panel_titles[2]) +
    panel(deciles_both, panel_titles[3]) +
    patchwork::plot_layout(nrow = 1, axis_titles = "collect") +
    patchwork::plot_annotation(
      title = wrap_to_width(parts$title, 20 * text_scale, 3 * panel_size + 1.5),
      subtitle = parts$step,
      theme = theme(plot.title = element_text(size = 20 * text_scale),
                    plot.subtitle = element_text(size = 17 * text_scale))
    )
}

# same plot, one panel (3.4 inches) per drug class that anchors the pairs (the
# index / matched patient's drug); the LOESS curve is fitted separately within
# each panel, from that drug's own group points. Each facet title is the drug
# class and its number of matched pairs, summing the per-group `n` returned by
# benefit_by_drug(). The number of groups is chosen per drug class.
plot_calibration_by_drug <- function(by_drug, title = NULL, min_n = 100) {
  if ("n_groups" %in% names(by_drug)) {
    by_drug <- by_drug %>%
      group_by(index_drug, n_groups) %>%
      mutate(.min_group_n = min(n, na.rm = TRUE)) %>%
      group_by(index_drug) %>%
      mutate(.chosen_groups = if (any(.min_group_n >= min_n)) {
        max(n_groups[.min_group_n >= min_n])
      } else {
        min(n_groups)
      }) %>%
      ungroup() %>%
      filter(n_groups == .chosen_groups) %>%
      select(-.min_group_n, -.chosen_groups, -n_groups)
  }

  panel_size <- 3.4
  panel_drugs <- levels(droplevels(factor(by_drug$index_drug)))
  n_panels   <- length(panel_drugs)
  n_cols     <- facet_ncol(n_panels)
  n_rows     <- ceiling(n_panels / n_cols)
  text_scale <- 1 + 0.2 * (n_cols - 1)

  # Facet titles: the drug class (named as in this figure) over its number of
  # matched pairs; as large as the longest line allows within a panel, but no
  # larger than the axis titles
  labels <- figure_labels(panel_drugs)
  facet_counts <- by_drug %>%
    group_by(index_drug) %>%
    summarise(n_pairs = sum(n, na.rm = TRUE), .groups = "drop")
  facet_labels <- setNames(
    paste0(labels[as.character(facet_counts$index_drug)],
           "\n(n = ", format(facet_counts$n_pairs, big.mark = ",", trim = TRUE), ")"),
    as.character(facet_counts$index_drug)
  )
  strip_size <- min(16 * text_scale, fit_font_size(unlist(strsplit(facet_labels, "\n")), panel_size - 0.25))

  # gap between panel columns for the outermost x axis numbers (the same axis in
  # every panel), plus 0.12 inches
  overhang  <- x_overhang(c(by_drug$pred.benefit, 0), 15 * text_scale, panel_size)
  panel_gap <- sum(overhang) + 0.12

  # the step goes in the subtitle, before the number of matched pairs
  parts <- split_title_step(title)
  subtitle <- paste0(if (!is.null(parts$step)) paste0(parts$step, "; "),
                     "n = ", format(sum(by_drug$n, na.rm = TRUE), big.mark = ","), " matched pairs")
  title_width <- n_cols * panel_size + 0.8

  plot_calibration(by_drug, text_scale = text_scale,
                   axis_title_size = axis_title_size("Observed HbA1c benefit (mmol/mol)", text_scale,
                                                     n_rows * (panel_size + 0.6))) +
    facet_wrap(~ index_drug, ncol = n_cols, labeller = labeller(index_drug = facet_labels)) +
    theme(
      # facet titles without a box
      strip.text = element_text(size = strip_size),
      strip.background = element_blank(),
      panel.spacing.x = grid::unit(panel_gap, "in"),
      plot.title.position = "plot"
    ) +
    labs(title = wrap_to_width(parts$title, 20 * text_scale, title_width),
         subtitle = wrap_to_width(subtitle, 17 * text_scale, title_width))
}

# helpers: save the overall benefit figures
save_overall <- function(fig, file) save_plot(fig, file, panel_width = 4.2, panel_height = 4.2)
save_overall_by_drug <- function(fig, file) save_plot(fig, file, panel_width = 3.4, panel_height = 3.4)


## Predicted vs observed figure ----

# One point per initiation: the observed HbA1c (at 12 months, or at 6 months when
# the 12-month value is missing) against the prediction for the drug class
# initiated; all in one panel (5.5 x 4.2 inches), or one panel per drug class
# (by_drug, 3.6 inches). pred_obs_data has one row per initiation: drugclass,
# pred_drug and posthba1cfinal. limits: list(x = , y = ) of the axis ranges,
# the same for every predicted vs observed figure of the script
# (pred_obs_limits()), so the figures of every step can be compared.
draw_pred_obs <- function(pred_obs_data, pred_drugs, title, step_title, by_drug = FALSE,
                          limits = pred_obs_limits(list(pred_obs_data))) {
  if (length(pred_drugs) == 0) return(NULL)

  panel_width  <- if (by_drug) 3.6 else 5.5
  panel_height <- if (by_drug) 3.6 else 4.2
  n_cols <- if (by_drug) facet_ncol(length(pred_drugs)) else 1
  n_rows <- if (by_drug) ceiling(length(pred_drugs) / n_cols) else 1
  text_scale <- if (by_drug) 1 + 0.2 * (n_cols - 1) else width_scale(panel_width)
  axis_text_size <- 15 * text_scale

  # Facet titles: the drug class (named as in this figure) and its number of
  # initiations, as large as the longest allows within a panel but no larger
  # than the axis titles; on two lines ("Drug" / "(n = ...)") when one line
  # would make them smaller than 75% of that
  labels <- figure_labels(pred_drugs)
  n_text <- format(as.integer(table(factor(pred_obs_data$drugclass, levels = pred_drugs))), big.mark = ",", trim = TRUE)
  panel_labels <- setNames(paste0(labels[pred_drugs], " (n = ", n_text, ")"), pred_drugs)
  strip_size <- min(16 * text_scale, fit_font_size(panel_labels, panel_width - 0.25))
  if (strip_size < 0.75 * 16 * text_scale) {
    panel_labels <- setNames(paste0(labels[pred_drugs], "\n(n = ", n_text, ")"), pred_drugs)
    strip_size <- min(16 * text_scale, fit_font_size(unlist(strsplit(panel_labels, "\n")), panel_width - 0.25))
  }
  pred_obs_data$panel <- factor(panel_labels[as.character(pred_obs_data$drugclass)], levels = panel_labels)

  # room between panels, and on the right, for the outermost x axis numbers
  overhang  <- x_overhang(limits$x, axis_text_size, panel_width)
  panel_gap <- sum(overhang) + 0.12
  title_width <- n_cols * panel_width + 0.8

  ggplot(pred_obs_data, aes(x = pred_drug, y = posthba1cfinal)) +
    geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
    # points and the fitted curve grow with the text
    geom_point(alpha = 0.4, size = 1.2 * text_scale) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, linewidth = 1 * text_scale) +
    { if (by_drug) facet_wrap(~ panel, ncol = n_cols) } +
    # the same axes in every predicted vs observed figure (zooms only: no point
    # or curve is dropped)
    coord_cartesian(xlim = limits$x, ylim = limits$y) +
    theme_classic() +
    theme(
      # axis titles on one line along the height of the panels
      axis.title = element_text(size = axis_title_size("Observed HbA1c (mmol/mol)", text_scale,
                                                       n_rows * (panel_height + 0.4))),
      axis.text = element_text(size = axis_text_size),
      # facet titles without a box
      strip.text = element_text(size = strip_size),
      strip.background = element_blank(),
      plot.title = element_text(size = 20 * text_scale),
      plot.subtitle = element_text(size = 17 * text_scale),
      plot.title.position = "plot",
      panel.spacing.x = grid::unit(panel_gap, "in"),
      plot.margin = margin(5.5, (overhang[["right"]] + 0.12) * 72, 5.5, 5.5, unit = "pt")
    ) +
    labs(
      x = expression("Predicted HbA"[1*c]*" (mmol/mol)"),
      y = expression("Observed HbA"[1*c]*" (mmol/mol)"),
      title = wrap_to_width(title, 20 * text_scale, title_width),
      subtitle = wrap_to_width(paste0(step_title, "; n = ", format(nrow(pred_obs_data), big.mark = ","), " initiations"),
                               17 * text_scale, title_width)
    )
}

# helper: the axis ranges covering every data frame in data_list (each with
# pred_drug and posthba1cfinal)
pred_obs_limits <- function(data_list) {
  all <- bind_rows(lapply(data_list, function(d) d[, c("pred_drug", "posthba1cfinal")]))
  list(x = range(all$pred_drug, na.rm = TRUE), y = range(all$posthba1cfinal, na.rm = TRUE))
}

# helper: save a predicted vs observed figure (nothing is saved when it was not made)
save_pred_obs <- function(fig, file, by_drug = FALSE) {
  if (!is.null(fig)) save_plot(fig, file, panel_width  = if (by_drug) 3.6 else 5.5,
                                          panel_height = if (by_drug) 3.6 else 4.2)
}


## Predicted optimal therapy figures ----

# width (inches) of a bar chart: 1.25 inches per bar, at least 3
bar_panel_width <- function(n_bars) max(3, 1.25 * n_bars)

# Distribution of the predicted best drug class (rank 1): vertical bars, largest
# first, without outlines on plain axes, with the percentage above each bar and
# the number of initiations under it in grey. rank1_breakdown has one row per
# drug class that is best for anyone: drug, n, proportion. analysis_drugs: the
# drug classes the best is chosen from (for their names in this figure).
draw_rank1 <- function(rank1_breakdown, plot_subtitle, analysis_drugs = rank1_breakdown$drug) {
  width      <- bar_panel_width(nrow(rank1_breakdown))
  text_scale <- width_scale(width)
  axis_text_size <- 15 * text_scale
  labels <- figure_labels(analysis_drugs)

  rank1_breakdown %>%
    mutate(drug = factor(drug, levels = drug[order(-n)]),
           pct_label = sprintf("%.1f%%", proportion),
           n_label   = paste0("n = ", format(n, big.mark = ",", trim = TRUE))) %>%
    ggplot(aes(x = drug, y = proportion, fill = drug)) +
    geom_col(width = 0.7) +
    # number of initiations (grey, smaller) just above the bar, percentage above it
    geom_text(aes(label = n_label), vjust = -0.5, size = 0.8 * axis_text_size / .pt, colour = "grey40") +
    geom_text(aes(label = pct_label), vjust = -1.75, size = axis_text_size / .pt, fontface = "bold") +
    scale_fill_manual(values = drug_colours, guide = "none") +
    # drug class names wrapped to the width of a bar, never inside a drug name
    scale_x_discrete(labels = function(x) vapply(labels[x], wrap_to_width, "", size = axis_text_size, width = 1.15)) +
    scale_y_continuous(limits = c(0, max(rank1_breakdown$proportion) * 1.3), expand = expansion(mult = c(0, 0)),
                       breaks = function(lim) { b <- pretty(lim); b[b <= 100] }) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 16 * text_scale),
      axis.text = element_text(size = axis_text_size, colour = "black"),
      # drug class names a little below the axis line
      axis.text.x = element_text(margin = margin(t = 8 * text_scale)),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 20 * text_scale),
      plot.subtitle = element_text(size = 17 * text_scale),
      plot.title.position = "plot"
    ) +
    labs(
      x = NULL,
      y = "Initiations (%)",
      title = wrap_to_width("Predicted best drug class", 20 * text_scale, width + 0.8),
      subtitle = wrap_to_width(plot_subtitle, 17 * text_scale, width + 0.8)
    )
}

# Drug classes added to and removed from the tolerance-3 sets by recalibration.
# changes has one row per change type ("added", "removed") and drug class:
# change_type, drug, n; n_total is the number of patients (for the percentage
# axis). analysis_drugs: the drug classes the sets are chosen from.
draw_tolerance3 <- function(changes, n_total, plot_subtitle, analysis_drugs = unique(changes$drug)) {
  width      <- bar_panel_width(n_distinct(changes$drug))
  text_scale <- width_scale(width)
  axis_text_size <- 15 * text_scale
  labels <- figure_labels(analysis_drugs)

  changes %>%
    mutate(
      drug = factor(drug, levels = intersect(model_drugs, drug)),
      signed_n = ifelse(change_type == "removed", -n, n)
    ) %>%
    ggplot(aes(x = drug, y = signed_n, fill = drug, alpha = change_type)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = drug_colours, guide = "none") +
    scale_alpha_manual(values = c("added" = 1, "removed" = 0.6)) +
    # drug class names wrapped to the width of a bar, never inside a drug name
    scale_x_discrete(labels = function(x) vapply(labels[x], wrap_to_width, "", size = axis_text_size, width = 1.15)) +
    scale_y_continuous(
      name = "Net count (n)",
      labels = scales::comma,
      sec.axis = sec_axis(
        ~ . * 100 / n_total,
        name = "Net percentage (%)",
        labels = function(x) sprintf("%.1f", x)
      )
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 16 * text_scale),
      axis.text = element_text(size = axis_text_size),
      plot.title = element_text(size = 20 * text_scale),
      plot.subtitle = element_text(size = 17 * text_scale),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.title = element_text(size = 16 * text_scale),
      legend.text = element_text(size = axis_text_size)
    ) +
    labs(
      x = "Drug Class",
      title = wrap_to_width("Tolerance-3 drug changes", 20 * text_scale, width + 1.6),
      subtitle = wrap_to_width(plot_subtitle, 17 * text_scale, width + 1.6),
      alpha = "Change"
    )
}

# helper: save a predicted optimal therapy figure (each bar gets the same room)
save_decision <- function(fig, file) {
  if (!is.null(fig)) save_plot(fig, file, panel_width = function(p) bar_panel_width(n_axis_levels(p)), panel_height = 3.6)
}


## Overall benefit estimates figure ----

# The overall benefit estimates (rows with key_result) in each step, one
# panel (5.5 x 1.6 inches) per analysis: observed, lower.ci, upper.ci, model
# (step) and analysis_label.
draw_benefit_estimates <- function(benefit_estimates, plot_subtitle = NULL,
                                   plot_title = "Overall benefit estimates") {
  panel_width <- 5.5
  text_scale  <- width_scale(panel_width)
  # first letter in capitals, the rest as written ("S.c. Semaglutide" stays)
  capitalise <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
  benefit_estimates <- benefit_estimates %>%
    filter(key_result) %>%
    mutate(
      model = factor(model, levels = c("Recalibration", "Recalibration-in-the-large",
                                       "Minimal recalibration-in-the-large", "Original")),
      analysis_label = capitalise(as.character(analysis_label)),
      analysis_label = factor(analysis_label, levels = unique(analysis_label))
    )

  # facet titles: as large as the longest allows within the panel, but no
  # larger than the axis titles
  strip_size  <- min(16 * text_scale, fit_font_size(levels(benefit_estimates$analysis_label), panel_width - 0.25))
  # titles can take the panel plus the step names on the left
  title_width <- panel_width + 2.6

  ggplot(benefit_estimates, aes(x = observed, xmin = lower.ci, xmax = upper.ci, y = model, colour = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    # points and lines grow with the text
    geom_pointrange(size = 0.7 * text_scale, linewidth = 0.5 * text_scale) +
    facet_wrap(~ analysis_label, ncol = 1) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 16 * text_scale),
      axis.text = element_text(size = 15 * text_scale),
      # facet titles without a box
      strip.text = element_text(size = strip_size),
      strip.background = element_blank(),
      plot.title = element_text(size = 20 * text_scale),
      plot.subtitle = element_text(size = 17 * text_scale),
      plot.title.position = "plot"
    ) +
    labs(
      x = expression("Observed HbA"[1*c]*" benefit estimate (mmol/mol)"),
      y = NULL,
      title = wrap_to_width(plot_title, 20 * text_scale, title_width),
      subtitle = if (!is.null(plot_subtitle)) wrap_to_width(plot_subtitle, 17 * text_scale, title_width)
    )
}

# helper: save the overall benefit estimates figure
save_benefit_estimates <- function(fig, file) save_plot(fig, file, panel_width = 5.5, panel_height = 1.6)


#############################################################################

# Part 5: Validation -----------------------------------------------------------

## The validation is done for each version of the predictions in turn (Parts
## 5.1 to 5.4), and then summarised (Part 5.5). Below, the analyses are set up
## and the code for each kind of figure is defined; Parts 5.1 to 5.5 run them.
##
## The prognostic factor (prognostic_vars) is a cross-cutting control applied
## within every version but the original when identified - adjusted for in the
## pairwise analysis and matched on in the overall-benefit analysis - not a
## version of its own.

## Pairwise drug calibration ----

## Interpretation of plots:
# - Points close to the red line (y=x) indicate good calibration
# - Points above the red line indicate underprediction (observed benefit greater than predicted)
# - Points below the red line indicate overprediction (observed benefit less than predicted)


# This methods has a few entries that can be adjusted based on your dataset:
# - cal_groups: Number of groups to stratify predictions into (e.g., 3, 5, 10)
#   Ensure that each group has sufficient patients for stable estimates (e.g., >= 100 per drug)
#   This can be checked in the output summary tables

cal_groups = c(3, 5, 10) # Number of groups to stratify predictions into

# - adjustment_var: baseline covariates to adjust for in the pairwise calibration
#   models. These should match the variables used in original model development.
#   Any additional prognostic factors go in prognostic_vars (above) and are added
#   on top (as an overall, treatment-independent effect) in the recalibrated
#   versions (adjustment_vars_prog), not the original.

adjustment_vars = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt",
                     "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth",
                     "sex", "smoke", "ncurrtx", "drugline")
adjustment_vars = c(adjustment_vars, "imd5") # SDRN run: each patient's own deprivation (SIMD quintile)
adjustment_vars_prog <- unique(c(adjustment_vars, prognostic_vars))

# - Pairwise figure: the drug classes analysed against each other

# helper: pairwise calibration plot, one panel per pair of drug classes
plot_pairwise <- function(pairwise_results, plot_title, plot_subtitle) {
  pairwise_results %>%
    group_by(drug1, drug2, n_groups) %>%
    mutate(min_val = min(n_drug1, n_drug2)) %>%
    group_by(drug1, drug2) %>%
    mutate(
      # the largest number of groups with at least 100 patients per drug in every group
      best_grouping = suppressWarnings(max(n_groups[min_val >= 100])),
      select_grouping = ifelse(is.finite(best_grouping), best_grouping, min(cal_groups))
    ) %>%
    ungroup() %>%
    filter(n_groups == select_grouping) %>%
    mutate(title = paste0(recode(drug1, !!!drug_labels), " vs ", recode(drug2, !!!drug_labels))) %>%
    pairwise_plot(plot_title, plot_subtitle)
}

# helper: pairwise calibration of every pair of the drug classes analysed, with the
# predictions with the given prefix
pairwise_results_of <- function(prefix) {
  cateval::pairwise_calibration(
    data = analysis_dataset,
    drug_var = "drugclass",
    drugs = drugs,
    prediction_vars = paste0(prefix, drugs),
    outcome_var = outcome_var,
    cal_groups = cal_groups,
    adjustment_var = if (prefix == "pred.") adjustment_vars else adjustment_vars_prog
  )
}


## Overall benefit calibration ----

## This compares patients who received the model's recommended drug (concordant)
## with clinically similar patients who received a different drug (discordant).
## Concordant patients are matched 1:1 to discordant patients, in both
## directions, and for each matched pair we compare the predicted benefit of the
## recommended drug against the benefit actually observed.
##
## Interpretation of plots:
# - Points close to the red line (y=x) indicate good calibration
# - Points above the red line indicate underprediction (observed benefit greater than predicted)
# - Points below the red line indicate overprediction (observed benefit less than predicted)

# The recommended (predicted best) drug class is chosen from the drug classes
# analysed (drugs), as every patient needs patients on their best drug class to be
# matched with.

# Covariates used to match concordant and discordant patients. Matching is 1:1
# with replacement and exact on the predicted-best drug (done inside the
# function), with a 0.4 SD caliper on baseline HbA1c. These are baseline
# characteristics only: unlike adjustment_vars above, hba1cmonth (time to
# outcome measurement) describes the outcome window rather than the patient at
# baseline, so it is left out of matching.
matching_vars <- c(
  "t2dmduration", "prebmi", "prehba1c", "agetx", "prealt",
  "preegfr", "pretotalcholesterol", "prehdl",
  "sex", "smoke", "ncurrtx", "drugline"
)
matching_vars <- c(matching_vars, "imd5") # SDRN run: each patient's own deprivation (SIMD quintile)
# The same covariates plus any prognostic factors in prognostic_vars. Used for
# matching in every version of the predictions but the original (the prognostic factor is a
# cross-cutting control, not a version of its own); equal to matching_vars until you populate
# prognostic_vars.
matching_vars_prog <- unique(c(matching_vars, prognostic_vars))

# Confidence intervals: FALSE uses fast analytical intervals; set to TRUE for
# pair-level bootstrap intervals (slower) that account for matching with
# replacement.
use_bootstrap <- FALSE
n_boot <- 1000
boot_seed <- 42

# Covariate names as they should read in the love plots (plotmath, so units
# keep their superscripts and subscripts)
love_labels <- c(
  "t2dmduration"        = "'Duration of diabetes (years)'",
  "prebmi"              = "'BMI (kg/m'^2*')'",
  "prehba1c"            = "'HbA'[1*c]*' (mmol/mol)'",
  "agetx"               = "'Current age (years)'",
  "prealt"              = "'ALT (IU/l)'",
  "preegfr"             = "'eGFR (ml/min per 1.73 m'^2*')'",
  "pretotalcholesterol" = "'Total cholesterol (mmol/l)'",
  "prehdl"              = "'HDL-cholesterol (mmol/l)'",
  "sexMale"             = "'Male'",
  "sexFemale"           = "'Female'",
  "smokeNon-smoker"     = "'Non-smoker'",
  "smokeActive smoker"  = "'Active smoker'",
  "smokeEx-smoker"      = "'Ex-smoker'",
  "smokeNot recorded"   = "'Not recorded'",
  "ncurrtx1"            = "'Number of other current glucose-lowering drugs = 1'",
  "ncurrtx2"            = "'Number of other current glucose-lowering drugs = 2'",
  "ncurrtx3"            = "'Number of other current glucose-lowering drugs = 3'",
  "ncurrtx4+"           = "'Number of other current glucose-lowering drugs = >3'",
  "drugline2"           = "'Number of previous glucose-lowering drugs = 2'",
  "drugline3"           = "'Number of previous glucose-lowering drugs = 3'",
  "drugline4"           = "'Number of previous glucose-lowering drugs = 4'",
  "drugline5+"          = "'Number of previous glucose-lowering drugs = >4'"
)
# SDRN run: deprivation (matched on above), the SIMD quintile
love_labels <- c(love_labels,
  "imd51 (least)" = "'Scottish Index of Multiple Deprivation = 1 (least deprived)'",
  "imd52"         = "'Scottish Index of Multiple Deprivation = 2'",
  "imd53"         = "'Scottish Index of Multiple Deprivation = 3'",
  "imd54"         = "'Scottish Index of Multiple Deprivation = 4'",
  "imd55 (most)"  = "'Scottish Index of Multiple Deprivation = 5 (most deprived)'"
)

# Balance (love) plot of one matching: absolute standardised mean
# difference per matching covariate, before and after matching, in each
# matching direction. Exact-matched variables are left out.
plot_loveplot <- function(match_obj) {

  # covariates with a single value cannot be summarised and are left out: e.g.
  # the predicted best drug class, when it is the same for every patient
  love_summary <- function(m) {
    m$X <- m$X[, vapply(m$X, function(x) length(unique(x)) > 1, logical(1)), drop = FALSE]
    summary(m)
  }
  love_conc <- love_summary(match_obj$matchit$concordant)
  love_disc <- love_summary(match_obj$matchit$discordant)

  love_n <- c(
    "Concordant as treated" = nrow(match_obj$concordant),
    "Discordant as treated" = nrow(match_obj$discordant)
  )
  love_strips <- paste0(
    names(love_n), "\nn = ", format(love_n, big.mark = ",", trim = TRUE), " pairs"
  )

  love_data <- bind_rows(
    as.data.frame(love_conc$sum.all) %>%
      rownames_to_column("covariate") %>%
      mutate(direction = "Concordant as treated", sample = "Before matching"),
    as.data.frame(love_conc$sum.matched) %>%
      rownames_to_column("covariate") %>%
      mutate(direction = "Concordant as treated", sample = "After matching"),
    as.data.frame(love_disc$sum.all) %>%
      rownames_to_column("covariate") %>%
      mutate(direction = "Discordant as treated", sample = "Before matching"),
    as.data.frame(love_disc$sum.matched) %>%
      rownames_to_column("covariate") %>%
      mutate(direction = "Discordant as treated", sample = "After matching")
  ) %>%
    filter(!str_detect(covariate, "^(cateval_best_drug|drugclass)")) %>%
    mutate(
      smd = abs(`Std. Mean Diff.`),
      covariate = factor(covariate, levels = rev(unique(covariate))),
      direction = factor(direction, levels = names(love_n), labels = love_strips),
      sample = factor(sample, levels = c("Before matching", "After matching"))
    )

  ggplot(love_data, aes(x = smd, y = covariate, shape = sample, colour = sample)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_vline(xintercept = c(0.05, 0.1), linetype = "dashed", colour = "grey60") +
    geom_point(size = 2.6, na.rm = TRUE) +
    facet_wrap(~ direction) +
    scale_shape_manual(values = c("Before matching" = 1, "After matching" = 16)) +
    scale_colour_manual(values = c("Before matching" = "grey45", "After matching" = "black")) +
    scale_y_discrete(labels = function(x) {
      out <- love_labels[x]
      out[is.na(out)] <- paste0("'", x[is.na(out)], "'")
      parse(text = out)
    }) +
    guides(
      colour = "none",
      shape = guide_legend(override.aes = list(
        shape = 22, size = 7, fill = c("white", "black"), colour = c("grey45", "black")
      ))
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 18),
      strip.text = element_text(size = 16),
      panel.spacing = grid::unit(1.5, "lines")
    ) +
    labs(
      x = "Absolute standardised mean difference",
      y = NULL, shape = NULL, colour = NULL
    )
}


# helper: overall benefit calibration of the drug classes analysed, with the
# predictions with the given prefix. Returns the matched pairs (and matchit
# objects, for the love plot), the mean observed benefit of following the model,
# the calibration slope / intercept and the figures.
overall_benefit <- function(prefix, step_title) {

  # match concordant vs discordant patients (returns the concordant, discordant
  # and combined matched populations, plus the matchit objects for love plots).
  # The matching can be changed here, e.g. if a love plot shows poor balance.
  matched <- cateval::match_benefit_pairs(
    data = analysis_dataset,
    drug_var = "drugclass",
    outcome_var = outcome_var,
    pred_cols = paste0(prefix, drugs),
    matching_var = if (prefix == "pred.") matching_vars else matching_vars_prog,
    caliper = c(prehba1c = 0.4),
    seed = 19840503
  )

  # group the matched pairs by predicted benefit (each cal_groups setting), for each population
  deciles_conc <- cateval::benefit_by_group(matched$concordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
  deciles_disc <- cateval::benefit_by_group(matched$discordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
  deciles_both <- cateval::benefit_by_group(matched$combined,   cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

  # calibration split by drug class (one panel per drug class analysed)
  by_drug <- cateval::benefit_by_drug(matched$combined, drugs = drugs, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

  title_end <- paste0(" (", tolower(step_title), ")")

  list(
    matched     = matched,
    # mean observed benefit of following the model, and calibration slope / intercept
    overall     = cateval::benefit_overall(matched$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed),
    calibration = cateval::benefit_citl_slope(matched$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed),
    # calibration plot for each matched population, side by side
    plot = overall_figure(deciles_conc, deciles_disc, deciles_both,
                          title = paste0("Overall Benefit Calibration", title_end)),
    plot_by_drug = plot_calibration_by_drug(by_drug, title = paste0("Overall Benefit Calibration by drug class", title_end))
  )
}


## Predicted vs Observed Calibration ----

## Predicted vs observed HbA1c (the outcome HbA1c: at 12 months, or at 6 months
## when the 12-month value is missing), with the prediction for the drug class
## initiated, for all initiations. The four figures have the same axes, so the
## versions of the predictions can be compared.
##
## Interpretation of plots:
# - Points close to the green line (y=x) indicate good calibration
# - Points above the green line indicate underprediction (observed HbA1c greater than predicted)
# - Points below the green line indicate overprediction (observed HbA1c less than predicted)

# helper: the initiations with an outcome, with the prediction (prefix) for the
# drug class initiated (pred_drug)
pred_obs_data_of <- function(prefix) {
  pred_obs_data <- analysis_dataset %>% filter(!is.na(posthba1cfinal))
  pred_obs_data$pred_drug <- pred_received(pred_obs_data, prefix)
  pred_obs_data
}

# the same axes in the four figures: the ranges in the four versions
pred_obs_axis_limits <- pred_obs_limits(lapply(c("pred.", "pred.recal_minimal.", "pred.recal_large.", "pred.recal."),
                                               pred_obs_data_of))

# helper: predicted vs observed figure, with the predictions with the given prefix
plot_pred_obs <- function(prefix, step_title) {
  draw_pred_obs(pred_obs_data_of(prefix), drugs, "Predicted vs Observed", step_title,
                limits = pred_obs_axis_limits)
}


## Predicted optimal therapy ----

## Clinical Decision Analysis: Does Recalibration Change Treatment Recommendations?
##
## This section compares the predicted optimal drug (rank 1 best predicted benefit) across
## different prediction approaches (original, minimal intercept recalibrated, intercept recalibrated,
## and intercept+slope recalibrated).
## Changes in optimal drug recommendations between models indicate that recalibration alters
## clinical decision-making for individual patients. This helps identify whether calibration
## adjustments have meaningful clinical implications beyond just improving prediction accuracy.
##
## Focus: Full cohort treatment recommendations

# The best drug class is chosen in two ways:
#   - all: from all drug classes the model predicts (model_drugs), as the model
#     would recommend them
#   - available: from the drug classes analysed (drugs) only
#   (only when some drug classes the model predicts are not analysed; otherwise
#   it is the same as all and is not made)

# helper: distribution of the predicted best drug class (rank 1), using the
# predictions with the given prefix for the drug classes in drug_set
plot_rank1 <- function(prefix, drug_set, step_title, set_title) {
  rank1_breakdown <- analysis_dataset %>%
    cateval::get_best_drugs(
      rank = 1,
      column_names = paste0(prefix, drug_set),
      final_var_name = prefix
    ) %>%
    mutate(drug = .data[[paste0(prefix, "rank1_drug_name")]]) %>%
    group_by(drug) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      total_n = sum(n),
      proportion = n / total_n * 100
    ) %>%
    arrange(desc(n))

  draw_rank1(rank1_breakdown, paste0(step_title, "; full cohort\n", str_to_sentence(set_title)), drug_set)
}

# helper: drug classes within 3 mmol/mol of the best, for each patient
tolerance3_sets <- function(prefix, drug_set) {
  analysis_dataset %>%
    cateval::get_best_drugs(
      tolerance = 3,
      column_names = paste0(prefix, drug_set),
      final_var_name = prefix
    ) %>%
    pull(paste0(prefix, "within_3_of_best_drug_name"))
}

# helper: drug classes added to and removed from each patient's tolerance-3 set
# by recalibration (original predictions vs the predictions with the given prefix)
plot_tolerance3_added_removed <- function(prefix, drug_set, step_title, set_title) {
  split_drugs <- function(x) str_split(ifelse(is.na(x), "", x), ",") %>% map(str_trim)
  original_list <- split_drugs(tolerance3_sets("pred.", drug_set))
  step_list     <- split_drugs(tolerance3_sets(prefix, drug_set))

  compare <- tibble(
    added   = map2_chr(step_list, original_list, ~ paste(setdiff(.x, .y), collapse = ",")),
    removed = map2_chr(original_list, step_list, ~ paste(setdiff(.x, .y), collapse = ","))
  ) %>%
    pivot_longer(everything(), names_to = "change_type", values_to = "drug") %>%
    separate_rows(drug, sep = ",") %>%
    mutate(drug = str_trim(drug)) %>%
    filter(drug != "") %>%
    count(change_type, drug) %>%
    complete(change_type = c("added", "removed"), drug, fill = list(n = 0))

  draw_tolerance3(compare, nrow(analysis_dataset),
                  paste0("Original vs ", step_title, "\n", str_to_sentence(set_title)), drug_set)
}

# available is only made when it differs from all
make_available <- !setequal(drugs, model_drugs)


## Best drug proportions and concordance ----

## For each set of predictions: how often each drug class is the predicted best
## (rank 1), and how often the drug initiated was that drug (concordance), with the
## best drug class chosen from all drug classes and from the drug classes analysed

# helper: one row per drug class (plus "All") for the predictions with the given
# prefix, with the best drug class chosen from the drug classes in drug_set
summarise_concordance <- function(prefix, drug_set, step_title, best_drug_from) {
  d <- analysis_dataset %>%
    cateval::get_best_drugs(
      rank = 1,
      column_names = paste0(prefix, drug_set),
      final_var_name = prefix
    ) %>%
    transmute(
      recommended_drug = factor(.data[[paste0(prefix, "rank1_drug_name")]], levels = drug_set),
      received_drug    = as.character(drugclass),
      concordant       = as.integer(as.character(recommended_drug) == received_drug)
    )

  total_n <- nrow(d)

  # .drop = FALSE so a drug the model never picks still gets a zero row
  by_drug <- d %>%
    group_by(recommended_drug, .drop = FALSE) %>%
    summarise(
      n_recommended = dplyr::n(),
      n_concordant  = sum(concordant),
      .groups = "drop"
    ) %>%
    mutate(recommended_drug = as.character(recommended_drug))

  bind_rows(
    by_drug,
    tibble(
      recommended_drug = "All",
      n_recommended    = total_n,
      n_concordant     = sum(d$concordant)
    )
  ) %>%
    mutate(
      model = step_title,
      best_drug_from = best_drug_from,
      total_n = total_n,
      # share of the cohort this drug is recommended for
      pct_recommended = n_recommended / total_n * 100,
      # of those recommended this drug, the share who received it
      pct_concordant_of_recommended = ifelse(
        n_recommended > 0, n_concordant / n_recommended * 100, NA_real_
      ),
      # share of the whole cohort both recommended and given this drug; the
      # "All" row of this column is the overall concordance rate
      pct_concordant_of_total = n_concordant / total_n * 100
    ) %>%
    select(
      model, best_drug_from, recommended_drug,
      n_recommended, total_n, pct_recommended,
      n_concordant, pct_concordant_of_recommended, pct_concordant_of_total
    )
}


## Validation of each version of the predictions ----

# helper: save a balance (love) plot, both matching directions in one figure
save_loveplot <- function(match_obj, file) {
  if (!is.null(match_obj)) save_plot(plot_loveplot(match_obj), file, panel_width = 3.4,
                                     panel_height = function(p) 0.3 * n_axis_levels(p, "y"))
}

# helper: the whole validation of one version of the predictions (prefix, e.g.
# "pred.recal."), with every figure saved as it is made. The files start with
# file_start (e.g. "1_original") and are numbered by figure, so they sort in this
# order:
#   01_pairwise, 02_overall, 03_overall_by_drug, 04_loveplot, 05_pred_obs
#   06_rank1: the predicted best drug class, of all drug classes (_all) and of the
#     drug classes analysed (_available, only when some are not analysed)
#   07_tolerance3: the same, for the drug classes within 3 mmol/mol of the best
#     (compared with the original predictions, so not made for them)
# Returns the pairwise results and overall benefit calibration, for the summaries.
validate_predictions <- function(prefix, step_title, file_start) {
  file <- function(number, figure) paste0(output_path, "/", file_start, "_", number, "_", figure, ".pdf")

  pairwise_results <- pairwise_results_of(prefix)
  save_pairwise(plot_pairwise(pairwise_results, "Pairwise Drug Calibration", step_title), file("01", "pairwise"))

  benefit <- overall_benefit(prefix, step_title)
  save_overall(benefit$plot, file("02", "overall"))
  save_overall_by_drug(benefit$plot_by_drug, file("03", "overall_by_drug"))
  save_loveplot(benefit$matched, file("04", "loveplot"))

  save_pred_obs(plot_pred_obs(prefix, step_title), file("05", "pred_obs"))

  save_decision(plot_rank1(prefix, model_drugs, step_title, "best of all drug classes"), file("06", "rank1_all"))
  if (make_available) {
    save_decision(plot_rank1(prefix, drugs, step_title, "best of the drug classes analysed"), file("06", "rank1_available"))
  }
  if (prefix != "pred.") {
    save_decision(plot_tolerance3_added_removed(prefix, model_drugs, step_title, "best of all drug classes"),
                  file("07", "tolerance3_all"))
    if (make_available) {
      save_decision(plot_tolerance3_added_removed(prefix, drugs, step_title, "best of the drug classes analysed"),
                    file("07", "tolerance3_available"))
    }
  }

  list(pairwise = pairwise_results, benefit = benefit)
}


### Part 5.1: Original predictions ----
# Does the model's predicted benefit calibrate as it is?

original <- validate_predictions("pred.", "Original predictions", "1_original")


### Part 5.2: Minimal recalibration-in-the-large ----
# Does it calibrate with the drug-specific intercepts beyond minimal_threshold?

recal_minimal <- validate_predictions("pred.recal_minimal.", "Minimal recalibration-in-the-large",
                                      "2_minimal_recalibration_in_the_large")


### Part 5.3: Recalibration-in-the-large ----
# Does it calibrate with a drug-specific intercept correction?

recal_large <- validate_predictions("pred.recal_large.", "Recalibration-in-the-large",
                                    "3_recalibration_in_the_large")


### Part 5.4: Recalibration ----
# Does it calibrate with a drug-specific intercept and slope correction?

recal <- validate_predictions("pred.recal.", "Recalibration", "4_recalibration")


### Part 5.5: Summaries ----
# Saved as 5_summary_...

#### Overall benefit estimate across the versions of the predictions

benefit_estimates <- bind_rows(
  original$benefit$overall %>% mutate(model = "Original"),
  recal_minimal$benefit$overall %>% mutate(model = "Minimal recalibration-in-the-large"),
  recal_large$benefit$overall %>% mutate(model = "Recalibration-in-the-large"),
  recal$benefit$overall %>% mutate(model = "Recalibration")
) %>%
  mutate(model = factor(
    model,
    levels = c("Original", "Minimal recalibration-in-the-large", "Recalibration-in-the-large", "Recalibration")
  ))

# (one analysis: the drug classes analysed, drawn as in v2)
save_benefit_estimates(
  draw_benefit_estimates(
    benefit_estimates %>% mutate(analysis_label = "drug classes analysed", key_result = TRUE)
  ),
  paste0(output_path, "/5_summary_overall_benefit_estimates.pdf")
)

# table of the average overall benefit per version of the predictions (estimate and 95% CI)
write.csv(
  benefit_estimates %>%
    select(model, estimate = observed, lower_ci = lower.ci, upper_ci = upper.ci, n_pairs = n),
  paste0(output_path, "/5_summary_overall_benefit_estimates.csv"),
  row.names = FALSE
)


#### Best drug proportions and concordance

concordance_summary <- bind_rows(
  summarise_concordance("pred.", model_drugs, "Original", "all"),
  if (make_available) summarise_concordance("pred.", drugs, "Original", "available"),
  summarise_concordance("pred.recal_minimal.", model_drugs, "Minimal recalibration-in-the-large", "all"),
  if (make_available) summarise_concordance("pred.recal_minimal.", drugs, "Minimal recalibration-in-the-large", "available"),
  summarise_concordance("pred.recal_large.", model_drugs, "Recalibration-in-the-large", "all"),
  if (make_available) summarise_concordance("pred.recal_large.", drugs, "Recalibration-in-the-large", "available"),
  summarise_concordance("pred.recal.", model_drugs, "Recalibration", "all"),
  if (make_available) summarise_concordance("pred.recal.", drugs, "Recalibration", "available")
)

print(concordance_summary)

write.csv(
  concordance_summary,
  paste0(output_path, "/5_summary_concordance_by_drug.csv"),
  row.names = FALSE
)

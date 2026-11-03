## Test run of example_validation_script_v1.R on CPRD: all initiations before
##  2020-10-14 (the cohort of the update scripts starts on that date). The
##  script is v1 unchanged apart from lines marked "CPRD test run": the CPRD
##  cohort is loaded in Step 1a, predictions come from Model_files/ and outputs
##  go to Scripts/00.V1_Validation/Outputs/Pre2020.
##
## This script demonstrates how to validate the 5-drug model in a new cohort
##  and how to visualise calibration results
##
## Model version: v1, the published model (Dennis et al. Lancet 2025) for the
##  5 drug classes.
##
## Model predictions come from the online model calculator or from the model
##  files (use_calculator, Step 0).
#
## There are several stages in this script:
# 0) Settings
# 1) Prepare your dataset with baseline features and outcomes
#    a) Formatting continuous and categorical variables
#    b) Creating a characteristics table
#    c) Imputing missing data if needed, and choosing the drug classes analysed
#    d) Writing the file to upload to the online calculator (calculator only)
# 2) Obtain model predictions for each drug class
#    a) Original predictions (online calculator or model files), made in each
#       imputed dataset and averaged
#    b) Recalibration-in-the-large
#    c) Recalibration
# 3) Perform calibration and clinical-decision analyses:
#    a) Pairwise drug calibration
#    b) Overall benefit calibration
#    c) Predicted vs observed calibration
#    d) Predicted optimal therapy (does recalibration change recommendations?)
# 4) Save plots


#############################################################################

# load required libraries
library(tidyverse)
library(tableone)
# library(remotes)
# remotes::install_github("PM-Cardoso/cateval")
library(cateval) # Custom R package developed by Pedro for calibration evaluation


#############################################################################

# Step 0: Settings ------------------------------------------------------------

# How the model predictions are made (Step 2):
#   TRUE:  with the online model calculator. Run the script twice: the first run
#          writes the file to upload (calculator_upload_file) and stops. Upload it in
#          the calculator's Validation Framework tab, save the file it returns as
#          calculator_results_file, and run the script again.
#   FALSE: from the model files, which need to be in the working directory.
use_calculator <- FALSE # CPRD test run: model files (Model_files/), CPRD data stays on the server

# where the outputs are saved
output_path <- "Scripts/00.V1_Validation/Outputs/Pre2020" # CPRD test run
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# files for the online calculator (only used when use_calculator is TRUE)
calculator_upload_file  <- paste0(output_path, "/calculator_upload.csv") # written in Step 1d
calculator_results_file <- "calculator_predictions.csv"                   # returned by the calculator

# The drug classes the model predicts (do not change these)
model_drugs <- c("SGLT2", "GLP1", "DPP4", "SU", "TZD")

# The drug classes to analyse. You do not need initiations of every drug class: the
# drug classes analysed are those of drugs_selected with initiations in your data
# (drugs, set in Step 1c). To leave out a drug class, e.g. one with too few
# initiations for stable estimates, remove it here, e.g.
#   drugs_selected <- c("SGLT2", "GLP1", "DPP4", "SU")
# Its initiations are dropped after imputation (Step 1c), so they still help impute
# missing values. The model still predicts it (original predictions, not
# recalibrated), so it can be the best drug class in the "all" outputs (Step 3).
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

# Step 1a: Prepare your dataset ------------------------------------------------
# Note: This portion of code can be adapted based on your data source and structure

# REQUIRED VARIABLES
# Your dataset (analysis_dataset) must include the following variables for validation:
#
# OUTCOME VARIABLES:
#   - posthba1cfinal: 12-month HbA1c outcome (numeric, mmol/mol)
#   - hba1cmonth: Months from treatment initiation to outcome measurement (numeric, 3-15 months)
#
# TREATMENT VARIABLES:
#   - drugclass: Drug class being initiated (categorical). Code the drug classes in
#     your data with these names (you do not need all of them):
#       "SGLT2", "GLP1", "DPP4", "SU", or "TZD"
#     Initiations with any other value are dropped.
#
# BASELINE CHARACTERISTICS:
#   - agetx: Age at treatment initiation (numeric, 18-79 years)
#   - sex: Patient sex (categorical: "Male" or "Female")
#   - t2dmduration: Duration of type 2 diabetes (numeric, 0-71 years)
#   - prehba1c: Baseline HbA1c (numeric, 53-110 mmol/mol)
#   - prebmi: Baseline BMI (numeric, 15-100 kg/m2)
#   - preegfr: Baseline eGFR (numeric, 45-300 mL/min per 1.73 m2)
#   - prealt: Baseline alanine aminotransferase (numeric, 0-200 U/L)
#   - pretotalcholesterol: Baseline total cholesterol (numeric, 0.5-20 mmol/L)
#   - prehdl: Baseline HDL cholesterol (numeric, 0.2-10 mmol/L)
#   - ethnicity: Patient ethnicity (categorical)
#     Valid values: "White", "South Asian", "Black", "Mixed", "Other", or "Missing"
#   - smoke: Smoking status (categorical)
#     Valid values: "Non-smoker", "Active smoker", "Ex-smoker", or "Not recorded"
#   - ncurrtx: Number of current antihyperglycaemic drugs (categorical: "1", "2", "3", or "4+")
#   - drugline: Number of previous drug classes tried (categorical: "2", "3", "4", or "5+")

# The ranges above are the ones the online calculator accepts. Initiations outside
# them are dropped after imputation (Step 1c), whichever way the predictions are made.
valid_ranges <- list(
  agetx = c(18, 79), t2dmduration = c(0, 71), prehba1c = c(53, 110),
  prebmi = c(15, 100), preegfr = c(45, 300), prealt = c(0, 200),
  pretotalcholesterol = c(0.5, 20), prehdl = c(0.2, 10), hba1cmonth = c(3, 15)
)


## CPRD test run: load the CPRD cohort (start) ---------------------------------
## Same as the cohort loading in 01.Temporal_Update/02.calibration_new_data_post2020.R,
## but initiations before 2020-10-14 and GLP1 kept as one class (as in v1).

# load required libraries
library(aurum)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

# Initiations before this date (first initiation date of the update scripts' cohort)
cohort_end <- as.Date("2020-10-14")

# Semaglutide substances (GLP1 drug_class in the cohort), left out as in the
# update scripts: semaglutide is not part of the published model's GLP1 class
sema_substances <- c("Low-dose semaglutide", "High-dose semaglutide",
                     "Oral semaglutide", "Semaglutide, dose unclear")

# set up dataset (cohort cached by 01.cohort_scripts.R)
analysis_raw <- analysis_raw %>%
  analysis$cached("ethnicity") %>%
  filter(dstartdate < cohort_end) %>%
  # select variables needed (to shorten the load time)
  select(
    patid, dstartdate,
    gender, dstartdate_age,
    t2dmduration, ethnicity_5cat,
    drug_substance, drug_class,
    smoking_cat,
    prebmi, prehba1c, preegfr, pretotalcholesterol, prehdl, prealt,
    drugline_all, MFN, SGLT2, GLP1, DPP4, TZD, SU, hba1cmonth,
    posthba1cfinal
  ) %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

analysis_dataset <- analysis_raw %>%
  mutate(
    pated = paste(patid, dstartdate, drug_class, sep = "."),
    agetx = dstartdate_age,

    # Sex
    sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),

    # Ethnicity
    ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
    ethnicity = factor(ethnicity, levels = c(0:5),
                       labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),

    # Smoking
    smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),

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
    "drugclass", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  ))) %>%
  as.data.frame()

cat("cohort: initiations from",
    format(min(analysis_dataset$dstartdate, na.rm = TRUE)), "to",
    format(max(analysis_dataset$dstartdate, na.rm = TRUE)), "\n")
table(analysis_dataset$drugclass, useNA = "ifany")

# v1 imputes with every column of analysis_dataset (Step 1c), so drop the date
analysis_dataset <- analysis_dataset %>% select(-dstartdate)

## CPRD test run: load the CPRD cohort (end) -----------------------------------


# initiations with a drugclass value that is not one of model_drugs (dropped below)
analysis_dataset %>%
  filter(!(as.character(drugclass) %in% model_drugs)) %>%
  count(drugclass, name = "n_dropped") %>%
  print()


# Caterogical variables must be properly formatted as factors with the specified levels.
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


# Step 1b: Create a characteristics table ---------------------------------

# This table describes the validation cohort at treatment initiation, stratified by
# initiated drug class.
# It includes baseline clinical characteristics, missingness indicators, outcome timing,
# observed 12-month HbA1c, and observed response (posthba1cfinal - prehba1c).


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

write.csv(
  characteristics_table_df,
  paste0(output_path, "/stage_1b_characteristics_by_drugclass.csv"),
  row.names = TRUE
)


# Step 1c: Prepare your dataset ------------------------------------------------

# Imputation of missing data can be performed using multiple imputation methods
# The following code generates m=10 imputed datasets. The model predictions are made in
# each imputed dataset and then averaged (Step 2). The imputed values are also averaged,
# to create a single complete dataset with mean imputed values for originally missing
# data, which gives the baseline characteristics used in the analyses (Step 3)

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
# predictions (Step 2). These have the same rows, in the same order, as analysis_dataset.
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

# Initiations outside the accepted ranges (valid_ranges, Step 1a) in any of the
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

# Initiations of drug classes not in drugs_selected (Step 0) are dropped
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
write.csv(drug_classes_analysed, paste0(output_path, "/stage_1c_drug_classes_analysed.csv"),
          row.names = FALSE)

if (length(drugs) < 2) {
  stop("At least two drug classes with initiations are needed for the validation; ",
       "found: ", paste(drugs, collapse = ", "), ".")
}


# Step 1d: File for the online calculator ------------------------------------

# Only when the predictions come from the online calculator (use_calculator, Step 0).
# The file holds every imputed dataset, one after the other (imputation says which),
# so it has mice_m rows per initiation, with numeric values rounded to 3 decimals. The calculator returns every column of this
# file in the same order, with the predictions added, so row_id links them back to
# these rows in Step 2.
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

# Step 2: Model predictions ------------------------------------------------

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

  ## Predictions from the model file, which needs to be in the working directory:
  ##   - fivedrugmodel_5knot_share_20230823.Rdata: the 5-drug model (m1.5.final)
  ## These give the same predictions as the online calculator.

  # load required libraries
  library(rms)

  # load model
  load("Model_files/fivedrugmodel_5knot_share_20230823.Rdata") # m1.5.final (CPRD test run: path)

  # The model also uses deprivation (imd5), which is set to the middle quintile as
  # in the online calculator. hba1cmonth is each initiation's own.
  # (the model needs a data.frame, not a tibble)
  prediction_list <- lapply(imputed_list, function(imputed) {
    as.data.frame(imputed) %>%
      mutate(imd5 = factor("3", levels = c("1 (least)", "2", "3", "4", "5 (most)")))
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
## and intercept + slope adjustments.

### Recalibration
# drugs (Step 1c) and model_drugs (Step 0)
outcome_var <- "posthba1cfinal"
pred_prefix <- "pred."

# Optional prognostic factors (leave empty if none). These are variables that
# affect the outcome irrespective of treatment, so they enter as a single overall
# (treatment-independent) effect: a shared main effect in the recalibration models
# below, plus covariate adjustment in the pairwise analysis and matching in the
# overall-benefit analysis (via matching_vars_prog).
prognostic_vars <- c()
prognostic_terms <- if (length(prognostic_vars) > 0) {
  paste0(" + ", paste(prognostic_vars, collapse = " + "))
} else {
  ""
}

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

#### Recalibration (drug-specific intercept & slope, shared prognostic effect)
formula_recal <- as.formula(
  paste0(outcome_var, " ~ 0 + drugclass + drugclass:pred_assigned", prognostic_terms)
)
fit_recal <- lm(formula_recal, data = recal_data)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal.", drug)]] <- predict_for_drug(fit_recal, drug)
}
recalib_recal_params <- summarize_lm(fit_recal, "recal")

#### Original predictions adjusted for prognostic factors
# Computed last, so the recalibrations above use the raw predictions as their base.
# The model's own prediction is kept (offset, coefficient 1) and only a shared
# coefficient per prognostic factor is added, shifting every drug equally: this
# moves the absolute predictions (and the predicted-vs-observed calibration) but
# leaves every between-drug benefit unchanged. The original pred.<drug> columns are
# overwritten so the rest of the script uses them directly. With prognostic_vars
# empty the predictions are returned unchanged.
formula_original <- as.formula(
  paste0(outcome_var, " ~ 0", prognostic_terms, " + offset(pred_assigned)")
)
fit_original <- lm(formula_original, data = recal_data)
for (drug in model_drugs) {
  analysis_dataset[[paste0(pred_prefix, drug)]] <- predict_for_drug(fit_original, drug)
}
recalib_original_params <- summarize_lm(fit_original, "original")

# drug classes not analysed (not in drugs) have no data to be recalibrated, so they
# keep their original predictions
for (drug in setdiff(model_drugs, drugs)) {
  analysis_dataset[[paste0("pred.recal_large.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
  analysis_dataset[[paste0("pred.recal.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
}

# Combine and save recalibration parameters
recalib_params_full <- bind_rows(recalib_original_params, recalib_recal_large_params, recalib_recal_params)

write_csv(recalib_params_full, paste0(output_path, "/recalibration_params_full.csv"))


#############################################################################

# ---- Saving figures ----
# Every figure is saved at a size worked out from its contents: each panel gets
# a fixed size in inches, and the file is exactly as large as the panels plus
# their axes, strips, titles and legend. A figure therefore comes out the same
# on any computer, with no white space around it and no text cut off (text is
# measured with the pdf device's built-in font metrics, not the system fonts).

# facet columns: one row for up to 4 panels, two rows above that
facet_ncol <- function(n_panels) if (n_panels <= 4) n_panels else ceiling(n_panels / 2)

# number of categories on a discrete axis ("x" or "y") of a plot
n_axis_levels <- function(plot, axis = "x") {
  length(ggplot2::ggplot_build(plot)$layout$panel_params[[1]][[axis]]$get_labels())
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

# helper: draw the pairwise calibration panels (one per pair of drug classes,
# labelled by `title`), each square and on one symmetric range
pairwise_plot <- function(plot_data, plot_title, plot_subtitle = NULL) {
  # a trailing "(step)" in the title goes in the subtitle, so the title stays on one line
  step <- regmatches(plot_title, regexec(" [(]([^()]*)[)]$", plot_title))[[1]]
  if (is.null(plot_subtitle) && length(step) == 2) {
    plot_title <- sub(" [(][^()]*[)]$", "", plot_title)
    plot_subtitle <- paste0(toupper(substr(step[2], 1, 1)), substring(step[2], 2))
  }
  n_panels <- n_distinct(plot_data$title)
  ggplot(plot_data, aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
    geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
    geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
    geom_point(colour = "#800020") +
    geom_errorbar(colour = "#800020") +
    facet_wrap(~title, ncol = facet_ncol(n_panels)) +
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
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 15),
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 17),
      plot.title.position = "plot",
      panel.spacing.x = grid::unit(1.6, "lines"), panel.spacing.y = grid::unit(0.6, "lines")
    ) +
    labs(
      x = expression("Predicted HbA"[1*c]*" benefit (mmol/mol)"),
      y = expression("Observed HbA"[1*c]*" benefit (mmol/mol)"),
      # wrapped to about the width of the panels
      title = str_wrap(plot_title, max(30, round(22 * facet_ncol(n_panels)))),
      subtitle = plot_subtitle
    )
}

# side of the square pairwise panels, in inches: the same in every script, and
# wide enough for the longest drug pair name on one line
pairwise_panel_size <- 3.4


#############################################################################

# Step 3: Evaluate performance ----------------------------------------------

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
#   on top (as an overall, treatment-independent effect) in the prognostic-adjusted
#   and recalibrated versions.

adjustment_vars = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt",
                     "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth",
                     "sex", "smoke", "ncurrtx", "drugline")

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

### Original predictions
pairwise_results_original <- cateval::pairwise_calibration(
  data = analysis_dataset,
  drug_var = "drugclass",
  drugs = drugs,
  prediction_vars = paste0("pred.", drugs),
  outcome_var = outcome_var,
  cal_groups = cal_groups,
  adjustment_var = unique(c(adjustment_vars, prognostic_vars))
)

plot_pairwise_original <- plot_pairwise(pairwise_results_original,
  "Pairwise Drug Calibration", "Original predictions")

### Recalibration-in-the-large predictions
pairwise_results_recal_large <- cateval::pairwise_calibration(
  data = analysis_dataset,
  drug_var = "drugclass",
  drugs = drugs,
  prediction_vars = paste0("pred.recal_large.", drugs),
  outcome_var = outcome_var,
  cal_groups = cal_groups,
  adjustment_var = unique(c(adjustment_vars, prognostic_vars))
)

plot_pairwise_recal_large <- plot_pairwise(pairwise_results_recal_large,
  "Pairwise Drug Calibration", "Recalibration-in-the-large")

### Recalibration predictions
pairwise_results_recal <- cateval::pairwise_calibration(
  data = analysis_dataset,
  drug_var = "drugclass",
  drugs = drugs,
  prediction_vars = paste0("pred.recal.", drugs),
  outcome_var = outcome_var,
  cal_groups = cal_groups,
  adjustment_var = unique(c(adjustment_vars, prognostic_vars))
)

plot_pairwise_recal <- plot_pairwise(pairwise_results_recal,
  "Pairwise Drug Calibration", "Recalibration")



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
# The same covariates plus any prognostic factors in prognostic_vars. Used for
# matching in every overall-benefit step (the prognostic factor is a cross-cutting
# control, not a step of its own); equal to matching_vars until you populate
# prognostic_vars.
matching_vars_prog <- unique(c(matching_vars, prognostic_vars))

# Confidence intervals: FALSE uses fast analytical intervals; set to TRUE for
# pair-level bootstrap intervals (slower) that account for matching with
# replacement.
use_bootstrap <- FALSE
n_boot <- 1000
boot_seed <- 42

# ---- Local plotting helpers ----
# Points are the mean observed benefit per group (from benefit_by_group() /
# benefit_by_drug()), with a confidence interval, against the mean predicted
# benefit. The dashed blue curve is a LOESS fit through these group point
# estimates themselves (not the individual matched pairs), so the curve
# summarises exactly the points shown, not the noisier patient-level data.
# The `n` column from benefit_by_group() is the number of matched pairs in each
# group, so its sum is the total number of pairs reported in the subtitle.
plot_calibration <- function(deciles, title = NULL, min_n = 100) {
  # When several cal_groups settings are present (an n_groups column from
  # benefit_by_group()), show the largest n_groups whose groups all contain at
  # least `min_n` matched pairs; if none does, fall back to the smallest n_groups.
  groups_label <- ""
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
    groups_label <- paste0(", ", chosen, ifelse(chosen == 1, " group", " groups"))
  }

  n_pairs <- sum(deciles$n, na.rm = TRUE)
  ggplot(deciles, aes(x = pred.benefit, y = obs.benefit)) +
    geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    # LOESS only where there are more than 5 group points (in each drug class);
    # through 5 or fewer its confidence band is unstable and can swamp the axis
    geom_smooth(data = function(d) d %>% group_by(across(any_of("index_drug"))) %>% filter(n() > 5) %>% ungroup(),
                method = "loess", span = 1, se = TRUE, colour = "blue", linetype = "dashed", linewidth = 0.75) +
    geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, colour = "black") +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 17)
    ) +
    labs(
      x = expression("Predicted HbA"[1*c]*" benefit (mmol/mol)"),
      y = expression("Observed HbA"[1*c]*" benefit (mmol/mol)"),
      title = title,
      subtitle = paste0("n = ", format(n_pairs, big.mark = ","), " matched pairs", groups_label)
    )
}

# same plot, faceted by the drug class that anchors each pair (the index /
# matched patient's drug); the LOESS curve is fitted separately within each
# facet, from that drug's own group points. Each facet strip is labelled with
# the number of matched pairs (index-drug initiations) for that drug class,
# summing the per-group `n` returned by benefit_by_drug().
plot_calibration_by_drug <- function(by_drug, title = NULL, min_n = 100) {
  # When several cal_groups settings are present (an n_groups column from
  # benefit_by_drug()), choose PER DRUG the largest n_groups whose groups all
  # contain at least `min_n` matched pairs. If even the coarsest grouping has a
  # group below min_n, fall back to that drug's smallest n_groups (fewest groups,
  # hence the largest per-group n).
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
      select(-.min_group_n, -.chosen_groups)
  }

  # Facet strip: number of matched pairs and the number of groups actually shown
  # for that drug (these can now differ between drug classes).
  facet_counts <- by_drug %>%
    group_by(index_drug) %>%
    summarise(n_pairs = sum(n, na.rm = TRUE), n_shown = dplyr::n(), .groups = "drop")
  facet_labels <- setNames(
    paste0(recode(as.character(facet_counts$index_drug), !!!drug_labels),
           "\n(n = ", format(facet_counts$n_pairs, big.mark = ",", trim = TRUE),
           ", ", facet_counts$n_shown, ifelse(facet_counts$n_shown == 1, " group)", " groups)")),
    as.character(facet_counts$index_drug)
  )
  # n_groups is dropped so plot_calibration() does not re-select across drugs
  plot_calibration(by_drug %>% select(-n_groups), title = title) +
    facet_wrap(~ index_drug, ncol = facet_ncol(nrow(facet_counts)), labeller = labeller(index_drug = facet_labels)) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 16),
      # room between panels so neighbouring x axis labels do not touch
      panel.spacing.x = grid::unit(2.2, "lines")
    )
}

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

# Balance (love) plot for one matching step: absolute standardised mean
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


### Original predictions ----

# match concordant vs discordant patients (returns the concordant, discordant
# and combined matched populations, plus the matchit objects for love plots).
# The matching can be changed for each step separately, e.g. if its love plot
# shows poor balance.
match_original <- cateval::match_benefit_pairs(
  data = analysis_dataset,
  drug_var = "drugclass",
  outcome_var = outcome_var,
  pred_cols = paste0("pred.", drugs),
  matching_var = matching_vars_prog,
  caliper = c(prehba1c = 0.4),
  seed = 19840503
)

# group the matched pairs by predicted benefit (each cal_groups setting), for each population
deciles_original_conc <- cateval::benefit_by_group(match_original$concordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_original_disc <- cateval::benefit_by_group(match_original$discordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_original_both <- cateval::benefit_by_group(match_original$combined,   cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

# mean observed benefit of following the model, and calibration slope / intercept
overall_original     <- cateval::benefit_overall(match_original$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
calibration_original <- cateval::benefit_citl_slope(match_original$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

# calibration plot for each matched population, then the three side by side
plot_original_conc <- plot_calibration(deciles_original_conc, title = "Concordant matched to discordant")
plot_original_disc <- plot_calibration(deciles_original_disc, title = "Discordant matched to concordant")
plot_original_both <- plot_calibration(deciles_original_both, title = "Combined")

plot_overall_original <- plot_original_conc + plot_original_disc + plot_original_both +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(title = "Overall Benefit Calibration (original predictions)", theme = theme(plot.title = element_text(size = 20)))

# calibration split by drug class (one panel per drug class analysed)
by_drug_original <- cateval::benefit_by_drug(match_original$combined, drugs = drugs, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
plot_overall_original_bydrug <- plot_calibration_by_drug(by_drug_original, title = "Overall Benefit Calibration by drug class (original predictions)")


### Recalibration-in-the-large predictions ----

match_recal_large <- cateval::match_benefit_pairs(
  data = analysis_dataset,
  drug_var = "drugclass",
  outcome_var = outcome_var,
  pred_cols = paste0("pred.recal_large.", drugs),
  matching_var = matching_vars_prog,
  caliper = c(prehba1c = 0.4),
  seed = 19840503
)

deciles_recal_large_conc <- cateval::benefit_by_group(match_recal_large$concordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_large_disc <- cateval::benefit_by_group(match_recal_large$discordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_large_both <- cateval::benefit_by_group(match_recal_large$combined,   cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

overall_recal_large     <- cateval::benefit_overall(match_recal_large$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
calibration_recal_large <- cateval::benefit_citl_slope(match_recal_large$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

plot_recal_large_conc <- plot_calibration(deciles_recal_large_conc, title = "Concordant matched to discordant")
plot_recal_large_disc <- plot_calibration(deciles_recal_large_disc, title = "Discordant matched to concordant")
plot_recal_large_both <- plot_calibration(deciles_recal_large_both, title = "Combined")

plot_overall_recal_large <- plot_recal_large_conc + plot_recal_large_disc + plot_recal_large_both +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(title = "Overall Benefit Calibration (recalibration-in-the-large)", theme = theme(plot.title = element_text(size = 20)))

by_drug_recal_large <- cateval::benefit_by_drug(match_recal_large$combined, drugs = drugs, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
plot_overall_recal_large_bydrug <- plot_calibration_by_drug(by_drug_recal_large, title = "Overall Benefit Calibration by drug class (recalibration-in-the-large)")


### Recalibration predictions ----

match_recal <- cateval::match_benefit_pairs(
  data = analysis_dataset,
  drug_var = "drugclass",
  outcome_var = outcome_var,
  pred_cols = paste0("pred.recal.", drugs),
  matching_var = matching_vars_prog,
  caliper = c(prehba1c = 0.4),
  seed = 19840503
)

deciles_recal_conc <- cateval::benefit_by_group(match_recal$concordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_disc <- cateval::benefit_by_group(match_recal$discordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_both <- cateval::benefit_by_group(match_recal$combined,   cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

overall_recal     <- cateval::benefit_overall(match_recal$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
calibration_recal <- cateval::benefit_citl_slope(match_recal$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

plot_recal_conc <- plot_calibration(deciles_recal_conc, title = "Concordant matched to discordant")
plot_recal_disc <- plot_calibration(deciles_recal_disc, title = "Discordant matched to concordant")
plot_recal_both <- plot_calibration(deciles_recal_both, title = "Combined")

plot_overall_recal <- plot_recal_conc + plot_recal_disc + plot_recal_both +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(title = "Overall Benefit Calibration (recalibration)", theme = theme(plot.title = element_text(size = 20)))

by_drug_recal <- cateval::benefit_by_drug(match_recal$combined, drugs = drugs, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
plot_overall_recal_bydrug <- plot_calibration_by_drug(by_drug_recal, title = "Overall Benefit Calibration by drug class (recalibration)")


### Compare the overall benefit estimate across model versions ----

benefit_estimates <- bind_rows(
  overall_original %>% mutate(model = "Original"),
  overall_recal_large %>% mutate(model = "Recalibration-in-the-large"),
  overall_recal %>% mutate(model = "Recalibration")
) %>%
  mutate(model = factor(
    model,
    levels = c("Original", "Recalibration-in-the-large", "Recalibration")
  ))

plot_benefit_estimates <- benefit_estimates %>%
  ggplot(aes(x = model, y = observed, ymin = lower.ci, ymax = upper.ci, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_pointrange(size = 0.7) +
  coord_flip() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none", axis.text = element_text(size = 16)) +
  labs(
    x = NULL,
    y = expression("Observed HbA"[1*c]*" benefit estimate (mmol/mol)"),
    title = "Overall benefit estimates"
  )


## Predicted vs Observed Calibration ----

## Interpretation of plots:
# - Points close to the green line (y=x) indicate good calibration
# - Points above the green line indicate underprediction (observed HbA1c greater than predicted)
# - Points below the green line indicate overprediction (observed HbA1c less than predicted)

### Original predictions
plot_pred_obs_original <- analysis_dataset %>%
  mutate(
    pred_drug = pred_received(analysis_dataset, "pred.")
  ) %>%
  select(pred_drug, posthba1cfinal) %>%
  ggplot(aes(x = pred_drug, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 17),
    # room for the last x axis label, which overhangs the panel
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  ) +
  labs(
    x = expression("Predicted 12-month HbA"[1*c]*" (mmol/mol)"),
    y = expression("Observed 12-month HbA"[1*c]*" (mmol/mol)"),
    title = "Predicted vs Observed",
    subtitle = paste0("Original predictions; n = ", format(sum(!is.na(analysis_dataset$posthba1cfinal)), big.mark = ","), " initiations")
  )

### Recalibration-in-the-large predictions
plot_pred_obs_recal_large <- analysis_dataset %>%
  mutate(
    pred_drug = pred_received(analysis_dataset, "pred.recal_large.")
  ) %>%
  select(pred_drug, posthba1cfinal) %>%
  ggplot(aes(x = pred_drug, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 17),
    # room for the last x axis label, which overhangs the panel
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  ) +
  labs(
    x = expression("Predicted 12-month HbA"[1*c]*" (mmol/mol)"),
    y = expression("Observed 12-month HbA"[1*c]*" (mmol/mol)"),
    title = "Predicted vs Observed",
    subtitle = paste0("Recalibration-in-the-large; n = ", format(sum(!is.na(analysis_dataset$posthba1cfinal)), big.mark = ","), " initiations")
  )

### Recalibration predictions
plot_pred_obs_recal <- analysis_dataset %>%
  mutate(
    pred_drug = pred_received(analysis_dataset, "pred.recal.")
  ) %>%
  select(pred_drug, posthba1cfinal) %>%
  ggplot(aes(x = pred_drug, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 17),
    # room for the last x axis label, which overhangs the panel
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  ) +
  labs(
    x = expression("Predicted 12-month HbA"[1*c]*" (mmol/mol)"),
    y = expression("Observed 12-month HbA"[1*c]*" (mmol/mol)"),
    title = "Predicted vs Observed",
    subtitle = paste0("Recalibration; n = ", format(sum(!is.na(analysis_dataset$posthba1cfinal)), big.mark = ","), " initiations")
  )


## Predicted optimal therapy ----

## Clinical Decision Analysis: Does Recalibration Change Treatment Recommendations?
##
## This section compares the predicted optimal drug (rank 1 best predicted benefit) across
## different prediction approaches (original, intercept recalibrated, and intercept+slope recalibrated).
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

  rank1_breakdown %>%
    mutate(drug = factor(drug, levels = model_drugs)) %>%
    ggplot(aes(x = reorder(drug, -n), y = proportion, fill = drug)) +
    geom_col(color = "black", size = 0.3) +
    geom_text(aes(label = paste0(format(n, big.mark = ","), "\n(", sprintf("%.1f", proportion), "%)")),
              vjust = -0.3, size = 5.5, lineheight = 0.9) +
    scale_fill_manual(
      values = drug_colours,
      guide = "none"
    ) +
    scale_x_discrete(labels = function(x) str_wrap(drug_labels[x], 12)) +
    scale_y_continuous(limits = c(0, max(rank1_breakdown$proportion) * 1.15),
                     breaks = function(lim) { b <- pretty(lim); b[b <= 100] }) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 17),
      plot.title.position = "plot",
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Drug Class",
      y = "Proportion of Cohort (%)",
      title = "Optimal Drug Distribution - Rank 1",
      subtitle = paste0(step_title, "; full cohort
", str_to_sentence(set_title))
    )
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
    complete(change_type = c("added", "removed"), drug, fill = list(n = 0)) %>%
    mutate(
      drug = factor(drug, levels = intersect(model_drugs, drug)),
      signed_n = ifelse(change_type == "removed", -n, n)
    )

  plot_title <- "Tolerance-3 drug changes"
  if (nrow(compare) == 0) {
    return(ggplot() +
             annotate("text", x = 0, y = 0, label = "No drug classes added or removed") +
             theme_void() +
             labs(title = plot_title, subtitle = paste0("Original vs ", step_title, "
", str_to_sentence(set_title))))
  }

  compare %>%
    ggplot(aes(x = drug, y = signed_n, fill = drug, alpha = change_type)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
    scale_fill_manual(
      values = drug_colours,
      guide = "none"
    ) +
    scale_alpha_manual(values = c("added" = 1, "removed" = 0.6)) +
    scale_x_discrete(labels = function(x) str_wrap(drug_labels[x], 12)) +
    scale_y_continuous(
      name = "Net count (n)",
      labels = scales::comma,
      sec.axis = sec_axis(
        ~ . * 100 / nrow(analysis_dataset),
        name = "Net percentage (%)",
        labels = function(x) sprintf("%.1f", x)
      )
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 17),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    ) +
    labs(
      x = "Drug Class",
      title = plot_title,
      subtitle = paste0("Original vs ", step_title, "
", str_to_sentence(set_title)),
      alpha = "Change"
    )
}

# available is only made when it differs from all
make_available <- !setequal(drugs, model_drugs)

### Original predictions
plot_rank1_original_all       <- plot_rank1("pred.", model_drugs, "Original predictions", "best of all drug classes")
plot_rank1_original_available <- if (make_available) plot_rank1("pred.", drugs, "Original predictions", "best of the drug classes analysed")

### Recalibration-in-the-large predictions
plot_rank1_recal_large_all       <- plot_rank1("pred.recal_large.", model_drugs, "Recalibration-in-the-large", "best of all drug classes")
plot_rank1_recal_large_available <- if (make_available) plot_rank1("pred.recal_large.", drugs, "Recalibration-in-the-large", "best of the drug classes analysed")

plot_tolerance3_recal_large_all       <- plot_tolerance3_added_removed("pred.recal_large.", model_drugs, "Recalibration-in-the-large", "best of all drug classes")
plot_tolerance3_recal_large_available <- if (make_available) plot_tolerance3_added_removed("pred.recal_large.", drugs, "Recalibration-in-the-large", "best of the drug classes analysed")

### Recalibration predictions
plot_rank1_recal_all       <- plot_rank1("pred.recal.", model_drugs, "Recalibration", "best of all drug classes")
plot_rank1_recal_available <- if (make_available) plot_rank1("pred.recal.", drugs, "Recalibration", "best of the drug classes analysed")

plot_tolerance3_recal_all       <- plot_tolerance3_added_removed("pred.recal.", model_drugs, "Recalibration", "best of all drug classes")
plot_tolerance3_recal_available <- if (make_available) plot_tolerance3_added_removed("pred.recal.", drugs, "Recalibration", "best of the drug classes analysed")


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

concordance_summary <- bind_rows(
  summarise_concordance("pred.", model_drugs, "Original", "all"),
  if (make_available) summarise_concordance("pred.", drugs, "Original", "available"),
  summarise_concordance("pred.recal_large.", model_drugs, "Recalibration-in-the-large", "all"),
  if (make_available) summarise_concordance("pred.recal_large.", drugs, "Recalibration-in-the-large", "available"),
  summarise_concordance("pred.recal.", model_drugs, "Recalibration", "all"),
  if (make_available) summarise_concordance("pred.recal.", drugs, "Recalibration", "available")
)

print(concordance_summary)

write.csv(
  concordance_summary,
  paste0(output_path, "/concordance_by_drug.csv"),
  row.names = FALSE
)


#############################################################################

# Step 4: Save Plots by Validation Step -----------------------------------

## This section saves plots for the validation procedure, which evaluates the
## differential (between-drug) treatment effect at three steps:
##
## STEP 1 — Original
##   Does the original model's differential effect calibrate as is?
##
## STEP 2 — Recalibration-in-the-large
##   Does a drug-specific intercept correction calibrate the differential effect?
##
## STEP 3 — Recalibration
##   Does a drug-specific intercept + slope correction calibrate it?
##
## The prognostic factor (prognostic_vars) is a cross-cutting control applied
## within every step when identified — adjusted for in the pairwise analysis and
## matched on in the overall-benefit analysis — not a step of its own.


## Save plots to separate PDFs by validation step ----

# helper: save a pairwise figure
save_pairwise <- function(fig, file) {
  save_plot(fig, file, panel_width = pairwise_panel_size, panel_height = pairwise_panel_size)
}


## STEP 1: Original ----

save_pairwise(plot_pairwise_original, paste0(output_path, "/step1_original_pairwise.pdf"))
save_plot(plot_overall_original, paste0(output_path, "/step1_original_overall.pdf"),
          panel_width = 4.2, panel_height = 4.2)
save_plot(plot_overall_original_bydrug, paste0(output_path, "/step1_original_overall_by_drug.pdf"),
          panel_width = 3.4, panel_height = 3.4)
save_plot(plot_pred_obs_original, paste0(output_path, "/step1_original_pred_obs.pdf"),
          panel_width = 5.5, panel_height = 4.2)

# balance (love) plot, both matching directions in one figure
save_plot(plot_loveplot(match_original), paste0(output_path, "/step1_original_loveplot.pdf"),
          panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))


## STEP 2: Recalibration-in-the-large ----

save_pairwise(plot_pairwise_recal_large, paste0(output_path, "/step2_recalibration_in_the_large_pairwise.pdf"))
save_plot(plot_overall_recal_large, paste0(output_path, "/step2_recalibration_in_the_large_overall.pdf"),
          panel_width = 4.2, panel_height = 4.2)
save_plot(plot_overall_recal_large_bydrug, paste0(output_path, "/step2_recalibration_in_the_large_overall_by_drug.pdf"),
          panel_width = 3.4, panel_height = 3.4)
save_plot(plot_pred_obs_recal_large, paste0(output_path, "/step2_recalibration_in_the_large_pred_obs.pdf"),
          panel_width = 5.5, panel_height = 4.2)

# balance (love) plot, both matching directions in one figure
save_plot(plot_loveplot(match_recal_large), paste0(output_path, "/step2_recalibration_in_the_large_loveplot.pdf"),
          panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))


## STEP 3: Recalibration ----

save_pairwise(plot_pairwise_recal, paste0(output_path, "/step3_recalibration_pairwise.pdf"))
save_plot(plot_overall_recal, paste0(output_path, "/step3_recalibration_overall.pdf"),
          panel_width = 4.2, panel_height = 4.2)
save_plot(plot_overall_recal_bydrug, paste0(output_path, "/step3_recalibration_overall_by_drug.pdf"),
          panel_width = 3.4, panel_height = 3.4)
save_plot(plot_pred_obs_recal, paste0(output_path, "/step3_recalibration_pred_obs.pdf"),
          panel_width = 5.5, panel_height = 4.2)

# balance (love) plot, both matching directions in one figure
save_plot(plot_loveplot(match_recal), paste0(output_path, "/step3_recalibration_loveplot.pdf"),
          panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))


## Additional Analyses: Clinical Decision Analysis ----

save_plot(plot_rank1_original_all, paste0(output_path, "/clinical_decision_rank1_original_all.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
if (make_available) save_plot(plot_rank1_original_available, paste0(output_path, "/clinical_decision_rank1_original_available.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)

save_plot(plot_rank1_recal_large_all, paste0(output_path, "/clinical_decision_rank1_recal_large_all.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
if (make_available) save_plot(plot_rank1_recal_large_available, paste0(output_path, "/clinical_decision_rank1_recal_large_available.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
save_plot(plot_tolerance3_recal_large_all, paste0(output_path, "/clinical_decision_tolerance3_added_removed_original_vs_recal_large_all.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
if (make_available) save_plot(plot_tolerance3_recal_large_available, paste0(output_path, "/clinical_decision_tolerance3_added_removed_original_vs_recal_large_available.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)

save_plot(plot_rank1_recal_all, paste0(output_path, "/clinical_decision_rank1_recal_all.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
if (make_available) save_plot(plot_rank1_recal_available, paste0(output_path, "/clinical_decision_rank1_recal_available.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
save_plot(plot_tolerance3_recal_all, paste0(output_path, "/clinical_decision_tolerance3_added_removed_original_vs_recal_all.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
if (make_available) save_plot(plot_tolerance3_recal_available, paste0(output_path, "/clinical_decision_tolerance3_added_removed_original_vs_recal_available.pdf"),
          panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)

save_plot(plot_benefit_estimates, paste0(output_path, "/overall_benefit_estimates_by_model_version.pdf"),
          panel_width = 5.5, panel_height = 2.4)

# table of the average overall benefit per step (estimate and 95% CI)
write.csv(
  benefit_estimates %>%
    select(model, estimate = observed, lower_ci = lower.ci, upper_ci = upper.ci, n_pairs = n),
  paste0(output_path, "/overall_benefit_estimates_by_model_version.csv"),
  row.names = FALSE
)

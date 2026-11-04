## This script demonstrates how to validate the 5-drug model in a new cohort
##  and how to visualise calibration results
##
## Model version: v2, the published model with its update. GLP1 is split into
##  two subtypes, dulaglutide and liraglutide (with exenatide), the main drug
##  classes are recalibrated, and oral semaglutide, injectable semaglutide and
##  tirzepatide are added. To validate the published model alone (v1), use
##  example_validation_script_v1.R instead.
##
## Model predictions come from the online model calculator or from the model
##  files (use_calculator, Step 0).
##
## The validation is done in analyses of a few drug classes each (Step 3):
##  - oral: the oral drug classes against each other (oral semaglutide, SGLT2i,
##    DPP4i, SU, TZD), with the predicted best drug among them
##  - one for each injectable drug class (dulaglutide, liraglutide, injectable
##    semaglutide, tirzepatide): that drug class against SGLT2i, DPP4i, SU and TZD
##  The GLP-1RA drug classes (the injectables and oral semaglutide) all start from
##  the published model's GLP1 prediction, so between them the predicted benefit is
##  (nearly) the same for every patient. They are compared with each other only as
##  one pooled estimate per pair.
##  An analysis is only done when your data has initiations of its drug classes
##  (a message says which analyses are left out). To leave out a drug class (e.g.
##  one with too few initiations), use drugs_selected (Step 0).
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
# 3) Perform calibration and clinical-decision analyses, for each set of
#    predictions (original, recalibration-in-the-large, recalibration):
#    a) Oral drug classes: pairwise calibration, overall benefit calibration,
#       predicted vs observed, predicted optimal therapy
#    b) Each injectable drug class against the oral drug classes: pairwise
#       calibration, overall benefit calibration
#    c) The injectable drug classes against each other (pooled pairwise estimates)
#    d) Predicted vs observed of every drug class validated (one panel each)
#    e) Summaries: overall benefit and concordance with the predicted best drug
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
# The online calculator does not yet return the GLP1 subtype predictions
# (pred.Dulaglutide, pred.Liraglutide), so use the model files for now.
use_calculator <- FALSE

# where the outputs are saved
output_path <- "Outputs/Validation_Procedure"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# files for the online calculator (only used when use_calculator is TRUE)
calculator_upload_file  <- paste0(output_path, "/calculator_upload.csv") # written in Step 1d
calculator_results_file <- "calculator_predictions.csv"                   # returned by the calculator

# Drug classes of this model version (do not change these):
#   oral_drugs: the oral drug classes, validated against each other (Step 3a).
#     Oral semaglutide is first, so it is drug 1 in the pairwise comparisons.
#   injectable_drugs: the injectable drug classes, each validated against the oral
#     drug classes (Step 3b). Dulaglutide is dulaglutide only; liraglutide includes
#     exenatide.
#   glp1_drugs: the GLP-1RA drug classes, whose predictions all start from the
#     published model's GLP1 prediction. Two of them differ by (nearly) the same
#     predicted benefit for every patient, so they are compared with each other only
#     as one pooled pairwise estimate, and an injectable drug class is compared with
#     comparator_drugs only in its overall benefit and best drug.
oral_drugs       <- c("Oral_Semaglutide", "SGLT2", "DPP4", "SU", "TZD")
injectable_drugs <- c("Dulaglutide", "Liraglutide", "Inje_Semaglutide", "TZP")
glp1_drugs       <- c("Oral_Semaglutide", injectable_drugs)
comparator_drugs <- setdiff(oral_drugs, glp1_drugs)  # SGLT2, DPP4, SU, TZD

# All drug classes the model predicts
model_drugs <- c(oral_drugs, injectable_drugs)

# The drug classes to analyse. You do not need initiations of every drug class: the
# drug classes analysed are those of drugs_selected with initiations in your data
# (drugs, set in Step 1c). To leave out a drug class, e.g. one with too few
# initiations for its own analysis or pairwise estimates, remove it here, e.g.
#   drugs_selected <- setdiff(model_drugs, "TZD")
# Its initiations are dropped after imputation (Step 1c), so they still help impute
# missing values.
drugs_selected <- model_drugs

if (!all(drugs_selected %in% model_drugs)) {
  stop("drugs_selected has drug classes the model does not predict: ",
       paste(setdiff(drugs_selected, model_drugs), collapse = ", "), call. = FALSE)
}

# drug class names as shown in the plots
drug_labels <- c("SGLT2" = "SGLT2i", "Dulaglutide" = "Dulaglutide", "Liraglutide" = "Liraglutide",
                 "DPP4" = "DPP4i", "SU" = "SU", "TZD" = "TZD",
                 "Oral_Semaglutide" = "Oral Semaglutide", "Inje_Semaglutide" = "S.c. Semaglutide",
                 "TZP" = "Tirzepatide")

# drug class colours (Okabe-Ito), as in the online calculator: the GLP1 subtypes and
# semaglutide take the GLP-1RA colour, and tirzepatide has its own
drug_colours <- c("SGLT2" = "#E69F00", "Dulaglutide" = "#56B4E9", "Liraglutide" = "#56B4E9",
                  "DPP4" = "#0072B2", "SU" = "#CC79A7", "TZD" = "#D55E00",
                  "Oral_Semaglutide" = "#56B4E9", "Inje_Semaglutide" = "#56B4E9",
                  "TZP" = "#009E73")


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
#       "SGLT2", "Dulaglutide", "Liraglutide" (liraglutide and exenatide), "DPP4", "SU",
#       "TZD", "Oral_Semaglutide", "Inje_Semaglutide", or "TZP" (tirzepatide)
#     GLP-1RA initiations other than semaglutide and tirzepatide are coded by
#     substance as "Dulaglutide" or "Liraglutide", not "GLP1".
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

# initiations per drug class analysed (also in recalibration_params_full.csv, Step 2)
print(analysis_dataset %>% count(drugclass, name = "n_initiations"))

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
#   - pred.Dulaglutide: Predicted 12-month HbA1c on dulaglutide (numeric, mmol/mol)
#   - pred.Liraglutide: Predicted 12-month HbA1c on liraglutide or exenatide (numeric, mmol/mol)
#   - pred.DPP4: Predicted 12-month HbA1c on DPP4 inhibitor (numeric, mmol/mol)
#   - pred.SU: Predicted 12-month HbA1c on Sulfonylurea (numeric, mmol/mol)
#   - pred.TZD: Predicted 12-month HbA1c on Thiazolidinedione (numeric, mmol/mol)
#   - pred.Oral_Semaglutide: Predicted 12-month HbA1c on oral semaglutide (numeric, mmol/mol)
#   - pred.Inje_Semaglutide: Predicted 12-month HbA1c on injectable semaglutide (numeric, mmol/mol)
#   - pred.TZP: Predicted 12-month HbA1c on tirzepatide (numeric, mmol/mol)

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

  ## Predictions from the model files, which need to be in the working directory:
  ##   - fivedrugmodel_5knot_share_20230823.Rdata: the 5-drug model (m1.5.final)
  ##   - model_update_recal_large.rds: the update, an intercept for each of the
  ##     main drug classes (the GLP1 subtypes included)
  ##   - model_update_recal_large_Oral_Semaglutide.rds and
  ##     model_update_recal_large_Inje_Semaglutide.rds: the semaglutide update, an
  ##     intercept for each semaglutide form
  ##   Each intercept is added to the published model's prediction of the file's
  ##   base_drugclass (GLP1 for the GLP1 subtypes and semaglutide).
  ## These give the same predictions as the online calculator.

  # load required libraries
  library(rms)

  # load model
  load("fivedrugmodel_5knot_share_20230823.Rdata") # m1.5.final

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

  # Update: an intercept for each main drug class (the GLP1 subtypes included) and
  # for each semaglutide form, on top of the published model's prediction of its
  # base drug class
  #   updated prediction = pred.<base_drugclass> (published) + intercept
  model_update <- bind_rows(lapply(
    c("model_update_recal_large.rds",
      "model_update_recal_large_Oral_Semaglutide.rds",
      "model_update_recal_large_Inje_Semaglutide.rds"),
    function(file) readRDS(file) %>% select(drugclass, base_drugclass, intercept)
  ))

  for (i in seq_len(nrow(model_update))) {
    analysis_dataset[[paste0("pred.", model_update$drugclass[i])]] <-
      analysis_dataset[[paste0("pred.", model_update$base_drugclass[i])]] + model_update$intercept[i]
  }

  # the published GLP1 prediction is not a drug class of this model version
  analysis_dataset$pred.GLP1 <- NULL

  # Tirzepatide: no intercept was estimated. As in the online calculator, it is the
  # injectable semaglutide prediction minus 5 mmol/mol, from a trial comparison of
  # tirzepatide 10 mg against semaglutide 1.0 mg (Karagiannis et al. Diabetologia 2024)
  tzp_offset <- -5
  analysis_dataset$pred.TZP <- analysis_dataset$pred.Inje_Semaglutide + tzp_offset
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

# drug classes not in your data have no data to be recalibrated, so they keep their
# original predictions
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
      subtitle = if (!is.null(plot_subtitle)) str_wrap(plot_subtitle, max(30, round(26 * facet_ncol(n_panels))))
    )
}

# side of the square pairwise panels, in inches: the same in every script, and
# wide enough for the longest drug pair name on one line
# ("Oral Semaglutide vs Dulaglutide")
pairwise_panel_size <- 3.4


#############################################################################

# Step 3: Evaluate performance ----------------------------------------------

## The analyses ----

## Each analysis validates a few drug classes against each other (see the top of
## this script). It is a list with:
##   - name, label: used in the file names and in the figure titles
##   - drugs: its drug classes in your data; the predicted best drug is chosen
##     from these
##   - focus: the drug classes it validates (drug 1 in its pairwise panels, and in
##     the predicted vs observed figures)
##   - pooled: GLP-1RA drug classes the focus drug class is compared with as one
##     pooled pairwise estimate
##   - data: the initiations of its drug classes
##   - done: FALSE when your data has no initiations of the drug classes it
##     validates, or fewer than two of its drug classes (a message says why); its
##     figures are then left out

# helper: an analysis of the drug classes of analysis_drugs in your data, which
# validates those of focus. It is done when your data has initiations of a focus
# drug class, and of at least two of its drug classes.
make_analysis <- function(name, label, analysis_drugs, focus, pooled = character(0)) {
  present <- intersect(analysis_drugs, drugs)
  focus   <- intersect(focus, present)
  pooled  <- intersect(pooled, drugs)

  reason <- if (length(focus) == 0) {
    "no initiations of the drug classes validated"
  } else if (length(present) < 2) {
    "fewer than two drug classes with initiations"
  } else {
    NA_character_
  }

  if (!is.na(reason)) cat("Analysis", name, "left out:", reason, "\n")

  list(
    name   = name,
    label  = label,
    drugs  = present,
    focus  = focus,
    pooled = pooled,
    data   = analysis_dataset %>%
      filter(drugclass %in% present) %>%
      mutate(drugclass = factor(drugclass, levels = present)),
    done   = is.na(reason)
  )
}

# a) The oral drug classes against each other (oral semaglutide is drug 1)
oral <- make_analysis("oral", "oral drug classes", oral_drugs, focus = oral_drugs)

# b) Each injectable drug class against the oral drug classes: its overall benefit
#    and best drug are against comparator_drugs, and oral semaglutide comes in as
#    one pooled pairwise estimate
dulaglutide      <- make_analysis("Dulaglutide", drug_labels[["Dulaglutide"]],
                                  c("Dulaglutide", comparator_drugs), focus = "Dulaglutide",
                                  pooled = "Oral_Semaglutide")
liraglutide      <- make_analysis("Liraglutide", drug_labels[["Liraglutide"]],
                                  c("Liraglutide", comparator_drugs), focus = "Liraglutide",
                                  pooled = "Oral_Semaglutide")
inje_semaglutide <- make_analysis("Inje_Semaglutide", drug_labels[["Inje_Semaglutide"]],
                                  c("Inje_Semaglutide", comparator_drugs), focus = "Inje_Semaglutide",
                                  pooled = "Oral_Semaglutide")
tzp              <- make_analysis("TZP", drug_labels[["TZP"]],
                                  c("TZP", comparator_drugs), focus = "TZP",
                                  pooled = "Oral_Semaglutide")

# c) The injectable drug classes against each other: one pooled pairwise estimate
#    for each pair of injectable drug classes in your data
injectable_pooled <- intersect(injectable_drugs, drugs)
injectable_pairs  <- if (length(injectable_pooled) >= 2) combn(injectable_pooled, 2, simplify = FALSE) else list()
if (length(injectable_pairs) == 0) {
  cat("Injectable drug classes against each other left out: fewer than two injectable drug classes with initiations\n")
}

# d) Predicted vs observed of every drug class validated by the analyses done (one
#    panel each, the oral drug classes first)
pred_obs_drugs <- intersect(model_drugs, unlist(lapply(
  list(oral, dulaglutide, liraglutide, inje_semaglutide, tzp),
  function(a) if (a$done) a$focus
)))


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

# - Pairwise figures:
#   - oral: the oral drug classes against each other
#   - one per injectable drug class: that drug class against the comparator drug
#     classes, and against oral semaglutide as one pooled estimate
#   - injectables: the injectable drug classes against each other, each pair as one
#     pooled estimate
#   A pooled estimate is a single group of all patients of the two drug classes (one
#   point in its panel). Two GLP-1RA drug classes have (nearly) the same predicted
#   benefit for every patient, so splitting them into groups by predicted benefit
#   means nothing.

# Any other pair with the same predicted benefit for every patient (the range of the
# group means below constant_benefit_tolerance) is also shown as a single group.
constant_benefit_tolerance <- 0.1 # mmol/mol

# helper: pairwise calibration of every pair of the drug classes of an analysis,
# with the predictions with the given prefix, grouped by predicted benefit
pairwise_grouped <- function(analysis, prefix) {
  cateval::pairwise_calibration(
    data = analysis$data,
    drug_var = "drugclass",
    drugs = analysis$drugs,
    prediction_vars = paste0(prefix, analysis$drugs),
    outcome_var = outcome_var,
    # 1 (a single group) is only used for pairs with the same predicted benefit for every patient
    cal_groups = c(1, cal_groups),
    adjustment_var = unique(c(adjustment_vars, prognostic_vars))
  )
}

# helper: pooled pairwise calibration (a single group of all patients) of drug1
# against drug2, with the predictions with the given prefix
pairwise_pooled <- function(drug1, drug2, prefix) {
  pair_drugs <- c(drug1, drug2)
  cateval::pairwise_calibration(
    data = analysis_dataset %>%
      filter(drugclass %in% pair_drugs) %>%
      mutate(drugclass = factor(drugclass, levels = pair_drugs)),
    drug_var = "drugclass",
    drugs = pair_drugs,
    prediction_vars = paste0(prefix, pair_drugs),
    outcome_var = outcome_var,
    cal_groups = 1,
    adjustment_var = unique(c(adjustment_vars, prognostic_vars))
  )
}

# helper: for each pair, the grouping shown: the largest number of groups with at
# least 100 patients per drug in every group, or a single group when the pair has
# the same predicted benefit for every patient
choose_grouping <- function(pairwise_results) {
  pairwise_results %>%
    group_by(drug1, drug2) %>%
    mutate(
      constant = diff(range(mean[n_groups == max(cal_groups)], na.rm = TRUE)) < constant_benefit_tolerance
    ) %>%
    group_by(drug1, drug2, n_groups) %>%
    mutate(min_val = min(n_drug1, n_drug2)) %>%
    group_by(drug1, drug2) %>%
    mutate(
      best_grouping = suppressWarnings(max(n_groups[n_groups %in% cal_groups & min_val >= 100])),
      select_grouping = case_when(
        constant ~ 1,
        is.finite(best_grouping) ~ best_grouping,
        TRUE ~ min(cal_groups)
      )
    ) %>%
    ungroup() %>%
    filter(n_groups == select_grouping)
}

# helper: panel title of each pair
pair_title <- function(drug1, drug2) {
  paste0(recode(drug1, !!!drug_labels), " vs ", recode(drug2, !!!drug_labels))
}

# helper: pairwise figure of an analysis: its focus drug classes (drug 1) against
# its other drug classes, then the pooled estimates against its pooled drug classes
pairwise_figure <- function(analysis, prefix, step_title) {
  if (!analysis$done) return(NULL)

  grouped <- pairwise_grouped(analysis, prefix) %>%
    filter(drug1 %in% analysis$focus) %>%
    choose_grouping() %>%
    arrange(match(drug1, analysis$drugs), match(drug2, analysis$drugs)) %>%
    mutate(title = pair_title(drug1, drug2))

  pooled <- bind_rows(lapply(analysis$pooled, function(other) {
    pairwise_pooled(analysis$focus, other, prefix)
  }))
  if (nrow(pooled) > 0) pooled <- pooled %>% mutate(title = pair_title(drug1, drug2))

  plot_data <- bind_rows(grouped, pooled)
  plot_data$title <- factor(plot_data$title, levels = unique(plot_data$title))

  pairwise_plot(plot_data, paste("Pairwise Drug Calibration:", analysis$label), step_title)
}

# helper: pairwise figure of the injectable drug classes in your data against each
# other: one pooled estimate (a single point) per pair, so one panel per pair
# (1 panel for 2 injectable drug classes, 3 for 3, 6 for 4)
pairwise_figure_injectables <- function(prefix, step_title) {
  if (length(injectable_pairs) == 0) return(NULL)

  plot_data <- bind_rows(lapply(injectable_pairs, function(pair) {
    pairwise_pooled(pair[1], pair[2], prefix)
  })) %>%
    mutate(title = pair_title(drug1, drug2))
  plot_data$title <- factor(plot_data$title, levels = unique(plot_data$title))

  pairwise_plot(plot_data, "Pairwise Drug Calibration: injectable drug classes", step_title)
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

# The recommended (predicted best) drug class is chosen from the drug classes of the
# analysis in your data, as every patient needs patients on their best drug class
# to be matched with.

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


# helper: overall benefit calibration of an analysis, with the predictions with the
# given prefix. Returns the matched pairs (and matchit objects, for the love
# plot), the mean observed benefit of following the model, the calibration slope /
# intercept and the figures.
# In an injectable drug class analysis most matched pairs are of the comparator
# drug classes, which the oral analysis covers already. Its key result is the
# matched pairs of the injectable's own initiations (its panel in the figure by drug
# class), whose mean observed benefit is overall_focus (NULL for the oral analysis,
# where these are all the pairs).
overall_benefit <- function(analysis, prefix, step_title) {
  if (!analysis$done) return(NULL)

  # match concordant vs discordant patients (returns the concordant, discordant
  # and combined matched populations, plus the matchit objects for love plots).
  # The matching can be changed here, e.g. if a love plot shows poor balance.
  matched <- cateval::match_benefit_pairs(
    data = analysis$data,
    drug_var = "drugclass",
    outcome_var = outcome_var,
    pred_cols = paste0(prefix, analysis$drugs),
    matching_var = matching_vars_prog,
    caliper = c(prehba1c = 0.4),
    seed = 19840503
  )

  # group the matched pairs by predicted benefit (each cal_groups setting), for each population
  deciles_conc <- cateval::benefit_by_group(matched$concordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
  deciles_disc <- cateval::benefit_by_group(matched$discordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
  deciles_both <- cateval::benefit_by_group(matched$combined,   cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

  # calibration split by drug class (one panel per drug class of the analysis)
  by_drug <- cateval::benefit_by_drug(matched$combined, drugs = analysis$drugs, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

  title_end <- paste0(": ", analysis$label, " (", tolower(step_title), ")")

  list(
    matched     = matched,
    # mean observed benefit of following the model, and calibration slope / intercept
    overall     = cateval::benefit_overall(matched$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed),
    # the same, in the matched pairs of the initiations of the drug classes validated
    # (index_drug: the drug class the pair is anchored on)
    overall_focus = if (!setequal(analysis$focus, analysis$drugs)) {
      cateval::benefit_overall(matched$combined %>% filter(index_drug %in% analysis$focus),
                               bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
    },
    calibration = cateval::benefit_citl_slope(matched$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed),
    # calibration plot for each matched population, side by side
    plot = plot_calibration(deciles_conc, title = "Concordant matched to discordant") +
      plot_calibration(deciles_disc, title = "Discordant matched to concordant") +
      plot_calibration(deciles_both, title = "Combined") +
      patchwork::plot_layout(nrow = 1) +
      patchwork::plot_annotation(title = paste0("Overall Benefit Calibration", title_end),
                                 theme = theme(plot.title = element_text(size = 20))),
    plot_by_drug = plot_calibration_by_drug(by_drug, title = paste0("Overall Benefit Calibration by drug class", title_end))
  )
}


## Predicted vs Observed Calibration ----

## Predicted vs observed 12-month HbA1c, with the prediction for the drug class
## initiated:
##   - the oral drug classes: all their initiations in one panel (_oral_pred_obs)
##   - every drug class validated, oral and injectable: one panel per drug class
##     (_pred_obs_by_drug), which shows each drug class's offset (what
##     recalibration-in-the-large corrects)
##
## Interpretation of plots:
# - Points close to the green line (y=x) indicate good calibration
# - Points above the green line indicate underprediction (observed HbA1c greater than predicted)
# - Points below the green line indicate overprediction (observed HbA1c less than predicted)

# helper: predicted vs observed figure of the initiations of pred_drugs, with the
# predictions with the given prefix: all in one panel, or one panel per drug class
# (by_drug)
plot_pred_obs <- function(pred_drugs, title, prefix, step_title, by_drug = FALSE) {
  if (length(pred_drugs) == 0) return(NULL)

  pred_obs_data <- analysis_dataset %>%
    filter(drugclass %in% pred_drugs, !is.na(posthba1cfinal))
  pred_obs_data$pred_drug <- pred_received(pred_obs_data, prefix)

  # panel labels: drug class and its number of initiations
  n_drug <- table(factor(pred_obs_data$drugclass, levels = pred_drugs))
  panel_labels <- setNames(
    paste0(drug_labels[pred_drugs], " (n = ", format(as.integer(n_drug), big.mark = ",", trim = TRUE), ")"),
    pred_drugs
  )
  pred_obs_data$panel <- factor(panel_labels[as.character(pred_obs_data$drugclass)], levels = panel_labels)

  ggplot(pred_obs_data, aes(x = pred_drug, y = posthba1cfinal)) +
    geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
    { if (by_drug) facet_wrap(~ panel, ncol = facet_ncol(length(pred_drugs))) } +
    theme_classic() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 15),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 17),
      plot.title.position = "plot",
      # room between panels and for the last x axis label, which overhangs its panel
      panel.spacing.x = grid::unit(1.6, "lines"),
      plot.margin = margin(5.5, 18, 5.5, 5.5)
    ) +
    labs(
      x = expression("Predicted 12-month HbA"[1*c]*" (mmol/mol)"),
      y = expression("Observed 12-month HbA"[1*c]*" (mmol/mol)"),
      title = title,
      subtitle = paste0(step_title, "; n = ", format(nrow(pred_obs_data), big.mark = ","), " initiations")
    )
}


## Predicted optimal therapy ----

## Clinical Decision Analysis: Does Recalibration Change Treatment Recommendations?
##
## This section compares the predicted optimal drug (rank 1 best predicted benefit) across
## different prediction approaches (original, intercept recalibrated, and intercept+slope recalibrated).
## Changes in optimal drug recommendations between models indicate that recalibration alters
## clinical decision-making for individual patients. This helps identify whether calibration
## adjustments have meaningful clinical implications beyond just improving prediction accuracy.
##
## Focus: the oral drug classes, with the best drug chosen from the oral drug classes
## in your data, for all their initiations. The concordance table (Step 3.4) also
## has a row for each injectable drug class: how often it is the predicted best
## drug over the comparator drug classes.

# helper: the predicted best drug class (rank 1) of each initiation of an
# analysis, chosen from its drug classes, with the predictions with the given prefix
best_drugs <- function(analysis, prefix) {
  analysis$data %>%
    cateval::get_best_drugs(
      rank = 1,
      column_names = paste0(prefix, analysis$drugs),
      final_var_name = prefix
    )
}

# helper: distribution of the predicted best drug class (rank 1)
plot_rank1 <- function(analysis, prefix, step_title) {
  if (!analysis$done) return(NULL)

  rank1_breakdown <- best_drugs(analysis, prefix) %>%
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
      subtitle = paste0(step_title, "\nBest of the ", analysis$label)
    )
}

# helper: drug classes within 3 mmol/mol of the best, for each initiation of an analysis
tolerance3_sets <- function(analysis, prefix) {
  analysis$data %>%
    cateval::get_best_drugs(
      tolerance = 3,
      column_names = paste0(prefix, analysis$drugs),
      final_var_name = prefix
    ) %>%
    pull(paste0(prefix, "within_3_of_best_drug_name"))
}

# helper: drug classes added to and removed from each patient's tolerance-3 set
# by recalibration (original predictions vs the predictions with the given prefix)
plot_tolerance3_added_removed <- function(analysis, prefix, step_title) {
  if (!analysis$done) return(NULL)

  split_drugs <- function(x) str_split(ifelse(is.na(x), "", x), ",") %>% map(str_trim)
  original_list <- split_drugs(tolerance3_sets(analysis, "pred."))
  step_list     <- split_drugs(tolerance3_sets(analysis, prefix))

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
  plot_subtitle <- paste0("Original vs ", step_title, "\nBest of the ", analysis$label)
  if (nrow(compare) == 0) {
    return(ggplot() +
             annotate("text", x = 0, y = 0, label = "No drug classes added or removed") +
             theme_void() +
             labs(title = plot_title, subtitle = plot_subtitle))
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
        ~ . * 100 / nrow(analysis$data),
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
      subtitle = plot_subtitle,
      alpha = "Change"
    )
}


## Validation of each analysis ----

# helper: the validation of an analysis with the predictions with the given
# prefix: pairwise figure and overall benefit calibration
validate <- function(analysis, prefix, step_title) {
  list(
    pairwise = pairwise_figure(analysis, prefix, step_title),
    benefit  = overall_benefit(analysis, prefix, step_title)
  )
}

# titles of the predicted vs observed figures
pred_obs_title         <- "Predicted vs Observed: oral drug classes"
pred_obs_by_drug_title <- "Predicted vs Observed by drug class"


### Step 3.1: Original predictions ----

# a) Oral drug classes
original_oral          <- validate(oral, "pred.", "Original predictions")
original_oral_pred_obs <- plot_pred_obs(if (oral$done) oral$focus, pred_obs_title, "pred.", "Original predictions")
original_oral_rank1    <- plot_rank1(oral, "pred.", "Original predictions")

# b) Each injectable drug class against the oral drug classes
original_dulaglutide      <- validate(dulaglutide, "pred.", "Original predictions")
original_liraglutide      <- validate(liraglutide, "pred.", "Original predictions")
original_inje_semaglutide <- validate(inje_semaglutide, "pred.", "Original predictions")
original_tzp              <- validate(tzp, "pred.", "Original predictions")

# c) The injectable drug classes against each other
original_injectables_pairwise <- pairwise_figure_injectables("pred.", "Original predictions")

# d) Predicted vs observed of every drug class validated
original_pred_obs_by_drug <- plot_pred_obs(pred_obs_drugs, pred_obs_by_drug_title, "pred.", "Original predictions", by_drug = TRUE)


### Step 3.2: Recalibration-in-the-large predictions ----

# a) Oral drug classes
recal_large_oral            <- validate(oral, "pred.recal_large.", "Recalibration-in-the-large")
recal_large_oral_pred_obs   <- plot_pred_obs(if (oral$done) oral$focus, pred_obs_title, "pred.recal_large.", "Recalibration-in-the-large")
recal_large_oral_rank1      <- plot_rank1(oral, "pred.recal_large.", "Recalibration-in-the-large")
recal_large_oral_tolerance3 <- plot_tolerance3_added_removed(oral, "pred.recal_large.", "Recalibration-in-the-large")

# b) Each injectable drug class against the oral drug classes
recal_large_dulaglutide      <- validate(dulaglutide, "pred.recal_large.", "Recalibration-in-the-large")
recal_large_liraglutide      <- validate(liraglutide, "pred.recal_large.", "Recalibration-in-the-large")
recal_large_inje_semaglutide <- validate(inje_semaglutide, "pred.recal_large.", "Recalibration-in-the-large")
recal_large_tzp              <- validate(tzp, "pred.recal_large.", "Recalibration-in-the-large")

# c) The injectable drug classes against each other
recal_large_injectables_pairwise <- pairwise_figure_injectables("pred.recal_large.", "Recalibration-in-the-large")

# d) Predicted vs observed of every drug class validated
recal_large_pred_obs_by_drug <- plot_pred_obs(pred_obs_drugs, pred_obs_by_drug_title, "pred.recal_large.", "Recalibration-in-the-large", by_drug = TRUE)


### Step 3.3: Recalibration predictions ----

# a) Oral drug classes
recal_oral            <- validate(oral, "pred.recal.", "Recalibration")
recal_oral_pred_obs   <- plot_pred_obs(if (oral$done) oral$focus, pred_obs_title, "pred.recal.", "Recalibration")
recal_oral_rank1      <- plot_rank1(oral, "pred.recal.", "Recalibration")
recal_oral_tolerance3 <- plot_tolerance3_added_removed(oral, "pred.recal.", "Recalibration")

# b) Each injectable drug class against the oral drug classes
recal_dulaglutide      <- validate(dulaglutide, "pred.recal.", "Recalibration")
recal_liraglutide      <- validate(liraglutide, "pred.recal.", "Recalibration")
recal_inje_semaglutide <- validate(inje_semaglutide, "pred.recal.", "Recalibration")
recal_tzp              <- validate(tzp, "pred.recal.", "Recalibration")

# c) The injectable drug classes against each other
recal_injectables_pairwise <- pairwise_figure_injectables("pred.recal.", "Recalibration")

# d) Predicted vs observed of every drug class validated
recal_pred_obs_by_drug <- plot_pred_obs(pred_obs_drugs, pred_obs_by_drug_title, "pred.recal.", "Recalibration", by_drug = TRUE)


### Step 3.4: Summaries ----

#### Overall benefit estimate of each analysis across model versions

# helper: the overall benefit estimates of an analysis in the three steps, of all its
# matched pairs ("all") and, for an injectable drug class analysis, of the pairs of
# the injectable's own initiations (overall_focus). key_result marks the estimate
# shown in the figure: all pairs for the oral analysis, the injectable's own pairs
# for an injectable analysis.
benefit_estimates_of <- function(analysis, original, recal_large, recal) {
  if (!analysis$done) return(NULL)
  steps <- list("Original" = original, "Recalibration-in-the-large" = recal_large, "Recalibration" = recal)
  bind_rows(lapply(names(steps), function(step) {
    benefit <- steps[[step]]$benefit
    bind_rows(
      benefit$overall %>% mutate(pairs = "all", key_result = is.null(benefit$overall_focus)),
      if (!is.null(benefit$overall_focus)) {
        benefit$overall_focus %>%
          mutate(pairs = paste(paste(analysis$focus, collapse = ", "), "initiations"), key_result = TRUE)
      }
    ) %>%
      mutate(model = step)
  })) %>%
    mutate(analysis_label = analysis$label, analysis = analysis$name, n_pairs = n)
}

benefit_estimates <- bind_rows(
  benefit_estimates_of(oral, original_oral, recal_large_oral, recal_oral),
  benefit_estimates_of(dulaglutide, original_dulaglutide, recal_large_dulaglutide, recal_dulaglutide),
  benefit_estimates_of(liraglutide, original_liraglutide, recal_large_liraglutide, recal_liraglutide),
  benefit_estimates_of(inje_semaglutide, original_inje_semaglutide, recal_large_inje_semaglutide, recal_inje_semaglutide),
  benefit_estimates_of(tzp, original_tzp, recal_large_tzp, recal_tzp)
)

# (no rows when no analysis was done)
if (nrow(benefit_estimates) > 0) {
benefit_estimates <- benefit_estimates %>%
  mutate(
    model = factor(model, levels = c("Recalibration", "Recalibration-in-the-large", "Original")),
    analysis_label = factor(str_to_sentence(analysis_label), levels = unique(str_to_sentence(analysis_label)))
  )

plot_benefit_estimates <- benefit_estimates %>%
  filter(key_result) %>%
  ggplot(aes(x = observed, xmin = lower.ci, xmax = upper.ci, y = model, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_pointrange(size = 0.7) +
  facet_wrap(~ analysis_label, ncol = 1) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none", axis.text = element_text(size = 16),
        strip.text = element_text(size = 16), plot.title.position = "plot",
        plot.subtitle = element_text(size = 15)) +
  labs(
    x = expression("Observed HbA"[1*c]*" benefit estimate (mmol/mol)"),
    y = NULL,
    title = "Overall benefit estimates",
    subtitle = paste0("Oral drug classes: all matched pairs\n",
                      "Injectable drug classes: matched pairs of their own initiations")
  )
}


#### Best drug proportions and concordance

## For each analysis and set of predictions: how often each drug class is the
## predicted best (rank 1, chosen from the drug classes of the analysis), and how
## often the drug initiated was that drug (concordance). The oral analysis has a
## row per drug class plus "All"; an injectable drug class analysis has only the
## row of the injectable: how often it is the predicted best drug over the
## comparator drug classes, and how often those patients initiated it.

# helper: the concordance rows of an analysis, for the predictions with the given
# prefix
summarise_concordance <- function(analysis, prefix, step_title) {
  if (!analysis$done) return(NULL)

  d <- best_drugs(analysis, prefix) %>%
    transmute(
      recommended_drug = factor(.data[[paste0(prefix, "rank1_drug_name")]], levels = analysis$drugs),
      received_drug    = as.character(drugclass),
      concordant       = as.integer(as.character(recommended_drug) == received_drug)
    )

  total_n <- nrow(d)

  # rows kept: every drug class and "All" for the oral analysis, only the injectable
  # for an injectable drug class analysis
  rows_kept <- if (setequal(analysis$focus, analysis$drugs)) c(analysis$drugs, "All") else analysis$focus

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
      analysis = analysis$name,
      model = step_title,
      total_n = total_n,
      # share of the analysis this drug is recommended for
      pct_recommended = n_recommended / total_n * 100,
      # of those recommended this drug, the share who received it
      pct_concordant_of_recommended = ifelse(
        n_recommended > 0, n_concordant / n_recommended * 100, NA_real_
      ),
      # share of the whole analysis both recommended and given this drug; the
      # "All" row of this column is the overall concordance rate
      pct_concordant_of_total = n_concordant / total_n * 100
    ) %>%
    select(
      analysis, model, recommended_drug,
      n_recommended, total_n, pct_recommended,
      n_concordant, pct_concordant_of_recommended, pct_concordant_of_total
    ) %>%
    filter(recommended_drug %in% rows_kept)
}

# helper: the concordance rows of an analysis in the three steps
concordance_of <- function(analysis) {
  bind_rows(
    summarise_concordance(analysis, "pred.", "Original"),
    summarise_concordance(analysis, "pred.recal_large.", "Recalibration-in-the-large"),
    summarise_concordance(analysis, "pred.recal.", "Recalibration")
  )
}

concordance_summary <- bind_rows(
  concordance_of(oral),
  concordance_of(dulaglutide),
  concordance_of(liraglutide),
  concordance_of(inje_semaglutide),
  concordance_of(tzp)
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
##
## Files are named <step>_<analysis>_<figure>.pdf, and
## <step>_pred_obs_by_drug.pdf for predicted vs observed of every drug class
## validated. Only the figures of the analyses done are saved.


## Save plots to separate PDFs by validation step ----

# helper: save a pairwise figure (nothing is saved when it was not made)
save_pairwise <- function(fig, file) {
  if (!is.null(fig)) save_plot(fig, file, panel_width = pairwise_panel_size, panel_height = pairwise_panel_size)
}

# helper: save the figures of the validation of an analysis (validate()), with
# file names starting with file_stem
save_validation <- function(validation, analysis, file_stem) {
  if (!analysis$done) return(invisible(NULL))
  file_start <- paste0(output_path, "/", file_stem, "_", analysis$name)

  save_pairwise(validation$pairwise, paste0(file_start, "_pairwise.pdf"))
  save_plot(validation$benefit$plot, paste0(file_start, "_overall.pdf"),
            panel_width = 4.2, panel_height = 4.2)
  save_plot(validation$benefit$plot_by_drug, paste0(file_start, "_overall_by_drug.pdf"),
            panel_width = 3.4, panel_height = 3.4)
  # balance (love) plot, both matching directions in one figure
  save_plot(plot_loveplot(validation$benefit$matched), paste0(file_start, "_loveplot.pdf"),
            panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))
}

# helper: save a predicted vs observed figure, one panel (oral drug classes) or one
# panel per drug class (nothing is saved when it was not made)
save_pred_obs <- function(fig, file, by_drug = FALSE) {
  if (!is.null(fig)) save_plot(fig, file, panel_width  = if (by_drug) 3.6 else 5.5,
                                          panel_height = if (by_drug) 3.6 else 4.2)
}

# helper: save a predicted optimal therapy figure (nothing is saved when it was not made)
save_decision <- function(fig, file) {
  if (!is.null(fig)) save_plot(fig, file, panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
}


## STEP 1: Original ----

save_validation(original_oral, oral, "step1_original")
save_pred_obs(original_oral_pred_obs, paste0(output_path, "/step1_original_oral_pred_obs.pdf"))
save_decision(original_oral_rank1, paste0(output_path, "/step1_original_oral_rank1.pdf"))

save_validation(original_dulaglutide, dulaglutide, "step1_original")
save_validation(original_liraglutide, liraglutide, "step1_original")
save_validation(original_inje_semaglutide, inje_semaglutide, "step1_original")
save_validation(original_tzp, tzp, "step1_original")

save_pairwise(original_injectables_pairwise, paste0(output_path, "/step1_original_injectables_pairwise.pdf"))

save_pred_obs(original_pred_obs_by_drug, paste0(output_path, "/step1_original_pred_obs_by_drug.pdf"), by_drug = TRUE)


## STEP 2: Recalibration-in-the-large ----

save_validation(recal_large_oral, oral, "step2_recalibration_in_the_large")
save_pred_obs(recal_large_oral_pred_obs, paste0(output_path, "/step2_recalibration_in_the_large_oral_pred_obs.pdf"))
save_decision(recal_large_oral_rank1, paste0(output_path, "/step2_recalibration_in_the_large_oral_rank1.pdf"))
save_decision(recal_large_oral_tolerance3, paste0(output_path, "/step2_recalibration_in_the_large_oral_tolerance3_added_removed.pdf"))

save_validation(recal_large_dulaglutide, dulaglutide, "step2_recalibration_in_the_large")
save_validation(recal_large_liraglutide, liraglutide, "step2_recalibration_in_the_large")
save_validation(recal_large_inje_semaglutide, inje_semaglutide, "step2_recalibration_in_the_large")
save_validation(recal_large_tzp, tzp, "step2_recalibration_in_the_large")

save_pairwise(recal_large_injectables_pairwise, paste0(output_path, "/step2_recalibration_in_the_large_injectables_pairwise.pdf"))

save_pred_obs(recal_large_pred_obs_by_drug, paste0(output_path, "/step2_recalibration_in_the_large_pred_obs_by_drug.pdf"), by_drug = TRUE)


## STEP 3: Recalibration ----

save_validation(recal_oral, oral, "step3_recalibration")
save_pred_obs(recal_oral_pred_obs, paste0(output_path, "/step3_recalibration_oral_pred_obs.pdf"))
save_decision(recal_oral_rank1, paste0(output_path, "/step3_recalibration_oral_rank1.pdf"))
save_decision(recal_oral_tolerance3, paste0(output_path, "/step3_recalibration_oral_tolerance3_added_removed.pdf"))

save_validation(recal_dulaglutide, dulaglutide, "step3_recalibration")
save_validation(recal_liraglutide, liraglutide, "step3_recalibration")
save_validation(recal_inje_semaglutide, inje_semaglutide, "step3_recalibration")
save_validation(recal_tzp, tzp, "step3_recalibration")

save_pairwise(recal_injectables_pairwise, paste0(output_path, "/step3_recalibration_injectables_pairwise.pdf"))

save_pred_obs(recal_pred_obs_by_drug, paste0(output_path, "/step3_recalibration_pred_obs_by_drug.pdf"), by_drug = TRUE)


## Summaries ----

if (nrow(benefit_estimates) > 0) {
  save_plot(plot_benefit_estimates, paste0(output_path, "/overall_benefit_estimates_by_model_version.pdf"),
            panel_width = 5.5, panel_height = 1.6)

  # table of the average overall benefit per analysis, matched pairs and step
  # (estimate and 95% CI); key_result marks the rows in the figure
  write.csv(
    benefit_estimates %>%
      mutate(model = factor(model, levels = c("Original", "Recalibration-in-the-large", "Recalibration"))) %>%
      arrange(analysis_label, desc(key_result), model) %>%
      select(analysis, pairs, key_result, model, estimate = observed, lower_ci = lower.ci, upper_ci = upper.ci, n_pairs),
    paste0(output_path, "/overall_benefit_estimates_by_model_version.csv"),
    row.names = FALSE
  )
}

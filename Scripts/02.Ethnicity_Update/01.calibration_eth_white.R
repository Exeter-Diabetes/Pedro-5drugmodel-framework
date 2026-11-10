# Validation Script for the Exeter 5-Drug Model
#
## Validation of the post-2020 updated model in one ethnic group: CPRD
##  initiations from 2020-10-14 onwards, restricted to White patients
##
## The model validated here is the published model plus the recalibration-in-
##  the-large intercepts estimated in 01.Temporal_Update/02.calibration_new_data_post2020.R
##  (read from Model_files). "Original" below therefore means the
##  updated model; recalibration-in-the-large and recalibration are then
##  estimated within the ethnic group on top of it.
##
## The model's GLP1 class is split into two subtypes (glp1_subtypes, Part 1):
##  dulaglutide, and liraglutide (with exenatide), each with its own update.
##  Within the ethnic group the GLP1 subtypes are recalibrated together, as one
##  GLP1 class: one intercept (and slope) for GLP1, applied to each patient's
##  prediction for their own subtype.
##
## The validation (Part 5) is done with one GLP1 subtype standing for
##  GLP1 next to the other four drug classes: dulaglutide in the main analysis
##  (glp1_main) and liraglutide in a sensitivity analysis, done right after the
##  main analysis each time.
#
## The script has six parts:
## 1) Drug classes: for the model and for the analysis
## 2) Prepare your data
##    a) Load the CPRD cohort and format the variables
##    b) Characteristics table
##    c) Impute missing data
## 3) Model predictions for each drug class, in three versions:
##    a) original predictions
##    b) recalibration-in-the-large (a drug-specific intercept)
##    c) recalibration (a drug-specific intercept and slope)
## 4) Figure code (needs to be run; no need to read it)
## 5) Validation, main and sensitivity analysis:
##    a) pairwise drug calibration, overall benefit calibration, predicted vs
##       observed calibration and predicted optimal therapy
##    b) save the figures. The files sort by version of the predictions:
##         1_original_..., 2_recalibration_in_the_large_..., 3_recalibration_...,
##       and within each, by figure: 01_pairwise, 02_overall, 03_overall_by_drug,
##       04_loveplot, 05_pred_obs, 06_rank1, 07_tolerance3. Then the summaries
##       (4_summary_...). The sensitivity analysis has the suffix _sens_liraglutide.
## 6) Supplementary: predicted optimal therapy with the GLP1 subtypes as separate
##    drug classes (5_supplement_...)


#############################################################################

# load required libraries
library(tidyverse)
library(tableone)
# library(remotes)
# remotes::install_github("PM-Cardoso/cateval")
library(cateval) # Custom R package developed by Pedro for calibration evaluation


#############################################################################

# Part 1: Drug classes ------------------------------------------------------

# glp1_subtypes: the model has one GLP1 class, split here into subtypes by
#   drug_substance. Each subtype starts from the model's GLP1 prediction plus its
#   own post-2020 update; within the ethnic group they are recalibrated together
#   as one GLP1 class. Exenatide is grouped with liraglutide.
glp1_subtypes <- list(
  "Dulaglutide" = c("Dulaglutide"),
  "Liraglutide" = c("Liraglutide", "Exenatide", "Exenatide prolonged-release")
)

# model_drugs: the drug classes the model makes predictions for; the predicted best
#   drug is always chosen from all of these (in Part 5, all but the GLP1
#   subtype not standing for GLP1).
# drugs: the drug classes analysed, which can be a subset of model_drugs (see
#   example_validation_script_v2.R). Only the drug classes in drugs are recalibrated
#   and have pairwise and per-drug plots.
model_drugs <- c("SGLT2", names(glp1_subtypes), "DPP4", "SU", "TZD")
drugs <- model_drugs

# glp1_main: the GLP1 subtype that stands for GLP1 next to the other drug classes
#   in the main analysis (Part 5); the other subtype stands for GLP1 in the
#   sensitivity analysis. Choose the most common subtype in your data (the
#   initiations per drug class are printed in Part 2a).
glp1_main <- "Dulaglutide"


#############################################################################

# Part 2: Prepare your data ------------------------------------------------------

# Part 2a: Load the CPRD cohort, and format the variables ----
# Note: This portion of code can be adapted based on your data source and structure

# REQUIRED VARIABLES
# Your dataset must include the following variables for validation:
#
# OUTCOME VARIABLES:
#   - posthba1cfinal: 12-month HbA1c outcome (numeric, mmol/mol)
#   - hba1cmonth: Months from treatment initiation to outcome measurement (numeric, 3-15 months)
#
# TREATMENT VARIABLES:
#   - drugclass: Drug class being initiated (categorical)
#     Valid values: "SGLT2", "Dulaglutide", "Liraglutide", "DPP4", "SU", or "TZD"
#     (GLP1 split into its subtypes, glp1_subtypes in Part 1)
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


# load required libraries
library(aurum)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

# First initiation date of the cohort: two years after the development cohort
# closed (2018-10-13), so no initiation here was used to fit the model
cohort_start <- as.Date("2020-10-14")

# Semaglutide substances (GLP1 drug_class in the cohort), left out here:
# semaglutide is added in 03.Semaglutide_Update
sema_substances <- c("Low-dose semaglutide", "High-dose semaglutide",
                     "Oral semaglutide", "Semaglutide, dose unclear")

# GLP1 subtype of each GLP1 drug_substance (glp1_subtypes, Part 1)
glp1_lookup <- setNames(rep(names(glp1_subtypes), lengths(glp1_subtypes)),
                        unlist(glp1_subtypes, use.names = FALSE))

# set up dataset (cohort cached by 01.cohort_scripts.R)
analysis_raw <- analysis_raw %>%
  analysis$cached("ethnicity") %>%
  filter(dstartdate >= cohort_start) %>%
  # select variables needed (to shorten the load time)
  select(
    patid, dstartdate,
    gender, dstartdate_age,
    t2dmduration, ethnicity_5cat,
    drug_substance, drug_class,
    imd_decile, smoking_cat,
    prebmi, prehba1c, preegfr, pretotalcholesterol, prehdl, prealt,
    drugline_all, MFN, SGLT2, GLP1, DPP4, TZD, SU, hba1cmonth,
    posthba1cfinal
  ) %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

# every GLP1 substance other than semaglutide must be in a GLP1 subtype
unassigned_glp1 <- setdiff(unique(analysis_raw$drug_substance[analysis_raw$drug_class == "GLP1"]),
                           c(names(glp1_lookup), sema_substances))
if (length(unassigned_glp1) > 0) {
  stop("GLP1 substances in no GLP1 subtype (glp1_subtypes, Part 1): ",
       paste(unassigned_glp1, collapse = ", "))
}

# Caterogical variables must be properly formatted as factors with the specified levels.
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
    smoke = factor(smoke,
                   levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),

    # Deprivation (IMD quintile)
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

    # Drug class: GLP1 split into its subtypes (semaglutide is in none, so it is
    # dropped below)
    drugclass = ifelse(drug_class == "GLP1", glp1_lookup[drug_substance], drug_class),
    drugclass = factor(drugclass,
                       levels = drugs)
  ) %>%
  group_by(pated) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  select(-row) %>%
  select(all_of(c(
    "pated", "dstartdate", "agetx", "sex", "t2dmduration", "ethnicity",
    "drug_substance", "drugclass",
    "imd5", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  )
  )) %>%
  as.data.frame() %>%
  # the model drug classes only; semaglutide is added in a later script
  filter(!is.na(drugclass)) %>%
  filter(!(drug_substance %in% sema_substances))

cat("cohort: initiations from", format(cohort_start), "to",
    format(max(analysis_dataset$dstartdate, na.rm = TRUE)), "\n")
table(analysis_dataset$ethnicity, analysis_dataset$drugclass)


# Part 2b: Characteristics table ----

# This table describes the validation cohort at treatment initiation, stratified by
# initiated drug class (the drug classes in drugs).
# It includes baseline clinical characteristics, missingness indicators, outcome timing,
# observed 12-month HbA1c, and observed response (posthba1cfinal - prehba1c).


# Add observed response and missingness indicators
analysis_dataset_drugclass <- analysis_dataset %>%
  filter(ethnicity == "White") %>% # White cohort only
  mutate(
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
  "sex", "ethnicity", "imd5", "smoke",
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
  "sex", "ethnicity", "imd5", "smoke", "drugline", "ncurrtx",
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
  addOverall = TRUE
)

## Export table
characteristics_table_df <- print(
  characteristics_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  includeNA = TRUE
)

output_path <- "Scripts/02.Ethnicity_Update/Outputs/White"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

write.csv(
  characteristics_table_df,
  paste0(output_path, "/characteristics_by_drugclass.csv"),
  row.names = TRUE
)



# Part 2c: Impute missing data ----

# Imputation of missing data can be performed using multiple imputation methods
# The following code generates m=20 imputed datasets and then averages the imputed values
# to create a single complete dataset with mean imputed values for originally missing data

# Example imputation code (using mice package):

# load require libraries
library(mice)

# Define variables to exclude from imputation (e.g., IDs)
# Modify this vector based on your dataset
vars_to_exclude <- c("pated", "dstartdate", "drug_substance", "drugclass",
                     "hba1cmonth", "posthba1cfinal")

# Filter to only those that actually exist in your dataset
vars_to_exclude <- vars_to_exclude[vars_to_exclude %in% names(analysis_dataset)]

# Create dataset for imputation (excluding specified variables)
vars_to_impute <- setdiff(names(analysis_dataset), vars_to_exclude)
data_for_imputation <- analysis_dataset[, vars_to_impute]

# Identify which values are missing in the dataset to be imputed
missing_indicator <- is.na(data_for_imputation)

# Perform multiple imputation (m=20 imputed datasets)
mice_m = 3; mice_maxit = 3
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

# The cohort is restricted to this ethnic group only now, after imputation on
# the full post-2020 cohort (so imputed values are the same as in the temporal
# update), and initiations with a missing outcome HbA1c are dropped
analysis_dataset <- analysis_dataset %>%
  filter(ethnicity == "White") %>%
  filter(!is.na(posthba1cfinal))

cat("White cohort:", nrow(analysis_dataset), "initiations\n")
table(analysis_dataset$drugclass)


#############################################################################

# Part 3: Model predictions ------------------------------------------------

## Predictions from the model can be obtained from the 5-drug model calculator (test version)
## available at: https://pm-cardoso.shinyapps.io/test_T2DTSM/)

# The model predictions should be in the following format:
#   - pred.SGLT2: Predicted 12-month HbA1c on SGLT2 inhibitor (numeric, mmol/mol)
#   - pred.GLP1: Predicted 12-month HbA1c on GLP-1 receptor agonist (numeric, mmol/mol)
#   - pred.DPP4: Predicted 12-month HbA1c on DPP4 inhibitor (numeric, mmol/mol)
#   - pred.SU: Predicted 12-month HbA1c on Sulfonylurea (numeric, mmol/mol)
#   - pred.TZD: Predicted 12-month HbA1c on Thiazolidinedione (numeric, mmol/mol)
# Each GLP1 subtype (pred.Dulaglutide, pred.Liraglutide) starts from pred.GLP1.


## When bringing the predictions back into this script, ensure that the variable names
## and formatting match those specified above.

# load required libraries
library(rms)

# load model
load("Model_files/fivedrugmodel_5knot_share_20230823.Rdata") # m1.5.final

# the model needs a data.frame (not a tibble) to make predictions
analysis_dataset <- as.data.frame(analysis_dataset)

# Generate predictions
analysis_dataset <- analysis_dataset %>%
  mutate(
    pred.SGLT2 = predict(m1.5.final, newdata = analysis_dataset %>% mutate(drugclass = "SGLT2")),
    pred.GLP1 = predict(m1.5.final, newdata = analysis_dataset %>% mutate(drugclass = "GLP1")),
    pred.DPP4 = predict(m1.5.final, newdata = analysis_dataset %>% mutate(drugclass = "DPP4")),
    pred.SU = predict(m1.5.final, newdata = analysis_dataset %>% mutate(drugclass = "SU")),
    pred.TZD = predict(m1.5.final, newdata = analysis_dataset %>% mutate(drugclass = "TZD"))
  )

# Apply the post-2020 update estimated in
# 01.Temporal_Update/02.calibration_new_data_post2020.R: an intercept for each drug
# class, the GLP1 subtypes included, on top of the published model's prediction
# of its base drug class (GLP1 for the GLP1 subtypes)
#   updated prediction = pred.<base_drugclass> + intercept
model_update <- readRDS("Model_files/model_update_recal_large.rds")
print(model_update)

for (i in seq_len(nrow(model_update))) {
  analysis_dataset[[paste0("pred.", model_update$drugclass[i])]] <-
    analysis_dataset[[paste0("pred.", model_update$base_drugclass[i])]] + model_update$intercept[i]
}


## We then generate versions of the predictions that are recalibrated using intercept-only
## and intercept + slope adjustments.

### Recalibration
# drugs and model_drugs are set at the top of the script (Part 1)
outcome_var <- "posthba1cfinal"
pred_prefix <- "pred."

# drug class names as shown in the plots
drug_labels <- c("SGLT2" = "SGLT2i", "Dulaglutide" = "Dulaglutide", "Liraglutide" = "Liraglutide",
                 "DPP4" = "DPP4i", "SU" = "SU", "TZD" = "TZD")

# drug class colours (Okabe-Ito), as in the online calculator: the GLP1 subtypes
# take the GLP-1RA colour
drug_colours <- c("SGLT2" = "#E69F00", "Dulaglutide" = "#56B4E9", "Liraglutide" = "#56B4E9",
                  "DPP4" = "#0072B2", "SU" = "#CC79A7", "TZD" = "#D55E00")

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

# The GLP1 subtypes are recalibrated together, as one GLP1 class: they share one
# intercept (and slope), applied to each patient's prediction for their own
# subtype, so the post-2020 difference between the subtypes is kept.
recal_class <- function(drug) ifelse(drug %in% names(glp1_subtypes), "GLP1", drug)
recal_classes <- unique(recal_class(drugs))

# Anchor each patient on the prediction for the drug they actually received (for
# the GLP1 subtypes, their own subtype); this is the offset (intercept model) /
# covariate (intercept + slope model).
recal_data <- analysis_dataset
recal_data$drugclass <- factor(recal_class(as.character(recal_data$drugclass)), levels = recal_classes)
pred_matrix <- as.matrix(analysis_dataset[, paste0(pred_prefix, drugs)])
recal_data$pred_assigned <- pred_matrix[
  cbind(seq_len(nrow(pred_matrix)), match(analysis_dataset$drugclass, drugs))
]

# helper: predict the recalibrated value for every patient under a given drug
predict_for_drug <- function(fit, drug) {
  nd <- analysis_dataset
  nd$drugclass <- factor(recal_class(drug), levels = recal_classes)
  nd$pred_assigned <- nd[[paste0(pred_prefix, drug)]]
  as.numeric(predict(fit, newdata = nd))
}

# helper: each patient's prediction (columns <prefix><drug>) for the drug they
# actually received
pred_received <- function(data, prefix) {
  pm <- as.matrix(data[, paste0(prefix, drugs)])
  pm[cbind(seq_len(nrow(pm)), match(data$drugclass, drugs))]
}

#### Recalibration-in-the-large (drug-specific intercept, GLP1 subtypes pooled, shared prognostic effect)
formula_recal_large <- as.formula(
  paste0(outcome_var, " ~ 0 + drugclass", prognostic_terms, " + offset(pred_assigned)")
)
fit_recal_large <- lm(formula_recal_large, data = recal_data)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal_large.", drug)]] <- predict_for_drug(fit_recal_large, drug)
}
recalib_recal_large_params <- summarize_lm(fit_recal_large, "recal_large")

#### Recalibration (drug-specific intercept & slope, GLP1 subtypes pooled, shared prognostic effect)
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

# drug classes not in drugs have no data to be recalibrated, so they keep their
# original predictions
for (drug in setdiff(model_drugs, drugs)) {
  analysis_dataset[[paste0("pred.recal_large.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
  analysis_dataset[[paste0("pred.recal.", drug)]] <- analysis_dataset[[paste0(pred_prefix, drug)]]
}

# Combine and save recalibration parameters
recalib_params_full <- bind_rows(recalib_original_params, recalib_recal_large_params, recalib_recal_params)

write_csv(recalib_params_full, paste0(output_path, "/recalibration_params_full.csv"))


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
# (pred_obs_limits()), so the figures of the three steps can be compared.
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

# The overall benefit estimates (rows with key_result) in the three steps, one
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
      model = factor(model, levels = c("Recalibration", "Recalibration-in-the-large", "Original")),
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

# Part 5: Validation, main analysis and GLP1 sensitivity analysis ----

# The model has one GLP1 class, so each analysis below is done with one GLP1
# subtype standing for GLP1 next to the other drug classes: initiations of the
# other GLP1 subtype are left out, and the predicted best drug is chosen without
# it. Dulaglutide (glp1_main, Part 1) is the main analysis; liraglutide
# (glp1_sens) is the sensitivity analysis, done right after the main analysis
# each time and saved with the suffix _sens_<subtype>. Predicted vs observed
# does not compare drug classes, so it is done once with all initiations.

glp1_sens <- setdiff(names(glp1_subtypes), glp1_main)

# helper: the initiations and drug classes of the analysis with the given GLP1
# subtype standing for GLP1 (suffix: added to the names of its output files)
glp1_analysis <- function(subtype, suffix) {
  other_subtypes <- setdiff(names(glp1_subtypes), subtype)
  analysis_drugs <- setdiff(drugs, other_subtypes)
  list(
    label       = paste(subtype, "for GLP1"),
    suffix      = suffix,
    data        = analysis_dataset %>%
      filter(drugclass %in% analysis_drugs) %>%
      mutate(drugclass = factor(drugclass, levels = analysis_drugs)),
    drugs       = analysis_drugs,
    model_drugs = setdiff(model_drugs, other_subtypes)
  )
}

main <- glp1_analysis(glp1_main, "")
sens <- glp1_analysis(glp1_sens, paste0("_sens_", tolower(glp1_sens)))

cat("Main analysis (", main$label, "): ", nrow(main$data), " initiations\n", sep = "")
cat("Sensitivity analysis (", sens$label, "): ", nrow(sens$data), " initiations\n", sep = "")


# Part 5a: Evaluate performance ----------------------------------------------

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
                     "sex", "smoke", "ncurrtx", "drugline", "imd5")

# - Pairwise figures: the drug classes of each analysis against each other

# helper: pairwise calibration of the drug classes of an analysis, with the
# predictions with the given prefix
pairwise_calibration_for <- function(analysis, prefix) {
  cateval::pairwise_calibration(
    data = analysis$data,
    drug_var = "drugclass",
    drugs = analysis$drugs,
    prediction_vars = paste0(prefix, analysis$drugs),
    outcome_var = outcome_var,
    cal_groups = cal_groups,  # Number of groups to stratify predictions into
    adjustment_var = unique(c(adjustment_vars, prognostic_vars))
  )
}

# helper: pairwise calibration plot, one panel per pair of drug classes
plot_pairwise <- function(pairwise_results, plot_title) {
  pairwise_results %>%
    mutate(drugcombo = paste(drug1, drug2)) %>%
    group_by(drugcombo, n_groups) %>%
    mutate(min_val = min(n_drug1, n_drug2)) %>%
    ungroup(n_groups) %>%
    mutate(
      select_grouping = ifelse(min_val >= 100, n_groups, NA),
      select_grouping = max(select_grouping, na.rm = TRUE),
      select_grouping = ifelse(is.infinite(abs(select_grouping)), min(cal_groups), select_grouping)
    ) %>%
    ungroup() %>%
    filter(select_grouping == n_groups) %>%
    select(-c(drugcombo, min_val, select_grouping)) %>%
    mutate(title = paste(recode(drug1, !!!drug_labels), "vs", recode(drug2, !!!drug_labels))) %>%
    pairwise_plot(plot_title)
}

### Original predictions
pairwise_results_original      <- pairwise_calibration_for(main, "pred.")
pairwise_results_original_sens <- pairwise_calibration_for(sens, "pred.")

plot_pairwise_original      <- plot_pairwise(pairwise_results_original,
  "Pairwise Drug Calibration (original predictions)")
plot_pairwise_original_sens <- plot_pairwise(pairwise_results_original_sens,
  paste0("Pairwise Drug Calibration, ", sens$label, " (original predictions)"))


### Recalibration-in-the-large predictions
pairwise_results_recal_large      <- pairwise_calibration_for(main, "pred.recal_large.")
pairwise_results_recal_large_sens <- pairwise_calibration_for(sens, "pred.recal_large.")

plot_pairwise_recal_large      <- plot_pairwise(pairwise_results_recal_large,
  "Pairwise Drug Calibration (recalibration-in-the-large)")
plot_pairwise_recal_large_sens <- plot_pairwise(pairwise_results_recal_large_sens,
  paste0("Pairwise Drug Calibration, ", sens$label, " (recalibration-in-the-large)"))


### Recalibration predictions
pairwise_results_recal      <- pairwise_calibration_for(main, "pred.recal.")
pairwise_results_recal_sens <- pairwise_calibration_for(sens, "pred.recal.")

plot_pairwise_recal      <- plot_pairwise(pairwise_results_recal,
  "Pairwise Drug Calibration (recalibration)")
plot_pairwise_recal_sens <- plot_pairwise(pairwise_results_recal_sens,
  paste0("Pairwise Drug Calibration, ", sens$label, " (recalibration)"))


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

# Covariates used to match concordant and discordant patients. Matching is 1:1
# with replacement and exact on the predicted-best drug (done inside the
# function; the best drug is chosen from all model_drugs), with a 0.4 SD caliper
# on baseline HbA1c. These are baseline
# characteristics only: unlike adjustment_vars above, hba1cmonth (time to
# outcome measurement) describes the outcome window rather than the patient at
# baseline, so it is left out of matching.
matching_vars <- c(
  "t2dmduration", "prebmi", "prehba1c", "agetx", "prealt",
  "preegfr", "pretotalcholesterol", "prehdl",
  "sex", "smoke", "ncurrtx", "drugline", "imd5"
)
# The same covariates plus any prognostic factors in prognostic_vars. Used for
# matching in every version of the predictions (the prognostic factor is a
# cross-cutting control, not a version of its own); equal to matching_vars until you populate
# prognostic_vars.
matching_vars_prog <- unique(c(matching_vars, prognostic_vars))

# Confidence intervals: FALSE uses fast analytical intervals; set to TRUE for
# pair-level bootstrap intervals (slower) that account for matching with
# replacement.
use_bootstrap <- TRUE
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
  "drugline5+"          = "'Number of previous glucose-lowering drugs = >4'",
  "imd51 (least)"       = "'Index of multiple deprivation = 1 (least deprived)'",
  "imd52"               = "'Index of multiple deprivation = 2'",
  "imd53"               = "'Index of multiple deprivation = 3'",
  "imd54"               = "'Index of multiple deprivation = 4'",
  "imd55 (most)"        = "'Index of multiple deprivation = 5 (most deprived)'"
)

# Balance (love) plot of one matching: absolute standardised mean
# difference per matching covariate, before and after matching, in each
# matching direction. Exact-matched variables are left out.
plot_loveplot <- function(match_obj) {

  love_conc <- summary(match_obj$matchit$concordant)
  love_disc <- summary(match_obj$matchit$discordant)

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


# helper: match concordant vs discordant patients of an analysis, with the
# predictions with the given prefix (returns the concordant, discordant and
# combined matched populations, plus the matchit objects for love plots).
# match_args_disc changes the matching in the discordant-as-treated direction.
match_pairs <- function(analysis, prefix, match_args_disc = list()) {
  cateval::match_benefit_pairs(
    data = analysis$data,
    drug_var = "drugclass",
    outcome_var = outcome_var,
    pred_cols = paste0(prefix, analysis$model_drugs),
    matching_var = matching_vars_prog,
    caliper = c(prehba1c = 0.4),
    match_args_disc = match_args_disc,
    seed = 19840503
  )
}

# helper: overall benefit calibration of an analysis, with the predictions with
# the given prefix (title_suffix ends the figure titles). Returns the matched
# pairs, the mean observed benefit of following the model, the calibration
# slope / intercept and the figures.
overall_benefit <- function(analysis, prefix, title_suffix, match_args_disc = list()) {
  matched <- match_pairs(analysis, prefix, match_args_disc)

  # group the matched pairs by predicted benefit (each cal_groups setting), for each population
  deciles_conc <- cateval::benefit_by_group(matched$concordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
  deciles_disc <- cateval::benefit_by_group(matched$discordant, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
  deciles_both <- cateval::benefit_by_group(matched$combined,   cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

  # calibration plot for each matched population, then the three side by side
  plot_conc <- plot_calibration(deciles_conc, title = "Concordant matched to discordant")
  plot_disc <- plot_calibration(deciles_disc, title = "Discordant matched to concordant")
  plot_both <- plot_calibration(deciles_both, title = "Combined")

  # calibration split by drug class (one panel per drug class of the analysis)
  by_drug <- cateval::benefit_by_drug(matched$combined, drugs = analysis$drugs, cal_groups = cal_groups, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

  list(
    matched      = matched,
    overall      = cateval::benefit_overall(matched$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed),
    calibration  = cateval::benefit_citl_slope(matched$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed),
    plot_conc    = plot_conc,
    plot_disc    = plot_disc,
    plot_both    = plot_both,
    plot         = overall_figure(deciles_conc, deciles_disc, deciles_both,
                                  paste0("Overall Benefit Calibration", title_suffix)),
    plot_by_drug = plot_calibration_by_drug(by_drug, title = paste0("Overall Benefit Calibration by drug class", title_suffix))
  )
}

# Matching for every version: stricter in the discordant-as-treated direction
# (extra calipers, exact on drug line and number of current drugs, each patient
# reused at most 3 times). The original predictions are those of the updated
# model, so they are matched in the same way as the two recalibrated versions.
match_args_disc_recal <- list(
  caliper = c(prehba1c = 0.4, prebmi = 0.5, t2dmduration = 0.5, agetx = 0.5),
  exact = c("drugline", "ncurrtx"),
  reuse.max = 3, m.order = "random")


### Original predictions ----

benefit_original      <- overall_benefit(main, "pred.", " (original predictions)",
                                         match_args_disc = match_args_disc_recal)
benefit_original_sens <- overall_benefit(sens, "pred.", paste0(", ", sens$label, " (original predictions)"),
                                         match_args_disc = match_args_disc_recal)


### Recalibration-in-the-large predictions ----

benefit_recal_large      <- overall_benefit(main, "pred.recal_large.", " (recalibration-in-the-large)",
                                            match_args_disc = match_args_disc_recal)
benefit_recal_large_sens <- overall_benefit(sens, "pred.recal_large.", paste0(", ", sens$label, " (recalibration-in-the-large)"),
                                            match_args_disc = match_args_disc_recal)


### Recalibration predictions ----

benefit_recal      <- overall_benefit(main, "pred.recal.", " (recalibration)",
                                      match_args_disc = match_args_disc_recal)
benefit_recal_sens <- overall_benefit(sens, "pred.recal.", paste0(", ", sens$label, " (recalibration)"),
                                      match_args_disc = match_args_disc_recal)


### Compare the overall benefit estimate across model versions ----

# helper: the overall benefit estimate of each version of the predictions, for one analysis
combine_benefit_estimates <- function(original, recal_large, recal) {
  bind_rows(
    original$overall %>% mutate(model = "Original"),
    recal_large$overall %>% mutate(model = "Recalibration-in-the-large"),
    recal$overall %>% mutate(model = "Recalibration")
  ) %>%
    mutate(model = factor(
      model,
      levels = c("Original", "Recalibration-in-the-large", "Recalibration")
    ))
}

benefit_estimates      <- combine_benefit_estimates(benefit_original, benefit_recal_large, benefit_recal)
benefit_estimates_sens <- combine_benefit_estimates(benefit_original_sens, benefit_recal_large_sens, benefit_recal_sens)

# helper: plot of the overall benefit estimates of one analysis (its label as
# the panel title)
plot_benefit_estimates <- function(estimates, plot_title, analysis_label) {
  draw_benefit_estimates(estimates %>% mutate(analysis_label = analysis_label, key_result = TRUE),
                         plot_title = plot_title)
}


## Predicted vs Observed Calibration ----

## Predicted vs observed HbA1c (the outcome HbA1c: at 12 months, or at 6 months
## when the 12-month value is missing) of every initiation (every drug class
## of the cohort), with the prediction for the drug class initiated. The three
## figures have the same axes, so the versions of the predictions can be compared.

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

# the same axes in the three figures: the ranges in the three versions
pred_obs_axis_limits <- pred_obs_limits(lapply(c("pred.", "pred.recal_large.", "pred.recal."), pred_obs_data_of))

plot_pred_obs_original    <- draw_pred_obs(pred_obs_data_of("pred."), drugs, "Predicted vs Observed",
                                           "Original predictions", limits = pred_obs_axis_limits)
plot_pred_obs_recal_large <- draw_pred_obs(pred_obs_data_of("pred.recal_large."), drugs, "Predicted vs Observed",
                                           "Recalibration-in-the-large", limits = pred_obs_axis_limits)
plot_pred_obs_recal       <- draw_pred_obs(pred_obs_data_of("pred.recal."), drugs, "Predicted vs Observed",
                                           "Recalibration", limits = pred_obs_axis_limits)


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


# helper: the predicted best drug class (rank 1) of each patient of an analysis,
# chosen from its model_drugs with the predictions with the given prefix (in the
# column <prefix>rank1_drug_name)
best_drug <- function(analysis, prefix) {
  analysis$data %>%
    cateval::get_best_drugs(
      rank = 1,
      column_names = paste0(prefix, analysis$model_drugs),
      final_var_name = prefix
    )
}

# helper: distribution of the predicted best drug class (rank 1)
plot_rank1 <- function(analysis, prefix, step_title) {
  rank1_breakdown <- best_drug(analysis, prefix) %>%
    count(drug = .data[[paste0(prefix, "rank1_drug_name")]]) %>%
    mutate(proportion = n / sum(n) * 100)

  draw_rank1(rank1_breakdown, paste0(step_title, "; full cohort, ", analysis$label), analysis$model_drugs)
}

# helper: drug classes within 3 mmol/mol of the best, for each patient
tolerance3_sets <- function(analysis, prefix) {
  analysis$data %>%
    cateval::get_best_drugs(
      tolerance = 3,
      column_names = paste0(prefix, analysis$model_drugs),
      final_var_name = prefix
    ) %>%
    pull(paste0(prefix, "within_3_of_best_drug_name"))
}

# helper: drug classes added to and removed from each patient's tolerance-3 set
# by recalibration (original predictions vs the predictions with the given prefix)
plot_tolerance3_added_removed <- function(analysis, prefix, step_title) {
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
    complete(change_type = c("added", "removed"), drug, fill = list(n = 0))

  draw_tolerance3(compare, nrow(analysis$data),
                  paste0("Original vs ", step_title, "; ", analysis$label), analysis$model_drugs)
}

### Original predictions
plot_rank1_original      <- plot_rank1(main, "pred.", "Original predictions")
plot_rank1_original_sens <- plot_rank1(sens, "pred.", "Original predictions")

### Recalibration-in-the-large predictions
plot_rank1_recal_large      <- plot_rank1(main, "pred.recal_large.", "Recalibration-in-the-large")
plot_rank1_recal_large_sens <- plot_rank1(sens, "pred.recal_large.", "Recalibration-in-the-large")

plot_tolerance3_recal_large      <- plot_tolerance3_added_removed(main, "pred.recal_large.", "Recalibration-in-the-large")
plot_tolerance3_recal_large_sens <- plot_tolerance3_added_removed(sens, "pred.recal_large.", "Recalibration-in-the-large")

### Recalibration predictions
plot_rank1_recal      <- plot_rank1(main, "pred.recal.", "Recalibration")
plot_rank1_recal_sens <- plot_rank1(sens, "pred.recal.", "Recalibration")

plot_tolerance3_recal      <- plot_tolerance3_added_removed(main, "pred.recal.", "Recalibration")
plot_tolerance3_recal_sens <- plot_tolerance3_added_removed(sens, "pred.recal.", "Recalibration")


## Best drug proportions and concordance ----

## For each set of predictions: how often each drug class is the predicted best
## (rank 1), and how often the drug initiated was that drug (concordance)

cohort_label <- "White"

# helper: one row per drug class (plus "All") for an analysis, with the
# predictions with the given prefix
summarise_concordance <- function(analysis, prefix, approach) {
  d <- best_drug(analysis, prefix) %>%
    transmute(
      recommended_drug = factor(as.character(.data[[paste0(prefix, "rank1_drug_name")]]), levels = analysis$model_drugs),
      received_drug    = factor(as.character(drugclass), levels = analysis$drugs),
      concordant       = as.integer(as.character(recommended_drug) == as.character(received_drug))
    )

  stopifnot(!anyNA(d$recommended_drug), !anyNA(d$received_drug))

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
      cohort = cohort_label,
      analysis = analysis$label,
      approach  = approach,
      total_n   = total_n,
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
      cohort, analysis, approach, recommended_drug,
      n_recommended, total_n, pct_recommended,
      n_concordant, pct_concordant_of_recommended, pct_concordant_of_total
    )
}

concordance_summary <- bind_rows(
  summarise_concordance(main, "pred.",             "original"),
  summarise_concordance(main, "pred.recal_large.", "recal_large"),
  summarise_concordance(main, "pred.recal.",       "recal"),
  summarise_concordance(sens, "pred.",             "original"),
  summarise_concordance(sens, "pred.recal_large.", "recal_large"),
  summarise_concordance(sens, "pred.recal.",       "recal")
)

print(concordance_summary)

write_csv(
  concordance_summary,
  paste0(output_path, "/4_summary_concordance_by_drug.csv")
)


#############################################################################

# Part 5b: Save the figures, by version of the predictions --------------

## The figures are saved for each version of the predictions in turn:
##   1) Original: does the model's predicted benefit calibrate as it is?
##   2) Recalibration-in-the-large: does it calibrate with a drug-specific
##      intercept correction?
##   3) Recalibration: does it calibrate with a drug-specific intercept and slope
##      correction?
## The files sort in the same order (1_original_..., 2_recalibration_in_the_large_...,
## 3_recalibration_...), then the summaries (4_summary_...).
##
## The prognostic factor (prognostic_vars) is a cross-cutting control applied
## within every version when identified - adjusted for in the pairwise analysis
## and matched on in the overall-benefit analysis - not a version of its own.

# The sensitivity analysis is saved with the same file names plus sens$suffix.

# helper: save the pairwise and overall benefit figures of one version of the predictions for one
# analysis (file names start with file_stem and end with suffix)
save_step_figures <- function(file_stem, plot_pairwise_step, benefit, suffix) {
  save_pairwise(plot_pairwise_step, paste0(output_path, "/", file_stem, "_01_pairwise", suffix, ".pdf"))
  save_overall(benefit$plot, paste0(output_path, "/", file_stem, "_02_overall", suffix, ".pdf"))
  save_overall_by_drug(benefit$plot_by_drug, paste0(output_path, "/", file_stem, "_03_overall_by_drug", suffix, ".pdf"))

  # balance (love) plot, both matching directions in one figure
  save_plot(plot_loveplot(benefit$matched), paste0(output_path, "/", file_stem, "_04_loveplot", suffix, ".pdf"),
            panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))
}


## 1) Original predictions ----

save_step_figures("1_original", plot_pairwise_original, benefit_original, main$suffix)
save_step_figures("1_original", plot_pairwise_original_sens, benefit_original_sens, sens$suffix)
save_pred_obs(plot_pred_obs_original, paste0(output_path, "/1_original_05_pred_obs.pdf"))


## 2) Recalibration-in-the-large ----

save_step_figures("2_recalibration_in_the_large", plot_pairwise_recal_large, benefit_recal_large, main$suffix)
save_step_figures("2_recalibration_in_the_large", plot_pairwise_recal_large_sens, benefit_recal_large_sens, sens$suffix)
save_pred_obs(plot_pred_obs_recal_large, paste0(output_path, "/2_recalibration_in_the_large_05_pred_obs.pdf"))


## 3) Recalibration ----

save_step_figures("3_recalibration", plot_pairwise_recal, benefit_recal, main$suffix)
save_step_figures("3_recalibration", plot_pairwise_recal_sens, benefit_recal_sens, sens$suffix)
save_pred_obs(plot_pred_obs_recal, paste0(output_path, "/3_recalibration_05_pred_obs.pdf"))


## Save the overall-benefit plot components, matched pairs and pairwise results ----
# The individual panels are saved rather than the assembled patchworks, so any
# combination of version of the predictions / matched population can be put together
# later without re-running this script. Of the matched pairs only the columns
# needed to regroup them are kept, so groups can be recomputed later without
# repeating any matching. The pairwise results are group-level summaries, so the
# drug-pair panels can be redrawn later.

slim_pairs <- function(pairs) {
  pairs %>%
    transmute(
      index_drug      = as.character(index_drug),
      concordant_drug = as.character(concordant_drug),
      discordant_drug = as.character(discordant_drug),
      pred.benefit,
      obs.benefit
    )
}

# helper: save these for one analysis (benefits and pairwise_results: lists
# with one element per version of the predictions)
save_benefit_results <- function(benefits, pairwise_results, suffix) {
  saveRDS(
    lapply(benefits, function(b) list(conc = b$plot_conc, disc = b$plot_disc, both = b$plot_both)),
    paste0(output_path, "/overall_benefit_plots", suffix, ".rds")
  )
  saveRDS(
    lapply(benefits, function(b) list(
      conc = slim_pairs(b$matched$concordant),
      disc = slim_pairs(b$matched$discordant),
      both = slim_pairs(b$matched$combined)
    )),
    paste0(output_path, "/overall_benefit_pairs", suffix, ".rds")
  )
  saveRDS(pairwise_results, paste0(output_path, "/pairwise_results", suffix, ".rds"))
}

save_benefit_results(
  list(original = benefit_original, recal_large = benefit_recal_large, recal = benefit_recal),
  list(original = pairwise_results_original, recal_large = pairwise_results_recal_large, recal = pairwise_results_recal),
  main$suffix
)
save_benefit_results(
  list(original = benefit_original_sens, recal_large = benefit_recal_large_sens, recal = benefit_recal_sens),
  list(original = pairwise_results_original_sens, recal_large = pairwise_results_recal_large_sens, recal = pairwise_results_recal_sens),
  sens$suffix
)


## Additional Analyses: Clinical Decision Analysis ----

# helper: save a clinical decision figure
save_decision_plot <- function(plot, file_stem, suffix) {
  save_decision(plot, paste0(output_path, "/", file_stem, suffix, ".pdf"))
}

save_decision_plot(plot_rank1_original, "1_original_06_rank1", main$suffix)
save_decision_plot(plot_rank1_original_sens, "1_original_06_rank1", sens$suffix)

save_decision_plot(plot_rank1_recal_large, "2_recalibration_in_the_large_06_rank1", main$suffix)
save_decision_plot(plot_rank1_recal_large_sens, "2_recalibration_in_the_large_06_rank1", sens$suffix)
save_decision_plot(plot_tolerance3_recal_large, "2_recalibration_in_the_large_07_tolerance3", main$suffix)
save_decision_plot(plot_tolerance3_recal_large_sens, "2_recalibration_in_the_large_07_tolerance3", sens$suffix)

save_decision_plot(plot_rank1_recal, "3_recalibration_06_rank1", main$suffix)
save_decision_plot(plot_rank1_recal_sens, "3_recalibration_06_rank1", sens$suffix)
save_decision_plot(plot_tolerance3_recal, "3_recalibration_07_tolerance3", main$suffix)
save_decision_plot(plot_tolerance3_recal_sens, "3_recalibration_07_tolerance3", sens$suffix)

save_benefit_estimates(plot_benefit_estimates(benefit_estimates, "Overall benefit estimates", main$label),
                       paste0(output_path, "/4_summary_overall_benefit_estimates", main$suffix, ".pdf"))
save_benefit_estimates(plot_benefit_estimates(benefit_estimates_sens, "Overall benefit estimates", sens$label),
                       paste0(output_path, "/4_summary_overall_benefit_estimates", sens$suffix, ".pdf"))

# table of the average overall benefit per version of the predictions (estimate and 95% CI)
write.csv(
  bind_rows(
    benefit_estimates %>% mutate(analysis = main$label),
    benefit_estimates_sens %>% mutate(analysis = sens$label)
  ) %>%
    select(analysis, model, estimate = observed, lower_ci = lower.ci, upper_ci = upper.ci, n_pairs = n),
  paste0(output_path, "/4_summary_overall_benefit_estimates.csv"),
  row.names = FALSE
)


#############################################################################

# Part 6: Supplementary: GLP1 subtypes as separate drug classes ---------------

## The predicted best drug (rank 1) in the full cohort, chosen from all drug
## classes with the GLP1 subtypes as separate classes (as the updated model would
## recommend them), and how often the drug initiated was that drug (a patient
## recommended one GLP1 subtype who initiated another is discordant).

glp1_separate <- list(
  label       = "GLP1 subtypes separate",
  suffix      = "_glp1_subtypes",
  data        = analysis_dataset,
  drugs       = drugs,
  model_drugs = model_drugs
)

save_decision_plot(plot_rank1(glp1_separate, "pred.", "Original predictions"),
                   "5_supplement_06_rank1_original", glp1_separate$suffix)
save_decision_plot(plot_rank1(glp1_separate, "pred.recal_large.", "Recalibration-in-the-large"),
                   "5_supplement_06_rank1_recalibration_in_the_large", glp1_separate$suffix)
save_decision_plot(plot_rank1(glp1_separate, "pred.recal.", "Recalibration"),
                   "5_supplement_06_rank1_recalibration", glp1_separate$suffix)

concordance_supplement <- bind_rows(
  summarise_concordance(glp1_separate, "pred.",             "original"),
  summarise_concordance(glp1_separate, "pred.recal_large.", "recal_large"),
  summarise_concordance(glp1_separate, "pred.recal.",       "recal")
)
print(concordance_supplement)

write_csv(
  concordance_supplement,
  paste0(output_path, "/5_supplement_concordance_by_drug_glp1_subtypes.csv")
)

# Validation Script for the Exeter 5-Drug Model (v1, developed XXXX)
#
## This script demonstrates how to validate the 5-drug model in a new cohort
##  and how to visualise calibration results
#
## There are several stages in this script:
# 1) Prepare your dataset with baseline features and outcomes
#    a) Formatting continuous and categorical variables
#    b) Creating a characteristics table
#    c) Imputing missing data if needed
# 2) Obtain model predictions for each drug class
#    a) Generate original predictions
#    b) Generate recalibration-in-the-large
#    c) Generate recalibration
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
# remotes::install_github("PM-Cardoso/cateval")
library(cateval) # Custom R package developed by Pedro for calibration evaluation


#############################################################################

# Step 1a: Prepare your dataset ------------------------------------------------
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
#     Valid values: "SGLT2", "GLP1", "DPP4", "SU", or "TZD"
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
    drugclass = factor(drug_class,
                       levels = c("SGLT2", "GLP1", "DPP4", "SU", "TZD"))
  )


# Step 1b: Create a characteristics table ---------------------------------

# This table describes the validation cohort at treatment initiation, stratified by
# initiated drug class (SGLT2, GLP1, DPP4, SU, TZD).
# It includes baseline clinical characteristics, missingness indicators, outcome timing,
# observed 12-month HbA1c, and observed response (posthba1cfinal - prehba1c).


# Add observed response and missingness indicators
analysis_dataset_drugclass <- analysis_dataset %>%
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
  "ethnicity", "smoke",
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
  "ethnicity", "smoke", "drugline", "ncurrtx",
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
  test = FALSE
)

## Export table
characteristics_table_df <- print(
  characteristics_table,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  includeNA = TRUE,
  includeOverall = TRUE
)

output_path <- "Outputs/Validation_Procedure"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

write.csv(
  characteristics_table_df,
  paste0(output_path, "/stage_1b_characteristics_by_drugclass.csv"),
  row.names = TRUE
)


# Step 1c: Prepare your dataset ------------------------------------------------

# Imputation of missing data can be performed using multiple imputation methods
# The following code generates m=20 imputed datasets and then averages the imputed values
# to create a single complete dataset with mean imputed values for originally missing data

# Example imputation code (using mice package):

# load require libraries
library(mice)

# Define variables to exclude from imputation (e.g., IDs)
# Modify this vector based on your dataset
vars_to_exclude <- c("patid")

# Filter to only those that actually exist in your dataset
vars_to_exclude <- vars_to_exclude[vars_to_exclude %in% names(analysis_dataset)]

# Create dataset for imputation (excluding specified variables)
vars_to_impute <- setdiff(names(analysis_dataset), vars_to_exclude)
data_for_imputation <- analysis_dataset[, vars_to_impute]

# Identify which values are missing in the dataset to be imputed
missing_indicator <- is.na(data_for_imputation)

# Perform multiple imputation (m=20 imputed datasets)
mice_m = 20; mice_maxit = 20
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


#############################################################################

# Step 2: Model predictions ------------------------------------------------

## Predictions from the model can be obtained from the 5-drug model calculator (test version)
## available at: https://pm-cardoso.shinyapps.io/test_T2DTSM/)

# The model predictions should be in the following format:
#   - pred.SGLT2: Predicted 12-month HbA1c on SGLT2 inhibitor (numeric, mmol/mol)
#   - pred.GLP1: Predicted 12-month HbA1c on GLP-1 receptor agonist (numeric, mmol/mol)
#   - pred.DPP4: Predicted 12-month HbA1c on DPP4 inhibitor (numeric, mmol/mol)
#   - pred.SU: Predicted 12-month HbA1c on Sulfonylurea (numeric, mmol/mol)
#   - pred.TZD: Predicted 12-month HbA1c on Thiazolidinedione (numeric, mmol/mol)


## When bringing the predictions back into this script, ensure that the variable names
## and formatting match those specified above.


## We then generate versions of the predictions that are recalibrated using intercept-only 
## and intercept + slope adjustments.

### Recalibration
drugs = c("SGLT2", "GLP1", "DPP4", "SU", "TZD")
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
for (drug in drugs) {
  analysis_dataset[[paste0(pred_prefix, drug)]] <- predict_for_drug(fit_original, drug)
}
recalib_original_params <- summarize_lm(fit_original, "original")

# Combine and save recalibration parameters
recalib_params_full <- bind_rows(recalib_original_params, recalib_recal_large_params, recalib_recal_params)

write_csv(recalib_params_full, paste0(output_path, "/recalibration_params_full.csv"))



# Step 3: Evaluate performance ----------------------------------------------

## Pairwise drug calibration ----

## Interpretation of plots:
# - Points close to the red line (y=x) indicate good calibration
# - Points above the red line indicate underprediction (observed benefit greater than predicted)
# - Points below the red line indicate overprediction (observed benefit less than predicted)


# This methods has a few entries that can be adjusted based on your dataset:
# - cal_groups: Number of groups to stratify predictions into (e.g., 3, 5, 10)
#   Ensure that each group has sufficient patients for stable estimates (e.g., >100 per drug)
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

plot_pairwise_original <- pairwise_results_original %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), min(cal_groups), select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point(colour = "#800020") +
  geom_errorbar(colour = "#800020") +
  facet_wrap(~title, nrow = 2) +
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
  labs(
    x = "Predicted HbA1c benefit (mmol/mol)",
    y = "Observed HbA1c benefit (mmol/mol)",
    title = "Pairwise Drug Calibration (original predictions)"
  )


### Recalibration-in-the-large predictions
pairwise_results_recal_large <- cateval::pairwise_calibration(
  data = analysis_dataset,
  drug_var = "drugclass",
  drugs = drugs,
  prediction_vars = paste0("pred.recal_large.", drugs),
  outcome_var = outcome_var,
  cal_groups = cal_groups,  # Number of groups to stratify predictions into
  adjustment_var = unique(c(adjustment_vars, prognostic_vars))
)

plot_pairwise_recal_large <- pairwise_results_recal_large %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), min(cal_groups), select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point(colour = "#800020") +
  geom_errorbar(colour = "#800020") +
  facet_wrap(~title, nrow = 2) +
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
  labs(
    x = "Predicted HbA1c benefit (mmol/mol)",
    y = "Observed HbA1c benefit (mmol/mol)",
    title = "Pairwise Drug Calibration (recalibration-in-the-large)"
  )


### Recalibration predictions
pairwise_results_recal <- cateval::pairwise_calibration(
  data = analysis_dataset,
  drug_var = "drugclass",
  drugs = drugs,
  prediction_vars = paste0("pred.recal.", drugs),
  outcome_var = outcome_var,
  cal_groups = cal_groups,  # Number of groups to stratify predictions into
  adjustment_var = unique(c(adjustment_vars, prognostic_vars))
)

plot_pairwise_recal <- pairwise_results_recal %>%
  mutate(drugcombo = paste(drug1, drug2)) %>%
  group_by(drugcombo, n_groups) %>%
  mutate(min_val = min(n_drug1, n_drug2)) %>%
  ungroup(n_groups) %>%
  mutate(
    select_grouping = ifelse(min_val > 100, n_groups, NA),
    select_grouping = max(select_grouping, na.rm = TRUE),
    select_grouping = ifelse(is.infinite(abs(select_grouping)), min(cal_groups), select_grouping)
  ) %>%
  ungroup() %>%
  filter(select_grouping == n_groups) %>%
  select(-c(drugcombo, min_val, select_grouping)) %>%
  mutate(title = paste(drug1, "vs", drug2)) %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point(colour = "#800020") +
  geom_errorbar(colour = "#800020") +
  facet_wrap(~title, nrow = 2) +
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
  labs(
    x = "Predicted HbA1c benefit (mmol/mol)",
    y = "Observed HbA1c benefit (mmol/mol)",
    title = "Pairwise Drug Calibration (recalibration)"
  )



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
plot_calibration <- function(deciles, title = NULL) {
  n_pairs <- sum(deciles$n, na.rm = TRUE)
  ggplot(deciles, aes(x = pred.benefit, y = obs.benefit)) +
    geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_smooth(method = "loess", span = 1, se = TRUE, colour = "blue", linetype = "dashed", linewidth = 0.75) +
    geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, colour = "black") +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    labs(
      x = "Predicted HbA1c benefit (mmol/mol)",
      y = "Observed HbA1c benefit (mmol/mol)",
      title = title,
      subtitle = paste0("n = ", format(n_pairs, big.mark = ","), " matched pairs")
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
    paste0(as.character(facet_counts$index_drug),
           " (n = ", format(facet_counts$n_pairs, big.mark = ","),
           ", ", facet_counts$n_shown, " groups)"),
    as.character(facet_counts$index_drug)
  )
  plot_calibration(by_drug, title = title) +
    facet_wrap(~ index_drug, labeller = labeller(index_drug = facet_labels))
}


### Original predictions ----

# match concordant vs discordant patients (returns the concordant, discordant
# and combined matched populations, plus the matchit objects for love plots)
match_original <- cateval::match_benefit_pairs(
  data = analysis_dataset,
  drug_var = "drugclass",
  outcome_var = outcome_var,
  pred_cols = paste0("pred.", drugs),
  matching_var = matching_vars_prog,
  caliper = c(prehba1c = 0.4),
  seed = 19840503
)

# group the matched pairs into deciles of predicted benefit, for each population
deciles_original_conc <- cateval::benefit_by_group(match_original$concordant, cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_original_disc <- cateval::benefit_by_group(match_original$discordant, cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_original_both <- cateval::benefit_by_group(match_original$combined,   cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

# mean observed benefit of following the model, and calibration slope / intercept
overall_original     <- cateval::benefit_overall(match_original$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
calibration_original <- cateval::benefit_citl_slope(match_original$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

# calibration plot for each matched population, then the three side by side
plot_original_conc <- plot_calibration(deciles_original_conc, title = "Concordant matched to discordant")
plot_original_disc <- plot_calibration(deciles_original_disc, title = "Discordant matched to concordant")
plot_original_both <- plot_calibration(deciles_original_both, title = "Combined")

plot_overall_original <- plot_original_conc + plot_original_disc + plot_original_both +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(title = "Overall Benefit Calibration (original predictions)")

# calibration split by recommended drug class (pass the full drug list so every
# recommended class is shown, even if some are never received)
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

deciles_recal_large_conc <- cateval::benefit_by_group(match_recal_large$concordant, cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_large_disc <- cateval::benefit_by_group(match_recal_large$discordant, cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_large_both <- cateval::benefit_by_group(match_recal_large$combined,   cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

overall_recal_large     <- cateval::benefit_overall(match_recal_large$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
calibration_recal_large <- cateval::benefit_citl_slope(match_recal_large$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

plot_recal_large_conc <- plot_calibration(deciles_recal_large_conc, title = "Concordant matched to discordant")
plot_recal_large_disc <- plot_calibration(deciles_recal_large_disc, title = "Discordant matched to concordant")
plot_recal_large_both <- plot_calibration(deciles_recal_large_both, title = "Combined")

plot_overall_recal_large <- plot_recal_large_conc + plot_recal_large_disc + plot_recal_large_both +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(title = "Overall Benefit Calibration (recalibration-in-the-large)")

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

deciles_recal_conc <- cateval::benefit_by_group(match_recal$concordant, cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_disc <- cateval::benefit_by_group(match_recal$discordant, cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
deciles_recal_both <- cateval::benefit_by_group(match_recal$combined,   cal_groups = 10, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

overall_recal     <- cateval::benefit_overall(match_recal$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)
calibration_recal <- cateval::benefit_citl_slope(match_recal$combined, bootstrap = use_bootstrap, n_boot = n_boot, boot_seed = boot_seed)

plot_recal_conc <- plot_calibration(deciles_recal_conc, title = "Concordant matched to discordant")
plot_recal_disc <- plot_calibration(deciles_recal_disc, title = "Discordant matched to concordant")
plot_recal_both <- plot_calibration(deciles_recal_both, title = "Combined")

plot_overall_recal <- plot_recal_conc + plot_recal_disc + plot_recal_both +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(title = "Overall Benefit Calibration (recalibration)")

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
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Model version",
    y = "Observed HbA1c benefit estimate (mmol/mol)",
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
    pred_drug = case_when(
      drugclass == "SGLT2" ~ pred.SGLT2,
      drugclass == "GLP1"  ~ pred.GLP1,
      drugclass == "DPP4"  ~ pred.DPP4,
      drugclass == "SU"    ~ pred.SU,
      drugclass == "TZD"   ~ pred.TZD,
      TRUE ~ NA_real_
    )
  ) %>%
  select(pred_drug, posthba1cfinal) %>%
  ggplot(aes(x = pred_drug, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_classic() +
  labs(
    x = "Predicted 12-month HbA1c (mmol/mol)",
    y = "Observed 12-month HbA1c (mmol/mol)",
    title = "Predicted vs Observed (original predictions)",
    subtitle = paste0("n = ", format(sum(!is.na(analysis_dataset$posthba1cfinal)), big.mark = ","), " initiations")
  )

### Recalibration-in-the-large predictions
plot_pred_obs_recal_large <- analysis_dataset %>%
  mutate(
    pred_drug = case_when(
      drugclass == "SGLT2" ~ pred.recal_large.SGLT2,
      drugclass == "GLP1"  ~ pred.recal_large.GLP1,
      drugclass == "DPP4"  ~ pred.recal_large.DPP4,
      drugclass == "SU"    ~ pred.recal_large.SU,
      drugclass == "TZD"   ~ pred.recal_large.TZD,
      TRUE ~ NA_real_
    )
  ) %>%
  select(pred_drug, posthba1cfinal) %>%
  ggplot(aes(x = pred_drug, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_classic() +
  labs(
    x = "Predicted 12-month HbA1c (mmol/mol)",
    y = "Observed 12-month HbA1c (mmol/mol)",
    title = "Predicted vs Observed (recalibration-in-the-large)",
    subtitle = paste0("n = ", format(sum(!is.na(analysis_dataset$posthba1cfinal)), big.mark = ","), " initiations")
  )

### Recalibration predictions
plot_pred_obs_recal <- analysis_dataset %>%
  mutate(
    pred_drug = case_when(
      drugclass == "SGLT2" ~ pred.recal.SGLT2,
      drugclass == "GLP1"  ~ pred.recal.GLP1,
      drugclass == "DPP4"  ~ pred.recal.DPP4,
      drugclass == "SU"    ~ pred.recal.SU,
      drugclass == "TZD"   ~ pred.recal.TZD,
      TRUE ~ NA_real_
    )
  ) %>%
  select(pred_drug, posthba1cfinal) %>%
  ggplot(aes(x = pred_drug, y = posthba1cfinal)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = "green") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_classic() +
  labs(
    x = "Predicted 12-month HbA1c (mmol/mol)",
    y = "Observed 12-month HbA1c (mmol/mol)",
    title = "Predicted vs Observed (recalibration)",
    subtitle = paste0("n = ", format(sum(!is.na(analysis_dataset$posthba1cfinal)), big.mark = ","), " initiations")
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
## Analysis includes both rank 1 (best) and rank 2 (second best) to visualize treatment alternatives


# Original predictions - Find rank 1 and rank 2 optimal drugs
tolerance3_original <- analysis_dataset %>%
  cateval::get_best_drugs(
    tolerance = 3,
    column_names = paste0("pred.", drugs),
    final_var_name = "pred."
  ) %>%
  select(matches("within_3"))


rank1_original <- analysis_dataset %>%
  cateval::get_best_drugs(
    rank = 1,
    column_names = paste0("pred.", drugs),
    final_var_name = "pred."
  ) %>%
  select(matches("rank1"))

rank2_original <- analysis_dataset %>%
  cateval::get_best_drugs(
    rank = 2,
    column_names = paste0("pred.", drugs),
    final_var_name = "pred."
  ) %>%
  select(matches("rank2"))

rank1_rank2_original <- cbind(rank1_original, rank2_original)

# Calculate difference between rank 1 and rank 2
rank1_rank2_original_filtered <- analysis_dataset %>%
  cbind(rank1_rank2_original) %>%
  mutate(
    rank1_benefit = pred.rank1_drug_value,
    rank2_benefit = pred.rank2_drug_value,
    benefit_difference = rank2_benefit - rank1_benefit
  )

rank1_original_breakdown <- rank1_rank2_original_filtered %>%
  mutate(drug = pred.rank1_drug_name) %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total_n = sum(n),
    proportion = n / total_n * 100
  ) %>%
  arrange(desc(n))

plot_rank1_original <- rank1_original_breakdown %>%
  mutate(drug = factor(drug, levels = c("SGLT2", "GLP1", "DPP4", "SU", "TZD"))) %>%
  ggplot(aes(x = reorder(drug, -n), y = proportion, fill = drug)) +
  geom_col(color = "black", size = 0.3) +
  geom_text(aes(label = paste0(n, "\n(", sprintf("%.1f", proportion), "%)")), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(
    values = c("SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "DPP4" = "#0072B2", 
               "SU" = "#CC79A7", "TZD" = "#D55E00"),
    guide = "none"
  ) +
  scale_x_discrete(labels = c("SGLT2" = "SGLT2i", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", 
                              "SU" = "SU", "TZD" = "TZD")) +
  scale_y_continuous(limits = c(0, max(rank1_original_breakdown$proportion) * 1.15)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Drug Class",
    y = "Proportion of Cohort (%)",
    title = "Optimal Drug Distribution - Rank 1 (Original Predictions)",
    subtitle = "Full cohort"
  )


# Recalibration-in-the-large predictions
tolerance3_recal_large <- analysis_dataset %>%
  cateval::get_best_drugs(
    tolerance = 3,
    column_names = paste0("pred.recal_large.", drugs),
    final_var_name = "pred.recal_large."
  ) %>%
  select(matches("within_3"))

rank1_recal_large <- analysis_dataset %>%
  cateval::get_best_drugs(
    rank = 1,
    column_names = paste0("pred.recal_large.", drugs),
    final_var_name = "pred.recal_large."
  ) %>%
  select(matches("rank1"))

rank2_recal_large <- analysis_dataset %>%
  cateval::get_best_drugs(
    rank = 2,
    column_names = paste0("pred.recal_large.", drugs),
    final_var_name = "pred.recal_large."
  ) %>%
  select(matches("rank2"))

rank1_rank2_recal_large <- cbind(rank1_recal_large, rank2_recal_large)

rank1_rank2_recal_large_filtered <- analysis_dataset %>%
  cbind(rank1_rank2_recal_large) %>%
  mutate(
    rank1_benefit = pred.recal_large.rank1_drug_value,
    rank2_benefit = pred.recal_large.rank2_drug_value,
    benefit_difference = rank2_benefit - rank1_benefit
  )

rank1_recal_large_breakdown <- rank1_rank2_recal_large_filtered %>%
  mutate(drug = pred.recal_large.rank1_drug_name) %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total_n = sum(n),
    proportion = n / total_n * 100
  ) %>%
  arrange(desc(n))

plot_rank1_recal_large <- rank1_recal_large_breakdown %>%
  mutate(drug = factor(drug, levels = c("SGLT2", "GLP1", "DPP4", "SU", "TZD"))) %>%
  ggplot(aes(x = reorder(drug, -n), y = proportion, fill = drug)) +
  geom_col(color = "black", size = 0.3) +
  geom_text(aes(label = paste0(n, "\n(", sprintf("%.1f", proportion), "%)")), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(
    values = c("SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "DPP4" = "#0072B2", 
               "SU" = "#CC79A7", "TZD" = "#D55E00"),
    guide = "none"
  ) +
  scale_x_discrete(labels = c("SGLT2" = "SGLT2i", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", 
                              "SU" = "SU", "TZD" = "TZD")) +
  scale_y_continuous(limits = c(0, max(rank1_recal_large_breakdown$proportion) * 1.15)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Drug Class",
    y = "Proportion of Cohort (%)",
    title = "Optimal Drug Distribution - Rank 1 (Recalibration-in-the-large)",
    subtitle = "Full cohort"
  )

## comparing tolerance3 combinations original vs intercept (drugs added and removed)
compare_tolerance3_orig_recal_large <- tolerance3_original %>%
  select(original = pred.within_3_of_best_drug_name) %>%
  cbind(
    tolerance3_recal_large %>%
      select(intercept = pred.recal_large.within_3_of_best_drug_name)
  ) %>%
  mutate(
    original_list  = map(original, ~ str_split(ifelse(is.na(.x), "", .x), ",")[[1]] %>% str_trim()),
    intercept_list = map(intercept, ~ str_split(ifelse(is.na(.x), "", .x), ",")[[1]] %>% str_trim()),
    added   = map2_chr(intercept_list, original_list,
                       ~ paste(setdiff(.x, .y), collapse = ",")),
    removed = map2_chr(original_list, intercept_list,
                       ~ paste(setdiff(.x, .y), collapse = ","))
  ) %>%
  select(added, removed) %>%
  pivot_longer(everything(), names_to = "change_type", values_to = "drugs") %>%
  separate_rows(drugs, sep = ",") %>%
  mutate(drugs = str_trim(drugs)) %>%
  filter(drugs != "") %>%
  count(change_type, drugs) %>%
  mutate(
    percent = 100 * n / nrow(analysis_dataset)
  )

plot_tolerance3_added_removed_orig_recal_large <- compare_tolerance3_orig_recal_large %>%
  complete(change_type = c("added", "removed"), drugs = drugs, fill = list(n = 0, percent = 0)) %>%
  mutate(
    drugs = factor(drugs, levels = unique(drugs)),
    signed_n = ifelse(change_type == "removed", -n, n)
  ) %>%
  ggplot(aes(x = drugs, y = signed_n, fill = drugs, alpha = change_type)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = c("SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "DPP4" = "#0072B2",
               "SU" = "#CC79A7", "TZD" = "#D55E00"),
    guide = "none"
  ) +
  scale_alpha_manual(values = c("added" = 1, "removed" = 0.6)) +
  scale_x_discrete(labels = c("SGLT2" = "SGLT2i", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i",
                              "SU" = "SU", "TZD" = "TZD")) +
  scale_y_continuous(
    name = "Net count (n)",
    breaks = function(x) {
      lower <- floor(min(x) / 5000) * 5000
      upper <- ceiling(max(x) / 5000) * 5000
      seq(lower, upper, by = 5000)
    },
    sec.axis = sec_axis(
      ~ . * 100 / nrow(analysis_dataset),
      name = "Net percentage (%)",
      labels = function(x) sprintf("%.1f", x),
      breaks = seq(-50, 50, by = 5)
    )
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  labs(
    x = "Drug Class",
    title = "Tolerance-3 drug changes: Original vs Recalibration-in-the-large",
    alpha = "Change"
  )



# Recalibration predictions
tolerance3_recal <- analysis_dataset %>%
  cateval::get_best_drugs(
    tolerance = 3,
    column_names = paste0("pred.recal.", drugs),
    final_var_name = "pred.recal."
  ) %>%
  select(matches("within_3"))

rank1_recal <- analysis_dataset %>%
  cateval::get_best_drugs(
    rank = 1,
    column_names = paste0("pred.recal.", drugs),
    final_var_name = "pred.recal."
  ) %>%
  select(matches("rank1"))

rank2_recal <- analysis_dataset %>%
  cateval::get_best_drugs(
    rank = 2,
    column_names = paste0("pred.recal.", drugs),
    final_var_name = "pred.recal."
  ) %>%
  select(matches("rank2"))

rank1_rank2_recal <- cbind(rank1_recal, rank2_recal)

rank1_rank2_recal_filtered <- analysis_dataset %>%
  cbind(rank1_rank2_recal) %>%
  mutate(
    rank1_benefit = pred.recal.rank1_drug_value,
    rank2_benefit = pred.recal.rank2_drug_value,
    benefit_difference = rank2_benefit - rank1_benefit
  )

rank1_recal_breakdown <- rank1_rank2_recal_filtered %>%
  mutate(drug = pred.recal.rank1_drug_name) %>%
  group_by(drug) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total_n = sum(n),
    proportion = n / total_n * 100
  ) %>%
  arrange(desc(n))

plot_rank1_recal <- rank1_recal_breakdown %>%
  mutate(drug = factor(drug, levels = c("SGLT2", "GLP1", "DPP4", "SU", "TZD"))) %>%
  ggplot(aes(x = reorder(drug, -n), y = proportion, fill = drug)) +
  geom_col(color = "black", size = 0.3) +
  geom_text(aes(label = paste0(n, "\n(", sprintf("%.1f", proportion), "%)")), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(
    values = c("SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "DPP4" = "#0072B2", 
               "SU" = "#CC79A7", "TZD" = "#D55E00"),
    guide = "none"
  ) +
  scale_x_discrete(labels = c("SGLT2" = "SGLT2i", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", 
                              "SU" = "SU", "TZD" = "TZD")) +
  scale_y_continuous(limits = c(0, max(rank1_recal_breakdown$proportion) * 1.15)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Drug Class",
    y = "Proportion of Cohort (%)",
    title = "Optimal Drug Distribution - Rank 1 (Recalibration)",
    subtitle = "Full cohort"
  )

## comparing tolerance3 combinations original vs intercept + slope (drugs added and removed)
compare_tolerance3_orig_recal <- tolerance3_original %>%
  select(original = pred.within_3_of_best_drug_name) %>%
  cbind(
    tolerance3_recal %>%
      select(recal = pred.recal.within_3_of_best_drug_name)
  ) %>%
  mutate(
    original_list = map(original, ~ str_split(ifelse(is.na(.x), "", .x), ",")[[1]] %>% str_trim()),
    recal_list = map(recal, ~ str_split(ifelse(is.na(.x), "", .x), ",")[[1]] %>% str_trim()),
    added = map2_chr(recal_list, original_list,
                     ~ paste(setdiff(.x, .y), collapse = ",")),
    removed = map2_chr(original_list, recal_list,
                       ~ paste(setdiff(.x, .y), collapse = ","))
  ) %>%
  select(added, removed) %>%
  pivot_longer(everything(), names_to = "change_type", values_to = "drugs") %>%
  separate_rows(drugs, sep = ",") %>%
  mutate(drugs = str_trim(drugs)) %>%
  filter(drugs != "") %>%
  count(change_type, drugs) %>%
  mutate(
    percent = 100 * n / nrow(analysis_dataset)
  )

plot_tolerance3_added_removed_orig_recal <- compare_tolerance3_orig_recal %>%
  complete(change_type = c("added", "removed"), drugs = drugs, fill = list(n = 0, percent = 0)) %>%
  mutate(
    drugs = factor(drugs, levels = unique(drugs)),
    signed_n = ifelse(change_type == "removed", -n, n)
  ) %>%
  ggplot(aes(x = drugs, y = signed_n, fill = drugs, alpha = change_type)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = c("SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "DPP4" = "#0072B2",
               "SU" = "#CC79A7", "TZD" = "#D55E00"),
    guide = "none"
  ) +
  scale_alpha_manual(values = c("added" = 1, "removed" = 0.6)) +
  scale_x_discrete(labels = c("SGLT2" = "SGLT2i", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i",
                              "SU" = "SU", "TZD" = "TZD")) +
  scale_y_continuous(
    name = "Net count (n)",
    breaks = function(x) {
      lower <- floor(min(x) / 5000) * 5000
      upper <- ceiling(max(x) / 5000) * 5000
      seq(lower, upper, by = 5000)
    },
    sec.axis = sec_axis(
      ~ . * 100 / nrow(analysis_dataset),
      name = "Net percentage (%)",
      labels = function(x) sprintf("%.1f", x),
      breaks = seq(-50, 50, by = 5)
    )
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  labs(
    x = "Drug Class",
    title = "Tolerance-3 drug changes: Original vs Recalibration",
    alpha = "Change"
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


## STEP 1: Original ----

ggsave(paste0(output_path, "/step1_original_pairwise.pdf"),
       plot_pairwise_original, width = 12, height = 8)
ggsave(paste0(output_path, "/step1_original_overall.pdf"),
       plot_overall_original, width = 18, height = 6)
ggsave(paste0(output_path, "/step1_original_overall_by_drug.pdf"),
       plot_overall_original_bydrug, width = 15, height = 10)
ggsave(paste0(output_path, "/step1_original_pred_obs.pdf"),
       plot_pred_obs_original, width = 10, height = 8)

# balance (love) plots for the two matching directions (one file each)
pdf(paste0(output_path, "/step1_original_loveplot_concordant.pdf"), width = 8, height = 7)
plot(summary(match_original$matchit$concordant), main = "Concordant as treated")
dev.off()
pdf(paste0(output_path, "/step1_original_loveplot_discordant.pdf"), width = 8, height = 7)
plot(summary(match_original$matchit$discordant), main = "Discordant as treated")
dev.off()


## STEP 2: Recalibration-in-the-large ----

ggsave(paste0(output_path, "/step2_recalibration_in_the_large_pairwise.pdf"),
       plot_pairwise_recal_large, width = 12, height = 8)
ggsave(paste0(output_path, "/step2_recalibration_in_the_large_overall.pdf"),
       plot_overall_recal_large, width = 18, height = 6)
ggsave(paste0(output_path, "/step2_recalibration_in_the_large_overall_by_drug.pdf"),
       plot_overall_recal_large_bydrug, width = 15, height = 10)
ggsave(paste0(output_path, "/step2_recalibration_in_the_large_pred_obs.pdf"),
       plot_pred_obs_recal_large, width = 10, height = 8)

# balance (love) plots for the two matching directions (one file each)
pdf(paste0(output_path, "/step2_recalibration_in_the_large_loveplot_concordant.pdf"), width = 8, height = 7)
plot(summary(match_recal_large$matchit$concordant), main = "Concordant as treated")
dev.off()
pdf(paste0(output_path, "/step2_recalibration_in_the_large_loveplot_discordant.pdf"), width = 8, height = 7)
plot(summary(match_recal_large$matchit$discordant), main = "Discordant as treated")
dev.off()


## STEP 3: Recalibration ----

ggsave(paste0(output_path, "/step3_recalibration_pairwise.pdf"),
       plot_pairwise_recal, width = 12, height = 8)
ggsave(paste0(output_path, "/step3_recalibration_overall.pdf"),
       plot_overall_recal, width = 18, height = 6)
ggsave(paste0(output_path, "/step3_recalibration_overall_by_drug.pdf"),
       plot_overall_recal_bydrug, width = 15, height = 10)
ggsave(paste0(output_path, "/step3_recalibration_pred_obs.pdf"),
       plot_pred_obs_recal, width = 10, height = 8)

# balance (love) plots for the two matching directions (one file each)
pdf(paste0(output_path, "/step3_recalibration_loveplot_concordant.pdf"), width = 8, height = 7)
plot(summary(match_recal$matchit$concordant), main = "Concordant as treated")
dev.off()
pdf(paste0(output_path, "/step3_recalibration_loveplot_discordant.pdf"), width = 8, height = 7)
plot(summary(match_recal$matchit$discordant), main = "Discordant as treated")
dev.off()


## Additional Analyses: Clinical Decision Analysis ----

ggsave(paste0(output_path, "/clinical_decision_rank1_original_predictions.pdf"), 
       plot_rank1_original, width = 8, height = 5)

ggsave(paste0(output_path, "/clinical_decision_rank1_recal_large_recalibrated.pdf"), 
       plot_rank1_recal_large, width = 8, height = 5)

ggsave(paste0(output_path, "/clinical_decision_tolerance3_added_removed_original_vs_recal_large.pdf"),
       plot_tolerance3_added_removed_orig_recal_large, width = 9, height = 6)

ggsave(paste0(output_path, "/clinical_decision_rank1_recal_recalibrated.pdf"), 
       plot_rank1_recal, width = 8, height = 5)

ggsave(paste0(output_path, "/clinical_decision_tolerance3_added_removed_original_vs_recal.pdf"),
       plot_tolerance3_added_removed_orig_recal, width = 9, height = 6)

ggsave(paste0(output_path, "/overall_benefit_estimates_by_model_version.pdf"),
       plot_benefit_estimates, width = 9, height = 5)

# table of the average overall benefit per step (estimate and 95% CI)
write.csv(
  benefit_estimates %>%
    select(model, estimate = observed, lower_ci = lower.ci, upper_ci = upper.ci, n_pairs = n),
  paste0(output_path, "/overall_benefit_estimates_by_model_version.csv"),
  row.names = FALSE
)



# Validation Script for the Exeter 5-Drug Model
#
## Validation and update of the post-2020 updated model for oral semaglutide:
##  CPRD initiations of the model drug classes from 2020-10-14 onwards, and all
##  semaglutide initiations (any date)
##
## The model has no semaglutide arm, so oral semaglutide is added as a new drug
##  class ("Oral_Semaglutide") next to the model drug classes, where GLP1 is split into
##  its subtypes (glp1_subtypes): dulaglutide, and liraglutide (with exenatide).
##  Injectable semaglutide is left out; it is updated in
##  02.calibration_sema_injectable.R. Its starting prediction is the published model's
##  GLP1 prediction.
##
## The model validated here is the published model plus the post-2020 update
##  estimated in 01.Temporal_Update/02.calibration_new_data_post2020.R (read
##  from Model_files), which has its own intercept for each GLP1 subtype.
##  "Original" below therefore means the updated model for the model drug
##  classes, with pred.Oral_Semaglutide = the published pred.GLP1.
##  Recalibration-in-the-large and recalibration then estimate an intercept (and
##  slope) for the new semaglutide class only: the model drug classes are kept
##  fixed at the post-2020 update (no further adjustment), so their recalibrated
##  predictions equal the original ones. The Oral_Semaglutide intercept, on top of the
##  published GLP1 prediction, is the oral semaglutide update, saved to
##  03.Semaglutide_Update/Data with its difference from each GLP1 subtype's
##  post-2020 update. Each step's semaglutide estimate is also compared, with its
##  95% CI, against the published GLP1 prediction and each GLP1 subtype.
##
## The validation (Steps 3 and 4) compares semaglutide with the drug classes
##  that are not GLP1 (comparator_drugs: SGLT2, DPP4, SU, TZD); the predicted best
##  drug is chosen from these 5. The GLP1 subtypes, built on the same published
##  GLP1 prediction as semaglutide, are compared with it only as one pooled
##  pairwise estimate each (cal_groups = 1). Predicted vs observed uses all
##  initiations, the GLP1 subtypes included.
##
## Drug class names follow the prediction columns of the model calculator
##  (pred.Oral_Semaglutide, pred.Inje_Semaglutide).
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
#    a) Pairwise drug calibration (and semaglutide vs each GLP1 subtype)
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
#     Valid values: "SGLT2", "Dulaglutide", "Liraglutide", "DPP4", "SU", or "TZD"
#     (GLP1 split into its subtypes, glp1_subtypes below), plus the new drug
#     class validated in this script (new_drug, below)
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

# First initiation date of the model drug classes: two years after the
# development cohort closed (2018-10-13), so no initiation here was used to fit
# the model, and the same cohort the post-2020 update was estimated in.
# Semaglutide was never in the development data, so all its initiations are
# used, whatever their date.
cohort_start <- as.Date("2020-10-14")

# Semaglutide forms, each its own drug class, and the drug_substance values that
# make them up (grouped as in 01.Temporal_Update/01.tx_response_year.R)
sema_classes <- list(
  "Oral_Semaglutide" = c("Oral semaglutide"),
  "Inje_Semaglutide" = c("Low-dose semaglutide", "High-dose semaglutide", "Semaglutide, dose unclear")
)

# glp1_subtypes: the model has one GLP1 class, split here into subtypes by
#   drug_substance. Each subtype starts from the model's GLP1 prediction plus its
#   own post-2020 update, and is not recalibrated here. Exenatide is grouped with
#   liraglutide.
glp1_subtypes <- list(
  "Dulaglutide" = c("Dulaglutide"),
  "Liraglutide" = c("Liraglutide", "Exenatide", "Exenatide prolonged-release")
)

# The new drug class validated and updated in this script, and the published
# model's drug class whose prediction it starts from
new_drug  <- "Oral_Semaglutide"
base_drug <- "GLP1"

# Drug classes of the cohort: the new one first (so it is drug 1 in the pairwise
# comparisons), then the model drug classes, with GLP1 split into its subtypes.
# All of them are used in the recalibration and predicted vs observed.
drugs <- c(new_drug, "SGLT2", names(glp1_subtypes), "DPP4", "SU", "TZD")
model_drugs <- drugs

# comparator_drugs: the drug classes semaglutide is validated against in Steps 3
#   and 4 (pairwise, overall benefit, predicted best drug). The GLP1 subtypes are
#   left out: they are built on the same published GLP1 prediction as semaglutide.
comparator_drugs <- c("SGLT2", "DPP4", "SU", "TZD")

# glp1_panel_min_n: initiations needed in both semaglutide and a GLP1 subtype for
#   their pooled pairwise estimate (cal_groups = 1) to be shown
glp1_panel_min_n <- 100

# all semaglutide substances, used to keep their initiations from any date
sema_substances <- unlist(sema_classes, use.names = FALSE)

# GLP1 subtype of each GLP1 drug_substance (glp1_subtypes, Step 0)
glp1_lookup <- setNames(rep(names(glp1_subtypes), lengths(glp1_subtypes)),
                        unlist(glp1_subtypes, use.names = FALSE))

# set up dataset (cohort cached by 01.cohort_scripts.R)
analysis_raw <- analysis_raw %>%
  analysis$cached("ethnicity") %>%
  filter(dstartdate >= cohort_start | drug_substance %in% sema_substances) %>%
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
  stop("GLP1 substances in no GLP1 subtype (glp1_subtypes, Step 0): ",
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

    # Drug class: semaglutide is split out of GLP1 into its own classes, and the
    # other GLP1 substances into the GLP1 subtypes
    drugclass = case_when(
      drug_substance %in% sema_classes$Oral_Semaglutide ~ "Oral_Semaglutide",
      drug_substance %in% sema_classes$Inje_Semaglutide ~ "Inje_Semaglutide",
      drug_class == "GLP1" ~ unname(glp1_lookup[drug_substance]),
      TRUE ~ drug_class
    ),
    drugclass = factor(drugclass,
                       levels = c("SGLT2", names(glp1_subtypes), "DPP4", "SU", "TZD", names(sema_classes)))
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
  # the model drug classes and both semaglutide forms; the cohort is
  # restricted to this script's drug classes (drugs) after imputation
  filter(!is.na(drugclass))

# initiations and date range per drug class (semaglutide may start before cohort_start)
analysis_dataset %>%
  group_by(drugclass) %>%
  summarise(n = n(), first = min(dstartdate), last = max(dstartdate))


# Step 1b: Create a characteristics table ---------------------------------

# This table describes the validation cohort at treatment initiation, stratified by
# initiated drug class (the new drug class and the model drug classes, with the
# GLP1 subtypes combined into GLP1).
# It includes baseline clinical characteristics, missingness indicators, outcome timing,
# observed 12-month HbA1c, and observed response (posthba1cfinal - prehba1c).


# Add observed response and missingness indicators
analysis_dataset_drugclass <- analysis_dataset %>%
  filter(drugclass %in% drugs) %>% # this script's drug classes only
  mutate(
    drugclass = as.character(drugclass),
    drugclass = ifelse(drugclass %in% names(glp1_subtypes), "GLP1", drugclass),
    drugclass = factor(drugclass, levels = c(new_drug, "SGLT2", "GLP1", "DPP4", "SU", "TZD")),
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

output_path <- "Scripts/03.Semaglutide_Update/Outputs/Oral_semaglutide"
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

# The cohort is restricted to this script's drug classes only now, after
# imputation on the full cohort with both semaglutide forms (so
# imputed values are the same in the oral and injectable semaglutide scripts),
# and initiations with a missing outcome HbA1c are dropped
analysis_dataset <- analysis_dataset %>%
  filter(drugclass %in% drugs) %>%
  mutate(drugclass = factor(drugclass, levels = drugs)) %>%
  filter(!is.na(posthba1cfinal))

cat(new_drug, "cohort:", nrow(analysis_dataset), "initiations\n")
table(analysis_dataset$drugclass)


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

# The new drug class starts from the published model's prediction of its base
# class (pred.GLP1, which the update leaves as it is)
analysis_dataset[[paste0("pred.", new_drug)]] <- analysis_dataset[[paste0("pred.", base_drug)]]


## We then generate versions of the predictions that are recalibrated using intercept-only
## and intercept + slope adjustments.

### Recalibration
# drugs and model_drugs (new_drug plus the model drug classes, with the GLP1
# subtypes) are set at the top of the script
outcome_var <- "posthba1cfinal"
pred_prefix <- "pred."

# drug class names as shown in the plots
drug_labels <- c("SGLT2" = "SGLT2i", "Dulaglutide" = "Dulaglutide", "Liraglutide" = "Liraglutide",
                 "DPP4" = "DPP4i", "SU" = "SU", "TZD" = "TZD",
                 "Oral_Semaglutide" = "Oral Semaglutide", "Inje_Semaglutide" = "S.c. Semaglutide")

# drug class colours (Okabe-Ito), as in the online calculator: the GLP1 subtypes
# and semaglutide take the GLP-1RA colour
drug_colours <- c("SGLT2" = "#E69F00", "Dulaglutide" = "#56B4E9", "Liraglutide" = "#56B4E9",
                  "DPP4" = "#0072B2", "SU" = "#CC79A7", "TZD" = "#D55E00",
                  "Oral_Semaglutide" = "#56B4E9", "Inje_Semaglutide" = "#56B4E9")

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

# Recalibration is done with a single pooled model per method. Only new_drug is
# adjusted (its own intercept, and slope): the model drug classes are already
# fixed by the post-2020 update, so they enter through the offset only (intercept
# 0, slope 1) and their recalibrated predictions equal the original ones. Any
# prognostic_vars enter as a single shared main effect, so their influence on the
# outcome is the same regardless of which treatment a patient received. With
# prognostic_vars empty these pooled fits reproduce a fit in new_drug alone.

# Anchor each patient on the prediction for the drug they actually received; this
# is the offset (intercept model) / covariate (intercept + slope model).
recal_data <- analysis_dataset
recal_data$drugclass <- factor(recal_data$drugclass, levels = drugs)
pred_matrix <- as.matrix(analysis_dataset[, paste0(pred_prefix, drugs)])
recal_data$pred_assigned <- pred_matrix[
  cbind(seq_len(nrow(pred_matrix)), match(recal_data$drugclass, drugs))
]

# new_drug design columns: indicator (intercept), its prediction (slope), and the
# prediction of the fixed classes (offset in the intercept + slope model)
add_new_drug_terms <- function(data) {
  data$new_ind    <- as.numeric(data$drugclass == new_drug)
  data$new_pred   <- data$new_ind * data$pred_assigned
  data$fixed_pred <- (1 - data$new_ind) * data$pred_assigned
  data
}
recal_data <- add_new_drug_terms(recal_data)

# helper: name the new_drug terms as drug-specific terms (drugclass<new_drug>,
# drugclass<new_drug>:pred_assigned), with n = patients on new_drug
label_new_drug_terms <- function(params) {
  n_new <- sum(recal_data$drugclass == new_drug)
  params %>%
    mutate(
      n    = ifelse(term %in% c("new_ind", "new_pred"), n_new, n),
      term = recode(term,
                    new_ind  = paste0("drugclass", new_drug),
                    new_pred = paste0("drugclass", new_drug, ":pred_assigned"))
    )
}

# helper: predict the recalibrated value for every patient under a given drug
predict_for_drug <- function(fit, drug) {
  nd <- analysis_dataset
  nd$drugclass <- factor(drug, levels = drugs)
  nd$pred_assigned <- nd[[paste0(pred_prefix, drug)]]
  nd <- add_new_drug_terms(nd)
  as.numeric(predict(fit, newdata = nd))
}

# helper: each patient's prediction (columns <prefix><drug>) for the drug they
# actually received
pred_received <- function(data, prefix) {
  pm <- as.matrix(data[, paste0(prefix, drugs)])
  pm[cbind(seq_len(nrow(pm)), match(data$drugclass, drugs))]
}

#### Recalibration-in-the-large (new_drug intercept, shared prognostic effect)
formula_recal_large <- as.formula(
  paste0(outcome_var, " ~ 0 + new_ind", prognostic_terms, " + offset(pred_assigned)")
)
fit_recal_large <- lm(formula_recal_large, data = recal_data)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal_large.", drug)]] <- predict_for_drug(fit_recal_large, drug)
}
recalib_recal_large_params <- label_new_drug_terms(summarize_lm(fit_recal_large, "recal_large"))

#### Recalibration (new_drug intercept & slope, shared prognostic effect)
formula_recal <- as.formula(
  paste0(outcome_var, " ~ 0 + new_ind + new_pred", prognostic_terms, " + offset(fixed_pred)")
)
fit_recal <- lm(formula_recal, data = recal_data)
for (drug in drugs) {
  analysis_dataset[[paste0("pred.recal.", drug)]] <- predict_for_drug(fit_recal, drug)
}
recalib_recal_params <- label_new_drug_terms(summarize_lm(fit_recal, "recal"))

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

# Compare new_drug's prediction at each step with the published GLP1 prediction
# (v1 model) and with each GLP1 subtype's prediction (post-2020 update): the
# mean difference over the new_drug initiations, with its 95% CI. With an
# intercept only this is the intercept itself (minus the subtype's intercept);
# with a slope it is intercept + (slope - 1) x the mean published GLP1 prediction,
# which least squares makes equal to the intercept-only estimate, so the
# recalibration step also reports the slope of new_drug's prediction on the
# published GLP1 prediction (1 for the published GLP1 prediction and the GLP1
# subtypes, so the same against every reference).
# The original step is the published GLP1 prediction, so it has no CI. The GLP1
# subtype intercepts are fixed, so they shift the difference and its CI. Shared
# prognostic factors (prognostic_vars) move every drug class equally and are
# left out.
compare_new_drug <- function(fit, step) {
  n_new <- sum(recal_data$new_ind == 1)
  references <- c(paste0("Published ", base_drug, " (v1 model)"), names(glp1_subtypes))
  reference_intercept <- c(0, model_update$intercept[match(names(glp1_subtypes), model_update$drugclass)])
  if (is.null(fit)) {
    est <- 0; ci <- c(NA_real_, NA_real_)
  } else {
    b <- coef(fit)
    L <- setNames(numeric(length(b)), names(b))
    L["new_ind"] <- 1
    shift <- 0
    if ("new_pred" %in% names(b)) {
      mean_pred <- mean(recal_data$pred_assigned[recal_data$new_ind == 1])
      L["new_pred"] <- mean_pred
      shift <- -mean_pred
    }
    est <- sum(L * b) + shift
    se  <- sqrt(drop(t(L) %*% vcov(fit) %*% L))
    ci  <- est + c(-1, 1) * qt(0.975, df.residual(fit)) * se
  }
  out <- tibble(
    recalibration = step,
    reference     = references,
    term          = "mean difference",
    estimate      = est - reference_intercept,
    ci_low        = ci[1] - reference_intercept,
    ci_high       = ci[2] - reference_intercept,
    n             = n_new
  )
  if (!is.null(fit) && "new_pred" %in% names(coef(fit))) {
    slope_ci <- confint(fit)["new_pred", ]
    out <- bind_rows(out, tibble(
      recalibration = step,
      reference     = references,
      term          = "slope",
      estimate      = unname(coef(fit)["new_pred"]),
      ci_low        = unname(slope_ci[1]),
      ci_high       = unname(slope_ci[2]),
      n             = n_new
    ))
  }
  out
}

new_drug_comparison <- bind_rows(
  compare_new_drug(NULL, "original"),
  compare_new_drug(fit_recal_large, "recal_large"),
  compare_new_drug(fit_recal, "recal")
)
print(new_drug_comparison)

write_csv(new_drug_comparison, paste0(output_path, "/recalibration_", new_drug, "_vs_glp1.csv"))

# Save the recalibration-in-the-large update for the new drug class so it can
# be applied to the model in other analyses. Only new_drug is saved: the
# model drug classes keep the post-2020 update (01.Temporal_Update/Data).
#   intercept:     on top of the published model's prediction of base_drug
#                  updated prediction = pred.<base_drugclass> (published) + intercept
#   intercept_vs_<GLP1 subtype>: for reporting, the same update relative to that
#                  GLP1 subtype's post-2020 update (intercept minus the subtype's
#                  intercept), so
#                  updated prediction = pred.<GLP1 subtype> (post-2020 update) + intercept_vs_<GLP1 subtype>
data_path <- "Scripts/03.Semaglutide_Update/Data"
dir.create(data_path, showWarnings = FALSE, recursive = TRUE)

model_update_new_drug <- recalib_recal_large_params %>%
  filter(term == paste0("drugclass", new_drug)) %>%
  transmute(
    drugclass      = new_drug,
    base_drugclass = base_drug,
    intercept      = round(estimate, 3)
  )

for (subtype in names(glp1_subtypes)) {
  model_update_new_drug[[paste0("intercept_vs_", subtype)]] <- round(
    model_update_new_drug$intercept - model_update$intercept[model_update$drugclass == subtype], 3
  )
}

print(model_update_new_drug)

saveRDS(model_update_new_drug,
        paste0(data_path, "/model_update_recal_large_", new_drug, ".rds"))


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
# ("Oral Semaglutide vs Dulaglutide")
pairwise_panel_size <- 3.4


#############################################################################

# Steps 3 and 4: semaglutide vs the drug classes that are not GLP1 -------------

# The analysis is done with new_drug and comparator_drugs (Step 0): initiations
# of the GLP1 subtypes are left out, and the predicted best drug is chosen from
# these drug classes only. The GLP1 subtypes come in only in the pooled pairwise
# estimates against new_drug. Predicted vs observed does not compare drug
# classes, so it is done with all initiations.

analysis_drugs <- c(new_drug, comparator_drugs)

main <- list(
  label       = paste(paste(drug_labels[head(analysis_drugs, -1)], collapse = ", "),
                      "and", drug_labels[[tail(analysis_drugs, 1)]]),
  data        = analysis_dataset %>%
    filter(drugclass %in% analysis_drugs) %>%
    mutate(drugclass = factor(drugclass, levels = analysis_drugs)),
  drugs       = analysis_drugs,
  model_drugs = analysis_drugs
)

cat("Analysis (", main$label, "): ", nrow(main$data), " initiations\n", sep = "")


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
                     "sex", "smoke", "ncurrtx", "drugline", "imd5")

# - Pairwise figures: the drug classes of the analysis against each other, and
#   new_drug against each of them plus one panel per GLP1 subtype. New_drug and a
#   GLP1 subtype are built on the same published GLP1 prediction, so their
#   predicted difference is (near) constant: each of these panels is a single
#   pooled estimate (cal_groups = 1), shown when both drug classes have at least
#   glp1_panel_min_n initiations (Step 0).

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

# GLP1 subtypes with enough initiations, next to new_drug, for a pooled estimate
glp1_counts <- table(factor(analysis_dataset$drugclass, levels = drugs))
glp1_panel_drugs <- names(glp1_subtypes)[
  glp1_counts[names(glp1_subtypes)] >= glp1_panel_min_n & glp1_counts[[new_drug]] >= glp1_panel_min_n
]
for (subtype in setdiff(names(glp1_subtypes), glp1_panel_drugs)) {
  cat("No pooled pairwise estimate for", new_drug, "vs", subtype, ": fewer than",
      glp1_panel_min_n, "initiations (", new_drug, glp1_counts[[new_drug]], ",",
      subtype, glp1_counts[[subtype]], ")\n")
}

# helper: pooled pairwise calibration (cal_groups = 1) of new_drug against each
# GLP1 subtype in glp1_panel_drugs, with the predictions with the given prefix
pairwise_calibration_glp1 <- function(prefix) {
  map_dfr(glp1_panel_drugs, function(subtype) {
    pair_drugs <- c(new_drug, subtype)
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
  })
}

# helper: pairwise calibration plot of the drug classes in drug1_set (drug 1)
# against those in drug2_set (drug 2), one panel per pair
plot_pairwise <- function(pairwise_results, plot_title, drug1_set, drug2_set = drug1_set) {
  pairwise_results %>%
    filter(drug1 %in% drug1_set, drug2 %in% drug2_set) %>%
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

# helper: the pairwise results of one step: the drug classes of the analysis
# against each other, and the pooled estimates of new_drug vs the GLP1 subtypes
pairwise_calibration_step <- function(prefix) {
  bind_rows(pairwise_calibration_for(main, prefix), pairwise_calibration_glp1(prefix))
}

# helper: the pairwise figures of one step, named by the end of their file
# names: the drug classes of the analysis ("") and new_drug against them and the
# GLP1 subtypes ("_<new_drug>_with_glp1")
plot_pairwise_figures <- function(pairwise_results, step_title) {
  plots <- setNames(list(plot_pairwise(pairwise_results,
    paste0("Pairwise Drug Calibration (", step_title, ")"), main$drugs)), "")
  plots[[paste0("_", new_drug, "_with_glp1")]] <- plot_pairwise(pairwise_results,
    paste0("Pairwise Drug Calibration: ", drug_labels[[new_drug]], " (", step_title, ")"),
    new_drug, c(comparator_drugs, glp1_panel_drugs))
  plots
}

### Original predictions
pairwise_results_original <- pairwise_calibration_step("pred.")
plot_pairwise_original    <- plot_pairwise_figures(pairwise_results_original, "original predictions")


### Recalibration-in-the-large predictions
pairwise_results_recal_large <- pairwise_calibration_step("pred.recal_large.")
plot_pairwise_recal_large    <- plot_pairwise_figures(pairwise_results_recal_large, "recalibration-in-the-large")


### Recalibration predictions
pairwise_results_recal <- pairwise_calibration_step("pred.recal.")
plot_pairwise_recal    <- plot_pairwise_figures(pairwise_results_recal, "recalibration")


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
# function; the best drug is chosen from new_drug and comparator_drugs), with
# a 0.4 SD caliper on baseline HbA1c. These are baseline
# characteristics only: unlike adjustment_vars above, hba1cmonth (time to
# outcome measurement) describes the outcome window rather than the patient at
# baseline, so it is left out of matching.
matching_vars <- c(
  "t2dmduration", "prebmi", "prehba1c", "agetx", "prealt",
  "preegfr", "pretotalcholesterol", "prehdl",
  "sex", "smoke", "ncurrtx", "drugline", "imd5"
)
# The same covariates plus any prognostic factors in prognostic_vars. Used for
# matching in every overall-benefit step (the prognostic factor is a cross-cutting
# control, not a step of its own); equal to matching_vars until you populate
# prognostic_vars.
matching_vars_prog <- unique(c(matching_vars, prognostic_vars))

# Confidence intervals: FALSE uses fast analytical intervals; set to TRUE for
# pair-level bootstrap intervals (slower) that account for matching with
# replacement.
use_bootstrap <- TRUE
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
  "drugline5+"          = "'Number of previous glucose-lowering drugs = >4'",
  "imd51 (least)"       = "'Index of multiple deprivation = 1 (least deprived)'",
  "imd52"               = "'Index of multiple deprivation = 2'",
  "imd53"               = "'Index of multiple deprivation = 3'",
  "imd54"               = "'Index of multiple deprivation = 4'",
  "imd55 (most)"        = "'Index of multiple deprivation = 5 (most deprived)'"
)

# Balance (love) plot for one matching step: absolute standardised mean
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
    plot         = plot_conc + plot_disc + plot_both +
      patchwork::plot_layout(nrow = 1) +
      patchwork::plot_annotation(title = paste0("Overall Benefit Calibration", title_suffix), theme = theme(plot.title = element_text(size = 20))),
    plot_by_drug = plot_calibration_by_drug(by_drug, title = paste0("Overall Benefit Calibration by drug class", title_suffix))
  )
}

# Matching for every step: stricter in the discordant-as-treated direction
# (extra calipers, exact on drug line and number of current drugs, each patient
# reused at most 3 times). The original predictions are those of the updated
# model, so they are matched in the same way as the two recalibration steps.
match_args_disc_recal <- list(
  caliper = c(prehba1c = 0.4, prebmi = 0.5, t2dmduration = 0.5, agetx = 0.5),
  exact = c("drugline", "ncurrtx"),
  reuse.max = 3, m.order = "random")


### Original predictions ----

benefit_original <- overall_benefit(main, "pred.", " (original predictions)",
                                    match_args_disc = match_args_disc_recal)


### Recalibration-in-the-large predictions ----

benefit_recal_large <- overall_benefit(main, "pred.recal_large.", " (recalibration-in-the-large)",
                                       match_args_disc = match_args_disc_recal)


### Recalibration predictions ----

benefit_recal <- overall_benefit(main, "pred.recal.", " (recalibration)",
                                 match_args_disc = match_args_disc_recal)


### Compare the overall benefit estimate across model versions ----

benefit_estimates <- bind_rows(
  benefit_original$overall %>% mutate(model = "Original"),
  benefit_recal_large$overall %>% mutate(model = "Recalibration-in-the-large"),
  benefit_recal$overall %>% mutate(model = "Recalibration")
) %>%
  mutate(model = factor(
    model,
    levels = c("Original", "Recalibration-in-the-large", "Recalibration")
  ))

# helper: plot of the overall benefit estimates
plot_benefit_estimates <- function(estimates, plot_title) {
  estimates %>%
    ggplot(aes(x = model, y = observed, ymin = lower.ci, ymax = upper.ci, colour = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_pointrange(size = 0.7) +
    coord_flip() +
    theme_classic(base_size = 18) +
    theme(legend.position = "none", axis.text = element_text(size = 16)) +
    labs(
      x = NULL,
      y = expression("Observed HbA"[1*c]*" benefit estimate (mmol/mol)"),
      title = plot_title
    )
}


## Predicted vs Observed Calibration ----

## Predicted vs observed HbA1c of every initiation (the GLP1 subtypes included),
## with the prediction for the drug class initiated

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

  rank1_breakdown %>%
    mutate(drug = factor(drug, levels = analysis$model_drugs)) %>%
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
      subtitle = paste0(step_title, "; full cohort, ", analysis$label)
    )
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
  split_drugs <- function(x) map(x, ~ str_split(ifelse(is.na(.x), "", .x), ",")[[1]] %>% str_trim())
  original_list <- split_drugs(tolerance3_sets(analysis, "pred."))
  step_list     <- split_drugs(tolerance3_sets(analysis, prefix))

  tibble(
    added   = map2_chr(step_list, original_list, ~ paste(setdiff(.x, .y), collapse = ",")),
    removed = map2_chr(original_list, step_list, ~ paste(setdiff(.x, .y), collapse = ","))
  ) %>%
    pivot_longer(everything(), names_to = "change_type", values_to = "drugs") %>%
    separate_rows(drugs, sep = ",") %>%
    mutate(drugs = str_trim(drugs)) %>%
    filter(drugs != "") %>%
    count(change_type, drugs) %>%
    mutate(
      percent = 100 * n / nrow(analysis$data)
    ) %>%
    complete(change_type = c("added", "removed"), drugs, fill = list(n = 0, percent = 0)) %>%
    mutate(
      drugs = factor(drugs, levels = unique(drugs)),
      signed_n = ifelse(change_type == "removed", -n, n)
    ) %>%
    ggplot(aes(x = drugs, y = signed_n, fill = drugs, alpha = change_type)) +
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
      title = "Tolerance-3 drug changes",
      subtitle = paste0("Original vs ", step_title, "; ", analysis$label),
      alpha = "Change"
    )
}

# Only new_drug's prediction changes between steps (the other drug classes are
# fixed at the post-2020 update), so the tolerance-3 changes are mostly new_drug
# entering or leaving a patient's set; another drug class changes only when
# new_drug becomes or stops being the best drug, moving the 3 mmol/mol threshold.

### Original predictions
plot_rank1_original <- plot_rank1(main, "pred.", "Original predictions")

### Recalibration-in-the-large predictions
plot_rank1_recal_large      <- plot_rank1(main, "pred.recal_large.", "Recalibration-in-the-large")
plot_tolerance3_recal_large <- plot_tolerance3_added_removed(main, "pred.recal_large.", "Recalibration-in-the-large")

### Recalibration predictions
plot_rank1_recal      <- plot_rank1(main, "pred.recal.", "Recalibration")
plot_tolerance3_recal <- plot_tolerance3_added_removed(main, "pred.recal.", "Recalibration")


## Best drug proportions and concordance ----

## For each set of predictions: how often each drug class is the predicted best
## (rank 1), and how often the drug initiated was that drug (concordance)

cohort_label <- unname(drug_labels[new_drug])

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
  summarise_concordance(main, "pred.recal.",       "recal")
)

print(concordance_summary)

write_csv(
  concordance_summary,
  paste0(output_path, "/concordance_by_drug.csv")
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

# helper: save the pairwise and overall benefit figures of one step (file names
# start with file_stem)
save_step_figures <- function(file_stem, plot_pairwise_step, benefit) {
  for (i in seq_along(plot_pairwise_step)) {
    save_plot(plot_pairwise_step[[i]], paste0(output_path, "/", file_stem, "_pairwise", names(plot_pairwise_step)[i], ".pdf"),
              panel_width = pairwise_panel_size, panel_height = pairwise_panel_size)
  }
  save_plot(benefit$plot, paste0(output_path, "/", file_stem, "_overall.pdf"),
            panel_width = 4.2, panel_height = 4.2)
  save_plot(benefit$plot_by_drug, paste0(output_path, "/", file_stem, "_overall_by_drug.pdf"),
            panel_width = 3.4, panel_height = 3.4)

  # balance (love) plot, both matching directions in one figure
  save_plot(plot_loveplot(benefit$matched), paste0(output_path, "/", file_stem, "_loveplot.pdf"),
            panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))
}


## STEP 1: Original ----

save_step_figures("step1_original", plot_pairwise_original, benefit_original)
save_plot(plot_pred_obs_original, paste0(output_path, "/step1_original_pred_obs.pdf"),
          panel_width = 5.5, panel_height = 4.2)


## STEP 2: Recalibration-in-the-large ----

save_step_figures("step2_recalibration_in_the_large", plot_pairwise_recal_large, benefit_recal_large)
save_plot(plot_pred_obs_recal_large, paste0(output_path, "/step2_recalibration_in_the_large_pred_obs.pdf"),
          panel_width = 5.5, panel_height = 4.2)


## STEP 3: Recalibration ----

save_step_figures("step3_recalibration", plot_pairwise_recal, benefit_recal)
save_plot(plot_pred_obs_recal, paste0(output_path, "/step3_recalibration_pred_obs.pdf"),
          panel_width = 5.5, panel_height = 4.2)


## Save the overall-benefit plot components, matched pairs and pairwise results ----
# The individual panels are saved rather than the assembled patchworks, so any
# combination of recalibration step / matched population can be put together
# later without re-running this script. Of the matched pairs only the columns
# needed to regroup them are kept, so groups can be recomputed later without
# repeating any matching. The pairwise results are group-level summaries (the
# pooled new_drug vs GLP1 subtype estimates included), so the drug-pair panels
# can be redrawn later.

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

benefits <- list(original = benefit_original, recal_large = benefit_recal_large, recal = benefit_recal)

saveRDS(
  lapply(benefits, function(b) list(conc = b$plot_conc, disc = b$plot_disc, both = b$plot_both)),
  paste0(output_path, "/overall_benefit_plots.rds")
)
saveRDS(
  lapply(benefits, function(b) list(
    conc = slim_pairs(b$matched$concordant),
    disc = slim_pairs(b$matched$discordant),
    both = slim_pairs(b$matched$combined)
  )),
  paste0(output_path, "/overall_benefit_pairs.rds")
)
saveRDS(
  list(original = pairwise_results_original, recal_large = pairwise_results_recal_large, recal = pairwise_results_recal),
  paste0(output_path, "/pairwise_results.rds")
)


## Additional Analyses: Clinical Decision Analysis ----

# helper: save a clinical decision figure
save_decision_plot <- function(plot, file_stem) {
  save_plot(plot, paste0(output_path, "/", file_stem, ".pdf"),
            panel_width = function(p) max(3, 1.25 * n_axis_levels(p)), panel_height = 3.6)
}

save_decision_plot(plot_rank1_original, "clinical_decision_rank1_original_predictions")

save_decision_plot(plot_rank1_recal_large, "clinical_decision_rank1_recal_large_recalibrated")
save_decision_plot(plot_tolerance3_recal_large, "clinical_decision_tolerance3_added_removed_original_vs_recal_large")

save_decision_plot(plot_rank1_recal, "clinical_decision_rank1_recal_recalibrated")
save_decision_plot(plot_tolerance3_recal, "clinical_decision_tolerance3_added_removed_original_vs_recal")

save_plot(plot_benefit_estimates(benefit_estimates, "Overall benefit estimates"),
          paste0(output_path, "/overall_benefit_estimates_by_model_version.pdf"),
          panel_width = 5.5, panel_height = 2.4)

# table of the average overall benefit per step (estimate and 95% CI)
write.csv(
  benefit_estimates %>%
    select(model, estimate = observed, lower_ci = lower.ci, upper_ci = upper.ci, n_pairs = n),
  paste0(output_path, "/overall_benefit_estimates_by_model_version.csv"),
  row.names = FALSE
)

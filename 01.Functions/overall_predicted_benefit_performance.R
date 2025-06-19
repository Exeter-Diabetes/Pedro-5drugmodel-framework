#' Estimate Observed vs Predicted Benefit Across Calibration Groups
#'
#' This function assesses how well predicted treatment benefits correspond to observed differences in outcomes. Patients are grouped by predicted benefit, and differences in observed outcomes between concordant and discordant treatments are regressed against the predicted differences. Patients are labeled as concordant if they received a treatment predicted to be optimal (either the top predicted or within a specified tolerance of the best). Matching is used to control for confounding between treatment groups.
#'
#' @param data A data frame containing treatment assignments, outcomes, and predicted benefits.
#' @param drug_var Character string. Column name for treatment assignment (actual treatment taken).
#' @param outcome_var Character string. Column name for the outcome variable.
#' @param cal_groups Integer. Number of calibration groups to divide the population into.
#' @param pred_cols Character vector. Column names of predicted outcomes or risks for each drug. These should have a common prefix (e.g., "pred_GLP1", "pred_SGLT2").
#' @param conc_tolerance Optional numeric value. If specified, concordance is defined as receiving a treatment within this absolute difference of the best predicted outcome. If NULL, concordance is based on rank-1 predicted drug.
#' @param matching_var Character vector. Covariates used for Mahalanobis distance matching.
#' @param match.exact Optional character vector. Variables for exact matching. The best predicted drug is automatically included.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching. The actual treatment variable is automatically included.
#'
#' @return A data frame with one row per calibration group, containing:
#' \describe{
#'   \item{mean}{Mean predicted benefit in the group (discordant minus concordant).}
#'   \item{coef}{Estimated observed difference in outcomes (regression coefficient).}
#'   \item{coef_low}{Lower bound of the 95\% confidence interval.}
#'   \item{coef_high}{Upper bound of the 95\% confidence interval.}
#'   \item{n_groups}{Number of calibration groups used.}
#'   \item{n}{Number of matched pairs in each group.}
#' }
#'
#' @details
#' This function uses nearest-neighbor matching (via the `MatchIt` package) on a set of covariates to create matched groups of concordant and discordant patients. Within each calibration group defined by predicted benefit, the observed difference in outcome is regressed on the predicted difference to assess calibration.
#'
#' @import tidyverse
#' @import MatchIt
#'
#' @export
overall_predicted_benefit_performance <- function(data, 
                                      drug_var, 
                                      outcome_var, 
                                      pred_cols = NULL, 
                                      conc_tolerance = NULL,
                                      matching_var = NULL, 
                                      match.exact = NULL, 
                                      match.antiexact = NULL) {
  
  # Load Required Libraries ----
  require(tidyverse)  # For data manipulation and piping
  require(MatchIt)    # For matching procedures
  
  
  # Input Validation ----
  # Check if treatment variable is present
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  # Check if outcome variable is present
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  # Check if all predicted columns are present
  if (!all(pred_cols %in% colnames(data))) stop("Some pred_cols not found in data")
  # Check conc_tolerance is numeric if provided
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("conc_tolerance must be numeric")
  # Check matching variables are present if provided
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  # Check exact matching variables are present if provided
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
  # Check antiexact matching variables are present if provided
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  
  # Prepare Dataset ----
  pre_data <- data %>%
    rename_with(~ "dataset_drug_var", all_of(drug_var)) %>%    # Rename treatment column to a consistent internal name
    rename_with(~ "dataset_outcome_var", all_of(outcome_var))  # Rename outcome column similarly
  
  # Identify Best Predicted Drug ----
  common_prefix <- function(strings) {
    min_len <- min(nchar(strings))
    prefix <- ""
    for (i in seq_len(min_len)) {
      current_chars <- substr(strings, i, i)
      if (length(unique(current_chars)) == 1) {
        prefix <- paste0(prefix, current_chars[1])
      } else {
        break
      }
    }
    return(prefix)
  }
  
  prefix <- common_prefix(pred_cols)            # Extract common prefix
  drug_names <- gsub(prefix, "", pred_cols)     # Get drug names by removing prefix
  
  # Use a helper function get_best_drugs to add best predicted drug info (rank 1)
  interim_dataset <- get_best_drugs(
    data = pre_data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    rename("function_rank1_drug_name" := paste0(prefix, "rank1_drug_name"))  # Rename new column for clarity
  
  
  # Define concordance label: 1 if patient received predicted best drug, else 0 ----
  
  if (is.null(conc_tolerance)) {
    # Without tolerance: concordant if received the top predicted drug exactly
    interim_dataset <- interim_dataset %>%
      mutate(conc_disc_label = ifelse(dataset_drug_var == function_rank1_drug_name, 1, 0))
  } else {
    # With tolerance: allow concordance if drug is within tolerance of best predicted
    
    n_drugs_required <- length(pred_cols)                      # Number of drugs (columns)
    tolerance_vars <- paste0("tolerance_drug_", 1:n_drugs_required)  # Prepare names for tolerance variables
    
    interim_dataset <- get_best_drugs(
      data = interim_dataset,
      tolerance = conc_tolerance,
      column_names = pred_cols,
      final_var_name = prefix
    ) %>%
      rename("function_tolerance_drug_name" := paste0(prefix, "within_", conc_tolerance, "_of_best_drug_name")) %>%
      mutate(
        # Concordant if actual drug is in the tolerance drug names list
        conc_disc_label = ifelse(str_detect(function_tolerance_drug_name, paste0("\\b", dataset_drug_var, "\\b")), 1, 0)
      ) %>%
      # Split the tolerance drug names into a list of drugs
      mutate(drug_list = str_split(function_tolerance_drug_name, "\\s*[;,\\s]\\s*")) %>%
      # Ensure each drug list has length equal to number of drugs by repeating
      mutate(drug_list = map(drug_list, ~ rep(.x, length.out = n_drugs_required))) %>%
      # Name each element of the list for unnesting later
      mutate(drug_list = map(drug_list, ~ set_names(.x, paste0("tolerance_drug_", seq_along(.x))))) %>%
      # Expand drug_list columns to separate columns
      unnest_wider(drug_list) %>%
      # For patients not concordant, overwrite tolerance drug columns with actual drug to avoid mismatch
      mutate(across(
        all_of(tolerance_vars),
        ~ if_else(conc_disc_label == 0, dataset_drug_var, .x)
      ))
  }
  
  # Matching Procedure ----
  interim_dataset <- interim_dataset %>% drop_na(all_of(matching_var))
  
  categorical_vars <- matching_var[sapply(interim_dataset[matching_var], \(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
  
  # Build matching formula dynamically for MatchIt
  matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))
  for (v in categorical_vars) {
    if (length(unique(interim_dataset[[v]])) > 1) {
      matching_formula <- paste(matching_formula, "+", v)
    }
  }
  
  # Define antiexact matching variables based on tolerance use
  match_model_antiexact_vars <- if (!is.null(conc_tolerance)) {
    unique(c(tolerance_vars, "dataset_drug_var", match.antiexact))
  } else {
    unique(c("dataset_drug_var", match.antiexact))
  }
  
  match_model <- MatchIt::matchit(
    formula = as.formula(matching_formula),     # Formula for propensity score or distance calculation
    data = interim_dataset,
    method = "nearest",                         # Nearest neighbor matching
    distance = "mahalanobis",                   # Mahalanobis distance metric
    replace = FALSE,                            # No replacement in matching
    exact = unique(c("function_rank1_drug_name", match.exact)),  # Exact matching variables
    antiexact = match_model_antiexact_vars                 # Antiexact matching variables
  )
  
  # Extract matched pairs data
  matched_data <- MatchIt::get_matches(match_model, data = interim_dataset)
  
  # Compute Observed and Predicted Differences ----
  processed_data <- matched_data %>%
    group_by(subclass) %>%                            # Group by matched pair subclass
    mutate(
      concordant_drugclass  = dataset_drug_var[1],   # Drug of concordant patient in pair
      discordant_drugclass  = dataset_drug_var[2],   # Drug of discordant patient in pair
      calibration_obs       = diff(dataset_outcome_var) # Difference in outcome between discordant - concordant
    ) %>%
    ungroup() %>%
    distinct(subclass, .keep_all = TRUE) %>%          # Keep one row per pair
    rowwise() %>%
    mutate(
      # Calculate predicted difference between discordant and concordant drugs using predicted outcome columns
      calibration_pred = get(paste0(prefix, discordant_drugclass)) - get(paste0(prefix, concordant_drugclass))
    ) %>%
    ungroup() %>%
    select(calibration_pred, calibration_obs)        # Keep only needed columns
  
  # summary metrics
  linear_model <- lm(calibration_obs ~ calibration_pred, data = processed_data)
  
  # Compile final result
  result <- list(
    calibration_intercept = data.frame(
      row.names = "Intercept",
      value = coef(linear_model)[1],
      lci = confint(linear_model)[1, 1],
      uci = confint(linear_model)[1, 2]
    ),
    calibration_slope = data.frame(
      row.names = "Slope",
      value = coef(linear_model)[2],
      lci = confint(linear_model)[2, 1],
      uci = confint(linear_model)[2, 2]
    )
  )
  
  # Return result
  return(result)
  
}

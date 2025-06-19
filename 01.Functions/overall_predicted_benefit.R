#' Estimate Observed vs. Predicted Treatment Benefit by Calibration Group
#'
#' This function evaluates how well predicted treatment benefits correspond to observed differences in outcomes.
#' Patients are divided into calibration groups based on predicted benefit, and observed outcome differences
#' between concordant and discordant treatments are regressed on the predicted benefit difference.
#' Concordance is defined either by receiving the top predicted treatment or any within a specified tolerance of the best.
#' Matching is used to adjust for confounding.
#'
#' @param data A data frame containing treatment assignments, outcomes, covariates, and predicted benefits.
#' @param drug_var Character. Name of the column indicating the actual treatment received.
#' @param outcome_var Character. Name of the outcome variable.
#' @param cal_groups Integer. Number of calibration groups to divide the population into (e.g., 5 for quintiles).
#' @param pred_cols Character vector. Column names of predicted outcomes (or risks) for each treatment. Should share a common prefix (e.g., "pred_GLP1", "pred_SGLT2").
#' @param conc_tolerance Optional numeric. If supplied, defines the absolute tolerance for determining concordance with the best predicted treatment. If NULL, concordance is based on the top predicted treatment only.
#' @param matching_var Character vector. Covariate names used for Mahalanobis distance matching.
#' @param match.exact Optional character vector. Variables required to match exactly across groups. The best predicted treatment is always included.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching. The actual treatment variable is always included.
#'
#' @return A data frame with one row per calibration group, containing:
#' \describe{
#'   \item{mean}{Average predicted benefit in the group (discordant minus concordant).}
#'   \item{coef}{Estimated observed difference in outcome (regression coefficient).}
#'   \item{coef_low}{Lower bound of the 95\% confidence interval.}
#'   \item{coef_high}{Upper bound of the 95\% confidence interval.}
#'   \item{n_groups}{Total number of calibration groups.}
#'   \item{n}{Number of matched pairs in each group.}
#' }
#'
#' @details
#' The function uses `MatchIt` for nearest-neighbor matching to control confounding. Concordant patients are those who received their predicted-optimal treatment, while discordant patients received a suboptimal one. Within each calibration group, observed outcomes are compared using regression.
#'
#' @import tidyverse
#' @import MatchIt
#' @export
overall_predicted_benefit <- function(data, 
                                      drug_var, 
                                      outcome_var, 
                                      cal_groups, 
                                      pred_cols = NULL, 
                                      conc_tolerance = NULL,
                                      matching_var = NULL, 
                                      match.exact = NULL, 
                                      match.antiexact = NULL) {
  
  # Load Required Libraries ----
  require(tidyverse)  # For data manipulation and piping (%>%)
  require(MatchIt)    # For matching procedures (nearest neighbor matching)
  
  # Input Validation ----
  # Ensure required columns exist in data
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!all(pred_cols %in% colnames(data))) stop("Some pred_cols not found in data")
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("conc_tolerance must be numeric")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  
  # Prepare Dataset ----
  # Standardize variable names internally for consistency
  pre_data <- data %>%
    rename_with(~ "dataset_drug_var", all_of(drug_var)) %>%      # Rename actual treatment variable
    rename_with(~ "dataset_outcome_var", all_of(outcome_var))    # Rename outcome variable
  
  # Identify Common Prefix of Predicted Columns ----
  # This assumes predicted columns have the same prefix + drug name suffix
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
  
  prefix <- common_prefix(pred_cols)                    # Extract prefix from predicted columns
  drug_names <- gsub(prefix, "", pred_cols)             # Extract drug names by removing prefix
  
  # Identify Best Predicted Drug per Patient ----
  # Uses helper function `get_best_drugs` to assign best predicted drug (rank 1)
  interim_dataset <- get_best_drugs(
    data = pre_data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    rename("function_rank1_drug_name" := paste0(prefix, "rank1_drug_name"))  # Rename for clarity
  
  # Define Concordance Label ----
  # 1 if patient received predicted best drug, else 0
  if (is.null(conc_tolerance)) {
    # Concordant if actual treatment exactly matches the top predicted drug
    interim_dataset <- interim_dataset %>%
      mutate(conc_disc_label = ifelse(dataset_drug_var == function_rank1_drug_name, 1, 0))
  } else {
    # With tolerance: concordant if drug is within tolerance of best predicted outcome
    n_drugs_required <- length(pred_cols)
    tolerance_vars <- paste0("tolerance_drug_", 1:n_drugs_required)  # Variable names for tolerance drugs
    
    # Re-run get_best_drugs with tolerance to get drugs within tolerance
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
      # Split tolerance drug names into list for separate columns
      mutate(drug_list = str_split(function_tolerance_drug_name, "\\s*[;,\\s]\\s*")) %>%
      # Ensure lists have length equal to number of drugs by repeating elements if needed
      mutate(drug_list = map(drug_list, ~ rep(.x, length.out = n_drugs_required))) %>%
      # Name elements for unnesting
      mutate(drug_list = map(drug_list, ~ set_names(.x, paste0("tolerance_drug_", seq_along(.x))))) %>%
      # Expand list columns to separate tolerance_drug_X columns
      unnest_wider(drug_list) %>%
      # For patients not concordant, replace tolerance drug columns with actual drug to avoid mismatch in matching
      mutate(across(
        all_of(tolerance_vars),
        ~ if_else(conc_disc_label == 0, dataset_drug_var, .x)
      ))
  }
  
  # Filter Out Missing Matching Variables ----
  interim_dataset <- interim_dataset %>% drop_na(all_of(matching_var))
  
  # Identify categorical and continuous matching variables ----
  categorical_vars <- matching_var[sapply(interim_dataset[matching_var], \(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
  
  # Construct Matching Formula for MatchIt ----
  # Use continuous variables as additive terms
  matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))
  # Add categorical variables if they have more than one unique level
  for (v in categorical_vars) {
    if (length(unique(interim_dataset[[v]])) > 1) {
      matching_formula <- paste(matching_formula, "+", v)
    }
  }
  
  # Define variables for anti-exact matching (e.g. treatment variable and tolerance vars if applicable)
  match_model_antiexact_vars <- if (!is.null(conc_tolerance)) {
    unique(c(tolerance_vars, "dataset_drug_var", match.antiexact))
  } else {
    unique(c("dataset_drug_var", match.antiexact))
  }
  
  # Run Nearest Neighbor Matching using MatchIt ----
  match_model <- MatchIt::matchit(
    formula = as.formula(matching_formula),            # Formula for matching
    data = interim_dataset,
    method = "nearest",                                # Nearest neighbor matching
    distance = "mahalanobis",                          # Mahalanobis distance metric
    replace = FALSE,                                   # Matching without replacement
    exact = unique(c("function_rank1_drug_name", match.exact)), # Exact matching variables
    antiexact = match_model_antiexact_vars             # Anti-exact matching variables
  )
  
  # Extract Matched Data ----
  matched_data <- MatchIt::get_matches(match_model, data = interim_dataset)
  
  # Calculate Observed and Predicted Differences Within Matched Pairs ----
  processed_data <- matched_data %>%
    group_by(subclass) %>%                           # Group by matched pair group
    mutate(
      concordant_drugclass  = dataset_drug_var[1],  # Concordant patient's drug
      discordant_drugclass  = dataset_drug_var[2],  # Discordant patient's drug
      calibration_obs       = diff(dataset_outcome_var) # Observed outcome difference (discordant - concordant)
    ) %>%
    ungroup() %>%
    distinct(subclass, .keep_all = TRUE) %>%         # Keep one row per matched pair
    rowwise() %>%
    mutate(
      # Calculate predicted outcome difference using predicted columns for discordant vs concordant drug
      calibration_pred = get(paste0(prefix, discordant_drugclass)) - get(paste0(prefix, concordant_drugclass))
    ) %>%
    ungroup() %>%
    # Assign patients to calibration groups based on predicted difference quantiles
    select(calibration_pred, calibration_obs) %>%
    mutate(grouping = ntile(calibration_pred, cal_groups))
  
  # Initialize Output Vectors for Results ----
  coef      <- rep(NA_real_, cal_groups)      # Estimated observed differences (coefficients)
  coef_low  <- rep(NA_real_, cal_groups)      # Lower bounds of 95% CI
  coef_high <- rep(NA_real_, cal_groups)      # Upper bounds of 95% CI
  mean_vals <- rep(NA_real_, cal_groups)      # Mean predicted differences in groups
  n_group   <- rep(0, cal_groups)              # Number of pairs in each group
  
  # Fit Regression Models Within Each Calibration Group ----
  for (cg in seq_len(cal_groups)) {
    group_data <- processed_data %>% filter(grouping == cg) # Subset group
    
    mean_vals[cg] <- mean(group_data$calibration_pred, na.rm = TRUE) # Average predicted benefit
    n_group[cg] <- nrow(group_data)                               # Number of matched pairs
    
    # Fit linear regression: observed difference ~ predicted difference
    model <- lm(calibration_obs ~ calibration_pred, data = group_data)
    
    # Predict outcome and confidence intervals at average predictor values
    predicted_value <- predict(
      model,
      newdata = group_data %>% colMeans() %>% t() %>% as.data.frame(),
      interval = "confidence"
    ) %>% unlist()
    
    # Store coefficient and confidence interval bounds
    coef[cg]      <- predicted_value[1]
    coef_low[cg]  <- predicted_value[2]
    coef_high[cg] <- predicted_value[3]
  }
  
  # Compile Results into a Data Frame ----
  result <- data.frame(
    mean     = mean_vals,
    coef     = coef,
    coef_low = coef_low,
    coef_high= coef_high,
    n_groups = cal_groups,
    n        = n_group
  )
  
  return(result)  # Return summarized results
}

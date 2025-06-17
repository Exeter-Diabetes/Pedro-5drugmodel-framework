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
overall_predicted_benefit <- function(data, 
                                      drug_var, 
                                      outcome_var, 
                                      cal_groups, 
                                      pred_cols = NULL, 
                                      conc_tolerance = NULL,
                                      matching_var = NULL, 
                                      match.exact = NULL, 
                                      match.antiexact = NULL) {
  
  # load libraries
  require(tidyverse)
  require(MatchIt)
  
  # ---------------------------
  # Validate inputs
  # ---------------------------
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!all(pred_cols %in% colnames(data))) stop("Some pred_cols not found in data")
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("conc_tolerance has to be numeric when supplied")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  
  # ---------------------------
  # Filter and prepare dataset
  # ---------------------------
  pre_data <- data %>%
    rename_with(~ "dataset_drug_var", all_of(drug_var)) %>%
    rename_with(~ "dataset_outcome_var", all_of(outcome_var))
  
  # ---------------------------
  # Identify best drug for each patient
  # ---------------------------
  
  # Function to extract common prefix in predicted column names
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
  
  # Extract common prefix and drug names
  prefix <- common_prefix(pred_cols)
  drug_names <- gsub(prefix, "", pred_cols)
  
  # Get best drug by prediction (rank 1 and tolerance)
  interim_dataset <- get_best_drugs(
    data = pre_data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    dplyr::rename("function_rank1_drug_name" := paste0(prefix, "rank1_drug_name"))
  
  if (is.null(conc_tolerance)) {
    interim_dataset <- interim_dataset %>%
      dplyr::mutate(
        conc_disc_label = ifelse(dataset_drug_var == function_rank1_drug_name, 1, 0)
      )
  } else {
    interim_dataset <- get_best_drugs(
      data = interim_dataset,
      tolerance = conc_tolerance,
      column_names = pred_cols,
      final_var_name = prefix
    ) %>%
      dplyr::rename("function_tolerance_drug_name" := paste0(prefix, "within_", conc_tolerance, "_of_best_drug_name")) %>%
      dplyr::mutate(
        conc_disc_label = ifelse(str_detect(function_tolerance_drug_name, paste0("\\b", dataset_drug_var, "\\b")), 1, 0)
      )
  }
  
  # ---------------------------
  # Matching step
  # ---------------------------
  
  # Drop rows with missing matching variables
  interim_dataset <- interim_dataset %>%
    drop_na(all_of(matching_var))
  
  # Build matching formula
  categorical_vars <- matching_var[sapply(interim_dataset[matching_var], \(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
  
  matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))
  for (v in categorical_vars) {
    if (length(unique(interim_dataset[[v]])) > 1) {
      matching_formula <- paste(matching_formula, "+", v)
    }
  }
  
  # Perform matching
  match_model <- MatchIt::matchit(
    formula = as.formula(matching_formula),
    data = interim_dataset,
    method = "nearest",
    distance = "mahalanobis",
    replace = FALSE,
    exact = unique(c("function_rank1_drug_name", match.exact)),
    antiexact = unique(c("dataset_drug_var", match.antiexact))
  )
  
  matched_data <- MatchIt::get_matches(match_model, data = interim_dataset)
  
  # ---------------------------
  # Calculate observed and predicted calibration differences
  # ---------------------------
  processed_data <- matched_data %>%
    group_by(subclass) %>%
    mutate(
      concordant_drugclass  = dataset_drug_var[1],
      discordant_drugclass  = dataset_drug_var[2],
      calibration_obs       = diff(dataset_outcome_var)
    ) %>%
    ungroup() %>%
    distinct(subclass, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(
      calibration_pred = get(paste0(prefix, discordant_drugclass)) - get(paste0(prefix, concordant_drugclass))
    ) %>%
    ungroup() %>%
    dplyr::select(calibration_pred, calibration_obs) %>%
    dplyr::mutate(grouping = dplyr::ntile(calibration_pred, cal_groups))
  
  # ---------------------------
  # Initialize result vectors
  # ---------------------------
  coef      <- rep(NA_real_, cal_groups)
  coef_low  <- rep(NA_real_, cal_groups)
  coef_high <- rep(NA_real_, cal_groups)
  mean_vals <- rep(NA_real_, cal_groups)
  n_group   <- rep(0, cal_groups)
  
  # ---------------------------
  # Loop through calibration groups
  # ---------------------------
  for (cg in seq_len(cal_groups)) {
    
    # Subset for current group
    group_data <- processed_data %>%
      dplyr::filter(grouping == cg)
    
    # Record mean predicted benefit and sample size
    mean_vals[cg] <- mean(group_data$calibration_pred, na.rm = TRUE)
    n_group[cg]   <- nrow(group_data)
    
    # ---------------------------
    # Fit regression model
    # ---------------------------
    formula_str <- "calibration_obs ~ calibration_pred"
    
    model <- lm(as.formula(formula_str), data = group_data)
    
    predicted_value <- predict(model, newdata = group_data %>% colMeans() %>% t() %>% as.data.frame(), interval = "confidence") %>% unlist()
    
    coef[cg]      <- predicted_value[1]
    coef_low[cg]  <- predicted_value[2]
    coef_high[cg] <- predicted_value[3]
    
  }
  
  # ---------------------------
  # Compile final results
  # ---------------------------
  result <- data.frame(
    mean     = mean_vals,
    coef     = coef,
    coef_low = coef_low,
    coef_high= coef_high,
    n_groups = cal_groups,
    n        = n_group
  )
  
  # Return result
  return(result)
  
}

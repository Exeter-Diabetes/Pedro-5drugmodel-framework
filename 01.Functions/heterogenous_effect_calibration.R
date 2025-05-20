#' Estimate Heterogeneous Treatment Effects Across Calibration Groups
#'
#' Splits the data into groups based on predicted benefit and estimates the treatment effect
#' (and confidence intervals) between two drugs within each group using regression models.
#'
#' @param data A data frame containing observed outcomes, treatment assignments, and predicted benefits.
#' @param drug_var Character string. Column name for treatment assignment (drug).
#' @param drugs Character vector of length 2. Names of the two drugs to compare.
#' @param benefit_var Character string. Column name containing predicted benefit scores.
#' @param outcome_var Character string. Column name for the outcome variable (e.g., clinical measurement).
#' @param cal_groups Integer. Number of calibration groups to split the population into.
#' @param adjustment_var Optional character vector. Names of covariates to adjust for in the regression.
#'
#' @return A data frame containing one row per calibration group with:
#' \describe{
#'   \item{mean}{Mean predicted benefit for the group.}
#'   \item{coef}{Estimated treatment effect (regression coefficient).}
#'   \item{coef_low}{Lower bound of 95% confidence interval.}
#'   \item{coef_high}{Upper bound of 95% confidence interval.}
#'   \item{grouping}{Calibration group identifier.}
#'   \item{drug1}{Name of the first drug.}
#'   \item{drug2}{Name of the second drug.}
#' }
#'
#' @examples
#' \dontrun{
#' heterogenous_effect_calibration(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 5,
#'   adjustment_var = c("age", "bmi")
#' )
#' }
#' @export
heterogenous_effect_calibration <- function(data,
                                            drug_var,
                                            drugs,
                                            benefit_var,
                                            outcome_var,
                                            cal_groups,
                                            matching = FALSE,
                                            adjustment_var = NULL,
                                            matching_var = adjustment_var) {
  
  # load libraries
  require(tidyverse)
  require(MatchIt)
  
  # ---------------------------
  # Validate inputs
  # ---------------------------
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(benefit_var %in% colnames(data))) stop("benefit_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!is.null(adjustment_var) && !all(adjustment_var %in% colnames(data))) stop("Some adjustment_var not in data")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (length(drugs) != 2) stop("Exactly two drugs must be specified")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some specified drugs not present in drug_var column")
  if (!is.numeric(cal_groups) || cal_groups <= 0) stop("cal_groups must be a positive number")
  if (!is.logical(matching)) stop("matching must be TRUE or FALSE")
  
  # ---------------------------
  # Prepare and group the data
  # ---------------------------
  initial_dataset <- data %>%
    dplyr::rename(
      dataset_benefit = !!benefit_var,
      dataset_drug_var = !!drug_var,
      dataset_outcome_var = !!outcome_var
    ) %>%
    dplyr::filter(dataset_drug_var %in% drugs)
  
  
  # ---------------------------
  # Optional matching
  # ---------------------------
  
  if (isTRUE(matching)) {
  
    # drop missing data
    initial_dataset <- initial_dataset %>%
      drop_na(matching_var)
    
    # matching formula
    categorical_vars <- matching_var[sapply(initial_dataset[matching_var], \(x) is.factor(x) || is.character(x))]
    cont_vars <- setdiff(matching_var, categorical_vars)
    
    matching_formula <- paste("dataset_drug_var ~", paste(cont_vars, collapse = " + "))
    for (v in categorical_vars) {
      if (length(unique(initial_dataset[[v]])) > 1) {
        matching_formula <- paste(matching_formula, "+", v)
      }
    }
    
    match_model <- MatchIt::matchit(
      formula = as.formula(matching_formula),
      data = initial_dataset,
      method = "nearest",
      distance = "mahalanobis",
      replace = FALSE
    )
    
    calibration_data <- MatchIt::get_matches(match_model, data = initial_dataset)
    
  } else {
    
    calibration_data <- initial_dataset
    
  }
  
  # ---------------------------
  # Initialize result vectors
  # ---------------------------
  coef      <- rep(NA_real_, cal_groups)
  coef_low  <- rep(NA_real_, cal_groups)
  coef_high <- rep(NA_real_, cal_groups)
  mean_vals <- rep(NA_real_, cal_groups)
  n_group   <- rep(0, cal_groups)
  
  
  # ---------------------------
  # Grouping patients
  # ---------------------------
  
  calibration_data <- calibration_data %>%
    dplyr::mutate(grouping = dplyr::ntile(dataset_benefit, cal_groups))
  
  
  # ---------------------------
  # Loop through calibration groups
  # ---------------------------
  for (g in seq_len(cal_groups)) {
    
    # Subset and factorize drug variable
    group_data <- calibration_data %>%
      dplyr::filter(grouping == g) %>%
      dplyr::mutate(dataset_drug_var = factor(dataset_drug_var, levels = rev(drugs)))
    
    # Record group size and mean benefit
    mean_vals[g] <- mean(group_data$dataset_benefit, na.rm = TRUE)
    n_group[g]   <- nrow(group_data)
    
    if (length(unique(group_data$dataset_drug_var)) < 2) next
    
    # ---------------------------
    # Build regression formula
    # ---------------------------
    formula_str <- "dataset_outcome_var ~ dataset_drug_var"
    
    if (!is.null(adjustment_var)) {
      cat_vars  <- adjustment_var[sapply(group_data[, adjustment_var], function(x) is.factor(x) || is.character(x))]
      cont_vars <- setdiff(adjustment_var, cat_vars)
      
      if (length(cont_vars) > 0) {
        formula_str <- paste(formula_str, paste(cont_vars, collapse = " + "), sep = " + ")
      }
      for (v in cat_vars) {
        if (length(unique(group_data[[v]])) > 1) {
          formula_str <- paste(formula_str, v, sep = " + ")
        }
      }
    }
    
    # ---------------------------
    # Fit model and extract estimates
    # ---------------------------
    model <- glm(as.formula(formula_str), data = group_data)
    ci <- suppressMessages(confint.default(model))
    
    coef[g]      <- coef(model)[2]
    coef_low[g]  <- ci[2, 1]
    coef_high[g] <- ci[2, 2]
  }
  
  # ---------------------------
  # Compile results into a data frame
  # ---------------------------
  result <- data.frame(
    mean     = mean_vals,
    coef     = coef,
    coef_low = coef_low,
    coef_high= coef_high,
    n        = n_group,
    grouping = seq_len(cal_groups),
    drug1    = drugs[1],
    drug2    = drugs[2]
  )
  
  # Return result table
  return(result)
  
}

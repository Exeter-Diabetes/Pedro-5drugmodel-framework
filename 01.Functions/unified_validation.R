#' Perform Pairwise Heterogeneous Treatment Effect Calibration Across Multiple Drugs
#'
#' This function computes heterogeneous treatment effect calibration curves for all pairwise
#' comparisons between a set of specified drugs. For each drug pair, it calculates the predicted
#' benefit (the difference between their predicted outcomes), partitions the data into calibration
#' groups based on predicted benefit, and estimates treatment effects within each group.
#'
#' The function internally calls `heterogenous_effect_calibration()` for each drug pair.
#'
#' @param data A data.frame containing observed outcomes, treatment assignments, and predicted outcomes for each drug.
#' @param drug_var Character string specifying the column name for the treatment variable (e.g., assigned drug).
#' @param drugs Character vector of two or more drug names to be compared pairwise.
#' @param prediction_vars Character vector of the same length as `drugs`, with the column names of the predicted outcomes for each corresponding drug.
#' @param outcome_var Character string specifying the column name for the observed outcome variable.
#' @param cal_groups Numeric scalar or vector. The number(s) of calibration groups (quantiles) to partition the population into for subgroup analysis.
#' @param matching Logical. If `TRUE`, performs covariate matching using the `MatchIt` package before estimating treatment effects.
#' @param adjustment_var Optional character vector of column names to include as covariates in the regression models.
#' @param matching_var Optional character vector specifying variables to match on. Defaults to `adjustment_var` if not specified.
#' @param match.exact Optional character vector. Variables for exact matching. Matching on best predicted drug automatically added.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching. drug_var automatically added.
#'
#' @return A data.frame with one row per calibration group and drug pair combination, containing:
#' \describe{
#'   \item{cal_groups}{Number of calibration groups used.}
#'   \item{grouping}{Calibration group identifier (e.g., 1 through `cal_groups`).}
#'   \item{mean}{Mean predicted benefit in the group.}
#'   \item{coef}{Estimated treatment effect (regression coefficient).}
#'   \item{coef_low}{Lower bound of the 95% confidence interval.}
#'   \item{coef_high}{Upper bound of the 95% confidence interval.}
#'   \item{drug1}{Name of the first drug (as provided).}
#'   \item{n_drug1}{Number of individuals in the group with the first drug.}
#'   \item{drug2}{Name of the second drug (as provided).}
#'   \item{n_drug2}{Number of individuals in the group with the second drug.}
#' }
#'
#' @examples
#' \dontrun{
#' unified_validation(
#'   data = mydata,
#'   drug_var = "treatment",
#'   drugs = c("DrugA", "DrugB", "DrugC"),
#'   prediction_vars = c("pred_drugA", "pred_drugB", "pred_drugC"),
#'   outcome_var = "outcome",
#'   cal_groups = c(3, 5),
#'   adjustment_var = c("age", "sex"),
#'   matching = TRUE
#' )
#' }
#'
#' @export
unified_validation <- function(data, 
                               drug_var, 
                               drugs, 
                               prediction_vars, 
                               outcome_var, 
                               cal_groups,
                               matching = FALSE,
                               adjustment_var = NULL,
                               matching_var = adjustment_var,
                               match.exact = NULL, 
                               match.antiexact = NULL) {
  
  # load libraries
  require(tidyverse)
  
  # ---------------------------
  # Input validation
  # ---------------------------
  if (!(drug_var %in% colnames(data))) stop("`drug_var` not found in data.")
  if (!all(prediction_vars %in% colnames(data))) stop("Some `prediction_vars` not found in data.")
  if (!(outcome_var %in% colnames(data))) stop("`outcome_var` not found in data.")
  if (!is.null(adjustment_var) && !all(adjustment_var %in% colnames(data))) stop("Some `adjustment_var` columns not found in data.")
  if (isTRUE(matching)) {
    if (length(matching_var) == 0) stop("Provide at least one matching_var")
    if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
    if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
    if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  }
  if (length(drugs) < 2) stop("At least two drugs must be specified.")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some `drugs` not present in the `drug_var` column.")
  if (length(drugs) != length(prediction_vars)) stop("`drugs` and `prediction_vars` must have the same length.")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  if (!is.logical(matching)) stop("matching must be TRUE or FALSE")
  
  # ---------------------------
  # Prepare data
  # ---------------------------
  calibration_data <- data %>%
    rename(
      dataset_drug_var = !!drug_var,
      dataset_outcome_var = !!outcome_var
    ) %>%
    filter(dataset_drug_var %in% drugs)
  
  prediction_map <- setNames(prediction_vars, drugs)
  output_table <- NULL
  
  # ---------------------------
  # Iterate over all drug pairs
  # ---------------------------
  drug_combinations <- combn(drugs, 2, simplify = FALSE)
  
  for (pair in drug_combinations) {
    
    # Subset and prepare data
    pair_data <- calibration_data %>%
      filter(dataset_drug_var %in% pair) %>%
      mutate(
        dataset_drug_var = factor(dataset_drug_var, levels = rev(pair)),
        drug_1_pred = .[[prediction_map[[pair[1]]]]],
        drug_2_pred = .[[prediction_map[[pair[2]]]]],
        benefit = drug_1_pred - drug_2_pred
      )
    
    # Run calibration
    calibration_result <- heterogenous_effect_calibration(
      data = pair_data,
      drug_var = "dataset_drug_var",
      drugs = pair,
      benefit_var = "benefit",
      outcome_var = "dataset_outcome_var",
      cal_groups = cal_groups,
      matching= matching,
      adjustment_var = adjustment_var,
      matching_var = matching_var,
      match.exact = match.exact, 
      match.antiexact = match.antiexact
    )
    
    # Store results
    output_table <- dplyr::bind_rows(output_table, calibration_result)
  }
  
  return(output_table)
}

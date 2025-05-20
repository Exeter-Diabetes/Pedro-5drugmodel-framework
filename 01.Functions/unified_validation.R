#' Perform Pairwise Heterogeneous Treatment Effect Calibration
#'
#' This function computes calibration curves for all pairwise comparisons between specified drugs,
#' based on their respective predicted benefits and observed outcomes. It internally calls 
#' `heterogenous_effect_calibration()` for each pair.
#'
#' @param data A data.frame containing the dataset.
#' @param drug_var A character string specifying the column name for the treatment variable.
#' @param drugs A character vector of exactly two or more drug names to compare.
#' @param prediction_vars A character vector of the same length as `drugs`, containing
#'   the column names of the predicted outcomes for each corresponding drug.
#' @param outcome_var A character string specifying the column name for the observed outcome.
#' @param cal_groups An integer specifying the number of calibration subgroups (quantiles).
#' @param adjustment_var Optional character vector of column names to include as adjustment variables in the models.
#'
#' @return A data.frame with calibration group statistics for each drug pair:
#'   mean predicted benefit, estimated treatment effect, confidence intervals, and drug pair labels.
#'
#' @examples
#' \dontrun{
#' unified_validation(
#'   data = mydata,
#'   drug_var = "treatment",
#'   drugs = c("DrugA", "DrugB"),
#'   prediction_vars = c("pred_drugA", "pred_drugB"),
#'   outcome_var = "outcome",
#'   cal_groups = 5,
#'   adjustment_var = c("age", "sex")
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
                               matching_var = adjustment_var) {
  
  # load libraries
  require(tidyverse)
  
  # ---------------------------
  # Input validation
  # ---------------------------
  if (!(drug_var %in% colnames(data))) stop("`drug_var` not found in data.")
  if (!all(prediction_vars %in% colnames(data))) stop("Some `prediction_vars` not found in data.")
  if (!(outcome_var %in% colnames(data))) stop("`outcome_var` not found in data.")
  if (!is.null(adjustment_var) && !all(adjustment_var %in% colnames(data))) stop("Some `adjustment_var` columns not found in data.")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (length(drugs) < 2) stop("At least two drugs must be specified.")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some `drugs` not present in the `drug_var` column.")
  if (length(drugs) != length(prediction_vars)) stop("`drugs` and `prediction_vars` must have the same length.")
  if (!is.numeric(cal_groups) || cal_groups <= 0) stop("`cal_groups` must be a positive number.")
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
      matching_var = matching_var
    )
    
    # Store results
    output_table <- dplyr::bind_rows(output_table, calibration_result)
  }
  
  return(output_table)
}
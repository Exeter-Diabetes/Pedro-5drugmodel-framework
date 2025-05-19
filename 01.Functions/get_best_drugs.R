#' Assign Best Drugs Based on Prediction Rankings or Tolerance
#'
#' Applies \code{\link{get_ranked_or_tolerant_drugs}} across all rows of a dataset to determine
#' the best (or nearly best) drugs based on model predictions. Adds the result as new columns to the data.
#'
#' @param data A data frame containing predicted values in specified columns.
#' @param rank Integer. Rank of the prediction to select (e.g., 1 = best, 2 = second-best). Ignored if `tolerance` is set.
#' @param column_names Character vector. Names of the columns containing predicted values.
#' @param final_var_name Character string. Base name for output columns.
#' @param tolerance Optional numeric. If set, returns all drugs within this tolerance of the best.
#'
#' @return A modified version of the original data frame with two new columns:
#' \describe{
#'   \item{<label>_drug_value}{Predicted value(s) of the selected drug(s).}
#'   \item{<label>_drug_name}{Name(s) of the selected drug(s).}
#' }
#'
#' @examples
#' \dontrun{
#' preds <- data.frame(pred_A = c(0.1, 0.3), pred_B = c(0.2, 0.2), pred_C = c(0.3, 0.1))
#' get_best_drugs(preds, rank = 1, column_names = names(preds), final_var_name = "best")
#' get_best_drugs(preds, tolerance = 0.05, column_names = names(preds), final_var_name = "tolerant")
#' }
#' @export
# Wrapper function to apply row-wise
get_best_drugs <- function(data, rank = 1, column_names = NULL, final_var_name = "", tolerance = NULL) {
  
  # ---------------------------
  # Validate inputs
  # ---------------------------
  if (is.null(column_names) || !all(column_names %in% colnames(data))) {stop("Invalid or missing column_names")}
  
  # ---------------------------
  # Apply ranking or tolerance logic to each row
  # ---------------------------
  # Calls 'get_ranked_or_tolerant_drugs' row-wise on the selected columns.
  # Depending on whether 'rank' or 'tolerance' is specified, it extracts the top-ranked or sufficiently good drug.
  # Returns a matrix with two columns: 'value' (e.g., predicted score) and 'name' (e.g., drug name).
  results <- t(apply(data[, column_names], 1, get_ranked_or_tolerant_drugs, 
                     rank = rank, tolerance = tolerance, column_names = column_names, prediction_column = final_var_name))
  
  # Assign standard column names to the result matrix
  colnames(results) <- c("value", "name")
  
  # ---------------------------
  # Generate label for new columns
  # ---------------------------
  # If a custom label is not provided, create one based on ranking or tolerance logic.
  # This label will prefix the output column names to describe what they contain.
  label <- if (!is.null(tolerance)) paste0(final_var_name, "within_", tolerance, "_of_best") else paste0(final_var_name, "rank", rank)
  
  # ---------------------------
  # Append new columns to original dataset
  # ---------------------------
  # Add the drug value column to the dataset. If tolerance is not used, coerce to numeric.
  if (is.null(tolerance)) {
    data[[paste0(label, "_drug_value")]] <- as.numeric(results[, "value"])
  } else {
    data[[paste0(label, "_drug_value")]] <- results[, "value"]
  }
  
  # Add the drug name column (e.g., best drug name by rank or tolerance)
  data[[paste0(label, "_drug_name")]] <- results[, "name"]
  
  # ---------------------------
  # Return modified dataset
  # ---------------------------
  return(data)
  
}

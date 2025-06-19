#' Assign Best Drugs Based on Prediction Rankings or Tolerance
#'
#' Applies \code{\link{get_ranked_or_tolerant_drugs}} across all rows of a dataset to determine the best (or nearly best) drugs based on model predictions. Adds the result as new columns to the data.
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
# Main wrapper function to apply drug selection logic row-wise
get_best_drugs <- function(data, rank = 1, column_names = NULL, final_var_name = "", tolerance = NULL) {
  
  # Input validation ----
  
  # Check that column_names are provided and exist in the data
  if (is.null(column_names) || !all(column_names %in% colnames(data))) {
    stop("Invalid or missing column_names")
  }
  
  # Apply selection logic row-wise ----
  
  # Apply get_ranked_or_tolerant_drugs to each row using only the specified columns
  # Results should be a two-column matrix: value (predicted score), name (drug name)
  results <- do.call(rbind, apply(
    data[, column_names], 
    1, 
    get_ranked_or_tolerant_drugs, 
    rank = rank, 
    tolerance = tolerance, 
    column_names = column_names, 
    prediction_column = final_var_name
  ))
  
  # Set column names of the results for clarity
  colnames(results) <- c("value", "name")
  
  # Create label for output columns ----
  
  # Generate a label that describes the selection logic used (rank or tolerance)
  label <- if (!is.null(tolerance)) {
    paste0(final_var_name, "within_", tolerance, "_of_best")
  } else {
    paste0(final_var_name, "rank", rank)
  }
  
  # Append new columns to dataset ----
  
  # Add predicted value(s) to the data
  # If using rank, coerce to numeric to ensure consistent column type
  if (is.null(tolerance)) {
    data[[paste0(label, "_drug_value")]] <- as.numeric(results[, "value"])
  } else {
    data[[paste0(label, "_drug_value")]] <- results[, "value"]
  }
  
  # Add selected drug name(s) to the data
  data[[paste0(label, "_drug_name")]] <- results[, "name"]
  
  # Return updated dataset ----
  
  # Return the modified dataset including the new columns
  return(data)
}
#' Get Ranked or Tolerance-Based Drug Recommendation
#'
#' This function identifies the drug(s) with the best (lowest) predicted outcome
#' from a set of model predictions. It can return the nth-best drug or all drugs
#' within a specified tolerance of the best drug.
#'
#' @param row A named vector or one-row data frame of predictions for a single patient.
#' @param rank Integer. Rank of drug to return (e.g., 1 = best, 2 = second-best). Ignored if `tolerance` is provided.
#' @param column_names Character vector. Names of the columns containing predicted values.
#' @param prediction_column Character string. Common prefix used in prediction column names (e.g., "pred_").
#' @param tolerance Optional numeric. If provided, return all drugs within this margin above the best predicted value.
#'
#' @return A character vector with two named elements:
#' \describe{
#'   \item{value}{Comma-separated string of predicted values (or single value if using `rank`).}
#'   \item{name}{Comma-separated string of corresponding drug names (or single name if using `rank`).}
#' }
#'
#' @examples
#' \dontrun{
#' row <- c(pred_A = 0.2, pred_B = 0.5, pred_C = 0.21)
#' get_ranked_or_tolerant_drugs(row, rank = 1, column_names = names(row), prediction_column = "pred_")
#' get_ranked_or_tolerant_drugs(row, tolerance = 0.05, column_names = names(row), prediction_column = "pred_")
#' }
#' @export
#' 
# Function to find nth-best drug OR drugs within a tolerance
get_ranked_or_tolerant_drugs <- function(row, rank = 1, column_names = NULL, prediction_column = NULL, tolerance = NULL) {
  
  # ---------------------------
  # Extract relevant predictions
  # ---------------------------
  # Subset the current row to only include the prediction columns for drugs
  values <- row[column_names]
  
  # Sort predictions in ascending order (lower = better; e.g., lower risk or error)
  sorted_values <- sort(values)
  
  if (!is.null(tolerance)) {
    # ---------------------------
    # Tolerance mode
    # ---------------------------
    # Identify all drugs whose predicted values are within 'tolerance' units of the best prediction
    best_val <- sorted_values[1]  # Best (lowest) prediction value
    close_vals <- sorted_values[sorted_values <= best_val + tolerance]  # Select predictions within range
    
    # Remove the prediction prefix to get just the drug names
    drug_names <- gsub(prediction_column, "", names(close_vals))
    
    # Return comma-separated lists of values and drug names
    return(c(
      value = paste(close_vals[order(drug_names)], collapse = ","),
      name  = paste(drug_names[order(drug_names)], collapse = ",")
    ))
    
  } else {
    # ---------------------------
    # Ranking mode
    # ---------------------------
    # Return the nth-best prediction value and corresponding drug
    
    # If requested rank exceeds number of available drugs, return NA
    if (rank > length(sorted_values)) {
      return(c(value = NA, name = NA))
    }
    
    # Extract the prediction value and corresponding drug name at the specified rank
    best_value <- sorted_values[rank]
    best_name <- gsub(prediction_column, "", names(sorted_values)[rank])
    
    return(c(value = best_value, name = best_name))
  }
}
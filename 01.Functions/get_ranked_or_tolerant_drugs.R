#' Get Ranked or Tolerance-Based Drug Recommendation
#'
#' This function identifies the drug(s) with the best (lowest) predicted outcome from a set of model predictions. It can return the nth-best drug or all drugs within a specified tolerance of the best drug.
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
  
  # Extract Relevant Predictions ----
  
  # Subset only the prediction values from the specified columns
  values <- row[column_names]
  
  # Sort the predictions in ascending order (lower values are better)
  sorted_values <- sort(values)
  
  # Tolerance Mode ----
  
  if (!is.null(tolerance)) {
    
    # Get the best (lowest) prediction value
    best_val <- sorted_values[1]
    
    # Select all drugs within the tolerance of the best value
    close_vals <- sorted_values[sorted_values <= best_val + tolerance]
    
    # Remove the prediction prefix (e.g., "pred_") to extract drug names
    drug_names <- gsub(prediction_column, "", names(close_vals))
    
    # Return comma-separated predicted values and drug names, sorted alphabetically by name
    return(c(
      value = paste(close_vals[order(drug_names)], collapse = ","),
      name  = paste(drug_names[order(drug_names)], collapse = ",")
    ))
    
  } else {
    
    # Ranking Mode ----
    
    # If the requested rank is beyond available predictions, return NAs
    if (rank > length(sorted_values)) {
      return(c(value = NA, name = NA))
    }
    
    # Get the value and drug name at the specified rank
    best_value <- sorted_values[rank]
    best_name <- gsub(prediction_column, "", names(sorted_values)[rank])
    
    # Return single best value and name
    return(c(value = best_value, name = best_name))
  }
}
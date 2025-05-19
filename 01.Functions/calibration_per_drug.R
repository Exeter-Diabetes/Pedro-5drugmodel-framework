#' Calibration Regression per Drug
#'
#' Performs calibration analysis for a specific drug by matching concordant and discordant treatment cases
#' and assessing the relationship between predicted and observed outcome differences across quantile groups.
#'
#' @param data A data frame containing the relevant variables.
#' @param drug_var Character string. Name of the column indicating the assigned treatment/drug.
#' @param drug Character string. The drug for which calibration is being evaluated.
#' @param pred_var Character string. Prefix of predicted probabilities column names (e.g., "pred_" if columns are named "pred_drugA", "pred_drugB", etc.).
#' @param outcome_var Character string. Name of the outcome variable column.
#' @param best_drug_var Character string. Column name indicating the model's best predicted treatment.
#' @param cal_groups Integer. Number of quantile groups to create for calibration analysis (default is 5).
#' @param match.vars Character vector. Variables to include in the matching process.
#' @param match.exact Optional character vector. Variables to match exactly.
#' @param match.antiexact Optional character vector. Variables to avoid matching exactly.
#'
#' @return A data frame summarizing calibration per quantile group:
#' \describe{
#'   \item{mean}{Median predicted difference in outcome between discordant and concordant treatment pairs.}
#'   \item{coef}{Estimated difference in observed outcomes regressed on predicted difference.}
#'   \item{coef_low}{Lower bound of 95% confidence interval.}
#'   \item{coef_high}{Upper bound of 95% confidence interval.}
#'   \item{n}{Number of matched pairs in the group.}
#'   \item{total_conc}{Total number of concordant cases across all matches.}
#'   \item{total_disconc}{Total number of discordant cases across all matches.}
#' }
#'
#' @import tidyverse
#' @import MatchIt
#' @export
#'
#' @examples
#' \dontrun{
#' calibration_per_drug(data = mydata,
#'                      drug_var = "drugclass",
#'                      drug = "DrugA",
#'                      pred_var = "pred.",
#'                      outcome_var = "posthba1c",
#'                      best_drug_var = "best_treatment",
#'                      match.vars = c("age", "sex", "comorbidity"))
#' }


calibration_per_drug <- function(data, drug_var, drug, pred_var, outcome_var, best_drug_var, cal_groups = 5, match.vars, match.exact = NULL, match.antiexact = NULL) {
  
  # load libraries
  require(tidyverse)
  require(MatchIt)
  
  # ---------------------------
  # Validate inputs
  # ---------------------------
  required_vars <- list(
    drug_var = drug_var,
    pred_var = pred_var,
    outcome_var = outcome_var,
    best_drug_var = best_drug_var,
    drug = drug,
    match.vars = match.vars
  )
  missing_vars <- names(Filter(is.null, required_vars))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required arguments:", paste(missing_vars, collapse = ", ")))
  }
  
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  
  var_check <- c(drug_var, best_drug_var, outcome_var)
  if (!all(var_check %in% colnames(data))) stop("Some specified variables not found in data")
  
  drugs_in_data <- unique(unlist(data[[drug_var]]))
  if (!all(drug %in% drugs_in_data)) stop("Specified drug not found in drug_var column")
  
  pred_columns <- paste0(pred_var, drugs_in_data)
  if (!all(pred_columns %in% colnames(data))) stop("Missing predicted probability columns")
  
  # ---------------------------
  # Filter and prepare dataset
  # ---------------------------
  pre_data <- data %>%
    rename_with(~ "dataset_drug_var", all_of(drug_var)) %>%
    rename_with(~ "dataset_best_drug_var", all_of(best_drug_var)) %>%
    rename_with(~ "dataset_outcome_var", all_of(outcome_var)) %>%
    filter(dataset_best_drug_var == drug) %>%
    mutate(conc_discordant_label = as.integer(dataset_drug_var == dataset_best_drug_var))
  
  # ---------------------------
  # Define matching formula
  # ---------------------------
  categorical_vars <- match.vars[sapply(pre_data[match.vars], \(x) is.factor(x) || is.character(x))]
  
  formula_terms <- setdiff(match.vars, c(categorical_vars, match.exact, match.antiexact))
  matching_formula <- paste0("conc_discordant_label ~ ", paste(formula_terms, collapse = " + "))
  
  # Avoids error if categorical variable only has one entry
  for (var in categorical_vars) {
    if (length(unique(pre_data[[var]])) > 1) {
      matching_formula <- paste(matching_formula, "+", var)
    }
  }
  
  # ---------------------------
  # Run matching
  # ---------------------------
  match_model <- MatchIt::matchit(
    formula = as.formula(matching_formula),
    data = pre_data,
    method = "nearest",
    distance = "mahalanobis",
    replace = FALSE, # changed from replace = TRUE
    exact = match.exact,
    antiexact = match.antiexact
  )
  
  matched_data <- MatchIt::get_matches(match_model, data = pre_data)
  
  # ---------------------------
  # Calculate observed and predicted outcomes
  # ---------------------------
  processed_data <- matched_data %>%
    group_by(subclass) %>%
    mutate(
      concordant_drugclass = dataset_drug_var[1],
      discordant_drugclass = dataset_drug_var[2],
      calibration_obs = diff(dataset_outcome_var)
    ) %>%
    ungroup() %>%
    distinct(subclass, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(
      calibration_pred = get(paste0(pred_var, discordant_drugclass)) - get(paste0(pred_var, concordant_drugclass))
    ) %>%
    ungroup() %>%
    select(calibration_pred, calibration_obs) %>%
    mutate(grouping = ntile(calibration_pred, cal_groups))
  
  # ---------------------------
  # Run calibration regression per group
  # ---------------------------
  calibration_summary <- processed_data %>%
    group_by(group_id = grouping) %>%
    summarise(
      mean = median(calibration_pred, na.rm = TRUE),
      n = n(),
      model = list(lm(calibration_obs ~ calibration_pred, data = cur_data())),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      prediction = list({
        pred <- predict(model, newdata = data.frame(calibration_pred = mean), interval = "confidence")
        # Ensure prediction has named columns
        colnames(pred) <- c("fit", "lwr", "upr")
        pred[1, ]  # return as named vector for unnesting
      })
    ) %>%
    unnest_wider(prediction, names_sep = "_") %>%
    rename(
      coef = prediction_fit,
      coef_low = prediction_lwr,
      coef_high = prediction_upr
    ) %>%
    select(mean, coef, coef_low, coef_high, n) %>%
    as.data.frame()
  
  # ---------------------------
  # Add total concordant/discordant counts
  # ---------------------------
  calibration_summary$total_conc <- matched_data %>% filter(conc_discordant_label == 1) %>% nrow()
  calibration_summary$total_disconc <- matched_data %>% filter(conc_discordant_label == 0) %>% nrow()
  
  return(calibration_summary)
  
  
} 

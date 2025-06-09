#' Estimate Heterogeneous Treatment Effects Across Calibration Groups
#'
#' Estimates heterogeneous treatment effects between two drugs by stratifying patients
#' into calibration groups based on predicted benefit scores. Optionally performs 
#' covariate matching before estimating treatment effects.
#'
#' If multiple values are provided to `cal_groups`, the function will repeat the process
#' for each value and return a combined result.
#'
#' @param data A data frame containing observed outcomes, treatment assignments, and predicted benefits.
#' @param drug_var Character string. Column name for the treatment assignment variable (e.g., drug).
#' @param drugs Character vector of length 2. Names of the two drugs to compare.
#' @param benefit_var Character string. Column name containing predicted benefit scores.
#' @param outcome_var Character string. Column name for the outcome variable (e.g., clinical measurement).
#' @param cal_groups Numeric or numeric vector. Number(s) of calibration groups (e.g., quantiles) to divide the data into based on predicted benefit.
#' @param matching Logical. Whether to perform covariate matching using the `MatchIt` package before estimating treatment effects.
#' @param adjustment_var Optional character vector. Names of covariates to include as adjustment variables in the regression model.
#' @param matching_var Optional character vector. Covariates to use for matching. Defaults to `adjustment_var` if not specified.
#' @param match.exact Optional character vector. Variables for exact matching. Matching on best predicted drug automatically added.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching. drug_var automatically added.
#'
#' @return A data frame with one row per calibration group and per `cal_groups` value, containing:
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
#' heterogenous_effect_calibration(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = c(3, 5),
#'   adjustment_var = c("age", "bmi"),
#'   matching = TRUE
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
                                            matching_var = adjustment_var,
                                            match.exact = NULL, 
                                            match.antiexact = NULL) {
  
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
  if (isTRUE(matching)) {
    if (length(matching_var) == 0) stop("Provide at least one matching_var")
    if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
    if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
    if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  }
  if (length(drugs) != 2) stop("Exactly two drugs must be specified")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some specified drugs not present in drug_var column")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  if (!is.logical(matching)) stop("matching must be TRUE or FALSE")
  
  # ---------------------------
  # Iterate through groups
  # ---------------------------
  
  result <- NULL
  
  for (cg in cal_groups) {
    if (cg <= 0) stop("Each value in cal_groups must be a positive number")
    
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
      
      # ---------------------------
      # Identify best drug for each patient
      # ---------------------------
      
      initial_dataset <- initial_dataset %>%
        dplyr::mutate(
          conc_disc_label = ifelse(dataset_benefit <= 0, 1, 0)
        )
      
      # drop missing data
      initial_dataset <- initial_dataset %>%
        drop_na(matching_var)
      
      # matching formula
      categorical_vars <- matching_var[sapply(initial_dataset[matching_var], \(x) is.factor(x) || is.character(x))]
      cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
      
      matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))  # changed from dataset_drug_var 
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
        replace = FALSE,
        exact = match.exact,
        antiexact = match.antiexact # c("dataset_drug_var", match.antiexact) 
      )
      
      calibration_data <- MatchIt::get_matches(match_model, data = initial_dataset)
      
    } else {
      
      calibration_data <- initial_dataset
      
    }
    
    # ---------------------------
    # Initialize result vectors
    # ---------------------------
    coef      <- rep(NA_real_, cg)
    coef_low  <- rep(NA_real_, cg)
    coef_high <- rep(NA_real_, cg)
    mean_vals <- rep(NA_real_, cg)
    n_drug1 <- rep(0, cg)
    n_drug2 <- rep(0, cg)
    
    # ---------------------------
    # Grouping patients
    # ---------------------------
    
    calibration_data <- calibration_data %>%
      dplyr::mutate(grouping = dplyr::ntile(dataset_benefit, cg))
    
    # ---------------------------
    # Loop through calibration groups
    # ---------------------------
    for (g in seq_len(cg)) {
      
      # Subset and factorize drug variable
      group_data <- calibration_data %>%
        dplyr::filter(grouping == g) %>%
        dplyr::mutate(dataset_drug_var = factor(dataset_drug_var, levels = rev(drugs)))
      
      # Record group size and mean benefit
      mean_vals[g] <- mean(group_data$dataset_benefit, na.rm = TRUE)
      n_drug1[g] <- nrow(group_data %>% filter(dataset_drug_var == drugs[1]))
      n_drug2[g] <- nrow(group_data %>% filter(dataset_drug_var == drugs[2]))

            
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
    result <- bind_rows(
      result,
      data.frame(
        mean     = mean_vals,
        coef     = coef,
        coef_low = coef_low,
        coef_high= coef_high,
        n_groups = cg,
        drug1    = drugs[1],
        n_drug1  = n_drug1,
        drug2    = drugs[2],
        n_drug2  = n_drug2
      )
    )
    
  }
  
  # Return result table
  return(result)
  
}

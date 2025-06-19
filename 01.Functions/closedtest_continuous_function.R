#
# Function: closed testing procedure -------------------------------------------
#

## Inputs
# cohort:         name of the cohort on which the recalibration procedure is applied
# dataset:        dataset used for recalibration
# observed:       observed outcome values
# predicted:      predicted values from original model in new data
# p_value:        nominal p_value 
# original_model: object with original model, requires to have "coefficients" information  

## Outputs
# testing_results: results with the test decision
# summary_models: list with all models ("Original", "Recalibrated intercept", "Recalibrated")
# name_chosen_model
# chosen_model


closedtest_continuous_function <- function(cohort, dataset, original_model, outcome_name, p_value){
  
  ## Some definitions ------------------------------------------------------------
  
  coefs              <- as.vector(original_model$coefficients)
  
  X_matrix_predict   <-  predict(original_model, dataset, type = "x")
  number_p_total     <- dim(X_matrix_predict)[2] - sum(as.vector(apply(X_matrix_predict, 2, sum)) == 0) - 1
  observed           <- as.vector(unlist(dataset[, outcome_name]))
  predicted          <- predict(original_model, dataset) # need to check if there is missing data (gives error throughout)
  n                  <- dim(dataset)[1]
  
  
  
  
  ## Original model --------------------------------------------------------------
  
  residuals_ormo          <- observed-predicted                                   
  sigma2_ormo             <- var(residuals_ormo)                                  
  
  # log likelihood of the original model (in new data)
  logLik_ormo             <- -n/2 * log(2 * pi * sigma2_ormo) - 1/(2 * sigma2_ormo) * sum(residuals_ormo^2)      
  
  
  
  ## Recalibration method 1: Calculate recalibration in the large (model intercept) -----
  
  model_upintercept       <- lm(observed - predicted ~ 1, data = dataset)     
  modelcoef_upintercept   <- cbind(model_upintercept$coefficients[1], confint(model_upintercept)[1], confint(model_upintercept)[2])
  
  residuals_upintmo       <- residuals(model_upintercept)
  sigma2_upintmo          <- var(residuals_upintmo)
  
  # log likelihood for the recalibrated intercept model (in new data)
  logLik_upintmo          <- -n/2 * log(2 * pi * sigma2_upintmo) - 1/(2 * sigma2_upintmo) * sum(residuals_upintmo^2)
  
  
  
  ## Recalibration method 2: Calculate coefficients after recalibration --------
  # recalibrate (overall) slope and intercept
  
  model_recalibrated           <- lm(observed ~ predicted, data = dataset)                                
  
  model_recalibrated_intercept <- cbind(model_recalibrated$coefficients[1], confint(model_recalibrated)[1,1], confint(model_recalibrated)[1,2])
  model_recalibrated_slope     <- cbind(model_recalibrated$coefficients[2], confint(model_recalibrated)[2,1], confint(model_recalibrated)[2,2])  # overall slope
  
  residuals_recalibrated       <- residuals(model_recalibrated)
  sigma2_recalibrated          <- var(residuals_recalibrated)
  
  # log likelihood for the recalibrated model (in new data)
  logLik_recalibmo             <- -n/2 * log(2 * pi * sigma2_recalibrated) - 1/(2 * sigma2_recalibrated) * sum(residuals_recalibrated^2)
  
  
  ## Test updating models  -------------------------------------------------------
  
  # Test 1: Recalibrated intercept model (recalibration in the large) against the original model using df = p (no extra coefficient estimated)
  # Test 2: If Test 1 is significant, test recalibrated model against recalibrated intercept model using df = p + 1 (1 extra coefficient is estimated)
  # --> If test 2 is significant, select re calibrated model as final model, 
  # --> if test 1 is significant but not test 2, select recalibrated intercept model
  # --> if test 1 and test 2 or not significant, select the original model 
  
  # calculation of the deviances
  deviance_test1    <- -2*logLik_ormo    + 2*logLik_upintmo
  deviance_test2    <- -2*logLik_upintmo + 2*logLik_recalibmo
  
  # decide on degrees of freedom
  df_test1          <- number_p_total
  df_test2          <- number_p_total + 1
  
  # perform testing
  desicion_test1    <- (1-pchisq(deviance_test1, df_test1)) < p_value
  desicion_test2    <- (1-pchisq(deviance_test2, df_test2)) < p_value
  
  # collect test p-values
  p_value_test1     <- (1-pchisq(deviance_test1, df_test1))
  p_value_test2     <- (1-pchisq(deviance_test2, df_test2))
  
  # choice of model based on closed testing procedure 
  choice_of_upintmo    <- 1 * (!desicion_test1)
  choice_of_recalibmo  <- 2 * ((!desicion_test1)&(!desicion_test2))
  
  index_tests          <- (choice_of_upintmo + choice_of_recalibmo)
  
  
  # collect all results from the testing procedure
  
  testing_results <- data.frame(cohort            = c(cohort,cohort,cohort),
                                n_population      = c(nrow(dataset), nrow(dataset), nrow(dataset)),
                                model             = c("Original", "Recalibrated intercept", "Recalibrated"),
                                loglikelihood     = c(logLik_ormo, logLik_upintmo, logLik_recalibmo),
                                intercept         = c(NA, modelcoef_upintercept[1], model_recalibrated_intercept[1]),
                                intercept_with_CI = c(NA,
                                                      paste0(round(modelcoef_upintercept[1],2), " (",paste0(round(modelcoef_upintercept[2],2),", ",paste0(round(modelcoef_upintercept[3],2)),")")),
                                                      paste0(round(model_recalibrated_intercept[1],2), " (",paste0(round(model_recalibrated_intercept[2],2),", ",paste0(round(model_recalibrated_intercept[3],2)),")"))),
                                
                                slope             = c(NA,NA,model_recalibrated_slope[1]),
                                slope_with_CI     = c(NA,
                                                      NA,
                                                      paste0(round(model_recalibrated_slope[1],2), " (",paste0(round(model_recalibrated_slope[2],2),", ",paste0(round(model_recalibrated_slope[3],2)),")"))),
                                df_test           = c(NA, df_test1, df_test2), 
                                test_pvalue       = c(NA, round(p_value_test1,5), round(p_value_test2,5)),
                                model_selected    = c(ifelse(desicion_test1 == FALSE & desicion_test2 == FALSE, "Yes", "No"),
                                                      ifelse(desicion_test1 == TRUE  & desicion_test2 == FALSE, "Yes", "No"),
                                                      ifelse(desicion_test2 == TRUE, "Yes", "No"))
  )
  
  
  summary_models       <- list("Original"               = original_model, 
                               "Recalibrated intercept" = model_upintercept, 
                               "Recalibrated"           = model_recalibrated)
  
  name_chosen_model    <- testing_results$model[testing_results$model_selected == "Yes"]
  chosen_model         <- summary_models[[name_chosen_model]]
  
  
  return(list(testing_results   = testing_results,
              summary_models    = summary_models, 
              name_chosen_model = name_chosen_model,
              chosen_model      = chosen_model))
  
}




#
# Function: Predictions with the model chosen by testing procedure -------------
#

## input: 
# test_results: object returned by from closedtest_continuous_function
# data: data frame for which outcome shouldbe predicted 

## output
# returns vector of predicted outcome values 


predict_with_modelchoice_function <- function(test_results, data){
  
  if (test_results$name_chosen_model == "Original") {
    
    Y_predicted <- as.vector(predict(test_results$summary_models$Original, data))
    
    Y_out <- Y_predicted
    
  } else if (test_results$name_chosen_model == "Recalibrated intercept") {
    
    Y_predicted <- as.vector(predict(test_results$summary_models$Original, data))
    
    Y_out <- Y_predicted + as.numeric(test_results$chosen_model$coefficients["(Intercept)"])
    
  } else if (test_results$name_chosen_model == "Recalibrated") {
    
    Y_predicted <- as.vector(predict(test_results$summary_models$Original, data))
    
    Y_out <- (as.numeric(test_results$chosen_model$coefficients["predicted"]) * Y_predicted) + as.numeric(test_results$chosen_model$coefficients["(Intercept)"])
    
  } else {
    
    stop("This hasn't been coded, the function should not be here.")
    
  }
  
  return(Y_out)
  
}



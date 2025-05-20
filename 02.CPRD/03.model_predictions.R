# # load libraries
# library(tidyverse)
# library(aurum)
# 
# 
# # set up aurum
# cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
# analysis = cprd$analysis("pedro_mm")
# 
# # set up functions
# is.integer64 <- function(x){
#   class(x)=="integer64"
# }
# source("02.CPRD/02.impute_missingness.R")
# 
# # set up dataset
# analysis_post_2020_raw <- analysis_post_2020_raw %>%
#   analysis$cached("analysis_post_2020") %>%
#   collect() %>%
#   mutate(patid=as.character(patid)) %>%
#   mutate_if(is.integer64, as.integer)
# 
# # load model
# load("fivedrugmodel_5knot_share_20230823.Rdata")
# 
# #:------------------------------------------------------
# # Pre-processing datasets
# #:-- Post 2020-10-14
# analysis_post_2020 <- analysis_post_2020_raw %>%
#   mutate(
#     pated = paste(patid, dstartdate, drug_class, sep = "."),
#     sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
#     agetx = dstartdate_age,
#     ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
#     ethnicity = factor(ethnicity, levels = c(0:5), labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),
# 
#     smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),
#     smoke = factor(smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),
#     imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
#     imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5), labels = c("1 (least)", "2", "3", "4", "5 (most)")),
# 
#     ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
#     ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
#     ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),
#     drugline = ifelse(drugline_all > 5, 5, drugline_all),
#     drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),
#     drugclass = drug_class
#   ) %>%
#   group_by(pated) %>%
#   mutate(row = 1:n()) %>%
#   ungroup() %>%
#   filter(row == 1) %>%
#   select(-row) %>%
#   select(all_of(c(
#       "pated", "agetx", "sex", "t2dmduration", "ethnicity", 
#       "drug_substance", "drug_class",
#       "imd5", "smoke",
#       "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
#       "drugline", "ncurrtx", "hba1cmonth",
#       "posthba1cfinal"
#     )
#   )) %>%
#   rename("drugclass" = "drug_class")
# 
# 
# 
# 
# #:------------------------------------------------------
# # Impute missing data
# source("02.CPRD/02.impute_missingness.R")
# analysis_post_2020_missingness <- imputation_methods(data = analysis_post_2020, 
#                                                      method = "mice", 
#                                                      mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))


#:------------------------------------------------------
# Predict from model
predict_5drugmodel <- function(data, model, drug_var = NULL, drugs = NULL, pred_col = "pred.") {
  
  # check whether objects is provided
  if (is.null(drugs)) {stop("drugs needs to be provided")}
  if (is.null(drug_var)) {stop("drug_var needs to be provided")}
  
  # check class of object (method used if class == "imputation method")
  if ("imputation method" %in% class(data)) {
    
    # check the imputation method
    if (data$type == "group") {
      # if group used
      
      # output object
      output_dataset <- data$imputation
      
      # check whether drugs provided are included in the dataset
      if (!all(drugs %in% unique(output_dataset %>% select(all_of(drug_var)) %>% unlist() %>% unique()))) {stop("some drugs provided are not in drug_var")}
      if (is.null(drugs)) {stop("drugs needs to be provided")}
      
      # name of variables imputated
      cols_imputed <- output_dataset %>% select(contains("_impute")) %>% colnames()
      orig_cols_imputed <- sub("_group_impute$", "", cols_imputed)
      
      # remove the original vars with missingness
      imputed_dataset <- output_dataset %>%
        select(-all_of(orig_cols_imputed))
      
      # rename imputed vars
      colnames(imputed_dataset) <- sub("_group_impute$", "", colnames(imputed_dataset))
      
      # predict from model
      for (drug_name in drugs) {
        
        # change drug_var to the specific drug being predicted
        current_prediction_dataset <- imputed_dataset %>%
          rename("target_var" = drug_var) %>%
          mutate(target_var = drug_name) %>%
          rename({{drug_var}} := "target_var")
        
        # predict from model
        prediction_vector <- predict(model, current_prediction_dataset)
        
        # attach this var to original dataset
        output_dataset <- cbind(output_dataset, prediction_vector) %>%
          as.data.frame() %>%
          rename(!!paste0("pred.group.", drug_name) := "prediction_vector")
        
      }
      
      
    } else if (data$type == "mice") {
      # if mice used
      
      # output object
      output_dataset <- data$imputation
      
      # check whether drugs provided are included in the dataset
      if (!all(drugs %in% unique(output_dataset %>% select(all_of(drug_var)) %>% unlist() %>% unique()))) {stop("some drugs provided are not in drug_var")}
      if (is.null(drugs)) {stop("drugs needs to be provided")}
      
      # name of variables imputated
      cols_imputed <- output_dataset %>% select(contains("_mice_impute")) %>% colnames()
      orig_cols_imputed <- sub("_mice_impute$", "", cols_imputed)
      
      # remove the original vars with missingness
      imputed_dataset <- output_dataset %>%
        select(-all_of(orig_cols_imputed))
      
      # rename imputed vars
      colnames(imputed_dataset) <- sub("_mice_impute$", "", colnames(imputed_dataset))
      
      # predict from model
      for (drug_name in drugs) {
        
        # change drug_var to the specific drug being predicted
        current_prediction_dataset <- imputed_dataset %>%
          rename("target_var" = drug_var) %>%
          mutate(target_var = drug_name) %>%
          rename({{drug_var}} := "target_var")
        
        # predict from model
        prediction_vector <- predict(model, current_prediction_dataset)
        
        # attach this var to original dataset
        output_dataset <- cbind(output_dataset, prediction_vector) %>%
          as.data.frame() %>%
          rename(!!paste0("pred.mice.", drug_name) := "prediction_vector")
        
      }
      
    } else if (data$type == "dpmm") {
      # if DPMM used
      stop("DPMM not codded yet")
    } else {stop("This is an error, you shouldn't be here")}
    
  } else {
    # no imputation method / prediction from model
    
    # output object
    output_dataset <- data
    
    # predict from model
    for (drug_name in drugs) {
      
      # change drug_var to the specific drug being predicted
      current_prediction_dataset <- output_dataset %>%
        rename("target_var" = drug_var) %>%
        mutate(target_var = drug_name) %>%
        rename({{drug_var}} := "target_var")
      
      # predict from model
      prediction_vector <- predict(model, current_prediction_dataset)
      
      # attach this var to original dataset
      output_dataset <- cbind(output_dataset, prediction_vector) %>%
        as.data.frame() %>%
        rename(!!paste0(pred_col, drug_name) := "prediction_vector")
      
    }
    
  }
  
  return(output_dataset)
  
}


# library(rms)
# analysis_post_2020_prediction <- predict_5drugmodel(analysis_post_2020_missingness, 
#                                                     model = m1.5.final, 
#                                                     drug_var = "drugclass",
#                                                     drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))
# 
# analysis_post_2020_prediction <- predict_5drugmodel(analysis_post_2020_prediction,
#                                                     model = m1.5.final,
#                                                     drug_var = "drugclass",
#                                                     drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"),
#                                                     pred_col = "pred.orig.")













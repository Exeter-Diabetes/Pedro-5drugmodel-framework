# # load libraries
# library(mice)
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
# 
# # set up dataset
# analysis_post_2020_raw <- analysis_post_2020_raw %>%
#   analysis$cached("analysis_post_2020") %>%
#   collect() %>%
#   mutate(patid=as.character(patid)) %>%
#   mutate_if(is.integer64, as.integer)
# 
# 
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
#   select(-row)


# Function for different imputation methods
imputation_methods <- function(data, method, mice.ignore.vars = NULL, mice.m = 5, mice.maxit = 10) {
  
  # load libraries
  require(tidyverse)
  
  # check whether method is one of the available options
  if (!(tolower(method) %in% c("group", "mice", "dpmm"))) {stop("The only methods available are: group / mice / DPMM")}
  
  
  #:------------------------------------------------------
  # Imputation methods
  
  if (tolower(method) == "group") {
    
    ## Group patients by selected variables
    ### sex, ethnicity, hba1c by 10 mmol/mol
    ## Impute with mean value
    
    imputation <- data %>%
      mutate(
        hba1c_groups = cut(prehba1c, breaks = seq(50, 110, 10), include.lowest = TRUE)
      ) %>%
      group_by(sex, ethnicity, hba1c_groups) %>% # reomved age (?)
      mutate(
        prebmi_impute = mean(prebmi, na.rm = TRUE),
        preegfr_impute = mean(preegfr, na.rm = TRUE),
        pretotalcholesterol_impute = mean(pretotalcholesterol, na.rm = TRUE),
        prehdl_impute = mean(prehdl, na.rm = TRUE),
        prealt_impute = mean(prealt, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(
        prebmi_group_impute = ifelse(is.na(prebmi), prebmi_impute, prebmi),
        preegfr_group_impute = ifelse(is.na(preegfr), preegfr_impute, preegfr),
        pretotalcholesterol_group_impute = ifelse(is.na(pretotalcholesterol), pretotalcholesterol_impute, pretotalcholesterol),
        prehdl_group_impute = ifelse(is.na(prehdl), prehdl_impute, prehdl),
        prealt_group_impute = ifelse(is.na(prealt), prealt_impute, prealt)
      ) %>%
      select(-c(prebmi_impute, preegfr_impute, pretotalcholesterol_impute, prehdl_impute, prealt_impute)) %>%
      as.data.frame()
    
    # output object
    output <- list(
      type = "group",
      imputation = imputation
    )
    
  } else if (tolower(method) == "mice") {
    
    # MICE

    # load libraries
    require(mice)
    
    if (!is.null(mice.ignore.vars)) {
      
      prediction_matrix <- mice::mice(
        data = data,
        maxit = 0
      )$predictorMatrix
      
      # change ignore variables to not be used
      prediction_matrix[mice.ignore.vars,] <- 0

    }
    
    # fit mice model and impute
    mice_output <- mice::mice(
      data = data,
      m = mice.m,
      maxit = mice.maxit,
      predictorMatrix = prediction_matrix
    )
    
    # datasets with imputations
    imputed_datasets <- complete(mice_output, action = "all")
    
    # list of imputed variables
    selected_list <- lapply(imputed_datasets, \(df) select(df, all_of(names(Filter(\(x) x != "", mice_output$method)))))
    
    # get average of imputed values
    imputed_variables <- Reduce("+", selected_list) / mice_output$m
    colnames(imputed_variables) <- paste0(colnames(imputed_variables), "_mice_impute")
    
    # final dataset
    imputation <- data %>%
      cbind(imputed_variables) %>%
      as.data.frame ()
    
    # output object
    output <- list(
      type = "mice",
      imputation = imputation
    )
    
  } else if (tolower(method) == "dpmm") {
    
    # DPMM
    ## Fit the model (complete data)
    
    ## Make the predictions
    
    stop("DPMM not coded yet")
    
  } else {stop("This is an error, the function shouldn't be here.")}
  
  class(output) = "imputation method"
  
  return(output)
  
}


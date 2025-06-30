
# Predict from model ----
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
      # if (!all(drugs %in% unique(output_dataset %>% select(all_of(drug_var)) %>% unlist() %>% unique()))) {stop("some drugs provided are not in drug_var")}
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

# 
# # Initial set-up ###########################################
# 
# # load libraries
# library(tidyverse)
# library(aurum)
# library(rms)
# 
# # set up aurum
# cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
# analysis = cprd$analysis("pedro_mm")
# 
# ## functions ----
# is.integer64 <- function(x){
#   class(x)=="integer64"
# }
# 
# sapply(
#   paste0("01.Functions/", list.files("01.Functions", pattern = "\\.R$")),
#   source
# )
# source("02.CPRD/03.impute_missingness.R")
# 
# # load model
# load("fivedrugmodel_5knot_share_20230823.Rdata")
# 
# # set up dataset
# analysis_post_2020_raw <- analysis_post_2020_raw %>%
#   analysis$cached("analysis_post_2020") %>%
#   collect() %>%
#   mutate(patid=as.character(patid)) %>%
#   mutate_if(is.integer64, as.integer)
# 
# 
# # Pre-processing datasets ########################################
# 
# ## Post 2020-10-14 ----
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
#     "pated", "agetx", "sex", "t2dmduration", "ethnicity",
#     "drug_substance", "drug_class",
#     "imd5", "smoke",
#     "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
#     "drugline", "ncurrtx", "hba1cmonth",
#     "posthba1cfinal"
#   )
#   )) %>%
#   rename("drugclass" = "drug_class") %>%
#   as.data.frame()
# 
# 
# # Missing data imputation ######################################
# 
# ## Post 2020-10-14 ----
# ### By subgroups
# analysis_post_2020_group_imputation <- imputation_methods(data = analysis_post_2020,
#                                                           method = "group")
# 
# ### By mice
# analysis_post_2020_mice_imputation <- imputation_methods(data = analysis_post_2020,
#                                                          method = "mice",
#                                                          mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))
# 
# 
# # Predictions from datasets #########################################
# 
# ## Post 2020-10-14 ----
# ### Original variables
# analysis_post_2020_prediction_orig <- predict_5drugmodel(analysis_post_2020,
#                                                          model = m1.5.final,
#                                                          drug_var = "drugclass",
#                                                          drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"),
#                                                          pred_col = "pred.orig.")
# 
# ### Imputation by groups
# analysis_post_2020_prediction_group <- predict_5drugmodel(analysis_post_2020_group_imputation,
#                                                           model = m1.5.final,
#                                                           drug_var = "drugclass",
#                                                           drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))
# 
# ### Imputation by mice
# analysis_post_2020_prediction_mice <- predict_5drugmodel(analysis_post_2020_mice_imputation,
#                                                          model = m1.5.final,
#                                                          drug_var = "drugclass",
#                                                          drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))
# 
# 
# ### merge impute columns into main dataset
# analysis_post_2020 <- analysis_post_2020 %>%
#   cbind(
#     analysis_post_2020_prediction_orig %>%
#       select(contains("pred.orig")),
#     analysis_post_2020_prediction_group %>%
#       select(contains("group_impute")),
#     analysis_post_2020_prediction_group %>%
#       select(contains("pred.group")),
#     analysis_post_2020_prediction_mice %>%
#       select(contains("mice_impute")),
#     analysis_post_2020_prediction_mice %>%
#       select(contains("pred.mice"))
#   )
# 
# 
# 
# # Closed loop test #################################################
# 
# ## Post 2020-10-14 ----
# 
# ### Original variables
# #### SGLT2
# closed_loop_test_results_SGLT2_post_2020_orig <- closedtest_continuous_function(
#   cohort = "SGLT2 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "SGLT2") %>% drop_na(),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### GLP1
# closed_loop_test_results_GLP1_post_2020_orig <- closedtest_continuous_function(
#   cohort = "GLP1 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "GLP1") %>% drop_na(),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### DPP4
# closed_loop_test_results_DPP4_post_2020_orig <- closedtest_continuous_function(
#   cohort = "DPP4 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "DPP4") %>% drop_na(),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### TZD
# closed_loop_test_results_TZD_post_2020_orig <- closedtest_continuous_function(
#   cohort = "TZD subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "TZD") %>% drop_na(),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### SU
# closed_loop_test_results_SU_post_2020_orig <- closedtest_continuous_function(
#   cohort = "SU subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "SU") %>% drop_na(),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# 
# ### group variables
# #### SGLT2
# closed_loop_test_results_SGLT2_post_2020_group <- closedtest_continuous_function(
#   cohort = "SGLT2 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "SGLT2") %>%
#     mutate(
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### GLP1
# closed_loop_test_results_GLP1_post_2020_group <- closedtest_continuous_function(
#   cohort = "GLP1 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "GLP1") %>%
#     mutate(
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### DPP4
# closed_loop_test_results_DPP4_post_2020_group <- closedtest_continuous_function(
#   cohort = "DPP4 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "DPP4") %>%
#     mutate(
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### TZD
# closed_loop_test_results_TZD_post_2020_group <- closedtest_continuous_function(
#   cohort = "TZD subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "TZD") %>%
#     mutate(
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### SU
# closed_loop_test_results_SU_post_2020_group <- closedtest_continuous_function(
#   cohort = "SU subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "SU") %>%
#     mutate(
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# 
# 
# ### mice variables
# #### SGLT2
# closed_loop_test_results_SGLT2_post_2020_mice <- closedtest_continuous_function(
#   cohort = "SGLT2 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "SGLT2") %>%
#     mutate(
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### GLP1
# closed_loop_test_results_GLP1_post_2020_mice <- closedtest_continuous_function(
#   cohort = "GLP1 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "GLP1") %>%
#     mutate(
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### DPP4
# closed_loop_test_results_DPP4_post_2020_mice <- closedtest_continuous_function(
#   cohort = "DPP4 subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "DPP4") %>%
#     mutate(
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### TZD
# closed_loop_test_results_TZD_post_2020_mice <- closedtest_continuous_function(
#   cohort = "TZD subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "TZD") %>%
#     mutate(
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# #### SU
# closed_loop_test_results_SU_post_2020_mice <- closedtest_continuous_function(
#   cohort = "SU subcohort",
#   dataset = analysis_post_2020 %>% filter(drugclass == "SU") %>%
#     mutate(
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     ),
#   original_model = m1.5.final,
#   outcome_name = "posthba1cfinal",
#   p_value = 0.05
# )
# 
# 
# 
# 
# 
# 
# #### Make predictions ----
# analysis_post_2020 <- analysis_post_2020 %>%
#   mutate(
#     # original variables
#     pred.orig.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "SGLT2")),
#     pred.orig.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "GLP1")),
#     pred.orig.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "DPP4")),
#     pred.orig.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_orig, analysis_post_2020 %>% mutate(drugclass = "TZD")),
#     pred.orig.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_orig, analysis_post_2020 %>% mutate(drugclass == "SU")),
#     # group variables
#     pred.group.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_group, analysis_post_2020 %>% mutate(
#       drugclass = "SGLT2",
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     )),
#     pred.group.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_group, analysis_post_2020 %>% mutate(
#       drugclass = "GLP1",
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     )),
#     pred.group.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_group, analysis_post_2020 %>% mutate(
#       drugclass = "DPP4",
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     )),
#     pred.group.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_group, analysis_post_2020 %>% mutate(
#       drugclass = "TZD",
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     )),
#     pred.group.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_group, analysis_post_2020 %>% mutate(
#       drugclass = "SU",
#       prebmi = prebmi_group_impute,
#       preegfr = preegfr_group_impute,
#       pretotalcholesterol = pretotalcholesterol_group_impute,
#       prehdl = prehdl_group_impute,
#       prealt = prealt_group_impute
#     )),
#     # mice variables
#     pred.mice.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_post_2020_mice, analysis_post_2020 %>% mutate(
#       drugclass = "SGLT2",
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     )),
#     pred.mice.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_post_2020_mice, analysis_post_2020 %>% mutate(
#       drugclass = "GLP1",
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     )),
#     pred.mice.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_post_2020_mice, analysis_post_2020 %>% mutate(
#       drugclass = "DPP4",
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     )),
#     pred.mice.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_post_2020_mice, analysis_post_2020 %>% mutate(
#       drugclass = "TZD",
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     )),
#     pred.mice.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_post_2020_mice, analysis_post_2020 %>% mutate(
#       drugclass = "SU",
#       prebmi = prebmi_mice_impute,
#       preegfr = preegfr_mice_impute,
#       pretotalcholesterol = pretotalcholesterol_mice_impute,
#       prehdl = prehdl_mice_impute,
#       prealt = prealt_mice_impute
#     ))
#   )
# 





# Author: pcardoso
###############################################################################

# Per drug pair benefit calibration

###############################################################################

# load libraries
library(tableone)
# remotes::install_github("harrelfe/rms")
library(rms)
library(qpdf)

# load function for generating dataset
source("/home/pcardoso/workspace/Pedro-5drugmodel-framework/03.SDRN/01.cohort_script.R")
# load functions needed
sapply(paste0("/home/pcardoso/workspace/Pedro-5drugmodel-framework/01.Functions/", list.files("/home/pcardoso/workspace/Pedro-5drugmodel-framework/01.Functions/", pattern = "\\.R$")),
		source
)
source("/home/pcardoso/workspace/Pedro-5drugmodel-framework/03.SDRN/03.impute_missingness.R")
source("/home/pcardoso/workspace/Pedro-5drugmodel-framework/03.SDRN/04.model_predictions.R")

# load dataset
analysis_cohort_raw <- set_up_data(data = "mm_20250506_t2d_1stinstance")

# load model
load("/home/pcardoso/workspace/Pedro-5drugmodel-framework/fivedrugmodel_5knot_share_20230823.Rdata")


# Pre-processing datasets #########################################

## analysis cohort
analysis_cohort <- analysis_cohort_raw %>%
		mutate(
				pated = paste(serialno, dstartdate, drug_class, sep = "."),
				sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
				agetx = dstartdate_age,
				ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
				ethnicity = factor(ethnicity, levels = c(0:5), labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),
				
				smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),
				smoke = factor(smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),
				imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
				imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5), labels = c("1 (least)", "2", "3", "4", "5 (most)")),
				
				ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
				ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
				ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),
				drugline = ifelse(drugline_all > 5, 5, drugline_all),
				drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),
				drugclass = drug_class
		) %>%
		group_by(pated) %>% 
		mutate(row = 1:n()) %>% 
		ungroup() %>% 
		filter(row == 1) %>%
		select(-row) %>%
		select(all_of(c(
								"pated", "agetx", "sex", "t2dmduration", "ethnicity",
								"drug_substance", "drug_class",
								"imd5", "smoke",
								"prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
								"drugline" ,"ncurrtx", "hba1cmonth",
								"posthba1cfinal"
						))) %>%
		rename("drugclass" = "drug_class") %>%
		as.data.frame()


# Missing data imputation ----

## Standard drugs ----
### By subgroups
analysis_cohort_mice_imputation <- imputation_methods(data = analysis_cohort,
		method = "mice",
		mice.ignore.vars = c("pated", "drug_substance", "drugclass", "hba1cmonth", "posthba1cfinal"))

# Predictions from datasets ----

## Standard drugs ----
analysis_cohort_prediction_mice <- predict_5drugmodel(analysis_cohort_mice_imputation,
		model = m1.5.final,
		drug_var = "drugclass",
		drugs = c("SGLT2", "GLP1", "DPP4", "TZD", "SU"))

### merge imputed columns into main dataset
analysis_cohort <- analysis_cohort %>%
		cbind(
				analysis_cohort_prediction_mice %>%
						select(contains("mice_impute")),
				analysis_cohort_prediction_mice %>%
						select(contains("pred.mice")) %>%
						rename_with(~ str_replace(., "pred\\.mice\\.", "pred.mice.preclosed."))
		)


# Closed loop test ----

## Standard drugs ----

### mice variables
#### SGLT2
closed_loop_test_results_SGLT2_standard_mice <- closedtest_continuous_function(
		cohort = "SGLT2 subcohort",
		dataset = analysis_cohort %>%
				filter(!(drug_substance %in% c("Low-dose semaglutide", "Oral semaglutide"))) %>% 
				filter(drugclass == "SGLT2") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)


#### GLP1
closed_loop_test_results_GLP1_standard_mice <- closedtest_continuous_function(
		cohort = "GLP1 subcohort",
		dataset = analysis_cohort %>%
				filter(!(drug_substance %in% c("Low-dose semaglutide", "Oral semaglutide"))) %>% 
				filter(drugclass == "GLP1") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)


#### DPP4
closed_loop_test_results_DPP4_standard_mice <- closedtest_continuous_function(
		cohort = "DPP4 subcohort",
		dataset = analysis_cohort %>%
				filter(!(drug_substance %in% c("Low-dose semaglutide", "Oral semaglutide"))) %>% 
				filter(drugclass == "DPP4") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)



#### TZD
closed_loop_test_results_TZD_standard_mice <- closedtest_continuous_function(
		cohort = "TZD subcohort",
		dataset = analysis_cohort %>%
				filter(!(drug_substance %in% c("Low-dose semaglutide", "Oral semaglutide"))) %>% 
				filter(drugclass == "TZD") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)


#### SU
closed_loop_test_results_SU_standard_mice <- closedtest_continuous_function(
		cohort = "SU subcohort",
		dataset = analysis_cohort %>%
				filter(!(drug_substance %in% c("Low-dose semaglutide", "Oral semaglutide"))) %>% 
				filter(drugclass == "SU") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)


## Injectable semaglutide ----

### mice variables
#### GLP1
closed_loop_test_results_GLP1_semaglutide_mice <- closedtest_continuous_function(
		cohort = "Semaglutide subcohort",
		dataset = analysis_cohort %>%
				filter(drug_substance == "Low-dose semaglutide") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)


## Oral semaglutide ----

### mice variables
#### GLP1
closed_loop_test_results_GLP1_oral_semaglutide_mice <- closedtest_continuous_function(
		cohort = "Oral semaglutide subcohort",
		dataset = analysis_cohort %>%
				filter(drug_substance == "Oral semaglutide") %>%
				mutate(
						prebmi = prebmi_mice_impute,
						preegfr = preegfr_mice_impute,
						pretotalcholesterol = pretotalcholesterol_mice_impute,
						prehdl = prehdl_mice_impute,
						prealt = prealt_mice_impute
				),
		original_model = m1.5.final,
		outcome_name = "posthba1cfinal",
		p_value = 0.05
)


## Make predictions ----
analysis_standard <- analysis_cohort %>%
		filter(!(drug_substance %in% c("Low-dose semaglutide", "Oral semaglutide"))) %>%
		mutate(
				prebmi = prebmi_mice_impute,
				preegfr = preegfr_mice_impute,
				pretotalcholesterol = pretotalcholesterol_mice_impute,
				prehdl = prehdl_mice_impute,
				prealt = prealt_mice_impute
		)

analysis_standard <- analysis_standard %>%
		mutate(
				pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_standard_mice, analysis_standard %>% 
								mutate(
										drugclass = "SGLT2"
								)
				),
				pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_standard_mice, analysis_standard %>% 
								mutate(
										drugclass = "GLP1"
								)
				),
				pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_standard_mice, analysis_standard %>% 
								mutate(
										drugclass = "DPP4"
								)
				),
				pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_standard_mice, analysis_standard %>% 
								mutate(
										drugclass = "TZD"
								)
				),
				pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_standard_mice, analysis_standard %>% 
								mutate(
										drugclass = "SU"
								)
				),
				pred.Sema = predict_with_modelchoice_function(closed_loop_test_results_GLP1_semaglutide_mice, analysis_standard %>% 
								mutate(
										drugclass = "GLP1"
								)
				),
				pred.Oral = predict_with_modelchoice_function(closed_loop_test_results_GLP1_oral_semaglutide_mice, analysis_standard %>% 
								mutate(
										drugclass = "GLP1"
								)
				)
		)


## Injectable semaglutide ----
analysis_semaglutide <- analysis_cohort %>%
		filter(drug_substance == "Low-dose semaglutide") %>%
		mutate(
			prebmi = prebmi_mice_impute,
			preegfr = preegfr_mice_impute,
			pretotalcholesterol = pretotalcholesterol_mice_impute,
			prehdl = prehdl_mice_impute,
			prealt = prealt_mice_impute
		)


analysis_semaglutide <- analysis_semaglutide %>%
		mutate(
				pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_standard_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "SGLT2"
								)
				),
				pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_standard_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "GLP1"
								)
				),
				pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_standard_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "DPP4"
								)
				),
				pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_standard_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "TZD"
								)
				),
				pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_standard_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "SU"
								)
				),
				pred.Sema = predict_with_modelchoice_function(closed_loop_test_results_GLP1_semaglutide_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "GLP1"
								)
				),
				pred.Oral = predict_with_modelchoice_function(closed_loop_test_results_GLP1_oral_semaglutide_mice, analysis_semaglutide %>% 
								mutate(
										drugclass = "GLP1"
								)
				)
		)

## Oral semaglutide ----
analysis_oral_semaglutide <- analysis_cohort %>%
		filter(drug_substance == "Oral semaglutide") %>%
		mutate(
				prebmi = prebmi_mice_impute,
				preegfr = preegfr_mice_impute,
				pretotalcholesterol = pretotalcholesterol_mice_impute,
				prehdl = prehdl_mice_impute,
				prealt = prealt_mice_impute
		)

analysis_oral_semaglutide <- analysis_oral_semaglutide %>%
		mutate(
				pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2_standard_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "SGLT2"
								)
				),
				pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1_standard_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "GLP1"
								)
				),
				pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4_standard_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "DPP4"
								)
				),
				pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD_standard_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "TZD"
								)
				),
				pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU_standard_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "SU"
								)
				),
				pred.Sema = predict_with_modelchoice_function(closed_loop_test_results_GLP1_semaglutide_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "GLP1"
								)
				),
				pred.Oral = predict_with_modelchoice_function(closed_loop_test_results_GLP1_oral_semaglutide_mice, analysis_oral_semaglutide %>% 
								mutate(
										drugclass = "GLP1"
								)
				)
		)

# Unified validation ----

## Standard & Oral
analysis_standard_oral_calibration_adj <- unified_validation(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide),
		drug_var = "drugclass",
		drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
		prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
		outcome_var = "posthba1cfinal",
		cal_groups = c(3, 5, 10),
		adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

analysis_standard_oral_calibration_matching <- unified_validation(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide) %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
		prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
		outcome_var = "posthba1cfinal",
		cal_groups = c(3, 5, 10),
		matching = TRUE,
		adjustment_var = NULL,
		matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

analysis_standard_oral_calibration_matching_adj <- unified_validation(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide) %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"),
		prediction_vars = paste0("pred.", c("SGLT2", "GLP1", "TZD", "SU", "DPP4")),
		outcome_var = "posthba1cfinal",
		cal_groups = c(3, 5, 10),
		matching = TRUE,
		adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)



interim_dataset <- analysis_standard %>%
		rbind(analysis_oral_semaglutide) %>%
		rename("drugclass_old" = "drugclass") %>%
		mutate(drugclass = drugclass_old) %>%
		rbind(
			analysis_semaglutide %>%
					rename("drugclass_old" = "drugclass") %>%
					mutate(drugclass = "Semaglutide")
		)



analysis_semaglutide_calibration_adj <- unified_validation(
		data = interim_dataset,
		drug_var = "drugclass",
		drugs = c("Semaglutide", "SGLT2", "GLP1", "TZD", "SU", "DPP4"),
		prediction_vars = paste0("pred.", c("Sema", "SGLT2", "GLP1", "TZD", "SU", "DPP4")),
		outcome_var = "posthba1cfinal",
		cal_groups = c(3, 5, 10),
		adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline")
)

analysis_semaglutide_calibration_matching <- unified_validation(
		data = interim_dataset %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		drugs = c("Semaglutide", "SGLT2", "GLP1", "TZD", "SU", "DPP4"),
		prediction_vars = paste0("pred.", c("Sema", "SGLT2", "GLP1", "TZD", "SU", "DPP4")),
		outcome_var = "posthba1cfinal",
		cal_groups = c(3, 5, 10),
		matching = TRUE,
		adjustment_var = NULL,
		matching_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

analysis_semaglutide_calibration_matching_adj <- unified_validation(
		data = interim_dataset %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		drugs = c("Semaglutide", "SGLT2", "GLP1", "TZD", "SU", "DPP4"),
		prediction_vars = paste0("pred.", c("Sema", "SGLT2", "GLP1", "TZD", "SU", "DPP4")),
		outcome_var = "posthba1cfinal",
		cal_groups = c(3, 5, 10),
		matching = TRUE,
		adjustment_var = c("t2dmduration", "prebmi", "prehba1c", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "sex", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)


# Plots ----

## Standard drugs & oral
pdf("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/06.standard_oral_overall_calibration.pdf", width = 12, height = 5)
analysis_standard_oral_calibration_adj %>%
		mutate(Method = "Adjustment") %>%
		rbind(
				analysis_standard_oral_calibration_matching %>%
						mutate(Method = "Matching"),
				analysis_standard_oral_calibration_matching_adj %>%
						mutate(Method = "Matching + Adj")
		) %>%
		mutate(drugcombo = paste(drug1, drug2)) %>%
		group_by(Method, drugcombo, n_groups) %>%
		mutate(min_val = min(n_drug1, n_drug2)) %>%
		ungroup(n_groups) %>%
		mutate(
				select_grouping = ifelse(min_val > 100, n_groups, NA),
				select_grouping = max(select_grouping, na.rm = TRUE),
				select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
		) %>%
		ungroup() %>%
		filter(select_grouping == n_groups) %>%
		select(-c(drugcombo, min_val, select_grouping)) %>%
		mutate(title = paste(drug1, "vs", drug2)) %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		facet_wrap(~title, nrow = 2) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")


analysis_standard_oral_calibration_adj %>%
		mutate(drugcombo = paste(drug1, drug2)) %>%
		group_by(drugcombo, n_groups) %>%
		mutate(min_val = min(n_drug1, n_drug2)) %>%
		ungroup(n_groups) %>%
		mutate(
				select_grouping = ifelse(min_val > 100, n_groups, NA),
				select_grouping = max(select_grouping, na.rm = TRUE),
				select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
		) %>%
		ungroup() %>%
		filter(select_grouping == n_groups) %>%
		select(-c(drugcombo, min_val, select_grouping)) %>%
		mutate(title = paste(drug1, "vs", drug2)) %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		facet_wrap(~title, nrow = 2) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")

dev.off()

## Semaglutide
pdf("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/06.semaglutide_overall_calibration.pdf", width = 12, height = 5)
analysis_semaglutide_calibration_adj %>%
		mutate(Method = "Adjustment") %>%
		rbind(
				analysis_semaglutide_calibration_matching %>%
						mutate(Method = "Matching"),
				analysis_semaglutide_calibration_matching_adj %>%
						mutate(Method = "Matching + Adj")
		) %>%
		mutate(drugcombo = paste(drug1, drug2)) %>%
		group_by(Method, drugcombo, n_groups) %>%
		mutate(min_val = min(n_drug1, n_drug2)) %>%
		ungroup(n_groups) %>%
		mutate(
				select_grouping = ifelse(min_val > 100, n_groups, NA),
				select_grouping = max(select_grouping, na.rm = TRUE),
				select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
		) %>%
		ungroup() %>%
		filter(select_grouping == n_groups) %>%
		select(-c(drugcombo, min_val, select_grouping)) %>%
		mutate(title = paste(drug1, "vs", drug2)) %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high, colour = Method)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		facet_wrap(~title, nrow = 2) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")


analysis_semaglutide_calibration_adj %>%
		mutate(drugcombo = paste(drug1, drug2)) %>%
		group_by(drugcombo, n_groups) %>%
		mutate(min_val = min(n_drug1, n_drug2)) %>%
		ungroup(n_groups) %>%
		mutate(
				select_grouping = ifelse(min_val > 100, n_groups, NA),
				select_grouping = max(select_grouping, na.rm = TRUE),
				select_grouping = ifelse(is.infinite(abs(select_grouping)), 3, select_grouping)
		) %>%
		ungroup() %>%
		filter(select_grouping == n_groups) %>%
		select(-c(drugcombo, min_val, select_grouping)) %>%
		mutate(title = paste(drug1, "vs", drug2)) %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		facet_wrap(~title, nrow = 2) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Groups have at least 100 individuals on either drug")


dev.off()




		
		
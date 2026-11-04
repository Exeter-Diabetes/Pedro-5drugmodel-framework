# Author: pcardoso
###############################################################################

# Overall benefit calibration

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

# Calibration of predictions ----

## Standard drugs ----
plot_pred_response_analysis_standard_preclosed <- analysis_standard %>%
		rbind(analysis_oral_semaglutide) %>%
		pivot_longer(cols = contains("pred.mice.preclosed")) %>%
		select(name, value, drugclass, posthba1cfinal) %>%
		mutate(name = gsub("pred\\.mice\\.preclosed\\.", "", name)) %>%
		filter(drugclass == name) %>%
		ggplot(aes(x = value, y = posthba1cfinal)) +
		geom_abline(aes(intercept = 0, slope = 1)) +
		stat_smooth(method = "lm", formula = y~poly(x, 2)) +
		facet_wrap(~name) +
		labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c (mmol/mol)", title = "Standard drugs")

plot_pred_response_analysis_standard_dpp4 <- analysis_standard %>%
		filter(drugclass == "DPP4") %>%
		select("obs" = posthba1cfinal, "Original" = pred.DPP4, "Adjusted" = pred.mice.preclosed.DPP4) %>%
		gather(key = "Predictions", value = "pred", -obs) %>%
		ggplot(aes(x = pred, y = obs, colour = Predictions)) +
		geom_abline(aes(intercept = 0, slope = 1)) +
		stat_smooth(method = "lm", formula = y~poly(x, 2)) +
		labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c (mmol/mol)", title = "DPP4") +
		theme_minimal() +
		theme(
				legend.position = "bottom"
		)

plot_pred_response_analysis_standard_tzd <- analysis_standard %>%
		filter(drugclass == "TZD") %>%
		select("obs" = posthba1cfinal, "Original" = pred.TZD, "Adjusted" = pred.mice.preclosed.TZD) %>%
		gather(key = "Predictions", value = "pred", -obs) %>%
		ggplot(aes(x = pred, y = obs, colour = Predictions)) +
		geom_abline(aes(intercept = 0, slope = 1)) +
		stat_smooth(method = "lm", formula = y~poly(x, 2)) +
		labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c (mmol/mol)", title = "TZD") +
		theme_minimal() +
		theme(
				legend.position = "bottom"
		)

plot_pred_response_analysis_semaglutide <- analysis_semaglutide %>%
		select("obs" = posthba1cfinal, "Original" = pred.GLP1, "Injectable semaglutide" = pred.Sema) %>%
		gather(key = "Predictions", value = "pred", -obs) %>%
		ggplot(aes(x = pred, y = obs, colour = Predictions)) +
		geom_abline(aes(intercept = 0, slope = 1)) +
		stat_smooth(method = "lm", formula = y~poly(x, 2)) +
		labs(x = "Predicted HbA1c (mmol/mol)", y = "Observed HbA1c (mmol/mol)", title = "Injectable semaglutide") +
		theme_minimal() +
		theme(
				legend.position = "bottom"
		)



# Overall calibration summary ----

## Tolerance 3 mmol/mol ----

### Standard drugs & oral ----
overall_benefit_calibration_tolerance3_standard_oral_mice_match_sex_hba1c <- overall_predicted_benefit_performance(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide) %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var ="posthba1cfinal",
		pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
		conc_tolerance = 3,
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

## Rank 1 ----

### Standard drugs & oral ----
overall_benefit_calibration_rank1_standard_oral_mice_match_sex_hba1c <- overall_predicted_benefit_performance(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide) %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var ="posthba1cfinal",
		pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)


## Summary ----
calibration_summary <- do.call(rbind, overall_benefit_calibration_tolerance3_standard_oral_mice_match_sex_hba1c) %>%
		mutate(Matching = "Best drug, sex, hba1c", Dataset = "Standard & Oral", Method = "3mmol") %>%
		rbind(
			do.call(rbind, overall_benefit_calibration_rank1_standard_oral_mice_match_sex_hba1c) %>%
					mutate(Matching = "Best drug, sex, hba1c", Dataset = "Standard & Oral", Method = "rank1")
		)

saveRDS(calibration_summary, "/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.overall_benefit_calibration_summary.rds")


# Find optimal drug ----
## Tolerance 3 mmol/mol ----
optimal_drug_tolerance3_standard_oral <- get_best_drugs(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide),
		tolerance = 3,
		column_names = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
		final_var_name = "pred."
	) %>%
	select(contains("pred.within"))

# Overall calibration (plot benefit) ----

## Tolerance 3 mmol/mol ----
overall_benefit_tolerance3_standard_oral_mice_match_sex_hba1c <- overall_predicted_benefit(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide) %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var ="posthba1cfinal",
		cal_groups = 10,
		pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
		conc_tolerance = 3,
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)

## Rank 1 ----
overall_benefit_rank1_standard_oral_mice_match_sex_hba1c <- overall_predicted_benefit(
		data = analysis_standard %>% rbind(analysis_oral_semaglutide) %>% mutate(hba1c_group = ntile(prehba1c, 10)),
		drug_var = "drugclass",
		outcome_var ="posthba1cfinal",
		cal_groups = 10,
		pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "TZD", "SU")),
		matching_var = c("t2dmduration", "prebmi", "agetx", "prealt", "preegfr", "pretotalcholesterol", "prehdl", "hba1cmonth", "smoke", "imd5", "ncurrtx", "drugline"),
		match.exact = c("sex", "hba1c_group")
)


# Plots ----

## Best drug combinations ----
# Define groupings
groups <- list(
		"1-drug combination" = 1,
		"2-drug combinations" = 2,
		"3-drug combinations" = 3,
		"4/5-drug combinations" = 4:5
)

list_tolerance3_analysis_standard_oral <- optimal_drug_comparison_list(optimal_drug_tolerance3_standard_oral$pred.within_3_of_best_drug_name, groups)

saveRDS(list_tolerance3_analysis_standard_oral, "/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.list_tolerance3_analysis_standard_oral.rds")

## Tolerance 3 mmol/mol
### Standard drugs ----
plot_overall_benefit_tolerance_3_standard_oral_mice_match_sex_hba1c <- overall_benefit_tolerance3_standard_oral_mice_match_sex_hba1c %>%
		ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
		geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
		geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
		geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
		geom_point() +
		geom_errorbar() +
		geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
		theme_minimal() +
		labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Standard + Oral (tolerance 3mmol) (matching best drug / sex / hba1c_10)")

## Rank 1 ----
### Standard drugs ----
plot_overall_benefit_rank1_standard_oral_mice_match_sex_hba1c <- overall_benefit_rank1_standard_oral_mice_match_sex_hba1c %>%
  ggplot(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_point() +
  geom_errorbar() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) +
  theme_minimal() +
  labs(x = "Predicted HbA1c benefit (mmol/mol)", y = "Observed HbA1c benefit* (mmol/mol)", title = "Standard + Oral (matching best drug / sex / hba1c_10)")


## PDF ----
pdf("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions_a.pdf", width = 8, height = 5)
plot_pred_response_analysis_standard_preclosed
dev.off()

pdf("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions_b.pdf", width = 6, height = 4)
plot_pred_response_analysis_standard_dpp4
plot_pred_response_analysis_standard_tzd
plot_pred_response_analysis_semaglutide
dev.off()

pdf_combine(
		input = c("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions_a.pdf", "/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions_b.pdf"),
		output = "/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions.pdf"
		)
file.remove(c("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions_a.pdf", "/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.calibration_predictions_b.pdf"))

pdf("/home/pcardoso/workspace/Pedro-5drugmodel-framework/Outputs/SDRN/05.standard_oral_overall_calibration.pdf", width = 7, height = 5)
plot_overall_benefit_tolerance_3_standard_oral_mice_match_sex_hba1c
plot_overall_benefit_rank1_standard_oral_mice_match_sex_hba1c
dev.off()





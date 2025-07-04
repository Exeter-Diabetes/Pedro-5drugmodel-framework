# Author: pcardoso
###############################################################################

# Generate analysis cohort

# This is a function so you can call the same dataset over and over again
## without generating a new datasets other than original

###############################################################################

# function
set_up_data <- function(data = "mm_20250415_t2d_1stinstance") {
	
	# load libraries
	require(tidyverse)
	
	# load dataset
	load(paste0("/home/pcardoso/workspace/SDRN-Cohort-scripts/Final_Datasets/", data, ".RData")) # t2d_1stinstance
	
	
	###############################################################################
	# Flow diagram of inclusions and exclusions
	
	analysis_cohort <- t2d_1stinstance %>%
			## Exclusions based on type of drugsubstance: 
			### SU - non-gliclazide
			### TZD - rosiglitazone
			### GLP1-RA - lixisenatide, exenatide slow-release, semaglutide
			### MFN
			filter(
					drug_class != "MFN" & drug_class != "INS" &
							drug_class != "Glinide" & drug_class != "Acarbose"
			) %>%
			filter(
					drug_substance != "Glimepiride" & drug_substance != "Lixisenatide" &
							drug_substance != "Glipizide" & drug_substance != "Ertugliflozin" &
							drug_substance != "Glibenclamide" & drug_substance != "Tolbutamide" &
							drug_substance != "Tolazamide" & drug_substance != "Troglitazone" &
							drug_substance != "Gliquidone" & drug_substance != "Albiglutide" &
							drug_substance != "Chlorpropamide"
			) %>%
			## Exclusions
			### Currently treated with insulin
			### Initiating as first-line therapy
			### End-stage kidney disease
			### Age <80
			filter(is.na(INS) | INS != 1) %>%
			filter(drugline_all != 1) %>%
			filter(!(preckdstage %in% c("stage 5"))) %>%
			filter(dstartdate_age > 18 & dstartdate_age < 80) %>%
			## Exclusions
			### Multiple treatments on same day
			### Start within 61 days since last therapy
			### missing HbA1c
			### baseline HbA1c<53
			### baseline HbA1c>120
			filter(multi_drug_start_class == 0) %>%
			filter(timeprevcombo_class >= 61) %>%
			filter(!is.na(prehba1c)) %>%
			filter(prehba1c >= 53) %>%
			filter(prehba1c <= 110) %>%
			## Exclusions
			### missing baseline (BMI, eGFR, ALT, HDL, total cholesterol)
			mutate(t2dmduration = as.numeric(difftime(dstartdate, dm_diag_date_all, units = "days"))/365.25) %>%
			filter(!is.na(t2dmduration)) %>%
#			filter(!is.na(prebmi)) %>%
#			filter(!is.na(preegfr)) %>%
#			filter(!is.na(prealt)) %>%
#			filter(!is.na(prehdl)) %>%
#			filter(!is.na(pretotalcholesterol)) %>%
			# ## Exclusions
			# ### Missing outcome HbA1c
			mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
			# generating HbA1c month
			mutate(hba1cmonth_12 = as.numeric(difftime(posthba1c12mdate, dstartdate, units = "days")) / 30) %>%
			mutate(hba1cmonth_6 = as.numeric(difftime(posthba1c6mdate, dstartdate, units = "days")) / 30) %>%
			mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
			filter(!is.na(posthba1cfinal))
	
	# return final cohort	
	return(analysis_cohort)
	
}



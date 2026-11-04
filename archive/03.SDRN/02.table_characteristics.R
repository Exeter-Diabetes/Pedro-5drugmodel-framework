# Author: pcardoso
###############################################################################

# Table of characteristics

###############################################################################

# load libraries
library(tableone)

# load function for generating dataset
source("/home/pcardoso/workspace/Pedro-5drugmodel-framework/03.SDRN/01.cohort_script.R")

# load dataset
analysis_cohort_raw <- set_up_data(data = "mm_20250506_t2d_1stinstance")

###############################################################################

# Pre-processing datasets
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
		select(-row)

###############################################################################
# Table of characteristics for each dataset
analysis_cohort <- analysis_cohort %>%
		mutate(
				prebmi_missing = ifelse(is.na(prebmi), "Yes", "No"),
				preegfr_missing = ifelse(is.na(preegfr), "Yes", "No"),
				pretotalcholesterol_missing = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
				prehdl_missing = ifelse(is.na(prehdl), "Yes", "No"),
				prealt_missing = ifelse(is.na(prealt), "Yes", "No")
		)

vars <- c(
		"agetx", "sex", "t2dmduration", "ethnicity", 
		"drug_substance",
		"imd5", "smoke",
		"prebmi", "prebmi_missing", "prehba1c", "preegfr", "preegfr_missing", 
		"pretotalcholesterol", "pretotalcholesterol_missing", "prehdl", "prehdl_missing",
		"prealt", "prealt_missing",
		"drugline", "ncurrtx", "hba1cmonth",
		"posthba1cfinal"
)

cat_vars <- c(
		"sex", "ethnicity", "drug_substance", "imd5", "smoke",
		"prebmi_missing", "preegfr_missing", "pretotalcholesterol_missing", "prehdl_missing", "prealt_missing",
		"drugline", "ncurrtx"
)


#:-- Pre 2020-10-14
table_characteristics_cohort <- CreateTableOne(
		vars = vars,
		factorVars = cat_vars,
		includeNA = TRUE,
		strata = c("drug_class"),
		data = analysis_cohort,
		test = FALSE
)

table_characteristics_cohort_print <- print(table_characteristics_cohort, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)


write.csv(table_characteristics_cohort_print, file = "Outputs/SDRN/02.characteristics_cohort.csv")







# Initial set-up ###########################################

# load libraries
library(tidyverse)
library(aurum)
library(rms)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("mm")

# set up dataset
t2d_1stinstance <- t2d_1stinstance %>%
  analysis$cached("20250327_t2d_1stinstance")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

# Investigating Semaglutide initiations ----

library(broom)
library(tidyr)
library(purrr)
library(dplyr)
library(ggpubr)
library(patchwork)

semaglutide_initiations <- t2d_1stinstance %>%
  filter(
    drug_substance == "Low-dose semaglutide" | drug_substance == "Oral semaglutide" |
      drug_substance == "Semaglutide, dose unclear"
  ) %>%
  select(drug_substance, prehba1c, hba1cresp6m, dstartdate) %>%
  collect() %>%
  drop_na(hba1cresp6m) %>%
  mutate(yrdrugstart = format(dstartdate, "%Y")) %>%
  group_by(yrdrugstart, drug_substance) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(hba1cresp6m ~ prehba1c, data = .x)),
    tidied = map(model, ~tidy(.x, conf.int = TRUE)),
    augmented = map(model, augment),
    mean = map_dbl(augmented, ~mean(.x$.fitted, na.rm = TRUE)),
    sd_adjusted_response = map_dbl(augmented, ~sd(.x$.fitted, na.rm = TRUE)),
    n_initiations = map_int(data, nrow),
    n = map_int(data, nrow),
    se = sd_adjusted_response / sqrt(n),
    ci_low = mean - qt(0.975, df = n - 2) * se,
    ci_high = mean + qt(0.975, df = n - 2) * se
  ) %>%
  select(yrdrugstart, drug_substance, n_initiations, mean, ci_low, ci_high)

plot_test <- semaglutide_initiations %>%
  ggplot(aes(x = yrdrugstart, y = mean, ymin = ci_low, ymax = ci_high, colour = drug_substance)) +
  geom_point(aes(size = n_initiations), position=position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position=position_dodge(width=0.5)) +
  labs(y = "6-month HbA1c response", x = "Year of initiation", size = "Number of initiations", colour = "") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  ) +
  guides(
    colour = guide_legend(ncol = 2),
    size = guide_legend()
  )


# Step 1: Make sure yrdrugstart is numeric
test_ordered <- semaglutide_initiations %>%
  distinct(drug_substance, yrdrugstart, n_initiations) %>%
  mutate(yrdrugstart = as.numeric(as.character(yrdrugstart)))

# Step 2: Pivot wider to get drug_substance as rows and years as columns
table_data <- test_ordered %>%
  pivot_wider(
    names_from = yrdrugstart,
    values_from = n_initiations,
    values_fill = 0
  )

# Step 3: Order columns: drug_substance + sorted years
year_cols <- sort(unique(test_ordered$yrdrugstart))
table_data <- table_data %>%
  select("Drug" = drug_substance, all_of(as.character(year_cols)))  # reorder year columns explicitly

# Step 4: Create table
table_plot <- ggtexttable(table_data, rows = NULL, theme = ttheme("light"))

# Step 5: Combine with plot
final_plot <- plot_test / table_plot + plot_layout(heights = c(3, 1))

pdf("Outputs/CPRD/01.semaglutide_response.pdf", width = 9, height = 5)
final_plot
dev.off()


# Flow diagram: inclusion/exclusions ----
analysis = cprd$analysis("pedro_mm")
## Datasets
#:-- Post 2020-10-14
#:-- Pre 2020-10-14 (development cohort)
#:-- Semaglutide initiations
#:-- Oral semaglutide initiations

## Post 2020-10-14 ----
analysis_post_2020 <- t2d_1stinstance %>%
  ## Exclusions based on time - initiations after 14th October 2020
  filter(dstartdate > as.Date("2020-10-14")) %>%
  ## Exclusions based on type of drugsubstance: 
  ### SU - non-gliclazide
  ### TZD - rosiglitazone
  ### GLP1-RA - lixisenatide, exenatide slow-release, semaglutide
  ### MFN
  filter(
    drug_class != "MFN" & drug_class != "INS" & 
      drug_class != "Glinide"& drug_class != "Acarbose"
  ) %>%
  filter(
    drug_substance != "Glimepiride" & drug_substance != "Lixisenatide" &
      drug_substance != "Glipizide" & drug_substance != "Ertugliflozin" &
      drug_substance != "Glibenclamide" & drug_substance != "Tolbutamide" &
      drug_substance != "Low-dose semaglutide" & drug_substance != "Oral semaglutide" &
      drug_substance != "Semaglutide, dose unclear"
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
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  # filter(!is.na(prebmi)) %>%
  # filter(!is.na(preegfr)) %>%
  # filter(!is.na(prealt)) %>%
  # filter(!is.na(prehdl)) %>%
  # filter(!is.na(pretotalcholesterol)) %>%
  # ## Exclusions
  # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  filter(!is.na(posthba1cfinal)) %>%
  analysis$cached("analysis_post_2020", indexes=c("patid", "dstartdate", "drug_substance"))


## Pre 2020-10-14 ----
analysis_pre_2020 <- t2d_1stinstance %>%
  ## Exclusions based on time - initiations after 14th October 2020
  filter(dstartdate <= as.Date("2020-10-14")) %>%
  ## Exclusions based on type of drugsubstance:
  ### SU - non-gliclazide
  ### TZD - rosiglitazone
  ### GLP1-RA - lixisenatide, exenatide slow-release, semaglutide
  ### MFN
  filter(
    drug_class != "MFN" & drug_class != "INS" &
      drug_class != "Glinide"& drug_class != "Acarbose"
  ) %>%
  filter(
    drug_substance != "Glimepiride" & drug_substance != "Lixisenatide" &
      drug_substance != "Glipizide" & drug_substance != "Ertugliflozin" &
      drug_substance != "Glibenclamide" & drug_substance != "Tolbutamide" &
      drug_substance != "Low-dose semaglutide" & drug_substance != "Oral semaglutide" &
      drug_substance != "Semaglutide, dose unclear"
  ) %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <80
  filter(is.na(INS) | INS != 1) %>%
  filter(drugline_all != 1) %>%
  filter(!(preckdstage %in% c("stage 5"))) %>%
  mutate(dstartdate_age = (datediff(dstartdate, dob))/365.25) %>% # not created?
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
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  filter(!is.na(prebmi)) %>%
  filter(!is.na(preegfr)) %>%
  filter(!is.na(prealt)) %>%
  filter(!is.na(prehdl)) %>%
  filter(!is.na(pretotalcholesterol)) %>%
  # ## Exclusions
  # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  filter(!is.na(posthba1cfinal)) %>%
  analysis$cached("analysis_pre_2020", indexes=c("patid", "dstartdate", "drug_substance"))


## Semaglutide initiations ----
analysis_semaglutide <- t2d_1stinstance %>%
  ## Inclusion based on type of drugsubstance:
  ### GLP1-RA - semaglutide
  filter(
    drug_substance == "Low-dose semaglutide"
  ) %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <80
  filter(is.na(INS) | INS != 1) %>%
  filter(drugline_all != 1) %>%
  filter(!(preckdstage %in% c("stage 5"))) %>%
  mutate(dstartdate_age = (datediff(dstartdate, dob))/365.25) %>% # not created?
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
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  # filter(!is.na(prebmi)) %>%
  # filter(!is.na(preegfr)) %>%
  # filter(!is.na(prealt)) %>%
  # filter(!is.na(prehdl)) %>%
  # filter(!is.na(pretotalcholesterol)) %>%
  # ## Exclusions
  # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  filter(!is.na(posthba1cfinal)) %>%
  analysis$cached("analysis_semaglutide", indexes=c("patid", "dstartdate", "drug_substance"))



## Oral semaglutide initiations ----
analysis_oral_semaglutide <- t2d_1stinstance %>%
  ## Inclusion based on type of drugsubstance:
  ### GLP1-RA - semaglutide
  filter(
    drug_substance == "Oral semaglutide"
  ) %>%
  ## Exclusions
  ### Currently treated with insulin
  ### Initiating as first-line therapy
  ### End-stage kidney disease
  ### Age <80
  filter(is.na(INS) | INS != 1) %>%
  filter(drugline_all != 1) %>%
  filter(!(preckdstage %in% c("stage 5"))) %>%
  mutate(dstartdate_age = (datediff(dstartdate, dob))/365.25) %>% # not created?
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
  mutate(t2dmduration = datediff(dstartdate, dm_diag_date_all)/365.25) %>%
  filter(!is.na(t2dmduration)) %>%
  # filter(!is.na(prebmi)) %>%
  # filter(!is.na(preegfr)) %>%
  # filter(!is.na(prealt)) %>%
  # filter(!is.na(prehdl)) %>%
  # filter(!is.na(pretotalcholesterol)) %>%
  # ## Exclusions
  # ### Missing outcome HbA1c
  mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%  # generating posthba1c
  # generating HbA1c month
  mutate(hba1cmonth_12 = datediff(posthba1c12mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth_6 = datediff(posthba1c6mdate, dstartdate) / 30) %>%
  mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
  filter(!is.na(posthba1cfinal)) %>%
  analysis$cached("analysis_oral_semaglutide", indexes=c("patid", "dstartdate", "drug_substance"))



# Framework: Testing the 5-drug model

This repository contains the code used in the study presenting a framework for validating the 5-drug model.

The files are divided into three folders:

-   01.Functions: Functions used during in the framework
-   02.CPRD: R code used to validate the CPRD cohort
-   03.SDRN: R code used to validate the SDRN cohort

## Functions

<br>

Functions used in this framework:

- `closedtest_continuous_function.R`: Closed-test procedure to assess whether predictions need to be recalibration (intercept or intercept + slope) is needed.
- `get_ranked_or_tolerant_drugs.R`: This function identifies the drug(s) with the best (lowest) predicted outcome from a set of model predictions. It can return the nth-best drug or all drugs within a specified tolerance of the best drug.
- `get_best_drugs.R`: Applies `get_ranked_or_tolerant_drugs` across all rows of a dataset to determine the best (or nearly best) drugs based on model predictions. Adds the result as new columns to the data.
- `heterogenous_effect_calibration.R`: Estimates heterogeneous treatment effects between two drugs by stratifying patients into calibration groups based on predicted benefit scores. Optionally performs covariate matching before estimating treatment effects.
- `unified_validation.R`: This function computes heterogeneous treatment effect calibration curves for all pairwise comparisons between a set of specified drugs. For each drug pair, it calculates the predicted benefit (the difference between their predicted outcomes), partitions the data into calibration groups based on predicted benefit, and estimates treatment effects within each group. The function internally calls `heterogenous_effect_calibration` for each drug pair.
- `overall_predicted_benefit.R`: This function assesses how well predicted treatment benefits correspond to observed differences in outcomes. Patients are grouped by predicted benefit, and differences in observed outcomes between concordant and discordant treatments are regressed against the predicted differences. Patients are labeled as concordant if they received a treatment predicted to be optimal (either the top predicted or within a specified tolerance of the best). Matching is used to control for confounding between treatment groups.


## CPRD

<br>

The `CPRD` folder contains the scripts used to generate and validate the CPRD cohorts. The CPRD cohorts were defined according to the GitHub repository: <https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts>

Two cohorts have been developed:

- Treatment initiations pre-2020-10-14
- Treatment initiations post-2020-10-14


## SDRN

<br>

The `SDRN` folder contains the scripts used to generate and validate the SDRN cohort. The SDRN cohort was defined according to the GitHub repository: <https://github.com/Exeter-Diabetes/SDRN-Cohort-scripts>


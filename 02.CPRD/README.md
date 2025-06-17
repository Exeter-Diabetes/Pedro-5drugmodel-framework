# Framework: Testing the 5-drug model

## CPRD cohort validation

The `CPRD` folder contains the scripts used to generate and validate the CPRD cohorts. The CPRD cohorts were defined according to the GitHub repository: <https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts>

Two cohorts have been developed:

- Treatment initiations pre-2020-10-14
- Treatment initiations post-2020-10-14


Files structure:

- `01.cohort_scripts.R`: defining cohorts based on date and initial data cleaning.
- `02.table_characteristics.R`: table of characteristics for the cohorts.
- `03.impute_missingness.R`: imputation methods for missing data with function to help deployment.
- `04.model_predictions.R`: prediction from 5-drug model with function to help deployment.
- `05.per_drug_calibration.R`: calibration of overall benefit with examples of deployment.
- `06.per_drug_pair_calibration.R`: calibration of drug pairs benefits with examples of deployment.
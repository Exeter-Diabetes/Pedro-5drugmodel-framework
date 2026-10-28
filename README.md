# Pedro-5drugmodel-updates
This repository contains all the developments from the original T2D treatment selection model.

Updates of the Exeter 5-drug model in CPRD Aurum: a temporal update, validation by
ethnicity, and new drug classes (oral and injectable semaglutide, tirzepatide).

- `example_validation_script_v1.R`: validate the published model in your own data
- `example_validation_script_v2.R`: validate the updated model in your own data
- `Scripts/`: the CPRD analyses for the paper

Run the scripts from the repository root. The model files go in `Model_files/`
(not included in the repository). Calibration uses the
[cateval](https://github.com/PM-Cardoso/cateval) package.

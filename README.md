# Validating and updating the Exeter 5-drug treatment-selection model

This repository has the framework to validate the Exeter 5-drug HbA1c
treatment-selection model (Dennis et al. Lancet 2025) in a new cohort, and the CPRD
Aurum analyses we used to update it: a temporal update, a check by ethnicity, and new
drug classes (oral and injectable semaglutide, tirzepatide).

The validation is of the **predicted treatment benefit**: do patients who received the
drug class the model recommends do better than similar patients who received another
one, and by as much as the model predicted? The predictions are checked as they are,
and after recalibration for each drug class.

## Structure

|                                                                     |                                                                                                                                                                                     |
| ------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`example_validation_script_v1.R`](example_validation_script_v1.R) | Example: validate the **published model** (SGLT2i, GLP-1RA, DPP4i, SU, TZD) in your own data.                                                                                 |
| [`example_validation_script_v2.R`](example_validation_script_v2.R) | Example: validate the **updated model** in your own data, with GLP-1RA split into dulaglutide and liraglutide, plus oral semaglutide, injectable semaglutide and tirzepatide. |
| [`Scripts/`](Scripts)                                              | The CPRD Aurum analyses for the model update.                                                                                                                                       |
| [`archive/`](archive)                                              | SDRN validation scripts from the first version of this repository, kept for future work. Not maintained.                                                                            |

## Example scripts

The two example scripts are the ones to adapt to your own data. Each runs from data
preparation to the saved figures and tables; the steps and required variables are
described at the top of the script and in Step 1a.

They use our R package [cateval](https://github.com/PM-Cardoso/cateval):

```r
# install.packages("remotes")
remotes::install_github("PM-Cardoso/cateval")
```

Model predictions come from the online calculator or from the model files
(`use_calculator`, Step 0). The model files are not included in this repository.

## License

See [`LICENSE`](LICENSE).

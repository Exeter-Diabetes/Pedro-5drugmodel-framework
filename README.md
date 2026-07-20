# A Framework for Validating the Exeter 5-Drug Treatment-Selection Model

This repository contains the **framework** for validating (and, where needed, recalibrating)
the Exeter 5-drug HbA1c treatment-selection model in a new cohort. The model predicts 12-month
HbA1c on each of five drug classes — **SGLT2, GLP1, DPP4, SU, TZD** — so that the class with the
lowest predicted HbA1c can be recommended.

The quantity being validated is the **predicted treatment benefit** (the difference in predicted
HbA1c between a patient's model-recommended drug and the alternatives), not the absolute HbA1c.
The framework asks: *do patients who received their model-recommended drug achieve a better observed
HbA1c than clinically similar patients who received a different drug, and is that advantage the size
the model predicted?*

The complete, runnable walk-through is [`example_validation_script.R`](example_validation_script.R);
the analysis rationale is in [`VALIDATION_PLAN.md`](VALIDATION_PLAN.md).

---

## How the framework works

Validation is organised as a **staged procedure**. Each stage re-runs the same calibration
analyses but changes what is being tested, so you only move to the next (more interventionist) stage
if the model is not yet adequate. This avoids recalibrating a model whose apparent miscalibration is
really a case-mix / confounding artefact.

| Stage | Question | Predictions used | What is added |
|-------|----------|------------------|---------------|
| **A — New population testing** | Does the original model work as is? | Original | Model covariates only |
| **B — Differences in treatment allocation** | Is the miscalibration explained by prognostic factors that differ between the drug groups? | Original | Prognostic-factor adjustment (pairwise) / matching (overall) |
| **C — Adjust baseline risk** | Is the baseline (intercept) effect different in this population? | **Intercept-recalibrated** | Per-drug intercept shift |
| **D — Adjust baseline + predictor effect** | Is the predictor (slope) effect different in this population? | **Intercept + slope recalibrated** | Per-drug intercept and slope |

**Decision logic.** If calibration is adequate at Stage A, stop — the model transfers as is. If not,
Stage B tests whether accounting for prognostic factors resolves it. If recalibration is required,
Stages C and D apply per-drug intercept (and then slope) recalibration. Recalibration is done
**per drug class**, not as one overall shift, because the five classes differ pharmacologically and
in prescribing/adherence, so population differences are expected to affect each class differently.

At every stage the model is judged adequate when the benefit-calibration curve lies close to the
line of identity within its confidence intervals, **calibration-in-the-large ≈ 0** and
**calibration slope ≈ 1**, and the mean observed benefit of following the model is positive and
consistent with what was predicted.

### The concordant vs discordant design

- **Best drug** = the class with the lowest predicted 12-month HbA1c.
- **Concordant** = received the model's best drug; **discordant** = received a different drug.

Concordant and discordant patients are made comparable with **1:1 nearest-neighbour matching**
(Mahalanobis distance, exact match on the predicted best drug, 0.4 SD caliper on baseline HbA1c,
matching with replacement). Matching is run **symmetrically** (concordant-as-treated and
discordant-as-treated) and the two matched sets are combined. For each matched pair the framework
computes the **predicted benefit** and the **observed benefit**, then summarises calibration.

---

## Outputs produced

For each stage the script writes a consistent set of figures (as PDFs to your `output_path`):

- **Pairwise drug calibration** — observed vs predicted benefit for every drug pair.
- **Overall benefit calibration** — concordant / discordant / combined matched sets, side by side.
- **Overall benefit by drug class** — the same, faceted by the drug that anchors each pair.
- **Predicted vs observed calibration** — per-drug predicted-vs-observed HbA1c.
- **Balance (love) plots** — matching balance for each matching direction (Stages A and B).

Supporting tables include the per-drug recalibration parameters and a baseline-characteristics
table. An **additional clinical-decision analysis** reports how often the model-recommended (best)
drug changes between the original and recalibrated predictions.

---

## The `cateval` package

The analytical engine is the **`cateval`** R package (developed for this framework). The example
script calls it directly; you do not need the individual function files to run the framework.

```r
# install.packages("remotes")
remotes::install_github("PM-Cardoso/cateval")
```

Key functions used by the framework:

- `match_benefit_pairs()` — labels each patient concordant/discordant, runs the symmetric 1:1
  matching, and returns the `concordant`, `discordant` and combined matched populations plus the
  underlying `matchit` objects (for love plots).
- `benefit_by_group()` — groups matched pairs by predicted benefit and returns mean observed vs
  predicted benefit per group, with analytical or bootstrap confidence intervals.
- `benefit_by_drug()` — the same, split by the index patient's drug class; accepts several
  `cal_groups` settings (e.g. `c(3, 5, 10)`) and stacks them with an `n_groups` column.
- `benefit_overall()` — the mean observed benefit of following the model.
- `benefit_citl_slope()` — calibration-in-the-large and calibration slope.
- `pairwise_calibration()` — pairwise benefit-calibration across all drug pairs.
- `get_best_drugs()` / `get_ranked_or_tolerant_drugs()` — identify each patient's best (or
  near-best) drug, used for the clinical-decision analysis.

---

## Getting started

1. **Install** R (≥ 4.3) and the packages used by the script:
   `tidyverse`, `tableone`, `cobalt`, `mice` (for imputation) and `cateval`.
2. **Prepare your dataset** as a single data frame with the required variables (see the top of
   [`example_validation_script.R`](example_validation_script.R) for the full specification): the
   12-month outcome `posthba1cfinal`, `hba1cmonth`, `drugclass`, and the baseline covariates (age,
   sex, diabetes duration, baseline HbA1c, BMI, eGFR, ALT, total/HDL cholesterol, smoking, drug
   line, number of current therapies; ethnicity where a stratified analysis is run). Categorical
   variables must be factors with the specified levels; impute missing baseline covariates with
   `mice`.
3. **Obtain model predictions** for each drug class (`pred.SGLT2` … `pred.TZD`), e.g. from the
   5-drug model calculator, and join them back into your dataset.
4. **Run the script**, adapting the data-preparation and `output_path` sections to your source.
   It generates the original, intercept-recalibrated and intercept + slope-recalibrated predictions
   and writes every stage's figures.

The example script is the reference implementation; adapt the data-loading and formatting to your
cohort and the rest of the pipeline is intended to run unchanged.

---

## License

See [`LICENSE`](LICENSE).

## Check of the intercept update of the GLP-1RA drug classes (prototype, to be
##  moved into example_validation_script_v2.R once the figures are settled).
##
## In v2 every GLP-1RA drug class is the published model's GLP1 prediction plus its
##  own intercept: the GLP1 arm's differential effects (its interactions with the
##  patient's characteristics) are assumed to hold for the new drug classes, shifted
##  by a constant. This script checks that assumption.
##
## Each initiation of a drug class (treated) is matched 1:1, with replacement, to an
##  initiation of a reference drug class, on a propensity score with a caliper on the
##  variable of the x axis, so both patients of a pair have about the same x. Each
##  matching is done once (with the original predictions) and each version of the
##  predictions (original, recalibration-in-the-large, recalibration) is a row of
##  the figure, with a column per treated drug class, all on the same axes. Points
##  are the means of groups of pairs ordered by x (95% CI), the blue curve a natural
##  spline over all pairs (95% CI), and each panel gives the slope per 10 mmol/mol of
##  x of the observed minus predicted difference (0 if the model's differential
##  effects hold).
##
## Two kinds of figure (why, benefits and problems of each: 02.intercept_check_plots.md):
##   1. The method (method_reference, method_treated): the injectable drug classes
##      matched to oral semaglutide, whose intercept (+0.24) leaves it at the
##      published GLP1 prediction. Within the GLP1 family the model predicts the same
##      difference for every patient (the intercepts differ), so the y axis is the
##      observed difference itself (treated minus oral semaglutide) and the red line
##      the predicted difference: flat for the original and recalibration-in-the-large
##      predictions, sloped with recalibration (a slope per drug class). X axis:
##      predicted HbA1c.
##   2. The check against drug classes outside the GLP1 family (check_references:
##      SGLT2 and DPP4, with their own modelled effects): does each GLP-1RA drug class
##      follow the GLP1 arm's differential effects? The predicted difference differs
##      between patients, so the y axis is the observed minus predicted difference
##      (the treated patient's prediction error minus the reference patient's), with
##      the red line at 0. X axis: predicted HbA1c, and baseline HbA1c (check_x_axes).
##
## Outputs (Scripts/04.V2_Validation/Outputs/Intercept_check):
##   - intercept_check_method_vs_Oral_Semaglutide.pdf (figure 1)
##   - intercept_check_<x axis>_vs_<reference>.pdf (figure 2)
##   - Loveplots/: the balance of each matching
##   - intercept_check_estimates.csv and intercept_check_groups.csv: every figure
##
## Run Scripts/04.V2_Validation/01.validation_v2_post2020.R first, in the same R
##  session, at least up to its Part 5.1: everything before that loads the data,
##  makes the predictions and defines the settings and figure code that this
##  script uses.
## Outputs go to Scripts/04.V2_Validation/Outputs/Intercept_check.


#############################################################################

needed <- c("analysis_dataset", "drugs", "injectable_drugs", "outcome_var", "drug_labels",
            "matching_vars_prog", "love_labels", "use_bootstrap", "n_boot", "boot_seed",
            "save_plot", "fit_font_size", "x_overhang", "axis_title_size", "wrap_to_width",
            "n_axis_levels")
not_found <- needed[!vapply(needed, exists, logical(1))]
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("The sandwich package is needed for the confidence intervals: install.packages(\"sandwich\")",
       call. = FALSE)
}
if (length(not_found) > 0) {
  stop("Run Scripts/04.V2_Validation/01.validation_v2_post2020.R up to Part 5.1 in this ",
       "R session first (not found: ", paste(not_found, collapse = ", "), ").", call. = FALSE)
}


# Part 1: Settings ------------------------------------------------------------

# the versions of the predictions, one row each (prefix of their columns). The
# original predictions are the updated model (the published model plus the
# intercepts of the model files). The matching uses the original predictions.
check_versions <- c("Original predictions"        = "pred.",
                    "Recalibration-in-the-large"  = "pred.recal_large.",
                    "Recalibration"               = "pred.recal.")

# the GLP-1RA drug classes added with an intercept (drug classes without
# initiations are left out)
glp1_checked <- c("Dulaglutide", "Liraglutide", "Inje_Semaglutide", "TZP")

# Figure 1, the method: the injectable drug classes (treated) matched to oral
# semaglutide (the control: at the published GLP1 prediction), on the observed
# difference scale, by predicted HbA1c. (Another GLP-1RA drug class can be the
# control, e.g. "Dulaglutide" with method_treated "Inje_Semaglutide".)
method_reference <- "Oral_Semaglutide"
method_treated   <- glp1_checked

# Figure 2: the reference drug classes outside the GLP1 family, each with the drug
# classes matched to it (treated), on the observed minus predicted scale
check_references <- list(
  SGLT2 = c(glp1_checked, "Oral_Semaglutide"),
  DPP4  = c(glp1_checked, "Oral_Semaglutide")
)

# the x axes of figure 2 (figure 1 uses the predicted HbA1c), each the variable
# plotted and matched on with a caliper (x_caliper):
#   predicted: the treated drug class's predicted HbA1c for the treated patient
#   baseline:  baseline HbA1c
check_x_axes <- c(predicted = "predicted HbA1c", baseline = "baseline HbA1c")

# Caliper (SD) on the x variable in the matching. Without it, a patient at either end
# of the x range is matched to someone nearer the middle, which tilts the plot (in
# simulated data with no trend, a slope of about +1 mmol/mol per 10 mmol/mol).
x_caliper <- 0.25

# How the pairs are matched (the settings to change if a love plot shows poor
# balance; matching_vars_prog, from the validation script, are the covariates):
#   match_distance: "glm" matches on the propensity score (the probability of the
#     treated drug class rather than the reference given the covariates), which
#     balances all covariates together and loses fewer matches than "mahalanobis".
#   ps_caliper: caliper (SD of the propensity score on the logit scale) for "glm"; 0.2
#     is the usual choice. NULL for none.
#   covariate_calipers: calipers (SD) on single covariates, e.g. c(prebmi = 0.5).
#   match_exact: covariates matched exactly, e.g. "drugline". Each costs matches.
#   extra_matching_vars: terms added to the propensity score model, e.g.
#     c("I(prehba1c^2)", "I(prebmi^2)").
match_distance      <- "glm"
ps_caliper          <- 0.2
covariate_calipers  <- c()
match_exact         <- c()
extra_matching_vars <- c()

# The pairs are ordered by x and split into groups of equal size, and each group's
# mean is plotted. The figure shows the largest number of groups in check_groups
# with at least min_group_pairs pairs in every group (or the smallest number).
check_groups    <- c(5, 10)
min_group_pairs <- 100

check_output_path <- "Scripts/04.V2_Validation/Outputs/Intercept_check"
loveplot_path     <- paste0(check_output_path, "/Loveplots")
dir.create(loveplot_path, showWarnings = FALSE, recursive = TRUE)


#############################################################################

# Part 2: Matching -------------------------------------------------------------

# helper: each initiation of treated matched 1:1, with replacement, to the most
# similar initiation of reference (propensity score on matching_vars_prog by
# default, Part 1), with a caliper on the x variable. Returns the matchit object
# and its data (for the love plot) and the matched pairs:
#   - treated_row, reference_row: the rows of analysis_dataset of the two
#     initiations (a reference initiation can be in several pairs)
#   - obs_diff: outcome HbA1c of the treated initiation minus that of the reference
match_pair <- function(treated, reference, x_axis) {
  match_data <- analysis_dataset %>%
    mutate(row_id = row_number()) %>%
    filter(drugclass %in% c(treated, reference)) %>%
    mutate(
      check_treated = as.integer(drugclass == treated),
      # the x variable (the original predictions for the predicted HbA1c)
      x_match = if (x_axis == "predicted") .data[[paste0(check_versions[[1]], treated)]] else prehba1c
    ) %>%
    as.data.frame()

  # covariates with a single value cannot be matched on
  match_vars <- matching_vars_prog[vapply(match_data[matching_vars_prog],
                                          function(x) length(unique(x)) > 1, logical(1))]

  # calipers: the x variable, any on single covariates and, for the propensity
  # score, an unnamed element
  calipers <- c(x_match = x_caliper, covariate_calipers)
  if (match_distance == "glm" && !is.null(ps_caliper)) calipers <- c(calipers, ps_caliper)

  set.seed(19840503)
  match_obj <- MatchIt::matchit(
    reformulate(c(match_vars, extra_matching_vars), response = "check_treated"),
    data = match_data,
    method = "nearest", distance = match_distance,
    link = if (match_distance == "glm") "linear.logit" else "logit",
    replace = TRUE, ratio = 1,
    exact = if (length(match_exact) > 0) reformulate(match_exact),
    caliper = calipers, std.caliper = TRUE
  )

  matched <- MatchIt::get_matches(match_obj, data = match_data)
  pairs <- matched %>%
    filter(check_treated == 1) %>%
    transmute(subclass, treated_row = row_id, outcome_treated = .data[[outcome_var]]) %>%
    inner_join(
      matched %>%
        filter(check_treated == 0) %>%
        transmute(subclass, reference_row = row_id, outcome_reference = .data[[outcome_var]]),
      by = "subclass"
    ) %>%
    mutate(drug = treated, reference = reference, x_axis = x_axis,
           obs_diff = outcome_treated - outcome_reference)

  cat(treated, "vs", reference, "(x:", x_axis, "):", nrow(pairs), "of", sum(match_data$check_treated),
      "initiations matched to", n_distinct(pairs$reference_row), "reference initiations\n")

  list(matchit = match_obj, data = match_data, pairs = pairs)
}

# helper: the pairs with the predictions of one version (prefix):
#   - x: the x variable of the treated patient (predicted HbA1c of the treated drug
#     class in this version, or baseline HbA1c)
#   - pred_diff: the predicted difference of the pair: each patient's prediction for
#     the drug class they received (treated patient on the treated drug class minus
#     reference patient on the reference). Using each patient's own covariates, any
#     imbalance left after matching is accounted for by the model's prediction (with
#     both predictions from the treated patient, it tilts the plot: in simulated
#     data with no trend, slopes of up to 4 mmol/mol per 10).
#   - gap: observed minus predicted difference, i.e. the treated patient's prediction
#     error minus the reference patient's (0 on average if the model is right)
#   - pred_diff_treated: the predicted difference for the treated patient (both
#     drug classes' predictions for that patient), the model's predicted treatment
#     effect, and gap_treated, observed minus this. Used in figure 1, where (within
#     the GLP1 family) it is the same for every patient but with recalibration.
pairs_for_version <- function(pairs, prefix, version) {
  pred_treated   <- analysis_dataset[[paste0(prefix, pairs$drug[1])]][pairs$treated_row]
  pred_reference <- analysis_dataset[[paste0(prefix, pairs$reference[1])]][pairs$reference_row]
  pred_reference_treated <- analysis_dataset[[paste0(prefix, pairs$reference[1])]][pairs$treated_row]
  pairs %>%
    mutate(version = version,
           x = if (pairs$x_axis[1] == "predicted") pred_treated else analysis_dataset$prehba1c[pairs$treated_row],
           pred_diff = pred_treated - pred_reference,
           gap = obs_diff - pred_diff,
           pred_diff_treated = pred_treated - pred_reference_treated,
           gap_treated = obs_diff - pred_diff_treated)
}


#############################################################################

# Part 3: Observed minus predicted difference against x --------------------------

# helper: for the matched pairs of one drug class, reference, x axis and version,
# with y the variable plotted ("gap": observed minus predicted difference, figure 2;
# "obs_diff": observed difference, figure 1) and the slope of slope_var (observed
# minus predicted: "gap", or "gap_treated" in figure 1),
#   - groups: the pairs ordered by x in groups of equal size: mean x and mean y
#     (95% CI)
#   - curve: a natural spline (3 df) of y on x, drawn over the middle 98% of x (few
#     pairs at the ends)
#   - estimate: the mean predicted and observed difference; the slope of the
#     observed minus predicted difference per 10 mmol/mol of x (0 if the model's
#     differential effects hold)
# 95% CIs account for the reference initiations being in several pairs (matching
# with replacement; a smaller reference drug class, e.g. oral semaglutide, has each
# initiation in several pairs). They are cluster-robust (sandwich, clustered on the
# reference initiation) or, when use_bootstrap is TRUE, from a bootstrap of the
# reference initiations (each with all its pairs). Treating the pairs as independent
# gives CIs about half as wide as they should be: in simulated data with no trend,
# slopes against oral semaglutide ranged from -2.0 to +1.5 across datasets, with a
# naive standard error of 0.5.
fit_gap <- function(pairs, y = "gap", slope_var = "gap") {
  pairs$y     <- pairs[[y]]
  pairs$slope <- pairs[[slope_var]]
  grid <- data.frame(x = seq(quantile(pairs$x, 0.01), quantile(pairs$x, 0.99), length.out = 100))
  spline_fit <- function(d) lm(y ~ splines::ns(x, df = 3), data = d)
  slope_fit  <- function(d) lm(slope ~ I(x / 10), data = d)

  # groups of x: the largest number in check_groups with at least min_group_pairs
  # pairs in each, or the smallest number
  ok_groups <- check_groups[floor(nrow(pairs) / check_groups) >= min_group_pairs]
  n_groups  <- if (length(ok_groups) > 0) max(ok_groups) else min(check_groups)
  pairs$group <- dplyr::ntile(pairs$x, n_groups)
  group_means <- function(d) as.numeric(tapply(d$y, factor(d$group, levels = seq_len(n_groups)), mean))

  curve_model <- spline_fit(pairs)
  slope_model <- slope_fit(pairs)
  gap_model   <- lm(slope ~ 1, data = pairs)

  # cluster-robust variance (clustered on the reference initiation)
  cluster_vcov <- function(model) sandwich::vcovCL(model, cluster = ~ reference_row)
  z <- qnorm(0.975)

  curve_matrix <- model.matrix(delete.response(terms(curve_model)), data = grid)
  curve_fit <- as.numeric(curve_matrix %*% coef(curve_model))
  curve_se  <- sqrt(rowSums((curve_matrix %*% cluster_vcov(curve_model)) * curve_matrix))
  curve <- grid %>% mutate(fit = curve_fit, lower = curve_fit - z * curve_se, upper = curve_fit + z * curve_se)

  slope    <- coef(slope_model)[[2]]
  slope_se <- sqrt(cluster_vcov(slope_model)[2, 2])
  slope_ci <- slope + c(-1, 1) * z * slope_se
  slope_p  <- 2 * pnorm(-abs(slope / slope_se))
  gap_se   <- sqrt(cluster_vcov(gap_model)[1, 1])
  gap_ci   <- mean(pairs$slope) + c(-1, 1) * z * gap_se

  # each group's mean, with a cluster-robust CI
  groups <- pairs %>%
    group_by(group) %>%
    group_modify(function(d, key) {
      m <- lm(y ~ 1, data = d)
      se <- if (n_distinct(d$reference_row) > 1) sqrt(cluster_vcov(m)[1, 1]) else NA_real_
      tibble(n = nrow(d), x = mean(d$x), y = mean(d$y), lower = y - z * se, upper = y + z * se)
    }) %>%
    ungroup()

  if (isTRUE(use_bootstrap)) {
    set.seed(boot_seed)
    # bootstrap of the reference initiations, each with all its pairs
    reference_pairs <- split(seq_len(nrow(pairs)), pairs$reference_row)
    boot <- replicate(n_boot, {
      d <- pairs[unlist(sample(reference_pairs, length(reference_pairs), replace = TRUE), use.names = FALSE), ]
      c(coef(slope_fit(d))[[2]], mean(d$slope), predict(spline_fit(d), newdata = grid), group_means(d))
    })
    slope_ci <- quantile(boot[1, ], c(0.025, 0.975))
    slope_p  <- min(1, 2 * min(mean(boot[1, ] <= 0), mean(boot[1, ] >= 0)))
    gap_ci   <- quantile(boot[2, ], c(0.025, 0.975))
    curve$lower <- apply(boot[2 + seq_len(nrow(grid)), ], 1, quantile, 0.025)
    curve$upper <- apply(boot[2 + seq_len(nrow(grid)), ], 1, quantile, 0.975)
    group_boot  <- boot[-(1:(2 + nrow(grid))), , drop = FALSE]
    groups$lower <- apply(group_boot, 1, quantile, 0.025, na.rm = TRUE)
    groups$upper <- apply(group_boot, 1, quantile, 0.975, na.rm = TRUE)
  }

  ids <- pairs[1, c("drug", "reference", "x_axis", "version")]
  list(
    curve  = bind_cols(curve, ids),
    groups = bind_cols(groups, ids) %>% mutate(n_groups = n_groups),
    estimate = bind_cols(ids, tibble(
      n_pairs = nrow(pairs),
      n_reference = n_distinct(pairs$reference_row),
      predicted_difference = mean(pairs$obs_diff - pairs$slope),
      observed_difference = mean(pairs$obs_diff),
      gap = mean(pairs$slope), gap_lower_ci = gap_ci[[1]], gap_upper_ci = gap_ci[[2]],
      slope_per_10 = slope, slope_lower_ci = slope_ci[[1]], slope_upper_ci = slope_ci[[2]],
      slope_p = slope_p
    ))
  )
}


#############################################################################

# Part 4: Figure code ----------------------------------------------------------
# (to go into the figure code block of the examples once settled)

## Intercept check figure ----

# One panel (3.6 x 3.0 inches) per treated drug class (columns) and version of the
# predictions (rows), all on the same axes: the mean of each group of pairs (points,
# 95% CI), the model's prediction (red line) and the spline over all pairs (blue
# curve, 95% CI). Column titles: the drug class and its number of matched pairs. In
# each panel, the slope per 10 mmol/mol of x (95% CI) of the observed minus
# predicted difference: 0 if the model's differential effects hold. scale:
#   "gap": y is the observed minus predicted difference; the red line is 0
#   "observed": y is the observed difference; the red line is the predicted
#     difference for the treated patient (pred_diff_treated), with a dashed line at 0
draw_intercept_check <- function(pairs, groups, curves, estimates, reference, x_axis, scale = "gap") {
  panel_width  <- 3.6
  panel_height <- 3.0
  panel_drugs <- intersect(names(drug_labels), unique(pairs$drug))
  versions   <- names(check_versions)
  n_cols     <- length(panel_drugs)
  n_rows     <- length(versions)
  text_scale <- 1 + 0.2 * (n_cols - 1)
  axis_text_size <- 15 * text_scale

  # column titles (the drug class and its number of pairs) and row titles, as large
  # as the longest allows within a panel, but no larger than the axis titles
  n_pairs <- estimates$n_pairs[match(panel_drugs, estimates$drug)]
  col_labels <- setNames(paste0(drug_labels[panel_drugs], " (n = ",
                                format(n_pairs, big.mark = ",", trim = TRUE), ")"),
                         panel_drugs)
  col_size <- min(16 * text_scale, fit_font_size(col_labels, panel_width - 0.25))
  row_size <- min(16 * text_scale, fit_font_size(versions, panel_height - 0.3))

  # the slope, written in the top left of each panel
  fmt <- function(x) formatC(x, format = "f", digits = 1)
  estimates$label <- paste0("Slope ", fmt(estimates$slope_per_10), " (", fmt(estimates$slope_lower_ci),
                            " to ", fmt(estimates$slope_upper_ci), ")")
  label_size <- min(13 * text_scale, fit_font_size(estimates$label, panel_width - 0.3))

  as_panels <- function(d) d %>% mutate(drug = factor(drug, levels = panel_drugs),
                                        version = factor(version, levels = versions))
  pairs <- as_panels(pairs); groups <- as_panels(groups)
  curves <- as_panels(curves); estimates <- as_panels(estimates)

  # the red line: 0, or the predicted difference along x
  observed  <- scale == "observed"
  red_lines <- if (observed) pairs %>% arrange(x) %>% transmute(drug, version, x, y = pred_diff_treated)

  # the same axes in every panel: x the pairs' x, y every group mean, interval, line
  # and band, with room at the top for the slope
  x_limits <- range(pairs$x)
  y_limits <- range(c(groups$lower, groups$upper, 0, curves$lower, curves$upper, red_lines$y))
  y_limits[2] <- y_limits[2] + 0.18 * diff(y_limits)

  # room between panels, and on the right, for the outermost x axis numbers
  overhang  <- x_overhang(x_limits, axis_text_size, panel_width)
  panel_gap <- sum(overhang) + 0.12
  title_width <- n_cols * panel_width + 0.8

  x_title <- if (x_axis == "predicted") {
    expression("Predicted HbA"[1*c]*" (mmol/mol)")
  } else {
    expression("Baseline HbA"[1*c]*" (mmol/mol)")
  }
  y_title <- if (observed) "Observed HbA1c difference (mmol/mol)" else "Observed minus predicted difference (mmol/mol)"
  title <- if (observed) {
    paste0("Matched to ", drug_labels[[reference]], ": observed HbA1c difference by ", check_x_axes[[x_axis]])
  } else {
    paste0("Matched to ", drug_labels[[reference]], ": observed minus predicted HbA1c difference by ",
           check_x_axes[[x_axis]])
  }

  ggplot(groups, aes(x = x, y = y)) +
    # lines and band grow with the text
    { if (observed) geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") } +
    { if (observed) {
      geom_line(data = red_lines, aes(x = x, y = y), inherit.aes = FALSE, colour = "red",
                linewidth = 0.75 * text_scale)
    } else {
      geom_hline(yintercept = 0, colour = "red", linewidth = 0.75 * text_scale)
    } } +
    geom_ribbon(data = curves, aes(x = x, ymin = lower, ymax = upper), inherit.aes = FALSE,
                fill = "blue", alpha = 0.2) +
    geom_line(data = curves, aes(x = x, y = fit), inherit.aes = FALSE, colour = "blue",
              linewidth = 1 * text_scale) +
    # group means: error bars without caps; points and error bars grow with the text
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, colour = "black",
                  linewidth = 0.5 * text_scale) +
    geom_point(size = 1.8 * text_scale) +
    geom_text(data = estimates, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
              size = label_size / ggplot2::.pt, inherit.aes = FALSE) +
    facet_grid(rows = vars(version), cols = vars(drug), labeller = labeller(drug = col_labels)) +
    # zooms only: nothing is dropped
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_classic() +
    theme(
      # axis titles on one line along the height of the panels
      axis.title = element_text(size = axis_title_size(y_title, text_scale, n_rows * (panel_height + 0.3))),
      axis.text = element_text(size = axis_text_size),
      # facet titles without a box
      strip.text.x = element_text(size = col_size),
      strip.text.y = element_text(size = row_size),
      strip.background = element_blank(),
      plot.title = element_text(size = 20 * text_scale),
      plot.title.position = "plot",
      panel.spacing.x = grid::unit(panel_gap, "in"), panel.spacing.y = grid::unit(0.6, "lines"),
      plot.margin = margin(5.5, (overhang[["right"]] + 0.12) * 72, 5.5, 5.5, unit = "pt")
    ) +
    labs(
      x = x_title,
      y = y_title,
      title = wrap_to_width(title, 20 * text_scale, title_width)
    )
}

# Balance (love) plot of the matchings of one figure: absolute standardised mean
# difference per matching covariate, the propensity score and the x variable,
# before and after matching; one panel per treated drug class.
draw_intercept_check_loveplot <- function(matches, reference, x_axis) {
  check_love_labels <- c(love_labels, "distance" = "'Propensity score (logit)'",
                         "x_match" = if (x_axis == "predicted") "'Predicted HbA'[1*c]*' (mmol/mol)'"
                                     else "'Baseline HbA'[1*c]*' (x axis)'")

  love_data <- bind_rows(lapply(matches, function(m) {
    drug <- m$pairs$drug[1]
    s <- summary(m$matchit, addlvariables = "x_match", data = m$data)
    strip <- paste0(drug_labels[[drug]], " as treated\nn = ",
                    format(nrow(m$pairs), big.mark = ",", trim = TRUE), " pairs")
    bind_rows(
      as.data.frame(s$sum.all) %>% rownames_to_column("covariate") %>% mutate(sample = "Before matching"),
      as.data.frame(s$sum.matched) %>% rownames_to_column("covariate") %>% mutate(sample = "After matching")
    ) %>%
      mutate(drug = strip)
  })) %>%
    mutate(
      smd = abs(`Std. Mean Diff.`),
      covariate = factor(covariate, levels = rev(unique(covariate))),
      drug = factor(drug, levels = unique(drug)),
      sample = factor(sample, levels = c("Before matching", "After matching"))
    )

  ggplot(love_data, aes(x = smd, y = covariate, shape = sample, colour = sample)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_vline(xintercept = c(0.05, 0.1), linetype = "dashed", colour = "grey60") +
    geom_point(size = 2.6, na.rm = TRUE) +
    facet_wrap(~ drug, nrow = 1) +
    scale_shape_manual(values = c("Before matching" = 1, "After matching" = 16)) +
    scale_colour_manual(values = c("Before matching" = "grey45", "After matching" = "black")) +
    scale_y_discrete(labels = function(x) {
      out <- check_love_labels[x]
      out[is.na(out)] <- paste0("'", x[is.na(out)], "'")
      parse(text = out)
    }) +
    guides(
      colour = "none",
      shape = guide_legend(override.aes = list(
        shape = 22, size = 7, fill = c("white", "black"), colour = c("grey45", "black")
      ))
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 18),
      strip.text = element_text(size = 16),
      panel.spacing = grid::unit(1.5, "lines")
    ) +
    labs(
      x = "Absolute standardised mean difference",
      y = NULL, shape = NULL, colour = NULL,
      title = paste0("Balance: matched to ", drug_labels[[reference]], " (caliper on ",
                     check_x_axes[[x_axis]], ")")
    )
}


#############################################################################

# Part 5: Run ------------------------------------------------------------------

# helper: one figure: match each treated drug class to reference once (with the
# original predictions), fit each version of the predictions on the same pairs, and
# save the figure (file_name) and the love plot (in Loveplots/). scale as in
# draw_intercept_check(). Returns the estimates and groups, for the csv files.
run_check <- function(reference, treated, x_axis, scale, file_name) {
  treated_drugs <- setdiff(intersect(treated, drugs), reference)
  if (!(reference %in% drugs) || length(treated_drugs) == 0) {
    cat("Figure", file_name, "left out: no initiations of", reference, "or of its treated drug classes\n")
    return(NULL)
  }

  # y plotted and the variable whose slope is reported (observed minus predicted)
  y_var     <- if (scale == "observed") "obs_diff" else "gap"
  slope_var <- if (scale == "observed") "gap_treated" else "gap"

  matches <- lapply(treated_drugs, match_pair, reference = reference, x_axis = x_axis)
  fitted <- lapply(names(check_versions), function(version) {
    lapply(matches, function(m) {
      version_pairs <- pairs_for_version(m$pairs, check_versions[[version]], version)
      c(fit_gap(version_pairs, y = y_var, slope_var = slope_var), list(pairs = version_pairs))
    })
  }) %>% unlist(recursive = FALSE)

  pairs     <- bind_rows(lapply(fitted, `[[`, "pairs"))
  groups    <- bind_rows(lapply(fitted, `[[`, "groups"))
  curves    <- bind_rows(lapply(fitted, `[[`, "curve"))
  estimates <- bind_rows(lapply(fitted, `[[`, "estimate"))

  save_plot(draw_intercept_check(pairs, groups, curves, estimates, reference, x_axis, scale),
            paste0(check_output_path, "/intercept_check_", file_name, ".pdf"),
            panel_width = 3.6, panel_height = 3.0)
  save_plot(draw_intercept_check_loveplot(matches, reference, x_axis),
            paste0(loveplot_path, "/loveplot_", file_name, ".pdf"),
            panel_width = 3.4, panel_height = function(p) 0.3 * n_axis_levels(p, "y"))

  scale_label <- if (scale == "observed") "observed difference" else "observed minus predicted"
  list(estimates = estimates %>% mutate(figure = file_name, scale = scale_label, .before = 1),
       groups    = groups %>% mutate(figure = file_name, scale = scale_label, .before = 1))
}

# Figure 1: the method (the injectable drug classes against oral semaglutide)
results <- list(run_check(method_reference, method_treated, "predicted", "observed",
                          paste0("method_vs_", method_reference)))

# Figure 2: against the drug classes outside the GLP1 family, by each x axis
for (reference in names(check_references)) {
  for (x_axis in names(check_x_axes)) {
    results[[length(results) + 1]] <- run_check(reference, check_references[[reference]], x_axis, "gap",
                                                paste0(x_axis, "_vs_", reference))
  }
}
all_estimates <- lapply(results, `[[`, "estimates")
all_groups    <- lapply(results, `[[`, "groups")

estimates <- as_tibble(bind_rows(all_estimates))
print(estimates %>% select(figure, version, drug, n_pairs, gap, slope_per_10,
                           slope_lower_ci, slope_upper_ci), n = Inf)

# for each figure, version and treated drug class: the mean predicted and observed
# difference (mmol/mol), the mean observed minus predicted difference and its slope
# per 10 mmol/mol of x, with 95% CIs
write.csv(estimates, paste0(check_output_path, "/intercept_check_estimates.csv"), row.names = FALSE)

# the groups of pairs plotted: mean x, mean y (the scale of the figure) and 95% CI
write.csv(bind_rows(all_groups) %>%
            select(figure, scale, reference, x_axis, version, drug, n_groups, group, n_pairs = n, x, y, lower, upper),
          paste0(check_output_path, "/intercept_check_groups.csv"), row.names = FALSE)

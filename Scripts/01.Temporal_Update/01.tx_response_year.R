
#############################################################################

# load required libraries
library(tidyverse)
library(patchwork) # aligning the response and initiation-count panels
library(tableone)
# library(remotes)
# remotes::install_github("PM-Cardoso/cateval")
library(cateval) # Custom R package developed by Pedro for calibration evaluation


# load required libraries
library(aurum)

# set up aurum
cprd = CPRDData$new(cprdEnv = "diabetes-jun2024", cprdConf = "~/.aurum.yaml")
analysis = cprd$analysis("pedro_mm")

## functions ----
is.integer64 <- function(x){
  class(x)=="integer64"
}

# set up dataset (cohort cached by 01.cohort_scripts.R)
analysis_raw <- analysis_raw %>%
  analysis$cached("ethnicity") %>%
  # select variables needed (to shorten the load time)
  select(
    patid, dstartdate,
    gender, dstartdate_age,
    t2dmduration, ethnicity_5cat,
    drug_substance, drug_class,
    imd_decile, smoking_cat,
    prebmi, prehba1c, preegfr, pretotalcholesterol, prehdl, prealt,
    drugline_all, MFN, SGLT2, GLP1, DPP4, TZD, SU, hba1cmonth,
    posthba1cfinal
  ) %>%
  collect() %>%
  mutate(patid=as.character(patid)) %>%
  mutate_if(is.integer64, as.integer)

# Caterogical variables must be properly formatted as factors with the specified levels.
analysis_dataset <- analysis_raw %>%
  mutate(
    pated = paste(patid, dstartdate, drug_class, sep = "."),
    agetx = dstartdate_age,

    # Sex
    sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),

    # Ethnicity
    ethnicity = ifelse(is.na(ethnicity_5cat), 5, ethnicity_5cat),
    ethnicity = factor(ethnicity, levels = c(0:5),
                       labels = c("White", "South Asian", "Black", "Other", "Mixed", "Missing")),

    # Smoking
    smoke = ifelse(is.na(smoking_cat), "Not recorded", smoking_cat),
    smoke = factor(smoke,
                   levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded")),

    # Deprivation (IMD quintile)
    imd5 = ifelse(is.na(imd_decile), 5, imd_decile),
    imd5 = factor(ceiling(imd5/2), levels = c(1, 2, 3, 4, 5),
                  labels = c("1 (least)", "2", "3", "4", "5 (most)")),

    # Number of current treatments
    ncurrtx = MFN + SGLT2 + GLP1 + DPP4 + TZD + SU,
    ncurrtx = ifelse(ncurrtx > 4, 4, ncurrtx),
    ncurrtx = factor(ncurrtx, levels = c(1:4), labels = c("1", "2", "3", "4+")),

    # Drug line
    drugline = ifelse(drugline_all > 5, 5, drugline_all),
    drugline = factor(drugline, levels = c(2:5), labels = c("2", "3", "4", "5+")),

    # Drug class
    drugclass = factor(drug_class,
                       levels = c("SGLT2", "GLP1", "DPP4", "SU", "TZD"))
  ) %>%
  group_by(pated) %>%
  mutate(row = 1:n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  select(-row) %>%
  select(all_of(c(
    "pated", "dstartdate", "agetx", "sex", "t2dmduration", "ethnicity",
    "drug_substance", "drugclass",
    "imd5", "smoke",
    "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
    "drugline", "ncurrtx", "hba1cmonth",
    "posthba1cfinal"
  )
  )) %>%
  as.data.frame()


#############################################################################

# Treatment response by drug group and calendar year ------------------------

# Initiations from this date onwards are analysed and plotted
plot_start <- as.Date("2013-01-01")

# Minimum number of initiations for a drug group / period to be reported
min_n <- 100

# Drug groups: edit the substances in each group to change how they are grouped.
# Substances not listed here are dropped (printed below).
drug_groups <- list(
  "SGLT2i"                 = c("Dapagliflozin", "Empagliflozin", "Canagliflozin", "Ertugliflozin"),
  "DPP4i"                  = c("Sitagliptin", "Saxagliptin", "Alogliptin", "Linagliptin", "Vildagliptin"),
  "SU"                     = c("Gliclazide"),
  "TZD"                    = c("Pioglitazone"),
  "GLP-1RA"                = c("Liraglutide", "Dulaglutide", "Exenatide", "Exenatide prolonged-release"),
  "Oral semaglutide"       = c("Oral semaglutide"),
  "Injectable semaglutide" = c("Low-dose semaglutide", "High-dose semaglutide", "Semaglutide, dose unclear"),
  "Tirzepatide"            = c("Tirzepatide")
)

# GLP-1RA broken down by substance (semaglutide and tirzepatide are their own groups)
glp1_breakdown <- list(
  "Liraglutide"                 = "Liraglutide",
  "Dulaglutide"                 = "Dulaglutide",
  "Exenatide"                   = "Exenatide",
  "Exenatide prolonged-release" = "Exenatide prolonged-release"
)

# First initiation date of the validation cohort, as in
# 02.calibration_new_data_post2020.R
cohort_start <- as.Date("2020-10-14")

# Last period with complete 12-month follow-up in the data. On annual periods
# this is 2022: the data end in mid-2024, so initiations in the second half of
# 2023 do not yet have their 12-month HbA1c and 2023 is not a complete year.
last_period <- "2022"

unassigned <- setdiff(unique(analysis_dataset$drug_substance), unlist(drug_groups))
if (length(unassigned) > 0) message("Substances not in drug_groups (dropped): ", paste(unassigned, collapse = ", "))

analysis_dataset_plot <- analysis_dataset %>%
  filter(dstartdate >= plot_start)

# reference baseline HbA1c at which the adjusted response is reported
ref_prehba1c <- mean(analysis_dataset_plot$prehba1c, na.rm = TRUE)

# 12-month HbA1c adjusted for baseline HbA1c (linear regression, fitted
# separately in each drug group / period), reported at ref_prehba1c;
# response = adjusted 12-month HbA1c - ref_prehba1c. Baseline HbA1c is the main
# determinant of the fall in HbA1c, so without it a drug or year with higher
# baseline values would look more effective. Drug line is not adjusted for here;
# the GLP-1RA comparison below checks that adding it changes nothing.
adjusted_response <- function(d) {
  d <- d %>% drop_na()
  if (nrow(d) < 3) {
    return(tibble(adj_response = NA_real_, adj_response_low = NA_real_, adj_response_high = NA_real_))
  }
  fit <- lm(posthba1cfinal ~ prehba1c, data = d)
  p <- predict(fit, newdata = data.frame(prehba1c = ref_prehba1c), se.fit = TRUE)
  tibble(
    adj_response = p$fit - ref_prehba1c,
    adj_response_low = p$fit - qt(0.975, p$df) * p$se.fit - ref_prehba1c,
    adj_response_high = p$fit + qt(0.975, p$df) * p$se.fit - ref_prehba1c
  )
}

# Summarise initiations and adjusted response by drug group and calendar year.
# groups: named list, group name -> substances
summarise_response_by_period <- function(data, groups) {
  # lookup: substance -> group
  lookup <- stack(groups) %>%
    transmute(drug_substance = as.character(values), drug_group = as.character(ind))
  data %>%
    # a substance may belong to several groups (e.g. pooled GLP-1RA and its own)
    inner_join(lookup, by = "drug_substance", relationship = "many-to-many") %>%
    mutate(
      drug_group = factor(drug_group, levels = names(groups)),
      # 12-month periods: one calendar year
      period_start = make_date(year(dstartdate), 1, 1),
      period = as.character(year(dstartdate))
    ) %>%
    group_by(drug_group, period_start, period) %>%
    summarise(
      n_initiations = n(),
      n_with_outcome = sum(!is.na(posthba1cfinal)),
      adj = list(adjusted_response(pick(prehba1c, posthba1cfinal))),
      .groups = "drop"
    ) %>%
    unnest(adj) %>%
    filter(n_initiations >= min_n) %>%
    arrange(drug_group, period_start)
}

response_by_period <- summarise_response_by_period(analysis_dataset_plot, drug_groups)
response_by_period_glp1_breakdown <- summarise_response_by_period(analysis_dataset_plot, glp1_breakdown)

print(as.data.frame(response_by_period))
print(as.data.frame(response_by_period_glp1_breakdown))

output_path <- "Scripts/01.Temporal_Update/Outputs/treatment_response_by_period"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# ref_prehba1c: the baseline HbA1c the adjusted response is predicted at
write_csv(response_by_period %>% mutate(ref_prehba1c = ref_prehba1c),
          paste0(output_path, "/01.tx_response_by_period.csv"))
write_csv(response_by_period_glp1_breakdown %>% mutate(ref_prehba1c = ref_prehba1c),
          paste0(output_path, "/01.tx_response_by_period_glp1_breakdown.csv"))


## Why GLP1 is split into subtypes ----

# The 5-drug model has a single GLP-1RA class, developed when most GLP-1RA
# initiations were liraglutide or exenatide. Prescribing has since moved to
# dulaglutide. If the substances differ in how much they lower HbA1c, the
# response of the pooled class changes with this mix, so an update estimated for
# the pooled class (an intercept) fits only the substance mix of the cohort it
# is estimated in, and is a property of neither drug. The analysis scripts
# therefore update dulaglutide and liraglutide (with exenatide) separately.
# Three calculations show whether the split is needed:
#   1) the share of each substance among GLP-1RA initiations, before and
#      in the validation cohort (glp1_mix)
#   2) the difference in 12-month HbA1c between substances, adjusted for
#      baseline HbA1c and calendar year, so it is not the change over time
#      (glp1_substance_difference); drug line is added as a check
#   3) the adjusted response of the pooled class at the start and the end of
#      the period (glp1_pooled_change)

complete_label   <- paste0(year(plot_start), " to ", last_period, " (complete follow-up)")
validation_label <- paste0("Validation cohort (from ", format(cohort_start), ")")

glp1_data <- analysis_dataset_plot %>%
  filter(drug_substance %in% unlist(glp1_breakdown)) %>%
  mutate(
    # liraglutide first: the reference in the comparison below
    drug_substance = factor(drug_substance, levels = unlist(glp1_breakdown, use.names = FALSE)),
    year = factor(year(dstartdate))
  )

# 1) substance mix
glp1_mix <- glp1_data %>%
  mutate(cohort = ifelse(dstartdate >= cohort_start, validation_label,
                         paste0(format(plot_start), " to ", format(cohort_start - 1)))) %>%
  count(cohort, drug_substance, name = "n_initiations") %>%
  group_by(cohort) %>%
  mutate(pct_of_glp1 = 100 * n_initiations / sum(n_initiations)) %>%
  ungroup()

# 2) 12-month HbA1c of each substance minus liraglutide (negative: a larger fall
# in HbA1c than liraglutide), from one linear regression with baseline HbA1c and
# calendar year (and, as a check, drug line). Substances with fewer than min_n
# initiations with an outcome are not reported.
substance_difference <- function(data, cohort_label, adjustment, terms) {
  data <- data %>% filter(!is.na(posthba1cfinal)) %>% droplevels()
  fit <- lm(as.formula(paste("posthba1cfinal ~ drug_substance + prehba1c + year", terms)), data = data)
  ci <- confint(fit)
  n_outcome <- table(data$drug_substance)
  tibble(term = names(coef(fit)), estimate = coef(fit), low = ci[, 1], high = ci[, 2]) %>%
    filter(str_starts(term, "drug_substance")) %>%
    transmute(
      cohort = cohort_label,
      adjustment = adjustment,
      substance = sub("^drug_substance", "", term),
      n_with_outcome = as.integer(n_outcome[substance]),
      n_liraglutide_with_outcome = as.integer(n_outcome[["Liraglutide"]]),
      difference_vs_liraglutide = estimate,
      low, high
    ) %>%
    filter(n_with_outcome >= min_n)
}

glp1_complete   <- glp1_data %>% filter(year(dstartdate) <= as.numeric(last_period))
glp1_validation <- glp1_data %>% filter(dstartdate >= cohort_start)

glp1_substance_difference <- bind_rows(
  substance_difference(glp1_complete,   complete_label,   "baseline HbA1c, calendar year", ""),
  substance_difference(glp1_complete,   complete_label,   "baseline HbA1c, calendar year, drug line", "+ drugline"),
  substance_difference(glp1_validation, validation_label, "baseline HbA1c, calendar year", ""),
  substance_difference(glp1_validation, validation_label, "baseline HbA1c, calendar year, drug line", "+ drugline")
)

# 3) adjusted response of the pooled class (all GLP-1RA together, as the
# 5-drug model treats them) in the first and the last two complete years
pooled_years <- response_by_period %>%
  filter(drug_group == "GLP-1RA", !is.na(adj_response), period <= last_period) %>%
  pull(period)
glp1_pooled_change <- response_by_period %>%
  filter(drug_group == "GLP-1RA") %>%
  mutate(years = case_when(period %in% head(pooled_years, 2) ~ paste(head(pooled_years, 2), collapse = "-"),
                           period %in% tail(pooled_years, 2) ~ paste(tail(pooled_years, 2), collapse = "-"))) %>%
  filter(!is.na(years)) %>%
  group_by(years) %>%
  summarise(n_initiations = sum(n_initiations), mean_adj_response = mean(adj_response), .groups = "drop")

print(as.data.frame(glp1_mix))
print(as.data.frame(glp1_substance_difference))
print(as.data.frame(glp1_pooled_change))

write_csv(glp1_mix, paste0(output_path, "/01.glp1_substance_mix.csv"))
write_csv(glp1_substance_difference, paste0(output_path, "/01.glp1_substance_difference.csv"))
write_csv(glp1_pooled_change, paste0(output_path, "/01.glp1_pooled_change.csv"))


## Plot ----

last_period_start <- response_by_period$period_start[response_by_period$period == last_period][1]

# Drug colours (Okabe-Ito). The two blocks below are the groups that ever share
# a single panel; each block has been checked for colourblind separation. The
# previous grey and black slots read as grey (no chroma), and grey sat too close
# to the SU pink under deuteranopia, so semaglutide and tirzepatide take hues.
drug_group_colours <- c(
  # 5-drug model classes
  "SGLT2i"                      = "#E69F00",
  "GLP-1RA"                     = "#56B4E9",
  "DPP4i"                       = "#0072B2",
  "SU"                          = "#CC79A7",
  "TZD"                         = "#D55E00",
  # GLP-1 family
  "Oral semaglutide"            = "#009E73",
  "Injectable semaglutide"      = "#D55E00",
  "Tirzepatide"                 = "#CC79A7",
  # GLP-1RA substances. Exenatide prolonged-release takes the orange that
  # SGLT2i uses in the drug-class figures: the two never appear together, and
  # the pink it used before sat too close to the dulaglutide green for
  # deuteranopes once these five shared a single panel.
  "Liraglutide"                 = "#0072B2",
  "Dulaglutide"                 = "#009E73",
  "Exenatide"                   = "#D55E00",
  "Exenatide prolonged-release" = "#E69F00"
)

# ---- Saving figures ----
# Every figure is saved at a size worked out from its contents: each panel gets
# a fixed size in inches, and the file is exactly as large as the panels plus
# their axes, strips, titles and legend. A figure therefore comes out the same
# on any computer, with no white space around it and no text cut off (text is
# measured with the pdf device's built-in font metrics, not the system fonts).

# facet columns: one row for up to 4 panels, two rows above that
facet_ncol <- function(n_panels) if (n_panels <= 4) n_panels else ceiling(n_panels / 2)

# number of categories on a discrete axis ("x" or "y") of a plot
n_axis_levels <- function(plot, axis = "x") {
  length(ggplot2::ggplot_build(plot)$layout$panel_params[[1]][[axis]]$get_labels())
}

# save a ggplot or patchwork with panels of panel_width x panel_height inches
# (pdf, or png at 300 dpi for slides). Either size can be a function of the
# plot, e.g. to give each bar or each covariate the same room.
save_plot <- function(plot, file, panel_width, panel_height) {
  if (is.function(panel_width))  panel_width  <- panel_width(plot)
  if (is.function(panel_height)) panel_height <- panel_height(plot)
  grDevices::pdf(NULL)
  g <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)

  # panels: the relative ("null") rows and columns get their size in inches
  null_w <- grid::unitType(g$widths) == "null"
  null_h <- grid::unitType(g$heights) == "null"
  g$widths[null_w]  <- grid::unit(panel_width  * as.numeric(g$widths[null_w]), "in")
  g$heights[null_h] <- grid::unit(panel_height * as.numeric(g$heights[null_h]), "in")

  inches <- function(u) grid::convertWidth(u, "in", valueOnly = TRUE)
  col_edges <- c(0, cumsum(vapply(seq_along(g$widths), function(i) inches(g$widths[i]), numeric(1))))
  width <- max(col_edges)

  # titles, captions and legends wider than the plot widen the page instead of
  # being cut off (titles start at their left edge, captions end at their
  # right edge, legends are centred)
  pad_left <- 0; pad_right <- 0
  for (i in seq_len(nrow(g$layout))) {
    name <- g$layout$name[i]
    if (!grepl("^(title|subtitle|caption|guide-box)", name)) next
    grob_width <- inches(grid::grobWidth(g$grobs[[i]]))
    left  <- col_edges[g$layout$l[i]]
    right <- col_edges[g$layout$r[i] + 1]
    if (grepl("^caption", name)) {
      pad_left <- max(pad_left, grob_width - right)
    } else if (grepl("^guide-box", name)) {
      centre <- (left + right) / 2
      pad_left  <- max(pad_left, grob_width / 2 - centre)
      pad_right <- max(pad_right, centre + grob_width / 2 - width)
    } else {
      pad_right <- max(pad_right, left + grob_width - width)
    }
  }
  # a little room on the right for the last x axis label, which overhangs its panel
  pad_right <- max(pad_right, 0.15)
  if (pad_left > 0)  g <- gtable::gtable_add_cols(g, grid::unit(pad_left, "in"), pos = 0)
  if (pad_right > 0) g <- gtable::gtable_add_cols(g, grid::unit(pad_right, "in"), pos = -1)

  width  <- width + pad_left + pad_right
  height <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  grDevices::dev.off()

  if (grepl("\\.png$", file)) {
    ragg::agg_png(file, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  }
  grid::grid.draw(g)
  invisible(grDevices::dev.off())
}

# Periods restricted to those with complete 12-month follow-up, from x_from on.
plot_ready <- function(data, groups, x_from = plot_start) {
  data %>%
    filter(!is.na(adj_response),
           period_start <= last_period_start,
           period_start >= as.Date(x_from),
           drug_group %in% groups) %>%
    mutate(drug_group = factor(drug_group, levels = groups)) %>%
    droplevels()
}

# Five long drug names do not fit on one legend row, so wrap them onto two.
legend_guide <- function(d) {
  guide_legend(nrow = if (nlevels(droplevels(d$drug_group)) > 4) 2 else 1, byrow = TRUE)
}

# Text and mark sizes for print, at journal page width (about 7 inches), scaled
# together so lines and points sit in proportion to the labels. Everything below
# derives from text_base, so raise it (to ~19) for slides and the set scales.
text_base   <- 10
size_axis_text  <- text_base * 0.9
size_axis_title <- text_base
size_legend     <- text_base * 0.9
size_strip      <- text_base
width_line      <- 0.6
size_point      <- 1.5
width_rule      <- 0.4   # the zero line and the cohort-start line

theme_response <- function() {
  theme_classic(base_size = text_base) +
    theme(
      # light horizontal guides only, for reading values off the y axis
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(size = size_axis_text, colour = "black"),
      axis.title       = element_text(size = size_axis_title),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = size_legend),
      legend.key.width = unit(1.6, "lines"),
      legend.margin    = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text       = element_text(size = size_strip, face = "bold", hjust = 0),
      plot.tag         = element_text(size = size_axis_title * 1.1, face = "bold"),
      plot.margin      = margin(4, 8, 4, 4)
    )
}

# The response is the 12-month change in HbA1c predicted at ref_prehba1c. That
# reference value is reported in the CSVs and the figure legend rather than on
# the figures, which carry only axis titles.
y_response_label <- bquote("Adjusted change in HbA"[1*c]*" (mmol/mol)")

# Number of initiations per period, drawn as a short strip sharing the x axis
# with the response panel above it. This replaces the point-size aesthetic:
# size is a weak visual channel, and it was on a scale shared across drugs
# spanning two orders of magnitude. It is also where the precision behind each
# response point is now read, the confidence bands having been dropped.
# count_legend: the faceted figures name each drug in the facet title but give
# the strip no key, so the strip carries the legend there; the overlaid figures
# already have one above the response panel.
count_strip <- function(d, x_scale, log_counts = FALSE, count_legend = FALSE) {
  p <- ggplot(d, aes(period_start, n_initiations, colour = drug_group)) +
    geom_line(linewidth = width_line * 0.85) +
    geom_point(size = size_point * 0.75) +
    x_scale +
    scale_colour_manual(values = drug_group_colours) +
    theme_response() +
    guides(colour = legend_guide(d)) +
    theme(legend.position = if (count_legend) "bottom" else "none") +
    labs(x = "Year of initiation", y = "Initiations")
  if (log_counts) {
    p + scale_y_log10(labels = scales::comma) +
      labs(y = "Initiations\n(log scale)")
  } else {
    p + scale_y_continuous(labels = scales::comma)
  }
}

# One panel per drug group, common y scale so the size of the change stays
# comparable between drugs. Points and lines only: the confidence bands
# overlapped heavily and obscured the trends, and the count strip below carries
# how much each point rests on.
# vline_groups: groups on which the validation cohort start line is drawn.
plot_response_faceted <- function(data, groups, x_from = plot_start,
                                  vline_groups = groups, log_counts = FALSE,
                                  date_breaks = "2 years") {
  d <- plot_ready(data, groups, x_from)
  x_scale <- scale_x_date(date_breaks = date_breaks, date_labels = "%Y",
                          limits = c(as.Date(x_from), max(d$period_start) + 120))
  # only draw the line on panels that have data (groups with < min_n
  # initiations in every period are not shown)
  vline_data <- d %>% distinct(drug_group) %>% filter(drug_group %in% vline_groups)

  p_resp <- ggplot(d, aes(period_start, adj_response, colour = drug_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
               linewidth = width_rule) +
    geom_vline(data = vline_data, aes(xintercept = cohort_start), colour = "black",
               linewidth = width_rule) +
    geom_line(linewidth = width_line) +
    geom_point(size = size_point) +
    x_scale +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_colour_manual(values = drug_group_colours) +
    facet_wrap(~drug_group, ncol = 1) +
    # drug colour is already given by the facet title
    guides(colour = "none") +
    theme_response() +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = y_response_label)

  p_resp / count_strip(d, x_scale, log_counts, count_legend = TRUE) +
    plot_layout(heights = c(4, 1))
}

# All groups overlaid in one panel, for figures where the comparison between
# drugs is the point rather than each drug's own trajectory.
plot_response_overlaid <- function(data, groups, x_from = plot_start,
                                   vline = TRUE, log_counts = FALSE,
                                   date_breaks = "1 year") {
  d <- plot_ready(data, groups, x_from)
  x_scale <- scale_x_date(date_breaks = date_breaks, date_labels = "%Y",
                          limits = c(as.Date(x_from), max(d$period_start) + 120))

  p_resp <- ggplot(d, aes(period_start, adj_response, colour = drug_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
               linewidth = width_rule) +
    geom_line(linewidth = width_line) +
    geom_point(size = size_point) +
    x_scale +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_colour_manual(values = drug_group_colours) +
    guides(colour = legend_guide(d)) +
    theme_response() +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = y_response_label)

  if (vline) {
    p_resp <- p_resp + geom_vline(xintercept = cohort_start, colour = "black",
                                  linewidth = width_rule)
  }

  p_resp / count_strip(d, x_scale, log_counts) + plot_layout(heights = c(3, 1))
}

# 1) Drug classes in the 5-drug model
drugclass_groups <- c("SGLT2i", "DPP4i", "SU", "TZD", "GLP-1RA")

plot_response_by_period_drugclasses <- plot_response_faceted(
  response_by_period,
  groups = drugclass_groups,
  # n runs from ~200 (TZD) to ~17,000 (SGLT2i), so the strip needs a log scale
  log_counts = TRUE
)

save_plot(plot_response_by_period_drugclasses, paste0(output_path, "/01.tx_response_by_period_drugclasses.pdf"),
          panel_width = 6, panel_height = 1.6)

# 1b) The same drug classes overlaid in one panel, for reading the drugs against
# each other rather than each drug's own trajectory
plot_response_by_period_drugclasses_overlaid <- plot_response_overlaid(
  response_by_period,
  groups = drugclass_groups,
  log_counts = TRUE,
  date_breaks = "2 years"
)

save_plot(plot_response_by_period_drugclasses_overlaid, paste0(output_path, "/01.tx_response_by_period_drugclasses_overlaid.pdf"),
          panel_width = 6, panel_height = 1.1)

# 2) GLP-1RA broken down by substance
plot_response_by_period_glp1_breakdown <- plot_response_faceted(
  response_by_period_glp1_breakdown,
  groups = names(glp1_breakdown)
)

save_plot(plot_response_by_period_glp1_breakdown, paste0(output_path, "/01.tx_response_by_period_glp1_breakdown.pdf"),
          panel_width = 6, panel_height = 1.6)

# 2b) Display item: why GLP1 is split into subtypes (see "Why GLP1 is split into
# subtypes" above). Panel A: adjusted response of each substance by year, with
# 95% confidence intervals, as these carry the comparison between substances.
# Panel B: each substance's share of GLP-1RA initiations, the change in
# prescribing that moves the pooled class response. Both stop at the last year
# with complete follow-up.
glp1_share_by_year <- glp1_complete %>%
  mutate(period_start = make_date(year(dstartdate), 1, 1)) %>%
  count(period_start, drug_group = drug_substance, name = "n_initiations") %>%
  group_by(period_start) %>%
  mutate(pct_of_glp1 = 100 * n_initiations / sum(n_initiations)) %>%
  ungroup()

write_csv(glp1_share_by_year, paste0(output_path, "/01.glp1_substance_share_by_year.csv"))

glp1_resp <- plot_ready(response_by_period_glp1_breakdown, names(glp1_breakdown))
glp1_x_scale <- scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                             limits = c(plot_start - 150, last_period_start + 150))
# substances side by side within a year, so the intervals do not overlap
glp1_dodge <- position_dodge(width = 160)

# no zero line: every response is a fall, and a zero line would stretch the
# y axis and shrink the gap between substances, which is what is compared
plot_glp1_split_resp <- ggplot(glp1_resp, aes(period_start, adj_response, colour = drug_group)) +
  geom_vline(xintercept = cohort_start, colour = "black", linewidth = width_rule) +
  annotate("text", x = cohort_start, y = Inf, label = " Validation\n cohort", hjust = 0, vjust = 1.1,
           size = size_axis_text / .pt, colour = "grey20", lineheight = 0.9) +
  geom_line(linewidth = width_line, position = glp1_dodge) +
  geom_errorbar(aes(ymin = adj_response_low, ymax = adj_response_high),
                width = 0, linewidth = width_line * 0.6, position = glp1_dodge) +
  geom_point(size = size_point, position = glp1_dodge) +
  glp1_x_scale +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_colour_manual(values = drug_group_colours) +
  guides(colour = legend_guide(glp1_resp)) +
  theme_response() +
  theme(axis.text.x = element_blank()) +
  labs(x = NULL, y = y_response_label)

plot_glp1_split_share <- ggplot(glp1_share_by_year, aes(period_start, pct_of_glp1, colour = drug_group)) +
  geom_vline(xintercept = cohort_start, colour = "black", linewidth = width_rule) +
  geom_line(linewidth = width_line) +
  geom_point(size = size_point * 0.85) +
  glp1_x_scale +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  scale_colour_manual(values = drug_group_colours) +
  theme_response() +
  theme(legend.position = "none") +
  labs(x = "Year of initiation", y = "Share of GLP-1RA\ninitiations (%)")

plot_response_by_period_glp1_breakdown_overlaid <- (plot_glp1_split_resp / plot_glp1_split_share) +
  plot_layout(heights = c(3, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = size_axis_title, face = "bold"))

save_plot(plot_response_by_period_glp1_breakdown_overlaid, paste0(output_path, "/01.tx_response_by_period_glp1_breakdown_overlaid.pdf"),
          panel_width = 6, panel_height = 1.1)

# Figure legend, with the numbers from this run, for the paper
fmt <- function(x) sprintf("%.1f", x)
difference_text <- function(cohort_label, adjustment_label) {
  d <- glp1_substance_difference %>%
    filter(cohort == cohort_label, adjustment == adjustment_label, substance == "Dulaglutide")
  paste0(fmt(d$difference_vs_liraglutide), " mmol/mol (95% CI ", fmt(d$low), " to ", fmt(d$high), ")")
}
mix_pct <- function(cohort_label, substances) {
  fmt(sum(glp1_mix$pct_of_glp1[glp1_mix$cohort == cohort_label & glp1_mix$drug_substance %in% substances]))
}
before_label <- setdiff(unique(glp1_mix$cohort), validation_label)
exenatide <- c("Exenatide", "Exenatide prolonged-release")

glp1_split_legend <- c(
  "# Figure legend: response to GLP-1 receptor agonists by substance and year",
  "",
  paste0(
    "**(A)** Change in HbA1c 12 months after starting each GLP-1 receptor agonist (GLP-1RA), ",
    "by calendar year of initiation, CPRD Aurum ", year(plot_start), " to ", last_period, ". ",
    "For each substance and year, 12-month HbA1c was regressed on baseline HbA1c and the change is ",
    "reported at the mean baseline HbA1c of ", fmt(ref_prehba1c), " mmol/mol; more negative values are larger falls. ",
    "Bars are 95% confidence intervals. Substance-years with fewer than ", min_n, " initiations are not shown. ",
    "**(B)** Each substance's share of GLP-1RA initiations. Semaglutide (added to the model as its own ",
    "drug classes) and tirzepatide are not included. The vertical line marks the start of the validation cohort (",
    format(cohort_start, "%d %B %Y"), "). ", last_period, " is the last year with complete 12-month follow-up."
  ),
  "",
  paste0(
    "**Why the GLP-1RA class is split.** The 5-drug model has one GLP-1RA class, developed when most ",
    "GLP-1RA initiations were liraglutide or exenatide. Within each substance, the response changed little ",
    "over time (A), but the substances differ: after adjustment for baseline HbA1c and calendar year, the ",
    "difference in 12-month HbA1c, dulaglutide minus liraglutide, was ",
    difference_text(complete_label, "baseline HbA1c, calendar year"), " (",
    difference_text(complete_label, "baseline HbA1c, calendar year, drug line"),
    " when also adjusted for drug line; ", difference_text(validation_label, "baseline HbA1c, calendar year"),
    " within the validation cohort). Prescribing moved from liraglutide and exenatide to dulaglutide (B): ",
    "dulaglutide was ", mix_pct(before_label, "Dulaglutide"), "% of GLP-1RA initiations from ",
    format(plot_start, "%Y"), " until the validation cohort and ", mix_pct(validation_label, "Dulaglutide"),
    "% within it. The adjusted response of the pooled class therefore moved from ",
    fmt(glp1_pooled_change$mean_adj_response[1]), " mmol/mol in ", glp1_pooled_change$years[1], " to ",
    fmt(glp1_pooled_change$mean_adj_response[2]), " mmol/mol in ", glp1_pooled_change$years[2],
    ", largely through the change in mix. An update estimated for the pooled class would describe the ",
    "substance mix of the cohort it was estimated in rather than either drug, so dulaglutide and liraglutide ",
    "are updated separately. Exenatide (", mix_pct(validation_label, exenatide),
    "% of GLP-1RA initiations in the validation cohort) is grouped with liraglutide."
  ),
  "",
  paste0(
    "The comparison between substances is observational and adjusted only for baseline HbA1c and calendar ",
    "year (and drug line as a check), so part of it may reflect differences in who was prescribed each drug."
  )
)

writeLines(glp1_split_legend, paste0(output_path, "/01.glp1_split_figure_legend.md"))

# 3) GLP-1RA vs semaglutide / tirzepatide, overlaid so the drugs can be
# compared directly. Starts in 2018 rather than 2013: semaglutide only appears
# from 2019 (injectable) and 2021 (oral), so a 2013 start left those series
# squeezed into the right-hand third of the panel.
plot_response_by_period_glp1 <- plot_response_overlaid(
  response_by_period,
  groups = c("GLP-1RA", "Oral semaglutide", "Injectable semaglutide", "Tirzepatide"),
  x_from = "2018-01-01"
)

save_plot(plot_response_by_period_glp1, paste0(output_path, "/01.tx_response_by_period_glp1.pdf"),
          panel_width = 6, panel_height = 1.1)

# 4) Change in adjusted response relative to each drug's own 2013-14 level.
# The absolute figures above share a ~22 mmol/mol y axis because the drugs sit
# at different response levels, which leaves the 3-9 mmol/mol within-drug change
# hard to see. Referencing each drug to its own starting level puts them on a
# ~10 mmol/mol axis in one panel, so the changes can be read and compared.
# Note that the baseline is itself estimated, so this understates uncertainty.
baseline_end <- as.Date("2015-01-01")

response_change_from_baseline <- plot_ready(response_by_period, drugclass_groups) %>%
  arrange(drug_group, period_start) %>%
  group_by(drug_group) %>%
  mutate(baseline_response = mean(adj_response[period_start < baseline_end]),
         change_from_baseline = adj_response - baseline_response) %>%
  ungroup()

write_csv(response_change_from_baseline,
          paste0(output_path, "/01.tx_response_change_from_baseline.csv"))

plot_response_change_from_baseline <- response_change_from_baseline %>%
  ggplot(aes(period_start, change_from_baseline, colour = drug_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = width_rule) +
  geom_vline(xintercept = cohort_start, colour = "black",
             linewidth = width_rule) +
  geom_line(linewidth = width_line) +
  geom_point(size = size_point) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = drug_group_colours) +
  guides(colour = legend_guide(response_change_from_baseline)) +
  theme_response() +
  labs(x = "Year of initiation",
       y = bquote("Change vs 2013-14 (mmol/mol)"))

save_plot(plot_response_change_from_baseline, paste0(output_path, "/01.tx_response_change_from_baseline.pdf"),
          panel_width = 6, panel_height = 3)

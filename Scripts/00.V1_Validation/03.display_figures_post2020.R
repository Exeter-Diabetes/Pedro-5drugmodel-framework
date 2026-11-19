#############################################################################
## Display figures of the framework paper (worked example)
##
## Runs 02.validation_v1_post2020.R, then saves the figures of the paper from its
## results, drawn with its figure code but without the figure titles (they are in
## the legends), as pdf and png (300 dpi):
##   Figures/Manuscript:    figure2 (average benefit of following the model),
##                          figure3 (concordant-discordant benefit calibration,
##                          selected approach), figure4 (head-to-head calibration,
##                          selected approach)
##   Figures/Supplementary: figureS2-S5 (balance of the matched pairs),
##                          figureS6, S8, S10 (concordant-discordant benefit
##                          calibration), figureS7, S9, S11 (head-to-head
##                          calibration), figureS12-S15 (predicted and observed
##                          12-month HbA1c)
## Figure 1 (the framework) and figure S1 (flow chart) are drawn separately.
## Run from the repository root.
#############################################################################

source("Scripts/00.V1_Validation/02.validation_v1_post2020.R")

fig_main <- "Figures/Manuscript"
fig_supp <- "Figures/Supplementary"
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)

# the selected approach of the worked example
selected <- recal_minimal

# the figure without its title and subtitle
no_title <- function(p) {
  if (inherits(p, "patchwork")) p + patchwork::plot_annotation(title = NULL, subtitle = NULL)
  else p + labs(title = NULL, subtitle = NULL)
}

# save a figure as pdf and png with one of the save_*() functions of the validation
save_both <- function(save_fun, fig, dir, name) {
  for (ext in c(".pdf", ".png")) save_fun(fig, file.path(dir, paste0(name, ext)))
}

# the gtable and size (inches) that save_plot() draws, for figures made of parts
sized_gtable <- function(plot, panel_width, panel_height) {
  grDevices::pdf(NULL)
  g <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)
  null_w <- grid::unitType(g$widths) == "null"
  null_h <- grid::unitType(g$heights) == "null"
  g$widths[null_w]  <- grid::unit(panel_width  * as.numeric(g$widths[null_w]), "in")
  g$heights[null_h] <- grid::unit(panel_height * as.numeric(g$heights[null_h]), "in")
  width  <- sum(vapply(seq_along(g$widths), function(i) grid::convertWidth(g$widths[i], "in", valueOnly = TRUE), numeric(1))) + 0.15
  g <- gtable::gtable_add_cols(g, grid::unit(0.15, "in"), pos = -1)
  height <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  grDevices::dev.off()
  list(g = g, width = width, height = height)
}

# save parts one above the other, labelled A, B, ... (pdf or png by extension)
save_stacked <- function(parts, file) {
  tag_height <- 0.35
  width  <- max(vapply(parts, function(p) p$width, numeric(1)))
  height <- sum(vapply(parts, function(p) p$height, numeric(1))) + tag_height * length(parts)
  if (grepl("[.]png$", file)) {
    ragg::agg_png(file, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  }
  grid::grid.newpage()
  top <- height
  for (i in seq_along(parts)) {
    grid::grid.text(LETTERS[i], x = grid::unit(0.1, "in"), y = grid::unit(top - tag_height / 2, "in"),
                    just = "left", gp = grid::gpar(fontsize = 22, fontface = "bold"))
    top <- top - tag_height
    grid::pushViewport(grid::viewport(x = grid::unit(0, "in"), y = grid::unit(top, "in"),
                                      width = grid::unit(parts[[i]]$width, "in"),
                                      height = grid::unit(parts[[i]]$height, "in"), just = c("left", "top")))
    grid::grid.draw(parts[[i]]$g)
    grid::popViewport()
    top <- top - parts[[i]]$height
  }
  invisible(grDevices::dev.off())
}

# concordant-discordant benefit calibration of one version: (A) the three matched
# populations, (B) by recommended drug class
save_benefit_figure <- function(version, dir, name) {
  parts <- list(sized_gtable(no_title(version$benefit$plot), 4.2, 4.2),
                sized_gtable(no_title(version$benefit$plot_by_drug), 3.4, 3.4))
  for (ext in c(".pdf", ".png")) save_stacked(parts, file.path(dir, paste0(name, ext)))
}


## Main figures ----

# Figure 2: average benefit of following the model under each approach
save_both(save_benefit_estimates,
          no_title(draw_benefit_estimates(benefit_estimates %>%
                                            mutate(analysis_label = "drug classes analysed", key_result = TRUE))),
          fig_main, "figure2")

# Figure 3: concordant-discordant benefit calibration, selected approach
save_benefit_figure(selected, fig_main, "figure3")

# Figure 4: head-to-head calibration, selected approach
save_both(save_pairwise, no_title(plot_pairwise(selected$pairwise, "Pairwise Drug Calibration", NULL)),
          fig_main, "figure4")


## Supplementary figures ----

versions <- list(original = original, recal_minimal = recal_minimal, recal_large = recal_large, recal = recal)
prefixes <- c("Original predictions" = "pred.", "Minimal recalibration-in-the-large" = "pred.recal_minimal.",
              "Recalibration-in-the-large" = "pred.recal_large.", "Recalibration" = "pred.recal.")

# Figures S2-S5: balance of the matched pairs
for (i in seq_along(versions)) {
  for (ext in c(".pdf", ".png")) {
    save_loveplot(versions[[i]]$benefit$matched, file.path(fig_supp, paste0("figureS", i + 1, ext)))
  }
}

# Figures S6-S11: concordant-discordant benefit calibration and head-to-head
# calibration of the approaches not shown in the main figures
for (v in c("original", "recal_large", "recal")) {
  number <- c(original = 6, recal_large = 8, recal = 10)[[v]]
  save_benefit_figure(versions[[v]], fig_supp, paste0("figureS", number))
  save_both(save_pairwise, no_title(plot_pairwise(versions[[v]]$pairwise, "Pairwise Drug Calibration", NULL)),
            fig_supp, paste0("figureS", number + 1))
}

# Figures S12-S15: predicted and observed 12-month HbA1c
for (i in seq_along(prefixes)) {
  save_both(save_pred_obs, no_title(plot_pred_obs(prefixes[[i]], names(prefixes)[i])), fig_supp, paste0("figureS", i + 11))
}

#' Plotting optimal drug combinations
#'
#'
#' This function generates a side-by-side bar plot visualizing:
#' 1. The proportion of individuals with 1 to 5 optimal drugs identified.
#' 2. The most frequent drug combinations by category (e.g., 1-drug, 2-drug).
#'
#' @param data A character vector where each element is a comma-separated string of drugs considered optimal for an individual.
#' @param groups A named list where each name represents a category (e.g., "1-drug combination") and each value is a vector of drug counts (e.g., `1`, `2`, `4:5`).
#'
#' @return A patchwork plot combining two ggplots.
#'
#' @examples
#' groups <- list(
#'   "1-drug combination" = 1,
#'   "2-drug combinations" = 2,
#'   "3-drug combinations" = 3,
#'   "3/4/5-drug combinations" = 4:5
#' )
#' optimal_drug_comparison_plot(c("SGLT2,DPP4,GLP1", "GLP1", "TZD,DPP4"), groups)
#'
optimal_drug_comparison_plot <- function(data, groups) {
  
  # Summarize number of drugs per individual ----
  
  summary_number <- data.frame(drugs = data) %>%
    # Count how many drugs each person has, based on commas
    mutate(drug_count = if_else(drugs == "", 0L, str_count(drugs, ",") + 1)) %>%
    # Filter out those with zero drugs
    filter(drug_count != 0) %>%
    # Count frequency of each drug count
    count(drug = drug_count, name = "Count") %>%
    # Calculate proportion relative to total data
    mutate(Percentage = round(Count / length(data), 4))
  
  # Summarize most common drug combinations by group ----
  
  summary_combinations <- purrr::map_dfr(
    names(groups),
    function(label) {
      tibble(drugs = data) %>%
        # Recalculate number of drugs per row
        mutate(drug_count = if_else(drugs == "", 0L, str_count(data, ",") + 1)) %>%
        # Keep only drugs matching group definition
        filter(drug_count %in% groups[[label]]) %>%
        # Count unique drug combinations
        count(drug = drugs, name = "Count") %>%
        # Add percentage and group label
        mutate(
          Percentage = round(Count / length(data), 4),
          Combinations = label
        )
    }
  ) %>%
    # Order within group by percentage descending
    group_by(Combinations) %>%
    arrange(Combinations, desc(Percentage), .by_group = TRUE) %>%
    # Calculate group-level percentage (for info only)
    mutate(Group_Percentage = round(sum(Count) / length(data), 4)) %>%
    ungroup() %>%
    # Filter out combinations with less than 0.1%
    filter(Percentage >= 0.001)
  
  # Plot number of drugs ----
  
  plot_a <- summary_number %>%
    ggplot(aes(x = reorder(drug, rev(drug)), y = Percentage)) +
    geom_text(aes(label = drug), 
              hjust = -0.5, size = 8) +
    geom_col(fill = "#076fa2") +
    geom_hline(aes(yintercept = 0)) +
    scale_y_continuous(labels = scales::percent, 
                       breaks = seq(0, max(summary_number$Percentage) + 0.05, by = 0.05)) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Number of predicted optimal drug",
         subtitle = "Proportion of individuals with optimal drugs (>3 mmol/mol benefit)") + 
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank()
    )
  
  # Prepare combination plot data ----
  
  # Define factor order based on original group names
  combi_order <- names(groups)
  
  df2 <- summary_combinations %>%
    # Convert Combinations to ordered factor
    mutate(Combinations = factor(Combinations, levels = combi_order)) %>%
    # Order drugs within each group
    group_by(Combinations) %>%
    arrange(desc(Percentage), .by_group = TRUE) %>%
    ungroup() %>%
    # Reorder drug levels for y-axis plotting
    mutate(drug = factor(drug, levels = rev(unique(drug)))) %>%
    arrange(Combinations, desc(Percentage)) %>%
    mutate(drug = factor(drug, levels = rev(drug)))
  
  # Plot combinations by group ----
  
  plot_b <- ggplot(df2, aes(x = Percentage, y = drug)) +
    geom_text(aes(label = scales::percent(Percentage, accuracy = 0.1)), 
              hjust = -0.1, size = 3) +
    geom_col(fill = "#076fa2") +
    labs(title = "Predicted HbA1c-optimal drug classes") +
    scale_x_continuous(
      labels = scales::percent,
      breaks = seq(0, max(df2$Percentage) + 0.05, by = 0.05),
      expand = expansion(mult = c(0, 0.1))
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank()
    )
  
  # Combine plots side-by-side ----
  
  final_plot <- patchwork::wrap_plots(plot_a, plot_b, ncol = 2, nrow = 1)
  
  # Return final plot
  return(final_plot)
}

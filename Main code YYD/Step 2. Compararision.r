# ============================================================
# Supplementary analysis: comparison of community structure
# between observed and simulated forest quadrats
# ============================================================
# This script summarises community-level structural metrics
# across 110 forest quadrats and produces the supplementary
# figures and tables used to compare observed and simulated
# community structure.
#
# Input:
#   110quadrat(Gini_LAC).xlsx
#
# Outputs:
#   Fig. S2  - Paired comparison of observed vs. simulated
#              spatial Gini coefficients
#   Fig. S3  - Correlation between delta_Gini_spatial and Gini_size
#   Table S1 - Descriptive statistics for Gini metrics
#   Table S2 - Descriptive statistics for LAC metrics
#
# Key metrics:
#   Gini_spatial_obs   = spatial Gini from observed tree positions
#   Gini_spatial_sim   = spatial Gini from simulated tree positions
#   delta_Gini_spatial = Gini_spatial_obs - Gini_spatial_sim
#   Gini_size          = size inequality based on individual size
#   LAC_spatial_obs    = Lorenz asymmetry coefficient for spatial structure
#   LAC_size           = Lorenz asymmetry coefficient for size structure
# ============================================================

rm(list = ls())

# ---------------------------
# 0. Packages and settings
# ---------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(writexl)
library(ggplot2)
library(scales)

save_dir <- getwd()

plot_colors <- list(
  main   = "#2E3440",
  accent = "#BF616A",
  paired = c("#1F77B4", "#E31A1C")
)


# ---------------------------
# 1. Read pre-computed metrics
# ---------------------------
# Each row represents one forest quadrat.
data_path <- "110quadrat(Gini_LAC).xlsx"
df_metrics <- read_excel(data_path)


# ---------------------------
# 2. Paired comparison plot
# ---------------------------
# This function compares two paired variables using a paired
# Wilcoxon signed-rank test and returns a violin-boxplot figure.
plot_paired_comparison <- function(data, col1, col2,
                                   name1, name2,
                                   y_label, colors) {

  df_pair <- data %>%
    dplyr::select(all_of(c(col1, col2))) %>%
    tidyr::drop_na()

  vec1 <- df_pair[[col1]]
  vec2 <- df_pair[[col2]]
  n_samples <- nrow(df_pair)

  wilcox_res <- wilcox.test(vec1, vec2, paired = TRUE)
  p_val <- wilcox_res$p.value

  sig_label <- dplyr::case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  )

  p_label <- dplyr::case_when(
    p_val < 0.001 ~ "p < 0.001",
    p_val < 0.01  ~ "p < 0.01",
    p_val < 0.05  ~ "p < 0.05",
    TRUE          ~ sprintf("p = %.3f", p_val)
  )

  df_long <- df_pair %>%
    rename(
      !!name1 := all_of(col1),
      !!name2 := all_of(col2)
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Type",
      values_to = "value"
    ) %>%
    mutate(Type = factor(Type, levels = c(name1, name2)))

  y_max   <- max(df_long$value, na.rm = TRUE)
  y_min   <- min(df_long$value, na.rm = TRUE)
  y_range <- y_max - y_min

  if (y_range == 0) y_range <- 0.05

  y_line  <- y_max + 0.05 * y_range
  y_star  <- y_max + 0.11 * y_range
  y_ptext <- y_max + 0.17 * y_range
  y_ntext <- y_max + 0.23 * y_range

  ggplot(df_long, aes(x = Type, y = value, fill = Type)) +
    geom_violin(
      width = 0.6, alpha = 0.7, scale = "width",
      trim = TRUE, color = "white"
    ) +
    geom_boxplot(
      width = 0.15, alpha = 0.9, fatten = 1.5,
      outlier.shape = NA, color = "black"
    ) +
    geom_jitter(
      width = 0.08, alpha = 0.5, size = 1.2,
      color = "grey30"
    ) +
    scale_fill_manual(values = setNames(colors, c(name1, name2))) +
    labs(x = NULL, y = y_label) +
    geom_segment(
      aes(x = 1, xend = 2, y = y_line, yend = y_line),
      inherit.aes = FALSE, color = "black", linewidth = 0.8
    ) +
    annotate(
      "text", x = 1.5, y = y_star,
      label = sig_label, size = 6, fontface = "bold"
    ) +
    annotate(
      "text", x = 1.5, y = y_ptext,
      label = p_label, size = 3.5, color = "grey30"
    ) +
    annotate(
      "text", x = 1.5, y = y_ntext,
      label = sprintf("n = %d", n_samples),
      size = 3.5, color = "grey30"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title.y = element_text(face = "bold", size = 13),
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      panel.border = element_rect(color = "black", linewidth = 1.2),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )
}


# ---------------------------
# 3. Fig. S2
# ---------------------------
# Paired comparison of observed and simulated spatial Gini coefficients.
plot_gini <- plot_paired_comparison(
  data    = df_metrics,
  col1    = "Gini_spatial_obs",
  col2    = "Gini_spatial_sim",
  name1   = "Observed",
  name2   = "Simulated",
  y_label = "Spatial Gini coefficient",
  colors  = plot_colors$paired
)

print(plot_gini)


# ---------------------------
# 4. Fig. S3
# ---------------------------
# Spearman correlation between delta_Gini_spatial and Gini_size.
df_corr <- df_metrics %>%
  dplyr::select(delta_Gini_spatial, Gini_size) %>%
  tidyr::drop_na()

spearman_test <- cor.test(
  df_corr$delta_Gini_spatial,
  df_corr$Gini_size,
  method = "spearman"
)

print(spearman_test)

is_significant <- spearman_test$p.value < 0.05

if (is_significant) {
  lm_model  <- lm(Gini_size ~ delta_Gini_spatial, data = df_corr)
  r_squared <- summary(lm_model)$r.squared
}

p_scatter <- ggplot(df_corr, aes(x = delta_Gini_spatial, y = Gini_size)) +
  geom_point(
    color = plot_colors$main,
    alpha = 0.5, size = 2.5, shape = 16
  ) +
  {if (is_significant) {
    geom_smooth(
      method = "lm",
      color = plot_colors$accent,
      fill = plot_colors$accent,
      linewidth = 1,
      alpha = 0.15,
      se = TRUE
    )
  }} +
  annotate(
    "text", x = Inf, y = Inf,
    label = if (is_significant) {
      sprintf(
        "\u03c1 = %.3f\np %s\nR\u00b2 = %.3f\nn = %d",
        unname(spearman_test$estimate),
        ifelse(spearman_test$p.value < 0.001, "< 0.001",
               sprintf("= %.3f", spearman_test$p.value)),
        r_squared,
        nrow(df_corr)
      )
    } else {
      sprintf(
        "\u03c1 = %.3f\np = %.3f\nn = %d",
        unname(spearman_test$estimate),
        spearman_test$p.value,
        nrow(df_corr)
      )
    },
    hjust = 1.08, vjust = 1.1, size = 3.5, lineheight = 0.9
  ) +
  labs(
    x = expression(bold(Delta*"Gini"[spatial])),
    y = expression(bold("Gini"[size]))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05)),
    breaks = pretty_breaks(n = 5)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05)),
    breaks = pretty_breaks(n = 5)
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(color = plot_colors$main, linewidth = 0.75),
    axis.ticks = element_line(color = plot_colors$main, linewidth = 0.75),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(10, 15, 10, 10, "pt"),
    legend.position = "none"
  )

print(p_scatter)


# ---------------------------
# 5. Table S1
# ---------------------------
# Descriptive statistics for the three Gini metrics.
table_s1_gini <- df_metrics %>%
  dplyr::select(Gini_spatial_obs, Gini_spatial_sim, Gini_size) %>%
  tidyr::drop_na() %>%
  summarise(
    across(
      everything(),
      list(
        Mean = ~ round(mean(.x), 3),
        SD   = ~ round(sd(.x), 3),
        Min  = ~ round(min(.x), 3),
        Max  = ~ round(max(.x), 3),
        N    = ~ n()
      ),
      .names = "{.col}__{.fn}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Variable", "Statistic"),
    names_sep = "__"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = value
  ) %>%
  mutate(
    Component = case_when(
      Variable == "Gini_spatial_obs" ~ "Observed Gini_spatial",
      Variable == "Gini_spatial_sim" ~ "Simulated Gini_spatial",
      Variable == "Gini_size"        ~ "Gini_size",
      TRUE                           ~ Variable
    ),
    `Mean ± SD` = sprintf("%.3f ± %.3f", Mean, SD)
  ) %>%
  dplyr::select(Component, Mean, SD, `Mean ± SD`, Min, Max, N)

cat("\n--- Table S1: Gini summary ---\n")
print(table_s1_gini)


# ---------------------------
# 6. Table S2
# ---------------------------
# Descriptive statistics for LAC metrics.
table_s2_lac <- df_metrics %>%
  dplyr::select(LAC_spatial_obs, LAC_size) %>%
  tidyr::drop_na() %>%
  summarise(
    Mean_spatial = round(mean(LAC_spatial_obs), 2),
    SD_spatial   = round(sd(LAC_spatial_obs), 2),
    Mean_size    = round(mean(LAC_size), 2),
    SD_size      = round(sd(LAC_size), 2)
  ) %>%
  {
    tibble(
      Component = c("LAC_spatial_obs", "LAC_size"),
      Mean = c(.$Mean_spatial, .$Mean_size),
      SD   = c(.$SD_spatial, .$SD_size),
      `Mean ± SD` = sprintf(
        "%.2f ± %.2f",
        c(.$Mean_spatial, .$Mean_size),
        c(.$SD_spatial, .$SD_size)
      )
    )
  }

cat("\n--- Table S2: LAC summary ---\n")
print(table_s2_lac)


# ---------------------------
# 7. Save outputs
# ---------------------------
ggsave(
  file.path(save_dir, "Fig. S2 Spatial_Gini_Comparison.pdf"),
  plot_gini,
  width = 89, height = 100,
  units = "mm", dpi = 600, device = cairo_pdf
)

ggsave(
  file.path(save_dir, "Fig. S3 Gini_Correlation.pdf"),
  p_scatter,
  width = 89, height = 89,
  units = "mm", dpi = 600, device = cairo_pdf
)

writexl::write_xlsx(
  list(
    "Table_S1_Gini_Summary" = table_s1_gini,
    "Table_S2_LAC_Summary"  = table_s2_lac
  ),
  path = file.path(save_dir, "TableS1_S2_Summary.xlsx")
)

cat("\nAll outputs saved to:", save_dir, "\n")

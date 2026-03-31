# ============================================================
#
#  Supplementary Analysis: Community Structure Comparison
#  Between Observed and Simulated Forest Plots
#
#  Scientific Context:
#  We test whether the spatial size inequality of individual
#  trees (measured by Spatial Gini Coefficient) significantly
#  differs between empirical observations (real forest data)
#  and a null model (complete spatial randomness simulation).
#  A significant difference indicates non-random spatial
#  aggregation of tree sizes — a key signal of ecological
#  structure in subtropical forests.
#
#  Outputs:
#    · Fig. S2  — Violin + boxplot: Observed vs. Simulated
#                 Spatial Gini Coefficient (n = 110 quadrats)
#    · Fig. S3  — Scatter plot: ΔGini_spatial vs. Gini_size
#                 Includes linear regression and Spearman correlation
#    · Table S1 — Mean ± SD of LAC (Spatial vs. Size-based)
#
#  Statistical tests:
#    - Paired Wilcoxon signed-rank test (Fig. S2)
#    - Spearman rank correlation + linear regression (Fig. S3)
#
#  Authors: [Your Name]
#  Last updated: 2026-03-30
#
# ============================================================


# ------------------------------------------------------------
# 0. Load required packages & global settings
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(writexl)
library(ggplot2)
library(scales)

# Output directory for all figures and tables
save_dir <- "F:/DYY-博士论文/1. QY_SAR/QY_SAR_PROJECT/5. 小论文工作/5.4 性状数据的获取/QYS_QYS_物种名称纠正/小论文结果数据1/"

# Shared color palette
nature_colors <- list(
  main   = "#2E3440",   # Dark blue-grey: axes, points
  accent = "#BF616A"    # Reddish: regression line / emphasis
)


# ------------------------------------------------------------
# 1. Load pre-computed community structure metrics
#
#    Each row = one 20×20 m forest quadrat (n = 110 total)
#    Key columns:
#      · Gini_spatial_obs   — Spatial Gini from real tree positions
#      · Gini_spatial_sim   — Spatial Gini under random simulation
#      · delta_Gini_spatial — Gini_spatial_obs − Gini_spatial_sim
#      · Gini_size          — Gini coefficient of DBH-based size
#      · LAC_spatial_obs    — Lorenz Asymmetry Coeff. (spatial)
#      · LAC_size           — Lorenz Asymmetry Coeff. (DBH-based)
# ------------------------------------------------------------
data_path <- "110quadrat(Gini_LAC).xlsx"

df_metrics <- read_excel(data_path)

# Quick sanity check
df_metrics %>% filter(plot_id == 74)


# ------------------------------------------------------------
# 2. Core plotting function — violin + boxplot comparison
#
#    Reusable for any paired two-group comparison.
#    Steps:
#      (a) Reshape two columns → long format for ggplot
#      (b) Run paired Wilcoxon test and extract p-value
#      (c) Draw violin + boxplot + jittered points
#      (d) Annotate significance bracket (line + stars + p + n)
#      (e) Return ggplot object (saving handled in Section 6)
#
#    Arguments:
#      data        — input data frame (df_metrics)
#      col1/col2   — column names of the two groups
#      name1/name2 — display labels on the X axis
#      y_label     — Y axis title
#      colors      — named color vector (length 2)
# ------------------------------------------------------------
plot_eco_comparison <- function(data, col1, col2,
                                name1, name2,
                                y_label, colors) {

  # --- 2a. Reshape to long format --------------------------
  df_long <- data %>%
    dplyr::select(all_of(c(col1, col2))) %>%
    drop_na() %>%
    rename(!!name1 := all_of(col1),
           !!name2 := all_of(col2)) %>%
    pivot_longer(cols      = everything(),
                 names_to  = "Type",
                 values_to = "value")

  # Fix factor order so ggplot respects left → right layout
  df_long$Type <- factor(df_long$Type, levels = c(name1, name2))

  # --- 2b. Paired Wilcoxon signed-rank test ----------------
  vec1 <- data[[col1]][!is.na(data[[col1]]) & !is.na(data[[col2]])]
  vec2 <- data[[col2]][!is.na(data[[col1]]) & !is.na(data[[col2]])]

  wilcox_res <- wilcox.test(vec1, vec2, paired = TRUE)
  p_val      <- wilcox_res$p.value
  n_samples  <- length(vec1)

  sig_symbol <- dplyr::case_when(
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

  # --- 2c. Y positions for significance bracket ------------
  y_max   <- max(df_long$value, na.rm = TRUE)
  y_range <- y_max - min(df_long$value, na.rm = TRUE)
  y_line  <- y_max + 0.05 * y_range
  y_star  <- y_max + 0.11 * y_range
  y_ptext <- y_max + 0.17 * y_range
  y_ntext <- y_max + 0.23 * y_range

  # --- 2d. Build figure ------------------------------------
  p <- ggplot(df_long, aes(x = Type, y = value, fill = Type)) +

    geom_violin(width = 0.6, alpha = 0.7, scale = "width",
                trim  = TRUE, color = "white") +

    geom_boxplot(width = 0.15, alpha = 0.9, fatten = 1.5,
                 outlier.shape = NA, color = "black") +

    geom_jitter(width = 0.08, alpha = 0.5, size = 1.2,
                color = "grey30") +

    scale_fill_manual(values = setNames(colors, c(name1, name2))) +

    labs(y = y_label, x = "") +

    # Significance bracket
    geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line),
                 color = "black", linewidth = 0.8, inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = y_star,
             label = sig_symbol, size = 6, fontface = "bold") +
    annotate("text", x = 1.5, y = y_ptext,
             label = p_label,   size = 3.5, color = "grey30") +
    annotate("text", x = 1.5, y = y_ntext,
             label = sprintf("n = %d", n_samples),
             size = 3.5, color = "grey30") +

    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid      = element_blank(),
      axis.title.y    = element_text(face = "bold", size = 13),
      axis.text.x     = element_text(face = "bold", size = 12,
                                     color = "black"),
      axis.text.y     = element_text(size = 11, color = "black"),
      panel.border    = element_rect(color = "black", linewidth = 1.2),
      axis.ticks      = element_line(color = "black", linewidth = 0.8)
    )

  return(p)
}


# ------------------------------------------------------------
# 3. Fig. S2 — Spatial Gini: Observed vs. Simulated
#
#    Blue = Observed (real tree positions)
#    Red  = Simulated (complete spatial randomness null model)
# ------------------------------------------------------------
color_gini <- c("#1F77B4",   # Blue → Observed
                "#E31A1C")   # Red  → Simulated

plot_gini <- plot_eco_comparison(
  data    = df_metrics,
  col1    = "Gini_spatial_obs",
  col2    = "Gini_spatial_sim",
  name1   = "Observed",
  name2   = "Simulated",
  y_label = "Spatial Gini Coefficient",
  colors  = color_gini
)

print(plot_gini)


# ------------------------------------------------------------
# 4. Fig. S3 — ΔGini_spatial vs. Gini_size
#
#    Spearman rank correlation to assess whether plots with
#    greater spatial size inequality (ΔGini_spatial) also
#    exhibit greater overall size inequality (Gini_size).
# ------------------------------------------------------------

# --- 4a. Spearman correlation test -----------------------
spearman_test <- cor.test(df_metrics$delta_Gini_spatial,
                          df_metrics$Gini_size,
                          method = "spearman")
print(spearman_test)

is_significant <- spearman_test$p.value < 0.05

# --- 4b. Linear regression (only if significant) ---------
if (is_significant) {
  lm_model  <- lm(Gini_size ~ delta_Gini_spatial, data = df_metrics)
  r_squared <- summary(lm_model)$r.squared
}

# --- 4c. Build scatter plot ------------------------------
p_scatter <- ggplot(df_metrics, aes(x = delta_Gini_spatial, y = Gini_size)) +

  geom_point(color = nature_colors$main,
             alpha = 0.5, size = 2.5, shape = 16) +

  {if (is_significant) {
    geom_smooth(method    = "lm",
                color     = nature_colors$accent,
                fill      = nature_colors$accent,
                linewidth = 1, alpha = 0.15, se = TRUE)
  }} +

  annotate("text", x = Inf, y = Inf,
           label = if (is_significant) {
             sprintf("\u03c1 = %.3f\nP %s\nR\u00b2 = %.3f\nn = %d",
                     spearman_test$estimate,
                     ifelse(spearman_test$p.value < 0.001, "< 0.001",
                            sprintf("= %.3f", spearman_test$p.value)),
                     r_squared,
                     nrow(df_metrics))
           } else {
             sprintf("\u03c1 = %.3f\nP = %.3f\nn = %d",
                     spearman_test$estimate,
                     spearman_test$p.value,
                     nrow(df_metrics))
           },
           hjust = 1.08, vjust = 1.1, size = 3.5, lineheight = 0.9) +

  labs(
    x = expression(bold(Delta*"Gini"[spatial])),
    y = expression(bold("Gini"[size]))
  ) +

  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     breaks = pretty_breaks(n = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     breaks = pretty_breaks(n = 5)) +

  theme_classic(base_size = 11) +
  theme(
    axis.line         = element_line(color = nature_colors$main, linewidth = 0.75),
    axis.ticks        = element_line(color = nature_colors$main, linewidth = 0.75),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 9),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.grid        = element_blank(),
    plot.margin       = margin(10, 15, 10, 10, "pt"),
    legend.position   = "none"
  )

print(p_scatter)


# ------------------------------------------------------------
# 5. Table S1 — LAC: Spatial vs. Size-based
#
#    Distinguishes whether inequality arises from the spatial
#    arrangement or from the underlying size hierarchy.
# ------------------------------------------------------------
table_s1 <- df_metrics %>%
  dplyr::select(LAC_spatial_obs, LAC_size) %>%
  drop_na() %>%
  summarise(
    Mean_spatial = round(mean(LAC_spatial_obs), 2),
    SD_spatial   = round(sd(LAC_spatial_obs),   2),
    Mean_size    = round(mean(LAC_size),         2),
    SD_size      = round(sd(LAC_size),           2)
  ) %>%
  {tibble::tibble(
    Component = c("LAC_spatial_obs", "LAC_size"),
    Mean      = c(.$Mean_spatial, .$Mean_size),
    SD        = c(.$SD_spatial,   .$SD_size)
  )}

print(table_s1)


# ------------------------------------------------------------
# 6. Save all outputs
# ------------------------------------------------------------

# Fig. S2
ggsave(file.path(save_dir, "FigS2_Spatial_Gini_Comparison.pdf"),
       plot_gini, width = 89, height = 100,
       units = "mm", dpi = 600, device = cairo_pdf)

# Fig. S3
ggsave(file.path(save_dir, "FigS3_Gini_Correlation.pdf"),
       p_scatter, width = 89, height = 89,
       units = "mm", dpi = 600, device = cairo_pdf)

# Table S1
writexl::write_xlsx(
  list("Table_S1_LAC" = table_s1),
  path = file.path(save_dir, "TableS1_LAC_Summary.xlsx")
)

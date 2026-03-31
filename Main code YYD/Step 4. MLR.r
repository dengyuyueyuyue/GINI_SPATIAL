# ============================================================
# Supplementary analysis:
# multiple regression of NPP against interpolated diversity metrics
# ============================================================
# This script evaluates the relationships between NPP and three
# interpolated diversity predictors at two spatial scales
# (100 ha and 25 ha).
#
# Response variable:
#   NPP
#
# Predictors:
#   SR_interpolated
#   delta_Gini_spatial_interpolated
#   Gini_size_interpolated
#
# Statistical framework:
# 1) standardized multiple linear regression
# 2) LMG variance decomposition of model R²
#
# Main outputs:
#   Fig. 4  - Forest plot of standardized coefficients and
#             bar plot of relative contributions
#   NPP_coefficients.csv
#   NPP_forest.pdf
#   NPP_variance.pdf
#   Fig. 4 Combined_NPP.pdf
#
# Key variables:
#   SR_interpolated                  = interpolated species richness
#   delta_Gini_spatial_interpolated  = interpolated spatial structural diversity
#   Gini_size_interpolated           = interpolated size structural diversity
#
# Note:
# All predictors and the response variable are z-score standardized
# before model fitting so that regression coefficients are directly
# comparable across variables and spatial scales.
# ============================================================

rm(list = ls())

# ---------------------------
# 0. Packages
# ---------------------------
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(readxl)
library(gridExtra)
library(relaimpo)

# ---------------------------
# 1. Standardized regression
# ---------------------------
# Fit a multiple linear regression after z-score standardization
# of the response and all predictors.
#
# Output:
# - fitted model
# - tidy coefficient table with confidence intervals
# - model summary statistics
# - filtered analysis dataset
run_mlr <- function(data, response_var, predictors) {

  analysis_vars <- c(response_var, predictors)

  data_model <- data %>%
    dplyr::select(all_of(analysis_vars)) %>%
    tidyr::drop_na()

  data_model[[paste0(response_var, "_norm")]] <- as.numeric(scale(data_model[[response_var]]))

  predictors_norm <- lapply(predictors, function(x) {
    as.numeric(scale(data_model[[x]]))
  }) %>%
    as.data.frame()

  colnames(predictors_norm) <- paste0(predictors, "_norm")
  data_model <- bind_cols(data_model, predictors_norm)

  formula <- as.formula(
    paste0(
      response_var, "_norm ~ ",
      paste(paste0(predictors, "_norm"), collapse = " + ")
    )
  )

  model <- lm(formula, data = data_model)

  return(list(
    model  = model,
    tidy   = broom::tidy(model, conf.int = TRUE),
    glance = broom::glance(model),
    data   = data_model
  ))
}

# ---------------------------
# 2. Forest plot
# ---------------------------
# Visualize standardized regression coefficients and 95% confidence
# intervals across the two spatial scales.
#
# Solid points indicate p < 0.05 and open points indicate p >= 0.05.
plot_forest <- function(res_100, res_25) {

  extract_coefs <- function(res, scale_label) {
    res$tidy %>%
      filter(term != "(Intercept)") %>%
      mutate(
        Scale = scale_label,
        Label = case_when(
          grepl("delta_Gini_spatial", term) ~ "ΔGini[spatial]",
          grepl("Gini_size", term)          ~ "Gini[size]",
          grepl("SR_interpolated", term)    ~ "SR",
          TRUE                              ~ term
        )
      )
  }

  combined <- bind_rows(
    extract_coefs(res_100, "100 ha"),
    extract_coefs(res_25,  "25 ha")
  )

  combined$Label <- factor(
    combined$Label,
    levels = rev(c("SR", "ΔGini[spatial]", "Gini[size]"))
  )

  combined$Scale <- factor(combined$Scale, levels = c("100 ha", "25 ha"))

  # Slight vertical offsets avoid overlap between scales
  combined <- combined %>%
    mutate(
      y_offset = case_when(
        Scale == "100 ha" ~  0.15,
        Scale == "25 ha"  ~ -0.15
      )
    )

  subtitle_text <- sprintf(
    "100 ha: R² = %.3f, Adj.R² = %.3f | 25 ha: R² = %.3f, Adj.R² = %.3f",
    res_100$glance$r.squared, res_100$glance$adj.r.squared,
    res_25$glance$r.squared,  res_25$glance$adj.r.squared
  )

  ggplot(combined,
         aes(x = estimate,
             y = as.numeric(Label) + y_offset,
             color = Scale)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray60", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.12, linewidth = 0.8) +
    geom_point(aes(shape = p.value < 0.05), size = 5) +
    scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 16),
      labels = c("FALSE" = "P ≥ 0.05", "TRUE" = "P < 0.05"),
      name = "Significance"
    ) +
    scale_color_manual(
      values = c("100 ha" = "#66C2A5", "25 ha" = "#FC8D62"),
      name = "Scale"
    ) +
    scale_y_continuous(
      breaks = 1:3,
      labels = parse(text = c("Gini[size]", "Delta*Gini[spatial]", "SR")),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    labs(
      x = "Standardized coefficient",
      y = NULL,
      subtitle = subtitle_text
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(size = 13, face = "bold"),
      legend.text      = element_text(size = 12),
      axis.title       = element_text(size = 14, face = "bold"),
      axis.text.y      = element_text(size = 13),
      axis.text.x      = element_text(size = 12),
      plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "gray30"),
      panel.grid.minor = element_blank()
    )
}

# ---------------------------
# 3. LMG variance decomposition
# ---------------------------
# Quantify the relative contribution of each predictor to the
# explained variance (R²) using the LMG method.
variance_decomposition <- function(results, scale_name) {
  model <- results$model
  imp   <- relaimpo::calc.relimp(model, type = "lmg", rela = TRUE)
  r2    <- summary(model)$r.squared
  adjr2 <- summary(model)$adj.r.squared

  data.frame(
    Variable   = names(imp$lmg),
    Importance = 100 * imp$lmg,
    Scale      = sprintf("%s (R² = %.3f, Adj.R² = %.3f)", scale_name, r2, adjr2),
    R2         = r2,
    Adj_R2     = adjr2
  ) %>%
    mutate(
      Label = case_when(
        grepl("delta_Gini_spatial", Variable) ~ "ΔGini[spatial]",
        grepl("Gini_size", Variable)          ~ "Gini[size]",
        grepl("SR_interpolated", Variable)    ~ "SR",
        TRUE                                  ~ Variable
      )
    )
}

# Visualize the relative importance of predictors at the two scales.
plot_variance <- function(imp_100, imp_25) {
  combined <- bind_rows(imp_100, imp_25)

  combined$Scale <- factor(combined$Scale, levels = unique(combined$Scale))
  combined$Label <- factor(
    combined$Label,
    levels = c("SR", "ΔGini[spatial]", "Gini[size]")
  )

  ggplot(combined, aes(x = Label, y = Importance, fill = Label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Importance)),
              vjust = -0.3, size = 3) +
    facet_wrap(~ Scale, ncol = 2) +
    scale_fill_manual(values = c(
      "SR"             = "#66C2A5",
      "ΔGini[spatial]" = "#4DBBD5",
      "Gini[size]"     = "#E64B35"
    )) +
    scale_x_discrete(labels = parse(text = c("SR", "Delta*Gini[spatial]", "Gini[size]"))) +
    labs(
      x = NULL,
      y = "Relative contribution to explained variance (%)"
    ) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
      strip.text      = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

# ---------------------------
# 4. Export coefficient table
# ---------------------------
# Save standardized regression coefficients and model-level R²
# statistics for both spatial scales.
save_coef_table <- function(res_100, res_25, response_var, output_folder) {
  bind_rows(
    res_100$tidy %>%
      mutate(
        Scale  = "100 ha",
        R2     = res_100$glance$r.squared,
        Adj_R2 = res_100$glance$adj.r.squared
      ),
    res_25$tidy %>%
      mutate(
        Scale  = "25 ha",
        R2     = res_25$glance$r.squared,
        Adj_R2 = res_25$glance$adj.r.squared
      )
  ) %>%
    filter(term != "(Intercept)") %>%
    dplyr::select(Scale, term, estimate, std.error, p.value,
                  conf.low, conf.high, R2, Adj_R2) %>%
    write.csv(
      file.path(output_folder, paste0(response_var, "_coefficients.csv")),
      row.names = FALSE
    )
}

# ---------------------------
# 5. Main analysis
# ---------------------------
# Read grid-level zonal statistics at the 100-ha and 25-ha scales,
# fit standardized regressions, export Fig. 4, and save coefficient tables.

data_path     <- "zonal_stats_all_scales.xlsx"
output_folder <- "zonal_statistics"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

predictors <- c(
  "SR_interpolated",
  "delta_Gini_spatial_interpolated",
  "Gini_size_interpolated"
)

cat("\n=== Analyzing NPP ===\n")

data_100 <- readxl::read_xlsx(data_path, sheet = "1000m")
data_25  <- readxl::read_xlsx(data_path, sheet = "500m")

npp_100 <- run_mlr(data_100, "NPP", predictors)
npp_25  <- run_mlr(data_25,  "NPP", predictors)

# ---------------------------
# 6. Model summary
# ---------------------------
cat("\n--- Model fit summary ---\n")
cat(sprintf(
  "NPP 100 ha: R² = %.4f, Adj.R² = %.4f, F-test p = %.2e\n",
  npp_100$glance$r.squared,
  npp_100$glance$adj.r.squared,
  npp_100$glance$p.value
))
cat(sprintf(
  "NPP  25 ha: R² = %.4f, Adj.R² = %.4f, F-test p = %.2e\n",
  npp_25$glance$r.squared,
  npp_25$glance$adj.r.squared,
  npp_25$glance$p.value
))

# ---------------------------
# 7. Fig. 4
# ---------------------------
# Left panel: standardized regression coefficients
# Right panel: LMG relative importance

p_npp_forest <- plot_forest(npp_100, npp_25)
ggsave(
  file.path(output_folder, "NPP_forest.pdf"),
  p_npp_forest,
  width = 8, height = 6
)

imp_npp_100 <- variance_decomposition(npp_100, "100 ha")
imp_npp_25  <- variance_decomposition(npp_25,  "25 ha")

p_npp_var <- plot_variance(imp_npp_100, imp_npp_25)
ggsave(
  file.path(output_folder, "NPP_variance.pdf"),
  p_npp_var,
  width = 8, height = 5
)

combined_npp <- gridExtra::arrangeGrob(
  p_npp_forest, p_npp_var, ncol = 2
)

ggsave(
  file.path(output_folder, "Fig. 4 Combined_NPP.pdf"),
  combined_npp,
  width = 14, height = 6
)

# ---------------------------
# 8. Export coefficient table
# ---------------------------
save_coef_table(npp_100, npp_25, "NPP", output_folder)

# ---------------------------
# 9. Console output
# ---------------------------
cat("\n--- 100 ha model summary ---\n")
print(summary(npp_100$model))

cat("\n--- 25 ha model summary ---\n")
print(summary(npp_25$model))

cat("\nNPP analysis completed.\n")
cat(sprintf("Results saved to: %s\n", output_folder))

# ============================================================
# Supplementary analysis:
# piecewise structural equation modelling of NPP
# ============================================================
# This script fits piecewise structural equation models (pSEM)
# at the 100-ha and 25-ha scales to evaluate how SR,
# spatial structural diversity, and size structural diversity
# jointly influence NPP.
#
# Analytical steps:
# 1) read grid-level data from zonal statistics;
# 2) fit candidate pSEM models at each spatial scale;
# 3) compare candidate models using AIC;
# 4) extract standardized path coefficients from the
#    best-supported model at each scale;
# 5) calculate direct, indirect, and total effects on NPP;
# 6) export effect plots and summary tables.
#
# Inputs:
#   zonal_stats_all_scales.xlsx
#
# Main outputs:
#   NPP_path_coefficients_all_scales.pdf
#   NPP_path_coefficients_all_scales.csv
#
# Key variables:
#   SR                  = species richness
#   delta_Gini_spatial  = spatial structural diversity
#   Gini_size           = size structural diversity
#   NPP                 = net primary productivity
# ============================================================

# ---------------------------
# 0. Packages
# ---------------------------
library(piecewiseSEM)
library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)

# ---------------------------
# 1. Read and prepare data
# ---------------------------
# Read grid-level datasets for the two spatial scales and
# rename variables to a common format for pSEM fitting.

file_path <- "zonal_stats_all_scales.xlsx"

data_100ha_raw <- read_excel(file_path, sheet = "1000m")
data_25ha_raw  <- read_excel(file_path, sheet = "500m")

prepare_psem_data <- function(dat) {
  data.frame(
    site               = seq_len(nrow(dat)),
    SR                 = dat$SR_interpolated,
    delta_Gini_spatial = dat$delta_Gini_spatial_interpolated,
    NPP                = dat$NPP,
    Gini_size          = dat$Gini_size_interpolated
  ) %>%
    na.omit()
}

data_100ha_mat <- prepare_psem_data(data_100ha_raw)
data_25ha_mat  <- prepare_psem_data(data_25ha_raw)

# ---------------------------
# 2. Candidate pSEM models
# ---------------------------
# Model A includes a direct path from SR to NPP.
# Model B removes the direct SR -> NPP path and represents
# the best-supported model used for effect decomposition.

# 100 ha
psem_100ha_H1_A <- psem(
  lm(NPP ~ delta_Gini_spatial + Gini_size + SR, data = data_100ha_mat),
  lm(delta_Gini_spatial ~ SR,                   data = data_100ha_mat),
  lm(Gini_size ~ SR,                            data = data_100ha_mat),
  delta_Gini_spatial %~~% Gini_size
)
summary(psem_100ha_H1_A)

psem_100ha_H1_B <- psem(
  lm(NPP ~ delta_Gini_spatial + Gini_size, data = data_100ha_mat),
  lm(delta_Gini_spatial ~ SR,              data = data_100ha_mat),
  lm(Gini_size ~ SR,                       data = data_100ha_mat)
)
summary(psem_100ha_H1_B)
plot(psem_100ha_H1_B)

AIC(psem_100ha_H1_A, psem_100ha_H1_B)

# 25 ha
psem_25ha_H1_A <- psem(
  lm(NPP ~ delta_Gini_spatial + Gini_size + SR, data = data_25ha_mat),
  lm(delta_Gini_spatial ~ SR,                   data = data_25ha_mat),
  lm(Gini_size ~ SR,                            data = data_25ha_mat),
  delta_Gini_spatial %~~% Gini_size
)
summary(psem_25ha_H1_A)

psem_25ha_H1_B <- psem(
  lm(NPP ~ delta_Gini_spatial + Gini_size, data = data_25ha_mat),
  lm(delta_Gini_spatial ~ SR,              data = data_25ha_mat),
  lm(Gini_size ~ SR,                       data = data_25ha_mat),
  delta_Gini_spatial %~~% Gini_size
)
summary(psem_25ha_H1_B)
plot(psem_25ha_H1_B)

AIC(psem_25ha_H1_A, psem_25ha_H1_B)

# ---------------------------
# 3. Extract standardized path coefficients
# ---------------------------
# Extract standardized coefficients and p values from a pSEM summary.

extract_sem_coefs <- function(psem_model) {
  sem_summary <- summary(psem_model)
  coef_table  <- sem_summary$coefficients

  data.frame(
    Response     = coef_table[, "Response"],
    Predictor    = coef_table[, "Predictor"],
    Std.Estimate = coef_table[, "Std.Estimate"],
    P.Value      = coef_table[, "P.Value"],
    stringsAsFactors = FALSE
  )
}

# ---------------------------
# 4. Effect decomposition
# ---------------------------
# Calculate direct, indirect, and total effects of each predictor
# on NPP based on standardized path coefficients.

calculate_effects_to_NPP <- function(coefs_df, all_vars = NULL) {
  effects_list    <- list()
  target_response <- "NPP"

  if (is.null(all_vars)) {
    all_predictors <- unique(coefs_df$Predictor[coefs_df$Predictor != target_response])
    all_predictors <- all_predictors[!grepl("~~", all_predictors)]
  } else {
    all_predictors <- all_vars
  }

  direct_to_NPP <- coefs_df[
    coefs_df$Response == target_response &
      !grepl("~~", coefs_df$Predictor), ]

  cat("\nDirect paths to NPP:\n")
  print(direct_to_NPP[, c("Predictor", "Std.Estimate", "P.Value")])

  for (predictor in all_predictors) {
    cat("\n", rep("-", 70), "\n", sep = "")
    cat("Predictor:", predictor, "\n")
    cat(rep("-", 70), "\n", sep = "")

    # Direct effect
    direct_effect <- 0
    direct_row <- direct_to_NPP[direct_to_NPP$Predictor == predictor, ]

    if (nrow(direct_row) > 0) {
      direct_effect <- direct_row$Std.Estimate
      cat(sprintf("Direct effect (%s -> NPP): %.4f (P = %.4f)\n",
                  predictor, direct_effect, direct_row$P.Value))
    } else {
      cat(sprintf("No direct path detected (%s -> NPP)\n", predictor))
    }

    # Indirect effect
    indirect_effect <- 0
    indirect_paths  <- list()
    mediators <- direct_to_NPP$Predictor[direct_to_NPP$Predictor != predictor]

    cat("\nIndirect paths:\n")
    for (mediator in mediators) {
      path_to_mediator <- coefs_df[
        coefs_df$Response == mediator &
          coefs_df$Predictor == predictor &
          !grepl("~~", coefs_df$Predictor), ]

      if (nrow(path_to_mediator) > 0) {
        coef1 <- path_to_mediator$Std.Estimate
        coef2_row <- direct_to_NPP[direct_to_NPP$Predictor == mediator, ]

        if (nrow(coef2_row) > 0) {
          coef2 <- coef2_row$Std.Estimate
          indirect_contrib <- coef1 * coef2
          indirect_effect <- indirect_effect + indirect_contrib

          path_desc <- sprintf(
            "  %s -> %s -> NPP: %.4f × %.4f = %.4f",
            predictor, mediator, coef1, coef2, indirect_contrib
          )
          cat(path_desc, "\n")
          indirect_paths[[length(indirect_paths) + 1]] <- path_desc
        }
      }
    }

    if (length(indirect_paths) == 0) cat("  No indirect path detected\n")

    total_effect <- direct_effect + indirect_effect
    cat(sprintf("\nTotal indirect effect: %.4f\n", indirect_effect))
    cat(sprintf("Total effect: %.4f (direct: %.4f + indirect: %.4f)\n",
                total_effect, direct_effect, indirect_effect))

    effects_list[[predictor]] <- list(
      direct   = direct_effect,
      indirect = indirect_effect,
      total    = total_effect
    )
  }

  return(effects_list)
}

# Convert the effect list to a wide summary table.
prepare_effects_df <- function(effects_list, var_order = NULL) {
  effects_df <- data.frame(
    predictor       = names(effects_list),
    direct_effect   = sapply(effects_list, `[[`, "direct"),
    indirect_effect = sapply(effects_list, `[[`, "indirect"),
    total_effect    = sapply(effects_list, `[[`, "total"),
    stringsAsFactors = FALSE
  )

  if (!is.null(var_order)) {
    effects_df$predictor <- factor(effects_df$predictor, levels = var_order)
    effects_df <- effects_df[order(effects_df$predictor), ]
    effects_df$predictor <- as.character(effects_df$predictor)
  }

  return(effects_df)
}

# Convert the wide effect table to long format for plotting.
effects_to_long_df <- function(effects_df, scale_name) {
  data.frame(
    scale       = scale_name,
    predictor   = rep(effects_df$predictor, 3),
    effect_type = rep(c("direct", "indirect", "total"), each = nrow(effects_df)),
    value       = c(effects_df$direct_effect,
                    effects_df$indirect_effect,
                    effects_df$total_effect),
    stringsAsFactors = FALSE
  )
}

# ---------------------------
# 5. Run effect decomposition
# ---------------------------
# Effect decomposition is performed using Model B at both scales.

all_variables <- c("SR", "delta_Gini_spatial", "Gini_size")

cat("\n", rep("=", 80), "\n", sep = "")
cat("100 ha effect decomposition\n")
cat(rep("=", 80), "\n", sep = "")
coefs_100ha      <- extract_sem_coefs(psem_100ha_H1_B)
effects_100ha    <- calculate_effects_to_NPP(coefs_100ha, all_vars = all_variables)
effects_df_100ha <- prepare_effects_df(effects_100ha, var_order = all_variables)

cat("\n", rep("=", 80), "\n", sep = "")
cat("25 ha effect decomposition\n")
cat(rep("=", 80), "\n", sep = "")
coefs_25ha      <- extract_sem_coefs(psem_25ha_H1_B)
effects_25ha    <- calculate_effects_to_NPP(coefs_25ha, all_vars = all_variables)
effects_df_25ha <- prepare_effects_df(effects_25ha, var_order = all_variables)

cat("\n", rep("=", 80), "\n", sep = "")
cat("Summary of decomposed effects\n")
cat(rep("=", 80), "\n\n", sep = "")
print("100 ha effects:"); print(effects_df_100ha)
print("25 ha effects:"); print(effects_df_25ha)

# ---------------------------
# 6. Plot direct, indirect, and total effects
# ---------------------------
# Create one panel for each spatial scale and compare the
# three effect components across predictors.

predictor_labels <- c(
  "SR"                 = expression(SR),
  "delta_Gini_spatial" = expression(SD[spatial]),
  "Gini_size"          = expression(SD[size])
)

var_levels <- c("SR", "delta_Gini_spatial", "Gini_size")

df_100ha <- effects_to_long_df(effects_df_100ha, "100 ha")
df_25ha  <- effects_to_long_df(effects_df_25ha,  "25 ha")

for (df_name in c("df_100ha", "df_25ha")) {
  df <- get(df_name)
  df$predictor   <- factor(df$predictor,   levels = var_levels)
  df$effect_type <- factor(df$effect_type, levels = c("direct", "indirect", "total"))
  assign(df_name, df)
}

plot_path_coefficients <- function(scale_data, scale_name) {
  scale_data$predictor <- factor(scale_data$predictor, levels = var_levels)

  colors   <- c("direct" = "#F9A7A7", "indirect" = "#83B5DD", "total" = "#8DC9A4")
  y_max    <- max(abs(scale_data$value), na.rm = TRUE) * 1.3
  y_limits <- c(-y_max, y_max)

  ggplot(scale_data, aes(x = predictor, y = value, fill = effect_type)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", value)),
              position = position_dodge(width = 0.8),
              vjust = ifelse(scale_data$value >= 0, -0.5, 1.5),
              size = 3.5, fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "solid",
               color = "black", linewidth = 0.6) +
    scale_fill_manual(values = colors,
                      labels = c("Direct", "Indirect", "Total")) +
    scale_x_discrete(labels = predictor_labels) +
    labs(title = scale_name,
         x     = NULL,
         y     = "Standardized effect on NPP",
         fill  = "Effect type") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title.y       = element_text(size = 13, face = "bold"),
      axis.text          = element_text(size = 11, color = "black"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5),
      legend.position    = "top",
      legend.text        = element_text(size = 11),
      legend.title       = element_text(size = 12, face = "bold"),
      plot.title         = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_y_continuous(limits = y_limits,
                       breaks = seq(-1, 1, 0.2),
                       expand = expansion(mult = c(0.02, 0.02)))
}

plot_100ha <- plot_path_coefficients(df_100ha, "100 ha")
plot_25ha  <- plot_path_coefficients(df_25ha,  "25 ha")

combined_plot <- plot_100ha / plot_25ha +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(legend.position = "bottom")
  )

print(combined_plot)

# ---------------------------
# 7. Export outputs
# ---------------------------
output_dir <- "pSEM"
if (!dir.exists(output_dir)) dir.create(output_dir)

ggsave(file.path(output_dir, "NPP_path_coefficients_all_scales.pdf"),
       combined_plot, width = 6, height = 10, dpi = 300)

effects_df_100ha$scale <- "100 ha"
effects_df_25ha$scale  <- "25 ha"

all_effects_table <- rbind(effects_df_100ha, effects_df_25ha)
all_effects_table <- all_effects_table[, c("scale", "predictor",
                                           "direct_effect", "indirect_effect", "total_effect")]

write.csv(all_effects_table,
          file.path(output_dir, "NPP_path_coefficients_all_scales.csv"),
          row.names = FALSE)

cat("\n", rep("=", 80), "\n", sep = "")
cat("Analysis completed. All outputs were saved to:", output_dir, "\n")
cat(rep("=", 80), "\n", sep = "")
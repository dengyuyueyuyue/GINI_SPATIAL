library(piecewiseSEM)
library(lme4)
library(fastDummies)
library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)

# ============================================
# 读取数据
# ============================================
file_path <- "zonal_stats_all_scales.xlsx"

data_100ha_raw <- read_excel(file_path, sheet = "1000m")
data_25ha_raw  <- read_excel(file_path, sheet = "500m")

# 100 ha 数据整理
data_100ha_mat <- data.frame(
  site               = 1:nrow(data_100ha_raw),
  SR                 = data_100ha_raw$SR_interpolated,
  delta_Gini_spatial = data_100ha_raw$delta_Gini_spatial_interpolated,
  NPP                = data_100ha_raw$NPP,
  Gini_size          = data_100ha_raw$Gini_size_interpolated
)

# 25 ha 数据整理（变量名与 100ha 统一）
data_25ha_mat <- data.frame(
  site               = 1:nrow(data_25ha_raw),
  SR                 = data_25ha_raw$SR_interpolated,
  delta_Gini_spatial = data_25ha_raw$delta_Gini_spatial_interpolated,
  NPP                = data_25ha_raw$NPP,
  Gini_size          = data_25ha_raw$Gini_size_interpolated
)

head(data_100ha_mat)
head(data_25ha_mat)
colnames(data_100ha_mat)
colnames(data_25ha_mat)

# ============================================
# 100 ha pSEM 模型
# ============================================
psem_100ha_H1_A <- psem(
  lm(NPP ~ delta_Gini_spatial + Gini_size + SR, data = data_100ha_mat),
  lm(delta_Gini_spatial ~ SR,                   data = data_100ha_mat),
  lm(Gini_size ~ SR,                            data = data_100ha_mat),
  delta_Gini_spatial %~~% Gini_size
)
summary(psem_100ha_H1_A)

#Best model for 100ha
psem_100ha_H1_B <- psem(
  lm(NPP ~ delta_Gini_spatial + Gini_size, data = data_100ha_mat),
  lm(delta_Gini_spatial ~ SR,              data = data_100ha_mat),
  lm(Gini_size ~ SR,                       data = data_100ha_mat)
)
summary(psem_100ha_H1_B)
plot(psem_100ha_H1_B)

AIC(psem_100ha_H1_A, psem_100ha_H1_B)

# ============================================
# 25 ha pSEM 模型
# ============================================
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

# ============================================
# 效应分解函数
# ============================================

# 函数1：提取标准化系数
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

# 函数2：计算直接、间接、总效应
calculate_effects_to_NPP <- function(coefs_df, all_vars = NULL) {
  effects_list    <- list()
  target_response <- "NPP"

  if (is.null(all_vars)) {
    all_predictors <- unique(coefs_df$Predictor[coefs_df$Predictor != target_response])
    all_predictors <- all_predictors[!grepl("~~", all_predictors)]
  } else {
    all_predictors <- all_vars
  }

  direct_to_NPP <- coefs_df[coefs_df$Response == target_response &
                               !grepl("~~", coefs_df$Predictor), ]

  cat("\n直接指向 NPP 的路径:\n")
  print(direct_to_NPP[, c("Predictor", "Std.Estimate", "P.Value")])

  for (predictor in all_predictors) {
    cat("\n", rep("-", 70), "\n", sep = "")
    cat("分析变量:", predictor, "\n")
    cat(rep("-", 70), "\n", sep = "")

    # 直接效应
    direct_effect <- 0
    direct_row    <- direct_to_NPP[direct_to_NPP$Predictor == predictor, ]
    if (nrow(direct_row) > 0) {
      direct_effect <- direct_row$Std.Estimate
      cat(sprintf("✓ 直接效应 (%s -> NPP): %.4f (P = %.4f)\n",
                  predictor, direct_effect, direct_row$P.Value))
    } else {
      cat(sprintf("✗ 没有直接路径 (%s -> NPP)\n", predictor))
    }

    # 间接效应
    indirect_effect <- 0
    indirect_paths  <- list()
    mediators       <- direct_to_NPP$Predictor[direct_to_NPP$Predictor != predictor]

    cat("\n间接路径分析:\n")
    for (mediator in mediators) {
      path_to_mediator <- coefs_df[coefs_df$Response == mediator &
                                     coefs_df$Predictor == predictor &
                                     !grepl("~~", coefs_df$Predictor), ]
      if (nrow(path_to_mediator) > 0) {
        coef1     <- path_to_mediator$Std.Estimate
        coef2_row <- direct_to_NPP[direct_to_NPP$Predictor == mediator, ]
        if (nrow(coef2_row) > 0) {
          coef2            <- coef2_row$Std.Estimate
          indirect_contrib <- coef1 * coef2
          indirect_effect  <- indirect_effect + indirect_contrib
          path_desc <- sprintf("  %s -> %s -> NPP: %.4f × %.4f = %.4f",
                               predictor, mediator, coef1, coef2, indirect_contrib)
          cat(path_desc, "\n")
          indirect_paths[[length(indirect_paths) + 1]] <- path_desc
        }
      }
    }
    if (length(indirect_paths) == 0) cat("  ✗ 没有间接路径\n")

    total_effect <- direct_effect + indirect_effect
    cat(sprintf("\n总间接效应: %.4f\n", indirect_effect))
    cat(sprintf("总效应: %.4f (直接: %.4f + 间接: %.4f)\n",
                total_effect, direct_effect, indirect_effect))

    effects_list[[predictor]] <- list(
      direct   = direct_effect,
      indirect = indirect_effect,
      total    = total_effect
    )
  }
  return(effects_list)
}

# 函数3：整理为宽格式数据框
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

# 函数4：转换为长格式
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

# ============================================
# 执行效应分解（均使用模型 B）
# ============================================
all_variables <- c("SR", "delta_Gini_spatial", "Gini_size")

cat("\n", rep("=", 80), "\n", sep = "")
cat("100 ha 模型效应分解\n")
cat(rep("=", 80), "\n", sep = "")
coefs_100ha      <- extract_sem_coefs(psem_100ha_H1_B)
effects_100ha    <- calculate_effects_to_NPP(coefs_100ha, all_vars = all_variables)
effects_df_100ha <- prepare_effects_df(effects_100ha, var_order = all_variables)

cat("\n", rep("=", 80), "\n", sep = "")
cat("25 ha 模型效应分解\n")
cat(rep("=", 80), "\n", sep = "")
coefs_25ha      <- extract_sem_coefs(psem_25ha_H1_B)
effects_25ha    <- calculate_effects_to_NPP(coefs_25ha, all_vars = all_variables)
effects_df_25ha <- prepare_effects_df(effects_25ha, var_order = all_variables)

cat("\n", rep("=", 80), "\n", sep = "")
cat("效应分解结果摘要\n")
cat(rep("=", 80), "\n\n", sep = "")
print("100 ha 效应值:"); print(effects_df_100ha)
print("25 ha 效应值:");  print(effects_df_25ha)

# ============================================
# 绘图
# ============================================

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

# 绘图函数
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

# ============================================
# 保存结果到 pSEM 文件夹
# ============================================
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
cat("分析完成！所有结果已保存到:", output_dir, "\n")
cat(rep("=", 80), "\n", sep = "")

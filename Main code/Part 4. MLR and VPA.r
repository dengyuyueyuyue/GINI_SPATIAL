###############################################
# 两个空间尺度（100ha / 25ha）的 NPP–多样性回归分析
# 响应变量  : NPP
# 预测变量  : SR_interpolated（物种丰富度插值）
#             delta_Gini_spatial_interpolated（空间 Gini 差异插值）
#             Gini_size_interpolated（个体大小 Gini 插值）
# 统计方法  : 标准化多元线性回归 + LMG 方差分解
# 对应图件  : Fig. 4（森林图 + 方差分解组合图）
###############################################

rm(list = ls())

# 加载必要的库
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(gridExtra)
library(readxl)
library(bestNormalize)
library(grid)
library(relaimpo)

#-----------------------
# 第 1 节：数据准备和建模
# 对响应变量与预测变量均进行 z-score 标准化，
# 使回归系数可直接跨变量比较（标准化系数）
#-----------------------
run_mlr <- function(data, response_var, predictors) {

  # 标准化响应变量
  data[[paste0(response_var, "_norm")]] <- scale(data[[response_var]])

  # 标准化所有预测变量
  predictors_norm <- sapply(predictors, function(x) scale(data[[x]]))
  predictors_norm <- as.data.frame(predictors_norm)
  colnames(predictors_norm) <- paste0(predictors, "_norm")
  data <- cbind(data, predictors_norm)

  # 构建标准化回归公式
  formula <- as.formula(paste0(
    response_var, "_norm", " ~ ",
    paste(paste0(predictors, "_norm"), collapse = " + ")
  ))
  model <- lm(formula, data = data)

  return(list(
    model  = model,
    tidy   = tidy(model, conf.int = TRUE),
    glance = glance(model),
    data   = data
  ))
}

#-----------------------
# 第 2 节：森林图（跨尺度比较）
# 展示两个空间尺度下各预测变量的标准化回归系数及 95% CI；
# 实心点表示 P < 0.05，副标题标注各尺度 R² / Adj.R²
#-----------------------
plot_forest <- function(res_100, res_25) {

  extract_coefs <- function(res, scale) {
    res$tidy %>%
      filter(term != "(Intercept)") %>%
      mutate(
        scale = scale,
        term_clean = case_when(
          grepl("delta_Gini_spatial", term) ~ "ΔGini[spatial]",
          grepl("Gini_size",          term) ~ "Gini[size]",
          grepl("SR_interpolated",    term) ~ "SR",
          TRUE ~ term
        )
      )
  }

  combined <- bind_rows(
    extract_coefs(res_100, "100ha"),
    extract_coefs(res_25,  "25ha")
  )

  # Y 轴顺序：SR > ΔGini_spatial > Gini_size
  combined$term_clean <- factor(combined$term_clean,
                                 levels = rev(c("SR", "ΔGini[spatial]", "Gini[size]")))
  combined$scale <- factor(combined$scale, levels = c("100ha", "25ha"))

  # 同一变量两个尺度点位上下轻微错开，避免重叠
  combined <- combined %>%
    mutate(y_offset = case_when(
      scale == "100ha" ~  0.15,
      scale == "25ha"  ~ -0.15
    ))

  # 构建含 R² 和 Adj.R² 的副标题
  r2_100    <- res_100$glance$r.squared
  adjr2_100 <- res_100$glance$adj.r.squared
  r2_25     <- res_25$glance$r.squared
  adjr2_25  <- res_25$glance$adj.r.squared

  subtitle_text <- sprintf(
    "100ha: R² = %.3f, Adj.R² = %.3f | 25ha: R² = %.3f, Adj.R² = %.3f",
    r2_100, adjr2_100, r2_25, adjr2_25
  )

  ggplot(combined, aes(x = estimate, y = as.numeric(term_clean) + y_offset, color = scale)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", size = 0.8) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.12, size = 0.8) +
    geom_point(aes(shape = p.value < 0.05), size = 5) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                       labels = c("FALSE" = "P ≥ 0.05", "TRUE" = "P < 0.05"),
                       name = "Significance") +
    scale_color_manual(values = c("100ha" = "#66C2A5", "25ha" = "#FC8D62"),
                       name = "Scale") +
    scale_y_continuous(
      breaks = 1:3,
      labels = parse(text = c("Gini[size]", "Delta*Gini[spatial]", "SR")),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    labs(
      x        = "Standardized Coefficient",
      y        = "",
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

#-----------------------
# 第 3 节：方差分解（LMG 相对重要性）
# 使用 relaimpo::calc.relimp 计算各预测变量对
# 模型 R² 的相对贡献（LMG 法，结果归一化为 100%）
#-----------------------
variance_decomposition <- function(results, scale_name) {
  model <- results$model
  imp   <- calc.relimp(model, type = "lmg", rela = TRUE)
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
        grepl("Gini_size",          Variable) ~ "Gini[size]",
        grepl("SR_interpolated",    Variable) ~ "SR",
        TRUE ~ Variable
      )
    )
}

plot_variance <- function(imp_100, imp_25) {
  combined <- bind_rows(imp_100, imp_25)

  # facet 标签保留含 R² 的字符串，按 100ha → 25ha 排序
  scale_levels   <- unique(combined$Scale)
  combined$Scale <- factor(combined$Scale, levels = scale_levels)
  combined$Label <- factor(combined$Label,
                           levels = c("SR", "ΔGini[spatial]", "Gini[size]"))

  ggplot(combined, aes(x = Label, y = Importance, fill = Label)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Importance)), vjust = -0.3, size = 3) +
    facet_wrap(~ Scale, ncol = 2) +
    scale_fill_manual(values = c(
      "SR"             = "#66C2A5",
      "ΔGini[spatial]" = "#4DBBD5",
      "Gini[size]"     = "#E64B35"
    )) +
    scale_x_discrete(labels = parse(text = c("SR", "Delta*Gini[spatial]", "Gini[size]"))) +
    labs(x = NULL, y = "Relative contribution to explained variance (%)") +
    theme_bw() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
      strip.text      = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )
}

#-----------------------
# 第 4 节：主程序 — NPP 分析
# 读取两个空间尺度的区域统计汇总表，分别拟合标准化多元线性回归，
# 输出森林图、方差分解图及系数 CSV
#-----------------------

# 输入 / 输出路径
data_path     <- "zonal_stats_all_scales.xlsx"
output_folder <- "zonal_statistics"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# 预测变量（均来自插值栅格的区域统计结果）
predictors <- c("SR_interpolated", "delta_Gini_spatial_interpolated", "Gini_size_interpolated")

# ---------- NPP 分析 ----------
cat("\n=== Analyzing NPP ===\n")
data_100 <- read_xlsx(data_path, sheet = "1000m")
data_25  <- read_xlsx(data_path, sheet = "500m")

npp_100 <- run_mlr(data_100, "NPP", predictors)
npp_25  <- run_mlr(data_25,  "NPP", predictors)

# 打印 R² 汇总
cat("\n--- R² Summary ---\n")
cat(sprintf("NPP 100ha: R² = %.4f, Adj.R² = %.4f, F-stat p = %.2e\n",
            npp_100$glance$r.squared, npp_100$glance$adj.r.squared, npp_100$glance$p.value))
cat(sprintf("NPP  25ha: R² = %.4f, Adj.R² = %.4f, F-stat p = %.2e\n",
            npp_25$glance$r.squared,  npp_25$glance$adj.r.squared,  npp_25$glance$p.value))

# 森林图（含 R² 副标题）
p_npp_forest <- plot_forest(npp_100, npp_25)
ggsave(file.path(output_folder, "NPP_forest.pdf"), p_npp_forest, width = 8, height = 6)

# 方差分解（facet 标签含 R²）
imp_npp_100 <- variance_decomposition(npp_100, "100ha")
imp_npp_25  <- variance_decomposition(npp_25,  "25ha")

p_npp_var <- plot_variance(imp_npp_100, imp_npp_25)
ggsave(file.path(output_folder, "NPP_variance.pdf"), p_npp_var, width = 8, height = 5)

# 组合图（森林图 + 方差分解并排）
combined_npp <- grid.arrange(p_npp_forest, p_npp_var, ncol = 2)
ggsave(file.path(output_folder, "Fig. 4 Combined_NPP.pdf"), combined_npp, width = 14, height = 6)

# 保存系数表（含 R² 和 Adj.R² 列）
save_coef_table <- function(res_100, res_25, response_var) {
  bind_rows(
    res_100$tidy %>% mutate(Scale  = "100ha",
                             R2     = res_100$glance$r.squared,
                             Adj_R2 = res_100$glance$adj.r.squared),
    res_25$tidy  %>% mutate(Scale  = "25ha",
                             R2     = res_25$glance$r.squared,
                             Adj_R2 = res_25$glance$adj.r.squared)
  ) %>%
    filter(term != "(Intercept)") %>%
    dplyr::select(Scale, term, estimate, std.error, p.value, conf.low, conf.high, R2, Adj_R2) %>%
    write.csv(file.path(output_folder, paste0(response_var, "_coefficients.csv")), row.names = FALSE)
}

save_coef_table(npp_100, npp_25, "NPP")

# 打印模型摘要
cat("\n--- 100ha Model Summary ---\n")
print(summary(npp_100$model))
cat("\n--- 25ha Model Summary ---\n")
print(summary(npp_25$model))

cat("\n✅ NPP analysis completed!\n")
cat(sprintf("📁 Results saved to: %s\n", output_folder))






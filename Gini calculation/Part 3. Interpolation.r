# ============================================================
#  Spatial Interpolation: SR / delta_Gini_spatial / Gini_size
#  Method : LOESS + 10-Fold Cross-Validation
# ============================================================
rm(list = ls())
#------------------------
#! 1. 加载所需的包
#------------------------
library(sf)
library(sp)
library(raster)
library(caret)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(readxl)
library(writexl)
library(scales)

#------------------------
#! 2. 基础函数定义
#------------------------

# --- 2.1 数据转换与逆变换 ---

transform_data <- function(data, attr_name) {
  constraints    <- value_constraints[[attr_name]]
  transform_type <- constraints$transform

  if (transform_type != "square" && min(data, na.rm = TRUE) < 0) {
    data <- data - min(data, na.rm = TRUE) + 0.001
  }

  switch(transform_type,
    "log"   = log(data + 1),
    "logit" = log(data / (1 - data)),
    "none"  = data
  )
}

inverse_transform <- function(data, attr_name) {
  constraints    <- value_constraints[[attr_name]]
  transform_type <- constraints$transform

  result <- switch(transform_type,
    "log"   = exp(data) - 1,
    "logit" = 1 / (1 + exp(-data)),
    "none"  = data
  )

  if (!is.null(constraints$integer) && constraints$integer) {
    result <- round(result)
  }
  return(result)
}

# --- 2.2 创建预测网格 ---

create_prediction_grid <- function(boundary_sp, cellsize_deg) {
  r   <- raster(extent(boundary_sp))
  res(r) <- cellsize_deg
  grd    <- rasterToPoints(r, spatial = TRUE)
  proj4string(grd) <- proj4string(boundary_sp)
  gridded(grd) <- TRUE
  return(grd)
}

# --- 2.3 空间插值主函数（10折CV） ---

simulate_and_loss <- function(boundary, points, attrName,
                               cellsize_deg = cellsize_deg_30,
                               k_fold = 10,
                               span   = 0.2) {

  message(sprintf("\n===== 开始处理: %s =====", attrName))

  # 统一坐标系
  if (st_crs(boundary) != st_crs(points)) {
    points <- st_transform(points, st_crs(boundary))
  }

  boundary_sp <- as(boundary, "Spatial")
  points_sp   <- as(points,   "Spatial")

  # 创建预测网格
  grd    <- create_prediction_grid(boundary_sp, cellsize_deg)
  grd_df <- as.data.frame(coordinates(grd))
  names(grd_df) <- c("x", "y")

  # 数据转换
  original_values             <- points_sp@data[[attrName]]
  points_sp$transformed_value <- transform_data(original_values, attrName)

  pts_coords        <- as.data.frame(coordinates(points_sp))
  names(pts_coords) <- c("x", "y")
  pts_df            <- cbind(pts_coords, value = points_sp$transformed_value)

  # 10 折交叉验证
  message("开始10折交叉验证...")
  set.seed(42)
  folds     <- createFolds(1:nrow(points_sp), k = k_fold, list = TRUE)
  fold_rmse <- numeric(k_fold)

  for (fold in 1:k_fold) {
    test_idx  <- folds[[fold]]
    train_idx <- setdiff(1:nrow(points_sp), test_idx)

    model     <- loess(value ~ x * y,
                       data      = pts_df[train_idx, ],
                       span      = span, degree = 2,
                       normalize = TRUE, family = "gaussian")
    test_pred <- predict(model, newdata = pts_df[test_idx, ])

    fold_rmse[fold] <- sqrt(mean((test_pred - pts_df[test_idx, "value"])^2,
                                  na.rm = TRUE))
    message(sprintf("  Fold %d/%d 完成, RMSE: %.4f", fold, k_fold, fold_rmse[fold]))
  }

  cv_rmse_mean <- mean(fold_rmse)
  cv_rmse_sd   <- sd(fold_rmse)
  message(sprintf("\n交叉验证完成: RMSE = %.4f ± %.4f", cv_rmse_mean, cv_rmse_sd))

  # 全数据最终预测
  message("使用全部数据进行最终预测...")
  full_model <- loess(value ~ x * y,
                      data      = pts_df,
                      span      = span, degree = 2,
                      normalize = TRUE, family = "gaussian")

  pred_values <- inverse_transform(
    predict(full_model, newdata = grd_df), attrName
  )

  # 应用固定约束
  if (attrName %in% names(value_constraints)) {
    fr <- value_constraints[[attrName]]$fixed_range
    if (!is.null(fr)) pred_values <- pmin(pmax(pred_values, fr[1]), fr[2])
  }

  # 创建并裁剪栅格
  r         <- raster(extent(grd))
  res(r)    <- cellsize_deg
  values(r) <- NA
  values(r)[cellFromXY(r, coordinates(grd))] <- pred_values

  full_r <- r
  r      <- mask(r, boundary_sp, updatevalue = NA, touches = TRUE)

  message(sprintf("===== %s 处理完成 =====\n", attrName))

  return(list(
    raster      = r,
    full_raster = full_r,
    model       = full_model,
    cv_results  = list(
      rmse_mean      = cv_rmse_mean,
      rmse_sd        = cv_rmse_sd,
      rmse_all_folds = fold_rmse,
      n_folds        = k_fold
    )
  ))
}

# --- 2.4 可视化函数（无标题纯净版 · Y轴垂直 · 大字号） ---

raster_to_ggplot <- function(raster, title, rmse_mean, rmse_sd, boundary) {

  df           <- as.data.frame(rasterToPoints(raster))
  names(df)[3] <- "value"

  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = boundary, fill = NA, color = "black", linewidth = 0.5) +
    scale_fill_gradientn(
      colors = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100)
    ) +
    scale_x_continuous(breaks = scales::breaks_width(0.04)) +
    scale_y_continuous(breaks = scales::breaks_width(0.04)) +
    labs(fill = title, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.title       = element_blank(),
      axis.title       = element_text(size = 17, face = "bold"),
      axis.text.x      = element_text(size = 15, color = "black"),
      axis.text.y      = element_text(size = 15, color = "black",
                                      angle = 90, hjust = 0.5),
      legend.position  = "right",
      legend.title     = element_text(size = 15, face = "bold"),
      legend.text      = element_text(size = 14)
    )
}

#------------------------
#! 3. 常量与配置
#------------------------

cellsize_deg_30 <- 30 / 111320   # ≈ 0.000269°（30 米）

value_constraints <- list(
  "SR" = list(
    fixed_range = c(12, 67),
    transform   = "log",
    integer     = TRUE
  ),
  "delta_Gini_spatial" = list(
    fixed_range = c(-0.02933, 0.62206),
    transform   = "logit",
    integer     = FALSE
  ),
  "Gini_size" = list(
    fixed_range = c(0.4760, 0.9091),
    transform   = "logit",
    integer     = FALSE
  )
)

#------------------------
#! 4. 数据准备
#------------------------

# 读取边界
boundary <- st_read("QYS_Boundary.shp")

# 读取原始点数据，提取几何并保存为独立 gpkg
#points_updated <- st_read(
#  "F:/DYY-博士论文/GPT_测试_geodata/QYS_ALL_POINT_updated.gpkg"
#)
#geom_data <- points_updated %>% dplyr::select(PlotID, geom)

# 保存 geom_data（后续可直接从此读取，无需重新依赖原始大文件）
st_write(geom_data, "geom_data.gpkg", delete_dsn = TRUE)
geom_data <- st_read("geom_data.gpkg")   # ← 已保存后直接从这里读取

# 读取样地指标数据
df_metrics <- read_excel("110quadrat(Gini_LAC).xlsx")

df_metrics_sf <- df_metrics %>%
  left_join(geom_data, by = c("plot_id" = "PlotID")) %>%
  st_as_sf()



#------------------------
#! 5A. 执行空间插值（3个指标，10折CV）
#------------------------

# 1. SR
result_SR <- simulate_and_loss(
  boundary, df_metrics_sf , "SR", cellsize_deg = cellsize_deg_30
)
saveRDS(result_SR, "result_SR_loess_10fold.rds")

# 2. delta_Gini_spatial
result_delta_Gini_spatial <- simulate_and_loss(
  boundary, df_metrics_sf , "delta_Gini_spatial", cellsize_deg = cellsize_deg_30
)
saveRDS(result_delta_Gini_spatial, "result_delta_Gini_spatial_loess_10fold.rds")

# 3. Gini_size
result_Gini_size <- simulate_and_loss(
  boundary, df_metrics_sf , "Gini_size", cellsize_deg = cellsize_deg_30
)
saveRDS(result_Gini_size, "result_Gini_size_loess_10fold.rds")

#------------------------
#! 5B. 读取已保存的插值结果（5A 已跑过则从此处开始）
#------------------------

result_SR                 <- readRDS("result_SR_loess_10fold.rds")
result_delta_Gini_spatial <- readRDS("result_delta_Gini_spatial_loess_10fold.rds")
result_Gini_size          <- readRDS("result_Gini_size_loess_10fold.rds")

#------------------------
#! 6. 可视化 Fig. S4
#------------------------

p_SR <- raster_to_ggplot(
  result_SR$raster,
  "SR",
  result_SR$cv_results$rmse_mean,
  result_SR$cv_results$rmse_sd,
  boundary
)

p_delta_Gini_spatial <- raster_to_ggplot(
  result_delta_Gini_spatial$raster,
  "\u03b4Gini",
  result_delta_Gini_spatial$cv_results$rmse_mean,
  result_delta_Gini_spatial$cv_results$rmse_sd,
  boundary
)

p_Gini_size <- raster_to_ggplot(
  result_Gini_size$raster,
  "Size Gini",
  result_Gini_size$cv_results$rmse_mean,
  result_Gini_size$cv_results$rmse_sd,
  boundary
)

# 组合图表（1×3 横向）
interpolation_plot <- (p_SR + p_delta_Gini_spatial + p_Gini_size) +
  plot_layout(guides = "keep") +
  plot_annotation(
    theme = theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      plot.title  = element_text(size = 20, face = "bold",
                                 hjust = 0.5,
                                 margin = ggplot2::margin(b = 20))
    )
  )

ggsave("interpolation_results_3vars_10fold.pdf",
       interpolation_plot, width = 22, height = 6, dpi = 300)


#------------------------
#! 7. RMSE 摘要输出 Table S2
#------------------------

rmse_summary <- data.frame(
  Variable  = c("Species Richness", "\u03b4Gini", "Size Gini"),
  RMSE_Mean = c(
    result_SR$cv_results$rmse_mean,
    result_delta_Gini_spatial$cv_results$rmse_mean,
    result_Gini_size$cv_results$rmse_mean
  ),
  RMSE_SD   = c(
    result_SR$cv_results$rmse_sd,
    result_delta_Gini_spatial$cv_results$rmse_sd,
    result_Gini_size$cv_results$rmse_sd
  ),
  N_Folds   = c(
    result_SR$cv_results$n_folds,
    result_delta_Gini_spatial$cv_results$n_folds,
    result_Gini_size$cv_results$n_folds
  )
) %>%
  mutate(`RMSE (Mean ± SD)` = sprintf("%.4f ± %.4f", RMSE_Mean, RMSE_SD)) %>%
  arrange(RMSE_Mean)

print("========== RMSE Summary (10-Fold Cross-Validation) ==========")
print(rmse_summary)
write_xlsx(rmse_summary, "rmse_summary_3vars_10fold.xlsx")

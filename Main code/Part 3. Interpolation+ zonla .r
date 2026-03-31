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
#! 5A. 执行空间插值（3个指标，10折CV） Table S3
#------------------------

# 确保目录存在
if (!dir.exists("raster_data")) dir.create("raster_data")

# 1. SR
result_SR <- simulate_and_loss(
  boundary, df_metrics_sf, "SR", cellsize_deg = cellsize_deg_30
)
saveRDS(result_SR, file.path("raster_data", "result_SR_loess_10fold.rds"))
writeRaster(result_SR$raster, file.path("raster_data", "SR_interpolated.tif"), format = "GTiff", overwrite = TRUE)

# 2. delta_Gini_spatial
result_delta_Gini_spatial <- simulate_and_loss(
  boundary, df_metrics_sf, "delta_Gini_spatial", cellsize_deg = cellsize_deg_30
)
saveRDS(result_delta_Gini_spatial, file.path("raster_data", "result_delta_Gini_spatial_loess_10fold.rds"))
writeRaster(result_delta_Gini_spatial$raster, file.path("raster_data", "delta_Gini_spatial_interpolated.tif"), format = "GTiff", overwrite = TRUE)

# 3. Gini_size
result_Gini_size <- simulate_and_loss(
  boundary, df_metrics_sf, "Gini_size", cellsize_deg = cellsize_deg_30
)
saveRDS(result_Gini_size, file.path("raster_data", "result_Gini_size_loess_10fold.rds"))
writeRaster(result_Gini_size$raster, file.path("raster_data", "Gini_size_interpolated.tif"), format = "GTiff", overwrite = TRUE)

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
  "ΔGini Spatial",
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




##---
## zonal statistics
#####


# -----------------------------
#! 1. 加载所需包
# -----------------------------
library(sf)         # 处理矢量数据
library(terra)      # 替代raster包，更高效
library(dplyr)      # 数据处理
library(writexl)    # 写出 Excel 文件
library(ggplot2)    # 可视化
library(pbapply)    # 添加进度条
library(parallel)   # 并行处理

# -----------------------------
#! 2. 读取数据并转换为投影（单位：米）
# -----------------------------
# 读取边界数据
boundary <- st_read("QYS_Boundary.shp")

# 转换为投影坐标系（单位：米）- UTM zone 50N
boundary_proj <- st_transform(boundary, crs = 32650)

# 获取目录下所有 .tif 文件的完整路径
tif_dir <- "raster_data"
tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)

cat("找到", length(tif_files), "个栅格文件\n")


# -----------------------------
#! 3. 构建完整网格并筛选有效单元
# -----------------------------
# 生成网格并设置ID的函数
generate_grid <- function(boundary, cellsize) {
  # 生成完整网格
  grid_full <- st_make_grid(boundary, cellsize = cellsize, square = TRUE)
  
  # 使用st_intersects筛选与边界相交的网格
  intersects <- st_intersects(grid_full, boundary, sparse = FALSE)[,1]
  grid_valid <- grid_full[intersects]
  
  # 返回带有ID的sf对象
  st_sf(id = 1:length(grid_valid), geometry = grid_valid)
}

# 预先生成并缓存所有分辨率的网格
cat("生成网格...\n")
grid_1000 <- generate_grid(boundary_proj, 1000)
grid_500  <- generate_grid(boundary_proj, 500)


cat("生成的网格单元数: 1000m =", nrow(grid_1000), 
    ", 500m =", nrow(grid_500), "\n")


####################################################
####  Fig.S 5
#############################################################
library(ggplot2)
library(patchwork)

# 创建三个独立的图形
p1 <- ggplot() + 
  geom_sf(data = boundary_proj, fill = NA, color = "red", size = 1) +
  geom_sf(data = grid_1000, fill = NA, color = "black", size = 0.3) +
  labs(title = paste("(a) 100 ha Grid (n =", nrow(grid_1000), ")")) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

p2 <- ggplot() + 
  geom_sf(data = boundary_proj, fill = NA, color = "red", size = 1) +
  geom_sf(data = grid_500, fill = NA, color = "black", size = 0.2) +
  labs(title = paste("(b) 25 ha Grid (n =", nrow(grid_500), ")")) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())



# 组合图形
combined_plot <- p1 + p2 

print(combined_plot)

ggsave(
  filename = "Fig.S5.pdf",
  plot = combined_plot,
  width = 15,
  height = 5,
  dpi = 300,
  bg = "white"
)

cat("图形已保存完成！\n")




# -----------------------------
#! 4. 优化的区域统计函数
# -----------------------------

process_zonal_stats <- function(grid_sf, resolution_label, tif_files, boundary_proj) {
  # 转换为terra格式的向量对象，比raster更高效
  grid_vect <- terra::vect(grid_sf)
  output_dir <- dirname(tif_files[1])
  
  # 初始化结果数据框
  result_all <- data.frame(id = grid_sf$id)
  
  # 报告处理起始时间
  start_time <- Sys.time()
  cat(sprintf("\n开始处理 %s 分辨率的区域统计，开始时间: %s\n", 
              resolution_label, format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # 设置terra内存使用和多线程
  terra::terraOptions(memfrac = 0.7, verbose = FALSE)  # 使用70%的RAM
  
  # 如果有多核心可用，设置多线程处理
  if(parallel::detectCores() > 1) {
    n_cores <- max(1, parallel::detectCores() - 1)  # 保留一个核心给系统
    terra::terraOptions(threads = n_cores)
    cat(sprintf("启用多线程处理，使用 %d 个核心\n", n_cores))
  }
  
  # 处理每个栅格文件
  for(i in seq_along(tif_files)) {
    file_start_time <- Sys.time()
    file <- tif_files[i]
    base_name <- sub("^best_loess_", "", sub("\\.tif$", "", basename(file)))
    out_csv <- file.path(output_dir, paste0(base_name, "_zonal_stats_", resolution_label, ".csv"))
    
    # 使用tryCatch处理可能的错误
    tryCatch({
      # 使用terra读取栅格，更高效
      cat(sprintf("[%d/%d] 开始读取: %s\n", i, length(tif_files), basename(file)))
      r <- terra::rast(file)
      
      # 仅在必要时进行坐标转换
      if(terra::crs(r) != st_crs(boundary_proj)$wkt) {
        cat("  进行坐标系转换...\n")
        r <- terra::project(r, st_crs(boundary_proj)$wkt)
      }
      
      # 使用terra的更高效extract，并设置向量化处理
      cat("  计算区域统计...\n")
      mean_vals <- terra::extract(r, grid_vect, fun = "mean", na.rm = TRUE, touches = FALSE)
      
      # 仅保留值列（第二列）
      mean_vals <- mean_vals[, 2]
      
      # 生成结果数据框
      df <- data.frame(id = grid_sf$id, mean = mean_vals)
      
      # 写出CSV
      write.csv(df, out_csv, row.names = FALSE)
      
      # 合并到总结果中
      result_all[[base_name]] <- mean_vals
      
      # 清理内存
      r <- NULL
      gc(reset = TRUE)
      
      # 计算并显示处理时间
      file_end_time <- Sys.time()
      elapsed <- difftime(file_end_time, file_start_time, units = "secs")
      cat(sprintf("[%d/%d] 已处理: %s (用时: %.2f秒)\n", 
                  i, length(tif_files), basename(file), as.numeric(elapsed)))
      
    }, error = function(e) {
      cat(sprintf("[%d/%d] 处理文件时出错: %s\n错误信息: %s\n", 
                 i, length(tif_files), basename(file), e$message))
    })
    
    # 每处理5个文件保存一次总结果
    if(i %% 5 == 0) {
      temp_file <- file.path(output_dir, paste0("temp_results_", resolution_label, ".rds"))
      saveRDS(result_all, temp_file)
      
      # 计算进度和预计剩余时间
      current_time <- Sys.time()
      elapsed_total <- difftime(current_time, start_time, units = "mins")
      avg_time_per_file <- as.numeric(elapsed_total) / i
      remaining_files <- length(tif_files) - i
      est_remaining_time <- avg_time_per_file * remaining_files
      
      cat(sprintf("保存临时结果到 %s\n", basename(temp_file)))
      cat(sprintf("进度: %.1f%% (已用时: %.1f分钟, 预计剩余: %.1f分钟)\n", 
                  i/length(tif_files)*100, 
                  as.numeric(elapsed_total),
                  est_remaining_time))
    }
  }
  
  # 处理完所有文件后保存最终结果
  final_file <- file.path(output_dir, paste0("final_results_", resolution_label, ".rds"))
  saveRDS(result_all, final_file)
  
  # 计算总处理时间
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("已完成所有处理，保存最终结果到 %s\n", basename(final_file)))
  cat(sprintf("总处理时间: %.2f分钟\n", as.numeric(total_time)))
  
  return(result_all)
}

# -----------------------------
#! 5. 逐个分辨率处理（使用优化函数）
# -----------------------------
# 为防止内存问题，分别处理每个分辨率

# 处理1000m网格
df_1000_all <- process_zonal_stats(grid_1000, "1000m", tif_files, boundary_proj)
# 处理完后清理一下内存
gc(reset = TRUE)
head(df_1000_all)

# 处理500m网格
df_500_all <- process_zonal_stats(grid_500, "500m", tif_files, boundary_proj)
gc(reset = TRUE)


# -----------------------------
#! 6. 导出结果为Excel和Shapefile  
#

# 先关联原始网格，恢复sf对象
cat("正在准备shapefile数据...\n")
df_1000_sf <- left_join(grid_1000, df_1000_all, by = "id")
df_500_sf  <- left_join(grid_500, df_500_all, by = "id")

# 然后删除包含NA的行
df_1000_clean_sf <- df_1000_sf[complete.cases(st_drop_geometry(df_1000_sf)), ]
df_500_clean_sf  <- df_500_sf[complete.cases(st_drop_geometry(df_500_sf)), ]


# 创建非空间数据框用于Excel导出
df_1000_clean <- st_drop_geometry(df_1000_clean_sf)
df_500_clean <- st_drop_geometry(df_500_clean_sf)


# 显示删除后的行数
cat("删除NA后的行数：\n")
cat("1000m数据：", nrow(df_1000_clean), "行\n")
cat("500m数据：", nrow(df_500_clean), "行\n")

library(writexl)
# 保存Excel
cat("正在保存Excel文件...\n")
write_xlsx(
  list(
    "1000m" = df_1000_clean,
    "500m"  = df_500_clean
  ),
  path = file.path("zonal_stats_all_scales.xlsx")
)

# 写出Shapefile/GeoPackage文件
cat("正在保存GeoPackage文件...\n")
st_write(df_1000_clean_sf, file.path( "zonal_stats_1000m.gpkg"), delete_layer = TRUE)
st_write(df_500_clean_sf, file.path( "zonal_stats_500m.gpkg"), delete_layer = TRUE)


cat("处理完成！所有结果已保存。\n")


##! Zonal exhibit###


#=========================================================
# Fig. S6   空间多边形数据双尺度对比批量制图
# 尺度: 100 ha (1000m) vs 25 ha (500m)
# 排版: 4行(指标) x 2列(尺度)
#=========================================================

#------------------------
# 1. 加载包
#------------------------
library(sf)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

#------------------------
# 2. 定义路径与加载数据
#------------------------
# 确保 output_shp_dir 已经在环境中定义好了

boundary_path <- "QYS_Boundary.shp"
output_pdf <- file.path( "Fig. S6 Grid_Metrics_MultiScale.pdf")

# 读取双尺度的 gpkg 数据
df_1000_sf <- st_read(file.path("zonal_stats_1000m.gpkg"))
df_500_sf  <- st_read(file.path("zonal_stats_500m.gpkg"))

# 读取边界数据
boundary_sf <- st_read(boundary_path)

#------------------------
# 3. 坐标系统一与转换
#------------------------
# 统一转换为 WGS84 经纬度 (EPSG:4326)，确保 0.04 度网格生效
df_1000_sf  <- st_transform(df_1000_sf, crs = 4326)
df_500_sf   <- st_transform(df_500_sf, crs = 4326)
boundary_sf <- st_transform(boundary_sf, crs = 4326)

#------------------------
# 4. 封装统一的学术绘图函数 (仅图例放大版)
#------------------------
plot_grid_metric <- function(sf_data, var_name, legend_title, color_theme = "eco", boundary = NULL, plot_title = NULL) {
  
  # 配置颜色方案
  if (color_theme == "eco") {
    fill_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
  } else if (color_theme == "prod") {
    fill_colors <- c("#d73027", "#fee08b", "#1a9850")
  }
  
  # 基础绘图
  p <- ggplot() +
    geom_sf(data = sf_data, aes(fill = .data[[var_name]]), color = "white", linewidth = 0.05)
  
  # 叠加边界
  if (!is.null(boundary)) {
    p <- p + geom_sf(data = boundary, fill = NA, color = "black", linewidth = 0.8)
  }
  
  # 格式设置
  p <- p + 
    # 【修改点 1】: 加入 guide_colorbar 稍微加宽加长色条，配合放大的文字
    scale_fill_gradientn(
      colors = fill_colors, 
      name = legend_title,
      guide = guide_colorbar(barwidth = unit(1.2, "cm"), barheight = unit(6, "cm"))
    ) +
    scale_x_continuous(breaks = scales::breaks_width(0.04)) +
    scale_y_continuous(breaks = scales::breaks_width(0.04)) +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black", angle = 90, hjust = 0.5),
      
      # 【修改点 2】: 图例文字放大
      legend.position = "right",
      legend.title = element_text(size = 22, face = "bold"), # 17 放大到 22
      legend.text = element_text(size = 18)                  # 16 放大到 18
    )
  
  # 如果传入了标题，则显示标题
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title) +
      theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  } else {
    p <- p + theme(plot.title = element_blank())
  }
  
  return(p)
}

#------------------------
# 5. 生成各个指标的图表 (双尺度)
#------------------------

# --- 1000m (100 ha) 尺度 ---
p_sr_1000     <- plot_grid_metric(df_1000_sf, "SR_interpolated", "SR", "eco", boundary_sf, plot_title = "100 ha")
p_t_gini_1000 <- plot_grid_metric(df_1000_sf, "delta_Gini_spatial_interpolated", "ΔGini spatial", "eco", boundary_sf)
p_size_g_1000 <- plot_grid_metric(df_1000_sf, "Gini_size_interpolated", "Gini size", "eco", boundary_sf)
p_npp_1000    <- plot_grid_metric(df_1000_sf, "NPP", bquote("NPP"), "prod", boundary_sf)

# --- 500m (25 ha) 尺度 ---
p_sr_500     <- plot_grid_metric(df_500_sf, "SR_interpolated", "SR", "eco", boundary_sf, plot_title = "25 ha")
p_t_gini_500 <- plot_grid_metric(df_500_sf, "delta_Gini_spatial_interpolated", "ΔGini spatial", "eco", boundary_sf)
p_size_g_500 <- plot_grid_metric(df_500_sf, "Gini_size_interpolated", "Gini size", "eco", boundary_sf)
p_npp_500    <- plot_grid_metric(df_500_sf, "NPP", bquote("NPP"), "prod", boundary_sf)

#------------------------
# 6. 图表拼接与导出 (4行 x 2列)
#------------------------
combined_plot <- (p_sr_1000     | p_sr_500) /
                 (p_t_gini_1000 | p_t_gini_500) /
                 (p_size_g_1000 | p_size_g_500) / 
                 (p_npp_1000    | p_npp_500) +
                 plot_layout(guides = "keep") # 保持各自独立的图例

# 保存为高质量 PDF
ggsave(output_pdf, 
       plot = combined_plot, 
       width = 22,    
       height = 36, 
       limitsize = FALSE, 
       dpi = 300)

message("===== 绘图完成！图例放大版双尺度对比文件已保存至: ", output_pdf, " =====")

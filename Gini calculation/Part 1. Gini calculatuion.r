############################################
##########这个Code主要介绍delta_Gini_spatial 和 Gini_size以及其对应的lorenz asymetric coefficient （S，在代码中用LAC表示）。由于政策原因，仅使用Gird 74作为workflow展示，对应原始过程数据中的grid_id ==74, 由于防止点重复而设置了抖动，
########## 因此，尽管设置了seed，您在本次计算得出的 Spatial structural diversity 相关指标均可能会与参与计算时使用的数据产生轻微的出入##################  
##################  



# =====================================================================
# 群落生态学空间格局与大小等级分析 — Demo Version: Quadrat 74
# (Spatial Voronoi Area, Gini Index, & Lorenz Asymmetry Coefficient)
# =====================================================================

library(dplyr)
library(readxl)
library(writexl)
library(spatstat.geom)
library(spatstat.random)
library(deldir)
library(ggplot2)
library(tidyr)

# =====================================================================
# 1. 核心统计算法定义
# =====================================================================

calc_gini <- function(x) {
  x <- sort(na.omit(x))
  n <- length(x)
  if (n < 2) return(NA_real_)
  if (any(x < 0)) stop("Gini requires non-negative values.")
  sx <- sum(x)
  if (sx <= 0) return(NA_real_)
  sum((2 * seq_len(n) - n - 1) * x) / ((n - 1) * sx)
}

calc_lac <- function(x, tol = sqrt(.Machine$double.eps)) {
  x_sorted <- sort(na.omit(x))
  n <- length(x_sorted)
  if (n < 2) return(NA_real_)
  if (any(x_sorted < 0)) stop("LAC requires non-negative values.")
  mu <- mean(x_sorted)
  Ln <- sum(x_sorted)
  if (Ln <= 0) return(NA_real_)
  m  <- sum(x_sorted < mu - tol)
  a  <- sum(abs(x_sorted - mu) <= tol)
  Lm <- if (m > 0) sum(x_sorted[1:m]) else 0
  if (a > 0) {
    L_m_plus_a <- sum(x_sorted[1:(m + a)])
    S_lower <- (m / n) + (Lm / Ln)
    S_upper <- ((m + a) / n) + (L_m_plus_a / Ln)
    return((S_lower + S_upper) / 2)
  } else {
    x_m        <- x_sorted[m]
    x_m_plus_1 <- x_sorted[m + 1]
    d    <- (mu - x_m) / (x_m_plus_1 - x_m)
    F_mu <- (m + d) / n
    L_mu <- (Lm + d * x_m_plus_1) / Ln
    return(F_mu + L_mu)
  }
}

# =====================================================================
# 2. 数据读取、筛选、重命名、写出、重新读入
# =====================================================================

data_path <- "F:/DYY-博士论文/1. QY_SAR/QY_SAR_PROJECT/5. 小论文工作/5.4 性状数据的获取/QYS_QYS_物种名称纠正/data_with_STD_NAME_final_AGB.xlsx"
data_raw <- read_excel(data_path)

data_grid74 <- data_raw %>%
  filter(网格号 == 74) %>%
  select(
    grid_id      = 网格号,
    sp_name      = SP_NAME,
    dbh_cm       = 胸径cm,
    height_m     = 高度m,
    wood_density = WoodDensity,
    agb_kg       = ABG_kg,        # ← 原 ABG_kg，后续全部用 agb_kg
    sim_x        = 模拟坐标X,
    sim_y        = 模拟坐标Y,
    abs_x        = 绝对坐标X,
    abs_y        = 绝对坐标Y
  ) %>%
  filter(!is.na(sp_name), dbh_cm > 0, height_m > 0, agb_kg >= 0) %>%
  mutate(unique_id = row_number())

cat(sprintf("Grid 74 | 保留个体: %d 株 | 物种数: %d 种\n",
            nrow(data_grid74), n_distinct(data_grid74$sp_name)))
glimpse(data_grid74)

path_grid74 <- "Grid74_renamed.xlsx"
write_xlsx(data_grid74, path_grid74)
message("已保存至: ", path_grid74)

data <- read_excel(path_grid74)
cat("重新读入成功，列名如下：\n")
print(names(data))

# =====================================================================
# 3. 初步清洗：边界检查
# =====================================================================

out_of_bounds <- data %>%
  filter(abs_x < 0 | abs_x > 25 | abs_y < 0 | abs_y > 25)

cat("超出边界的点（共", nrow(out_of_bounds), "条）：\n")
print(out_of_bounds)

data_clean <- data %>%
  filter(abs_x >= 0 & abs_x <= 25 & abs_y >= 0 & abs_y <= 25) %>%
  mutate(unique_id = row_number())

# =====================================================================
# 4. 空间 Voronoi 多边形拟合与面积提取
# =====================================================================

window <- owin(c(0, 25), c(0, 25))

convert_to_ppp <- function(data, x_col, y_col, id_col, grid_col) {
  grids <- unique(data[[grid_col]])
  ppp_list <- lapply(grids, function(grid) {
    grid_data <- data %>% filter(!!sym(grid_col) == grid)
    ppp(grid_data[[x_col]], grid_data[[y_col]],
        window = window, marks = grid_data[[id_col]])
  })
  names(ppp_list) <- grids
  return(ppp_list)
}
set.seed(42)  
ppp_absolute  <- convert_to_ppp(data_clean, "abs_x", "abs_y", "unique_id", "grid_id")
ppp_simulated <- convert_to_ppp(data_clean, "sim_x", "sim_y", "unique_id", "grid_id")

apply_jitter_to_ppp <- function(ppp_objects, jitter_amount = 0.2) {
  lapply(ppp_objects, function(ppp_obj)
    rjitter(ppp_obj, radius = jitter_amount, retry = TRUE))
}

ppp_absolute_jittered <- apply_jitter_to_ppp(ppp_absolute, jitter_amount = 0.2)

generate_voronoi_with_ppp <- function(ppp_objects) {
  voronoi_areas <- lapply(ppp_objects, function(ppp_obj) {
    x_coords <- ppp_obj$x
    y_coords <- ppp_obj$y
    ids      <- marks(ppp_obj)
    voronoi_result <- deldir(x_coords, y_coords,
                             rw = c(0, 25, 0, 25), suppressMsge = TRUE)
    tiles <- tile.list(voronoi_result)
    areas <- sapply(tiles, function(tile) {
      area.owin(owin(poly = list(x = tile$x, y = tile$y)))
    })
    data.frame(unique_id = ids, clipped_area = areas)
  })
  do.call(rbind, voronoi_areas)
}

voronoi_absolute_areas  <- generate_voronoi_with_ppp(ppp_absolute_jittered)
voronoi_simulated_areas <- generate_voronoi_with_ppp(ppp_simulated)

final_data <- data_clean %>%
  left_join(voronoi_absolute_areas  %>% rename(observed_area  = clipped_area), by = "unique_id") %>%
  left_join(voronoi_simulated_areas %>% rename(simulated_area = clipped_area), by = "unique_id")

# =====================================================================
# 5. 核心汇总：Plot-level Gini 和 LAC
# ✅ 修正: ABG_kg → agb_kg（重命名后的正确列名）
# =====================================================================

final_metrics_per_grid <- final_data %>%
  group_by(grid_id) %>%
  summarise(
    N_individuals        = n(),
    total_observed_area  = sum(observed_area,  na.rm = TRUE),
    total_simulated_area = sum(simulated_area, na.rm = TRUE),

    Gini_spatial_obs  = calc_gini(observed_area),
    LAC_spatial_obs   = calc_lac(observed_area),

    Gini_spatial_sim  = calc_gini(simulated_area),
    LAC_spatial_sim   = calc_lac(simulated_area),

    delta_Gini_spatial = Gini_spatial_obs - Gini_spatial_sim,

    Gini_size = calc_gini(agb_kg),   
    LAC_size  = calc_lac(agb_kg),   

    .groups = 'drop'
  )

print("计算完成！汇总数据如下：")
print(final_metrics_per_grid)

# =====================================================================
# 6. 结果导出
# =====================================================================

write_xlsx(
  final_data,
  "Grid74_Voronoi_Area_Individual.xlsx"
)
write_xlsx(
  final_metrics_per_grid,
  "Grid74_Gini_LAC_Summary.xlsx"
)

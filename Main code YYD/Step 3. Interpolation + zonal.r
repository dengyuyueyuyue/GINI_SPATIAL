# ============================================================
# Supplementary spatial workflow:
# interpolation, zonal statistics, and multi-scale mapping
# ============================================================
# This script performs three linked spatial analyses:
# 1) LOESS interpolation of SR, delta_Gini_spatial, and Gini_size
#    using 10-fold cross-validation;
# 2) zonal statistics to aggregate interpolated rasters to
#    100-ha and 25-ha analytical grids;
# 3) figure export for interpolated surfaces and grid-level maps.
#
# Inputs:
#   1) QYS_Boundary.shp
#   2) geom_data.gpkg
#   3) 110quadrat(Gini_LAC).xlsx
#   4) raster layers in raster_data/
#
# Main outputs:
#   Fig. S4  - Interpolated maps of SR, delta_Gini_spatial, and Gini_size
#   Fig. S5  - Analytical grids at 100-ha and 25-ha scales
#   Fig. S6  - Grid-level maps of interpolated predictors and NPP
#   Table S3 - Cross-validation performance of LOESS interpolation
#   zonal_stats_all_scales.xlsx
#   zonal_stats_1000m.gpkg
#   zonal_stats_500m.gpkg
#
# Key variables:
#   SR                  = species richness
#   delta_Gini_spatial  = observed minus simulated spatial Gini
#   Gini_size           = size inequality within quadrats
#
# Note:
# Interpolation is based on confidentiality-safe quadrat-level inputs.
# The full analytical loop across all variables and scales was run
# separately; this script presents the core spatial workflow.
# ============================================================

rm(list = ls())

# ---------------------------
# 0. Packages
# ---------------------------
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

# ---------------------------
# 1. Core helper functions
# ---------------------------

# Apply variable-specific transformation before interpolation.
# Transformations are defined in `value_constraints`.
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

# Back-transform LOESS predictions to the original data scale.
# Integer variables are rounded after inverse transformation.
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

# Create a regular prediction grid covering the study boundary.
# Grid resolution is controlled by `cellsize_deg`.
create_prediction_grid <- function(boundary_sp, cellsize_deg) {
  r <- raster(extent(boundary_sp))
  res(r) <- cellsize_deg
  grd <- rasterToPoints(r, spatial = TRUE)
  proj4string(grd) <- proj4string(boundary_sp)
  gridded(grd) <- TRUE
  return(grd)
}

# Interpolate one structural metric using LOESS and evaluate
# predictive performance with k-fold cross-validation.
#
# Steps:
# 1) harmonize coordinate systems;
# 2) transform the response variable if required;
# 3) fit LOESS models under k-fold cross-validation;
# 4) calculate fold-specific and mean RMSE;
# 5) refit the final model using all observations;
# 6) predict to a regular grid and clip the raster to the boundary.
#
# Output:
# - interpolated raster
# - full raster before masking
# - fitted LOESS model
# - cross-validation summary
simulate_and_loss <- function(boundary, points, attrName,
                              cellsize_deg = cellsize_deg_30,
                              k_fold = 10,
                              span   = 0.2) {

  message(sprintf("\n===== Processing: %s =====", attrName))

  if (st_crs(boundary) != st_crs(points)) {
    points <- st_transform(points, st_crs(boundary))
  }

  boundary_sp <- as(boundary, "Spatial")
  points_sp   <- as(points,   "Spatial")

  grd    <- create_prediction_grid(boundary_sp, cellsize_deg)
  grd_df <- as.data.frame(coordinates(grd))
  names(grd_df) <- c("x", "y")

  original_values             <- points_sp@data[[attrName]]
  points_sp$transformed_value <- transform_data(original_values, attrName)

  pts_coords        <- as.data.frame(coordinates(points_sp))
  names(pts_coords) <- c("x", "y")
  pts_df            <- cbind(pts_coords, value = points_sp$transformed_value)

  message("Running 10-fold cross-validation...")
  set.seed(42)
  folds     <- createFolds(1:nrow(points_sp), k = k_fold, list = TRUE)
  fold_rmse <- numeric(k_fold)

  for (fold in 1:k_fold) {
    test_idx  <- folds[[fold]]
    train_idx <- setdiff(1:nrow(points_sp), test_idx)

    model <- loess(value ~ x * y,
                   data      = pts_df[train_idx, ],
                   span      = span,
                   degree    = 2,
                   normalize = TRUE,
                   family    = "gaussian")

    test_pred <- predict(model, newdata = pts_df[test_idx, ])

    fold_rmse[fold] <- sqrt(mean((test_pred - pts_df[test_idx, "value"])^2,
                                 na.rm = TRUE))
    message(sprintf("  Fold %d/%d completed, RMSE: %.4f",
                    fold, k_fold, fold_rmse[fold]))
  }

  cv_rmse_mean <- mean(fold_rmse)
  cv_rmse_sd   <- sd(fold_rmse)
  message(sprintf("\nCross-validation completed: RMSE = %.4f ± %.4f",
                  cv_rmse_mean, cv_rmse_sd))

  message("Fitting the final model using all observations...")
  full_model <- loess(value ~ x * y,
                      data      = pts_df,
                      span      = span,
                      degree    = 2,
                      normalize = TRUE,
                      family    = "gaussian")

  pred_values <- inverse_transform(
    predict(full_model, newdata = grd_df), attrName
  )

  if (attrName %in% names(value_constraints)) {
    fr <- value_constraints[[attrName]]$fixed_range
    if (!is.null(fr)) pred_values <- pmin(pmax(pred_values, fr[1]), fr[2])
  }

  r         <- raster(extent(grd))
  res(r)    <- cellsize_deg
  values(r) <- NA
  values(r)[cellFromXY(r, coordinates(grd))] <- pred_values

  full_r <- r
  r      <- mask(r, boundary_sp, updatevalue = NA, touches = TRUE)

  message(sprintf("===== %s completed =====\n", attrName))

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

# Convert an interpolated raster into a publication-style map.
# Used for Fig. S4.
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

# ---------------------------
# 2. Analysis settings
# ---------------------------
# Define interpolation resolution and variable-specific constraints.

cellsize_deg_30 <- 30 / 111320

# Variable-specific ranges and transformation rules used during interpolation.
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

# ---------------------------
# 3. Input data preparation
# ---------------------------
# Read the study boundary, quadrat geometry, and quadrat-level
# structural metrics used for interpolation.

boundary <- st_read("QYS_Boundary.shp")

# Quadrat geometry is stored as a lightweight GeoPackage for routine reruns.
geom_data <- st_read("geom_data.gpkg")

df_metrics <- read_excel("110quadrat(Gini_LAC).xlsx")

df_metrics_sf <- df_metrics %>%
  left_join(geom_data, by = c("plot_id" = "PlotID")) %>%
  st_as_sf()

# ---------------------------
# 4. LOESS interpolation with 10-fold cross-validation
# ---------------------------
# This section interpolates SR, delta_Gini_spatial, and Gini_size
# across the study boundary and exports both raster outputs and
# cross-validation results.

if (!dir.exists("raster_data")) dir.create("raster_data")

# Interpolate species richness (SR)
result_SR <- simulate_and_loss(
  boundary, df_metrics_sf, "SR", cellsize_deg = cellsize_deg_30
)
#saveRDS(result_SR, file.path("raster_data", "result_SR_loess_10fold.rds"))
#writeRaster(result_SR$raster,
#            file.path("raster_data", "SR_interpolated.tif"),
#            format = "GTiff", overwrite = TRUE)

# Interpolate spatial structural diversity (delta_Gini_spatial)
result_delta_Gini_spatial <- simulate_and_loss(
  boundary, df_metrics_sf, "delta_Gini_spatial", cellsize_deg = cellsize_deg_30
)
#saveRDS(result_delta_Gini_spatial,
#        file.path("raster_data", "result_delta_Gini_spatial_loess_10fold.rds"))
#writeRaster(result_delta_Gini_spatial$raster,
#            file.path("raster_data", "delta_Gini_spatial_interpolated.tif"),
#            format = "GTiff", overwrite = TRUE)

# Interpolate size structural diversity (Gini_size)
result_Gini_size <- simulate_and_loss(
  boundary, df_metrics_sf, "Gini_size", cellsize_deg = cellsize_deg_30
)
#saveRDS(result_Gini_size,
#        file.path("raster_data", "result_Gini_size_loess_10fold.rds"))
#writeRaster(result_Gini_size$raster,
#           file.path("raster_data", "Gini_size_interpolated.tif"),
#           format = "GTiff", overwrite = TRUE)

# ---------------------------
# 5. Reload saved interpolation results
# ---------------------------
# Use this section when interpolation has already been completed
# and only figure generation or downstream analyses are needed.

result_SR <- readRDS(file.path("raster_data", "result_SR_loess_10fold.rds"))
result_delta_Gini_spatial <- readRDS(file.path("raster_data", "result_delta_Gini_spatial_loess_10fold.rds"))
result_Gini_size <- readRDS(file.path("raster_data", "result_Gini_size_loess_10fold.rds"))

# ---------------------------
# 6. Fig. S4
# ---------------------------
# Visualize interpolated surfaces for SR, delta_Gini_spatial,
# and Gini_size.

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

# Combine the three interpolated maps into a single horizontal panel.
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

ggsave("Fig. S4 interpolation_results_3vars_10fold.pdf",
       interpolation_plot, width = 22, height = 6, dpi = 300)

# ============================================================
# Zonal statistics and multi-scale grid aggregation
# ============================================================
# This section aggregates interpolated raster layers to
# 100-ha and 25-ha analytical grids for downstream analyses.
# ============================================================

# ---------------------------
# 7. Packages for zonal statistics
# ---------------------------
library(sf)
library(terra)
library(dplyr)
library(writexl)
library(ggplot2)
library(pbapply)
library(parallel)

# ---------------------------
# 8. Boundary projection and raster discovery
# ---------------------------
# Convert the study boundary to a projected coordinate system
# so that grid cell size is defined in meters.

boundary <- st_read("QYS_Boundary.shp")
boundary_proj <- st_transform(boundary, crs = 32650)

tif_dir <- "raster_data"
tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)

cat("Found", length(tif_files), "raster files\n")

# ---------------------------
# 9. Fig. S5
# ---------------------------
# Generate analytical grids and visualize the 100-ha and 25-ha
# grid systems used for zonal statistics.

# Generate square analytical grids and retain only cells that
# intersect the study boundary.
generate_grid <- function(boundary, cellsize) {
  grid_full <- st_make_grid(boundary, cellsize = cellsize, square = TRUE)
  intersects <- st_intersects(grid_full, boundary, sparse = FALSE)[, 1]
  grid_valid <- grid_full[intersects]
  st_sf(id = 1:length(grid_valid), geometry = grid_valid)
}

# Construct analytical grids corresponding to 100 ha and 25 ha.
cat("Generating analytical grids...\n")
grid_1000 <- generate_grid(boundary_proj, 1000)
grid_500  <- generate_grid(boundary_proj, 500)

cat("Number of grid cells: 1000 m =", nrow(grid_1000),
    ", 500 m =", nrow(grid_500), "\n")

library(ggplot2)
library(patchwork)

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

combined_plot <- p1 + p2
print(combined_plot)

ggsave(
  filename = "Fig. S5 zonal.pdf",
  plot = combined_plot,
  width = 15,
  height = 5,
  dpi = 300,
  bg = "white"
)

cat("Fig. S5 saved\n")

# ---------------------------
# 10. Zonal statistics at two spatial scales
# ---------------------------
# Extract mean raster values for the 100-ha and 25-ha grids.

# Extract mean raster values for each analytical grid cell.
#
# For each raster layer, this function:
# 1) reads the raster;
# 2) reprojects it if needed;
# 3) calculates mean values within grid cells;
# 4) exports layer-specific CSV files;
# 5) combines all extracted values into one grid-level table.
process_zonal_stats <- function(grid_sf, resolution_label, tif_files, boundary_proj) {
  grid_vect <- terra::vect(grid_sf)
  output_dir <- dirname(tif_files[1])

  result_all <- data.frame(id = grid_sf$id)

  start_time <- Sys.time()
  cat(sprintf("\nProcessing %s zonal statistics, start time: %s\n",
              resolution_label, format(start_time, "%Y-%m-%d %H:%M:%S")))

  terra::terraOptions(memfrac = 0.7, verbose = FALSE)

  if (parallel::detectCores() > 1) {
    n_cores <- max(1, parallel::detectCores() - 1)
    terra::terraOptions(threads = n_cores)
    cat(sprintf("Using %d CPU cores\n", n_cores))
  }

  for (i in seq_along(tif_files)) {
    file_start_time <- Sys.time()
    file <- tif_files[i]
    base_name <- sub("^best_loess_", "", sub("\\.tif$", "", basename(file)))
    out_csv <- file.path(output_dir,
                         paste0(base_name, "_zonal_stats_", resolution_label, ".csv"))

    tryCatch({
      cat(sprintf("[%d/%d] Reading: %s\n", i, length(tif_files), basename(file)))
      r <- terra::rast(file)

      if (terra::crs(r) != st_crs(boundary_proj)$wkt) {
        cat("  Reprojecting raster...\n")
        r <- terra::project(r, st_crs(boundary_proj)$wkt)
      }

      cat("  Extracting mean values...\n")
      mean_vals <- terra::extract(r, grid_vect, fun = "mean",
                                  na.rm = TRUE, touches = FALSE)

      mean_vals <- mean_vals[, 2]
      df <- data.frame(id = grid_sf$id, mean = mean_vals)

      write.csv(df, out_csv, row.names = FALSE)
      result_all[[base_name]] <- mean_vals

      r <- NULL
      gc(reset = TRUE)

      file_end_time <- Sys.time()
      elapsed <- difftime(file_end_time, file_start_time, units = "secs")
      cat(sprintf("[%d/%d] Finished: %s (%.2f s)\n",
                  i, length(tif_files), basename(file), as.numeric(elapsed)))

    }, error = function(e) {
      cat(sprintf("[%d/%d] Error in %s\nMessage: %s\n",
                  i, length(tif_files), basename(file), e$message))
    })

    if (i %% 5 == 0) {
      temp_file <- file.path(output_dir,
                             paste0("temp_results_", resolution_label, ".rds"))
      saveRDS(result_all, temp_file)

      current_time <- Sys.time()
      elapsed_total <- difftime(current_time, start_time, units = "mins")
      avg_time_per_file <- as.numeric(elapsed_total) / i
      remaining_files <- length(tif_files) - i
      est_remaining_time <- avg_time_per_file * remaining_files

      cat(sprintf("Temporary results saved to %s\n", basename(temp_file)))
      cat(sprintf("Progress: %.1f%% (elapsed: %.1f min, remaining: %.1f min)\n",
                  i / length(tif_files) * 100,
                  as.numeric(elapsed_total),
                  est_remaining_time))
    }
  }

  final_file <- file.path(output_dir,
                          paste0("final_results_", resolution_label, ".rds"))
  saveRDS(result_all, final_file)

  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("All processing completed, final results saved to %s\n",
              basename(final_file)))
  cat(sprintf("Total processing time: %.2f min\n", as.numeric(total_time)))

  return(result_all)
}

df_1000_all <- process_zonal_stats(grid_1000, "1000m", tif_files, boundary_proj)
gc(reset = TRUE)
head(df_1000_all)

df_500_all <- process_zonal_stats(grid_500, "500m", tif_files, boundary_proj)
gc(reset = TRUE)

# ---------------------------
# 11. Export grid-level outputs
# ---------------------------
# Save grid-level summaries as Excel tables and spatial layers
# for downstream modelling and visualization.
#
# Grid cells containing missing values in any extracted layer are removed
# to produce complete analytical datasets for subsequent analyses.

cat("Preparing grid-level spatial outputs...\n")
df_1000_sf <- left_join(grid_1000, df_1000_all, by = "id")
df_500_sf  <- left_join(grid_500,  df_500_all,  by = "id")

df_1000_clean_sf <- df_1000_sf[complete.cases(st_drop_geometry(df_1000_sf)), ]
df_500_clean_sf  <- df_500_sf[complete.cases(st_drop_geometry(df_500_sf)), ]

df_1000_clean <- st_drop_geometry(df_1000_clean_sf)
df_500_clean  <- st_drop_geometry(df_500_clean_sf)

cat("Rows retained after removing missing values:\n")
cat("1000 m:", nrow(df_1000_clean), "rows\n")
cat("500 m :", nrow(df_500_clean), "rows\n")

library(writexl)

cat("Saving Excel output...\n")
write_xlsx(
  list(
    "1000m" = df_1000_clean,
    "500m"  = df_500_clean
  ),
  path = file.path("zonal_stats_all_scales.xlsx") # PROCESS DATA
)

cat("Saving GeoPackage outputs...\n")
st_write(df_1000_clean_sf, file.path("zonal_stats_1000m.gpkg"), delete_layer = TRUE)
st_write(df_500_clean_sf,  file.path("zonal_stats_500m.gpkg"),  delete_layer = TRUE)

cat("Grid-level outputs saved\n")

# ============================================================
# Fig. S6
# Multi-scale maps of interpolated predictors and NPP
# ============================================================

# ---------------------------
# 12. Packages for grid-level mapping
# ---------------------------
library(sf)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

# ---------------------------
# 13. Read spatial layers for mapping
# ---------------------------
boundary_path <- "QYS_Boundary.shp"
output_pdf <- file.path("Fig. S6 Grid_Metrics_MultiScale.pdf")

df_1000_sf <- st_read(file.path("zonal_stats_1000m.gpkg"))
df_500_sf  <- st_read(file.path("zonal_stats_500m.gpkg"))
boundary_sf <- st_read(boundary_path)

# ---------------------------
# 14. Coordinate harmonization
# ---------------------------
# Convert all spatial layers to WGS84 for consistent map display.
df_1000_sf  <- st_transform(df_1000_sf, crs = 4326)
df_500_sf   <- st_transform(df_500_sf,  crs = 4326)
boundary_sf <- st_transform(boundary_sf, crs = 4326)

# ---------------------------
# 15. Build maps for each variable and scale
# ---------------------------

# Plot one grid-level variable for one spatial scale.
# Used to assemble the multi-panel map in Fig. S6.
plot_grid_metric <- function(sf_data, var_name, legend_title,
                             color_theme = "eco",
                             boundary = NULL,
                             plot_title = NULL) {

  if (color_theme == "eco") {
    fill_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
  } else if (color_theme == "prod") {
    fill_colors <- c("#d73027", "#fee08b", "#1a9850")
  }

  p <- ggplot() +
    geom_sf(data = sf_data,
            aes(fill = .data[[var_name]]),
            color = "white", linewidth = 0.05)

  if (!is.null(boundary)) {
    p <- p + geom_sf(data = boundary, fill = NA, color = "black", linewidth = 0.8)
  }

  p <- p +
    scale_fill_gradientn(
      colors = fill_colors,
      name = legend_title,
      guide = guide_colorbar(barwidth = unit(1.2, "cm"),
                             barheight = unit(6, "cm"))
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
      axis.text.y = element_text(size = 16, color = "black",
                                 angle = 90, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 18)
    )

  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title) +
      theme(plot.title = element_text(size = 24, face = "bold",
                                      hjust = 0.5,
                                      margin = margin(b = 15)))
  } else {
    p <- p + theme(plot.title = element_blank())
  }

  return(p)
}

# 1000 m (100 ha)
p_sr_1000     <- plot_grid_metric(df_1000_sf, "SR_interpolated", "SR", "eco", boundary_sf, plot_title = "100 ha")
p_t_gini_1000 <- plot_grid_metric(df_1000_sf, "delta_Gini_spatial_interpolated", "ΔGini spatial", "eco", boundary_sf)
p_size_g_1000 <- plot_grid_metric(df_1000_sf, "Gini_size_interpolated", "Gini size", "eco", boundary_sf)
p_npp_1000    <- plot_grid_metric(df_1000_sf, "NPP", bquote("NPP"), "prod", boundary_sf)

# 500 m (25 ha)
p_sr_500     <- plot_grid_metric(df_500_sf, "SR_interpolated", "SR", "eco", boundary_sf, plot_title = "25 ha")
p_t_gini_500 <- plot_grid_metric(df_500_sf, "delta_Gini_spatial_interpolated", "ΔGini spatial", "eco", boundary_sf)
p_size_g_500 <- plot_grid_metric(df_500_sf, "Gini_size_interpolated", "Gini size", "eco", boundary_sf)
p_npp_500    <- plot_grid_metric(df_500_sf, "NPP", bquote("NPP"), "prod", boundary_sf)

# ---------------------------
# 16. Assemble and export Fig. S6
# ---------------------------
# Arrange all maps into a 4 × 2 panel layout:
# rows = variables, columns = spatial scales.
combined_plot <- (p_sr_1000     | p_sr_500) /
                 (p_t_gini_1000 | p_t_gini_500) /
                 (p_size_g_1000 | p_size_g_500) /
                 (p_npp_1000    | p_npp_500) +
                 plot_layout(guides = "keep")

ggsave(output_pdf,
       plot = combined_plot,
       width = 22,
       height = 36,
       limitsize = FALSE,
       dpi = 300)

message("===== Fig. S6 saved to: ", output_pdf, " =====")
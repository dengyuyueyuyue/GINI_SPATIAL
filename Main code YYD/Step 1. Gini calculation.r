# ============================================================
# Demo workflow for calculating spatial and size structural diversity
# ============================================================
# This script demonstrates the workflow used to calculate
# delta_Gini_spatial, Gini_size, and their corresponding
# Lorenz asymmetry coefficients (S; denoted here as LAC).
#
# For confidentiality reasons, this demo uses only Grid 74
# (corresponding to plot/grid_id == 74 in the original dataset)
# as an example of the full analytical workflow.
#
# Spatial structural diversity is quantified from Voronoi area
# derived from observed and simulated individual coordinates.
# Size structural diversity is quantified from individual AGB.
#
# To avoid computational issues caused by duplicated point locations,
# a small random jitter is applied when necessary before Voronoi
# tessellation. Therefore, although a random seed is set, the
# spatial structural diversity metrics shown here may differ slightly
# from those reported in the full analysis.
#
# Input:
#   Grid74_renamed.xlsx
#
# Outputs:
#   1) Grid74_Voronoi_Area_Individual.xlsx
#   2) Grid74_Gini_LAC_Summary.xlsx
#
# Key metrics:
#   observed_area       = Voronoi area from observed coordinates
#   simulated_area      = Voronoi area from simulated coordinates
#   Gini_spatial_obs    = spatial inequality in observed area allocation
#   Gini_spatial_sim    = spatial inequality in simulated area allocation
#   delta_Gini_spatial  = observed minus simulated spatial inequality
#   Gini_size           = size inequality based on individual AGB
#   LAC                 = Lorenz asymmetry coefficient (S)
# ============================================================
rm(list = ls())
library(dplyr)
library(readxl)
library(writexl)
library(spatstat.geom)
library(spatstat.random)
library(deldir)

# ---------------------------
# 1. Core functions
# ---------------------------
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
  x <- sort(na.omit(x))
  n <- length(x)
  if (n < 2) return(NA_real_)
  if (any(x < 0)) stop("LAC requires non-negative values.")

  mu <- mean(x)
  Ln <- sum(x)
  if (Ln <= 0) return(NA_real_)
  if (all(abs(x - mu) <= tol)) return(1)

  m  <- sum(x < mu - tol)
  a  <- sum(abs(x - mu) <= tol)
  Lm <- if (m > 0) sum(x[1:m]) else 0

  if (a > 0) {
    L_m_plus_a <- sum(x[1:(m + a)])
    S_lower <- (m / n) + (Lm / Ln)
    S_upper <- ((m + a) / n) + (L_m_plus_a / Ln)
    return((S_lower + S_upper) / 2)
  } else {
    x_m <- x[m]
    x_m1 <- x[m + 1]
    d <- (mu - x_m) / (x_m1 - x_m)
    F_mu <- (m + d) / n
    L_mu <- (Lm + d * x_m1) / Ln
    return(F_mu + L_mu)
  }
}

# ---------------------------
# 2. Read and screen data
# ---------------------------
# The input file is a confidentiality-safe demo dataset for Grid 74.
path_grid74 <- "Grid74_renamed.xlsx"

dat <- read_excel(path_grid74) %>%
  mutate(source_row_id = row_number()) %>%
  filter(
    !is.na(sp_name),
    !is.na(agb_kg), agb_kg >= 0,
    !is.na(abs_x), !is.na(abs_y),
    !is.na(sim_x), !is.na(sim_y),
    abs_x >= 0, abs_x <= 25,
    abs_y >= 0, abs_y <= 25,
    sim_x >= 0, sim_x <= 25,
    sim_y >= 0, sim_y <= 25
  )

# ---------------------------
# 3. Voronoi tessellation
# ---------------------------
window <- owin(c(0, 25), c(0, 25))

make_ppp <- function(df, x_col, y_col, id_col) {
  ppp(df[[x_col]], df[[y_col]], window = window, marks = df[[id_col]])
}

jitter_if_needed <- function(ppp_obj, radius = 0.2) {
  xy <- data.frame(x = ppp_obj$x, y = ppp_obj$y)
  if (any(duplicated(xy))) {
    rjitter(ppp_obj, radius = radius, retry = TRUE)
  } else {
    ppp_obj
  }
}

extract_voronoi_area <- function(ppp_obj, area_name) {
  vor <- deldir(ppp_obj$x, ppp_obj$y, rw = c(0, 25, 0, 25), suppressMsge = TRUE)
  tiles <- tile.list(vor)

  areas <- sapply(tiles, function(tile) {
    area.owin(owin(poly = list(x = tile$x, y = tile$y)))
  })

  out <- data.frame(
    source_row_id = marks(ppp_obj),
    area = areas
  )
  names(out)[2] <- area_name
  out
}

set.seed(42)

ppp_obs <- dat %>%
  make_ppp("abs_x", "abs_y", "source_row_id") %>%
  jitter_if_needed(radius = 0.2)

ppp_sim <- dat %>%
  make_ppp("sim_x", "sim_y", "source_row_id") %>%
  jitter_if_needed(radius = 0.2)

area_obs <- extract_voronoi_area(ppp_obs, "observed_area")
area_sim <- extract_voronoi_area(ppp_sim, "simulated_area")

final_data <- dat %>%
  left_join(area_obs, by = "source_row_id") %>%
  left_join(area_sim, by = "source_row_id")

# ---------------------------
# 4. Grid-level summary
# ---------------------------
final_metrics_per_grid <- final_data %>%
  group_by(grid_id) %>%
  summarise(
    N_individuals = n(),
    N_species = n_distinct(sp_name),

    total_observed_area = sum(observed_area, na.rm = TRUE),
    total_simulated_area = sum(simulated_area, na.rm = TRUE),

    Gini_spatial_obs = calc_gini(observed_area),
    LAC_spatial_obs = calc_lac(observed_area),

    Gini_spatial_sim = calc_gini(simulated_area),
    LAC_spatial_sim = calc_lac(simulated_area),

    delta_Gini_spatial = Gini_spatial_obs - Gini_spatial_sim,

    Gini_size = calc_gini(agb_kg),
    LAC_size = calc_lac(agb_kg),

    .groups = "drop"
  )

# ---------------------------
# 5. Export outputs
# ---------------------------
write_xlsx(final_data, "Grid74_Voronoi_Area_Individual.xlsx")
write_xlsx(final_metrics_per_grid, "Grid74_Gini_LAC_Summary.xlsx")

print(final_metrics_per_grid)
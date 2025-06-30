#' Prepare Stan Data from Real Cougar Tracking Dataset
#'
#' @param df A data.frame containing GPS points with columns: x, y, animal_id, and any user-specified covariates
#' @param include_hab_cov Logical; whether to include habitat covariates
#' @param hab_cov_names Character vector of habitat covariate column names
#' @param include_ind_cov Logical; whether to include individual-level covariates
#' @param ind_cov_names Character vector of individual-level covariate column names
#' @param include_cell_effect Logical; whether to include spatial cell effect
#' @param grid_bounds Named list with xmin, xmax, ymin, ymax (only required if include_cell_effect = TRUE)
#' @param grid_res Grid resolution in meters (required if include_cell_effect = TRUE), please make sure that the grid bounds are divisible by the grid resolution
#' @return A named list ready to be passed to a CmdStan model
#' @export
prepare_stan_data <- function(df, include_hab_cov = TRUE, include_ind_cov = TRUE,
                              hab_cov_names = c("prey_avail", "hunter"),
                              ind_cov_names = c("age", "sex_m"),
                              include_cell_effect = TRUE,
                              grid_bounds = NULL,
                              grid_res = NULL) {
  library(dplyr)
  library(FNN)
  library(lubridate)
  library(tidyr)

  n_animals <- length(unique(df$animal_id))
  n_fixes <- df %>% count(animal_id) %>% pull(n) %>% max()

  # Fill in missing locations (assumes complete tracks or user has preprocessed)
  df <- df %>% arrange(animal_id, timestamp)
  df$fix_id <- ave(df$animal_id, df$animal_id, FUN = seq_along)

  # Time step calculation
  df <- df %>% group_by(animal_id) %>% mutate(
    lead_time = lead(timestamp),
    dt = as.numeric(difftime(lead_time, timestamp, units = "hours"))
  ) %>% ungroup()

  # Replace NA in last time step with median dt
  median_dt <- median(df$dt, na.rm = TRUE)
  df$dt[is.na(df$dt)] <- median_dt

  # Reshape to time_step matrix
  time_step_df <- df[, c("animal_id", "fix_id", "dt")] %>% pivot_wider(
    names_from = fix_id, values_from = dt)
  time_step <- as.matrix(time_step_df[,-1])
  time_step[is.na(time_step)] <- 0 # Replace NAs with 0

  n_locs <- df %>% count(animal_id) %>% pull(n)
  delta <- tapply(df$delta, df$animal_id, mean)

  # Cell effect
  if (include_cell_effect) {
    stopifnot(!is.null(grid_bounds), !is.null(grid_res))
    x_seq <- seq(grid_bounds$xmin + grid_res / 2, grid_bounds$xmax, by = grid_res)
    y_seq <- seq(grid_bounds$ymin + grid_res / 2, grid_bounds$ymax, by = grid_res)
    grid_centers <- as.matrix(expand.grid(x_seq, y_seq))
    knots_ce <- grid_centers

    # Assign cells
    cell_idx <- get.knnx(knots_ce, df[, c("x", "y")], k = 1)$nn.index[,1]
    df$cell_id <- cell_idx
    n_knots <- nrow(knots_ce)
  } else {
    knots_ce <- matrix(0, 1, 2)
    df$cell_id <- 1
    n_knots <- 1
  }

  cell_mat <- df[, c("animal_id", "fix_id", "cell_id")] %>% pivot_wider(values_from = cell_id, names_from = fix_id)

  cell_mat <- as.matrix(cell_mat[,-1])
  cell_mat[is.na(cell_mat)] <- 0  # Replace NAs with 0

  # Habitat covariates (3D array)
  num_hab_covs <- length(hab_cov_names)
  hab_cov <- array(NA, dim = c(n_animals, n_fixes, num_hab_covs))
  for (i in seq_along(hab_cov_names)) {
    temp <- df[, c("animal_id", "fix_id", hab_cov_names[i])] %>% pivot_wider(
      names_from = fix_id, values_from = hab_cov_names[i])
    hab_cov[,,i] <- as.matrix(temp[,-1])
  }
  hab_cov[is.na(hab_cov)] <- 0  # Replace NAs with 0

  # Individual covariates
  ind_data <- df %>% group_by(animal_id) %>% slice(1) %>% ungroup()
  z <- as.matrix(ind_data[, ind_cov_names])
  num_indv_covs <- ncol(z)

  list(
    n = n_animals,
    max_locs = n_fixes,
    n_locs = n_locs,
    time_step = time_step,
    delta = delta,
    n_knots = n_knots,
    knots_ce = knots_ce,
    sigma = 1.0,
    rho = 100.0,
    cell_mat = cell_mat,
    ind_cell_effect = as.integer(include_cell_effect),
    beta_prior = c(0, 1),
    llambda_prior = c(0, 2),
    alpha_prior = c(0, 1),
    num_hab_covs = num_hab_covs,
    hab_cov = hab_cov,
    num_indv_covs = num_indv_covs,
    z = z
  )
}

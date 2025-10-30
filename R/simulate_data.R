#' Simulate Animal Movement and Mortality Data
#'
#' @param n_animals Number of animals to simulate
#' @param n_fixes Maximum GPS fixes per animal
#' @param n_dead Number of dead animals (must be a squared number)
#' @param n_knots Number of spatial grid cells
#' @param landscape_size Size of landscape (in meters)
#' @return A list containing all inputs for the Stan model
#' @export
simulate_data <- function(n_animals = 200, n_fixes = 1000, n_dead = 40,
                                 n_knots = 100, landscape_size = 1000) {
  set.seed(123)

  # Create spatial grid
  grid_res <- sqrt(n_knots)
  length <- landscape_size / grid_res
  grid_x <- seq(length/2, landscape_size - length / 2, length.out = grid_res)
  grid_y <- seq(length/2, landscape_size - length / 2, length.out = grid_res)
  grid_centers <- as.matrix(expand.grid(grid_x, grid_y))
  knots_ce <- grid_centers

  # Create hunter density (centered peak) and prey availability (random)
  xx <- matrix(rep(grid_x, each = length(grid_y)), nrow = length(grid_y))
  yy <- matrix(rep(grid_y, times = length(grid_x)), nrow = length(grid_y))
  hunter <- 10*exp(-((xx - 750)^2 + (yy - 750)^2) / 1000)
  hunter[hunter<0.01] <- 0
  prey_avail <- matrix(rnorm(length(xx), mean = 500, sd = 50), nrow = nrow(xx))/100

  # Sample animal-level data
  dead_ids <- sample(n_animals, n_dead)
  age <- sample(1:20, n_animals, replace = TRUE)
  age[dead_ids] <- pmin(20, age[dead_ids] + sample(3:8, n_dead, replace = TRUE))
  sex <- sample(c("F", "M"), n_animals, replace = TRUE)
  sex_m <- as.integer(sex == "M")

  # Decide n_locs for each individual
  n_locs_dead <- sample(1:(n_fixes - 1), n_dead, replace = TRUE)
  n_locs_alive <- c(sample(1:n_fixes, n_animals - n_dead - 1, replace = TRUE), n_fixes)

  n_locs <- c()
  d <- 1
  a <- 1
  for(i in 1:n_animals) {
    if (i %in% dead_ids) {
      n_locs[i] <- n_locs_dead[d]
      d <- d + 1
    } else {
      n_locs[i] <- n_locs_alive[a]
      a <- a + 1
    }
  }
  # Prepare time step matrix
  time_step <- matrix(1.0, n_animals, n_fixes)
  for (i in 1:n_animals){
    # Set time steps for dead animals to 0 after their last location
    if (n_locs[i] == n_fixes) {
      next
    }
    time_step[i, c((n_locs[i] + 1) : n_fixes)] <- 0.0
  }

  # Helper for movement simulation
  simulate_track <- function(dead) {
    if (dead) {
      x <- rnorm(n_fixes, mean = 750, sd = 100)
      y <- rnorm(n_fixes, mean = 750, sd = 100)
    } else {
      x <- rnorm(n_fixes, mean = 250, sd = 100)
      y <- rnorm(n_fixes, mean = 250, sd = 100)
    }
    x <- pmin(pmax(x, 0), landscape_size)
    y <- pmin(pmax(y, 0), landscape_size)
    data.frame(x = x, y = y)
  }

  # Simulate data
  library(FNN)
  all_data <- list()
  for (i in 1:n_animals) {
    is_dead <- i %in% dead_ids
    coords <- simulate_track(is_dead)
    coords$animal_id <- i
    coords$fix_id <- 1:n_fixes
    coords$delta <- as.integer(is_dead)
    coords <- coords[1:n_locs[i],]
    coords$timestamp <- seq(1, n_locs[i], length.out = n_locs[i])
    all_data[[i]] <- coords
  }
  df <- do.call(rbind, all_data)


  # Assign cells
  cell_idx <- get.knnx(knots_ce, df[, c("x", "y")], k = 1)$nn.index[,1]
  df$cell_id <- cell_idx

  # Get raster indices for prey/hunter
  for (i in 1:nrow(df)) {
    xi <- max(1, min(grid_res, floor(df$x[i] / landscape_size * grid_res) + 1))
    yi <- max(1, min(grid_res, floor(df$y[i] / landscape_size * grid_res) + 1))
    df$prey_avail[i] <- prey_avail[yi, xi]
    df$hunter[i] <- hunter[yi, xi]
  }


  # Shape for Stan
  cell_mat <- reshape(df[, c("animal_id", "fix_id", "cell_id")],
                      timevar = "fix_id", idvar = "animal_id", direction = "wide")
  cell_mat <- as.matrix(cell_mat[,-1])
  cell_mat[is.na(cell_mat)] <- 0

  hab_cov <- array(NA, dim = c(n_animals, n_fixes, 2))
  for (i in 1:2) {
    vname <- c("prey_avail", "hunter")[i]
    temp <- reshape(df[, c("animal_id", "fix_id", vname)],
                    timevar = "fix_id", idvar = "animal_id", direction = "wide")
    hab_cov[,,i] <- as.matrix(temp[,-1])
  }

  hab_cov[is.na(hab_cov)] <- 0


  list(
    n = n_animals,
    max_locs = n_fixes,
    n_locs = n_locs,
    time_step = time_step,
    delta = tapply(df$delta, df$animal_id, mean),
    n_knots = n_knots,
    knots_ce = knots_ce,
    sigma = 1.0,
    rho = 100.0,
    cell_mat = cell_mat,
    ind_cell_effect = 1,
    beta_prior = c(0, 1),
    llambda_prior = c(0, 2),
    alpha_prior = c(0, 1),
    num_hab_covs = 2,
    hab_cov = hab_cov,
    num_indv_covs = 2,
    z = cbind(age, sex_m),
    raw_data = df
  )
}

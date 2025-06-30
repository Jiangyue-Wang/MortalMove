library(MortalMove)
library(sf)
library(dplyr)
# Load the example data
df <- readRDS("D:/OneDrive - UW/Research/OP_cougars/data/collar_dist2hs.rds")
summary <- readRDS("D:/OneDrive - UW/Research/OP_cougars/data/collar_meta_with_hs.rds")
df <- left_join(df, summary[,c("name", "death_human")], by = "name")
df$dist2hs <- df$dist2hs/1000
df$x <- st_coordinates(df$geometry)[, 1]
df$y <- st_coordinates(df$geometry)[, 2]
df <- st_drop_geometry(df)
df <- df %>% rename(
  animal_id = name,
  timestamp = time,
  delta = death_human
)
df$delta <- as.numeric(df$delta)
df$delta <- ifelse(df$delta==1, 0, 1)

# Prepare data

df_stan_data <- prepare_stan_data(
  df = df,
  include_hab_cov = TRUE, include_ind_cov = FALSE,
  hab_cov_names = c("dist2hs"),
  ind_cov_names = c(),
  include_cell_effect = TRUE,
  grid_bounds = list(xmin = 370000, xmax = 570000, ymin = 5110000, ymax = 5370000),
  grid_res = 20000
)

# Fit the model
set_cmdstan_path("C:/Users/jyuewang/.cmdstan/cmdstan-2.36.0/cmdstan-2.36.0")
library(cmdstanr)
model_file <- paste("D:/OneDrive - UW/Research/MortalMove/inst/stan", "mortality_model.stan", sep = "/")
mod <- cmdstan_model(model_file, exe_file = "D:/OneDrive - UW/Research/MortalMove/inst/stan/mortality_model.exe")
fit <- mod$sample(
  data = df_stan_data,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 1500,
  seed = 123
)
fit$summary()

fit$draws("llambda") |> bayesplot::mcmc_trace()
# Or for all parameters:
bayesplot::mcmc_trace(fit$draws())



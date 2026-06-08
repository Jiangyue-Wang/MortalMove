# Get movement tracking datasets ready for modeling, fit the model using cmdstanr or rstan

``` r

library(MortalMove)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Required data structure

You probably have two data frames: one with the GPS fixes with
timestamps and one with the individual information (e.g., age, sex,
fate, etc.). The first one should have the following columns:

- `animal_id`: unique identifier for each individual
- `timestamp`: date-time of the GPS fix (converted to
  [`lubridate::ymd_hms()`](https://lubridate.tidyverse.org/reference/ymd_hms.html)
  format, or numeric format in seconds)
- `x`: x-coordinate (in meters)
- `y`: y-coordinate (in meters)
- Other habitat covariates you want to include. Please note that if you
  want to use categorical variables (e.g., 4 season), you should input
  the variables as independent columns (e.g., `season_spring`,
  `season_summer`, `season_autumn`, with winter as reference level),
  where each column is a binary variable indicating the presence of that
  season. All values should be numeric.

Here, I included prey availability (`prey_avail`) and hunter density
(`hunter`) as habitat covariates.

``` r

head(df)
#>           x        y animal_id timestamp prey_avail   hunter
#> 1 495.44004 320.3625         1         1   5.179907 5.976312
#> 2  85.06545 219.6892         1         2   4.884911 4.535864
#> 3 329.58052 333.4154         1         3   5.230458 5.291962
#> 4 215.87068 260.6982         1         4   5.230458 5.291962
#> 5 109.47462 271.1017         1         5   4.884911 4.535864
#> 6 157.70392 252.8365         1         6   4.884911 4.535864
```

The second one should have the following columns: - `animal_id`: unique
identifier for each individual - `fate`: fate of the individual (e.g., 0
for alive, 1 for dead) - Other habitat covariates you want to include.
Please note that if you want to use categorical variables (e.g., sex),
you should input the variables as independent columns (e.g.,
\`sex_male’, with 1 for male, 0 for female). All values should be
numeric.

Here, I included age and sex as individual covariates.

``` r

head(indv)
#>   age sex_m animal_id fate
#> 1  16     0         1    0
#> 2  20     0         2    0
#> 3   6     0         3    0
#> 4  11     0         4    0
#> 5   8     1         5    0
#> 6   7     0         6    0
```

## Prepare data for modeling

First, we need to integrate these two tables into one data frame, which
will be used for modeling. The
[`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
function does this for you. Further, we need to change column name
`fate` to `delta`.

``` r

df <- left_join(df, indv, by = "animal_id")
df <- df %>% 
  dplyr::rename(delta = fate)
head(df)
#>           x        y animal_id timestamp prey_avail   hunter age sex_m delta
#> 1 495.44004 320.3625         1         1   5.179907 5.976312  16     0     0
#> 2  85.06545 219.6892         1         2   4.884911 4.535864  16     0     0
#> 3 329.58052 333.4154         1         3   5.230458 5.291962  16     0     0
#> 4 215.87068 260.6982         1         4   5.230458 5.291962  16     0     0
#> 5 109.47462 271.1017         1         5   4.884911 4.535864  16     0     0
#> 6 157.70392 252.8365         1         6   4.884911 4.535864  16     0     0
```

Next, we prepare the data for modeling. The `prepare_data_for_stan()`
function does this for you. You can specify which habitat covariates and
individual covariates you want to include in the model, as well as
whether you want to include a cell effect (i.e., a random effect for
each grid cell). You can also specify the grid bounds and resolution if
you want to include a cell effect. The grid bounds should be specified
as a list with `xmin`, `xmax`, `ymin`, and `ymax` values, and the grid
resolution should be specified in meters. Please note that the grid
bounds should cover the entire area where your data is collected, and
the difference between `xmax` and `xmin`, `ymax` and `ymin` should be
divisible by the grid resolution.

``` r

df_stan_data <- prepare_stan_data(
  df = df,
  include_hab_cov = TRUE, include_ind_cov = TRUE,
  hab_cov_names = c("hunter", "prey_avail"),
  ind_cov_names = c("age", "sex_m"),
  include_cell_effect = TRUE,
  grid_bounds = list(xmin = 0, xmax = 1000, ymin = 0, ymax = 1000),
  grid_res = 100
)
#> 
#> Attaching package: 'lubridate'
#> The following objects are masked from 'package:base':
#> 
#>     date, intersect, setdiff, union
names(df_stan_data)
#>  [1] "n"               "n_causes"        "max_locs"        "n_locs"         
#>  [5] "time_step"       "delta"           "n_knots"         "knots_ce"       
#>  [9] "sigma"           "rho"             "cell_mat"        "ind_cell_effect"
#> [13] "beta_prior"      "llambda_prior"   "alpha_prior"     "num_hab_covs"   
#> [17] "hab_cov"         "num_indv_covs"   "z"
```

## Fit model using cmdstanr

It is generally faster to fit the model using `cmdstanr` than `rstan`,
but requires a bit compilation before running. You can install
`cmdstanr` from [GitHub](https://github.com/stan-dev/cmdstanr), and then
set the path to the CmdStan installation directory using
[`set_cmdstan_path()`](https://mc-stan.org/cmdstanr/reference/set_cmdstan_path.html).
You can specify the number of chains, iterations, and other parameters
as needed.

``` r

library(cmdstanr)
set_cmdstan_path("YOUR_CMD_PATH")
model_file <- system.file("stan/mortality_model.stan", package = "MortalMove")
mod <- cmdstan_model(model_file, exe_file = system.file("stan/mortality_model.exe", package = "MortalMove"))
fit <- mod$sample(
  data = df_stan_data,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 50,
  iter_sampling = 150,
  seed = 123
)
```

## Fit model using rstan

If you prefer to use `rstan`, you can do so as well. It’s much simpler
to [install](https://mc-stan.org/install/), but may take longer to fit
than `cmdstanr`, especially for complicated models with spatial effects.

``` r

library(rstan)
model_file <- system.file("stan/mortality_model.stan", package = "MortalMove")
fit <- stan(
  file = model_file,
  data = df_stan_data,
  chains = 2,
  iter = 200,
  warmup = 50,
  seed = 123
)
```

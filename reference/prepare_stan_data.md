# Prepare Stan Data from Real Tracking Dataset

Prepare Stan Data from Real Tracking Dataset

## Usage

``` r
prepare_stan_data(
  df,
  include_hab_cov = TRUE,
  include_ind_cov = TRUE,
  hab_cov_names = c("prey_avail", "hunter"),
  ind_cov_names = c("age", "sex_m"),
  include_cell_effect = TRUE,
  grid_bounds = NULL,
  grid_res = NULL,
  cause_col = "delta",
  n_causes = NULL,
  max_time_gap = Inf
)
```

## Arguments

- df:

  A data.frame containing GPS points with columns: x, y, animal_id, and
  any user-specified covariates

- include_hab_cov:

  Logical; whether to include habitat covariates

- include_ind_cov:

  Logical; whether to include individual-level covariates

- hab_cov_names:

  Character vector of habitat covariate column names

- ind_cov_names:

  Character vector of individual-level covariate column names

- include_cell_effect:

  Logical; whether to include spatial cell effect

- grid_bounds:

  Named list with xmin, xmax, ymin, ymax (only required if
  include_cell_effect = TRUE)

- grid_res:

  Grid resolution in meters (required if include_cell_effect = TRUE),
  please make sure that the grid bounds are divisible by the grid
  resolution

- cause_col:

  Character; name of the column in df that contains cause of death
  information

- n_causes:

  Integer; number of unique causes of death (if applicable)

- max_time_gap:

  Numeric; maximum time interval, in hours, allowed between consecutive
  GPS fixes. Intervals larger than this value are truncated to this
  value. The default is Inf, meaning that intervals are not truncated.

## Value

A named list ready to be passed to a CmdStan model

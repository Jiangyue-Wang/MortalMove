# Perform Diagnostics on CmdStan Fit

Perform Diagnostics on CmdStan Fit

## Usage

``` r
diagnostic_fit(fit, plot = TRUE, observed_y = NULL, delta = NULL)
```

## Arguments

- fit:

  A CmdStanMCMC object returned from cmdstanr::sample()

- plot:

  Logical; whether to show diagnostic plots

- observed_y:

  Optional vector of observed values if predictive p-value is to be
  computed, the survey/survival length of each individual

- delta:

  Optional vector of censoring indicators (0 = censored, 1 =
  dead/failed) if predictive p-value is to be computed

## Value

Invisible list of diagnostics

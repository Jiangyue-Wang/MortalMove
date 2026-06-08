# Simulate Animal Movement and Mortality Data

Simulate Animal Movement and Mortality Data

## Usage

``` r
simulate_data(
  n_animals = 200,
  n_fixes = 1000,
  n_dead = 40,
  n_knots = 100,
  landscape_size = 1000
)
```

## Arguments

- n_animals:

  Number of animals to simulate

- n_fixes:

  Maximum GPS fixes per animal

- n_dead:

  Number of dead animals (must be a squared number)

- n_knots:

  Number of spatial grid cells

- landscape_size:

  Size of landscape (in meters)

## Value

A list containing all inputs for the Stan model

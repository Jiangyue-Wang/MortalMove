test_that("simulate data works", {
  sim <- simulate_data(n_animals = 20, n_fixes = 200, n_dead = 5, n_knots = 25)
  expect_s3_class(sim, "list")
})

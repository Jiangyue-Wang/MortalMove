#' Perform Diagnostics on CmdStan Fit
#'
#' @param fit A CmdStanMCMC object returned from cmdstanr::sample()
#' @param plot Logical; whether to show diagnostic plots
#' @param observed_y Optional vector of observed values if predictive p-value is to be computed, the survey/survival length of each individual
#' @param delta Optional vector of censoring indicators (0 = censored, 1 = dead/failed) if predictive p-value is to be computed
#' @return Invisible list of diagnostics
#' @export
diagnostic_fit <- function(fit, plot = TRUE, observed_y = NULL, delta = NULL) {
  library(posterior)
  library(bayesplot)
  library(loo)
  library(ggplot2)

  if(!inherits(fit, "CmdStanMCMC")) {
    stop("fit must be a CmdStanMCMC object from cmdstanr::sample()")
  }

  subset_draws <- c("llambda")

  # Check if some coefficients (alpha, beta) are included in the fit
  if ("alpha" %in% fit$metadata()$stan_variables) {
    subset_draws <- c("alpha", subset_draws)
  }
  if ("beta" %in% fit$metadata()$stan_variables) {
    subset_draws <- c("beta", subset_draws)
  }
  draws <- fit$draws(subset_draws)
  draws_summary <- summarise_draws(draws)
  print(draws_summary[,c("variable", "mean", "sd", "q5", "q95")])

  cat("\n--- Rhat and ESS Summary ---\n")
  print(draws_summary[, c("variable", "rhat", "ess_bulk", "ess_tail")])


  # LOO diagnostics
  if ("log_lik" %in% fit$metadata()$stan_variables) {
    log_lik_array <- fit$draws("log_lik")
    log_lik_matrix <- posterior::as_draws_matrix(log_lik_array)
    loo_result <- loo(log_lik_matrix)
    print(loo_result)
  } else {
    cat("\nNote: 'log_lik' not found in Stan output, skipping LOO diagnostics.\n")
  }

  # Posterior predictive p-values
  if (!is.null(observed_y)) {
    log_S <- posterior::as_draws_matrix(fit$draws(c("log_S")))
    failed_y <- observed_y[delta == 1]
    censored_y <- observed_y[delta == 0]

    expected_y <- log_S[,delta==1]
    expected_y <- failed_y/(-expected_y)

    expected_censored <- log_S[,delta==0]
    expected_censored <- censored_y/(-expected_censored)

    for(i in 1:nrow(expected_censored)) {
      expected_censored[i,] <- runif(length(censored_y), min = 0, max = min(max(observed_y),expected_censored[i,]))
    }

    p_failed <- colMeans(expected_y, na.rm = TRUE) < failed_y
    p_censored <- colMeans(expected_censored, na.rm = TRUE) < censored_y

    p_vals <- data.frame(
      variable = c("Failed", "Censored"),
      p_value = c(mean(p_failed, na.rm = TRUE), mean(p_censored, na.rm = TRUE))
    )
    cat("\nPosterior predictive p-value summary:\n")
    print(p_vals)

    if (plot) {

      ppc_dens_overlay(y = failed_y, yrep = expected_y) + ggtitle("Posterior Predictive Check: Failure times")
      ppc_dens_overlay(y = censored_y, yrep = expected_censored) + ggtitle("Posterior Predictive Check: Censored times")
    }
  } else {
    cat("\nNote: Observed y not available for posterior predictive checks.\n")
  }

  invisible(list(summary = draws_summary, loo = if (exists("loo_result")) loo_result else NULL,
                 p_values = if (exists("p_vals")) p_vals else NULL))
}

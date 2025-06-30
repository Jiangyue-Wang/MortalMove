#' Perform Diagnostics on CmdStan Fit
#'
#' @param fit A CmdStanMCMC object returned from cmdstanr::sample()
#' @param plot Logical; whether to show diagnostic plots
#' @param observed_y Optional vector of observed values if predictive p-value is to be computed
#' @return Invisible list of diagnostics
#' @export
diagnostic_fit <- function(fit, plot = TRUE, observed_y = NULL) {
  library(posterior)
  library(bayesplot)
  library(loo)

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
    expected_y <- matrix(data = NA, nrow = nrow(log_S), ncol = ncol(log_S))
    for(i in 1:nrow(log_S)){
      for (j in 1:ncol(log_S)) {
        expected_y[i,j] <- rbinom(1, 1, exp(log_S[i,j]))
      }
    }
    p_vals <- c()

    for(i in 1:nrow(log_S)) {
      p_vals <- c(p_vals, sum(expected_y[i,] <= observed_y)/length(observed_y))
      }


    cat("\nPosterior predictive p-value summary:\n")
    print(summary(p_vals))

    if (plot) {
      hist(p_vals, breaks = 20, main = "Posterior Predictive p-values",
           xlab = "p-value", col = "lightblue", border = "white", xlim = c(0, 1))
      abline(v = c(0.1, 0.9), col = c("red"), lty = 2)
    }
  } else {
    cat("\nNote: Observed y not available for posterior predictive checks.\n")
  }

  invisible(list(summary = draws_summary, loo = if (exists("loo_result")) loo_result else NULL,
                 p_values = if (exists("p_vals")) p_vals else NULL))
}

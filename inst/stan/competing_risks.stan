data{
  int<lower=1> n;
  int<lower=1> max_locs;
  array[n] int<lower=1> n_locs;
  matrix[n, max_locs] time_step;

  // Competing risks indicator:
  // delta[i,k] = 1 if individual i dies of cause k at its last location;
  // row is all 0 if censored.
  int<lower=1> n_causes;
  matrix<lower=0, upper=1>[n, n_causes] delta;

  // Spatial random effect
  int<lower=1> n_knots;
  array[n_knots] vector[2] knots_ce;
  array[n, max_locs] int cell_mat;
  real<lower=0> sigma;
  real<lower=0> rho;
  int ind_cell_effect;

  // Priors (shared across causes)
  vector[2] beta_prior;
  vector[2] llambda_prior;
  vector[2] alpha_prior;

  // Covariates
  int<lower=0> num_hab_covs;
  array[n, max_locs] vector[num_hab_covs] hab_cov;
  int<lower=0> num_indv_covs;
  matrix[n, num_indv_covs] z;
}

transformed data {
  matrix[n_knots, n_knots] L_cov;
  matrix[n_knots, n_knots] cov;
  if(ind_cell_effect == 1){
    for (i in 1:n_knots) {
      for (j in 1:n_knots) {
        cov[i, j] = square(sigma) * exp(-0.5 * squared_distance(knots_ce[i], knots_ce[j]) / square(rho));
      }
    }
    cov += diag_matrix(rep_vector(1e-10, n_knots));
    L_cov = cholesky_decompose(cov);
  }
}

parameters{
  // Cause-specific baselines and covariate effects
  vector[n_causes] llambda;
  matrix[num_hab_covs, n_causes] alpha;
  matrix[num_indv_covs, n_causes] beta;

  // Spatial effect (shared across causes)
  vector[n_knots*ind_cell_effect] eta;

  // --- If you want cause-specific spatial surfaces instead, use:
  // matrix[n_knots*ind_cell_effect, n_causes] eta;
}

transformed parameters{
  matrix[n, n_causes] nonspat;
  vector[n_knots] cell_effect;
  matrix[n, n_causes] log_h;
  matrix[n, n_causes] log_S;
  vector[n] log_lik;

  vector[max_locs] haz; // temporary per-cause hazard along the track
  real hab_spat;

  if(ind_cell_effect == 1){
    cell_effect = L_cov * eta;
    cell_effect = cell_effect - mean(cell_effect);
  }

  // nonspat[i,k] = llambda[k] + z_i * beta[,k]
  for(k in 1:n_causes){
    if(num_indv_covs > 0){
      nonspat[, k] = llambda[k] + z * beta[, k];
    } else {
      nonspat[, k] = rep_vector(llambda[k], n);
    }
  }

  // Per-cause cumulative hazards and event hazards
  for(k in 1:n_causes){
    for(i in 1:n){
      haz = rep_vector(0, max_locs);

      for(j in 1:n_locs[i]){
        if(num_hab_covs > 0){
          hab_spat = dot_product(alpha[, k], hab_cov[i, j]);
        } else {
          hab_spat = 0;
        }

        if(ind_cell_effect == 1){
          haz[j] = exp(nonspat[i, k] + hab_spat + cell_effect[cell_mat[i, j]]);
        } else {
          haz[j] = exp(nonspat[i, k] + hab_spat);
        }
      }

      // log_S[i,k] = -\int h_{ik}(t) dt (piecewise-constant over steps)
      log_S[i, k] = -dot_product(haz[1:n_locs[i]], time_step[i, 1:n_locs[i]]);
      log_h[i, k] = log(haz[n_locs[i]]);
    }
  }

  // Competing risks likelihood:
  // log p_i = [log h_{i,cause}(t_i)] + sum_k log S_{ik}(t_i)
  // where sum_k log S_{ik} = -\int \sum_k h_{ik}(t) dt
  for(i in 1:n){
    log_lik[i] = dot_product(log_h[i, ], delta[i, ]) + sum(log_S[i, ]);
  }
}

model{
  // Priors (same prior hyperparameters reused for each cause)
  llambda ~ normal(llambda_prior[1], llambda_prior[2]);
  if(num_hab_covs > 0){
    to_vector(alpha) ~ normal(alpha_prior[1], alpha_prior[2]);
  }
  if(num_indv_covs > 0){
    to_vector(beta) ~ normal(beta_prior[1], beta_prior[2]);
  }

  if(ind_cell_effect == 1){
    eta ~ std_normal();
  }

  target += sum(log_lik);
}


data{

		int<lower=1> n;
		int<lower=1> max_locs;
		array[n] int<lower=1> n_locs;
		matrix[n, max_locs] time_step;
		vector<lower=0, upper=1>[n] delta;
		int<lower=1> n_knots; // Number of grids in assessing spatial autocorrelation
		array[n_knots] vector[2] knots_ce; // grid cell center coordinates (in utm)
		array[n, max_locs] int cell_mat; // location -> cell ID mapping
		real<lower=0> sigma; // marginal sd of spatial effect
		real<lower=0> rho; // range parameter
		int ind_cell_effect;
		vector[2] beta_prior;
		vector[2] llambda_prior;
		vector[2] alpha_prior;
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
    cov += diag_matrix(rep_vector(1e-10, n_knots));  // Add jitter
    L_cov = cholesky_decompose(cov);
  }

  }



		parameters{

		real llambda;
		vector[num_hab_covs] alpha;
		vector[num_indv_covs] beta;
    vector[n_knots*ind_cell_effect] eta;

}

		transformed parameters{


		  vector[n] nonspat;
		  real hab_spat;
		  vector[max_locs] haz;
		  vector[n] log_h;
		  vector[n] log_S;
		  vector[n] log_lik;
      vector[n_knots] cell_effect;
      if(ind_cell_effect == 1){
        cell_effect = L_cov * eta;
		    cell_effect = cell_effect - mean(cell_effect);
      }


	  if(num_indv_covs > 0){
		  nonspat = llambda + z * beta;
	  }else{
	    nonspat = rep_vector(llambda, n);
	  }



		for(i in 1:n){

		  haz = rep_vector(0, max_locs);

		  for(j in 1:n_locs[i]){

		    if(num_hab_covs > 0){
		      hab_spat = dot_product(alpha, hab_cov[i,j]);
		    }else{
		      hab_spat = 0;
		    }

        if(ind_cell_effect == 1){
          haz[j] = exp(nonspat[i] + hab_spat + cell_effect[cell_mat[i,j]]);
        }else{
          haz[j] = exp(nonspat[i] + hab_spat);
        }


		  }

      log_S[i] = -dot_product(haz[1:n_locs[i]], time_step[i, 1:n_locs[i]]);
			log_h[i] = log(haz[n_locs[i]]);
		}

		log_lik = delta.*(log_h+log_S)	+ (1-delta).*log_S;
}

		model{

		llambda ~ normal(llambda_prior[1], llambda_prior[2]);
		if(num_hab_covs > 0){
				alpha ~ normal(alpha_prior[1], alpha_prior[2]);
		}
		if(num_indv_covs > 0){
		    beta ~ normal(beta_prior[1], beta_prior[2]);
		}
    if(ind_cell_effect == 1){
        eta ~ std_normal();
    }

    target += sum(log_lik);

}


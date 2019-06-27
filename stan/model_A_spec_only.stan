// -------------------------------------------------------------------------
// ------ Nonlinear hierarchical model for radial sap flux profiles --------
// ------ Model C: including random tree effects                    --------
// ------          and fixed effects for mu and k                   --------
// -------------------------------------------------------------------------
data{
    int<lower=0> I;                 // number of observations
    vector<lower=0>[I] flux;        // daily averages of sap flux rates
    vector<lower=0>[I] depth;       // sensor depths
    int<lower=0> J;                 // number of tree:date combinations
    int<lower=0> K;                 // number of trees
    int<lower=0> L;                 // number of species
    int<lower=1, upper=I> treedate[I]; // date ID
    int<lower=1, upper=I> tree[I];  // tree ID
    int<lower=1, upper=I> spec[I];  // species ID
  }

transformed data{
  // scaled response variable
  real<lower=0> sdflux;
  vector<lower=0>[I] scaleflux;
  // calculate scaled response variable
  sdflux = sd(flux);
  scaleflux = flux/sdflux;
}

parameters{
  // regression parameter vectors
  real mu0;       // mean RSFP depth
  real k0;        // concentration parameter
  real mult0_sc;  // average of multiplicator parameter (for scaled response)
  // variance and variance covariates
  real<lower=0> sigma_sc; // residual SD
  // hierarchical variance parameters
  // (multivariate correlations between effects are modeled after cholesky decomposition
      // of correlation matrix)
  vector[J] u_raw;
  vector[K] v_raw;
  matrix[2, L] z_spec;
  cholesky_factor_corr[2] L_Omega_spec;   // cholesky factor for species effects
  real<lower=0> tau_u;      // standard deviation for day effects
  real<lower=0> tau_v;      // tree level standard deviation
  vector<lower=0>[2] tau_w; // species level standard deviations
  // hierarchical scale parameter
  real<lower=0> hier_scale;
}

transformed parameters{
  // rescaled raw random effects
  matrix[2, L] z_spec_1;
  vector[J] u;
  vector[K] v;
  // -- variables from cholesky decomposition --
  matrix[L, 2] w;            // random species effects
  // cauchy scales from uniform transformed scales --
  // -- definitions of temporary variables --
  vector<lower=0>[I] exp_profile; 
  vector<lower=1e-10, upper=1-1e-10>[I] mu;
  vector<lower=1e-10>[I] k;
  vector<lower=0>[I] mult;
  
 // recenter and rescale raw random effects to aid convergence
  z_spec_1[1, ] = (z_spec[1, ] - mean(z_spec[1, ])) / sd(z_spec[1, ]);
  z_spec_1[2, ] = (z_spec[2, ] - mean(z_spec[2, ])) / sd(z_spec[2, ]);

  // calculate random effects from cholesky decomposition
  u = ((u_raw - mean(u_raw)) / sd(u_raw))* tau_u;
  v = ((v_raw - mean(v_raw)) / sd(v_raw))* tau_v;
  w = (diag_pre_multiply(tau_w, L_Omega_spec) * z_spec_1)';

  // regression equations for parameters
  mu   = inv_logit(mu0 +  w[spec, 1]);
  k    =       exp(k0  +  w[spec, 2]);
  mult =       exp(mult0_sc + u[treedate] + v[tree]);   
  
  // compute expected value of profile (from reparameterized beta pdf)
    for (i in 1:I){  exp_profile[i] = mult[i] * exp(beta_lpdf(depth[i] | mu[i] * k[i], (1 - mu[i]) * k[i]));
    }
}

model{
  // -- model --
  scaleflux ~ normal(exp_profile, sigma_sc);
  
  // -- priors on parameters --
  mu0      ~ cauchy(0, 10); // 
  k0       ~ cauchy(0, 10); //
  mult0_sc ~ cauchy(0, 10); //

  // -- priors on variance components --
  // -- residual SD --
  sigma_sc ~ cauchy(0, 2);

  // --- priors on random effects --
  // -- scaled random effects --
  u_raw             ~ normal(0, 1);
  v_raw             ~ normal(0, 1);
  to_vector(z_spec) ~ normal(0, 1);
  // -- cholesky factors of RE correlation matrices --
  L_Omega_spec ~ lkj_corr_cholesky(2.5); //
  // -- SDs of random effects --
  tau_u ~ cauchy(0, hier_scale);
  tau_v ~ cauchy(0, hier_scale);
  tau_w ~ cauchy(0, hier_scale);
  // hierarchical cauchy scale parameter
  hier_scale ~ uniform(0, 10);
}

generated quantities{
  // log likelihood 
  vector[I] log_lik;
  // parameters and parameter vectors
  real<lower=0> mult0;
  // correlation matrices
  corr_matrix[2] Omega_spec;
  // retransformed residual SD 
  real<lower=0> sigma;
  // variance decomposition for pseudo r-squared
  vector[3] y_vardecomp;

  // log likelihood
    for ( i in 1:I)  log_lik[i] = normal_lpdf(scaleflux[i] | exp_profile[i], sigma_sc);
  // transformed parameters and parameter vectors
  mult0      = mult0_sc + log(sdflux);
  // correlation matrices
  Omega_spec = multiply_lower_tri_self_transpose(L_Omega_spec);
  // variance parameters on untransformed scale
  sigma = sigma_sc * sdflux;
  // pseudo rsquared 
  { // tree level predictions (locally defined to save memory)
    vector[I] exp_profile_tree; 
    exp_profile_tree = (exp_profile ./ exp(u[treedate]));
    y_vardecomp  = [variance(exp_profile_tree), variance(exp_profile - exp_profile_tree), variance(scaleflux - exp_profile)]';
    y_vardecomp = y_vardecomp/sum(y_vardecomp);  
  } 
}

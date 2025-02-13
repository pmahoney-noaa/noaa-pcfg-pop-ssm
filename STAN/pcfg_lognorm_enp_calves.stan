// Kalman filter type state-space model (linear with guassian errors in log-space)
//  - Monte Carlo projections are in the generated quantities block
//  - Model fit to PCFG abundance estimates during 2002-2022 (Harris et al. 2024)
//  - Lambda [t] is conditioned on estimates of ENP calves from year t
//  - Note: The model fits in this version ignore the covariance between abundance estimates 
// Authors: John Brandon & Peter J. Mahoney
// Oct '24

data {
  int<lower=0> n_dat_yrs;              // number of abundance estimates to model
  int<lower=0> n_proj_yrs;             // number of projection years
  vector[n_dat_yrs] mu_logN_hat;       // mean of abundance estimates in log-space
  vector[n_dat_yrs] sigma_logN_hat;    // sd of abundance estimates in log-space
  vector[n_dat_yrs] mu_logC_hat;       // mean of ENP calf estimates in log-space
  vector[n_dat_yrs] sigma_logC_hat;    // sd of ENP calf estimates in log-space
  vector[n_proj_yrs] mu_logC_proj;     // mean of ENP calf estimates for projected years in log-space
  vector[n_proj_yrs] sigma_logC_proj;  // sd of ENP calf estimates for projected years in log-space
  vector[n_proj_yrs] N_harvest;        // number of known harvested animals during projected years
}

parameters {  // --------------------------------------------------------------------
  real<lower = 0> logN_init;
  vector[n_dat_yrs] logLambda;
  vector[n_dat_yrs] logC;
  real beta;
  real mu_logLambda;
  real<lower = 0> sigma_logLambda;
}

transformed parameters{  // ----------------------------------------------------
  //vector<lower = 0>[n_dat_yrs] logN;
  vector[n_dat_yrs] logN;      //removed the constraint on logN due to model fitting errors.
  logN[1] = logN_init;         // Initial abundance treated as parameter (but fit to data in model block)
  for(t in 2:(n_dat_yrs)){
    logN[t] = logN[t - 1] + logLambda[t - 1];
  }
}

model {  // --------------------------------------------------------------------
  // Priors
  mu_logLambda ~ normal(0, 1);        // hyper-prior for mean lambda
  sigma_logLambda ~ lognormal(0, 1);  // hyper-prior for sigma
  beta ~ normal(0, 1);                // Priors for beta coefs
  
  // Process error (lambda subsumes births, deaths, immigration, and emmigration)
  logLambda ~ normal(mu_logLambda + logC * beta, sigma_logLambda);
  
  // Likelihood of observations (observation error)
  logN ~ normal(mu_logN_hat, sigma_logN_hat);  // PCFG abundance
  logC ~ normal(mu_logC_hat, sigma_logC_hat);  // ENP calf abundance based on SWFSC modle
}

generated quantities{  // ------------------------------------------------------
 vector[n_dat_yrs] log_lik;     // save pointwise log-likelihood
 vector[n_proj_yrs] logN_proj;  // projected N in log-space some years into the future
 vector[n_proj_yrs] logC_proj;  // ENP calf estimates for projected years in log-space
 
 // following could be set up as an array, but stuck with vectors for simplicity in aggregating across sims.
 vector[n_dat_yrs - 1] logN_pyear1; // year 1 projected N in log-space post data year (starting in year t_start + 1); used in determining model performance
 vector[n_dat_yrs - 1] logN_pyear2; // year 2 projected N in log-space post data year (starting in year t_start + 2); used in determing model performance
 //vector[n_dat_yrs - 1] logN_pyear3; // year 3 projected N in log-space post data year (starting in year t_start + 3); used in determing model performance

 // save pointwise log-likelihood (e.g., use for leave-one-out [LOO] cross-validation) from model fits to abundance estimates
 for(t in 1:n_dat_yrs){
   log_lik[t] = normal_lpdf(logN[t] | mu_logN_hat[t], sigma_logN_hat[t]);  
 }
 
 // estimate ENP calf abundance in log-space during projected years
 for (t_proj in 1:n_proj_yrs) {
   logC_proj[t_proj] = normal_rng(mu_logC_proj[t_proj], sigma_logC_proj[t_proj]);
 }
 
 // project abundance in log-space
 logN_proj[1] = log(exp(normal_rng(logN[n_dat_yrs] + (mu_logLambda + logC[n_dat_yrs] * beta), sigma_logLambda)) - N_harvest[1]);  
 for(t_proj in 2:n_proj_yrs){
   logN_proj[t_proj] = log(exp(normal_rng(logN_proj[t_proj - 1] + (mu_logLambda + logC_proj[t_proj - 1] * beta), sigma_logLambda)) - N_harvest[t_proj]);
 }
 
 // project abundance in log-space for data year t; used in determing model performance
 for (t in 1:(n_dat_yrs - 1)){
   logN_pyear1[t] = normal_rng(logN[t] + (mu_logLambda + logC[t] * beta), sigma_logLambda);
   logN_pyear2[t] = normal_rng(logN_pyear1[t] + (mu_logLambda + logC[t + 1] * beta), sigma_logLambda);
   //logN_pyear3[t] = normal_rng(logN_pyear2[t] + mu_logLambda + logC[t + 2] * beta, sigma_logLambda);
 }
}

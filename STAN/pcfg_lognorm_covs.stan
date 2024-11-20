// Kalman filter type state-space model (linear with guassian errors in log-space)
//  - Monte Carlo projections are in the generated quantities block
//  - Model fit to PCFG abundance estimates during 2002-2022 (Harris et al. 2024)
//  - Lambda [t] is conditioned on PCFG calves and/or ENP Strandings from year t
//  - Note: The model fits in this version ignore the covariance between abundance estimates 
// Authors: John Brandon & Peter J. Mahoney
// Oct '24

data {
  int<lower=0> n_dat_yrs;                      // number of abundance estimates to model
  int<lower=0> n_proj_yrs;                     // number of projection years
  int<lower = 1> n_betas;                      // number of betas
  vector[n_dat_yrs] mu_logN_hat;               // mean of abundance estimates in log-space
  vector[n_dat_yrs] sigma_logN_hat;            // sd of abundance estimates in log-space
  matrix[n_dat_yrs, n_betas] x_lambda_dat;     // matrix of predictors on lambda (abundance data years)
  matrix[n_proj_yrs, n_betas] x_lambda_proj;   // matrix of predictors on lambda (projected years)
  vector[n_proj_yrs] N_harvest;                // number of known harvested animals during projected years
}

parameters {  // --------------------------------------------------------------------
  real<lower = 0> logN_init;
  vector[n_dat_yrs] logLambda;
  vector[n_betas] beta;
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
  logLambda ~ normal(x_lambda_dat * beta + mu_logLambda, sigma_logLambda);
  
  // Likelihood of observations (observation error)
  logN ~ normal(mu_logN_hat, sigma_logN_hat);
}

generated quantities{  // ------------------------------------------------------
 vector[n_dat_yrs] log_lik;     // save pointwise log-likelihood
 vector[n_proj_yrs] logN_proj;  // projected N in log-space some years into the future
 
 // following could be set up as an array, but stuck with vectors for simplicity in aggregating across sims.
 vector[n_dat_yrs - 1] logN_pyear1; // year 1 projected N in log-space post data year (starting in year t_start + 1); used in determining model performance
 vector[n_dat_yrs - 1] logN_pyear2; // year 2 projected N in log-space post data year (starting in year t_start + 2); used in determing model performance
 //vector[n_dat_yrs - 1] logN_pyear3; // year 3 projected N in log-space post data year (starting in year t_start + 3); used in determing model performance

 // save pointwise log-likelihood (e.g., use for leave-one-out [LOO] cross-validation) from model fits to abundance estimates
 for(t in 1:n_dat_yrs){
   log_lik[t] = normal_lpdf(logN[t] | mu_logN_hat[t], sigma_logN_hat[t]);  
 }
 
 // project abundance in log-space
 logN_proj[1] = log(exp(normal_rng(logN[n_dat_yrs] + (x_lambda_proj[1] * beta + mu_logLambda), sigma_logLambda)) - N_harvest[1]);  
 for(t_proj in 2:n_proj_yrs){
   logN_proj[t_proj] = log(exp(normal_rng(logN_proj[t_proj - 1] + (x_lambda_proj[t_proj] * beta + mu_logLambda), sigma_logLambda)) - N_harvest[t_proj]);
 }
 
 // project abundance in log-space for data year t; used in determing model performance
 for (t in 1:(n_dat_yrs - 1)){
   logN_pyear1[t] = normal_rng(logN[t] + (x_lambda_dat[t] * beta + mu_logLambda), sigma_logLambda);
   logN_pyear2[t] = normal_rng(logN_pyear1[t] + (x_lambda_dat[t + 1] * beta + mu_logLambda), sigma_logLambda);
   //logN_pyear3[t] = normal_rng(logN_pyear2[t] + (x_lambda_dat[t + 2] * beta + mu_logLambda), sigma_logLambda);
 }
}

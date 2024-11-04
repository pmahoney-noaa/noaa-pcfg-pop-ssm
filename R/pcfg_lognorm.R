#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Run PCFG state-space projection model(s)
# 
# Notes:
# 
# Authors: 
#  John R. Brandon, PhD [john dot brandon at icf dot com]
#  Peter J. Mahoney, PhD [peter dot mahoney at noaa dot gov] 
# Date: Oct '24
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

library(pacman)
p_load(cmdstanr, tidyverse, here, tidybayes, brms)
#set_cmdstan_path("~/cmdstan")  # your `cmdstan` path may be different

# Initialize -------------------------------------------------------------------

# Data import
source(here("R", "read_abun_retro.R"))     #; Ndata; tail(Ndata)  # Check
source(here("R", "read_calf_enp_retro.R")) #; Cdata; tail(Cdata)  # Check
source(here("R", "read_strand_retro.R"))   #; Sdata; tail(Sdata)  # Check

# Abundance Estimates
Ndata_input = filter(Ndata, year >= 2002) #; tail(Ndata_input)  # Check

# SWFSC ENP calf estimates
Cdata_input = filter(Cdata, year >= 2002) #; tail(Cdata_input)  # Check
Cdata_2020 <- Cdata_input %>%             # Impute mean of shoulder years for missing 2020 data
  slice(1) %>%
  mutate_all(~ NA) %>%
  mutate(
    year = 2020, method = "Mean of 2019 & 2021",
    mean_log = Cdata_input %>% filter(year >= 2019 & year <= 2021) %>% pull(mean_log) %>% mean(),
    sd_log = Cdata_input %>% filter(year >= 2019 & year <= 2021) %>% pull(sd_log) %>% mean()
  )

Cdata_input <- Cdata_input %>%
  add_row(Cdata_2020) %>%
  arrange(year)

# Covariate data (point estimates only) on lambda; at present, 
# not the precise time series and should be viewed as a filler for model development
x_lambda = data.frame(
  N_pcfg_calves = c(4,5,3,0,3,2,1,4,6,12,14,17,12,7,4,4,2,4,5,2, # abundance data years
                   0,0,0), # projected years, real values exist.
  N_strandings = Sdata %>% filter(state == "All" & year >= 2002) %>% pull(N_strand) %>% c(., 40) # This last value is not real!
)
x_lambda[,1] = scale(x_lambda[,1])

# Input data -------------------------------------------------------------------
init_pcfg_data = list(
  n_dat_yrs = nrow(Ndata_input),  
  n_proj_yrs = 2,                      # number of years to project into the future
  n_betas = ncol(x_lambda),            # number of coefficients on lambda
  mu_logN_hat = Ndata_input$mean_log,  # estimated pop abundance, log space
  sigma_logN_hat = Ndata_input$sd_log,
  mu_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(mean_log),       # estimated ENP calf abundance, log space
  sigma_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(sd_log),
  mu_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(mean_log),  # estimated contemporary ENP calf abundance, log space
  sigma_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(sd_log),
  x_lambda_dat = as.matrix(x_lambda %>% slice(1:nrow(Ndata_input))),
  x_lambda_proj = as.matrix(x_lambda %>% slice(-c(1:nrow(Ndata_input)))),
  N_harvest = c(0,0) #,0)
)
# init_pcfg_data  # Check

# STAN model definitions -------------------------------------------------------
f_pcfg_lognorm <- f_pcfg_base <- here::here('STAN', 'pcfg_lognorm_base.stan')
# f_pcfg_lognorm <- f_pcfg_base <- here::here('STAN', 'pcfg_lognorm_ar1_v1.stan')
# f_pcfg_lognorm <- f_pcfg_base <- here::here('STAN', 'pcfg_lognorm_ar1_v2.stan')
# f_pcfg_lognorm <- f_pcfg_base <- here::here('STAN', 'pcfg_lognorm_ar1_v1_delta.stan')
# f_pcfg_lognorm <- f_pcfg_enp <- here::here('STAN', 'pcfg_lognorm_enp_calves.stan')
# f_pcfg_lognorm <- f_pcfg_covs <- here::here('STAN', 'pcfg_lognorm_covs.stan')


# MCMC -------------------------------------------------------------------------
init_mod_cmdstanr = cmdstanr::cmdstan_model(f_pcfg_lognorm)  # Compile the Stan code for fitting model 
mcmc_pcfg = init_mod_cmdstanr$sample(data = init_pcfg_data, 
                               output_dir = here("out"),
                               seed = 42,
                               chains = 3,
                               parallel_chains = 3,
                               #iter_warmup = 3000, iter_sampling = 6000,
                               adapt_delta = 0.99)

##
## Model/chain performance ------------------------------------------------------------
##

mcmc_pcfg$diagnostic_summary()
# mcmc_pcfg$summary()
# bayesplot::mcmc_trace(tidy_mcmc, pars = c("logN_init", "mu0_logLambda", "sigma_logLambda"), 
#                       facet_args = list(ncol = 1, strip.position = "left"))

# LOO --------------------------------------------------------------------------
# Can use this for cross-validation to compare predictive performance against alternative projection-models
# In order to use this method you must compute and save the pointwise log-likelihood 
#   in your Stan code:
# https://mc-stan.org/loo/articles/loo2-with-rstan.html
loo_pcfg = mcmc_pcfg$loo(cores = 4)

# Wrangle Posterior ------------------------------------------------------------
# Tidy draws
tidy_mcmc = tidy_draws(mcmc_pcfg) #; names(tidy_mcmc)


# Plot -------------------------------------------------------------------------
threshold_N = 192    # Threshold on abundance below which a hunt is closed
threshold_Nmin = 171 # Threshold on minimum abundance below which a hunt is closed

tidy_plot_traj(Ndata_input, tidy_mcmc, threshold_N, threshold_Nmin)
# ggsave(filename = here("img", "pcfg_lognorm_ar1_v2_Oct-19-24.png"),
#        width = 8, height = 5, units = "in", dpi = 300)

# Table for hunt closures in projected years 1 through 2 or 3 from data years (t_start + 1) to (t_end - 1)
N_eval_table <- pred_summary_tbl(Ndata_input, tidy_mcmc, threshold_N, threshold_Nmin)
tidy_plot_propBelowThresh(N_eval_table)
tidy_plot_retroPred(N_eval_table)

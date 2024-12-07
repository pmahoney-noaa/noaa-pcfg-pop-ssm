#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Model fitting and diagnostics
# 
# Notes:
# 
# Authors: 
#  John R. Brandon, PhD [john dot brandon at icf dot com]
#  Peter J. Mahoney, PhD [peter dot mahoney at noaa dot gov] 
# Date: Oct '24
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

# Load required packages -------------------------------------------------------

library(pacman)
p_load(cmdstanr, tidyverse, here, tidybayes, brms, bayesplot)
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

# Define harvest data
N_harvest = c(0,0,0)

# Input data sets---------------------------------------------------------------
init_pcfg_data_3yr = list(
  n_dat_yrs = nrow(Ndata_input),  
  n_proj_yrs = 3,                      # number of years to project into the future
  n_betas = ncol(x_lambda),            # number of coefficients on lambda
  mu_logN_hat = Ndata_input$mean_log,  # estimated pop abundance, log space
  sigma_logN_hat = Ndata_input$sd_log,
  mu_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(mean_log),       # estimated ENP calf abundance, log space
  sigma_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(sd_log),
  mu_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(mean_log),  # estimated contemporary ENP calf abundance, log space
  sigma_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(sd_log),
  x_lambda_dat = as.matrix(x_lambda %>% slice(1:nrow(Ndata_input))),
  x_lambda_proj = as.matrix(x_lambda %>% slice(-c(1:nrow(Ndata_input)))),
  N_harvest = N_harvest
)

init_pcfg_data_2yr = list(
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
  N_harvest = N_harvest[1:2]
)

init_pcfg_data_2yr_calves = list(
  n_dat_yrs = nrow(Ndata_input),  
  n_proj_yrs = 2,                      # number of years to project into the future
  n_betas = 1,            # number of coefficients on lambda
  mu_logN_hat = Ndata_input$mean_log,  # estimated pop abundance, log space
  sigma_logN_hat = Ndata_input$sd_log,
  mu_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(mean_log),       # estimated ENP calf abundance, log space
  sigma_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(sd_log),
  mu_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(mean_log),  # estimated contemporary ENP calf abundance, log space
  sigma_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(sd_log),
  x_lambda_dat = as.matrix(x_lambda %>% dplyr::select(N_pcfg_calves) %>% dplyr::slice(1:nrow(Ndata_input))),
  x_lambda_proj = as.matrix(x_lambda %>% dplyr::select(N_pcfg_calves) %>% slice(-c(1:nrow(Ndata_input)))),
  N_harvest = N_harvest[1:2]
)

init_pcfg_data_2yr_strands = list(
  n_dat_yrs = nrow(Ndata_input),  
  n_proj_yrs = 2,                      # number of years to project into the future
  n_betas = 1,            # number of coefficients on lambda
  mu_logN_hat = Ndata_input$mean_log,  # estimated pop abundance, log space
  sigma_logN_hat = Ndata_input$sd_log,
  mu_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(mean_log),       # estimated ENP calf abundance, log space
  sigma_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input)) %>% pull(sd_log),
  mu_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(mean_log),  # estimated contemporary ENP calf abundance, log space
  sigma_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input))) %>% pull(sd_log),
  x_lambda_dat = as.matrix(x_lambda %>% dplyr::select(N_strandings) %>% slice(1:nrow(Ndata_input))),
  x_lambda_proj = as.matrix(x_lambda %>% dplyr::select(N_strandings) %>% slice(-c(1:nrow(Ndata_input)))),
  N_harvest = N_harvest[1:2]
)

# STAN model definitions -------------------------------------------------------
# models <- list.files("./STAN", "*.stan$")
f_pcfg_base <- here::here('STAN', 'pcfg_lognorm_base.stan')
f_pcfg_ar1v1 <- here::here('STAN', 'pcfg_lognorm_ar1_v1.stan')
f_pcfg_ar1v2 <- here::here('STAN', 'pcfg_lognorm_ar1_v2.stan')
f_pcfg_enp <- here::here('STAN', 'pcfg_lognorm_enp_calves.stan')
f_pcfg_covs <- here::here('STAN', 'pcfg_lognorm_covs.stan')
model_names <- factor(c("Base", "AR1v1", "AR1v2", "ENP Calves", "Calves/Strandings", "Calves only", "Strandings only"),
                      levels = c("Base", "AR1v1", "AR1v2", "ENP Calves", "Calves/Strandings", "Calves only", "Strandings only"))

# Model file pointers
models <- list(f_pcfg_base, f_pcfg_ar1v1, f_pcfg_ar1v2, f_pcfg_enp, 
               f_pcfg_covs, f_pcfg_covs, f_pcfg_covs)

# Specify input data
init_data <- list(init_pcfg_data_3yr, init_pcfg_data_3yr, init_pcfg_data_3yr,
                  init_pcfg_data_2yr, init_pcfg_data_2yr, 
                  init_pcfg_data_2yr_calves, init_pcfg_data_2yr_strands)

# Compile models ---------------------------------------------------------------
cmodels <- purrr::map(models, cmdstanr::cmdstan_model, .progress = T)

# MCMC -------------------------------------------------------------------------
mfit <- purrr::map2(cmodels, init_data, \(x, i) x$sample(
  data = i,
  output_dir = here("out"),
  seed = 42,
  chains = 3,
  parallel_chains = 3,
  #iter_warmup = 6000, iter_sampling = 9000,
  adapt_delta = 0.99
))

# Diagnostic summary -----------------------------------------------------------
tab_diag_summ <- purrr::map_dfr(mfit, \(x) {
  i <- x$diagnostic_summary()
  do.call(cbind, i) %>% 
    as.data.frame() %>% 
    add_column(chain = 1:3, .before = 1)
}) %>%
  add_column(
    Model = rep(model_names, each = 3),
    .before = 1
  )
tab_diag_summ

# LOO --------------------------------------------------------------------------
# Can use this for cross-validation to compare predictive performance against alternative projection-models
# In order to use this method you must compute and save the pointwise log-likelihood 
#   in your Stan code:
# https://mc-stan.org/loo/articles/loo2-with-rstan.html

loo_pcfg <- purrr::map(mfit, \(x) x$loo(cores = 4))
purrr::map(loo_pcfg, plot)

tab_loo <- purrr::map_dfr(loo_pcfg, \(x) t(x$estimates)[1, ]) %>%
  add_column(
    Model = model_names,
    .before = 1
  ) %>% 
  arrange(looic) %>%
  mutate(
    deltaLooic = looic - min(looic)
  )
tab_loo


# Prepping data for plotting ---------------------------------------------------
# Tidy draws
tfit = purrr::map(mfit, tidy_draws, .progress = T)
np <- purrr::map(mfit, nuts_params, .progress = T)

threshold_N = 192    # Threshold on abundance below which a hunt is closed
threshold_Nmin = 171 # Threshold on minimum abundance below which a hunt is closed


# Population trajectories ------------------------------------------------------
tidy_plot_traj_multimodel(Ndata_input, tfit, model_names, threshold_N, threshold_Nmin, ncols = 1)


# Retrospective prediction -----------------------------------------------------
N_eval_table <- pred_summary_tbl_multimodel(Ndata_input, tfit, model_names, threshold_N, threshold_Nmin)

# Summary stats
N_eval_summary_proj1yr <- N_eval_table %>%
  filter(proj_set == "pyear1") %>%
  group_by(model) %>%
  summarize(
    mnRSS = mean(rss),
    mnPercentile_abundEstN = mean(percentile_abundEstN),
    mnProp_below_threshold = mean(prop_below_threshold),
    mnProp_below_minThreshold = mean(prop_below_minThreshold),
    Nclosures = sum(closure)
  ) %>%
  arrange(mnRSS)

N_eval_summary_proj2yr <- N_eval_table %>%
  filter(proj_set == "pyear2") %>%
  group_by(model) %>%
  summarize(
    mnRSS = mean(rss),
    mnPercentile_abundEstN = mean(percentile_abundEstN),
    mnProp_below_threshold = mean(prop_below_threshold),
    mnProp_below_minThreshold = mean(prop_below_minThreshold),
    Nclosures = sum(closure)
  ) %>%
  arrange(mnRSS)

N_eval_summary_proj1yr; N_eval_summary_proj2yr;

# Figure
N_eval_table <- N_eval_table %>%
  mutate(model = factor(model, 
                              levels = c("Base", "AR1v1", "AR1v2", "Calves/Strandings", 
                                         "Calves only", "Strandings only", "ENP Calves")))
tidy_plot_retroPred(N_eval_table)



# Projected estimated final abundance + 2 years --------------------------------
#lastYr <- paste0("logN[", nrow(Ndata_input), "]")
lastYr <- "logN_proj[2]"
postLastYr <- purrr::imap(tfit, ~.x[lastYr] %>% add_column(Model = model_names[.y]))
postLastYr <- do.call(rbind, postLastYr)

postLastYr_mns <- postLastYr %>%
  group_by(Model) %>%
  summarize(
    Nmean = exp(mean(`logN_proj[2]`))
  )

ggplot(postLastYr, aes(group = Model, fill = Model, color = Model)) +
  geom_density(aes(x = exp(`logN_proj[2]`)), alpha = 0.4) +
  geom_rug(aes(x = threshold_N), color = "black", linewidth = 1.5) +
  geom_rug(aes(x = threshold_Nmin), color = "red", linewidth = 1.5) +
  # vline_at(threshold_N, linetype = 2) +
  # vline_at(threshold_Nmin, color = "red", linetype = 2) +
  geom_vline(data = postLastYr_mns, aes(xintercept = Nmean, group = Model, color = Model),
             alpha = 1, linewidth = 0.8, linetype = 2) +
  scale_x_continuous(limits = c(0, 350)) +
  theme_bw()
#bayesplot::mcmc_dens(exp(mfit[[1]]$draws(lastYr)))



# Divergence diagnostics -------------------------------------------------------
# See https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html

# Trace plots
# https://mc-stan.org/bayesplot/
mfit[[1]]$cmdstan_summary()
post <- mfit[[1]]$draws()
bayesplot::mcmc_trace(post, pars = "sigma_logLambda", n_warmup = 300)

# AR1v2  
mcmc_pairs(tfit[[3]], np = np[[3]], pars = c("mu0_logLambda", "mu_logLambda", "sigma0_logLambda", "sigma_logLambda"),
           off_diag_args = list(size = 0.75))
mcmc_pairs(tfit[[3]], np = np[[3]], pars = c("mu0_logLambda", "sigma0_logLambda", "logLambda[1]"),
           off_diag_args = list(size = 0.75))

mcmc_scatter(tfit[[3]], np = np[[3]], pars = c("mu0_logLambda", "sigma0_logLambda"), 
  #transform = list(logLambda = "exp"), 
  size = 1
)
mcmc_scatter(tfit[[3]], np = np[[3]], pars = c("logN[1]", "sigma0_logLambda"), 
             #transform = list(logN = "exp"), 
             size = 1
)

mcmc_trace(mfit[[3]]$draws(), np = np[[3]], pars = "sigma0_logLambda") #, window = c(300,500))

# ENP Calves
mcmc_pairs(tfit[[4]], np = np[[4]], pars = c("mu_logLambda", "sigma_logLambda", "beta"),
           off_diag_args = list(size = 0.75))

mcmc_scatter(tfit[[4]], np = np[[4]], pars = c("mu_logLambda", "sigma_logLambda"), 
             #transform = list(logLambda = "exp"), 
             size = 1
)

mcmc_trace(mfit[[4]]$draws(), np = np[[4]], pars = "sigma_logLambda") #, window = c(300,500))

# Strandings
mcmc_pairs(tfit[[7]], np = np[[7]], pars = c("mu_logLambda", "sigma_logLambda", "beta[1]"),
           off_diag_args = list(size = 0.75))












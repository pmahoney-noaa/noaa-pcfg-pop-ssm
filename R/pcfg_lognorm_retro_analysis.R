#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: 
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
  arrange(year) %>%
  
  # adding filler row for ENP calf data base on last year's
  # should be adjusted to some modeled estimate
  add_row(Cdata_input %>% filter(year == max(year)))

# Covariate data (point estimates only) on lambda; at present, 
# not the precise time series and should be viewed as a filler for model development
x_lambda = data.frame(
  N_pcfg_calves = c(9,4,5,3,0,3,2,1,4,6,12,14,17,12,7,4,4,2,4,5,2,11, # abundance data years
                    13,10), # projected years, real values exist.
  N_strandings = Sdata %>% filter(state == "All" & year >= 2002) %>% pull(N_strand) %>% c(., 40, 40) # This last value is not real!
)
x_lambda[,1] = scale(x_lambda[,1])

# Define harvest data
N_harvest = c(0,0)

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



# Retrospective evaluation -----------------------------------------------------
# Years to retrospectively evaluate
Y_retro = c(2013:2021)

purrr::map(Y_retro, function(y) {
  Ndata_input_y <- Ndata_input %>% filter(year <= y)
  
  # Input data set
  init_pcfg_data_calves <- init_pcfg_data_strandings <- init_pcfg_data <- list(
    n_dat_yrs = nrow(Ndata_input_y),  
    n_proj_yrs = 2,                      # number of years to project into the future
    n_betas = ncol(x_lambda),            # number of coefficients on lambda
    mu_logN_hat = Ndata_input_y$mean_log,  # estimated pop abundance, log space
    sigma_logN_hat = Ndata_input_y$sd_log,
    mu_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input_y)) %>% pull(mean_log),       # estimated ENP calf abundance, log space
    sigma_logC_hat = Cdata_input %>% slice(1:nrow(Ndata_input_y)) %>% pull(sd_log),
    mu_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input_y))) %>% slice(1:2) %>% pull(mean_log),  # estimated contemporary ENP calf abundance, log space
    sigma_logC_proj = Cdata_input %>% slice(-c(1:nrow(Ndata_input_y))) %>% slice(1:2) %>% pull(sd_log),
    x_lambda_dat = as.matrix(x_lambda %>% slice(1:nrow(Ndata_input_y))),
    x_lambda_proj = as.matrix(x_lambda %>% slice(-c(1:nrow(Ndata_input_y))) %>% slice(1:2)),
    N_harvest = N_harvest[1:2]
  )
  
  # Subset full dataset for models with single covariates
  init_pcfg_data_calves$n_betas <- 1
  init_pcfg_data_calves$x_lambda_dat <- as.matrix(x_lambda %>% dplyr::select(N_pcfg_calves) %>% slice(1:nrow(Ndata_input_y)))
  init_pcfg_data_calves$x_lambda_proj <- as.matrix(x_lambda %>% dplyr::select(N_pcfg_calves) %>% slice(-c(1:nrow(Ndata_input_y))) %>% slice(1:2))
  init_pcfg_data_strandings$n_betas <- 1
  init_pcfg_data_strandings$x_lambda_dat <- as.matrix(x_lambda %>% dplyr::select(N_strandings) %>% slice(1:nrow(Ndata_input_y)))
  init_pcfg_data_strandings$x_lambda_proj <- as.matrix(x_lambda %>% dplyr::select(N_strandings) %>% slice(-c(1:nrow(Ndata_input_y))) %>% slice(1:2))
  
  init_data <- list(
    init_pcfg_data, init_pcfg_data, init_pcfg_data, init_pcfg_data, init_pcfg_data,
    init_pcfg_data_calves, init_pcfg_data_strandings
  )
  
  # Compile models (if they haven't been)
  cmodels <- purrr::map(models, cmdstanr::cmdstan_model, .progress = T)
  
  # MCMC 
  mfit <- purrr::map2(cmodels, init_data, \(x, i) x$sample(
    data = i,
    output_dir = here("out"),
    seed = 42,
    chains = 3,
    parallel_chains = 3,
    #iter_warmup = 6000, iter_sampling = 9000,
    adapt_delta = 0.99
  ))
  
  save(mfit, file = here("out", paste0("Harris_2025_retro_y", y, ".dat")))
})


# Prepping data for plotting ---------------------------------------------------

# Years to retrospectively evaluate
Y_retro = c(2013:2021)

# Closure thresholds
threshold_N = 192    # Threshold on abundance below which a hunt is closed
threshold_Nmin = 171 # Threshold on minimum abundance below which a hunt is closed

N_eval_retro <- purrr::map(Y_retro, function(y) {
  # Load model results
  load(file = here("out", paste0("Harris_2025_retro_y", y, ".dat")))
  
  # Tidy draws
  tfit = purrr::map(mfit, tidy_draws, .progress = T)
  
  # Build tables
  if (y == max(Y_retro)) {
    input_data <- Ndata_input %>% filter(year <= (y + 1))
  } else{
    input_data <- Ndata_input %>% filter(year <= (y + 2))
  }
  
  N_eval_table <- purrr::imap_dfr(tfit, ~ .x %>% 
                                    gather_draws(logN_proj[year]) %>% #, logN_pyear3[year]) %>%
                                    select(year, N = .value) %>% 
                                    mutate(
                                      proj_set = paste0(year, " yr"),
                                      year = year + y,
                                      N = exp(N) 
                                    ) %>%
                                    left_join(input_data %>% dplyr::select(year, abundEstN = N), by = "year") %>% 
                                    #group_by(proj_set, year) %>% 
                                    summarize(
                                      proj_set = proj_set[1],
                                      model = model_names[.y],
                                      abundEst = mean(abundEstN),
                                      meanN = mean(N),
                                      medianN = median(N),
                                      loN_ci = quantile(N, 0.025),
                                      hiN_ci = quantile(N, 0.975),
                                      percentile_20 = quantile(N, 0.2),
                                      percentile_80 = quantile(N, 0.8),
                                      rss = (meanN - abundEst)^2, #sum((N - abundEst)^2),
                                      closure_Nmin = percentile_20 < threshold_Nmin,
                                      closure_N = meanN < threshold_N,
                                      closure = closure_Nmin | closure_N,
                                      percentile_abundEstN = ecdf(N)(abundEst),
                                      prop_below_threshold = mean(N < threshold_N)
                                    ) %>%
                                    ungroup())
  
  # Build population trajectories
  # At present, doesn't plot the last projected year in data year + 1 (e.g., 2023 if last abundance estimate in 2022)
  pl <- tidy_plot_traj_multimodel(input_data, tfit, model_names, 
                            threshold_N, threshold_Nmin, ncols = 1,
                            ylims = c(0, 350), truncated_retro = T)
  
  
  return(list(summary = N_eval_table, plot = pl))
})

# Summary stats
N_eval_retro_summ <- purrr::imap_dfr(N_eval_retro, ~ .x$summary %>% mutate(start_year = .y))

# save output for use in Quarto docs
save(N_eval_retro_summ, file = here("out", paste0("TruncatedRetroSummary_2025", ".dat")))

(N_eval_summary_proj1yr <- N_eval_retro_summ %>%
  filter(proj_set == "1 yr" & year != 2023) %>%
  filter(year > 2014) %>%
  group_by(model) %>%
  summarize(
    mnRSS = mean(rss),
    mdRSS = median(rss),
    mnPercentile_abundEstN = mean(percentile_abundEstN),
    mdPercentile_abundEstN = median(percentile_abundEstN),
    mnProp_below_threshold = mean(prop_below_threshold),
    mdProp_below_threshold = median(prop_below_threshold),
    Nclosures = sum(closure),
    closure_years = paste(cur_data()$year[closure == T], collapse = ",")
  ) %>%
  arrange(mnRSS))

(N_eval_summary_proj2yr <- N_eval_retro_summ %>%
  filter(proj_set == "2 yr" & year != 2023) %>%
  filter(year > 2014) %>%
  group_by(model) %>%
  summarize(
    mnRSS = mean(rss),
    mdRSS = median(rss),
    mnPercentile_abundEstN = mean(percentile_abundEstN),
    mdPercentile_abundEstN = median(percentile_abundEstN),
    mnProp_below_threshold = mean(prop_below_threshold),
    mdProp_below_threshold = median(prop_below_threshold),
    Nclosures = sum(closure),
    closure_years = paste(cur_data()$year[closure == T], collapse = ",")
  ) %>%
  arrange(mnRSS))


# Figure
N_eval_table <- N_eval_retro_summ %>%
  filter(year >= 2007 & year != 2023) %>%
  mutate(model = factor(model, 
                        levels = c("Base", "AR1v1", "AR1v2", "Calves/Strandings", 
                                   "Calves only", "Strandings only", "ENP Calves")))
tidy_plot_retroPred(N_eval_table, ylims = c(0, 500), truncated_retro = T)


# Relative model performance x time series length

(N_eval_summary_startYear <- N_eval_retro_summ %>%
    filter(proj_set == "1 yr" & year != 2023 & start_year > 1) %>%
    group_by(model, start_year) %>%
    summarize(
      mnRSS = mean(rss),
      #mdRSS = median(rss),
      mnPercentile_abundEstN = mean(percentile_abundEstN),
      mdPercentile_abundEstN = median(percentile_abundEstN),
      mnProp_below_threshold = mean(prop_below_threshold),
      mdProp_below_threshold = median(prop_below_threshold),
      Nclosures = sum(closure),
      closure_years = paste(cur_data()$year[closure == T], collapse = ",")
    ))

Retro_eval <- N_eval_summary_startYear %>%
  ungroup() %>%
  mutate(
    scRSS = (mnRSS - mean(mnRSS)) / sd(mnRSS)
  )

ggplot(Retro_eval, aes(x = start_year, y = mnRSS, group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 5000), oob = scales::squish) +
  #scale_y_continuous(limits = c(-0.7, 3), oob = scales::squish) +
  theme_bw()

(N_eval_summary_startYear <- N_eval_retro_summ %>%
    filter(proj_set == "2 yr" & year != 2023 & start_year > 1) %>%
    group_by(model, start_year) %>%
    summarize(
      mnRSS = mean(rss),
      #mdRSS = median(rss),
      mnPercentile_abundEstN = mean(percentile_abundEstN),
      mdPercentile_abundEstN = median(percentile_abundEstN),
      mnProp_below_threshold = mean(prop_below_threshold),
      mdProp_below_threshold = median(prop_below_threshold),
      Nclosures = sum(closure),
      closure_years = paste(cur_data()$year[closure == T], collapse = ",")
    ))

Retro_eval <- N_eval_summary_startYear %>%
  ungroup() %>%
  mutate(
    scRSS = (mnRSS - mean(mnRSS)) / sd(mnRSS)
  )

ggplot(Retro_eval, aes(x = start_year, y = mnRSS, group = model, color = model)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(-1, 10000), oob = scales::squish) +
  #scale_y_continuous(limits = c(-0.4, 0), oob = scales::squish) +
  theme_bw()

# Diagnostic summary -----------------------------------------------------------


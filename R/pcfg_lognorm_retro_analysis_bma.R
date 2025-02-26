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

# Additional functions ---------------------------------------------------------

tidy_traj_multimodel = function(input_data, tidy_mcmc, model_names, 
                                     truncated_retro = F) { # arg. to specify a truncated retrospective analysis
  
  N_out = purrr::imap_dfr(tidy_mcmc, ~ .x %>% 
                            spread_draws(logN[year]) %>% 
                            select(year, logN) %>% 
                            ungroup() %>% 
                            mutate(
                              model = model_names[.y],
                              year = year + input_data$year[1] - 1))
  
  if (truncated_retro) {
    N_proj_out = purrr::imap_dfr(tidy_mcmc, ~ .x %>% 
                                   spread_draws(logN_proj[year]) %>% 
                                   select(year, logN = logN_proj) %>% 
                                   ungroup() %>% 
                                   mutate(
                                     model = model_names[.y],
                                     year = year + max(input_data$year) - 2))
  } else {
    N_proj_out = purrr::imap_dfr(tidy_mcmc, ~ .x %>% 
                                   spread_draws(logN_proj[year]) %>% 
                                   select(year, logN = logN_proj) %>% 
                                   ungroup() %>% 
                                   mutate(
                                     model = model_names[.y],
                                     year = year + max(input_data$year)))
  }
  
  N_out_table = N_out %>% 
    bind_rows(N_proj_out) %>% 
    group_by(model, year) %>% 
    summarize(mean = mean(logN),
              median = median(logN),
              lo_ci = quantile(logN, 0.025),
              hi_ci = quantile(logN, 0.975),
              percentile_20 = quantile(logN, 0.2),
              percentile_80 = quantile(logN, 0.8)) %>% 
    ungroup() %>% 
    mutate(across(.cols = c(-model, -year), .fns = exp)) %>%
    arrange(model, year)

  return(N_out_table)
}

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
f_pcfg_enp <- here::here('STAN', 'pcfg_lognorm_enp_calves.stan')
f_pcfg_covs <- here::here('STAN', 'pcfg_lognorm_covs.stan')
model_names <- factor(c("Base", "AR1v1", "ENP Calves", "Calves only", "Strandings only", "Calves/Strandings"),
                      levels = c("Base", "AR1v1", "ENP Calves", "Calves only", "Strandings only", "Calves/Strandings"),
                      labels = c("Base", "AR1", "ENP Calves","PCFG Calves only", "ENP Strandings only", "PCFG Calves + ENP Strandings"))

# Model file pointers
models <- list(f_pcfg_base, f_pcfg_ar1v1, #f_pcfg_ar1v2, 
               f_pcfg_enp, f_pcfg_covs, f_pcfg_covs, f_pcfg_covs)



# Retrospective evaluation -----------------------------------------------------
# Years to retrospectively evaluate
Y_retro = c(2017:2022)

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
    init_pcfg_data, init_pcfg_data, init_pcfg_data, #init_pcfg_data,
    init_pcfg_data_calves, init_pcfg_data_strandings, init_pcfg_data
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
  
  save(mfit, file = here("out", paste0("Harris_2025_retro_y", y, "_cbma.dat")))
})


# Prepping data for plotting ---------------------------------------------------

# Years to retrospectively evaluate
Y_retro = c(2017:2022)

# Closure thresholds
threshold_N = 192    # Threshold on abundance below which a hunt is closed
threshold_Nmin = 171 # Threshold on minimum abundance below which a hunt is closed

N_eval_retro_ma <- purrr::map(Y_retro, function(y) {
  
  # Input_data subset
  if (y == max(Y_retro)) {
    input_data <- Ndata_input %>% filter(year <= (y + 1))
  } else{
    input_data <- Ndata_input %>% filter(year <= (y + 2))
  }
  
  # Load model results
  load(file = here("out", paste0("Harris_2025_retro_y", y, "_bma.dat")))
  
  # LOO
  loo_pcfg <- purrr::map(mfit, \(x) x$loo(cores = 4, k_threshold = 0.7))
  
  # PseudoBMA+ weights
  #wgts <- loo_model_weights(loo_pcfg)
  wgts <- loo_model_weights(loo_pcfg, method = "pseudobma", BB = TRUE)
  
  # Tidy draws
  tfit = purrr::map(mfit, tidy_draws, .progress = T)
  
  # Model-specific estimates
  mo <- tidy_traj_multimodel(input_data %>% slice(-((n()-1):n())), tfit, model_names, truncated_retro = F)
  
  # Model averaged estimates
  ma <- tidy_model_avg(tfit, wgts, iters = 1000, input_data %>% slice(-((n()-1):n())), seed = 101)
  
  ma_sub <- ma %>%
    slice_tail(n = 2) %>%
    add_column(
      proj_set = c("1 yr", "2 yr"),
      model = "BMA",.before = 1
    ) %>%
    left_join(input_data %>% dplyr::select(year, abundEstN = N), by = "year") %>%
    mutate(
      rss = (mean - abundEstN)^2, #sum((N - abundEst)^2),
      closure_Nmin = percentile_20 < threshold_Nmin,
      closure_N = mean < threshold_N,
      closure = closure_Nmin | closure_N
    )

  # Build population trajectories
  # At present, doesn't plot the last projected year in data year + 1 (e.g., 2023 if last abundance estimate in 2022)

  # mma <- ma %>%
  #   add_column(model = "BMA", .before = 1) %>%
  #   add_row(mo)
    
  pl <- input_data %>% 
    ggplot(aes(x = year, y = N)) +  
    geom_ribbon(data = ma, aes(y = mean, ymin = percentile_20, ymax = percentile_80), fill = "red", alpha = 0.2) +
    #geom_ribbon(data = ma, aes(y = mean, ymin = lo_ci, ymax = hi_ci), fill = "red", alpha = 0.2) +
    
    geom_ribbon(data = mo, aes(y = mean, ymin = percentile_20, ymax = percentile_80), alpha = 0.2) +
    geom_ribbon(data = mo, aes(y = mean, ymin = lo_ci, ymax = hi_ci), alpha = 0.2) +
    geom_line(data = mo, aes(y = mean)) +
    geom_point(size = 3, color = "white", fill = "black", shape = 21) +
    scale_x_continuous(limits = c(2002, y+2), 
                       breaks = seq(2002, y+2, by = 1), minor_breaks = NULL) +
    scale_y_continuous(limits = c(0, 350), oob = scales::squish) +
    geom_hline(yintercept = threshold_N) + 
    geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
    geom_errorbar(aes(ymin = low_95CI, ymax = high_95CI)) +  
    geom_errorbar(aes(ymin = low_60CI, ymax = N), linewidth = 0) +
    geom_point(aes(y = low_60CI), shape = 23, fill = "red", size = 2) +
    theme_bw(base_size = 16) +
    facet_wrap(~ model, ncol = 1) +
    labs(x = "Year", y = "PCFG Abundance") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    NULL
  
  return(list(summary = ma_sub, plot = pl, weights = wgts))
})

# Summary stats
N_eval_retro_ma_weights <- purrr::imap_dfr(N_eval_retro_ma, ~ .x$weights)
N_eval_retro_summ_ma <- purrr::imap_dfr(N_eval_retro_ma, ~ .x$summary %>% mutate(start_year = .y))

# save output for use in Quarto docs
save(N_eval_retro_summ_ma, file = here("out", paste0("TruncatedRetroSummary_BMA_2025", ".dat")))

(N_eval_summary_proj1yr <- N_eval_retro_summ_ma %>%
  filter(proj_set == "1 yr" & year != 2023) %>%
  group_by(model) %>%
  summarize(
    mnRSS = mean(rss),
    mdRSS = median(rss),
    # mnPercentile_abundEstN = mean(percentile_abundEstN),
    # mdPercentile_abundEstN = median(percentile_abundEstN),
    # mnProp_below_threshold = mean(prop_below_threshold),
    # mdProp_below_threshold = median(prop_below_threshold),
    Nclosures = sum(closure),
    closure_years = paste(cur_data()$year[closure == T], collapse = ",")
  ) %>%
  arrange(mnRSS))

(N_eval_summary_proj2yr <- N_eval_retro_summ_ma %>%
  filter(proj_set == "2 yr" & year != 2023) %>%
  group_by(model) %>%
  summarize(
    mnRSS = mean(rss),
    mdRSS = median(rss),
    # mnPercentile_abundEstN = mean(percentile_abundEstN),
    # mdPercentile_abundEstN = median(percentile_abundEstN),
    # mnProp_below_threshold = mean(prop_below_threshold),
    # mdProp_below_threshold = median(prop_below_threshold),
    Nclosures = sum(closure),
    closure_years = paste(cur_data()$year[closure == T], collapse = ",")
  ) %>%
  arrange(mnRSS))


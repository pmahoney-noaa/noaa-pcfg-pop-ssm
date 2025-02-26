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
p_load(cmdstanr, tidyverse, here, tidybayes, brms, bayesplot, furrr)
plan(multisession, workers = 4)

# Initialize -------------------------------------------------------------------

# Data import
source(here("R", "read_abun_retro.R"))     #; Ndata; tail(Ndata)  # Check

# Abundance Estimates
# Range of declining lambda
peakYr <- Ndata %>% filter(N == max(N)) %>% pull(year)

(lowLambda <- Ndata %>% filter(year > peakYr & year <= 2018) %>% pull(lambda) %>% mean())
(hiLambda <- Ndata %>% filter(year > 2018 & year <= 2022) %>% pull(lambda) %>% mean())
Ndata %>% filter(year > 2004 & year <= 2011) %>% pull(SE) %>% mean()

# Define harvest data
N_harvest = c(0,0)

# STAN model definitions -------------------------------------------------------
# models <- list.files("./STAN", "*.stan$")
f_pcfg_base <- here::here('STAN', 'pcfg_lognorm_base.stan')
f_pcfg_ar1v1 <- here::here('STAN', 'pcfg_lognorm_ar1_v1.stan')
model_names <- factor(c("Base", "AR1v1"), 
                      levels = c("Base", "AR1v1"), 
                      labels = c("Base", "AR1"))

# Model file pointers
models <- list(f_pcfg_base, f_pcfg_ar1v1)



# Simulation evaluation -----------------------------------------------------
# Constant decline

# Params
i_lambdas <- seq(0.92, 0.98, 0.005)
n_years <- 15
n_start <- 250
se_start <- 10
n_sim <- 10

# Simulate data, stable 10 years, decline 10 years
base <- data.frame(N = rep(n_start, n_years), SE = rep(se_start, n_years))

Ndata_input_list <- purrr::map(i_lambdas, function(y) {
  for (x in 1:n_sim) {
    if (x == 1) {
      N1 = base %>% slice(n_years) %>% pull(N)
      Nsim = N1*y
    } else {
      Nsim[x] = Nsim[x-1]*y
    }
  }
  
  dy <- data.frame(N = Nsim, SE = se_start) %>%
    add_row(base, .before = 1) %>%
    mutate(
      mean_log = calc_mu_log(mu = N, sd = SE),
      sd_log = calc_sd_log(mu = N, sd = SE),
      sd_log = sd_log[1]
    )
  
  do <- list()
  do[[as.character(y)]] <- dy
  return(do)
})

# Run simulations
furrr::future_map(Ndata_input_list, function (i) {
  nam <- names(i)
  Ndata_input <- i[[1]]
  
  Y_retro = c(n_years:(n_years + n_sim - 1))
  sims_y <- purrr::map(Y_retro, function(y) {
    Ndata_input_y <- Ndata_input %>% slice(1:y)
    
    # Input data set
    init_pcfg_data <- list(
      n_dat_yrs = nrow(Ndata_input_y),  
      n_proj_yrs = 2,                        # number of years to project into the future
      mu_logN_hat = Ndata_input_y$mean_log,  # estimated pop abundance, log space
      sigma_logN_hat = Ndata_input_y$sd_log,
      N_harvest = N_harvest[1:2]
    )
    
    init_data <- list(
      init_pcfg_data, init_pcfg_data
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
    
    return(mfit)
  })
  
  save(sims_y, file = here("out", paste0("Harris_2025_simDecline_lambda_", nam, ".dat")))
})


# Evaluating simulations
# Closure thresholds
threshold_N = 192    # Threshold on abundance below which a hunt is closed
threshold_Nmin = 171 # Threshold on minimum abundance below which a hunt is closed


# Difference between actual and projections
summ_proj <- furrr::map_dfr(1:length(i_lambdas), function(x) {
  lambda <- i_lambdas[x]
  load(here("out", paste0("Harris_2025_simDecline_lambda_", lambda, ".dat")))
  
  input_data <- Ndata_input_list[[x]][[1]] %>%
    mutate(year = 1:(n_years + n_sim))
  
  out <- purrr::map_dfr(1:length(sims_y), function(y) {
    tfit <- purrr::map(sims_y[[y]], tidy_draws, .progress = T)
    

    tbl <- purrr::imap_dfr(tfit, ~ .x %>% 
                             spread_draws(logN_proj[year]) %>% #, logN_pyear3[year]) %>%
                             mutate(
                               proj_set = if_else(year == 1, "1 yr", "2 yr")
                             ) %>%
                             select(proj_set, year, logN = logN_proj) %>% 
                             mutate(
                               year = year + y + (n_years-1),
                               N = exp(logN) 
                             ) %>%
                             filter(year <= max(input_data$year)) %>%
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
                               percentile_abundEstN = ecdf(N)(abundEst)
                             ) %>% 
                             ungroup() %>%
                             mutate(
                               sim_lambda = lambda
                             )
    )
    
    return(tbl)
  })

  return(out)
})
summary(lm((meanN - abundEst) ~ model + sim_lambda + year + proj_set, data = summ_proj))
summary(lm((meanN - abundEst) ~ model * sim_lambda + year + proj_set, data = summ_proj))
save(summ_proj, file = here("out", "SimulatedDecline_projections.dat"))

summ_proj <- summ_proj %>%
  filter(sim_lambda %in% c(seq(0.92, 0.98, by = 0.01))) %>%
  mutate(
    lambda = factor(sim_lambda, 
                     levels = unique(sim_lambda), 
                     labels = paste("lambda:", unique(sim_lambda)))
  )
ggplot(summ_proj %>% filter(proj_set == "1 yr"), aes(x = year, group = model)) +
  geom_hline(yintercept = threshold_N) + 
  geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
  geom_ribbon(aes(ymin = percentile_20, ymax = percentile_80, fill = model), alpha = 0.2) +
  geom_line(aes(y = meanN, color = model), alpha = 0.5) +
  geom_point(aes(y = abundEst), color = 'black') +
  scale_fill_manual(values = c("red", "black")) +
  scale_color_manual(values = c("red", "black")) +
  facet_wrap(~ lambda, nrow = 4, labeller = label_parsed) +
  labs(x = "Year", "Abundance (N)") +s
  scale_x_continuous(limits = c(16, 25), breaks = 16:25) +
  theme_bw()

ggplot(summ_proj %>% filter(proj_set == "2 yr"), aes(x = year, group = model)) +
  geom_hline(yintercept = threshold_N) + 
  geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
  geom_ribbon(aes(ymin = percentile_20, ymax = percentile_80, fill = model), alpha = 0.2) +
  geom_line(aes(y = meanN, color = model), alpha = 0.5) +
  geom_point(aes(y = abundEst), color = 'black') +
  scale_fill_manual(values = c("red", "black")) +
  scale_color_manual(values = c("red", "black")) +
  facet_wrap(~ lambda, nrow = 4, labeller = label_parsed) +
  labs(x = "Year", "Abundance (N)") +
  scale_x_continuous(limits = c(16, 25), breaks = 16:25) +
  theme_bw()

# Difference from actual and predicted trend
summ_trend <- purrr::map(1:length(i_lambdas), function(x) {
  lambda <- i_lambdas[x]
  load(here("out", paste0("Harris_2025_simDecline_lambda_", lambda, ".dat")))
  
  input_data <- Ndata_input_list[[x]][[1]] %>%
    mutate(year = 1:(n_years + n_sim))
  
  tfit <- purrr::map(sims_y[[length(sims_y)]], tidy_draws, .progress = T)
  
  out <- purrr::map_dfr(1:length(tfit), function(y) {
    N_out = tfit[[y]] %>% 
      spread_draws(logN[year]) %>% 
      select(year, logN) %>% 
      ungroup() %>% 
      mutate(
        model = model_names[y],
        year = year + input_data$year[1] - 1,
        ) 
  }) %>% 
    group_by(model, year) %>% 
    summarize(mean = mean(logN),
              median = median(logN),
              lo_ci = quantile(logN, 0.025),
              hi_ci = quantile(logN, 0.975),
              percentile_20 = quantile(logN, 0.2),
              percentile_80 = quantile(logN, 0.8)) %>% 
    ungroup() %>% 
    mutate(across(.cols = c(-year, -model), .fns = exp)) 
  
  out <- out %>% filter(year >= n_years + 1)
  
  pl <- input_data %>% 
    filter(year >= n_years + 1) %>%
    ggplot(aes(x = year, y = N)) +  
    geom_ribbon(data = out, aes(y = mean, ymin = percentile_20, ymax = percentile_80, group = model, fill = model), alpha = 0.2) +
    #geom_ribbon(data = out, aes(y = mean, ymin = lo_ci, ymax = hi_ci, group = model, fill = model), alpha = 0.2) +
    geom_line(data = out, aes(y = mean, group = model, color = model)) +
    geom_point(size = 3, color = "white", fill = "black", shape = 21) +
    scale_fill_manual(values = c("red", "black")) +
    scale_color_manual(values = c("red", "black")) +
    scale_x_continuous(limits = c(min(out$year), max(out$year)), 
                       breaks = seq(min(out$year), max(out$year), by = 1), minor_breaks = NULL) +
    scale_y_continuous(limits = c(0, 350), breaks = seq(0, 350, by = 100)) +
    geom_hline(yintercept = threshold_N) + 
    geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
    theme_bw(base_size = 16) +
    labs(x = "Year", y = "Simulated Abundance") +
    # theme(
    #   axis.text.x = element_text(angle = 45, hjust = 1)
    # ) +
    NULL
  
  return(pl)
})

# Mean difference between actual and observed threshold cross (N and Nmin)
summ_diff <- summ_proj %>%
  group_by(sim_lambda, proj_set, model) %>%
  summarize(
    mnDiffN = mean(meanN - abundEst),
    sdDiffN = sd(meanN - abundEst),
    lo_ci = quantile(meanN - abundEst, 0.025),
    hi_ci = quantile(meanN - abundEst, 0.975),
    mnPercentileN = mean(percentile_abundEstN)
  )

ggplot(summ_diff, aes(x = sim_lambda, y = mnDiffN, group = model, color = model)) + 
  geom_errorbar(aes(ymin = lo_ci, ymax = hi_ci), width = 0.001, position = position_dodge(width = 0.001)) +
  geom_line(position = position_dodge(width = 0.001)) +
  geom_point(position = position_dodge(width = 0.001)) +
  scale_color_manual(values = c("red", "black")) +
  labs(y = "Mean difference in Abundance", x = "Simulated Lambda") +
  facet_wrap(~proj_set, ncol = 1) +
  theme_bw()





# Simulation evaluation -----------------------------------------------------
# Stable

# Params
i_lambdas <- seq(1, 1, 0.005)
n_years <- 20
n_start <- 215
se_start <- 20
n_sim <- 10







# Simulation evaluation -----------------------------------------------------
# Constant increase

# Params
i_lambdas <- seq(1.02, 1.08, 0.005)
n_years <- 15
n_start <- 125
se_start <- 10
n_sim <- 10

# Simulate data, stable 10 years, decline 10 years
base <- data.frame(N = rep(n_start, n_years), SE = rep(se_start, n_years))

Ndata_input_list <- purrr::map(i_lambdas, function(y) {
  for (x in 1:n_sim) {
    if (x == 1) {
      N1 = base %>% slice(n_years) %>% pull(N)
      Nsim = N1*y
    } else {
      Nsim[x] = Nsim[x-1]*y
    }
  }
  
  dy <- data.frame(N = Nsim, SE = se_start) %>%
    add_row(base, .before = 1) %>%
    mutate(
      mean_log = calc_mu_log(mu = N, sd = SE),
      sd_log = calc_sd_log(mu = N, sd = SE),
      sd_log = sd_log[1]
    )
  
  do <- list()
  do[[as.character(y)]] <- dy
  return(do)
})

# Run simulations
furrr::future_map(Ndata_input_list, function (i) {
  nam <- names(i)
  Ndata_input <- i[[1]]
  
  Y_retro = c(n_years:(n_years + n_sim - 1))
  sims_y <- purrr::map(Y_retro, function(y) {
    Ndata_input_y <- Ndata_input %>% slice(1:y)
    
    # Input data set
    init_pcfg_data <- list(
      n_dat_yrs = nrow(Ndata_input_y),  
      n_proj_yrs = 2,                        # number of years to project into the future
      mu_logN_hat = Ndata_input_y$mean_log,  # estimated pop abundance, log space
      sigma_logN_hat = Ndata_input_y$sd_log,
      N_harvest = N_harvest[1:2]
    )
    
    init_data <- list(
      init_pcfg_data, init_pcfg_data
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
    
    return(mfit)
  })
  
  save(sims_y, file = here("out", paste0("Harris_2025_simIncline_lambda_", nam, ".dat")))
})


# Evaluating simulations
# Closure thresholds
threshold_N = 192    # Threshold on abundance below which a hunt is closed
threshold_Nmin = 171 # Threshold on minimum abundance below which a hunt is closed


# Difference between actual and projections
summ_proj <- furrr::future_map_dfr(1:length(i_lambdas), function(x) {
  lambda <- i_lambdas[x]
  load(here("out", paste0("Harris_2025_simIncline_lambda_", lambda, ".dat")))
  
  input_data <- Ndata_input_list[[x]][[1]] %>%
    mutate(year = 1:(n_years + n_sim))
  
  out <- purrr::map_dfr(1:length(sims_y), function(y) {
    tfit <- purrr::map(sims_y[[y]], tidy_draws, .progress = T)
    
    
    tbl <- purrr::imap_dfr(tfit, ~ .x %>% 
                             spread_draws(logN_proj[year]) %>% #, logN_pyear3[year]) %>%
                             mutate(
                               proj_set = if_else(year == 1, "1 yr", "2 yr")
                             ) %>%
                             select(proj_set, year, logN = logN_proj) %>% 
                             mutate(
                               year = year + y + (n_years-1),
                               N = exp(logN) 
                             ) %>%
                             filter(year <= max(input_data$year)) %>%
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
                               percentile_abundEstN = ecdf(N)(abundEst)
                             ) %>% 
                             ungroup() %>%
                             mutate(
                               sim_lambda = lambda
                             )
    )
    
    return(tbl)
  })
  
  return(out)
})
summary(lm((meanN - abundEst) ~ model + sim_lambda + year + proj_set, data = summ_proj))
summary(lm((meanN - abundEst) ~ model * sim_lambda + year + proj_set, data = summ_proj))
save(summ_proj, file = here("out", "SimulatedIncline_projections.dat"))

summ_proj <- summ_proj %>%
  filter(sim_lambda %in% c(seq(1.02, 1.08, by = 0.01))) %>%
  mutate(
    lambda = factor(sim_lambda, 
                    levels = unique(sim_lambda), 
                    labels = paste("lambda:", unique(sim_lambda)))
  )
ggplot(summ_proj %>% filter(proj_set == "1 yr"), aes(x = year, group = model)) +
  geom_hline(yintercept = threshold_N) + 
  geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
  geom_ribbon(aes(ymin = percentile_20, ymax = percentile_80, fill = model), alpha = 0.2) +
  geom_line(aes(y = meanN, color = model), alpha = 0.5) +
  geom_point(aes(y = abundEst), color = 'black') +
  scale_fill_manual(values = c("red", "black")) +
  scale_color_manual(values = c("red", "black")) +
  facet_wrap(~ lambda, nrow = 4, labeller = label_parsed) +
  labs(x = "Year", "Abundance (N)") +
  scale_x_continuous(limits = c(16, 25), breaks = 16:25) +
  theme_bw()

ggplot(summ_proj %>% filter(proj_set == "2 yr"), aes(x = year, group = model)) +
  geom_hline(yintercept = threshold_N) + 
  geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
  geom_ribbon(aes(ymin = percentile_20, ymax = percentile_80, fill = model), alpha = 0.2) +
  geom_line(aes(y = meanN, color = model), alpha = 0.5) +
  geom_point(aes(y = abundEst), color = 'black') +
  scale_fill_manual(values = c("red", "black")) +
  scale_color_manual(values = c("red", "black")) +
  facet_wrap(~ lambda, nrow = 4, labeller = label_parsed) +
  labs(x = "Year", "Abundance (N)") +
  scale_x_continuous(limits = c(16, 25), breaks = 16:25) +
  theme_bw()

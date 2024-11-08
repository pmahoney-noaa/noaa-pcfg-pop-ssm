#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Helper functions
# 
# Notes:
# 
# Author: John R. Brandon, PhD [john dot brandon at icf dot com]
# Date: 2024-Aug
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

# Mean of logarithms -----------------------------------------------------------
calc_mu_log = function(mu, sd){
  # Calculate: mean of the distribution on the log scale (aka "meanlog")
  # mu and sd are the expectation and standard deviation in untransformed space
  cv = sd / mu
  log(mu / sqrt(1 + cv * cv))
}
# SD of logarithms -------------------------------------------------------------
calc_sd_log = function(mu, sd){
  # Calculate: standard deviation of the distribution on the log scale (aka "sdlog")
  # mu and sd are the expectation and standard deviation in untransformed space
  cv = sd / mu
  var_log = log(1 + cv * cv)
  sqrt(var_log)  # sd of logarithms
}
# Log-normal RV ----------------------------------------------------------------
gen_lognorm = function(n = 1, mu = 0, sd = 1){
  # Generate a log-normal random deviate given untransformed mean and SD
  # If argument n > 1, return vector of variates
  cv = sd / mu
  sigma_tmp = 1 + cv * cv
  N_rand = log(mu / (sqrt(sigma_tmp))) 
  N_rand = N_rand + rnorm(n = n, mean = 0, sd = 1) * sqrt(log(sigma_tmp))
  exp(N_rand)
}
# tmp = gen_lognorm(n = 1e3, mu = 100, sd = 50)
# mean(tmp); median(tmp); sd(tmp)

# Confidence Intervals Log-Normal ----------------------------------------------
calc_log_ci = function(mu, sd, zz = 1.96){
  # Uses Wade 1998 formulation (Eqn 4)
  cv = sd / mu
  var = sd * sd
  var_log = log(1 + (var / (mu ^ 2)))
  cc = exp(zz * sqrt(var_log))
  data.frame(lo_ci = mu / cc, mu, hi_ci = cc * mu, sd_log = sqrt(var_log))
}
# calc_log_ci(mu = 100, sd = 10)

# Abundance plots --------------------------------------------------------------
tidy_plot_traj = function(input_data, tidy_mcmc, threshold_N, threshold_Nmin){
  N_out = tidy_mcmc %>% 
    spread_draws(logN[year]) %>% 
    select(year, logN) %>% 
    ungroup() %>% 
    mutate(year = year + Ndata_input$year[1] - 1) 
  
  N_proj_out = tidy_mcmc %>% 
    spread_draws(logN_proj[year]) %>% 
    select(year, logN = logN_proj) %>% 
    ungroup() %>% 
    mutate(year = year + max(Ndata_input$year)) 
  
  N_out_table = N_out %>% 
    bind_rows(N_proj_out) %>% 
    group_by(year) %>% 
    summarize(mean = mean(logN),
              median = median(logN),
              lo_ci = quantile(logN, 0.025),
              hi_ci = quantile(logN, 0.975),
              percentile_20 = quantile(logN, 0.2),
              percentile_80 = quantile(logN, 0.8)) %>% 
    ungroup() %>% 
    mutate(across(.cols = -year, .fns = exp)) 
  
  pl <- input_data %>% 
    ggplot(aes(x = year, y = N)) +  
    geom_ribbon(data = N_out_table, aes(y = mean, ymin = percentile_20, ymax = percentile_80), alpha = 0.2) +
    geom_ribbon(data = N_out_table, aes(y = mean, ymin = lo_ci, ymax = hi_ci), alpha = 0.2) +
    geom_line(data = N_out_table, aes(y = mean)) +
    geom_point(size = 3, color = "white", fill = "black", shape = 21) +
    ylim(c(0, 350)) +
    geom_hline(yintercept = threshold_N) + 
    geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
    geom_errorbar(aes(ymin = low_95CI, ymax = high_95CI)) +  
    geom_errorbar(aes(ymin = low_60CI, ymax = N), linewidth = 0) +
    geom_point(aes(y = low_60CI), shape = 23, fill = "red", size = 2) +
    theme_bw(base_size = 16) +
    labs(x = "Year", y = "PCFG")
  
  return(print(pl))
}

# Prediction plots --------------------------------------------------------------
tidy_plot_propBelowThresh = function(N_eval_table){
  N_eval_table <- N_eval_table %>%
    mutate(
      proj_set = factor(proj_set, levels = rep(paste0("pyear", 1:3)), labels = rep(paste("Projected year", 1:3)))
    )
  
  p1 <- N_eval_table %>% 
    ggplot(aes(x = year, y = prop_below_threshold, group = proj_set, color = proj_set)) + 
    geom_line() +
    geom_point() +
    scale_color_brewer(palette = "Reds") +
    labs(x = "Year", y = "Prop. of sims (< N threshold)") +
    guides(color = guide_legend(title = NULL, position = "top", direction = "horizontal")) +
    theme_bw(base_size = 16)
  
  p2 <- N_eval_table %>% 
    ggplot(aes(x = year, y = prop_below_minThreshold, group = proj_set, color = proj_set)) + 
    geom_line() +
    geom_point() +
    scale_color_brewer(palette = "Reds") +
    labs(x = "Year", y = "Prop. of sims (< Nmin threshold)") +
    guides(color = "none") +
    theme_bw(base_size = 16)
  
  return(cowplot::plot_grid(p1, p2, ncol = 1, align = "hv"))
}

tidy_plot_retroPred = function(N_eval_table){
  N_eval_table <- N_eval_table %>%
    mutate(
      proj_set = factor(proj_set, levels = rep(paste0("pyear", 1:3)), labels = rep(paste("Projected year", 1:3)))
    )
  
  width <- 0.5
  pl <- N_eval_table %>% 
    ggplot(aes(x = year, y = abundEst, group = proj_set, color = proj_set)) + 
    geom_errorbar(aes(ymin = loN_ci, ymax = hiN_ci), linetype = 2, width = width, position = position_dodge(width = width)) + 
    geom_errorbar(aes(ymin = percentile_20, ymax = percentile_80), width = width, position = position_dodge(width = width)) +
    geom_point(aes(x = year, y = meanN), position = position_dodge(width = width)) +
    geom_point(size = 2, color = "black", fill = "black", shape = 23) +
    #viridis::scale_color_viridis(discrete = T, option = "D") +
    scale_color_brewer(palette = "Reds") +
    ylim(c(0, 350)) +
    labs(x = "Year", y = "PCFG Abundance") +
    guides(color = guide_legend(title = NULL, position = "top", direction = "horizontal")) +
    theme_bw()
  
  return(print(pl))
}

# Prediction summary table -----------------------------------------------------
pred_summary_tbl = function(Ndata_input, tidy_mcmc, threshold_N, threshold_Nmin){
  tbl <- tidy_mcmc %>% 
    gather_draws(logN_pyear1[year], logN_pyear2[year]) %>% #, logN_pyear3[year]) %>%
    select(year, proj_set = .variable, N = .value) %>% 
    mutate(
      year = year + Ndata_input$year[1] - 1,
      proj_set = gsub("logN_", "", proj_set),
      N = exp(N),
      year = ifelse(proj_set == "pyear1", year + 1, 
                    ifelse(proj_set == "pyear2", year + 2, year + 3)) 
    ) %>%
    filter(year <= max(Ndata_input$year)) %>%
    left_join(Ndata_input %>% dplyr::select(year, abundEstN = N), by = "year") %>%
    #group_by(proj_set, year) %>% 
    summarize(abundEst = mean(abundEstN),
              meanN = mean(N),
              medianN = median(N),
              loN_ci = quantile(N, 0.025),
              hiN_ci = quantile(N, 0.975),
              percentile_20 = quantile(N, 0.2),
              percentile_80 = quantile(N, 0.8),
              rss = sum((N - abundEst)^2),
              closure = meanN < threshold_N | percentile_20 < threshold_Nmin,
              percentile_abundEstN = ecdf(N)(abundEst),
              prop_below_threshold = mean(N < threshold_N),
              prop_below_minThreshold = mean(N < threshold_Nmin)
    ) %>% 
    ungroup()
  
  return(tbl)
}

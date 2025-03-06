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

# Coefficient estimates --------------------------------------------------------
tidy_coefs = function(mfit, model_names){
  
  o <- purrr::imap(mfit, ~ .x$summary() %>%
                     mutate(model = model_names[.y]) %>%
                     filter(grepl("beta", variable)))
  
  do.call('rbind', o) %>%
    dplyr::select(model, variable, mean, median, sd, 
                  lo_ci = q5, hi_ci = q95, rhat, ess_bulk) %>%
    mutate(
      variable = '$\\beta$'
    )
}

# Abundance plots --------------------------------------------------------------
tidy_plot_traj_multimodel = function(input_data, tidy_mcmc, model_names, 
                                     threshold_N, threshold_Nmin, ncols = 1,
                                     ylims = c(0, 350), 
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
  
  pl <- input_data %>% 
    ggplot(aes(x = year, y = N)) +  
    geom_ribbon(data = N_out_table, aes(y = mean, ymin = percentile_20, ymax = percentile_80), alpha = 0.2) +
    geom_ribbon(data = N_out_table, aes(y = mean, ymin = lo_ci, ymax = hi_ci), alpha = 0.2) +
    geom_line(data = N_out_table, aes(y = mean)) +
    geom_point(size = 3, color = "white", fill = "black", shape = 21) +
    scale_x_continuous(limits = c(min(N_out_table$year), max(N_out_table$year)), 
                       breaks = seq(min(N_out_table$year), max(N_out_table$year), by = 1), minor_breaks = NULL) +
    scale_y_continuous(limits = ylims, oob = scales::squish) +
    geom_hline(yintercept = threshold_N, color = "#4d4d4d") + 
    geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "#4d4d4d") +
    geom_errorbar(aes(ymin = low_95CI, ymax = high_95CI)) +  
    geom_errorbar(aes(ymin = low_60CI, ymax = N), linewidth = 0) +
    geom_point(aes(y = low_60CI), shape = 23, fill = "red", size = 2) +
    theme_bw(base_size = 16) +
    facet_wrap(~ model, ncol = ncols) +
    labs(x = "Year", y = "PCFG Abundance") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    NULL
  
  return(print(pl))
}

tidy_plot_traj = function(input_data, tidy_mcmc, threshold_N, threshold_Nmin){
  N_out = tidy_mcmc %>% 
    spread_draws(logN[year]) %>% 
    select(year, logN) %>% 
    ungroup() %>% 
    mutate(year = year + input_data$year[1] - 1) 
  
  N_proj_out = tidy_mcmc %>% 
    spread_draws(logN_proj[year]) %>% 
    select(year, logN = logN_proj) %>% 
    ungroup() %>% 
    mutate(year = year + max(input_data$year)) 
  
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
    scale_x_continuous(limits = c(min(N_out_table$year), max(N_out_table$year)), 
                       breaks = seq(min(N_out_table$year), max(N_out_table$year), by = 1), minor_breaks = NULL) +
    scale_y_continuous(limits = c(0, 350), breaks = seq(0, 350, by = 100)) +
    geom_hline(yintercept = threshold_N) + 
    geom_hline(yintercept = threshold_Nmin, linetype = 2, color = "red") +
    geom_errorbar(aes(ymin = low_95CI, ymax = high_95CI)) +  
    geom_errorbar(aes(ymin = low_60CI, ymax = N), linewidth = 0) +
    geom_point(aes(y = low_60CI), shape = 23, fill = "red", size = 2) +
    theme_bw(base_size = 16) +
    labs(x = "Year", y = "PCFG Abundance") +
    scale_x_continuous(breaks = 1990:2100, minor_breaks = NULL) +
    # theme(
    #   axis.text.x = element_text(angle = 45, hjust = 1)
    # ) +
    NULL
  
  return(print(pl))
}

# Prediction plots -------------------------------------------------------------
tidy_plot_propBelowThresh = function(N_eval_table){
  N_eval_table <- N_eval_table %>%
    mutate(
      proj_set = factor(proj_set, levels = rep(paste0("pyear", 1:3)), labels = rep(paste("Projected year", 1:3)))
    )
  
  p1 <- N_eval_table %>% 
    ggplot(aes(x = year, y = prop_below_threshold, group = proj_set, color = proj_set)) + 
    geom_line() +
    geom_point() +
    scale_color_paletteer_d("Manu::Kotare") +
    labs(x = "Year", y = "Prop. of sims (< N threshold)") +
    guides(color = guide_legend(title = NULL, position = "top", direction = "horizontal")) +
    theme_bw(base_size = 16) +
    scale_x_continuous(breaks = 1990:2100, minor_breaks = NULL) +
    NULL
  
  p2 <- N_eval_table %>% 
    ggplot(aes(x = year, y = prop_below_minThreshold, group = proj_set, color = proj_set)) + 
    geom_line() +
    geom_point() +
    scale_color_paletteer_d("Manu::Kotare") +
    labs(x = "Year", y = "Prop. of sims (< Nmin threshold)") +
    guides(color = "none") +
    theme_bw(base_size = 16) +
    scale_x_continuous(breaks = 1990:2100, minor_breaks = NULL) +
    NULL
  
  return(cowplot::plot_grid(p1, p2, ncol = 1, align = "hv"))
}

tidy_plot_retroPred = function(N_eval_table, ylims = c(0, 350), truncated_retro = F){
  
  if(truncated_retro) {
    N_eval_table <- N_eval_table %>%
      mutate(
        proj_set = factor(proj_set, levels = c("1 yr", "2 yr"), 
                          labels = rep(paste("Projected year", 1:2)))
      )
  } else {
    N_eval_table <- N_eval_table %>%
      mutate(
        proj_set = factor(proj_set, levels = rep(paste0("pyear", 1:3)), 
                          labels = rep(paste("Projected year", 1:3)))
      )
  }
  
  width <- 0.5
  
  if("model" %in% names(N_eval_table)) {
    pl <- N_eval_table %>% 
      ggplot(aes(x = year, y = abundEst, group = model, color = model)) +
      facet_wrap(~ proj_set, ncol = 1)
  } else {
    pl <- N_eval_table %>% 
      ggplot(aes(x = year, y = abundEst, group = proj_set, color = proj_set))
  }
  
  pl <- pl + 
    geom_errorbar(aes(ymin = loN_ci, ymax = hiN_ci), linetype = 2, width = width, 
                  position = position_dodge(width = width)) + 
    geom_errorbar(aes(ymin = percentile_20, ymax = percentile_80), width = width, 
                  position = position_dodge(width = width)) +
    #geom_point(aes(x = year, y = meanN), position = position_dodge(width = width)) +
    geom_point(aes(x = year, y = medianN), position = position_dodge(width = width)) +
    geom_point(size = 2, color = "black", fill = "black", shape = 23) +
    scale_color_paletteer_d("Manu::Kotare") +
    coord_cartesian(ylim = ylims) +
    #scale_y_continuous(limits = ylims, oob = scales::squish) +
    labs(x = "Year", y = "PCFG Abundance") +
    scale_x_continuous(limits = c(min(N_eval_table$year) - 1.5, max(N_eval_table$year)) + 1, 
                       breaks = seq(min(N_eval_table$year), max(N_eval_table$year), 1),
                       expand = c(0,0)) +
    guides(color = guide_legend(title = NULL, position = "top", direction = "horizontal", nrow = 1)) +
    theme_bw()
  
  return(print(pl))
}

# Prediction summary table -----------------------------------------------------
pred_summary_tbl_multimodel = function(input_data, tidy_mcmc, model_names, 
                                       threshold_N, threshold_Nmin){
  tbl = purrr::imap_dfr(tidy_mcmc, ~ .x %>% 
                          gather_draws(logN_pyear1[year], logN_pyear2[year]) %>% #, logN_pyear3[year]) %>%
                          select(year, proj_set = .variable, N = .value) %>% 
                          mutate(
                            year = year + input_data$year[1] - 1,
                            proj_set = gsub("logN_", "", proj_set),
                            N = exp(N),
                            year = ifelse(proj_set == "pyear1", year + 1, 
                                          ifelse(proj_set == "pyear2", year + 2, year + 3)) 
                          ) %>%
                          filter(year <= max(input_data$year)) %>%
                          left_join(input_data %>% dplyr::select(year, abundEstN = N), by = "year") %>%
                          #group_by(proj_set, year) %>% 
                          summarize(
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
                            prop_below_threshold = mean(N < threshold_N),
                            prop_below_minThreshold = mean(N < threshold_Nmin)
                          ) %>% 
                          ungroup()
                        )
  return(tbl)
}

pred_summary_tbl = function(input_data, tidy_mcmc, threshold_N, threshold_Nmin){
  tbl <- tidy_mcmc %>% 
    gather_draws(logN_pyear1[year], logN_pyear2[year]) %>% #, logN_pyear3[year]) %>%
    select(year, proj_set = .variable, N = .value) %>% 
    mutate(
      year = year + input_data$year[1] - 1,
      proj_set = gsub("logN_", "", proj_set),
      N = exp(N),
      year = ifelse(proj_set == "pyear1", year + 1, 
                    ifelse(proj_set == "pyear2", year + 2, year + 3)) 
    ) %>%
    filter(year <= max(input_data$year)) %>%
    left_join(input_data %>% dplyr::select(year, abundEstN = N), by = "year") %>%
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

# Model averaged abundance -----------------------------------------------------
tidy_model_avg <- function (tfit, wgts, iters, input_data, seed = 1001) {
  
  set.seed(seed)
  
  N_out = map2(tfit, wgts, ~ .x %>%
                     slice(sample(1:n(), size = round(.y*iters), replace = F)) %>% 
                     spread_draws(logN[year]) %>% 
                     select(year, logN, .draw) %>% 
                     ungroup() %>% 
                     mutate(year = year + input_data$year[1] - 1))
  
  N_proj_out = map2_dfr(tfit, N_out, ~ .x %>%
                          spread_draws(logN_proj[year]) %>% 
                          filter(.draw %in% unique(.y$.draw)) %>%
                          select(year, logN = logN_proj, .draw) %>% 
                          ungroup() %>% 
                          mutate(year = year + max(input_data$year)))
  
  N_out = do.call("rbind", N_out) %>%
    add_row(N_proj_out)
  
  # N_out = map2_dfr(tfit, wgts, ~ .x %>%
  #                    slice(sample(1:n(), size = round(.y*iters), replace = F)) %>% 
  #                    spread_draws(logN[year]) %>% 
  #                    select(year, logN, .draw) %>% 
  #                    ungroup() %>% 
  #                    mutate(year = year + input_data$year[1] - 1))
  
  # N_proj_out = map2_dfr(tfit, wgts, ~ .x %>%
  #                         slice(sample(1:n(), size = round(.y*iters), replace = F)) %>% 
  #                         spread_draws(logN_proj[year]) %>% 
  #                         select(year, logN = logN_proj) %>% 
  #                         ungroup() %>% 
  #                         mutate(year = year + max(input_data$year)))
  
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
  
  return(list(draws = N_out, summary = N_out_table))
}

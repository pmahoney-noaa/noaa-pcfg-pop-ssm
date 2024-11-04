#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Read in retrospective time series of ENP calf estimates
# 
# Notes: Contains estimates of calf abundance from 1994 - 2024, updates are 
#        provided in the Spring of the current year.
# 
# Author: Peter J. Mahoney, PhD [peter dot mahoney at noaa dot gov]
#  Modified from read_abun_retro for consistency
# Date: Oct '24
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

library(pacman)
p_load(tidyverse, here)
source(here("R", "helper_functions.R"))

# Read Data --------------------------------------------------------------------
retro_dat_raw <- read_csv(here("data", "Calf_Estimates_v3_Mv1_2024-09-16.csv"))

# Wrangle Data -----------------------------------------------------------------
retro_dat <- retro_dat_raw %>%
  rename_with(tolower) %>%
  rename(N = mean) %>% select(-median) %>%
  mutate(lambda = N / lag(N)) %>% 
  #filter(year >= 2002) %>%   # Drop 1998-2001 as some early estimates not considered reliable (discovery curve etc.)
  mutate(year_id = row_number() - 1) %>% 
  mutate(init_N = first(na.omit(N))) %>% 
  mutate(N_delta_perc = ifelse(!is.na(N), (N - init_N) / init_N, NA),
         estimator = method[1],
         low_60CI = calc_log_ci(mu = N, sd = se, zz = qnorm(p = 0.8))$lo_ci,
         high_60CI = calc_log_ci(mu = N, sd = se, zz = qnorm(p = 0.8))$hi_ci,
         low_95CI = calc_log_ci(mu = N, sd = se, zz = 1.96)$lo_ci,
         high_95CI = calc_log_ci(mu = N, sd = se, zz = 1.96)$hi_ci,
         mean_log = calc_mu_log(mu = N, sd = se),
         sd_log = calc_sd_log(mu = N, sd = se),
         cv = se / N) 

retro_dat
# View(retro_dat)

# Simplifying
Cdata = retro_dat %>% 
  #filter(estimator == "Stewart&Weller") %>% 
  select(-N_delta_perc, -estimator, -init_N) %>% 
  select(year:se, lambda:sd_log, cv, method, year_id)

# Garbage collection
rm(list = c("retro_dat", "retro_dat_raw"))

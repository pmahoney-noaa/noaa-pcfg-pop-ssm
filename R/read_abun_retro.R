#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Read in retrospective time series of PCFG abundance estimates
# 
# Notes:
# 
# Author: John R. Brandon, PhD [john dot brandon at icf dot com]
# Date: Oct '24
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

library(pacman)
p_load(tidyverse, readxl, here)
source(here("R", "helper_functions.R"))

# Read Data --------------------------------------------------------------------
retro_dat_raw = read_excel(path = here("data", "PCFG_retrospective_abundance_time_series_JB.xlsx"),
                           sheet = "N_retro_long")  # sheet = "All time series together"

# Key-Val table for refs -------------------------------------------------------
ref_table = retro_dat_raw %>% 
  distinct(ref) %>% 
  mutate(ref_id = row_number())

ref_table

# Wrangle Data -----------------------------------------------------------------
retro_dat = retro_dat_raw %>% 
  filter(ref != "Calambokidis et al. 2010") %>%
  left_join(ref_table, by = "ref") %>% 
  group_by(ref) %>% 
  mutate(lambda = N / lag(N)) %>% 
  filter(year >= 2002) %>%   # Drop 1998-2001 as some early estimates not considered reliable (discovery curve etc.)
  mutate(year_id = row_number() - 1) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mutate(init_N = first(na.omit(N))) %>% 
  ungroup() %>% 
  mutate(N_delta_perc = ifelse(!is.na(N), (N - init_N) / init_N, NA),
         estimator = ifelse(ref == "Calambokidis et al. 2010", "POPAN", "JS1"),
         low_60CI = calc_log_ci(mu = N, sd = SE, zz = qnorm(p = 0.8))$lo_ci,
         high_60CI = calc_log_ci(mu = N, sd = SE, zz = qnorm(p = 0.8))$hi_ci,
         low_95CI = calc_log_ci(mu = N, sd = SE, zz = 1.96)$lo_ci,
         high_95CI = calc_log_ci(mu = N, sd = SE, zz = 1.96)$hi_ci,
         ref_yr = substr(ref, start = nchar(ref) - 3, nchar(ref)),
         mean_log = calc_mu_log(mu = N, sd = SE),
         sd_log = calc_sd_log(mu = N, sd = SE),
         cv = SE / N) 

retro_dat
# View(retro_dat)

# Most recent time series of abundance estimates
Ndata = retro_dat %>% 
  filter(ref == "Harris et al. 2024") %>% 
  select(-N_delta_perc, -estimator, -ref_yr, -init_N) %>% 
  select(year:SE, lambda:sd_log, cv, ref_id, ref, year_id)

# Garbage collection
rm(list = c("ref_table", "retro_dat", "retro_dat_raw"))

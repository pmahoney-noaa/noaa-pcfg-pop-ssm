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
retro_dat_raw <- read_csv(here(
  "data",
  "PCFG_calf_counts.csv"
))

# Wrangle Data -----------------------------------------------------------------
N_pcfg_calves <- retro_dat_raw %>%
  rename(
    year = Year
  )

# Garbage collection
rm(list = c("retro_dat_raw"))

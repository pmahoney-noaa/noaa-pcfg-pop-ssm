#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Read in retrospective time series of PCFG abundance estimates
# 
# Notes: Contains stranding data from 1998 through 2023, updates are provided 
#        in the Spring of the following year.
# 
# Author: Peter J. Mahoney, PhD [peter dot mahoney at noaa dot gov]
#  Modified from read_abun_retro for consistency
# Date: Oct '24
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

library(pacman)
p_load(tidyverse, readxl, here)

# Read Data --------------------------------------------------------------------
retro_orwa_strand = read_excel(path = here("data", "WCR_Er_Strandings_1998-2023.xlsx"),
                             sheet = "OR-WA Er 1998-2023") %>%
  select(
    national_database_number = `National Database Number`,
    species = `Common Name`,
    examiner = `Examiner Name`,
    county = County, state = State,
    longitude = Longitude, latitude = Latitude,
    obs_date = `Observation Date`,
    exam_date = `Date of Examination`,
    obs_status = `Observation Status`,
    type = `Report Type`,
    anthro = `human interaction contributed`
  )

retro_ca_strand_late = read_excel(path = here("data", "WCR_Er_Strandings_1998-2023.xlsx"),
                                  sheet = "CA Er 2006-2023") %>%
  select(
    national_database_number = `National Database Number`,
    species = `Common Name`,
    examiner = `Examiner Name`,
    county = County, state = State,
    longitude = Longitude, latitude = Latitude,
    obs_date = `Observation Date`,
    exam_date = `Date of Examination`,
    obs_status = `Observation Status`,
    type = `Report Type`,
    anthro = `human interaction contributed`
  )

retro_ca_strand_early = read_excel(path = here("data", "WCR_Er_Strandings_1998-2023.xlsx"),
                                  sheet = "CA Er 1998-2005") %>%
  mutate(
    national_database_number = NA, state = "CA", anthro = NA,
    species = "Whale, gray"
  ) %>%
  select(
    national_database_number,
    species,
    examiner = `Examiner`,
    county = County, state,
    longitude = Longitude, latitude = Latitude,
    obs_date = `Date of Occurance`,
    exam_date = `Date of Examination`,
    obs_status = `Initial Condition`,
    type = `Type of Occurance`,
    anthro
  )

retro_strand <- rbind(retro_orwa_strand, retro_ca_strand_early, retro_ca_strand_late)

# Wrangle Data -----------------------------------------------------------------
strand_dat <- retro_strand %>%
  mutate(year = year(obs_date)) %>%
  group_by(species, state, year) %>%
  summarize(
    N_strand = n()
  ) %>%
  ungroup()

strand_dat_all <- retro_strand %>%
  mutate(year = year(obs_date)) %>%
  group_by(species, year) %>%
  summarize(
    N_strand = n()
  ) %>%
  ungroup() %>%
  add_column(state = "All", .before = "year")

Sdata <- rbind(strand_dat, strand_dat_all)
# View(Sdata)


# Garbage collection
rm(list = c("retro_orwa_strand", "retro_ca_strand_late", "retro_ca_strand_early", "strand_dat_all", "retro_strand", "strand_dat"))

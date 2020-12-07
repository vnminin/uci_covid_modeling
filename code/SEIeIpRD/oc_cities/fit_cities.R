library(tidyverse)
library(stemr)
library(extraDistr)
library(fs)
library(here)
library(foreach)
library(doRNG)
library(doParallel)
registerDoParallel(cores = future::availableCores())

source('code/SEIeIpRD/oc_cities/oc_cities_functions.R')
model_name <- "SEIeIpRD"
loc_name <- "oc"

last_2_folder_names <- tibble(dirs = list.dirs(here("code", model_name, loc_name), recursive = F, full.names = F)) %>%
  filter(str_detect(dirs, "^\\d{4}-\\d{2}-\\d{2}_\\d{4}-\\d{2}-\\d{2}$")) %>%
  separate(dirs, c("start", "end"), "_") %>%
  arrange(end) %>%
  tail(2) %>%
  unite(col = dir, sep = "_") %>%
  pull(dir)

folder_name <- last_2_folder_names[2]
prev_folder_name <- last_2_folder_names[1]

fit_and_save_city_model(city_name = "Santa Ana", folder_name = folder_name, prev_folder_name = prev_folder_name, time_interval_in_days = 3)
fit_and_save_city_model(city_name = "Irvine", folder_name = folder_name, prev_folder_name = prev_folder_name, time_interval_in_days = 3)

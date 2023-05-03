# General workflow for testing with data

# 1) Wrangle the survey data for trawl to start
# 2a) Fit test models to single species
# 2b) Convert to function using try() to fit all models
# 3a) Get index from each fitted model;
# 3b) Convert to function using try() to apply to multiple species
# 4) 

# 00) Describe the differences in the species, particularly with respect to 
#     how characteristics of the different species compare or relate to the 
#     simulated data examples. 
# E.g., + range of autocorrelation
#       + gap between surveys
#       + expectations in variability from year to year

# Questions
# - "Survey domain year" --> this has to do with when the grid was adjusted? 
#   I.e., the last time it was adjusted to account for dropped cells/areas?
# - Can we walk through the different model scenarios to make sure I understand? 
#   it might help me if I could make a table describing what each is doing
# - Why is offset area divided by 100 000 instead of 1 000 000? m2 to km2 should 
#   be 1e-6?
#   A: Use 1e5 because this gets the log(offset) closer to zero
# - Should priors be used? Should they be standardised like they were in the 
# - sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L); default? nlminb_loops?
# - index_args; an argument in sdmTMB; these should be a default in fitting the 
#   model fitting? 
# - How to decide what k should be when using s(year)
# - I need to remind myself what this bias_correct is about for the index fitting
# - If an offset is included in the original model; does this follow through in
#   when we use predict.sdmTMB

# Should 2021 WCVI be left out? Should the analysis be compared with and without 
# this year? To see if it has any affect? 
# distinct(dat, year, survey_abbrev) %>% arrange(survey_abbrev, year)
# dat %>% filter(year == 2021, survey_abbrev == "SYN WCVI")


library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)  #install.packages('beepr'); for the cat gifs: https://www.sumsar.net/blog/2014/01/announcing-pingr/

source(here::here('stitch', 'utils.R'))

mytheme <- function() ggsidekick::theme_sleek()  # sometimes I add more layers to themes
theme_set(mytheme())

# Data preparations
# ------------------------------------------------------------------------------
# Specify model order and labels for plotting ----------------------------------
model_lookup <- 
  tibble(id = 1:8, 
         desc = c("st = 'rw'", # 6
                  "st IID covariate", # 3
                  "st IID s(year)", # 5
                  "st IID no covariate as.factor year", # 2????
                  "st time_varying RW", # 4
                  "st (1|year)", # 1
                  "spatial only",  # 7
                  "st (1|region)"),  # 8
         order = c(6, 3, 5, 2, 4, 1, 7, 8))  # try matching order of the simulated plots

# Prepare grids ----------------------------------------------------------------
syn_grid <- 
  gfplot::synoptic_grid %>%
  tibble() %>%  # because I accidentally print the full df too often
  dplyr::select(-survey_series_name, -utm_zone, -survey_domain_year) %>% 
  mutate(log_depth = log(depth), region = as.factor(survey))

# Clean survey data ------------------------------------------------------------
dat <- 
  readRDS(here::here('data/all_surv_catch.rds')) %>%  # What was the code that actually made this (from SOPO)
  mutate(density_kgkm2 = density_kgpm2 * 1000000, 
         log_depth = log(depth_m), 
         area_swept1 = doorspread_m * (speed_mpm * duration_min), 
         area_swept2 = tow_length_m * doorspread_m, 
         area_swept = ifelse(!is.na(area_swept2), area_swept2, area_swept1)) %>% 
  mutate(offset = log(area_swept / 1e5)) %>%  # Value used for offset
  #filter(!(year == 2021 & survey_abbrev == "SYN WCVI")) %>%  # this region not usually surveyed in odd years
  sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 32609) %>% 
  # simplify df columns
  select(survey_id, trip_id, fishing_event_id, survey_abbrev,
         year, month, day, latitude, longitude, X, Y,
         depth_m, log_depth,
         species_code, species_common_name, 
         catch_weight, catch_count,
         density_kgkm2, density_pcpm2, density_ppkm2, 
         area_swept, offset, hook_count, time_deployed) %>%
  # specify factor variables that will be used in models
  mutate(fyear = as.factor(year), 
         region = as.factor(survey_abbrev)) %>%
  # use complete datasets
  drop_na(area_swept) %>%   # drop empty area swept (no doorspread given)
  drop_na(depth_m)          # drop rows without depths

#hist(log(dat$area_swept / 1e5))  # check that log offset is close to 0
#source(here::here('stitch', 'explore-grids.R'))
survey_region_baseplot + 
  geom_point(data = dat %>% filter(str_detect(survey_abbrev, "SYN")), aes(colour = survey_abbrev), shape = 21, alpha = 0.5) + 
  facet_wrap(~ year) + 
  guides(colour = guide_legend(override.aes = list(shape = 1, alpha = 1)))

# Start with arrowtooth from synoptic trawl for now
arrow <- 
  filter(dat, str_detect(survey_abbrev, "SYN")) %>% 
  filter(., species_common_name == 'arrowtooth flounder')
# What data are missing and are there patterns in the missing data?
# - For arrowtooth there are no patterns
# - Area swept: One trip in one year had missing data (15 rows)
# - Depth: 17 rows are missing; 12 of these are from 2010, the rest are scattered throughout

arrow_no2014 <- 
  filter(dat, str_detect(survey_abbrev, "SYN")) %>% 
  filter(., species_common_name == 'arrowtooth flounder') %>% 
  filter(year != 2014)

bocaccio <- 
  filter(dat, str_detect(survey_abbrev, "SYN")) %>% 
  filter(., species_common_name == 'bocaccio')

# Inputs for predictions and index calcluations --------------------------------
fitted_yrs <- sort(unique(arrow$year))

fitted_yrs <- sort(unique(bocaccio$year))
fitted_yrs <- sort(unique(arrow_no2014$year))
nd <- 
  make_grid(syn_grid, years = fitted_yrs) %>% 
  mutate(fyear = as.factor(year))

# Fit models -------------------------------------------------------------------
ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)

bocaccio_mesh <- make_mesh(bocaccio, c("X", "Y"), cutoff = 20)
plot(bocaccio_mesh)
bocaccio_fits <- fit_models(bocaccio)
beep()

# Is it of interest to consider what ranges are estimated for the different models
# given that the gap size is 
sdmTMB:::print_range(fit1)
sdmTMB:::print_range(fit2)
sdmTMB:::print_range(fit3)
arrow_no2014_fits <- fit_models(dat = arrow_no2014, data_subset = "arrowtooth - 2014 removed")
beep()

# Predict on grid and calculate indices ----------------------------------------
preds <- get_pred_list(fit_list = arrow_no2014_fits, newdata = nd)
indices <- get_index_list(pred_list = preds)
beep()

preds <-  
  purrr::map(fits, function(.x) {
    if (inherits(.x, "sdmTMB")) {
      out <- predict(.x, newdata = nd, return_tmb_object = TRUE)
    } else {
      out <- NA
    }
    out
  })
beep()

indices <- 
  purrr::map(preds, function(.x) {
    if (length(.x) > 1) {
      get_index(.x, bias_correct = TRUE)
    }
  })
beep()

index_df <- 
  enframe(indices) %>% 
  unnest(col = "value") %>% 
  separate(col = 'name', into = c('id', 'species')) %>% 
  mutate(id = as.numeric(id)) %>% 
  as_tibble() %>% 
  right_join(., model_lookup) %>% 
  # Add this for quick and dirty north/south, but note that WCVI data from 2021
  # was in this arrow dataset used. 
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd"))

# Plot index over time ---------------------------------------------------------
ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_wrap(~ fct_reorder(desc, order), nrow = 1L, scales = "free_y")

# What about MRE, RMSE, see-saw, coverage etc. ? --------------------------
seesaw_metrics <- 
  index_df %>% 
    group_by(id, species, desc, order, odd_even) %>% 
    summarise(mean_est = mean(est)) %>%  # should this be log estimate? but then log ratio of log???
    pivot_wider(names_from = odd_even, values_from = mean_est) %>% 
    select(-`NA`) %>% 
    mutate(log_ratio = log(even / odd)) %>% 
    mutate(seesaw_index = abs(log_ratio)) %>% 
    left_join(., 
      index_df %>%   # add mean_se per time-series
      group_by(id, species, desc, order) %>% 
      summarise(mean_se = mean(se))) %>% 
    ungroup()

seesaw_metrics
       # mre = mean(log_ratio),  # I don't know if this metric makes sense to keep/use?
       # rmse = sqrt(mean(log_ratio^2)),  # Nor this one?
       #coverage = mean(total < upr & total > lwr))  # or this one?

seesaw_metrics %>% 
  arrange(seesaw_index) |>
    mutate(desc = forcats::fct_reorder(as.factor(desc), rev(seesaw_index), .na_rm = FALSE)) |>
    tidyr::pivot_longer(cols = c(log_ratio, seesaw_index, mean_se), names_to = "metric") |>
    mutate(metric = factor(metric,
      levels = c('seesaw_index', 'mean_se', 'log_ratio'))) |>
      #levels = c("seesaw_index", "rmse", "mre", "mean_se", "coverage"))) |>
    ggplot(aes(value, desc)) +
    geom_point(pch = 21, size = 1.6) +
    facet_wrap(~metric, scales = "free_x", nrow = 1L) +
    ggsidekick::theme_sleek() +
    theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y.left = element_blank()) +
    xlab("Metric value")


# TODO: RERUN arrow and bocaccio with above and check with and without 2014
# TODO: setup functions to produce plots
# TODO: get best HBLL grids; start with inside HBLL and use quillback and see 
#       what output looks likeindexes_df |>
  left_join(actual) |>
  left_join(select(d, year, sampled_region) %>% distinct()) |>
  group_by(model) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[sampled_region == "north"]) - mean(log_residual[sampled_region == "south"])),
    mre = mean(log_residual),
    rmse = sqrt(mean(log_residual^2)),
    mean_se = mean(se),
    coverage = mean(total < upr & total > lwr)
  )
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

# Species-specific mesh --------------------------------------------------------
mesh <- make_mesh(arrow, c("X", "Y"), cutoff = 20)

# Inputs for predictions and index calcluations --------------------------------
fitted_yrs <- sort(unique(arrow$year))


# Fit models -------------------------------------------------------------------
ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)

cli::cli_inform("Fitting st = 'rw'")
fit1 <- 
  sdmTMB(
    catch_weight ~ 1,
    family = tweedie(),
    data = arrow, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = mesh,
    offset = arrow$offset,
    control = ctrl
  )
  beep()

cli::cli_inform("Fitting st IID covariate")
fit2 <- try(sdmTMB(
    catch_weight ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
    family = tweedie(),
    data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
    silent = TRUE, mesh = mesh,
    offset = arrow$offset,
    control = ctrl
  ))

cli::cli_inform("Fitting st IID s(year)")
fit3 <- sdmTMB(
  catch_weight ~ s(year),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  offset = arrow$offset,
  control = ctrl
)
beepr::beep()

cli::cli_inform("Fitting st IID no covariate as.factor year")
fit4 <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = arrow$offset,
  control = ctrl
)
beepr::beep()

cli::cli_inform("Fitting st time_varying RW")
fit5 <- sdmTMB(
  catch_weight ~ 0,
  family = tweedie(),
  time_varying = ~1, time_varying_type = "rw",
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = arrow$offset,
  control = ctrl
)

cli::cli_inform("Fitting st (1|year)")
fit6 <- sdmTMB(
  catch_weight ~ 1 + (1 | fyear),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = arrow$offset,
  control = ctrl
)

cli::cli_inform("Fitting spatial only")
fit7 <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "off", spatial = "on",
  mesh = mesh,
  offset = arrow$offset,
  control = ctrl
)
beepr::beep()

cli::cli_inform("Fitting st (1 | region")
fit8 <- sdmTMB(
  catch_weight ~ 0 + fyear + region,
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = arrow$offset,
  control = ctrl
)
beep()

# Is it of interest to consider what ranges are estimated for the different models
# given that the gap size is 
sdmTMB:::print_range(fit1)
sdmTMB:::print_range(fit2)
sdmTMB:::print_range(fit3)

# Predict on grid and calculate indices ----------------------------------------
fits <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8)


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
  bind_rows(indices, .id = "id") %>% 
  mutate(id = as.numeric(id)) %>% 
  as_tibble() %>% 
  left_join(., model_lookup) %>% 
  # Add this for quick and dirty north/south, but note that WCVI data from 2021
  # was in this arrow dataset used. 
  mutate(sampled_region = ifelse(year %% 2 == 1, "north", "south"))

# Plot index over time ---------------------------------------------------------
ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = sampled_region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62")) + 
  labs(colour = "Sampled region") + 
  facet_wrap(~ fct_reorder(desc, order), nrow = 1L)


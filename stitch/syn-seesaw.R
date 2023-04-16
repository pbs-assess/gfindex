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
# - Should priors be used? Should they be standardised like they were in the 
# - sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L); default? nlminb_loops?

# Interesting that including random effect of survey block did not fix 
# see-saw like Patrick says it did in his HMSC paper. I wonder why not?


library(tidyverse)
library(gfplot)
library(sdmTMB)

source(here::here('stitch', 'utils.R'))

dat <- 
# What was the code that actually made this 
  readRDS(here::here('data/all_surv_catch.rds')) %>% 
  mutate(density_kgkm2 = density_kgpm2 * 1000000, 
         log_depth = log(depth_m), 
         area_swept1 = doorspread_m * (speed_mpm * duration_min), 
         area_swept2 = tow_length_m * doorspread_m, 
         area_swept = ifelse(!is.na(area_swept2), area_swept2, area_swept1)) %>% 
  mutate(area_swept_km2 = area_swept2 / 1000000) %>%  # This should be used for offset???
  mutate(log_area = log(area_swept_km2)) %>%  # Value used for offset
  sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 32609)

# Focus on synoptic trawl data to start

# Start with arrowtooth from synoptic trawl to start
arrow <- 
  filter(dat, str_detect(survey_abbrev, "SYN")) %>% 
  filter(., species_common_name == 'arrowtooth flounder') %>%
  select(survey_id, trip_id, fishing_event_id, 
         year, month, day, latitude, longitude, X, Y,
         depth_m, log_depth,
         species_code, species_common_name, 
         catch_weight, catch_count,
         density_kgkm2, density_pcpm2, density_ppkm2, 
         # area_swept1, area_swept2, doorspread_m, speed_mpm, duration_min, 
         # tow_length_m,
         area_swept, hook_count, time_deployed) %>%
  mutate(fyear = as.factor(year)) %>% 
  drop_na(area_swept) %>%   # drop empty area swept (no doorspread given)
  drop_na(depth_m)   # drop rows without depths

# Only one year (one trip specifically) had missing data on area swept
  # For arrowtooth this was only 15 rows of data
# How many rows do not have depth, and is there a pattern in which ones are 
# missing depth? 
# For arrowtooth only 17 rows are missing depth data, 12 of these are from 2010; 
  # the rest are scattered throughout. 

mesh <- make_mesh(arrow, c("X", "Y"), cutoff = 20)

model_lookup <- 
  tibble(id = 1:7, 
         desc = c("st = 'rw'", 
                  "st IID covariate", 
                  "st IID s(year)", 
                  "st IID no covariate as.factor year", 
                  "st time_varying RW",
                  "st (1|year)",   # I don't really understand this model
                  "spatial only"))

# grid for whole coast synoptic indices

syn_grid <- gfplot::synoptic_grid %>% tibble()

grid_locs <- gfplot::synoptic_grid %>%
      dplyr::select(X, Y, depth, survey)

fit_models <- function(dat, mesh,ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L))
  fits <- list()
  nms <- c()
  i <- 1

cli::cli_inform("Fitting st = 'rw'")
fit1 <- try(
  sdmTMB(
    catch_weight ~ 1,
    family = tweedie(),
    data = arrow, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = mesh,
    offset = log(arrow$area_swept / 100000), # why divided by this?
    control = sdmTMBcontrol(newton_loops = 1)
  )
)


cli::cli_inform("Fitting st IID covariate")
fit2 <- sdmTMB(
    catch_weight ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
    family = tweedie(),
    data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
    silent = TRUE, mesh = mesh,
    #priors = priors,
    offset = log(arrow$area_swept / 100000),
    control = sdmTMBcontrol(newton_loops = 1)
  )
sanity(fit2)
beepr::beep()

cli::cli_inform("Fitting st IID s(year)")
fit3 <- sdmTMB(
  catch_weight ~ s(year),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  #priors = priors,
  offset = log(arrow$area_swept / 100000),
  control = sdmTMBcontrol(newton_loops = 1)
)
beepr::beep()

cli::cli_inform("Fitting st IID no covariate as.factor year")
fit4 <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = log(arrow$area_swept / 100000),
  control = sdmTMBcontrol(newton_loops = 1)
)
beepr::beep()

cli::cli_inform("Fitting st time_varying RW")
fit5 <- sdmTMB(
  catch_weight ~ 0,
  family = tweedie(),
  time_varying = ~1,
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = log(arrow$area_swept / 100000),
  control = sdmTMBcontrol(newton_loops = 1)
)

cli::cli_inform("Fitting st (1|year)")
fit6 <- sdmTMB(
  catch_weight ~ 1 + (1 | fyear),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = log(arrow$area_swept / 100000),
  control = sdmTMBcontrol(newton_loops = 1)
)

cli::cli_inform("Fitting spatial only")
fit7 <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  family = tweedie(),
  data = arrow, time = "year", spatiotemporal = "off", spatial = "on",
  mesh = mesh,
  offset = log(arrow$area_swept / 100000),
  control = sdmTMBcontrol(newton_loops = 1)
)

sdmTMB:::print_range(fit1)
sdmTMB:::print_range(fit2)
sdmTMB:::print_range(fit3)



# model_name <- "Spatial only, no covariate"
# fit <- sdmTMB(
#   catch_weight ~ 0 + fyear, 
#   family = tweedie(), 
#   data = arrow, 
#   time = "year", spatiotemporal = "off", spatial = "on",
#   mesh = mesh, 
#   offset = log(arrow$area_swept / 100000), # why divided by this?
#   control = sdmTMBcontrol(newton_loops = 1)
# )

# Predict on grid and calculate indexes -----------------------------------
nd <- select(predictor_dat, X, Y, year, fyear, region)



out <- predict(.x, newdata = nd, return_tmb_object = TRUE)

get_index(.x, bias_correct = TRUE)



select(arrow, area_swept1, area_swept2, area_swept) %>% 
filter(is.na(area_swept))

    observed ~ 1,
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors,
    control = ctl

# FIX ME: this model won't run
# model_type <- "IID st (1|fyear) cov = FALSE"  # IID spatiotemporal (years); with spatial rf; no covariate
# fit <- try(
#   sdmTMB(
#        formula = catch_weight ~ 1 + (1 | fyear),
#        family = tweedie(link = "log"),
#        data = arrow,
#        time = "year", spatiotemporal = "iid", spatial = "on",
#        mesh = mesh, 
#        #offset = log(arrow$area_swept / 100000), # why divided by this?
#        # need priors?; why priors used in simulation fit?
#        sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L))
#   ) 
# )

arrow %>% 
filter(catch_weight > 0  | catch_count > 0 | 
         density_kgkm2 > 0 | density_pcpm2 > 0 | density_ppkm2 > 0) %>% 
ggplot(aes(x = X, y = Y)) + 
  geom_point(aes(colour = log(density_kgkm2))) + 
  scale_colour_viridis_c() 



distinct(trawl, species_common_name) %>% view

trawl %>%
  filter(catch_weight > 0  | catch_count > 0 | 
         density_kgpm2 > 0 | density_pcpm2 > 0 | density_ppkm2 > 0) %>%  
  count(species_common_name, name = 'n_samps') %>% view

trawl %>% 
  filter(catch_weight > 0  | catch_count > 0 | 
         density_kgpm2 > 0 | density_pcpm2 > 0 | density_ppkm2 > 0) %>%  
filter(species_common_name == 'bocaccio') %>% view

trawl %>% 
filter(species_common_name == 'arrowtooth flounder') %>%
filter(catch_weight > 0  | catch_count > 0 | 
         density_kgpm2 > 0 | density_pcpm2 > 0 | density_ppkm2 > 0) %>% 
ggplot(aes(x = X, y = Y)) + 
  geom_point(aes(colour = log(catch_weight))) + 
  scale_colour_viridis_c() 

trawl_spp <- c('bocaccio')

synoptic_grid %>% glimpse

bocaccio_trawl <- trawl %>% filter(species_common_name == "bocaccio")


bocaccio_trawl %>% names



# Attemps to plot the grid out of curiosity
# guessing 'crs = 3156'
# syn_sf <- st_as_sf(synoptic_grid %>% mutate(X = X, Y = Y), coords = c('X', 'Y'), crs = st_crs(3156))

# syn_bbox <- st_bbox(syn_sf)
# syn_grid <- st_make_grid(syn_sf, cellsize = c(2, 2))

# ggplot(syn_grid) +
#   geom_sf() + 
#   geom_sf(data = syn_sf)
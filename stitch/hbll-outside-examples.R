library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here('stitch', '00_prep-example-data.R'))
source(here::here('stitch', 'utils.R'))
source(here::here('stitch', 'fit_model_func.R'))

mytheme <- function() ggsidekick::theme_sleek()  # sometimes I add more layers to themes
theme_set(mytheme())

# Try looking at quillback and yelloweye:
# ------------------------------------------------------------------------------
outside_test_dat <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  filter(species_common_name %in% c('quillback rockfish', 'yelloweye rockfish')) %>%
  #split(f = .$species_common_name) %>%
  split(f = list(.$species_common_name))
  # map(~.x %>% mutate(data_subset = paste(species_common_name, survey_abbrev, sep = "-")))

tictoc::tic()
#future::plan(future::multisession, workers = 5) # or whatever number
fits <- outside_test_dat %>%
  map(fit_models, catch = "density_ppkm2") %>%

  list_flatten(name_spec = "{inner}")
#future::plan(future::sequential)
tictoc::toc()
beep()

fits_cleaned <- fits %>%
  purrr::map(check_sanity)

# fits_cleaned <- fits %>%
#   purrr::map(check_sanity)

fitted_yrs <- sort(unique(outside_test_dat[[1]]$year))
#fitted_yrs_extra <- min(outside_test_dat[[1]]$year):max(outside_test_dat[[1]]$year)
fitted_yrs_extra <- sort(c(fitted_yrs, sdmTMB:::find_missing_time(outside_test_dat[[1]]$year)))

nd <-
  hbll_outside %>% # filter(survey == "HBLL OUT N") %>%
  make_grid(years = fitted_yrs) %>%
  mutate(fyear = as.factor(year))

nd_extra_time <-
  hbll_outside %>% # filter(survey == "HBLL OUT N") %>%
  make_grid(years = fitted_yrs_extra) %>%
  mutate(fyear = as.factor(year))

preds <- get_pred_list(fits_cleaned, newdata = nd, newdata_extra_time = nd_extra_time)
indices <- get_index_list(pred_list = preds)
beep()

index_df <- mk_index_df(indices)

p3 <-
ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 1L, scales = "free_y") +
  ggtitle("North")

# p1 / p2

# p3 / p4

# can fix the standard deviation parameter in the random walk to try and get the
# st time_varying RW; could base this off of a known a similar species like quillback or
# plus or minus of plus 0.2; using species life history characteristic, you can set a prior
# sd of 0.2 means that 95% of the jumps from year to year should be smaller than 20%
# which is probably generous for a species like yelloweye that lives ~120 years

fit_mle <- update(
  fit,
  control = sdmTMBcontrol(
    start = list(
      ln_tau_V = log(0.2)  # with only one random walk; 1 sd of that walk; that parameter is called tau_v; or ln_tau_v is the log version of that
    ),
    map = list(
      ln_tau_V = factor(NA)  # map with the NA says start and stay at 0.2
    )
  )
)

# Go with the counts and use log(hook_count) as the numbers

# Examples for the talk:
# - try and get the veraious configurations working with yelloweye;
# - this one looks like the
# - include an AR(1) model for comparison
# - change colours to match N/S
# - if the random fields :
# - time_varying random walk with AR(1) random fields (`"spatiotemporal = "iid"`)
# - try using the zero-centred random walk ('rw0'); make sure that you have a 1;
# - Drop 2021; drop 2020 (since no data);

# Look into the extra_time and send an example (on slack)


# If you have two different surveys and years;
# if you have years of factor and survey as factor; you need to have rw start at mean of zero with first year as deviation; otherwise overparameterised
#

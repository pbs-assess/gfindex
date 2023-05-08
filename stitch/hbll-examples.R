library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)

source(here::here('stitch', '00_prep-example-data.R'))
source(here::here('stitch', 'utils.R'))
source(here::here('stitch', 'fit_model_func.R'))

mytheme <- function() ggsidekick::theme_sleek()  # sometimes I add more layers to themes
theme_set(mytheme())

# family <- nbinom2(link = 'log') # could use this instead, but would need to change units and add offset

# Try looking at quillback and yelloweye:
# ------------------------------------------------------------------------------
inside_test_dat <- dat %>% 
  filter(str_detect(survey_abbrev, "HBLL INS")) %>% 
  filter(species_common_name %in% c('quillback rockfish', 'yelloweye rockfish')) %>% 
  split(., f = .$species_common_name)

fits <- inside_test_dat %>% 
  map(., fit_models, catch = "density_ppkm2") %>% 
  list_flatten(name_spec = "{inner}")
beep()

fits_cleaned <- fits %>%
  modify(., check_sanity)  # omit plots made from models that did not pass sanity check

fitted_yrs <- sort(unique(inside_test_dat[[1]]$year))
fitted_yrs_extra <- min(inside_test_dat[[1]]$year):max(inside_test_dat[[1]]$year)
nd <- 
  make_grid(hbll_inside, years = fitted_yrs) %>% 
  mutate(fyear = as.factor(year))
nd_extra_time <- 
  make_grid(hbll_inside, years = fitted_yrs_extra) %>% 
  mutate(fyear = as.factor(year))

preds <- get_pred_list(fits_cleaned, newdata = nd, newdata_extra_time = nd_extra_time)
indices <- get_index_list(pred_list = preds)
beep()

index_df <- mk_index_df(indices)

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y")

# without cleaned fits 
preds <- get_pred_list(fits, newdata = nd, newdata_extra_time = nd_extra_time)
indices <- get_index_list(pred_list = preds)
beep()

index_df <- mk_index_df(indices)

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_grid(species ~ fct_reorder(desc, order), scales = "free_y") + 
  expand_limits(y = 0)

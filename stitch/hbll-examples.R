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

names(fits)
fits[[1]]

fits_cleaned <- fits
# fits_cleaned <- fits %>%
#   modify(., check_sanity)  # omit plots made from models that did not pass sanity check
# beep()

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

index_df <-   enframe(indices) %>% 
  unnest(col = "value") %>% 
  separate(col = 'name', into = c('id', 'species', 'model'), sep = ":") %>% 
  mutate(id = as.numeric(id)) %>% 
  right_join(., model_lookup) %>% 
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd"))

select(index_df, id, model, desc, order) %>% distinct

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y")

# without cleaned fits 
preds <- get_pred_list(fits, newdata = nd)
indices <- get_index_list(pred_list = preds, area = nd$area)
beep()

index_df <- mk_index_df(indices)

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_grid(species ~ fct_reorder(desc, order), scales = "free_y") + 
  expand_limits(y = 0)
beep()

# ------------------------------------------------------------------------------
# Investigate why the index values drop at the end. They do not look like any 
# other plots I have seen elsewhere. 
qb <- dat %>% 
 filter(str_detect(survey_abbrev, "HBLL INS"), species_common_name == "quillback rockfish")

# Look at raw data
ggplot(mapping = aes(x = X, y = Y)) + 
  geom_point(data = qb %>% filter(catch_count != 0), aes(size = catch_weight, colour = catch_count), shape = 21) + 
  geom_point(data = filter(qb, catch_weight == 0), colour = "grey", shape = 21) + 
  facet_wrap(~ year)

qb_mesh <- make_mesh(qb, xy_cols = c('X', 'Y'), cutoff = 20)

qb_test_fit <- sdmTMB(
      density_ppkm2 ~ 1,
      data = qb,
      mesh = qb_mesh,
      time = "year",
      extra_time = sdmTMB:::find_missing_time(qb$year),
      family = tweedie(link = "log"),
      spatial = "on",
      spatiotemporal = "RW",
      control = sdmTMBcontrol(newton_loops = 1L)
)

nd <- 
  make_grid(hbll_inside, years = min(qb$year):max(qb$year)) %>% 
  #make_grid(hbll_inside, years = sort(unique(qb$year))) %>% 
  mutate(fyear = as.factor(year)) %>% 
  select(X, Y, area, year)

qb_test_p <- predict(qb_test_fit, newdata = nd, return_tmb_object = TRUE) 
  # predict(qb_test_fit) %>% 
  #mutate(across(c(X, Y), \(x) round(x, digits = 0)))

qb_test_p$data %>% glimpse()

ggplot(qb_test_p$data, aes(X, Y, colour = est + omega_s + epsilon_st)) +
  geom_point() +
  facet_wrap(~year) +
  coord_fixed()

qb_test_i <- get_index(qb_test_p, bias_correct = TRUE, area = nd$area)
beep()
ggplot(data = qb_test_i, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_ribbon(alpha = 0.20, colour = NA)


ggplot(data = index_df %>% filter(order == 6, species == "quillback rockfish"), aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_ribbon(alpha = 0.20, colour = NA)

   +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region")
beep()

fitted_yrs <- sort(unique(qb$year))
nd <- 
  make_grid(hbll_inside, years = fitted_yrs) %>% 
  mutate(fyear = as.factor(year))

fit_weight <- fit_models(qb, catch = 'catch_weight')
fit_count  <- fit_models(qb, catch = 'catch_count')
beep()

pred_wt <- get_pred_list(fit_weight, newdata = nd)
index_wt <- get_index_list(pred_wt, area = nd$area)

pred_ct <- get_pred_list(fit_count, newdata = nd)
index_ct <- get_index_list(pred_ct, area = nd$area)

qb_wt_ct_index_list <- list("catch_weight" = index_wt, "catch_count" = index_ct)

index_df <- bind_rows(mk_index_df(index_wt) %>% bind_cols(catch = "catch_weight"), 
  mk_index_df(index_ct) %>% bind_cols(catch = "catch_count")
)


ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_wrap(catch ~ desc, nrow = 2L, scales = "free_y")

pred_1an <- predict(fit1_an, newdata = nd, return_tmb_object = TRUE)
qb_index_1an <- get_index(pred_1an, bias_correct = TRUE, area = nd$area)



qb_mesh <- make_mesh(qb, c("X", "Y"), n_knots = 250)

fitted_yrs <- sort(unique(qb$year))
nd <- 
  make_grid(hbll_inside, years = fitted_yrs) %>% 
  mutate(fyear = as.factor(year))

fit1 <- try(
  sdmTMB(
    catch_count ~ 1,
    family = tweedie(),
    data = qb, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = qb_mesh,
    offset = NULL
  )
)

pred_1 <- predict(fit1, newdata = nd, return_tmb_object = TRUE)
qb_index1 <- get_index(pred_1, bias_correct = TRUE, area = nd$area)

fit1_an <- try(
  sdmTMB(
    catch_count ~ 1,
    family = tweedie(),
    data = qb, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = qb_mesh,
    offset = NULL, 
    anisotropy = TRUE
  )
)
beep()

pred_1an <- predict(fit1_an, newdata = nd, , return_tmb_object = TRUE)
qb_index_1an <- get_index(pred_1an, bias_correct = TRUE, area = nd$area)

qb_index_list <- list("1 no aniso:quillback" = qb_index1, "1 aniso:quillback" = qb_index_1an)

index_df <- qb_index_list %>% 
  enframe() %>%
  unnest(col = "value") %>% 
  separate(col = 'name', into = c('id', 'species'), sep = ":") %>% 
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd"))

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = odd_even)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_wrap(~ id, nrow = 1L, scales = "free_y")

beep()

filter(qb, is.na(year))

qb %>% view

qb_mesh <- make_mesh(data = qb, xy_cols = c("X", "Y"), cutoff = 20)

plot(qb_mesh)

qb_fits <- fit_models(dat = qb, mesh = qb_mesh, family = tweedie())
beep()




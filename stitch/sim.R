# Simulation testing survey stitching with various models -----------------

# Observations ------------------------------------------------------------

# - Try lower or increasing observation error ('phi') (e.g., 2 vs. 7)
# - Try lowering or increasing numbers of samples (e.g., 250 vs. 400)

PHI <- 7
SAMPLE_N <- 400
SEED <- 123

# Simulate data -----------------------------------------------------------

predictor_dat <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
  year = 1:10
)
mesh_sim <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

# Define north and south regions ------------------------------------------

predictor_dat$region <- NA
predictor_dat$region[predictor_dat$Y < 0.5] <- "south"
predictor_dat$region[predictor_dat$Y >= 0.5] <- "north"
predictor_dat$region <- as.factor(predictor_dat$region)

set.seed(SEED * 889)
yrs <- rnorm(max(predictor_dat$year), mean = 1, sd = 0.2)
sim_dat <- sdmTMB_simulate(
  formula = ~ 0 + as.factor(year) + region,
  data = predictor_dat,
  time = "year",
  mesh = mesh_sim,
  family = tweedie(),
  range = 0.5,
  sigma_E = 0.2,
  rho = 0.5,
  phi = PHI,
  tweedie_p = 1.5,
  sigma_O = 0.3,
  seed = SEED * 1029,
  B = c(yrs, 1.3) # last value is multiplicative value for south
)
sim_dat$region <- predictor_dat$region
sim_dat$year <- predictor_dat$year

# Visualize what we just did ----------------------------------------------

ggplot(sim_dat, aes(X, Y, fill = eta)) +
  geom_raster() +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c()

# Zoom in on a year -------------------------------------------------------

# filter(sim_dat, year == 1) %>%
#   ggplot(aes(X, Y, fill = eta)) +
#   geom_raster() +
#   facet_wrap(vars(year)) +
#   scale_fill_viridis_c()

# Sample 400 per year -----------------------------------------------------

set.seed(SEED * 9283)
obs_dat <- sim_dat %>%
  group_by(year) %>%
  sample_n(SAMPLE_N)

# Lose half the survey most years -----------------------------------------

d <- obs_dat
d <- d[!(d$year %in% seq(3, 9, 2) & d$region == "north"), ]
d <- d[!(d$year %in% seq(2, 10, 2) & d$region != "north"), ]
d$sampled_region <- d$region
d$sampled_region <- as.character(d$sampled_region)
d$sampled_region[d$year == 1] <- "both"

# Visualize it ------------------------------------------------------------

ggplot(d, aes(X, Y, colour = log(observed))) +
  geom_point() +
  facet_wrap(vars(year)) +
  scale_colour_viridis_c()

# Calculate known true biomass/abundance ----------------------------------

actual <- group_by(sim_dat, year) %>%
  summarise(total = sum(mu))

# Fit models --------------------------------------------------------------

mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.1)

priors <- sdmTMBpriors(
  matern_s = pc_matern(range_gt = 0.3, sigma_lt = 0.4),
  matern_st = pc_matern(range_gt = 0.3, sigma_lt = 0.3)
)

fit_rw <- sdmTMB(
  observed ~ 1 + as.factor(region), family = tweedie(),
  data = d, time = "year", spatiotemporal = "rw", spatial = "off",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_rw_no_region <- sdmTMB(
  observed ~ 1, family = tweedie(),
  data = d, time = "year", spatiotemporal = "rw", spatial = "off",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_iid <- sdmTMB(
  observed ~ 0 + as.factor(year) + as.factor(region), family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_iid_no_region <- sdmTMB(
  observed ~ 0 + as.factor(year), family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

# Predict on grid and calculate indexes -----------------------------------

nd <- select(predictor_dat, X, Y, year, region)

pred_rw <- predict(fit_rw, newdata = nd, return_tmb_object = TRUE)
geo_index_rw <- get_index(pred_rw, bias_correct = TRUE)

pred_rw_no_region <- predict(fit_rw_no_region, newdata = nd, return_tmb_object = TRUE)
geo_index_rw_no_region <- get_index(pred_rw_no_region, bias_correct = TRUE)

pred <- predict(fit_iid, newdata = nd, return_tmb_object = TRUE)
geo_index <- get_index(pred, bias_correct = TRUE)

pred_no_region <- predict(fit_iid_no_region, newdata = nd, return_tmb_object = TRUE)
geo_index_no_region <- get_index(pred_no_region, bias_correct = TRUE)

# Plot it -----------------------------------------------------------------

mult <- 1
g <- mutate(geo_index_rw, type = "RW", with_region = "region = TRUE") %>%
  bind_rows(mutate(geo_index, type = "IID", with_region = "region = TRUE")) %>%
  bind_rows(mutate(geo_index_no_region, type = "IID", with_region = "region = FALSE")) %>%
  bind_rows(mutate(geo_index_rw_no_region, type = "RW", with_region = "region = FALSE")) %>%
  left_join(select(d, year, sampled_region) %>% distinct()) %>%
  ggplot(aes(year, est / mult,
  ymin = lwr / mult, ymax = upr / mult)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(data = actual, mapping = aes(year, total),
    inherit.aes = FALSE, lty = 2) +
  facet_grid(with_region~type) +
  ggtitle("IID vs. RW model; with and without region",
  subtitle = paste0("Dashed = true; dots/lines = estimated\n", "phi = ", PHI, ", N = ", SAMPLE_N)) +
  ylab("Abundance estimate") + xlab("Year") +
  labs(colour = "Sampled region")
print(g)

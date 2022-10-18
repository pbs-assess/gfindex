# To look at:
# - gap matters! how much; related to range I assume
# - **why does RW 'fix' things!?**
# - getting the depth covariate right seems to fix seesaw; *but* bigger CIs on
#   poorly sampled region
# - is coverage of RW still OK!?
# - happens with spatial-only model I think (from Quang)
# - dig in... need to articulate simply what's going wrong!
# - how do I know this is happening in reality (besides the seesaw)?

library(sdmTMB)
library(ggplot2)
library(dplyr)
# Simulation testing survey stitching with various models -----------------

PHI <- 8
SAMPLE_N <- 300
SEED <- 12346

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
yrs <- rnorm(max(predictor_dat$year), mean = 1.2, sd = 0.15)
predictor_dat$depth_cov <- predictor_dat$Y
sim_dat <- sdmTMB_simulate(
  formula = ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2),
  data = predictor_dat,
  time = "year",
  mesh = mesh_sim,
  family = tweedie(),
  range = 0.8,
  sigma_E = 0.4,
  rho = 0.5,
  phi = PHI,
  tweedie_p = 1.6,
  sigma_O = 1.2,
  seed = SEED * 1029,
  B = c(yrs, 2.2, 3.8)
)
xx <- seq(0, 1, length.out = 100)
plot(xx, xx * 2.2 + 3.8 * xx^2, ylab = "Depth effect", xlab = "Depth value")
sim_dat$region <- predictor_dat$region
sim_dat$year <- predictor_dat$year

# Visualize what we just did ----------------------------------------------

ggplot(sim_dat, aes(X, Y, fill = eta)) +
  geom_raster() +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c()

# Zoom in on a year -------------------------------------------------------

filter(sim_dat, year == 1) %>%
  ggplot(aes(X, Y, fill = eta)) +
  geom_raster() +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c()

# Sample 400 per year -----------------------------------------------------

set.seed(SEED * 9283)
obs_dat <- sim_dat %>%
  group_by(year) %>%
  sample_n(SAMPLE_N)

# Lose half the survey most years -----------------------------------------

d <- obs_dat
# d <- d[!(d$year %in% seq(3, 9, 2) & d$region == "north"), ]
d <- d[!(d$year %in% seq(1, 9, 2) & d$region == "north"), ]
d <- d[!(d$year %in% seq(2, 10, 2) & d$region != "north"), ]
d$sampled_region <- d$region
d$sampled_region <- as.character(d$sampled_region)
# d$sampled_region[d$year == 1] <- "both"

# Remove strip in middle to increase gap?
d <- d[!(d$Y > 0.3 & d$Y < 0.7),]

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
  observed ~ 1 + depth_cov + I(depth_cov^2), family = tweedie(),
  data = d, time = "year", spatiotemporal = "rw", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_rw_no_cov <- sdmTMB(
  observed ~ 1, family = tweedie(),
  data = d, time = "year", spatiotemporal = "rw", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_iid <- sdmTMB(
  observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2), family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_iid_no_cov <- sdmTMB(
  observed ~ 0 + as.factor(year), family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_sp <- sdmTMB(
  observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2),
  family = tweedie(),
  data = d, time = "year", spatiotemporal = "off", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

fit_sp_no_cov <- sdmTMB(
  observed ~ 0 + as.factor(year),
  family = tweedie(),
  data = d, time = "year", spatiotemporal = "off", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)

# Predict on grid and calculate indexes -----------------------------------

nd <- select(predictor_dat, X, Y, year, region)
nd$depth_cov <- nd$Y

pred_rw <- predict(fit_rw, newdata = nd, return_tmb_object = TRUE)
geo_index_rw <- get_index(pred_rw, bias_correct = TRUE)

pred_rw_no_cov <- predict(fit_rw_no_cov, newdata = nd, return_tmb_object = TRUE)
geo_index_rw_no_cov <- get_index(pred_rw_no_cov, bias_correct = TRUE)

pred <- predict(fit_iid, newdata = nd, return_tmb_object = TRUE)
geo_index <- get_index(pred, bias_correct = TRUE)

pred_no_cov <- predict(fit_iid_no_cov, newdata = nd, return_tmb_object = TRUE)
geo_index_no_cov <- get_index(pred_no_cov, bias_correct = TRUE)

pred_sp <- predict(fit_sp, newdata = nd, return_tmb_object = TRUE)
geo_index_sp <- get_index(pred_sp, bias_correct = TRUE)

pred_sp_no_cov <- predict(fit_sp_no_cov, newdata = nd, return_tmb_object = TRUE)
geo_index_sp_no_cov <- get_index(pred_sp_no_cov, bias_correct = TRUE)

# Plot it -----------------------------------------------------------------

mult <- 1
g <- mutate(geo_index_rw, type = "RW", with_depth = "depth = TRUE") %>%
  bind_rows(mutate(geo_index, type = "IID", with_depth = "depth = TRUE")) %>%
  bind_rows(mutate(geo_index_no_cov, type = "IID", with_depth = "depth = FALSE")) %>%
  bind_rows(mutate(geo_index_rw_no_cov, type = "RW", with_depth = "depth = FALSE")) %>%
  bind_rows(mutate(geo_index_sp, type = "Spatial only", with_depth = "depth = TRUE")) %>%
  bind_rows(mutate(geo_index_sp_no_cov, type = "Spatial only", with_depth = "depth = FALSE")) %>%
  left_join(select(d, year, sampled_region) %>% distinct()) %>%
  ggplot(aes(year, est / mult,
  ymin = lwr / mult, ymax = upr / mult)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(data = actual, mapping = aes(year, total),
    inherit.aes = FALSE, lty = 2) +
  facet_grid(with_depth~type) +
  ggtitle("IID vs. RW model; with and without cov",
  subtitle = paste0("Dashed = true; dots/lines = estimated\n", "phi = ", PHI, ", N = ", SAMPLE_N)) +
  ylab("Abundance estimate") + xlab("Year") +
  labs(colour = "Sampled region")
print(g)

# look at one point in space...

get_eg_cell <- function(obj, x, y) {
  obj$data[
    round(obj$data$X, 3) == round(x, 3) &
    round(obj$data$Y, 3) == round(y, 3), ]
}

objs <- list(pred_rw, pred_rw_no_cov, pred, pred_no_cov, pred_sp, pred_sp_no_cov)
names(objs) <- c("RW depth", "RW", "IID depth", "IID", "Spatial depth", "Spatial")

p1 <- purrr::map_dfr(objs, get_eg_cell, x = 0.50505051, y = 0.81818182,
  .id = "model") |>
  mutate(with_depth = grepl("depth", model)) |>
  mutate(type = gsub(" depth", "", model))

ggplot(p1, aes(year, est)) + geom_line() +
  facet_grid(with_depth~type)

# Question: can an extra latent year at the beginning fix the 'burnin' problem!?


# Observations ------------------------------------------------------------

# - Try lower or increasing observation error ('phi') (e.g., 2 vs. 7)
# - Try lowering or increasing numbers of samples (e.g., 250 vs. 400)

# - Is it important to have at least one year with full coverage?
# - If so, does it matter if it's near the beginning or end?
# - Is it/when is it important to have some kind of fixed effect covariate (region/depth)?
# - What form (IID, RW, AR1) works best in what conditions?
# - Does an AR1 see-saw towards mean with uneven coverage?
# - When is RW at risk of over smoothing process?

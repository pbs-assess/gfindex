# question: same thing happen with missing chunks from one index!?

# To look at:
# - gap matters! how much; related to range I assume
# - **why does RW 'fix' things!?**
# - getting the depth covariate right seems to fix seesaw; *but* bigger CIs on
#   poorly sampled region
# - is coverage of RW still OK!?
# - happens with spatial-only model I think (from Quang)
# - dig in... need to articulate simply what's going wrong!
# - how do I know this is happening in reality (besides the seesaw)?
# - does including a region/north/south covariate fix it: No! build in to
#   routine code?

# What is an index of seesawness?
#
# Things to show:
# - spatial random field messing up
# - RW deviations adapting to fit?
# - fix depth covariates correctly - works?


# quadratic, linear, or breakpoint covariate effect?

library(sdmTMB)
library(ggplot2)
library(dplyr)
source("stitch/funcs.R")
dir.create("stitch/figs", showWarnings = FALSE)
# Simulation testing survey stitching with various models -----------------

# n_year <- 12 # even number
# .seed <- 28817
# .seed <- 288171

# predictor_grid <- expand.grid(
#   X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
#   year = seq_len(n_year)
# )
# mesh_sim <- make_mesh(predictor_grid, xy_cols = c("X", "Y"), cutoff = 0.1)

# xx <- seq(0, 1, length.out = 100)
# plot(xx, xx * 2 + 5 * xx^2, ylab = "Depth effect", xlab = "Depth value")

sim_fit_and_index <- function(n_year,
                              .seed,
                              gap_size = 0.3,
                              obs_sampled_size = 400L,
                              year_marginal_sd = 0.3,
                              obs_yrs = list(
                                north_yrs = seq(1, n_year - 1, 2),
                                south_yrs = seq(2, n_year, 2)
                              ),
                              phi = 8,
                              region_cutoff = 0.5,
                              range = 0.8,
                              sigma_O = 1.6,
                              year_arima.sim = list(ar = 0.6),
                              sim_coefs = c(2, 5)) {
  is_even <- function(x) x %% 2 == 0
  if (!is_even(n_year)) cli::cli_abort("Number of years must be even.")

  cli::cli_alert_info(glue::glue("Running simulation for seed ", .seed))

  cli::cli_alert_success("Creating mesh...")
  predictor_grid <- expand.grid(
    X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
    year = seq_len(n_year)
  )
  mesh_sim <- make_mesh(predictor_grid, xy_cols = c("X", "Y"), cutoff = 0.1)

  cli::cli_alert_success("Simulating...")
  x <- sim(predictor_grid, mesh_sim,
    seed = .seed, phi = phi, range = range,
    region_cutoff = region_cutoff,
    year_arima.sim = year_arima.sim, year_marginal_sd = year_marginal_sd,
    coefs = sim_coefs,
    north_effect = 0, sigma_O = sigma_O
  )
  sim_dat <- x$sim_dat
  predictor_dat <- x$predictor_dat

  # Visualize what we just did ----------------------------------------------

  # # Year effects:
  # ggplot(data.frame(x = seq_len(n_year), y = x$year_effects), aes(x, y)) +
  #   geom_line()
  # ggplot(data.frame(x = seq_len(n_year), y = x$year_effects), aes(x, exp(y))) +
  #   geom_line()
  #
  # blank_theme_elements <- theme(panel.grid.major = element_line(colour = "grey90"),
  #   # panel.spacing.x = unit(20, "pt"),
  #   axis.text = element_blank(), axis.ticks = element_blank(),
  #   axis.title = element_blank(), legend.position = "right")
  #
  # # Spatio-temporal truth
  # ggplot(sim_dat, aes(X, Y, fill = eta)) +
  #   geom_raster() +
  #   facet_wrap(vars(year)) +
  #   scale_fill_viridis_c() +
  #   ggsidekick::theme_sleek() +
  #   coord_equal(expand = FALSE) +
  #   blank_theme_elements +
  #   labs(fill = "True\nsimulated\nlog abundance")
  # # ggsave("stitch/figs/spatio-temporal-truth.pdf", width = 8, height = 5)

  # Sample N per year -----------------------------------------------------

  cli::cli_alert_success("Observing...")
  d <- observe(
    sim_dat,
    sample_n = obs_sampled_size,
    region_cutoff = region_cutoff,
    seed = .seed,
    north_yrs = obs_yrs$north_yrs,
    south_yrs = obs_yrs$south_yrs, gap = gap_size
  )

  # Visualize it ------------------------------------------------------------

  # ggplot(d, aes(X, Y, colour = log(observed))) +
  #   geom_point() +
  #   facet_wrap(vars(year)) +
  #   scale_colour_viridis_c() +
  #   ggsidekick::theme_sleek() +
  #   blank_theme_elements +
  #   coord_equal(expand = FALSE) +
  #   theme(panel.grid.major = element_line(colour = "grey90"))
  # # panel.spacing.x = unit(20, "pt")
  # # ggsave("stitch/figs/spatio-temporal-observed.pdf", width = 8, height = 5)

  # Calculate known true biomass/abundance ----------------------------------

  regions_samp <- select(d, year, sampled_region) %>% distinct()

  actual <- group_by(sim_dat, year) %>%
    summarise(total = sum(mu)) |>
    left_join(regions_samp, by = "year")

  # ggplot(actual, aes(year, total)) +
  #   geom_line()

  # Fit models --------------------------------------------------------------

  cli::cli_alert_success("Fitting models...")
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.1)
  priors <- sdmTMBpriors(
    matern_s = pc_matern(range_gt = 0.3, sigma_lt = 0.4),
    matern_st = pc_matern(range_gt = 0.3, sigma_lt = 0.3)
  )

  fits <- list()
  nms <- c()
  i <- 1

  fits[[i]] <- sdmTMB(
    observed ~ 1,
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors
  )
  nms <- c(nms, "RW")
  i <- i + 1

  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors
  )
  nms <- c(nms, "IID covariate")
  i <- i + 1

  fits[[i]] <- sdmTMB(
    observed ~ s(year),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors
  )
  nms <- c(nms, "IID s(year)")
  i <- i + 1

  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = mesh,
    priors = priors
  )
  i <- i + 1
  nms <- c(nms, "IID")

  fits[[i]] <- sdmTMB(
    observed ~ 0,
    family = tweedie(),
    time_varying = ~1,
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = mesh,
    priors = priors
  )
  i <- i + 1
  nms <- c(nms, "IID RW year")

  fits[[i]] <- sdmTMB(
    observed ~ 0,
    family = tweedie(),
    time_varying = ~1, time_varying_type = "ar1",
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = mesh,
    priors = priors
  )
  i <- i + 1
  nms <- c(nms, "IID AR1 year")

  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "off", spatial = "on",
    mesh = mesh,
    priors = priors
  )
  i <- i + 1
  nms <- c(nms, "Spatial only")

  names(fits) <- nms

  # Predict on grid and calculate indexes -----------------------------------

  cli::cli_alert_success("Calculating indices...")
  nd <- select(predictor_dat, X, Y, year, region)
  nd$depth_cov <- nd$Y
  preds <- purrr::map(fits, predict, newdata = nd, return_tmb_object = TRUE)
  indexes <- purrr::map(preds, get_index, bias_correct = TRUE)

  indexes_df <- dplyr::bind_rows(indexes, .id = "model") |>
    mutate(with_depth = paste0("covariate = ", grepl("covariate", model))) |>
    mutate(type = gsub(" covariate", "", model))

  indexes_df <- left_join(indexes_df, actual, by = "year") |>
    mutate(seed = .seed)
  indexes_df
}

# Plot it -----------------------------------------------------------------

# out <- purrr::map_dfr(seq_len(1L), ~ sim_fit_and_index(n_year = 12L, .seed = .x))
# out <- purrr::map_dfr(seq_len(10L), ~ sim_fit_and_index(n_year = 12L, .seed = .x))
#

# study design:
# - [ ] no gap, small gap, large gap `gap_size`
#       gap_size = 0.00
#       gap_size = 0.25
#       gap_size = 0.50

# - [ ] no covariate effect; strong covariate effect
#       sim_coefs = c(2, 5)
#       sim_coefs = c(0, 0)

# - [ ] lower or increasing observation error ('phi')
#       phi = 2
#       phi = 7

# - [ ] try lowering or increasing numbers of samples
#       obs_sampled_size = 200L
#       obs_sampled_size = 400L

# - [ ] one year of complete coverage at beginning or end or none TODO
#       obs_yrs = list(north_yrs = seq(1, n_year - 1, 2), south_yrs = c(seq(2, n_year, 2))
#       obs_yrs = list(north_yrs = c(seq(1, n_year - 1, 2), n_year), south_yrs = c(seq(2, n_year, 2))

# - [ ] impact of uneven north south vs. even?
#       region_cutoff = 0.50
#       region_cutoff = 0.25

n_year <- 12L

base <- list(
  n_year = n_year,
  label = "empty",
  gap_size = 0.25,
  sim_coefs = c(2, 5),
  phi = 8,
  obs_sampled_size = 400L,
  obs_yrs = list(
    north_yrs = seq(1, n_year - 1, 2), south_yrs = c(seq(2, n_year, 2)),
    region_cutoff = 0.50,
    range = 0.8,
    sigma_O = 1.6,
    year_marginal_sd = 0.2
  )
)

sc <- purrr::map(seq_len(14), ~base)

i <- 1

sc[[i]]$label <- "Base"
i <- i + 1

sc[[i]]$label <- "No gap"
sc[[i]]$gap_size <- 0
i <- i + 1

sc[[i]]$label <- "Large gap"
sc[[i]]$gap_size <- 0.5
i <- i + 1

sc[[i]]$label <- "No covariate"
sc[[i]]$sim_coefs <- c(0, 0)
i <- i + 1

sc[[i]]$label <- "Low observation error"
sc[[i]]$phi <- 3
i <- i + 1

sc[[i]]$label <- "High observation error"
sc[[i]]$phi <- 14
i <- i + 1

sc[[i]]$label <- "Low sample size"
sc[[i]]$obs_sampled_size <- 200L
i <- i + 1

sc[[i]]$label <- "High sample size"
sc[[i]]$obs_sampled_size <- 800L
i <- i + 1

sc[[i]]$label <- "Year of overlap"
sc[[i]]$obs_yrs <- list(north_yrs = c(seq(1, n_year - 1, 2), n_year), south_yrs = c(seq(2, n_year, 2)))
i <- i + 1

sc[[i]]$label <- "Unequal regions"
sc[[i]]$region_cutoff <- 0.25
i <- i + 1

sc[[i]]$label <- "Low range"
sc[[i]]$range <- 0.2
i <- i + 1

sc[[i]]$label <- "Low sigma O"
sc[[i]]$sigma_O <- 0.2
i <- i + 1

sc[[i]]$label <- "High year SD"
sc[[i]]$year_marginal_sd <- 0.5
i <- i + 1

sc[[i]]$label <- "Low year SD"
sc[[i]]$year_marginal_sd <- 0.05
i <- i + 1

if (any(grepl("empty", purrr::map_chr(sc, "label")))) {
  stop("Too many slots")
}

names(sc) <- purrr::map_chr(sc, "label")
sc <- purrr::map(sc, ~ {
  .x$label <- NULL
  .x
})

# out <- purrr::pmap_dfr(sc, sim_fit_and_index, .id = "label")

seeds <- seq_len(20L)
out_df <- purrr:::map_dfr(seeds, function(seed_i) {
  for (i in seq_along(sc)) {
    sc[[i]]$.seed <- seed_i
    x[[i]] <- do.call(sim_fit_and_index, sc[[i]])
  }
  names(x) <- names(sc)
  bind_rows(x, .id = "label")
})

# ---------------------------------------------
# iterate above ....

actual <- select(out_df, label, year, total, seed, sampled_region) |>
  distinct()

# out_df$label <- forcats::fct_inorder(out_df$label)

cols <- RColorBrewer::brewer.pal(3L, "Set2")
names(cols) <- c("north", "south", "both")
actual2 <- mutate(actual, label = gsub("obs", "\\\nobs", label))

g <- out_df |>
  mutate(with_depth = gsub("covariate =", "cov =", with_depth)) |>
  mutate(label = gsub("obs", "\\\nobs", label)) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  # geom_line(colour = "grey50") +
  geom_ribbon(alpha = 0.20, colour = NA) +
  geom_line(
    data = actual2, mapping = aes(year, total),
    inherit.aes = FALSE, lty = 2
  ) +
  facet_grid(forcats::fct_inorder(label) ~ paste(type, with_depth),
    scales = "free_y"
  ) +
  ylab("Abundance estimate") +
  xlab("Year") +
  labs(colour = "Sampled region") +
  scale_colour_manual(values = cols[c(2, 1, 3)]) +
  scale_y_log10() +
  scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2))
print(g)
ggsave("stitch/figs/saw-tooth-scenarios.pdf", width = 14, height = 17)

# Look at one point in space... -------------------------------------------

# get_eg_cell <- function(obj, x, y) {
#   obj$data[
#     round(obj$data$X, 3) == round(x, 3) &
#       round(obj$data$Y, 3) == round(y, 3),
#   ]
# }
#
# p1 <- purrr::map_dfr(preds, get_eg_cell,
#   x = 0.50505051, y = 0.81818182,
#   .id = "model"
# ) |>
#   mutate(with_depth = grepl("covariate", model)) |>
#   mutate(type = gsub(" covariate", "", model))
#
# p1 |>
#   left_join(select(d, year, sampled_region) %>% distinct()) |>
#   ggplot(aes(year, est)) +
#   geom_line() +
#   facet_grid(with_depth ~ type) +
#   geom_point(aes(colour = sampled_region)) +
#   ggsidekick::theme_sleek()

# What about MRE, RMSE, see-saw, coverage etc. ? --------------------------

out |>
  # left_join(actual) |>
  # left_join(select(d, year, sampled_region) %>% distinct()) |>
  group_by(model) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[sampled_region == "north"]) -
      mean(log_residual[sampled_region == "south"])),
    mre = mean(log_residual),
    rmse = sqrt(mean(log_residual^2)),
    mean_se = mean(se),
    coverage = mean(total < upr & total > lwr)
  ) |>
  arrange(seesaw_index, rmse) |>
  mutate(model = forcats::fct_reorder(model, rev(seesaw_index))) |>
  tidyr::pivot_longer(cols = -model, names_to = "metric") |>
  mutate(metric = factor(metric,
    levels = c("seesaw_index", "rmse", "mre", "mean_se", "coverage")
  )) |>
  ggplot(aes(value, model)) +
  geom_point(pch = 21, size = 1.6) +
  facet_wrap(~metric, scales = "free_x", nrow = 1L) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y.left = element_blank()) +
  xlab("Metric value")

# observation: even with covariate correctly fixed, seems to want to put
# variance into year factors and shrink random fields
# nothing stopping it from see-sawing fixed effects and shrinking random field
# variance

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

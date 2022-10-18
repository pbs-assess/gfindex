sim <- function(predictor_grid, mesh, phi = 8, seed = 123,
  region_cutoff = 0.5, rho = 0, sigma_E = 0.4, range = 0.8, tweedie_p = 1.6,
  sigma_O = 1.6, coefs = c(2.2, 3.8), year_mean = 1.2, year_sd = 0.15
) {
  SEED <- seed
  predictor_grid$region <- NA
  predictor_grid$region[predictor_grid$Y < region_cutoff] <- "south"
  predictor_grid$region[predictor_grid$Y >= region_cutoff] <- "north"
  predictor_grid$region <- as.factor(predictor_grid$region)
  set.seed(SEED * 889)
  yrs <- rnorm(max(predictor_grid$year), mean = year_mean, sd = year_sd)
  predictor_grid$depth_cov <- predictor_grid$Y
  sim_dat <- sdmTMB_simulate(
    formula = ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2),
    data = predictor_grid,
    time = "year",
    mesh = mesh_sim,
    family = tweedie(),
    range = range,
    sigma_E = sigma_E,
    rho = rho,
    phi = phi,
    tweedie_p = tweedie_p,
    sigma_O = sigma_O,
    seed = SEED * 1029,
    B = c(yrs, coefs[1], coefs[2])
  )
  sim_dat$region <- predictor_grid$region
  sim_dat$year <- predictor_grid$year
  list(sim_dat = sim_dat, predictor_dat = predictor_grid)
}

observe <- function(sim_dat, sample_n = 300L, seed = 10282, gap = 0,
  north_yrs = seq(1, 9, 2), south_yrs = seq(2, 10, 2), both_yrs = NULL) {

  assertthat::assert_that(gap >= 0 && gap <= 0.99)
  assertthat::assert_that(sample_n > 1L)
  assertthat::assert_that(all(north_yrs %in% sim_dat$year) &&
      all(south_yrs %in% sim_dat$year))

  set.seed(seed * 9283)
  d <- sim_dat %>%
    group_by(year) %>%
    sample_n(sample_n)
  # Lose half the survey most years -----------------------------------------
  d <- d
  # d <- d[!(d$year %in% seq(3, 9, 2) & d$region == "north"), ]

  north <- d[(d$year %in% union(north_yrs, both_yrs) & d$region == "north"), ]
  south <- d[(d$year %in% union(south_yrs, both_yrs) & d$region == "south"), ]

  d <- bind_rows(north, south)
  d$sampled_region <- as.character(d$region)
  d$sampled_region[d$year %in% both_yrs] <- "both"
  # Remove strip in middle to increase gap?
  d <- d[!(d$Y > (0.5 - gap/2) & d$Y < (0.5 + gap/2)),]
}

library(dplyr)
library(ggplot2)
library(here)
library(sdmTMB)

region <- c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG")
# region <- c("SYN WCVI")
species <- "arrowtooth flounder"

dat <- readRDS(here("data/all-survey-sets-2021-select.rds")) %>%
  dplyr::filter(species_common_name == tolower(species)) %>%
  dplyr::filter(survey_abbrev %in% region)
# change from per m2 to per km2:
dat$density <- dat[["density_kgpm2"]] * 1000000
dat$log_depth <- log(dat$depth_m)
dat$area_swept1 <- dat$doorspread_m * (dat$speed_mpm * dat$duration_min)
dat$area_swept2 <- dat$tow_length_m * dat$doorspread_m
dat$area_swept <- ifelse(!is.na(dat$area_swept2), dat$area_swept2, dat$area_swept1)
dat <- dplyr::filter(dat, !is.na(area_swept))
dat <- sdmTMB::add_utm_columns(dat, c("longitude", "latitude"), utm_crs = 32609)

formula <- catch_weight ~ s(log_depth, k = 5)
# formula <- catch_weight ~ 0 + as.factor(year)

set_priors <- sdmTMBpriors(
  matern_s = pc_matern(range_gt = 10, sigma_lt = 4),
  matern_st = pc_matern(range_gt = 10, sigma_lt = 3)
)

# set_priors <- sdmTMBpriors(
#   matern_s = pc_matern(range_gt = 5, sigma_lt = 3),
#   matern_st = pc_matern(range_gt = 5, sigma_lt = 3)
# )

spde <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 20)

# bnd <- INLA::inla.nonconvex.hull(cbind(dat$X, dat$Y), convex = -0.04)
# inla_mesh <- INLA::inla.mesh.2d(
#   boundary = bnd,
#   max.edge = c(50, 120),
#   offset = c(-0.01, -0.02),
#   cutoff = c(20, 100),
#   min.angle = c(21, 21)
# )
# mesh <- make_mesh(dat, c("X", "Y"), mesh = inla_mesh)
# plot(mesh$mesh, asp = 1)
# points(dat$X, dat$Y, pch = ".")
# mesh$mesh$n

ctrl <- sdmTMBcontrol(newton_loops = 2L)

mesh <- make_mesh(dat, c("X", "Y"), cutoff = 20)
plot(mesh$mesh, asp = 1)
points(dat$X, dat$Y, pch = ".")
mesh$mesh$n

# yrs <- unique(dat$year)
# all_yrs <- seq(min(dat$year), max(dat$year))
# base::setdiff(yrs, all_yrs)

nd <- readRDS(here("grids/synoptic_grid.rds")) %>%
  dplyr::filter(survey %in% region)
make_grid <- function(.x, years) {
  years <- sort(unique(years))
  .nd <- do.call(
    "rbind",
    replicate(length(years), .x, simplify = FALSE)
  )
  .nd$year <- rep(years, each = nrow(.x))
  .nd
}
fitted_yrs <- sort(unique(dat$year))
nd <- make_grid(nd, years = fitted_yrs)
nd <- na.omit(nd)
nd$year <- as.integer(nd$year)
nd$log_depth <- log(nd$depth)
nd$cell_area <- 4e+6

TMB::openmp(n = 3, DLL = "sdmTMB")
fit <- try(
  sdmTMB(
    formula,
    dat,
    mesh = spde,
    time = "year",
    family = delta_gamma(),
    # family = tweedie(),
    spatial = "on",
    spatiotemporal = "rw",
    offset = log(dat$area_swept / 100000),
    share_range = TRUE,
    anisotropy = FALSE,
    silent = FALSE,
    control = ctrl,
    do_index = TRUE,
    index_args = list(area = nd$cell_area / 100000),
    predict_args = list(newdata = nd, return_tmb_object = TRUE),
  )
)

s <- sanity(fit)
fit

# p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
# ind <- get_index(p, bias_correct = TRUE, area = nd$cell_area / 100000)
ind <- get_index(fit, bias_correct = TRUE)

mult <- 1000
g <- ggplot(ind, aes(year, est / mult)) +
  ggsidekick::theme_sleek()
g <- g + geom_pointrange(aes(ymin = lwr / mult, ymax = upr / mult),
  alpha = 1, colour = "grey20"
)
g <- g + ylab("Biomass (tonnes)") +
  ggtitle(stringr::str_to_title(species), subtitle = paste(region, collapse = ", ")) +
  coord_cartesian(
    ylim = c(0, max(ind$upr / mult) * 1.03),
    expand = FALSE, xlim = range(ind$year) + c(-0.25, 0.25)
  ) +
  theme(axis.title.x = element_blank())
g

# g <- ggplot(ind, aes(year, est / 1000)) +
#   ggsidekick::theme_sleek()
# g <- g +
#   geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000),
#     alpha = 0.9, fill = "grey90"
#   ) +
#   geom_line(alpha = 0.9, lty = 1, lwd = 0.5)
# g

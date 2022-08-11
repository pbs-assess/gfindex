library(dplyr)
library(ggplot2)
library(here)
library(sdmTMB)

source("stitch/utils.R")

dat <- prep_data()
list_species <- read_species_list("stitch/spp.txt")

out <- purrr::map(list_species, ~
    fit_index(
      dat,
      species = .x
    )
) %>% setNames(list_species)
index <- purrr::map_dfr(out, "index", .id = "species")

out_nodepth <- purrr::map(list_species, ~
    fit_index(
      dat,
      species = .x,
      formula = catch_weight ~ 1,
    )
) %>% setNames(list_species)
index_nodepth <- purrr::map_dfr(out_nodepth, "index", .id = "species")

# out_ar1 <- purrr::map(list_species, ~
#     fit_index(
#       dat,
#       species = .x,
#       spatiotemporal = "ar1"
#     )
# ) %>% setNames(list_species)
#
# out_iid <- purrr::map(list_species, ~
#     fit_index(
#       dat,
#       species = .x,
#       formula = catch_weight ~ 0 + as.factor(year),
#       spatiotemporal = "iid"
#     )
# ) %>% setNames(list_species)

out_dg_nodepth <- purrr::map(list_species, ~
    fit_index(
      dat,
      species = .x,
      formula = catch_weight ~ 1,
      family = delta_gamma()
    )
) %>% setNames(list_species)
index_dg_nodepth <- purrr::map_dfr(out_dg_nodepth, "index", .id = "species")

out_dg <- purrr::map(list_species, ~
    fit_index(
      dat,
      species = .x,
      family = delta_gamma()
    )
) %>% setNames(list_species)
index_dg <- purrr::map_dfr(out_dg, "index", .id = "species")


# out_pg <- purrr::map(list_species, ~
#     fit_index(
#       dat,
#       formula = catch_weight ~ 1,
#       species = .x,
#       family = delta_poisson_link_gamma()
#     )
# ) %>% setNames(list_species)

# out_dg <- purrr::map(list_species, ~
#     fit_index(
#       dat,
#       species = .x,
#       spatiotemporal = list("off", "rw"),
#       family = delta_gamma()
#     )
# )

index <- index %>% mutate(type = "Tweedie s(depth)")
index_nodepth <- index_nodepth %>% mutate(type = "Tweedie")
index_dg <- index_dg %>% mutate(type = "Delta-Gamma s(depth)")
index_dg_nodepth <- index_dg_nodepth %>% mutate(type = "Delta-Gamma")

mult <- 1000

g <- index %>%
  bind_rows(index_nodepth) %>%
  bind_rows(index_dg) %>%
  bind_rows(index_dg_nodepth) %>%
  group_by(species, type) %>%
  mutate(max_log = max(log_est)) %>%
  filter(max_log < 20) %>%
  mutate(species = stringr::str_to_title(species)) %>%
  ggplot(aes(year, est / mult, colour = type, fill = type,
    ymin = lwr / mult, ymax = upr / mult)) +
  ggsidekick::theme_sleek() +
  facet_wrap(vars(species), scales = "free_y")
g <- g +
  # geom_pointrange(aes(ymin = lwr / mult, ymax = upr / mult))
  geom_line() +
  geom_ribbon(alpha = 0.2, colour = NA)
g <- g + ylab("Biomass (tonnes)") +
  # ggtitle(stringr::str_to_title(species), subtitle = paste(region, collapse = ", ")) +
  coord_cartesian(
    ylim = c(0, 70000),
    expand = FALSE, xlim = range(index$year) + c(-0.25, 0.25)
  ) +
  theme(axis.title.x = element_blank()) +
  facet_wrap(~type)
g

# AIC(out$`arrowtooth flounder`$fit)
# AIC(out_nodepth$`arrowtooth flounder`$fit)
# AIC(out_dg$`arrowtooth flounder`$fit)
# AIC(out_dg_nodepth$`arrowtooth flounder`$fit)

ind <- index %>%
  bind_rows(index_nodepth) %>%
  bind_rows(index_dg) %>%
  bind_rows(index_dg_nodepth)

group_by(ind, type) %>%
  summarise(mean_se = mean(se))

simulate(out$`arrowtooth flounder`$fit, 200) %>%
  dharma_residuals(out$`arrowtooth flounder`$fit)

simulate(out_dg$`arrowtooth flounder`$fit, 200) %>%
  dharma_residuals(out$`arrowtooth flounder`$fit)

simulate(out_dg_nodepth$`arrowtooth flounder`$fit, 200) %>%
  dharma_residuals(out$`arrowtooth flounder`$fit)

r0 <- residuals(out$`arrowtooth flounder`$fit, "mle-mcmc",
  model = 1L, mcmc_iter = 101, mcmc_warmup = 100)

r1 <- residuals(out_dg$`arrowtooth flounder`$fit, "mle-mcmc",
model = 1L, mcmc_iter = 101, mcmc_warmup = 100)

r2 <- residuals(out_dg$`arrowtooth flounder`$fit, "mle-mcmc",
  model = 2L, mcmc_iter = 101, mcmc_warmup = 100)

save(r0, r1, r2, file = "~/src/arrowtooth/arrowtooth-nongit/geostat-figs/mcmc-resids.rda")



m <- out_dg$`arrowtooth flounder`$fit
nd <- readRDS(here("grids/synoptic_grid.rds"))
dat <- out_dg$`arrowtooth flounder`$fit$data
fitted_yrs <- sort(unique(dat$year))
nd <- make_grid(nd, years = fitted_yrs)
nd <- na.omit(nd)
nd$year <- as.integer(nd$year)
nd$log_depth <- log(nd$depth)
p <- predict(m, newdata = nd)

data <- dat
coast <- gfplot:::load_coastline(
  range(data$longitude) + c(-0.5, 0.5),
  range(data$latitude) + c(-0.5, 0.5),
  utm_zone = 9
)
coords <- coord_equal(
  expand = FALSE, xlim = range(data$X) + c(-10, 10),
  ylim = range(data$Y) + c(-10, 10)
)
utm_labs <- labs(x = "Easting", y = "Northing")

plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column, colour = column)) +
    geom_polygon(
      data = coast, aes_string(x = "X", y = "Y", group = "PID"),
      fill = "grey87", col = "grey70", lwd = 0.2, inherit.aes = FALSE
    ) +
    geom_tile(width = 2, height = 2) +
    scale_colour_viridis_c(trans = "sqrt") +
    scale_fill_viridis_c(trans = "sqrt") +
    coord_fixed() +
    coords +
    utm_labs
}

p %>%
  mutate(est_total = plogis(est1) * exp(est2)) %>%
  mutate(est_total = ifelse(est_total > quantile(est_total, probs = 0.995), quantile(est_total, probs = 0.995), est_total)) %>%
  # filter(year %in% 2007:2009) %>%
  plot_map("est_total") +
  facet_wrap(vars(year)) +
  labs(fill = "Biomass density", colour = "Biomass density") +
  gfplot::theme_pbs()

ggsave("")

# g <- ggplot(ind, aes(year, est / 1000)) +
#   ggsidekick::theme_sleek()
# g <- g +
#   geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000),
#     alpha = 0.9, fill = "grey90"
#   ) +
#   geom_line(alpha = 0.9, lty = 1, lwd = 0.5)
# g

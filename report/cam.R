library(dplyr)
library(gfplot)
source(here::here("report/trawl/functions.R"))

filter_matures <- function(spp_hyphen = "pacific-cod", maturity_type = "Mature") {
  f <- paste0("../gfsynopsis-2021/report/data-cache-feb-2023/", spp_hyphen, ".rds")
  dat <- readRDS(f)
  surveys <- "SYN WCVI"
  d_survey_sets <- dat$survey_sets |> filter(survey_abbrev %in% surveys)
  d_survey_samples <- dat$survey_samples |> filter(survey_abbrev %in% surveys)
  d_survey_sets <- d_survey_sets[!duplicated(d_survey_sets$fishing_event_id), ]
  d_by_maturity <- gfplot::split_catch_by_sex(
    d_survey_sets, d_survey_samples,
    plot = TRUE,
    sample_id_re = FALSE,
    immatures_pooled = TRUE
  )
  maturity_assigned <- !is.null(d_by_maturity$maturity)
  if (maturity_assigned) {
    print(d_by_maturity$maturity_plot)
    print(d_by_maturity$weight_plot)
    d <- d_by_maturity$data
    d <- filter(d, grepl(paste0("^", maturity_type), group_name))
    # sum over males and females if needed (needed for matures)
    d <- group_by(
      d, year, fishing_event_id, survey_id,
      survey_abbrev, species_common_name
    ) |>
      summarise(
        catch_weight = sum(group_catch_est), doorspread_m = doorspread_m[1],
        speed_mpm = speed_mpm[1], tow_length_m = tow_length_m[1],
        area_swept = area_swept[1], latitude = latitude[1],
        longitude = longitude[1], depth_m = depth_m[1], duration_min = duration_min[1],
        .groups = "drop"
      )
    d$density_kgpm2 <- d$catch_weight / d$area_swept
  } else {
    d <- d_survey_sets
  }
  d <- filter(d, !is.na(catch_weight))
  list(dat = d, maturity_assigned = maturity_assigned)
}

# matures:
spp_mat <- c(
  "pacific-cod", "north-pacific-spiny-dogfish", "big-skate",
  "longnose-skate", "yellowtail-rockfish", "canary-rockfish",
  "arrowtooth-flounder"
)

purrr::walk(
  spp_mat, function(x) {
    d <- filter_matures(x)
    fit_index(
      region = "SYN WCVI",
      species = gsub("-", " ", x),
      covariate_type = "no depth",
      data = d$dat,
      extra_title = if (d$maturity_assigned) "(mature only)" else ""
    )
  },
  .progress = TRUE
)

# immatures:
spp_imm <- c(
  "lingcod"
)

purrr::walk(
  spp_imm, function(x) {
    d <- filter_matures(x, maturity_type = "Immature")
    fit_index(
      region = "SYN WCVI",
      species = gsub("-", " ", x),
      covariate_type = "no depth",
      data = d$dat,
      extra_title = if (d$maturity_assigned) "(immature only)" else ""
    )
  },
  .progress = TRUE
)

get_design_index <- function(spp_hyphen) {
  f <- paste0("../gfsynopsis-2021/report/data-cache-feb-2023/", spp_hyphen, ".rds")
  d <- readRDS(f)$survey_index
  d <- d |> dplyr::filter(survey_abbrev %in% "HBLL OUT S")
  d
}

spp_hbll <- c("big-skate", "lingcod", "longnose-skate")
dat <- purrr::map_dfr(spp_hbll, get_design_index)
table(dat$species_common_name)

saveRDS(dat, "report/hbll-indexes.rds")

library(gfiphc)
sp <- "pacific halibut"
# cache_pbs_data_iphc(sp, path = "report")
sp_set_counts <- readRDS("report/pacific-halibut.rds")
series_ABCD_full <- calc_iphc_full_res(sp_set_counts$set_counts)
# series_ABCD_full
plot_IPHC_ser_four_panels_ABCD(series_ABCD_full, sp = sp)
series_ABCD_full$type
series_ABCD_full$full_coast
series_ABCD_full$ser_longest
# series_ABCD_full$ser_all$ser_B
# series_ABCD_full$ser_all$ser_C
# series_ABCD_full$ser_all$ser_C
saveRDS(series_ABCD_full$ser_longest, "report/iphc-halibut-20.rds")

# gfiphc::plot_iphc_map(sp_set_counts$set_counts, years = 2002)

sp <- "big skate"
# cache_pbs_data_iphc(sp, path = "report")
sp_set_counts <- readRDS("report/big-skate.rds")
series_ABCD_full <- calc_iphc_full_res(sp_set_counts$set_counts)
# series_ABCD_full
plot_IPHC_ser_four_panels_ABCD(series_ABCD_full, sp = sp)
series_ABCD_full$type
series_ABCD_full$full_coast
series_ABCD_full$ser_longest
saveRDS(series_ABCD_full$ser_longest, "report/iphc-big-skate-20.rds")

sp <- "longnose skate"
# cache_pbs_data_iphc(sp, path = "report")
sp_set_counts <- readRDS("report/longnose-skate.rds")
series_ABCD_full <- calc_iphc_full_res(sp_set_counts$set_counts)
# series_ABCD_full
plot_IPHC_ser_four_panels_ABCD(series_ABCD_full, sp = sp)
series_ABCD_full$type
series_ABCD_full$full_coast
series_ABCD_full$ser_longest
saveRDS(series_ABCD_full$ser_longest, "report/iphc-longnose-skate-20.rds")

i <- list()
i[[1]] <- readRDS("report/iphc-big-skate-20.rds") |>
  select(year, est = I_t20SampleMean, lwr = I_t20BootLow, upr = I_t20BootHigh) |>
  mutate(species = "big skate", survey_abbrev = "IPHC 20", maturity = "all")
i[[2]] <- readRDS("report/iphc-longnose-skate-20.rds") |>
  select(year, est = I_t20SampleMean, lwr = I_t20BootLow, upr = I_t20BootHigh) |>
  mutate(species = "longnose skate", survey_abbrev = "IPHC 20", maturity = "all")
i[[3]] <- readRDS("report/iphc-halibut-20.rds") |>
  select(year, est = I_t20SampleMean, lwr = I_t20BootLow, upr = I_t20BootHigh) |>
  mutate(species = "pacific halibut", survey_abbrev = "IPHC 20", maturity = "all")

spp1 <- c("pacific-cod", "north-pacific-spiny-dogfish", "yellowtail-rockfish", "canary-rockfish", "arrowtooth-flounder"
)
for (j in seq_along(spp1)) {
  i[[length(i) + j]] <- readRDS(paste0("indices/index-", spp1[j], "--no-depth-(mature only).rds")) |> select(year, est, lwr, upr, species, survey_abbrev = survey) |>
    mutate(maturity = "mature only")
}

spp2 <- c("big-skate", "longnose-skate")
for (j in seq_along(spp2)) {
  i[[length(i) + j]] <- readRDS(paste0("indices/index-", spp2[j], "--no-depth-.rds")) |> select(year, est, lwr, upr, species, survey_abbrev = survey) |>
    mutate(maturity = "all")
}

spp3 <- c("lingcod")
for (j in seq_along(spp3)) {
  i[[length(i) + j]] <- readRDS(paste0("indices/index-", spp3[j], "--no-depth-(immature only).rds")) |> select(year, est, lwr, upr, species, survey_abbrev = survey) |>
    mutate(maturity = "immature only")
}

temp <- readRDS("report/hbll-indexes.rds") |>
  select(year, est = biomass, lwr = lowerci, upr = upperci, survey_abbrev, species = species_common_name) |>
  mutate(maturity = "all")

i[[length(i) + 1]] <- temp

ind <- bind_rows(i)

library(ggplot2)
ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) + geom_line() +
  geom_ribbon(alpha = 0.5) +
  facet_wrap(~paste(species, "\n", survey_abbrev, "\n", maturity), scales = "free_y") +
  theme_light()
ggsave("indices-cam.pdf", width = 10, height = 7)
saveRDS(ind, file = "indices-cam.rds")




# f <- c(
#   "index-arrowtooth-flounder--no-depth-(mature only).rds",
#   "index-big-skate--no-depth-.rds",
#   "index-canary-rockfish--no-depth-(mature only).rds",
#   "index-lingcod--no-depth-(immature only).rds",
#   "index-longnose-skate--no-depth-.rds",
#   "index-north-pacific-spiny-dogfish--no-depth-(mature only).rds",
#   "index-pacific-cod--no-depth-(mature only).rds",
#   "index-yellowtail-rockfish--no-depth-(mature only).rds")




# fake hbll as synoptic trawl:
#
# fake_hbll_as_trawl <- function(spp_hyphen) {
#   f <- paste0("../gfsynopsis-2021/report/data-cache-feb-2023/", spp_hyphen, ".rds")
#   dat <- readRDS(f)
#   surveys <- "HBLL OUT S"
#   d_survey_sets <- dat$survey_sets |> filter(survey_abbrev %in% surveys)
#   d <- d |> select(
#     year, fishing_event_id, survey_id,
#     survey_abbrev, species_common_name,
#     catch_count,
#     latitude,
#     longitude,
#     depth_m, hook_count
#   ) |>
#     rename(catch_weight = catch_count, area_swept = hook_count) |>
#     mutate(density_kgpm2 = catch_weight / area_swept)
#   d
# }
#
# spp <- c(
#   "big-skate",
#   "longnose-skate",
#   "lingcod"
# )
#
# purrr::walk(
#   spp, function(x) {
#     d <- fake_hbll_as_trawl(x)
#     browser()
#     fit_index(
#       region = "HBLL OUT S",
#       species = gsub("-", " ", x),
#       covariate_type = "no depth",
#       data = d
#     )
#   },
#   .progress = TRUE
# )
#
# fake_hbll_as_trawl(spp[1])
#
# # halibut - iphc - big and longnose skate

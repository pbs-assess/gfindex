# from SOPO repo:
data <- readRDS(here::here("data/all-survey-sets-2021.rds")) %>%
  rename(survey_series_id = survey_series_id.x) %>%
  filter(depth_m > 0) %>%
  filter(survey_abbrev %in% c("SYN WCHG", "SYN QCS", "SYN HS", "SYN WCVI"))

data <- select(
  data,
  survey_series_id,
  year,
  # fishing_event_id,
  latitude,
  longitude,
  grouping_code,
  depth_m,
  duration_min,
  doorspread_m,
  speed_mpm,
  tow_length_m,
  catch_weight,
  density_kgpm2,
  # catch_count,
  survey_abbrev,
  species_common_name,
  species_science_name
  # sample_id
  # area_km2
)
saveRDS(data, here::here("data/all-survey-sets-2021-select.rds"))

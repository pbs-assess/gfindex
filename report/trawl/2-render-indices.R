PARALLEL <- FALSE
CORES <- floor(parallel::detectCores() / 2 - 1)

library(dplyr)
source(here::here("report/trawl/functions.R"))
if (PARALLEL) setup_parallel(CORES)

# 1. adjust species in report/trawl/spp.txt ---------------------------

list_species <- read_species_list(here::here("report/trawl/spp.txt"))

# 2. select regions: --------------------------------------------------

list_regions <- c(
  # "SYN QCS"
  "SYN WCVI"
  # "SYN HS",
  # "SYN WCHG"
)

# 3. select covariate type: -------------------------------------------

# list_covariate_type <- c("with depth", "no depth")
list_covariate_type <- c("with depth")

# 4. grab the data if needed: -----------------------------------------

f <- here::here("data/all-survey-sets-2021-select.rds")
if (!file.exists(f)) {
  dropbox_file <- paste0("https://www.dropbox.com/s/uduck48z9jnkq9f/",
    "all-survey-sets-2021-select.rds?dl=1")
  download.file(dropbox_file, destfile = f)
}

# 5. run the following: -----------------------------------------------

to_fit <- make_dat(list_regions, list_species, list_covariate_type)
if (PARALLEL) {
  furrr::future_pwalk(to_fit, fit_index)
} else {
  purrr::pwalk(to_fit, fit_index)
}

# avoid crashing RStudio:
if (PARALLEL) future::plan(future::sequential)

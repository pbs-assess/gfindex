# update or call a previous version of sdmTMB
# remotes::install_github("pbs-assess/sdmTMB", ref = "delta")
library(dplyr)

make_dat <- function(r, s) {
  expand.grid(
    species = s,
    region = r,
    stringsAsFactors = FALSE
  )
}

list_regions <- c(
      "SYN QCS"
    , "SYN WCVI"
    , "SYN HS"
    , "SYN WCHG"
    )

list_species <- c(
  "Canary Rockfish"
)

to_fit <- make_dat(list_regions, list_species)

# https://github.com/rstudio/rmarkdown/issues/1673
render_separately <- function(...) callr::r(
  function(...) rmarkdown::render(..., envir = globalenv()),
  args = list(...), show = TRUE)

fit_index <- function(region, species) {
  spp <- gsub(" ", "-", gsub("\\/", "-", tolower(species)))
  # name <- "with depth" # describe model covariates
  name <- "no depth" # describe model covariates
  region_name <- region
  try({
    render_separately("report/trawl/index-standardization.Rmd",
      params = list(
        species = species,
        region = region,
        name = name,
        update_model = TRUE,
        update_index = TRUE,
        silent = TRUE
      ),
      output_file = paste0(spp, "-",
        gsub(" ", "-", gsub("\\/", "-", name)), "-",
        gsub(" ", "-", gsub("\\/", "-", region_name)), ".html")
    )
  })
}

# purrr::pwalk(to_fit, fit_index)

is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (is_unix && !is_rstudio) {
  future::plan(future::multicore, workers = 6L)
} else {
  future::plan(future::multisession, workers = 2L)
}
options(future.rng.onMisuse = "ignore")

furrr::future_pwalk(to_fit, fit_index)

future::plan(future::sequential)


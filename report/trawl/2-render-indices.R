library(dplyr)

make_dat <- function(r, s, .c) {
  expand.grid(
    species = s,
    region = r,
    covariate_type = .c,
    stringsAsFactors = FALSE
  )
}

list_regions <- c(
  "SYN QCS",
  "SYN WCVI",
  "SYN HS",
  "SYN WCHG"
)

list_species <- c(
  "Canary Rockfish"
)

list_covariate_type <- c("with depth", "no depth")

to_fit <- make_dat(list_regions, list_species, list_covariate_type)

# https://github.com/rstudio/rmarkdown/issues/1673
render_separately <- function(...) {
  callr::r(
    function(...) rmarkdown::render(..., envir = globalenv()),
    args = list(...), show = TRUE
  )
}

fit_index <- function(region, species, covariate_type) {
  spp <- gsub(" ", "-", gsub("\\/", "-", tolower(species)))
  region_name <- region
  try({
    render_separately("report/trawl/index-standardization.Rmd",
      params = list(
        species = species,
        region = region,
        covariate_type = covariate_type,
        update_model = TRUE,
        update_index = TRUE,
        silent = TRUE
      ),
      output_file = paste0(
        spp, "-",
        gsub(" ", "-", gsub("\\/", "-", covariate_type)), "-",
        gsub(" ", "-", gsub("\\/", "-", region_name)), ".html"
      )
    )
  })
}

# purrr::pwalk(to_fit, fit_index)

is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (is_unix && !is_rstudio) {
  future::plan(future::multicore, workers = 4L)
} else {
  future::plan(future::multisession, workers = 4L)
}
options(future.rng.onMisuse = "ignore")

furrr::future_pwalk(to_fit, fit_index)

future::plan(future::sequential)

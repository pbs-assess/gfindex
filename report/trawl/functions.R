read_species_list <- function(f) {
  x <- readr::read_delim(f, comment = "#", delim = "\n", col_types = "c")
  x$species
}

make_dat <- function(r, s, .c) {
  expand.grid(
    species = s,
    region = r,
    covariate_type = .c,
    stringsAsFactors = FALSE
  )
}

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

setup_parallel <- function(cores) {
  is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
  is_unix <- .Platform$OS.type == "unix"
  if (is_unix && !is_rstudio) {
    future::plan(future::multicore, workers = cores)
  } else {
    future::plan(future::multisession, workers = cores)
  }
  options(future.rng.onMisuse = "ignore")
}

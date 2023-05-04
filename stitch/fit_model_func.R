# Test model on survey data
check_sanity <- function(x) {
  if (!all(unlist(sanity(x)))) {
    return(NA)
  } else {
    return(x)
  }
}

fit_models <- function(
    dat, catch, data_subset = NULL, mesh = NULL, cutoff = 20, family = tweedie(), 
    offset = NULL, 
    ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)) {

  if (is.null(data_subset)) {
    data_subset <- unique(dat$species_common_name)
  }

  message(cat("\tFitting models for data subset:", data_subset, "\n"))

  if (is.null(mesh)) {
    message(cat("\n\tNo mesh provided, making mesh with cutoff:", cutoff, "\n"))
    mesh <- make_mesh(dat, c("X", "Y"), cutoff = 20)
  }

  dat <- droplevels(dat)  # drop extra factor levels before running models
  fits <- list()
  i <- 1

  cli::cli_inform("\tFitting st = 'rw'")
  fit1 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1,
      family = family,
      data = dat, time = "year", spatiotemporal = "rw", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit1
  model_ids <- i
  i <- i + 1

  cli::cli_inform("\tFitting st IID covariate")
  fit2 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit2
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st IID s(year)")
  fit3 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ s(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit3
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st IID no covariate as.factor year")
  fit4 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit4
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st time_varying RW")
  fit5 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0,
      family = family,
      time_varying = ~1, time_varying_type = "rw",
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit5
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st (1|year)")
  fit6 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1 + (1 | fyear),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit6
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting spatial only")
  fit7 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "off", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit7
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st (1 | region)")
  fit8 <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + fyear + region,
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit8
  model_ids <- c(model_ids, i)

  names(fits) <- paste(model_ids, data_subset, sep = ":")
  return(fits)
}

# is_even <- function(column) {
#   ifelse({{column}} %% 2 == 0, TRUE, FALSE)
# }

# -------
get_pred_list <- function(fit_list, newdata) 
  fit_list %>% 
  purrr::map(., function(.x) {
    if (inherits(.x, "sdmTMB")) {
      out <- predict(.x, newdata = newdata, return_tmb_object = TRUE)
    } else {
      out <- NA
    }
    out
  })

get_index_list <- function(pred_list, area = 1) {
  purrr::map(pred_list, function(.x) {
    if (length(.x) > 1) {
      get_index(.x, bias_correct = TRUE, area = area)
    } else {
      out <- NA  # keep empty fits as visual cue that these did not fit when plotting
    }
  })
}

mk_index_df <- function(index_list) {
  enframe(index_list) %>% 
  unnest(col = "value") %>% 
  separate(col = 'name', into = c('id', 'species'), sep = ":") %>% 
  mutate(id = as.numeric(id)) %>% 
  right_join(., model_lookup) %>% 
  mutate(odd_even = ifelse(year %% 2 == 0, "even", "odd"))
}
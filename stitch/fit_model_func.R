# Test model on survey data
fit_models <- function(dat, data_subset = NULL, mesh, family = tweedie(), ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)) {
  if (is.null(data_subset)) {
    data_subset <- unique(dat$species_common_name)
  }
  fits <- list()
  i <- 1
  
  ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)

  cli::cli_inform("\tFitting st = 'rw'")
  fit1 <- try(
    sdmTMB(
      catch_weight ~ 1,
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "rw", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit1
  model_ids <- i
  i <- i + 1

  cli::cli_inform("\tFitting st IID covariate")
  fit2 <- try(
    sdmTMB(
      catch_weight ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit2
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st IID s(year)")
  fit3 <- try(
    sdmTMB(
      catch_weight ~ s(year),
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit3
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st IID no covariate as.factor year")
  fit4 <- try(
    sdmTMB(
      catch_weight ~ 0 + as.factor(year),
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit4
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st time_varying RW")
  fit5 <- try(
    sdmTMB(
      catch_weight ~ 0,
      family = tweedie(),
      time_varying = ~1, time_varying_type = "rw",
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit5
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st (1|year)")
  fit6 <- try(
    sdmTMB(
      catch_weight ~ 1 + (1 | fyear),
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit6
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting spatial only")
  fit7 <- try(
    sdmTMB(
      catch_weight ~ 0 + as.factor(year),
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "off", spatial = "on",
      mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit7
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting st (1 | region)")
  fit8 <- try(
    sdmTMB(
      catch_weight ~ 0 + fyear + region,
      family = tweedie(),
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = dat$offset,
      control = ctrl
    )
  )
  fits[[i]] <- fit8
  model_ids <- c(model_ids, i)

  names(fits) <- paste(model_ids, data_subset, sep = ":")
  fits <- enframe(fits)
  return(fits)
}
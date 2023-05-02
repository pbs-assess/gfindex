# Prepare grids
syn_grid <- 
  gfplot::synoptic_grid %>%
  tibble() %>%  # because I accidentally print the full df too often
  dplyr::select(-survey_series_name, -utm_zone, -survey_domain_year)

# 
hbll_n_grid  <- 
  gfplot::hbll_n_grid$grid %>% tibble() %>% 
  rename(longitude = "X", latitude = "Y") %>% 
  add_utm_columns(c("longitude", "latitude"), utm_crs = 32609)
hbll_s_grid  <- 
  gfplot::hbll_s_grid$grid %>% tibble() %>% 
  rename(longitude = "X", latitude = "Y") %>% 
  add_utm_columns(c("longitude", "latitude"), utm_crs = 32609)

# Inside grids (between VI and mainland)
hbll_inside_n_grid <- gfplot::hbll_inside_n_grid$grid  %>% tibble() %>% 
  rename(longitude = "X", latitude = "Y") %>% 
  add_utm_columns(c("longitude", "latitude"), utm_crs = 32609)

hbll_inside_s_grid <- gfplot::hbll_inside_s_grid$grid  %>% tibble() %>% 
  rename(longitude = "X", latitude = "Y") %>% 
  add_utm_columns(c("longitude", "latitude"), utm_crs = 32609)

hbll_inside <- bind_rows(hbll_inside_n_grid, hbll_inside_s_grid)

survey_region_baseplot <- 
ggplot(mapping = aes(x = X, y = Y)) + 
  geom_tile(data = syn_grid, aes(width = 1.5, height = 1.5), colour = "grey90", alpha = 0.2) +
  geom_tile(data = hbll_n_grid, aes(width = 0.9, height = 0.9), colour = "black", alpha = 0.2) +
  geom_tile(data = hbll_s_grid, aes(width = 0.9, height = 0.9), colour = "darkslateblue", alpha = 0.2) +
  geom_tile(data = hbll_inside, aes(width = 0.9, height = 0.9), colour = "darkslateblue", alpha = 0.2)


formulae <- list(
  
  baseline = bf(
    n_npev | trials(n_npafp) ~ 
      # IID random intercepts on district
      (1 | guid) + 
      # Periodic spline on 1-12 for seasonality
      s(month_of_year, bs = "cc", k = 12) 
    ),
  
  m1 = bf(
    n_npev | trials(n_npafp) ~ 
      # Fixed effect on population density
      log(guid_pop_dens_c) +
      # 1st order autocorrelation by month
      ar(t, p = 1, gr = guid) +
      # BYM2 random effect across districts
      car(W, gr = guid, type = "bym2") +   
      # Periodic spline on 1-12 for seasonality
      s(month_of_year, bs = "cc", k = 12) 
  ),
  
  m2 = bf(
    n_npev | trials(n_npafp) ~ 
      # Fixed effect on population density
      log(guid_pop_dens_c) + 
      # 1st order random walk by month
      (1|t) +
      # BYM2 random effect across districts
      car(W, gr = guid, type = "bym2") +   
      # Global seasonal spline
      s(month_of_year, bs = "cc", k = 12) +
      # Seasonal spline by year
      s(month_of_year, bs = "cc", k = 12, by = year) 
  ),
  
  m3 = bf(
    n_npev | trials(n_npafp) ~ 
      # Fixed effect on population density 
      log(guid_pop_dens_c) +
      # 1st order autocorrelation by month
      ar(t, p = 1, gr = guid) +
      # BYM2 random effect across districts
      car(W, gr = guid, type = "bym2") +   
      # Global seasonal spline
      s(month_of_year, bs = "cc", k = 12) +
      # Seasonal spline by province
      s(month_of_year, bs = "cc", k = 12, by = adm1_name) 

   kroiss = bf(
    n_npev | trials(n_npafp) ~
      # Random intercepts by province
      (1 | adm1_name) +
      # 1st order autocorrelation by month
      ar(t, p = 1, gr = guid) +
      # Spatially-structured random intercepts by district
      car(W, gr = guid, type = "bym2") +
      # Global cyclic spline
      s(month_of_year, bs = "cc", k = 12) +
      # Cyclic spline by province
      s(month_of_year, bs = "cc", k = 12, by = adm1_name) 
    )
)

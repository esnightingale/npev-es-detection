check_model_fit <- function(fit){
  
  check_plots <- list(pp_check(fit, type = "bars"), 
                      pp_check(fit, type = "scatter_avg"),
                      pp_check(fit, "rootogram"))  
  return(check_plots)
}


check_residuals <- function(fit, data, shapefile){
  
  res <- residuals(fit) |> as.data.frame()
  data <- data |> 
    mutate(residual = res$Estimate, 
           residual_err = res$Est.Error,
           residual_std = residual/residual_err)
  
  # Check if residuals still show temporal autocorrelation
  ggplot(data, aes(x = t, y = residual_std, group = guid)) +
    geom_point(alpha = 0.2) +
    geom_smooth(se = F, lwd = 0.1) +
    labs(x = "Time", y = "Std. residuals (estimate/error)", title = "Temporal Autocorrelation in Residuals") -> p_time
  
  # Check if residuals still show spatial autocorrelation
  data |> 
    group_by(guid) |> 
    summarise(resid_mean = mean(residual),
              resid_sd = sd(residual),
              resid_std = resid_mean/resid_sd,
              resid_mean_err = mean(residual_err),
              resid_mean_std = mean(residual_std)) |> 
    right_join(shapefile) -> tmp
  
  tmp |> 
    ggplot(aes(geometry = SHAPE, fill = resid_std)) +
    geom_sf() +
    scale_fill_gradient2() +
    theme_void() + 
    labs(subtitle = "Standardised average residual (mean(residual)/sd(residual))") -> p_space1
  
  tmp |> 
    ggplot(aes(geometry = SHAPE, fill = resid_mean_err)) +
    geom_sf() +
    scale_fill_gradient2() +
    theme_void() + 
    labs(subtitle = "Average error of residuals (mean(error))") -> p_space2
  
  tmp |> 
    ggplot(aes(geometry = SHAPE, fill = resid_mean_std)) +
    geom_sf() +
    scale_fill_gradient2() +
    theme_void() + 
    labs(subtitle = "Average standardised residual (mean(residual/error))") -> p_space3
  
  print(list(p_time, p_space1, p_space2, p_space3))
  
}


 
# plot_reff_time <- function(fit){
#   
#   plot(conditional_effects(fit, effects = "month_of_year"), points = TRUE)
#   
#   # Extract smooth month estimates
#   ce <- conditional_effects(fit, effects = "month_of_year")
#   month_effects <- ce$month_of_year  # Extract month effect data
#   
#   # Plot with ggplot2
#   ggplot(month_effects, aes(x = month_of_year, y = estimate__)) +
#     geom_line(size = 1, color = "blue") + 
#     geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3) + 
#     labs(x = "Month", y = "Estimated Effect", title = "Cyclical Effect of Month") +
#     theme_minimal()
#   
#   posterior_samples(fit, pars = "sds_smonth_1") %>% 
#     ggplot(aes(x = .value)) +
#     geom_histogram(fill = "skyblue", bins = 30) +
#     labs(x = "Effect Size", y = "Density", title = "Posterior Distribution of Monthly Effect")
#   
# }
# 
# plot_reff_space <- function(fit, shapefile){
#   
#   spatial_effects <- ranef(fit)$region  # Extract spatial random effects
#   
#   # Decompose BYM2 spatial effect
#   structured_effect <- spatial_effects[, , "gr(guid, by = 'bym2')"] # grepl("bym2", rownames(spatial_effects))
#   iid_effect <- spatial_effects[, , "(Intercept)"]
#   total_spatial_effect <- structured_effect + iid_effect
#   
#   # Convert spatial effects to a data frame
#   spatial_df <- data.frame(
#     region = rownames(total_spatial_effect),
#     effect = total_spatial_effect[, "Estimate"]
#   )
#   
#   # Merge with spatial data
#   shapefile <- shapefile %>%
#     left_join(spatial_df, by = "region")
#   
#   ggplot(shapefile) +
#     geom_sf(aes(fill = effect), color = "black") +  # Plot spatial effect
#     scale_fill_viridis_c(option = "magma", name = "Spatial Effect") +  # Color scale
#     theme_void() +
#     labs(title = "Estimated Spatial Random Effect (BYM2)")
# }

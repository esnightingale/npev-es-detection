fit_model <- function(data, W, formula, prior){

  # Fit brms model
  fit <- brms::brm(formula = formula,
                       data = data, 
                       data2 = list(W = W),
                       prior = prior,
                       family = binomial(link = "logit"),
                       iter = 2000, chains = 4, cores = 4,
                       control = list(adapt_delta = 0.95),  # To avoid divergences
                       save_pars = save_pars(all = T)
                       )
  
  return(fit)
}

plot_pp_checks <- function(fit){
  
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

plot_intercepts_adm1 <- function(fit){
  
  fit %>%
    spread_draws(`b_Intercept`, r_adm1_name[adm1_name,]) %>%
    mutate(adm1_mean = b_Intercept + r_adm1_name) %>%
    median_qi(adm1_mean, .width = c(.95, .8, .5)) -> intercepts_rand
  
  # Plot estimates
  ggplot(intercepts_rand, 
         aes(y = factor(adm1_name, 
                        levels = rev(levels(factor(adm1_name)))), 
             x = adm1_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Intercept", y = "Province") |> return()

}

plot_intercepts_adm2 <- function(fit){
  
  fit$data |> 
    mutate(guid = as.numeric(as.factor(guid))) |>  
    select(guid, adm1_name) |> 
    distinct() -> adm1_adm2
  
  fit |> 
    spread_draws(`b_Intercept`, r_adm1_name[adm1_name,], rcar[guid]) |> 
    left_join(adm1_adm2, by = "guid") |> 
    mutate(flag = ifelse(adm1_name.x == adm1_name.y, 1, 0), 
           adm2_mean = b_Intercept + r_adm1_name*flag + rcar) |> 
    group_by(adm1_name.y, guid) |> 
    median_qi(adm2_mean) |> 
    mutate(did = row_number()) -> tmp

  # Plot estimates
  ggplot(tmp, 
         aes(y = did, x = adm2_mean, col = adm1_name.y,
             xmin = .lower, xmax = .upper)) +
    geom_pointinterval(cex = 0.5, lwd = 0.1) +
    scale_y_discrete(drop = TRUE, expand = c(0, 0)) +
    facet_grid(adm1_name.y~., scales = "free", space = "free_y", switch = "y") +
    theme(strip.placement = "outside",
          panel.spacing = unit(0, "in"),
          strip.background.y = element_rect(fill = "white", color = "gray75")) +
    guides(col = "none") +
    labs(x = "Intercept", y = NULL) |> 
    return()
  
}

plot_season_prov <- function(fit){
  
  fit |> 
    conditional_smooths() |>
    plot(plot = F) -> plots
  
  ggplot_build(plots[[1]]) -> p1
  ggplot_build(plots[[2]]) -> p2

p2$plot$data |> 
  rename(Province = adm1_name) |> 
  ggplot(mapping = p2$plot$mapping) + 
  geom_ribbon(aes(fill = Province), alpha = 0.1) +
  geom_line(aes(col = Province)) +
  geom_line(data = p1$plot$data, lwd = 0.2) + 
  geom_ribbon(data = p1$plot$data, 
              fill = NA, col = "black", lty = 2, lwd = 0.1) + 
    labs(x = NULL, y = "Provincial seasonal effect") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = month.abb) |> return()
  
}

plot_season_year <- function(fit){
  
  fit |> 
    conditional_smooths() |>
    plot() -> plots
  
  ggplot_build(plots[[1]]) -> p1
  ggplot_build(plots[[2]]) -> p2
  
  p2$plot$data |> 
    rename(Year = year) |> 
    ggplot(mapping = p2$plot$mapping) + 
    geom_ribbon(aes(fill = Year), alpha = 0.1) +
    geom_line(aes(col = Year)) +
    geom_line(data = p1$plot$data, lwd = 0.2) + 
    geom_ribbon(data = p1$plot$data, 
                fill = NA, col = "black", lty = 2, lwd = 0.1) + 
    labs(x = NULL, y = "Seasonal effect by year") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = month.abb) |> return()
  
}

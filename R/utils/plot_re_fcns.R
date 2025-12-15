plot_re_adm1 <- function(fit){
  
fit %>%
  gather_draws(r_adm1_name[adm1_name, Intercept]) %>%
  median_qi() -> tmp

tmp |> 
  ggplot(aes(adm1_name, .value,
             ymin = .lower, ymax = .upper,
             col = adm1_name)) + 
  geom_hline(yintercept = 0, col = "grey") +
  geom_errorbar() + 
  geom_point() +
  labs(x = NULL, y = "Posterior estimate", col = "Province",
       subtitle = "Province-level random effect") +
  theme(axis.text.x = element_blank()) + 
  scale_colour_viridis_d(option = "turbo") -> p

  return(p)

}
plot_re_guid <- function(fit){
  
  fit %>%
    gather_draws(r_guid[guid, Intercept]) %>%
    median_qi() |> 
    left_join(select(fitdata, guid, adm1_name) |> distinct()) -> tmp
  
  tmp |> 
    arrange(adm1_name, guid) |> 
    mutate(did = row_number()) |> 
    ggplot(aes(did, .value,
               ymin = .lower, ymax = .upper,
               col = adm1_name)) + 
    geom_hline(yintercept = 0, col = "grey") +
    geom_errorbar() + 
    geom_point() +
    labs(x = "District", y = "Posterior estimate", col = "Province",
         subtitle = "District-level random effect") +
    scale_colour_viridis_d(option = "turbo") +
    theme(axis.text.x = element_blank()) -> p
  
  return(p)
  
}

plot_re_adm1_guid <- function(fit, prov_guid, prov_lines = T){

  fit %>%
    spread_draws(r_adm1_name[adm1_name, Intercept],r_guid[guid, Intercept]) %>% 
    left_join(prov_guid, by = "guid") |> 
    filter(adm1_name.x == adm1_name.y) |> 
    median_qi(r_guid = r_adm1_name + r_guid) -> tmp
  
  tmp |> 
    left_join(prov_guid) |> 
    arrange(adm1_name, guid) |> 
    mutate(did = row_number()) |> 
    ggplot(aes(did, r_guid,
               ymin = .lower, ymax = .upper,
               col = adm1_name)) + 
    geom_hline(yintercept = 0, col = "grey") +
    geom_errorbar() + 
    geom_point() +
    labs(x = "District", y = "Posterior estimate", col = "Province",
         subtitle = "Total district deviation (province + district random effects)") +
    scale_colour_viridis_d(option = "turbo") +
    theme(axis.text.x = element_blank()) -> p
  
  if(prov_lines){
    
    fit %>%
      gather_draws(r_adm1_name[adm1_name, Intercept]) %>%
      median_qi() -> tmp2
    
    p <- p + 
      geom_hline(data = tmp2,
                 aes(yintercept = .value, col = adm1_name),
                 lwd = 0.2)
  }
  
  return(p)
}

plot_re_t <- function(fit){
  
  fit %>%
    gather_draws(r_t[t, Intercept]) %>%
    median_qi() |> 
    mutate(time = 1:n()) |> 
    left_join(select(fitdata, time, month) |> distinct()) -> tmp
  
  tmp |>  
    ggplot(aes(month, .value,
               ymin = .lower, ymax = .upper)) + 
    geom_hline(yintercept = 0, col = "grey") +
    geom_errorbar() + 
    geom_point() +
    labs(x = "Month", y = "Posterior estimate",
         subtitle = "Temporal random effect") -> p
  
    return(p)
  
}

plot_re_t <- function(fit){
  
  fit %>%
    gather_draws(r_t[t, Intercept]) %>%
    median_qi() |> 
    left_join(select(fitdata, t, month) |> distinct()) -> tmp
  
  tmp |>  
    ggplot(aes(month, .value,
               ymin = .lower, ymax = .upper)) + 
    geom_hline(yintercept = 0, col = "grey") +
    geom_errorbar() + 
    geom_point() +
    labs(x = "Month", y = "Posterior estimate") |> 
    return()
  
}

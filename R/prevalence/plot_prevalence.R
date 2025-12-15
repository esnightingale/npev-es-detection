################################################################################
################################################################################
# Plot predicted prevalence per district/month
################################################################################

pacman::p_load(tidyverse, here, gtsummary, sf, brms, tidybayes, bayesplot,  
               units, loo, posterior, patchwork)
theme_set(theme_minimal())

dir <- "data/Pakistan/analysis/prevalence"
outdir <- "output/prevalence/final"
figdir <- "figures/prevalence/prediction"

shape2 <- read_rds(here(dir,"../shape2.rds"))

pred <- read_rds(here(outdir,"pred_prev_npev.rds"))
pred_draws <- read_rds(here(outdir,"draws_prev_npev.rds"))

# Set colour palette for provinces
pal_prov <- viridis::turbo(n_distinct(shape2$adm1_name))  
names(pal_prov) <- unique(shape2$adm1_name)

# Observed vs predicted ---------------------------------------------------

pred |> 
  arrange(n_npafp) |> 
  ggplot(aes(n_npev, col = n_npafp, size = n_npafp,
             y = pred_y, ymin = low_y, ymax = hi_y)) + 
  geom_abline() +
  geom_jitter(alpha = 0.7) +  
  guides(size = "none") +
  scale_colour_viridis_b(option = "mako", direction= -1) +
  labs(x = "Observed NPEV positives (per district/month)",
       y = "Posterior prediction",
       col = "No. tested") + 
  theme(legend.position = c(0.85,0.2), legend.background = element_rect(fill = "white"))
ggsave(here(figdir,"pred_y_scatter.png"),
       height = 6, width = 8, bg = "white")

pred |> 
  filter(n_npafp > 5) |> 
  arrange(n_npafp) |> 
  ggplot(aes(npev_npafp_p, size = n_npafp, col = n_npafp, 
             y = pred_p, ymin = low_p, ymax = hi_p)) + 
  geom_abline() +
  geom_point(alpha = 0.7) + 
  guides(size = "none") +
  scale_colour_viridis_b(option = "mako", direction= -1) +
  # scale_alpha_continuous(range = c(0.01,1)) + #c(0.01,0.1,0.5,0.75,1)
  labs(x = "Observed NPEV prevalence (per district/month)",
       y = "Posterior prediction",
       col = "No. tested") + 
  theme(legend.position = c(0.85,0.2), legend.background = element_rect(fill = "white"))
ggsave(here(figdir,"pred_prev_scatter.png"),
       height = 6, width = 8, bg = "white")

pred |> 
  filter(n_npafp >=10) |> 
  arrange(n_npafp) |> 
  ggplot(aes(npev_npafp_p, size = n_npafp, col = n_npafp, 
             y = pred_p, ymin = low_p, ymax = hi_p)) + 
  geom_abline() +
  geom_jitter(alpha = 0.7) + 
  guides(size = "none") +
  scale_colour_viridis_b(option = "mako", direction= -1) +
  labs(x = "Observed NPEV prevalence (per district/month)",
       y = "Posterior prediction",
       col = "No. tested",
       caption = "Data subset to district:months with at least 10 NP-AFP notifications") + 
  theme(legend.position = c(0.85,0.2), legend.background = element_rect(fill = "white"))
ggsave(here(figdir,"pred_prev_scatter_gt10.png"),
       height = 6, width = 8, bg = "white")

# Many month:district with a single NPAFP notification. 
# - These still add some information due to pooling


# By month ----------------------------------------------------------------

pred |> 
  ggplot(aes(month, group = guid)) + 
  geom_jitter(aes(y = npev_npafp_p), 
              col = "grey", alpha = 0.5, cex = 0.2) +
  geom_ribbon(aes(ymin = low_p, ymax = hi_p), 
              alpha = 0.2, fill = "lightsteelblue") + 
  geom_line(aes(y = pred_p), alpha = 0.2, col = "steelblue4") + 
  # geom_line(aes(y = pred_mean), col = "red", alpha = 0.1) +
  labs(y = "Estimated prevalence", x = "Month", title= "Predicted mean & 95% quantile interval")

pred |> 
  ggplot(aes(month, group = guid)) + 
  geom_jitter(aes(y = npev_npafp_p), 
              col = "grey", alpha = 0.5, cex = 0.2) +
  geom_ribbon(aes(ymin = pred_p - 1.96*sd_p, ymax = pred_p + 1.96*sd_p), 
              alpha = 0.2, fill = "lightsteelblue") + 
  geom_line(aes(y = pred_p), alpha = 0.2, col = "steelblue4") + 
  # geom_line(aes(y = pred_mean), col = "red", alpha = 0.1) +
  labs(y = "Estimated prevalence", x = "Month", title= "Predicted mean & 95% gaussian interval")


# By district -------------------------------------------------------------

pred |> 
  ggplot(aes(month, pred_p, col = guid, fill = guid)) +
  geom_ribbon(aes(ymin = low_p, ymax = hi_p), alpha = 0.2, col = NA) +
  geom_line(lwd = 0.2) + 
  geom_line(aes(y = npev_npafp_p), lty = "dashed") +
  scale_colour_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  facet_wrap(~guid) +
  guides(col = "none", fill = "none") +
  labs(x = NULL, y = "Estimated prevalence") -> p_pred_guid

p_pred_guid
pdf(here(figdir,"pred_prev_guid.pdf"), 
    height = 20, width = 30)
p_pred_guid
dev.off()



# Aggregated --------------------------------------------------------------

## Monthly total
summdata <- pred |> 
  group_by(month) |> 
  summarise(p_obs = sum(n_npev)/sum(n_npafp)) |> ungroup()

pred_draws |> 
  group_by(month) |>
  summarise(pred = mean(.linpred),
            lo = quantile(.linpred, 0.025),
            hi = quantile(.linpred, 0.975)) |> 
  ungroup() |> 
  ggplot(aes(month, pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.4, fill = "grey60") + 
  geom_line(col = "grey60") + 
  # stat_lineribbon(aes(y = .linpred), .width = c(.95, .80, .50), alpha = 0.2, col = NA) +
  # geom_point(data = summdata, aes(y = p_obs)) +
  geom_line(data = summdata, aes(y = p_obs), lty = "dashed") +
  scale_fill_brewer() +
  labs(x = NULL, y = "Estimated prevalence") -> p_pred_total

p_pred_total
ggsave(here(figdir,"pred_mth_tot.png"),
       p_pred_total,
       height = 6, width = 10, bg = "white")


## Monthly by province 

summdata_prov <- pred |> 
  group_by(month, adm1_name) |> 
  summarise(p_obs = sum(n_npev)/sum(n_npafp)) |> ungroup()

pred_draws |> 
  group_by(month, adm1_name) |>
  summarise(pred = mean(.linpred),
            lo = quantile(.linpred, 0.025),
            hi = quantile(.linpred, 0.975)) |> 
  ungroup() |> 
  ggplot(aes(month, pred, fill = adm1_name, col = adm1_name)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, col = NA), alpha = 0.4) + 
  geom_line(aes(col = adm1_name)) + 
  # stat_lineribbon(aes(y = .linpred), .width = c(.95, .80, .50), alpha = 0.2, col = NA) +
  geom_line(data = summdata_prov, aes(y = p_obs, col = adm1_name), lty = "dashed") +
  # geom_point(data = summdata_prov, aes(y = p_obs), col = "black") +
  facet_wrap(~adm1_name) +
  scale_colour_viridis_d(option = "turbo", na.translate = F) +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill = "none") + 
  labs(x = NULL, y = "Estimated prevalence", col = "Province") -> p_pred_prov

p_pred_prov
ggsave(here(figdir,"pred_mth_prov.png"),
       p_pred_prov,
       height = 6, width = 10, bg = "white")


# Manuscript figure -------------------------------------------------------

p_pred_total + p_pred_prov + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths = c(4,5))

ggsave(here(figdir,"pred_mth_tot_prov.png"),
       height = 6, width = 14, bg = "white")


# Spatial trends ----------------------------------------------------------
# Here just averaging pred mean rather than aggregating draws

## By month
pred |> 
  filter(period %% 0.125 == 0) |>
  right_join(shape2) |> 
  filter(!is.na(pred_p)) -> tmp

tmp |> 
  ggplot(aes(fill = pred_p, geometry = SHAPE)) + 
  geom_sf(col = NA) +
  facet_wrap(~month) +
  scale_fill_viridis_c(name = "Predicted NPEV+\nPrevalence") + 
  theme_void() -> p_pred_mth

p_pred_mth
ggsave(here(figdir,"pred_dist_mth.png"),
       height = 6, width = 10, bg = "white")

tmp |> 
  ggplot(aes(fill = sd_p, geometry = SHAPE)) + 
  geom_sf(col = NA) +
  facet_wrap(~month) +
  scale_fill_viridis_c(name = "Predicted NPEV+\nPrevalence (SD)") + 
  theme_void() -> p_pred_mth_sd

p_pred_mth_sd
ggsave(here(figdir,"pred_sd_dist_mth.png"),
       height = 6, width = 10, bg = "white")

## By year

pred |>
  group_by(guid, year) |>
  summarise(pred_mean = mean(pred_p),
            pred_sd = mean(sd_p)) |>
  ungroup() |>
  right_join(shape2) |> 
  filter(!is.na(pred_mean)) -> tmp

tmp |>
  ggplot(aes(fill = pred_mean, geometry = SHAPE)) +
  geom_sf(col = NA) +
  facet_wrap(~year) +
  scale_fill_viridis_c(name = "Predicted NPEV+\nPrevalence") +
  theme_void() -> p_pred_yr
p_pred_yr

ggsave(here(figdir,"pred_dist_yr.png"),
       height = 6, width = 10, bg = "white")

tmp |>
  ggplot(aes(fill = pred_sd, geometry = SHAPE)) +
  geom_sf(col = NA) +
  facet_wrap(~year) +
  scale_fill_viridis_c(name = "Predicted NPEV+\nPrevalence (SD)") +
  theme_void() -> p_pred_yr_sd
p_pred_yr_sd

################################################################################
################################################################################
################################################################################
################################################################################
# Extract predicted prevalence from fitted model
# 
# Draw and summarise model-estimated prevalence for each district/month, 
# including months with zero NPAFP.
# 
# Note: 
# Currently exclude two entire districts absent from fitting data (no recorded 
# NPAFP in this period). These can't be predicted for due to the district-specific
# spline component, which BRMS doesn't let you exclude from prediction...
# Ideally want to predict for these districts based on neighbouring prevalence 
# (through BYM component), province prevalence, and overall time trends. 
################################################################################

pacman::p_load(tidyverse, here, brms)
rm(list = ls(all = TRUE)) 

dir <- "data/Pakistan/analysis/prevalence"
outdir <- "output/prevalence/final"

# Fitting data
fitdata <- read_rds(here(dir,"fitdata.rds"))

# Prediction data (including district:months with zero NPAFP)
preddata <- read_rds(here(dir,"preddata.rds"))

# Shapefiles
shape2 <- read_rds(here(dir,"../shape2.rds"))
W <- read_rds(here(dir, "../W.rds"))

# Exclude districts absent from fitting -------------------------------

# Split out 2 districts missing from fitting data
sub1 <- filter(preddata, !guid %in% unique(fitdata$guid))
sub2 <- filter(preddata, guid %in% unique(fitdata$guid)) |> 
  # For now, match guid levels to fitdata and don't try to predict for zero districts
  mutate(guid = factor(guid, levels = unique(fitdata$guid)))

# Posterior samples of linear predictor -----------------------------------

# Expected counts (relative to NPAFP)
y <- posterior_epred(fit, 
                     newdata = sub2,
                     # All components
                     re_formula = NULL)

## For districts absent from fitting data, predict from model without guid splines
# y1 <- posterior_epred(fit, 
#                       newdata = sub1,
#                       # Exclude s(t, by = guid)
#                       re_formula = ~ (1 | adm1_name) + (1 | t) + s(t, bs = "tp", k = 10) + car(W, gr = guid, type = "bym2") + s(month_of_year, bs = "cc", k = 12),
#                       allow_new_levels = TRUE, sample_new_levels = "uncertainty")
#                       
# y <- epred_draws(fit,
#                  newdata = preddata,
# newdata = filter(preddata, !grepl("0061E338-05BA", guid)) |> mutate(guid = as.character(guid)),
# allow_new_levels = TRUE, sample_new_levels = "gaussian")

# Predicted prevalence
prev <- posterior_linpred(fit, newdata = sub2, transform = T)

# prev1 <- posterior_linpred(fit, newdata = sub1, transform = T,
#                            re_formula = ~ (1 | adm1_name) + (1 | t) +
#                              car(W, gr = guid, type = "bym2") + 
#                              s(month_of_year, bs = "cc", k = 12),
#                            allow_new_levels = T, sample_new_levels = "gaussian")


# Summarise predictions ---------------------------------------------------

y_summ <- as.data.frame(t(apply(y, 2, 
                                function(x) {c(mean(x), 
                                               sd(x), 
                                               quantile(x, 0.025), 
                                               quantile(x, 0.975)[1])}))) |> 
  setNames(c("pred_y", "sd_y","low_y", "hi_y"))


prev_summ <- as.data.frame(t(apply(prev, 2, 
                                   function(x) {
                                     c(mean(x), 
                                       sd(x), 
                                       quantile(x, 0.025), 
                                       quantile(x, 0.975))}))) |> 
  setNames(c("pred_p", "sd_p", "low_p", "hi_p"))

# Bind to main dataset
pred <- cbind(cbind(sub2, y_summ), prev_summ)

# Add draws ---------------------------------------------------------------

pred_draws <- pred |> 
  add_linpred_draws(fit, transform = T) 

# Save predictions --------------------------------------------------------

saveRDS(pred, here(outdir,"pred_prev_npev.rds"))
saveRDS(pred_draws, here(outdir,"draws_prev_npev.rds"))

################################################################################
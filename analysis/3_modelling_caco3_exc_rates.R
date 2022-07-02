#### Aims

# Investigate the factors explaining variation in total carbonate excretion rate

# Potential factors:
# - body mass
# - temperature
# - caudal fin aspect ratio
# - salinity
# - relative intestinal length
# - pCO2: assumed to be at atmospheric equilibrium at all sites (~400 µatm)
# - diet: direct effect not considered as fish were starved

# Also two nuisance variables related to sampling:
# - titration methodology
# - sampling time

#########################################################################################
#### Load data
data_caco3_traits <- readr::read_csv(here::here("data", "derived_data", "data_caco3_traits.csv"))

#########################################################################################
#### EDA

# Check number of replicates per family
length(unique(data_caco3_traits$family)) # 35 families
print(count(data_caco3_traits, family), n = 35)
# Several families have only one or two replicates
# Keep only families with at least 3 replicates
data_caco3_f <- data_caco3_traits %>%
  group_by(family) %>%
  mutate(n = n()) %>%
  filter(n > 2) %>%
  ungroup()

nrow(data_caco3_f) # 175; removed 18 observations
length(unique(data_caco3_f$family)) # 21 families; removed 14

# check for error values
cols <- c(8:12, 14, 15, 18, 22)
lapply(data_caco3_f[, cols], function(x){
  list(nas = data.frame(table(is.na(x))),
       nulls = data.frame(table(x == '')),
       infs = data.frame(table(is.infinite(x))))}
)

cols <- c(9:12, 14, 18, 22)

# Check data distributions
ggpubr::ggarrange(plotlist = lapply(names(data_caco3_f)[cols], function(col) {
  ggplot(data_caco3_f, aes_string(x = col)) + geom_histogram()
}))
# natural-log transform body mass, excretion rate, and RIL
# square-root transform aspect ratio and sampling interval
data_caco3_f <- data_caco3_f %>%
  mutate(log_weight = log(mean_weight_kg),
         log_exc_rate = log(exc_rate_umol_h),
         log_ril = log(ril),
         sqrt_asp_ratio = sqrt(aspect_ratio),
         sqrt_sampling = sqrt(sampling_interval_h),
         # Change family for parrotfishes to Scaridae
         family = if_else(genus == "Scarus", "Scaridae", family))

cols <- c(11, 12, 15, 24:28)

# Check correlations among variables
GGally::ggpairs(data_caco3_f, columns = cols)

# Plot relationships between carbonate excretion rate and potential drivers
ggpubr::ggarrange(plotlist = lapply(names(data_caco3_f)[cols[-5]], function(col) {
  ggplot(data_caco3_f, aes_string(x = col, y = "log_exc_rate")) +
    geom_point() + geom_smooth(method = "lm")
}))

# As expected there is a strong relationship between body mass and carbonate excretion rate (r=0.9).
# All variables are somehow related to carbonate excretion rate, but they are all correlated with body mass
# Check whether the observed effects remain after accounting for body mass
res_er_w <- residuals(lm(log_exc_rate ~ log_weight, data = data_caco3_f))
cols <- c(11, 12, 15, 26, 27, 28)
GGally::ggpairs(cbind(data_caco3_f, res_er_w), columns = c(cols, 29))
ggpubr::ggarrange(plotlist = lapply(names(data_caco3_f)[cols], function(col) {
  ggplot(data_caco3_f, aes_string(x = col, y = "res_er_w")) +
    geom_point() + geom_smooth(method = "lm")
}))
rm(res_er_w, cols)
# Negative effect of RIL and positive effect of salinity, however they are correlated (r=-0.6).
# Moreover, salinity has the highest correlation with body mass (r=0.53).
# We will fit several models that include either salinity or RIL, and then compare them.
# No relationship between residuals and temperature or aspect ratio.
# The methodology used to quantify carbonate excretion rates is correlated to several variables,
# thus will not be included in the model to avoid multicollinearity.
# We will check the correlation with model residuals.

# Fish of the same species and of similar size were sometimes kept in group.
# Estimates of excretion rates obtained from these fish represent average rates
# for an individual of average body mass.
# Thus they should have a higher weight than estimates from individual fishes.
# Check distribution of number of fish per tank
ggplot(data_caco3_f, aes(x = fish_per_tank)) + geom_histogram(binwidth = 1)
summary(data_caco3_f$fish_per_tank)
# Min. 1st Qu.  Median    Mean  3rd Qu.    Max.
# 1.000   1.000   1.000   2.011   2.000  13.000
# The majority of observations have been made on 1 to 3 individuals with an average of 2.
# To avoid overweighting observations made on more than 2 individuals,
# give a weight of 2 to all observations made on multiple individuals.
data_caco3_f_weight <- data_caco3_f %>%
  mutate(weights = if_else(fish_per_tank > 1, 2, 1))

#########################################################################################
#### Analysis

# The following line will fit 36 Bayesian linear and multilevel models to the data and will return a tibble
fits_caco3 <- run_caco3_models(data = data_caco3_f_weight,
                               iter = 4000, warmup = 1000, chains = 4,
                               cores = 4, seed = 98765)

# Add LOO and R2 to each model
fits_caco3$fit <- lapply(fits_caco3$fit, function(x) add_criterion(x, c("loo", "bayes_R2")))

# Compare models
loo_all <- as_tibble(loo_compare(lapply(fits_caco3$fit, function(x) {loo <- x$criteria$loo})), rownames = "model")[, -c(4:7)] %>%
  mutate(across(2:5, round, 2))
fit_vars <- tibble(model = names(fits_caco3$fit)) %>%
  bind_cols(lapply(fits_caco3$fit, extract_fix_ran) %>%
              bind_rows())
loo_all <- left_join(loo_all, fit_vars) %>% select(1, 6, 7, 2:5)
loo_all
readr::write_csv(loo_all, here::here("outputs", "tables", "loo_comparison_all.csv"))

bind_rows(lapply(fits_caco3$fit, function(x) as_tibble(bayes_R2(x))), .id = "model")[,-3] %>%
  arrange(desc(Estimate)) %>%
  mutate(across(2:4, round, 3)) %>% print(n = 36)

# Models with random intercept per family are consistently better than those without.
# Models with RIL are consistently better than those with salinity.
# There is no noteworthy difference in the model fit when adding temperature and aspect ratio,
# but there is also no harm in including them, indeed the model with both has the lowest looic.
# Model fit does not vary when adding sampling interval for all models.
# Use m9_fam_caco3

# Further model checking for "m9_fam_caco3"
m_caco3 <- fits_caco3$fit$m9_fam_caco3
# Summary - Rhat - ESS
summary(m_caco3) # Rhat and ESS ok
# Negative effect of RIL
# Positive effect of temperature
# Slight positive effect of aspect ratio

# Check convergence
plot(m_caco3, N = 4, ask = FALSE) # ok

# Posterior predictive checks
# Overall
pp_check(m_caco3, nsamples = 100) +
  xlim(min(data_caco3_f_weight$log_exc_rate), max(data_caco3_f_weight$log_exc_rate)) # ok
# Mean for each family
pp_check(m_caco3, nsamples = 1000, type = "stat_grouped", stat = "mean", group = "family") # ok

# Get unscaled slopes
round(fixef(m_caco3)[2, ] / sd(data_caco3_f_weight$log_weight), 2) # 0.80 [0.74, 0.87] --> negative allometry
round(fixef(m_caco3)[3, ] / sd(data_caco3_f_weight$log_ril), 2) # -0.57 [-0.89, -0.24]
round(fixef(m_caco3)[4, ] / sd(data_caco3_f_weight$mean_T), 2) # 0.06 [0.01, 0.11]
round(fixef(m_caco3)[5, ] / sd(data_caco3_f_weight$sqrt_asp_ratio), 2) # 0.52 [0.01, 1.03]

# Calculate % variation within the range of RIL (0.33-7.56)
min <- tibble(log_weight = mean(data_caco3_f_weight$log_weight),
              log_ril = min(data_caco3_f_weight$log_ril),
              mean_T = mean(data_caco3_f_weight$mean_T),
              sqrt_asp_ratio = mean(data_caco3_f_weight$sqrt_asp_ratio))
max <- tibble(log_weight = mean(data_caco3_f_weight$log_weight),
              log_ril = max(data_caco3_f_weight$log_ril),
              mean_T = mean(data_caco3_f_weight$mean_T),
              sqrt_asp_ratio = mean(data_caco3_f_weight$sqrt_asp_ratio))

percent_variation(m_caco3, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#   -80.60   -93.77   -52.63

# Calculate % variation within the range of temperature (23-30°C)
min <- min %>% mutate(log_ril = mean(data_caco3_f_weight$log_ril),
                      mean_T = min(data_caco3_f_weight$mean_T))
max <- max %>% mutate(log_ril = mean(data_caco3_f_weight$log_ril),
                      mean_T = max(data_caco3_f_weight$mean_T))

percent_variation(m_caco3, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#    54.74     5.98   116.61

# Q10
Q10(r1 = exp(fitted(m_caco3, newdata = min, re_formula = NA)[1]),
    r2 = exp(fitted(m_caco3, newdata = max, re_formula = NA)[1]),
    t1 = min(data_caco3_f_weight$mean_T),
    t2 = max(data_caco3_f_weight$mean_T))
# 1.79

# Calculate % variation within the range of aspect ratio (0.76-3.3)
min <- min %>% mutate(mean_T = mean(data_caco3_f_weight$mean_T),
                      sqrt_asp_ratio = min(data_caco3_f_weight$sqrt_asp_ratio))
max <- max %>% mutate(mean_T = mean(data_caco3_f_weight$mean_T),
                      sqrt_asp_ratio = max(data_caco3_f_weight$sqrt_asp_ratio))

percent_variation(m_caco3, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#    69.05     0.78   163.85
rm(min, max)

# Correlation between residuals and the methodology, pH, and salinity
res_caco3 <- cbind(data_caco3_f_weight, residuals(m_caco3))

ggpubr::ggarrange(p1 = ggplot(res_caco3, aes(x = log(sampling_interval_h), y = Estimate)) +
                    geom_point() + xlab(italic(ln)~"sampling period (h)") + ylab("residuals") +
                    ggpubr::stat_cor(family = "serif", r.digits = 1, label.x.npc = "left", label.y.npc = 0.05),
                  p2 = ggplot(res_caco3, aes(x = method, y = Estimate)) +
                    geom_boxplot() + xlab("Methodology") + ylab("") +
                    ggpubr::stat_compare_means(family = "serif", label.x = 0.7, label.y.npc = 0.04),
                  p3 = ggplot(res_caco3, aes(x = mean_S, y = Estimate)) +
                    geom_point() + xlab("Salinity") + ylab("") +
                    ggpubr::stat_cor(family = "serif", label.x.npc = "left", label.y.npc = 0.05),
                  labels = letters[1:3], nrow = 1)

ggsave(here::here("outputs", "figures", "residuals_vs_other_vars.png"),
       width = 20, height = 8, unit = "cm", dpi = 600, type = "cairo")

# Test difference in residual variance between methods
car::leveneTest(Estimate ~ as.factor(method), data = res_caco3)
#        Df F value    Pr(>F)
# group   1  22.819 3.779e-06 ***
#       173

rm(res_caco3)

# We can model sigma as a function of titration method (distributional regression)
m_caco3_sigma <- update(m_caco3, bf(~ . , sigma ~ method),
                        newdata = data_caco3_f_weight, family = student(),
                        iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 98765,
                        file = here::here("outputs", "models", "m_caco3_sigma"))

m_caco3_sigma <- add_criterion(m_caco3_sigma, c("loo", "bayes_R2"))
loo_sigma <- as_tibble(loo_compare(m_caco3, m_caco3_sigma), rownames = "model")[, -c(4:7)] %>%
  mutate(across(2:5, round, 2))
fit_vars <- tibble(model = c("m_caco3", "m_caco3_sigma")) %>%
  bind_cols(lapply(list(m_caco3, m_caco3_sigma), extract_fix_ran) %>%
              bind_rows())
loo_sigma <- left_join(loo_sigma, fit_vars) %>% select(1, 6:8, 2:5)
loo_sigma
# model         fixed                            random sigma  elpd_diff  se_diff    looic      se_looic
# m_caco3_sigma log(M) + log(RIL) + T + sqrt(AR) family method   0.00     0.00       451.27     36.70
# m_caco3       log(M) + log(RIL) + T + sqrt(AR) family -      -16.92     7.61       485.11     36.03
readr::write_csv(loo_sigma, here::here("outputs", "tables", "loo_comparison_sigma.csv"))

rm(fit_vars)

summary(m_caco3_sigma)
plot(m_caco3_sigma, N = 4, ask = FALSE)
pp_check(m_caco3_sigma, nsamples = 1000) +
  xlim(min(data_caco3_f_weight$log_exc_rate), max(data_caco3_f_weight$log_exc_rate)) # ok
pp_check(m_caco3_sigma, nsamples = 1000, type = "stat_grouped", stat = "mean", group = "family") # ok

# Test if sigma differ from 0
hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_methodsingle) = 0")
hypothesis(m_caco3_sigma, hyp)
# 0.39 [0.31, 0.47]
# 0.80 [0.63, 0.98]

# Test if sigma of direct titration is larger than sigma of double titration
hyp <- "exp(sigma_Intercept + sigma_methodsingle) > exp(sigma_Intercept)"
(hyp <- hypothesis(m_caco3_sigma, hyp))
plot(hyp, chars = NULL)
# Yes: 0.4 [0.27, 0.55]
rm(hyp)

# Get unscaled slopes
round(fixef(m_caco3_sigma)[3, ] / sd(data_caco3_f_weight$log_weight), 2) # 0.78[0.72, 0.83]
round(fixef(m_caco3_sigma)[4, ] / sd(data_caco3_f_weight$log_ril), 2) # -0.59[-0.92, -0.27]
round(fixef(m_caco3_sigma)[5, ] / sd(data_caco3_f_weight$mean_T), 2) # 0.05[0.01, 0.10]
round(fixef(m_caco3_sigma)[6, ] / sd(data_caco3_f_weight$sqrt_asp_ratio), 2) # 0.71[0.25, 1.17]

# Calculate % variation within the range of RIL (0.33-7.56)
min <- tibble(log_weight = mean(data_caco3_f_weight$log_weight),
              log_ril = min(data_caco3_f_weight$log_ril),
              mean_T = mean(data_caco3_f_weight$mean_T),
              sqrt_asp_ratio = mean(data_caco3_f_weight$sqrt_asp_ratio),
              method = "double")
max <- tibble(log_weight = mean(data_caco3_f_weight$log_weight),
              log_ril = max(data_caco3_f_weight$log_ril),
              mean_T = mean(data_caco3_f_weight$mean_T),
              sqrt_asp_ratio = mean(data_caco3_f_weight$sqrt_asp_ratio),
              method = "double")

percent_variation(m_caco3_sigma, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#   -82.05   -94.44   -57.16

# Calculate % variation within the range of temperature (23-30°C)
min <- min %>% mutate(log_ril = mean(data_caco3_f_weight$log_ril),
                      mean_T = min(data_caco3_f_weight$mean_T))
max <- max %>% mutate(log_ril = mean(data_caco3_f_weight$log_ril),
                      mean_T = max(data_caco3_f_weight$mean_T))

percent_variation(m_caco3_sigma, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#    48.46     4.01   106.50

# Q10
q10 <- tibble(r1 = exp(fitted(m_caco3_sigma, newdata = min, re_formula = NA, summary = FALSE)),
              r2 = exp(fitted(m_caco3_sigma, newdata = max, re_formula = NA, summary = FALSE)),
              t1 = min(data_caco3_f_weight$mean_T),
              t2 = max(data_caco3_f_weight$mean_T))

q10 <- Q10(q10$r1, q10$r2, q10$t1, q10$t2)
round(mean(q10), 2)
# 1.74
quantile(q10, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>% round(2)
# 2.5%   25%   50%   75% 97.5%
# 1.06  1.43  1.69  2.00  2.73
rm(q10)

# Calculate % variation within the range of aspect ratio (0.76-3.3)
min <- min %>% mutate(mean_T = mean(data_caco3_f_weight$mean_T),
                      sqrt_asp_ratio = min(data_caco3_f_weight$sqrt_asp_ratio))
max <- max %>% mutate(mean_T = mean(data_caco3_f_weight$mean_T),
                      sqrt_asp_ratio = max(data_caco3_f_weight$sqrt_asp_ratio))

percent_variation(m_caco3_sigma, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#   100.47    26.76   201.64

rm(min, max)

# Check number of observations that fall within the 95% CI of the predictions
predict(m_caco3_sigma) %>%
  as_tibble() %>%
  bind_cols(data_caco3_f_weight) %>%
  mutate(check = ifelse(log_exc_rate >= Q2.5 & log_exc_rate <= Q97.5, 1, 0)) %>%
  count(check) %>%
  mutate("%" = n/sum(n)*100)
#  check     n   `%`
#     0     7     4
#     1   168    96

rm(fits_caco3)
gc()

#########################################################################################

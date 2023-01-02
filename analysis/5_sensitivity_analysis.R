#### Aims

# Sensitivity analysis
# Repeat analysis after correcting excretion rates from Palau
# Use average difference between titration methods (1.08)

#########################################################################################
#### Prepare data

# Total carbonate excretion rate
data_caco3_corrected <- data_caco3_f_weight %>%
  mutate(exc_rate_umol_h = if_else(region == "Palau",
                                   exc_rate_umol_h / 1.08,
                                   exc_rate_umol_h),
         log_exc_rate = log(exc_rate_umol_h))

# Carbonate composition
caco3_comp_corrected <- caco3_comp %>%
  mutate(exc_rate_umol_h = if_else(region == "Palau",
                                   exc_rate_umol_h / 1.08,
                                   exc_rate_umol_h),
         L_umol_h = exc_rate_umol_h * (L / 100),
         H_umol_h = exc_rate_umol_h * (H / 100),
         AR_umol_h = exc_rate_umol_h * (AR / 100),
         M_umol_h = exc_rate_umol_h * (M / 100),
         AC_umol_h = exc_rate_umol_h * (AC / 100),
         B_umol_h = exc_rate_umol_h * (B / 100),
         N_umol_h = exc_rate_umol_h * (N / 100))

#########################################################################################
#### Models

# Total carbonate excretion rate
m_caco3_sigma_corrected <- update(m_caco3_sigma,
                                  newdata = data_caco3_corrected,
                                  family = student(),
                                  iter = 4000, warmup = 1000,
                                  chains = 4, cores = 4, seed = 98765)

summary(m_caco3_sigma_corrected)

plot(m_caco3_sigma_corrected, N = 4, ask = FALSE)
pp_check(m_caco3_sigma_corrected, nsamples = 1000) +
  xlim(min(data_caco3_corrected$log_exc_rate),
       max(data_caco3_corrected$log_exc_rate))
pp_check(m_caco3_sigma_corrected, nsamples = 1000,
         type = "stat_grouped", stat = "mean", group = "family")

# Test if sigma differ from 0
hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_methodsingle) = 0")
hypothesis(m_caco3_sigma_corrected, hyp)
# 0.4 [0.32, 0.47]
# 0.80 [0.63, 0.98]

# Test if sigma of direct titration is larger than sigma of double titration
hyp <- "exp(sigma_Intercept + sigma_methodsingle) > exp(sigma_Intercept)"
(hyp <- hypothesis(m_caco3_sigma_corrected, hyp))
plot(hyp, chars = NULL)
# Yes: 0.4 [0.27, 0.55]
rm(hyp)

# Get unscaled slopes
round(fixef(m_caco3_sigma_corrected)[3, ] / sd(data_caco3_corrected$log_weight), 2) # 0.78[0.72, 0.84]
round(fixef(m_caco3_sigma_corrected)[4, ] / sd(data_caco3_corrected$log_ril), 2) # -0.61[-0.94, -0.27]
round(fixef(m_caco3_sigma_corrected)[5, ] / sd(data_caco3_corrected$mean_T), 2) # 0.05[0.01, 0.10]
round(fixef(m_caco3_sigma_corrected)[6, ] / sd(data_caco3_corrected$sqrt_asp_ratio), 2) # 0.71[0.24, 1.17]

# Calculate % variation within the range of RIL (0.33-7.56)
min <- tibble(log_weight = mean(data_caco3_corrected$log_weight),
              log_ril = min(data_caco3_corrected$log_ril),
              mean_T = mean(data_caco3_corrected$mean_T),
              sqrt_asp_ratio = mean(data_caco3_corrected$sqrt_asp_ratio),
              method = "double")
max <- tibble(log_weight = mean(data_caco3_corrected$log_weight),
              log_ril = max(data_caco3_corrected$log_ril),
              mean_T = mean(data_caco3_corrected$mean_T),
              sqrt_asp_ratio = mean(data_caco3_corrected$sqrt_asp_ratio),
              method = "double")

percent_variation(m_caco3_sigma_corrected, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#   -83.01   -94.74   -57.70

# Calculate % variation within the range of temperature (23-30°C)
min <- min %>% mutate(log_ril = mean(data_caco3_corrected$log_ril),
                      mean_T = min(data_caco3_corrected$mean_T))
max <- max %>% mutate(log_ril = mean(data_caco3_corrected$log_ril),
                      mean_T = max(data_caco3_corrected$mean_T))

percent_variation(m_caco3_sigma_corrected, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#    43.65     0.30    99.80

# Q10
q10 <- tibble(r1 = exp(fitted(m_caco3_sigma_corrected, newdata = min, re_formula = NA, summary = FALSE)),
              r2 = exp(fitted(m_caco3_sigma_corrected, newdata = max, re_formula = NA, summary = FALSE)),
              t1 = min(data_caco3_corrected$mean_T),
              t2 = max(data_caco3_corrected$mean_T))

q10 <- Q10(q10$r1, q10$r2, q10$t1, q10$t2)
round(mean(q10), 2)
# 1.66
quantile(q10, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>% round(2)
# 2.5%   25%   50%   75% 97.5%
# 1.00  1.37  1.62  1.89  2.60
rm(q10)

# Calculate % variation within the range of aspect ratio (0.76-3.3)
min <- min %>% mutate(mean_T = mean(data_caco3_corrected$mean_T),
                      sqrt_asp_ratio = min(data_caco3_corrected$sqrt_asp_ratio))
max <- max %>% mutate(mean_T = mean(data_caco3_corrected$mean_T),
                      sqrt_asp_ratio = max(data_caco3_corrected$sqrt_asp_ratio))

percent_variation(m_caco3_sigma_corrected, x1 = min, x2 = max, exp = TRUE, robust = FALSE) %>% round(2)
# Estimate     Q2.5    Q97.5
#   100.18    25.60   201.56

rm(min, max)

# Carbonate composition
m3_caco3_comp_corrected <- brm(bf_3,
                               data = caco3_comp_corrected,
                               family = hurdle_lognormal(),
                               prior = prior_multivariate,
                               chains = 3, cores = 3,
                               iter = 4000, warmup = 2000,
                               seed = 98765,
                               control = list(adapt_delta = 0.99, max_treedepth = 15))

# Summary - Rhat - ESS
summary(m3_caco3_comp_corrected)

# Posterior predictive checks
ggpubr::ggarrange(pp_check_mult(m3_caco3_comp_corrected, "Humolh", "HMC", "H_umol_h"),
                  pp_check_mult(m3_caco3_comp_corrected, "Lumolh", "LMC", "L_umol_h"),
                  pp_check_mult(m3_caco3_comp_corrected, "Mumolh", "MHC", "M_umol_h"),
                  pp_check_mult(m3_caco3_comp_corrected, "ARumolh", "Aragonite", "AR_umol_h"),
                  pp_check_mult(m3_caco3_comp_corrected, "ACumolh", "ACMC", "AC_umol_h"),
                  common.legend = TRUE)

# Check if the model retrieves the correct proportion of zeroes
ggpubr::ggarrange(plotlist = Map(function(x, y) {
  pp_check(m3_caco3_comp_corrected, resp = x, nsamples = 1000,
           type = "stat", stat = "prop_zero", binwidth = 0.005) +
    ggtitle(y)
  },
  x = c("Lumolh", "ARumolh", "Humolh", "Mumolh", "ACumolh"),
  y = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")),
  ncol = 3, nrow = 2, legend = "none")

# Get unscaled effect of body mass
round(fixef(m3_caco3_comp_corrected) %>%
        as.data.frame() %>%
        filter(endsWith(rownames(.), "_scalelog_weight")) / sd(caco3_comp_corrected$log_weight), 2)
# No difference

rm(bf_1, bf_2, bf_3, bf_AC, bf_AR, bf_H, bf_L, bf_M)
gc()

#########################################################################################

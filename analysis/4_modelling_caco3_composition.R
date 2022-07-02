#### Aims

# Investigate the drivers of carbonate composition

#########################################################################################
#### Load data
caco3_comp <- readr::read_csv(here::here("data", "raw_data", "caco3_composition.csv"))

#########################################################################################
# Check names
fishflux::name_errors(unique(caco3_comp$species))
# Inaccurate species names found:
# "Apogon limenus"   "Blenny sp."     "Flounder sp."    "Eucinostomus sp."
# "Haemulon sp."    "Pelates sexlineatus"

# Correct names
caco3_comp <- caco3_comp %>%
  mutate(species = recode(species,
                          "Apogon limenus" = "Ostorhinchus limenus",
                          "Pelates sexlineatus" = "Helotes sexlineatus"),
         genus = recode(genus,
                        "Apogon" = "Ostorhinchus",
                        "Pelates" = "Helotes"))

# Change parrotfish family to Scaride to join the datasets
caco3_comp <- caco3_comp %>%
  mutate(family = if_else(genus == "Scarus", "Scaridae", family))

setdiff(unique(data_caco3_f_weight$species), unique(caco3_comp$species))
# Six species don't have compositional data

# Join compositional data to the excretion rates
caco3_comp <- left_join(data_caco3_f_weight, caco3_comp[, -11])

# Compute excretion rate for each polymorph and remove NAs
caco3_comp <- caco3_comp %>%
  mutate(L_umol_h = exc_rate_umol_h * (L / 100),
         H_umol_h = exc_rate_umol_h * (H / 100),
         AR_umol_h = exc_rate_umol_h * (AR / 100),
         M_umol_h = exc_rate_umol_h * (M / 100),
         AC_umol_h = exc_rate_umol_h * (AC / 100),
         B_umol_h = exc_rate_umol_h * (B / 100),
         N_umol_h = exc_rate_umol_h * (N / 100)) %>%
  filter(!is.na(H))
# Removed 8 observations without compositional data

# Check number of replicates per family
length(unique(caco3_comp$family)) # 22 families
print(count(caco3_comp, family), n = 22) # all families have 3 or more replicates

# Check data distribution
ggpubr::ggarrange(plotlist = lapply(names(caco3_comp)[37:43], function(col) {
  ggplot(caco3_comp, aes_string(x = col)) + geom_histogram()
}))
# Only one fish in the dataset produces nesquehonite (N)
# Most values are 0
# % of 0 for each polymorph
sapply(names(caco3_comp)[37:43], function(col) {
  n0 = round(sum(caco3_comp[,col] == 0) / nrow(caco3_comp) * 100, 1)
  })
# L_umol_h  H_umol_h AR_umol_h  M_umol_h AC_umol_h  B_umol_h  N_umol_h
# 81.4      13.8      74.3      83.2      53.3      63.5      99.4

# Do not consider nesquehonite in the analysis as it is a minor component (5%) of one species only
# Do not consider brucite as it is magnesium hydroxide (Mg(OH)2) and a minor component (<5%) when present

#########################################################################################
#### Analysis

# Bayesian multivariate hurdle model
resp <- c("ACumolh", "Humolh", "Lumolh", "Mumolh", "ARumolh")
prior_multivariate <- c(
  set_prior("normal(0, 5)", class = "Intercept", resp = resp),
  set_prior("normal(0, 5)", class = "b", resp = resp),
  set_prior("normal(0, 5)", class = "sd", resp = resp),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = resp),
  set_prior("normal(0, 5)", class = "b", dpar = "hu", resp = resp),
  set_prior("normal(0, 5)", class = "sd", dpar = "hu", resp = resp),
  set_prior("student_t(3, 0, 2.5)", class = "Intercept", dpar = "sigma", resp = resp),
  set_prior("normal(0, 5)", class = "b", dpar = "sigma", resp = resp)
  )

# Full model:
# - apply to each response variable the model selected for the total excretion rates
# - model the hurdle probability (hu) with the group-level effect of family,
#   and the fixed effects of temperature and relative intestinal length,
#   to test whether they influence the probability of producing certain polymorphs
# - to account for the correlation among response variables,
#   the group-level effects are modelled correlated
bf_1 <- bf(mvbind(L_umol_h, H_umol_h, AR_umol_h, M_umol_h, AC_umol_h)|weights(weights) ~
             scale(log_weight) + scale(log_ril) + scale(mean_T) + scale(sqrt_asp_ratio) + (1|p|family),
           sigma ~ method,
           hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))
m1_caco3_comp <- brm(bf_1, data = caco3_comp, family = hurdle_lognormal(),
                     prior = prior_multivariate, chains = 3, cores = 3,
                     iter = 4000, warmup = 2000, seed = 98765,
                     control = list(adapt_delta = 0.999, max_treedepth = 15),
                     file = here::here("outputs", "models", "m1_caco3_comp"))
# Took about 35 minutes on this machine

# Add LOO and R2
m1_caco3_comp <- add_criterion(m1_caco3_comp, c("loo", "bayes_R2"))

# Summary - Rhat - ESS
summary(m1_caco3_comp)

# Bayesian R2
round(bayes_R2(m1_caco3_comp), 2)
#           Estimate Est.Error Q2.5 Q97.5
# R2Lumolh      0.63      0.26 0.11  0.96
# R2Humolh      0.54      0.16 0.28  0.86
# R2ARumolh     0.88      0.14 0.50  0.99
# R2Mumolh      0.58      0.14 0.41  0.89
# R2ACumolh     0.76      0.24 0.16  0.99

# Posterior predictive checks
pp_check_mult <- function(model, resp, polymorph, colname) {
  pp_check(model, resp = resp, nsamples = 100) +
    xlab(paste(polymorph, "exc. rate (&mu;mol h<sup>-1</sup>)")) +
    xlim(0, max(caco3_comp[, colname])) +
    theme(axis.title.x = element_markdown())
}
ggpubr::ggarrange(pp_check_mult(m1_caco3_comp, "Humolh", "HMC", "H_umol_h"),
                  pp_check_mult(m1_caco3_comp, "Lumolh", "LMC", "L_umol_h"),
                  pp_check_mult(m1_caco3_comp, "Mumolh", "MHC", "M_umol_h"),
                  pp_check_mult(m1_caco3_comp, "ARumolh", "Aragonite", "AR_umol_h"),
                  pp_check_mult(m1_caco3_comp, "ACumolh", "ACMC", "AC_umol_h"),
                  common.legend = TRUE)

# Check if the model retrieves the correct proportion of zeroes
prop_zero <- function(x) mean(x == 0)
ggpubr::ggarrange(plotlist = Map(function(x, y) {
  pp_check(m1_caco3_comp, resp = x, nsamples = 1000,
           type = "stat", stat = "prop_zero", binwidth = 0.005) +
    ggtitle(y)
}, x = resp, y = c("ACMC", "HMC", "LMC", "MHC", "Aragonite")),
ncol = 3, nrow = 2, legend = "none")

# Different formula for each response variable
# Remove variables with large errors on lognormal part (error > estimate),
# particularly T, AR, RIL for LMC and T for MHC
bf_L <- bf(L_umol_h|weights(weights) ~ scale(log_weight) + (1|p|family),
           sigma ~ method,
           hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))

bf_H <- bf(H_umol_h|weights(weights) ~ scale(log_weight) + scale(log_ril) +
             scale(mean_T) + scale(sqrt_asp_ratio) + (1|p|family),
           sigma ~ method,
           hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))

bf_AR <- bf(AR_umol_h|weights(weights) ~ scale(log_weight) + scale(log_ril) +
              scale(mean_T) + scale(sqrt_asp_ratio) + (1|p|family),
            sigma ~ method,
            hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))

bf_M <- bf(M_umol_h|weights(weights) ~ scale(log_weight) + scale(log_ril) +
             scale(sqrt_asp_ratio) + (1|p|family),
           sigma ~ method,
           hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))

bf_AC <- bf(AC_umol_h|weights(weights) ~ scale(log_weight) + scale(log_ril) +
              scale(mean_T) + scale(sqrt_asp_ratio) + (1|p|family),
            sigma ~ method,
            hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))

bf_2 <- bf_L + bf_H + bf_AR + bf_M + bf_AC
m2_caco3_comp <- brm(bf_2, data = caco3_comp, family = hurdle_lognormal(),
                     prior = prior_multivariate, chains = 3, cores = 3,
                     iter = 4000, warmup = 2000, seed = 98765,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     file = here::here("outputs", "models", "m2_caco3_comp"))

# Add LOO and R2
m2_caco3_comp <- add_criterion(m2_caco3_comp, c("loo", "bayes_R2"))

# Summary - Rhat - ESS
summary(m2_caco3_comp)

# Compare models
loo(m1_caco3_comp, m2_caco3_comp) # no difference

# Bayesian R2
round(bayes_R2(m2_caco3_comp), 2)
#           Estimate Est.Error Q2.5 Q97.5
# R2Lumolh      0.68      0.25 0.14  0.96
# R2Humolh      0.54      0.15 0.27  0.85
# R2ARumolh     0.88      0.14 0.50  0.99
# R2Mumolh      0.62      0.16 0.38  0.90
# R2ACumolh     0.76      0.25 0.16  0.99

# Posterior predictive checks
ggpubr::ggarrange(pp_check_mult(m2_caco3_comp, "Humolh", "HMC", "H_umol_h"),
                  pp_check_mult(m2_caco3_comp, "Lumolh", "LMC", "L_umol_h"),
                  pp_check_mult(m2_caco3_comp, "Mumolh", "MHC", "M_umol_h"),
                  pp_check_mult(m2_caco3_comp, "ARumolh", "Aragonite", "AR_umol_h"),
                  pp_check_mult(m2_caco3_comp, "ACumolh", "ACMC", "AC_umol_h"),
                  common.legend = TRUE)

# Check if the model retrieves the correct proportion of zeroes
ggpubr::ggarrange(plotlist = Map(function(x, y) {
  pp_check(m2_caco3_comp, resp = x, nsamples = 1000,
           type = "stat", stat = "prop_zero", binwidth = 0.005) +
    ggtitle(y)
  }, x = resp, y = c("ACMC", "HMC", "LMC", "MHC", "Aragonite")),
 ncol = 3, nrow = 2, legend = "none")

# Remove the distributional part of the model from MHC
bf_M <- bf(M_umol_h|weights(weights) ~ scale(log_weight) + scale(log_ril) +
             scale(sqrt_asp_ratio) + (1|p|family),
           hu ~ scale(log_ril) + scale(mean_T) + (1|p|family))

prior_multivariate <- prior_multivariate[-c(34, 39),]

bf_3 <- bf_L + bf_H + bf_AR + bf_M + bf_AC
m3_caco3_comp <- brm(bf_3, data = caco3_comp, family = hurdle_lognormal(),
                     prior = prior_multivariate, chains = 3, cores = 3,
                     iter = 4000, warmup = 2000, seed = 98765,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     file = here::here("outputs", "models", "m3_caco3_comp"))

# Add LOO and R2
m3_caco3_comp <- add_criterion(m3_caco3_comp, c("loo", "bayes_R2"))

# Summary - Rhat - ESS
summary(m3_caco3_comp)

# Compare models
loo(m1_caco3_comp, m2_caco3_comp, m3_caco3_comp) # no difference

# Bayesian R2
round(bayes_R2(m3_caco3_comp), 2)
#           Estimate Est.Error Q2.5 Q97.5
# R2Lumolh      0.69      0.24 0.14  0.96
# R2Humolh      0.54      0.16 0.26  0.85
# R2ARumolh     0.88      0.14 0.49  0.99
# R2Mumolh      0.66      0.16 0.35  0.90
# R2ACumolh     0.76      0.24 0.18  0.99

# Posterior predictive checks
ggpubr::ggarrange(pp_check_mult(m3_caco3_comp, "Humolh", "HMC", "H_umol_h"),
                  pp_check_mult(m3_caco3_comp, "Lumolh", "LMC", "L_umol_h"),
                  pp_check_mult(m3_caco3_comp, "Mumolh", "MHC", "M_umol_h"),
                  pp_check_mult(m3_caco3_comp, "ARumolh", "Aragonite", "AR_umol_h"),
                  pp_check_mult(m3_caco3_comp, "ACumolh", "ACMC", "AC_umol_h"),
                  common.legend = TRUE)

# Check if the model retrieves the correct proportion of zeroes
ggpubr::ggarrange(plotlist = Map(function(x, y) {
  pp_check(m3_caco3_comp, resp = x, nsamples = 1000,
           type = "stat", stat = "prop_zero", binwidth = 0.005) +
    ggtitle(y)
  },
  x = c("Lumolh", "ARumolh", "Humolh", "Mumolh", "ACumolh"),
  y = c("LMC", "Aragonite", "HMC", "MHC", "ACMC")),
  ncol = 3, nrow = 2, legend = "none")

# Check number of observations that fall within the 50% CI of the predictions
Map(function(x, y) {
  predict(m3_caco3_comp, robust = TRUE, resp = x, probs = c(0.25, 0.75)) %>%
    as_tibble() %>%
    bind_cols(caco3_comp) %>%
    mutate(check = ifelse(.[,y] >= Q25 & .[,y] <= Q75, 1, 0)) %>%
    count(check) %>%
    mutate("%" = n/sum(n)*100)
}, x = resp, y = c("AC_umol_h", "H_umol_h", "L_umol_h", "M_umol_h", "AR_umol_h")) %>%
  bind_rows(.id = "mineral")
#   mineral check[,1]     n   `%`
# ACumolh         0    39 23.4
# ACumolh         1   128 76.6
# Humolh          0    61 36.5
# Humolh          1   106 63.5
# Lumolh          0    13  7.78
# Lumolh          1   154 92.2
# Mumolh          0    12  7.19
# Mumolh          1   155 92.8
# ARumolh         0    21 12.6
# ARumolh         1   146 87.4

# Repeat with the 90% CI
Map(function(x, y) {
  predict(m3_caco3_comp, robust = TRUE, resp = x, probs = c(0.05, 0.95)) %>%
    as_tibble() %>%
    bind_cols(caco3_comp) %>%
    mutate(check = ifelse(.[,y] >= Q5 & .[,y] <= Q95, 1, 0)) %>%
    count(check) %>%
    mutate("%" = n/sum(n)*100)
}, x = resp, y = c("AC_umol_h", "H_umol_h", "L_umol_h", "M_umol_h", "AR_umol_h")) %>%
  bind_rows(.id = "mineral")
#mineral check[,1]     n    `%`
# ACumolh         0     5  2.99
# ACumolh         1   162 97.0
# Humolh          0     6  3.59
# Humolh          1   161 96.4
# Lumolh          0     2  1.20
# Lumolh          1   165 98.8
# Mumolh          0     1  0.599
# Mumolh          1   166 99.4
# ARumolh         0     1  0.599
# ARumolh         1   166 99.4

# Get unscaled effect of body mass
round(fixef(m3_caco3_comp) %>%
         as.data.frame() %>%
         filter(endsWith(rownames(.), "_scalelog_weight")) / sd(caco3_comp$log_weight), 2)
#                         Estimate Est.Error Q2.5 Q97.5
# Lumolh_scalelog_weight      0.89      0.16 0.56  1.19
# Humolh_scalelog_weight      0.74      0.03 0.68  0.80
# ARumolh_scalelog_weight     1.14      0.14 0.84  1.39
# Mumolh_scalelog_weight      1.40      0.19 1.04  1.78
# ACumolh_scalelog_weight     0.89      0.13 0.64  1.14

rm(m1_caco3_comp, m2_caco3_comp)
gc()

#########################################################################################

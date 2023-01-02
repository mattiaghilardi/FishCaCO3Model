#### Aims

# Extrapolation of intestinal length of fishes in the carbonate dataset

# 1 - Fit Bayesian phylogenetic model with data published in
#     Ghilardi et al. (2021) Phylogeny, body morphology, and trophic level
#     shape intestinal traits in coral reef fishes. Ecol. Evol.
# 2 - Validate model
# 3 - Extrapolate intestinal length for taxa in the carbonate dataset,
#     both fishes identified at the species or genus level, using
#     ancestral state reconstruction of the phylogenetic effect

#########################################################################################
#### Load data
int_moor <- readr::read_csv(here::here("data", "raw_data", "intestine_fish_Moorea.csv"))
data_caco3_FB_traits <- readr::read_csv(here::here("data", "derived_data", "data_caco3_FB_traits.csv"))

#########################################################################################
#### Retrieve species traits from FishBase

# Check if species names match those in FishBase
fishflux::name_errors(unique(int_moor$species))
# Inaccurate species names found:
# "Ostracion cubicus" "Myripristis speA" "Centropyge loricula" "Acanthurus nigros"

# We can drop the second but we retain the other three
int_moor <- filter(int_moor, species != "Myripristis speA")

# Correct name for "Ostracion cubicus", "Centropyge loricula", and "Acanthurus nigros"
rfishbase::validate_names(c("Ostracion cubicus", "Centropyge loricula", "Acanthurus nigros"))
# "Ostracion cubicum"  "Centropyge loriculus"  NA
# "Acanthurus nigros" is considered different from "A. nigroris" following
# Randall et al. (2011), but it is not present in FishBase
# Rename as "A. nigroris"

# Change these names in the dataset
int_moor <- int_moor %>%
  mutate(species = recode(species,
                          "Ostracion cubicus" = "Ostracion cubicum",
                          "Centropyge loricula" = "Centropyge loriculus",
                          "Acanthurus nigros" = "Acanthurus nigroris"))

# Check
fishflux::name_errors(unique(int_moor$species)) # ok

# TROPHIC LEVEL
troph_moor <- trophic_level_FB(unique(int_moor$species), type = "food items") %>%
  dplyr::rename(level_troph = taxon_level)
filter(troph_moor, is.na(trophic_level)) # no NAs

# ELONGATION
elon_moor <- elongation_FB(unique(int_moor$species)) %>%
  dplyr::rename(level_elon = taxon_level)
filter(elon_moor, is.na(elongation)) # no NAs

# Join traits to the dataset
int_moor <- left_join(int_moor, troph_moor) %>% left_join(elon_moor)

rm(troph_moor, elon_moor)

#########################################################################################
#### Prepare data and phylogenetic reletedness matrix

# Data transformation
int_moor <- int_moor %>%
  mutate(il_log = log(intestine_length),
         sl_log = log(sl),
         elon_log = log(elongation))

# Select species with at least 3 replicates
int_moor_f <- int_moor %>%
  group_by(species) %>%
  mutate(n = n()) %>%
  filter(n > 2) %>%
  ungroup()

## Download the phylogeny from The Fish Tree of Life,
## and create a covariance matrix for Bayesian phylogenetic model in brms

# Check species names
check_name_fishtree(unique(int_moor_f$species), sampled = TRUE)
# Some species name are incorrect or not present in the Fish Tree of Life
# Some species do not have genetic data
# $name_error
# [1] "Ostracion cubicum"
# $not_sampled
# [1] "Epinephelus hexagonatus" "Ostracion cubicum" "Hemiramphus depauperatus" "Exallias brevis"
check_name_fishtree("Ostracion cubicus") #ok
# Rename again this species as Ostracion cubicus
int_moor_f <- mutate(int_moor_f, species = recode(species,
                                                  "Ostracion cubicum" = "Ostracion cubicus"))

# Use function "fishtree_complete_phylogeny" (100 trees)
tree_int <- fishtree::fishtree_complete_phylogeny(species = unique(int_moor_f$species))

# Phylogenetic correlation matrix
corphy_int <- phylo_cor(tree_int)
corphy_int_m <- corphy_int$mean

# Add column "phylo" to the dataset
int_moor_f$phylo <- gsub(" ", "_", int_moor_f$species)

#########################################################################################
#### Build predictive model for intestinal length

## Bayesian phylogenetic hierarchical linear model in brms

# Priors
prior_int <- c(prior(normal(0, 10), class = Intercept),
               prior(normal(0, 5), class = b),
               prior(normal(0, 5), class = sd),
               prior(gamma(2, 0.1), class = nu))

# Model
m_intestine <- brm(il_log ~ scale(sl_log) + scale(trophic_level) + scale(elon_log) + (1 |gr(phylo, cov = corphy_int_m)),
                   data = int_moor_f, family = student(), data2 = list(corphy_int_m = corphy_int_m),
                   prior = prior_int, warmup = 1000, iter = 4000, chains = 4, cores = 4,
                   seed = 98765, file = "outputs/models/m_intestine")

# Summary and checks
summary(m_intestine) # ok
pp_check(m_intestine, type = "dens_overlay", nsamples = 100) +
  xlim(min(int_moor_f$il_log), max(int_moor_f$il_log)) # good
plot(m_intestine, N = 7) # good
bayes_R2(m_intestine) # 0.92 [0.92, 0.93]

# Get unscaled slopes
# standard length
round(fixef(m_intestine)[2,]/sd(int_moor_f$sl_log), 2) # -0.92 [-0.85, -0.99]
# trophic level
round(fixef(m_intestine)[3,]/sd(int_moor_f$trophic_level), 2) # -0.35 [-0.49, -0.22]
# elongation
round(fixef(m_intestine)[4,]/sd(int_moor_f$elon_log), 2) # -0.76 [-1.02, -0.51]

# Plot observed vs predicted
predict(m_intestine) %>%
  as_tibble() %>%
  bind_cols(int_moor_f) %>%
  plot_obs_vs_pred(pred = "Estimate", obs = "il_log",
                   xlab = "Predicted ln int. length (mm)",
                   ylab = "Observed ln int. length (mm)",
                   point_size = 1.5,
                   text_size = 4)

# The model is very accurate in predicting the observed data used to train the model (R2=0.93, b=1, a=-0.056)

ggsave(here::here("outputs", "figures", "obs_vs_pred_IL.png"),
       width = 10, height = 7, units = "cm", dpi = 600, type = "cairo")

# Check number of observations that fall within the 95% CI of the predictions
predict(m_intestine) %>%
  as_tibble() %>%
  bind_cols(int_moor_f) %>%
  mutate(check = ifelse(il_log >= Q2.5 & il_log <= Q97.5, 1, 0)) %>%
  count(check) %>%
  mutate("%" = n/sum(n)*100)
# check     n    `%`
#     0    43  3.56
#     1  1165  96.4

rm(corphy_int, corphy_int_m, tree_int, prior_int)

#########################################################################################
#### Model validation

# Exclude one species per time from the model and predict intestinal length for all individuals
# of this species in the dataset.

# To account for the phylogeny we first need to estimate the phylogenetic effect for the unobserved taxa,
# which can be done through ancestral state reconstruction (see function "phylo_effect()",
# following Parravicini et al. 2020 https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000702),
# then the intestinal length of unobserved taxa can be predicted (see function "intestinal_length()").

# Predict intestinal length
sp_list <- unique(int_moor_f$species)

pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total [:percent]",
  total = length(sp_list))

test_extrapolation <- function(){
  pb$tick(0)
  Sys.sleep(1 / 100)
  lapply(1:length(sp_list), function(x) {
    df <- filter(int_moor_f, species == sp_list[x])
    pred <- test_intestinal_length(data = df, ndraws = 2000)
    pb$tick(1)
    Sys.sleep(1 / 100)
    suppressMessages(left_join(df, pred))
  }) %>%
    bind_rows()
}

int_pred <- test_extrapolation()

rm(pb, test_extrapolation, sp_list)

plot_obs_vs_pred(data = int_pred,
                 pred = "int_length", obs = "il_log",
                 xlab = "Predicted ln int. length (mm)",
                 ylab = "Observed ln int. length (mm)",
                 point_size = 1.5,
                 text_size = 4)

# Predictions for unobserved species are quite accurate (R2=0.81, b=1, a=-0.16)

ggsave(here::here("outputs", "figures", "int_validation_species.png"),
       width = 10, height = 7, units = "cm", dpi = 600, type = "cairo")

# Check number of observations that fall within the 95% CI of the predictions
int_pred %>%
  mutate(check = ifelse(il_log >= .lower & il_log <= .upper, 1, 0)) %>%
  count(check) %>%
  mutate("%" = n/sum(n)*100)
# check     n   `%`
#     0   903  74.8
#     1   305  25.2
# ~25% of all predictions fall within the 95% CI
# Here, we are not incorporating the residual error in the prediction,
# which is instead done in "brms::predict.brmsfit()",
# thus our predictions are equivalent to "brms::fitted.brmsfit()"

#########################################################################################
#### Validate prediction for species identified at genus level

# In the carbonate dataset there are a couple of species identified at the genus level.
# Can we accurately predict intestinal length for them?
# The "intestinal_length_sp()" function predicts intestinal length at the genus-level
# by averaging predictions of all species with genetic data in the Fish Tree of Life belonging to the same genus.

# Validate these predictions on a subset of 200 fishes
# Sample 200 observations from the dataset used to train the model
# First, remove species without genetic data (n=3)
int_sp <- filter(int_moor_f, !species %in% check_name_fishtree(unique(int_moor_f$species), sampled = TRUE)$not_sampled)
# This step is not really necessary, but avoids potential NAs when there are no species with genetic data within a genus
# For example, for the genus Alloblennius (Blenniidae) there are:
fishtree::fishtree_taxonomy("Blenniidae")$Blenniidae$species  %>%
  stringr::str_subset(stringr::fixed("Alloblennius"))
# "Alloblennius anuchalis" "Alloblennius frondiculus" "Alloblennius jugularis" "Alloblennius parvus" "Alloblennius pictus"
# 5 species in the Fish Tree of Life, but
fishtree::fishtree_taxonomy("Blenniidae")$Blenniidae$sampled_species %>%
  stringr::str_subset(stringr::fixed("Alloblennius"))
# none have genetic data
# If we try to predict the genus-level intestinal length for Alloblennius sp. we will get NA
intestinal_length_sp(data = tibble(family = "Blenniidae",
                                   genus = "Alloblennius",
                                   species = "Alloblennius sp.",
                                   sl = 50,
                                   elongation = elongation_FB_sp(species)$elongation,
                                   trophic_level = trophic_level_FB_sp(species, type = "food items")$trophic_level),
                     ndraws = 100)
# family     genus           species              sl int_length .lower .upper .width
# Blenniidae Alloblennius    Alloblennius sp.     50         NA     NA     NA   0.95

# Now the sampling
set.seed(98765)
int_sp <- sample_n(int_sp, 200)
length(unique(int_sp$species)) # 100 different species
length(unique(int_sp$genus)) # 51 genera

# Change species names to "Genus sp."
int_sp <- mutate(int_sp, species = paste(genus, "sp."))

# Keep required columns plus observed values of intestinal length
int_sp <- select(int_sp, family, genus, species, sl, il_log)

# Get fish traits
# TROPHIC LEVEL
troph_int_sp <- trophic_level_FB_sp(unique(int_sp$species), type = "food items")
filter(troph_int_sp, is.na(trophic_level)) # no NAs

# ELONGATION
elon_int_sp <- elongation_FB_sp(unique(int_sp$species))
filter(elon_int_sp, is.na(elongation)) # no NAs

# Join traits to the dataset
int_sp <- left_join(int_sp, troph_int_sp[, -3]) %>% left_join(elon_int_sp[, -3])

rm(troph_int_sp, elon_int_sp)

# Predict intestinal length
int_sp_pred <- intestinal_length_sp(data = int_sp, ndraws = 2000)

# Plot observed vs predicted
left_join(int_sp, int_sp_pred) %>%
  plot_obs_vs_pred(pred = "int_length", obs = "il_log",
                   xlab = "Predicted ln int. length (mm)",
                   ylab = "Observed ln int. length (mm)",
                   point_size = 1.5,
                   text_size = 4)

# Genus-level predictions are quite accurate (R2=0.84, b=1.1, a=-0.44)

ggsave(here::here("outputs", "figures", "int_validation_genus.png"),
       width = 10, height = 7, units = "cm", dpi = 600, type = "cairo")

# We could also compare them with species-level predictions
set.seed(98765)
int_moor_f %>%
  filter(!species %in% check_name_fishtree(unique(int_moor_f$species), sampled = TRUE)$not_sampled) %>%
  sample_n(200) %>%
  predict(m_intestine, newdata = .) %>%
  cbind(left_join(int_sp, int_sp_pred)) %>%
  ggplot(aes(x = int_length, y = Estimate)) +
  geom_point(size = 1.5, alpha = 0.5, stroke = 0.1, color = "#0D0887FF") +
  geom_abline(linetype = 2, color = "red", size = 0.7) +
  geom_smooth(method = "lm", size = 0.7, alpha = 0.3, color = "#0D0887FF", fill = "#0D0887FF") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")), family = "serif", size = 3) +
  labs(x = "Genus-level predicted ln int. length (mm)",
       y = "Species-level predicted ln int. length (mm)") +
  theme(text = element_text(size = 10))

ggsave(here::here("outputs", "figures", "species_vs_genus_IL_pred.png"),
       width = 10, height = 8, units = "cm", dpi = 600, type = "cairo")

#########################################################################################
#### Predict intestinal length for fishes in the carbonate dataset
# Since intestinal length is related to body size, which will be included in the carbonate predictive model,
# we need to standardise intestinal length.
# To do this, we predict intestinal length at a common standard length for all species.
# We use 20 cm, but the value doesn't really matter as we don't have random slopes (results would be the same with a different size).

# Prepare dataset
data_caco3_il <- data_caco3_FB_traits %>%
  select(species, trophic_level, elongation) %>%
  unique() %>%
  mutate(sl = 200,
         id = species)

# Remove the two unidentified species
data_caco3_il <- filter(data_caco3_il, stringr::str_detect(species, stringr::fixed(" sp."), negate = TRUE))

# Check if species names match those in the Fish Tree of Life
check_name_fishtree(unique(data_caco3_il$species))
# These species names are incorrect or not present in the Fish Tree of Life:
# [1] "Pagrus auratus" "Ostorhinchus limenus"
# Family names of these species
data_caco3_FB_traits %>%
  filter(species == c("Pagrus auratus", "Ostorhinchus limenus")) %>%
  select(family, species) %>%
  unique()

fishtree::fishtree_taxonomy("Apogonidae")$Apogonidae[[1]]
# Changed "Apogon limenus" with "Ostorhinchus limenus" before to retrieve traits from FishBase,
# but "Ostorhinchus limenus" in the Fish Tree of Life is named "Apogon limenus"
fishtree::fishtree_taxonomy("Sparidae")$Sparidae[[1]]
# "Pagrus auratus" in the Fish Tree of Life is named "Chrysophrys auratus"
# Change these names
data_caco3_il <- data_caco3_il %>%
  mutate(species = recode(species,
                          "Ostorhinchus limenus" = "Apogon limenus",
                          "Pagrus auratus" = "Chrysophrys auratus"))

# Predict intestinal length
il_pred_caco3 <- intestinal_length(data = data_caco3_il, ndraw = 2000)

# Relative intestinal length
il_pred_caco3 <- il_pred_caco3 %>%
  mutate(ril = round(exp(int_length)/sl, 2))

# Change again the names to join the predictions to the dataset
il_pred_caco3 <- il_pred_caco3 %>%
  mutate(species = recode(species,
                          "Apogon limenus" = "Ostorhinchus limenus",
                          "Chrysophrys auratus" = "Pagrus auratus"))

# Predict genus-level intestinal length for the two unidentified species
data_caco3_il_sp <- data_caco3_FB_traits %>%
  filter(stringr::str_detect(species, stringr::fixed(" sp."))) %>%
  select(family, genus, species, trophic_level, elongation) %>%
  unique() %>%
  mutate(sl = 200)

il_pred_caco3_sp <- intestinal_length_sp(data = data_caco3_il_sp, ndraw = 2000)

# Relative intestinal length
il_pred_caco3_sp <- il_pred_caco3_sp %>%
  mutate(ril = round(exp(int_length)/sl, 2))

# Join the relative intestinal length to the dataset
data_caco3_traits <- left_join(data_caco3_FB_traits,
                               select(bind_rows(il_pred_caco3, il_pred_caco3_sp), species, ril))

# Save csv
#readr::write_csv(data_caco3_traits, here::here("data", "derived_data", "data_caco3_traits.csv"))

rm(data_caco3_il, data_caco3_il_sp, data_caco3_FB_traits, il_pred_caco3, il_pred_caco3_sp)

#########################################################################################

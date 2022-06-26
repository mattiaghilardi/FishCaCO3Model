#### Aims

# Retrieve species traits from Fishbase:
# trophic level, caudal fin aspect ratio, and body elongation

#########################################################################################
#### Load data
data_caco3 <- readr::read_csv(here::here("data", "raw_data", "caco3_excretion.csv"))

#########################################################################################
# Check if species names match those in FishBase
fishflux::name_errors(unique(data_caco3$species))
# Inaccurate species names found:
# "Apogon limenus" "Pelates sexlineatus" "Haemulon sp." "Mugil sp."

# The last two taxa are identified at the genus level
# Correct name for "Apogon limenus" and "Pelates sexlineatus"
rfishbase::validate_names(c("Apogon limenus", "Pelates sexlineatus"))
# "Ostorhinchus limenus" "Helotes sexlineatus"

# Change these names in the dataset
data_caco3 <- data_caco3 %>%
  mutate(species = dplyr::recode(species,
                                 "Apogon limenus" = "Ostorhinchus limenus",
                                 "Pelates sexlineatus" = "Helotes sexlineatus"),
         genus = dplyr::recode(genus,
                               "Apogon" = "Ostorhinchus",
                               "Pelates" = "Helotes"))

# Taxa identified at the species level
species <- unique(data_caco3$species[which(stringr::str_detect(data_caco3$species, stringr::fixed(" sp."), negate = TRUE))])
# Taxa identified at the genus level
species_sp <- unique(data_caco3$species[! (data_caco3$species %in% species)])

# TROPHIC LEVEL
troph_caco3 <- dplyr::bind_rows(trophic_level_FB(species, type = "food items"),
                                trophic_level_FB_sp(species_sp, type = "food items")) %>%
  dplyr::rename(level_troph = taxon_level)

dplyr::filter(troph_caco3, is.na(trophic_level)) # no NAs

# ASPECT RATIO
ar_caco3 <- dplyr::bind_rows(aspect_ratio_FB(species),
                             aspect_ratio_FB_sp(species_sp)) %>%
  dplyr::rename(level_ar = taxon_level)

dplyr::filter(ar_caco3, is.na(aspect_ratio)) # 1 NA - Gymnothorax javanicus

# Moray eels have no caudal fin, set aspect ratio to zero (following Villeger et al. 2010, Toussaint et al. 2016, Su et al. 2019)
ar_caco3$aspect_ratio <- tidyr::replace_na(ar_caco3$aspect_ratio, 0)

# ELONGATION
elon_caco3 <- dplyr::bind_rows(elongation_FB(species),
                               elongation_FB_sp(species_sp)) %>%
  dplyr::rename(level_elon = taxon_level)

dplyr::filter(elon_caco3, is.na(elongation)) # no NAs

# Join traits to the dataset
data_caco3_FB_traits <- dplyr::left_join(data_caco3, troph_caco3) %>%
  dplyr::left_join(ar_caco3) %>%
  dplyr::left_join(elon_caco3)

# Check if both traits are available at the family level for any species
nrow(dplyr::filter(data_caco3_FB_traits,
                   level_troph == "family" & level_ar == "family")) # none

rm(troph_caco3, ar_caco3, elon_caco3, data_caco3)

# Save csv
#readr::write_csv(data_caco3_FB_traits,
#                 here::here("data", "derived_data", "data_caco3_FB_traits.csv"))

#########################################################################################

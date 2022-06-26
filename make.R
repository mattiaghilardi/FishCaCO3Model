# Script to reproduce this project

# Install dependencies
renv::restore()

# Load packages and functions
devtools::load_all()

# Global options
options(FISHBASE_VERSION = "21.04") # set FishBase version
options(mc.cores = 4) # set number of cores

# Set plot theme
theme_set(theme_bw(base_size = 12, base_family = "serif") +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black")))

# Create output directory for models
dir.create(here::here("outputs", "models"))

# Load the scripts

# Analysis:

# 1 - FishBase traits
source(here::here("analysis", "1_FB_traits.R"))

# 2 - intestine extrapolation
source(here::here("analysis", "2_extrapolation_intestine_length.R"))

# 3 - CaCO3 excretion
source(here::here("analysis", "3_modelling_caco3_exc_rates.R"))

# 4 - CaCO3 composition
source(here::here("analysis", "4_modelling_caco3_composition.R"))

# 5 - sensitivity analysis
source(here::here("analysis", "5_sensitivity_analysis.R"))

# Figures:

# Fig. 1
source(here::here("analysis", "figure1.R"))

# Fig. 2
source(here::here("analysis", "figure2.R"))

# Fig. 3
source(here::here("analysis", "figure3.R"))

# Fig. 4
source(here::here("analysis", "figure4.R"))

# Fig. 5
source(here::here("analysis", "figure5.R"))

# Fig. 6
source(here::here("analysis", "figure6.R"))

# Supplementary figures
source(here::here("analysis", "supp_figures.R"))

# As they are in order can also simply run:
#scripts <- list.files(here::here("analysis"))
#lapply(scripts, function(x) source(here::here("analysis", x)))

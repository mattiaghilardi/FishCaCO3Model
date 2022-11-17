#' Test extrapolation of intestinal length
#'
#' @inheritParams intestinal_length
#' @inherit intestinal_length return
test_intestinal_length <- function(data, method = "fitted", ndraws = 1000, seed = 1,
                                   mc.cores = getOption("mc.cores", 1L),
                                   summary = TRUE, .point = mean, .width = 0.95) {

  # Check species names in the Fish Tree of Life
  # Not all species have the same name in FishBase and the Fish Tree of Life
  # For these species traits cannot be estimated
  suppressMessages(check_name_fishtree(unique(data$species)))

  # Complete phylogeny: species included in the model + new species
  sp   <- unique(c(as.character(int_moor_f$species), as.character(data$species)))
  tree <- fishtree::fishtree_complete_phylogeny(species = sp)

  # Estimate phylogenetic effect for new species
  unsampled <- gsub(" ", "_", unique(data$species))
  phylo_effect <- test_phylo_effect(ndraws = ndraws, seed = seed,
                                    mc.cores = mc.cores, phy = tree,
                                    unsampled = unsampled)

  # Extract model parameters
  model_param <- tidybayes::spread_draws(m_intestine,
                                         b_Intercept, b_scalesl_log, b_scaletrophic_level, b_scaleelon_log,
                                         sigma, nu, ndraws = ndraws, seed = seed)

  # Retain "ndraw" random draws
  model_param <- dplyr::select(model_param, c(-.chain, -.iteration))

  # Join phylogenetic effect to the other model parameters
  all_param <- suppressMessages(dplyr::left_join(phylo_effect, model_param))

  # Prepare dataset for trait's prediction
  all_param$species <- gsub("_", " ", all_param$species)
  newdata <- suppressMessages(dplyr::left_join(data, all_param))
  newdata$sl_log <- log(newdata$sl)
  newdata$elon_log <- log(newdata$elongation)

  # Predict trait for each draw
  sl_scale      <- stats::sd(int_moor_f$sl_log)
  sl_center     <- mean(int_moor_f$sl_log)
  troph_scale   <- stats::sd(int_moor_f$trophic_level)
  troph_center  <- mean(int_moor_f$trophic_level)
  elon_scale    <- stats::sd(int_moor_f$elon_log)
  elon_center   <- mean(int_moor_f$elon_log)

  newdata <- dplyr::mutate(newdata,
                           pred = b_Intercept + phy_eff +
                             (b_scalesl_log/sl_scale)*(sl_log - sl_center) +
                             (b_scaletrophic_level/troph_scale)*(trophic_level - troph_center) +
                             (b_scaleelon_log/elon_scale)*(elon_log - elon_center))

  if (method == "fitted") {
    newdata <- dplyr::select(newdata, c(id, species, sl, pred, .draw))
  } else {
    newdata <- dplyr::select(newdata, c(id, species, sl, pred, sigma, nu, .draw))
    newdata <- dplyr::group_by(newdata, id)
    newdata <- dplyr::mutate(newdata, pred = brms::rstudent_t(ndraws, nu, pred, sigma))
    newdata <- dplyr::select(newdata, c(id, species, sl, pred, .draw))
  }

  # If summary=TRUE return median or mean and the chosen CIs
  # If summary=FALSE return all draws
  if (summary) {
    newdata <- dplyr::group_by(newdata, id, species, sl)
    newdata <- tidybayes::point_interval(newdata, pred, .width = .width, .point = .point)
    newdata <- dplyr::ungroup(newdata)
    newdata <- dplyr::select(newdata, -c(.point, .interval))
  } else {
    newdata
  }

  colnames(newdata)[4] <- "int_length"

  newdata
}

#' Test estimation of phylogenetic effect for intestinal length
#' using observed or unobserved taxa
#'
#' @inheritParams phylo_effect
#' @param unsampled Vector of species (observed or unobserved) for
#'                  which the phylogenetic effect has to be predicted
#' @inherit phylo_effect return
test_phylo_effect <- function(ndraws = NULL, seed = NULL,
                              mc.cores = getOption("mc.cores", 1L), phy,
                              unsampled = NULL) {

  # Check tree
  if (!inherits(phy, "phylo") & !inherits(phy, "multiPhylo")) {
    stop("'phy' must be of class 'phylo' or 'multiPhylo'")
  }

  # Extract draws of phylogenetic effect
  phy_eff <- tidybayes::spread_draws(model = m_intestine, r_phylo[species,], ndraws = ndraws, seed = seed)
  phy_eff <- dplyr::filter(phy_eff, ! species %in% unsampled)
  phy_eff <- dplyr::select(phy_eff, c(-.chain, -.iteration))
  phy_eff <- .named_group_split(phy_eff, .draw, .keep = FALSE)

  # Predict traits for new species using phyEstimate() from "picante"
  # If 'phy' is of class 'multiPhylo' use one random tree for each draw
  phy_eff_pred <- parallel::mclapply(phy_eff, function(x){
    x <- tibble::column_to_rownames(x, colnames(x)[1])
    colnames(x)[1] <- "phy_eff"

    # If tree is of class "multiPhylo" sample one random tree
    if (inherits(phy, "phylo")) {
      tree <- phy
    } else {
      tree <- sample(phy, 1)[[1]]
    }

    # Predict trait
    trait_pred <- picante::phyEstimate(phy = tree, trait = x, method = "pic")
    colnames(trait_pred)[1] <- "phy_eff"
    trait_pred <- dplyr::select(trait_pred, -se)
    x <- rbind(x, trait_pred)
    x <- tibble::rownames_to_column(x, "species")
  }, mc.cores = mc.cores)

  # Convert list in data frame
  phy_eff_pred <- dplyr::bind_rows(phy_eff_pred, .id = ".draw")
  phy_eff_pred$.draw <- as.integer(phy_eff_pred$.draw)

  phy_eff_pred
}

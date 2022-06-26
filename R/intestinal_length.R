#' Estimate intestinal length of fishes
#'
#' Predict intestinal length of fish based on phylogeny, body morphology, and trophic level
#'
#' @param data A data frame with five columns:
#' \itemize{
#' \item{"id":} an identifier for each observation
#' \item{"species":} the taxonomic name in the format "Genus species"
#' \item{"sl":} standard length in "mm"
#' \item{"elongation":} body elongation (i.e. standard length/max body depth)
#' \item{"trophic_level":} trophic level (from 1 to 5, preferably based on food items)
#' }
#' @param method Either 'fitted' or 'predict'. Predictions performed with 'predict'
#'               have higher variance as they incorporate the residual error.
#'               However, the estimated means of both methods should be very similar
#' @param ndraws The number of draws to be used for the predictions
#' @param seed The seed for random number generation to make results reproducible.
#'             Must be different from NULL
#' @param mc.cores The number of cores to use for
#'                 \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#' @param summary Logical, if predictions should be summarised.
#'                If TRUE returns mean and sd, if FALSE returns all draws
#' @param .point Only if 'summary=TRUE', point summary function,
#'               e.g. 'mean' or 'median', default to 'mean'
#' @param .width Only if 'summary=TRUE', according to
#'               \code{\link[tidybayes:point_interval]{tidybayes::point_interval()}}
#'
#' @return A data frame with natural-log intestinal length in "mm"
#'
#' @importFrom fishtree fishtree_complete_phylogeny
#' @importFrom dplyr select left_join mutate group_by ungroup
#' @importFrom tidybayes spread_draws point_interval
#' @importFrom brms rstudent_t
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' dummy <- tibble(id = letters[1:4],
#'                 species = c("Gymnothorax javanicus",
#'                             "Lutjanus kasmira",
#'                             "Siganus doliatus",
#'                             "Sphyraena barracuda"),
#'                 sl = c(1200, 200, 150, 625))
#' d_elon <- elongation_FB(dummy$species)
#' d_troph <- trophic_level_FB(dummy$species, type = "food items")
#' dummy <- left_join(dummy, d_elon[, -3]) %>%
#'   left_join(d_troph[, -3])
#' intestinal_length(dummy, ndraw = 10)
#' intestinal_length(dummy, ndraw = 10, summary = FALSE)
#' intestinal_length(dummy, ndraw = 10, .point = median, .width = c(0.5, 0.95))
#' }
#'
#' @export
intestinal_length <- function(data, method = "fitted", ndraws = 1000, seed = 1,
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
  phylo_effect <- phylo_effect(model = m_intestine, r_phylo[species,], ndraws = ndraws,
                               seed = seed, mc.cores = mc.cores, phy = tree)

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

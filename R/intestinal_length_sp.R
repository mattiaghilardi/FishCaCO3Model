#' Estimate intestinal length for fishes identified to genus level
#'
#' Predict genus-level intestinal length of fish based on phylogeny, body morphology, and trophic level
#'
#' @param data A data frame with five columns:
#' \itemize{
#' \item{"family":} the taxonomic family
#' \item{"genus":} the taxonomic genus
#' \item{"species":} the taxonomic name in the form "Genus sp." or "Genus spp."
#' \item{"sl":} standard length in "mm"
#' \item{"elongation":} body elongation (i.e. standard length/max body depth)
#' \item{"trophic_level":} trophic level (from 1 to 5, preferably based on food items)
#' }
#' @param method Either 'fitted' or 'predict'. Predictions performed with 'predict'
#'               have higher variance as they incorporate the residual error.
#'               However, the estimated means of both methods should be very similar
#' @param ndraws The number of draws to be used for the predictions (defaults to 1000)
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
#' @importFrom fishtree fishtree_taxonomy
#' @importFrom dplyr left_join rename filter summarise_at group_by select ungroup bind_rows mutate
#' @importFrom tidybayes point_interval
#' @importFrom tidyr expand_grid
#'
#' @examples
#' \dontrun{
#' dummy <- data.frame(family = c("Siganidae", "Sparidae"),
#'                     genus = c("Siganus", "Sparus"),
#'                     species = c("Siganus sp.", "Sparus sp."),
#'                     sl = c(120, 160))
#' d_elon <- elongation_FB_sp(dummy$species)
#' d_troph <- trophic_level_FB_sp(dummy$species, type = "food items")
#' dummy <- left_join(dummy, d_elon[, -3]) %>%
#'   left_join(d_troph[, -3])
#' intestinal_length_sp(dummy, ndraws = 10)
#' intestinal_length_sp(dummy, ndraws = 10, summary = FALSE)
#' intestinal_length_sp(dummy, ndraws = 10, .point = median, .width = c(0.5, 0.95))
#' }
#'
#' @export
intestinal_length_sp <- function(data, method = "fitted", ndraws = 1000, seed = 1,
                                 mc.cores = getOption("mc.cores", 1L),
                                 summary = TRUE, .point = mean, .width = 0.95) {

  # Get all species with genetic data in the Fish Tree of Life for the genera in data
  sp <- fishtree::fishtree_taxonomy("Actinopteri")[[1]]$sampled_species
  sp <- data.frame(genus = gsub(" .*", "", sp), species = sp)

  # Merge to data
  df <- suppressMessages(dplyr::left_join(dplyr::rename(data, species_in = species), sp))

  # Filter NAs
  df_na <- dplyr::filter(df, is.na(species))
  df <- dplyr::filter(df, !is.na(species))

  if (nrow(df)) {
    # Add id for intestinal_length()
    df$id <- 1:nrow(df)

    # Predict intestinal length
    all_il <- intestinal_length(df, method = method, ndraws = ndraws, seed = seed,
                                mc.cores = mc.cores, summary = FALSE)

    # Add family and genus
    all_il <- suppressMessages(dplyr::left_join(dplyr::select(df, -elongation, -trophic_level), all_il))

    # Compute average for each draw
    mean_il <- dplyr::summarise_at(dplyr::group_by(all_il, family, genus, species_in, sl, .draw),
                                   "int_length", mean)
    mean_il <- dplyr::ungroup(mean_il)
  }

  if (exists("mean_il")) {
    # Add NAs
    if (nrow(df_na)) {
      df_na <- dplyr::select(df_na, family, genus, species_in, sl)
      df_na <- tidyr::expand_grid(df_na, .draw = unique(mean_il$.draw), int_length = NA)
      mean_il <- dplyr::bind_rows(mean_il, df_na)
    }
    mean_il
  } else {
      df_na <- dplyr::select(df_na, family, genus, species_in, sl)
      df_na <- dplyr::mutate(df_na, .draw = NA, int_length = NA)
      mean_il <- df_na
  }

  # Rename species
  mean_il <- dplyr::rename(mean_il, species = species_in)

  # If summary=TRUE return mean or median and the chosen CIs
  # If summary=FALSE return all draws
  if (summary) {
    mean_il <- dplyr::group_by(mean_il, family, genus, species, sl)
    mean_il <- tidybayes::point_interval(mean_il, int_length, .width = .width, .point = .point)
    mean_il <- dplyr::ungroup(mean_il)
    mean_il <- dplyr::select(mean_il, -c(.point, .interval))
  }

  mean_il
}

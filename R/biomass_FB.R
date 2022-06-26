#' Obtain fish biomass
#'
#' Compute fish biomass using a and b parameters from FishBase knowing total length
#'
#' @param data A data frame with two columns:
#' \itemize{
#' \item{"species":} taxonomic names in the form "Genus species"
#' \item{"length":} total length in "cm"
#' }
#' @param check_names Logical, check whether species names are correct
#'                    according to the accepted taxonomy in FishBase
#'
#' @return The input data frame with an additional column for biomass
#'
#' @importFrom rfishbase estimate
#' @importFrom dplyr mutate select left_join
#'
#' @examples
#' \dontrun{
#' biomass_FB(data.frame(species = c("Dentex dentex", "Chaetodon lunula"),
#'                       length = c(20, 12)),
#'            check_names = TRUE)
#' }
#'
#' @return
biomass_FB <- function(data, check_names = FALSE) {

  species_list <- unique(data$species)

  # Check species names
  if (check_names) {
    check_names_FB(species_list)
  }

  # Get a and b parameters for each species
  param <- rfishbase::estimate(species_list, fields = c("Species", "a", "b"))

  # Join parameters to the dataset and compute biomass
  dat <- dplyr::left_join(data, param, by = c("species" = "Species"))
  dat <- dplyr::mutate(dat, biomass = a*length^b)
  dat <- dplyr::select(dat, -a, -b)
  dat
}

#' Obtain fish biomass for taxa identified to genus level
#'
#' Compute fish biomass using genus-average a and b parameters from FishBase knowing total length
#'
#' @param data A data frame with two columns:
#' \itemize{
#' \item{"species":} taxonomic names in the form "Genus sp." or "Genus spp."
#' \item{"length":} total length in "cm"
#' }
#'
#' @return The input data frame with an additional column for biomass
#'
#' @importFrom rfishbase species_list estimate
#' @importFrom dplyr mutate summarise select left_join group_by
#'
#' @examples
#' \dontrun{
#' biomass_FB_sp(data.frame(species = c("Dentex sp.", "Chaetodon spp."),
#'                          length = c(20, 12)))
#' }
#'
#' @return
biomass_FB_sp <- function(data) {

  # Get all genera
  genus_list <- unique(gsub(" sp.*", "", unique(data$species)))

  # Get all valid species within these genera
  all_species <- rfishbase::species_list(Genus = genus_list)

  # Get a and b parameters for each species
  param <- rfishbase::estimate(all_species, fields = c("Species", "a", "b"))

  # Compute average for each genus
  param <- dplyr::mutate(param, genus = gsub(" .*", "", Species))
  param <- dplyr::summarise(dplyr::group_by(param, genus),
                            a = mean(a, na.rm = TRUE),
                            b = mean(b, na.rm = TRUE))
  # Join parameters to the dataset and compute biomass
  dat <- dplyr::mutate(data, genus = gsub(" .*", "", species))
  dat <- suppressMessages(dplyr::left_join(dat, param))
  dat <- dplyr::mutate(dat, biomass = a*length^b)
  dat <- dplyr::select(dat, -genus, -a, -b)
  dat
}

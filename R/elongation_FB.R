#' Obtain body elongation from FishBase
#'
#' Obtain average species elongation (i.e., standard length/maximum body depth)
#' from FishBase at the lowest available taxonomic level,
#' either species, genus or family level.
#'
#' @param species_list A vector of taxonomic names in the form "Genus species"
#' @param check_names Logical, check whether species names are correct
#'                    according to the accepted taxonomy in FishBase
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @return A tibble with species name, body elongation, and taxonomic level.
#'
#' @importFrom rfishbase load_taxa morphometrics
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup bind_rows tibble
#'
#' @examples
#' \dontrun{
#' elongation_FB("Sarpa salpa")
#' elongation_FB(c("Aulostomus chinensis", "Chaetodon lunulatus"))
#' library(rfishbase)
#' scombrids <- species_list(Family = "Scombridae")
#' elongation_FB(scombrids)
#' }
#'
#' @export
elongation_FB <- function(species_list, check_names = FALSE, mc.cores = getOption("mc.cores", 1L)) {

  # Check species names
  if (check_names) {
    check_names_FB(species_list)
  }

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]

  # All morphometrics
  morph <- rfishbase::morphometrics(taxo$Species, fields = c("Species", "SL", "BD"))
  morph$SL <- as.numeric(morph$SL)

  # Get elongation at the species, genus or family level for each species
  all_elon <- parallel::mclapply(species_list, function(x){

    # Get elongation for all species within the family
    fam <- taxo[taxo$Species == x, "Family"]
    fam_all <- taxo[taxo$Family == fam, "Species"]
    elon <- dplyr::filter(morph, Species %in% fam_all)
    elon <- dplyr::mutate(elon, elon = SL/BD)

    # Join taxonomic info to the aspect ratios
    elon <- suppressMessages(dplyr::left_join(elon, taxo[, c("Family", "Genus", "Species")]))

    # Compute species, genus and family level aspect ratio
    elon <- dplyr::group_by(elon, Species)
    elon <- dplyr::mutate(elon, elon_sp = mean(elon, na.rm = TRUE))
    elon <- dplyr::ungroup(elon)
    elon <- dplyr::group_by(elon, Genus)
    elon <- dplyr::mutate(elon, elon_gen = mean(elon_sp, na.rm = TRUE))
    elon <- dplyr::ungroup(elon)
    elon$elon_fam <- mean(elon$elon_sp, na.rm = TRUE)

    # Retain only the desired species
    elon <- elon[elon$Species == x,]

    # Retain only the first row
    elon <- elon[1,]

    # Retain the elongation at the lowest possible taxonomic level
    # Add taxonomic level of elongation's value
    elon <- if (!is.na(elon$elon_sp)) {
      dplyr::mutate(elon,
                    elongation = elon_sp,
                    level = "species")
    } else if (is.na(elon$elon_sp) & !is.na(elon$elon_gen)) {
      dplyr::mutate(elon,
                    elongation = elon_gen,
                    level = "genus")
    } else {
      dplyr::mutate(elon,
                    elongation = elon_fam,
                    level = "family")
    }

    # Round elongation
    elon$elongation <- round(elon$elongation, 2)

    # Tibble for each species
    dplyr::tibble(species = x,
                  elongation = elon$elongation,
                  taxon_level = elon$level)

  }, mc.cores = mc.cores)

  # List to tibble
  all_elon <- dplyr::bind_rows(all_elon)
  all_elon
}

#' Obtain body elongation from FishBase for taxa identified to genus level
#'
#' Obtain average species elongation (i.e., standard length/maximum body depth)
#' from FishBase at the lowest available taxonomic level,
#' either genus or family level.
#'
#' @param species_list A vector of taxonomic names in the form "Genus sp." or "Genus spp."
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @return A tibble with the provided name, body elongation, and taxonomic level.
#'
#' @importFrom rfishbase load_taxa morphometrics
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup bind_rows tibble
#'
#' @examples
#' \dontrun{
#' elongation_FB_sp("Sarpa sp.")
#' elongation_FB_sp(c("Epinephelus sp.", "Sillago spp."))
#' }
#'
#' @export
elongation_FB_sp <- function(species_list, mc.cores = getOption("mc.cores", 1L)) {

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]

  # All morphometrics
  morph <- rfishbase::morphometrics(taxo$Species, fields = c("Species", "SL", "BD"))
  morph$SL <- as.numeric(morph$SL)

  # Get elongation at the genus or family level for each species
  all_elon <- parallel::mclapply(species_list, function(x){

    # Genus
    genus <- gsub(" sp.*", "", x)

    # Get aspect ratio for all species within the family
    fam <- taxo[taxo$Genus == genus, "Family"][1]
    fam_all <- taxo[taxo$Family == fam, "Species"]
    elon <- dplyr::filter(morph, Species %in% fam_all)
    elon <- dplyr::mutate(elon, elon = SL/BD)

    # Join taxonomic info to the aspect ratios
    elon <- suppressMessages(dplyr::left_join(elon, taxo[, c("Family", "Genus", "Species")]))

    # Compute species, genus and family level aspect ratio
    elon <- dplyr::group_by(elon, Species)
    elon <- dplyr::mutate(elon, elon_sp = mean(elon, na.rm = TRUE))
    elon <- dplyr::ungroup(elon)
    elon <- dplyr::group_by(elon, Genus)
    elon <- dplyr::mutate(elon, elon_gen = mean(elon_sp, na.rm = TRUE))
    elon <- dplyr::ungroup(elon)
    elon$elon_fam <- mean(elon$elon_sp, na.rm = TRUE)

    # Retain only the desired species
    elon <- elon[elon$Genus == genus,]

    # Retain only the first row
    elon <- elon[1,]

    # Retain the elongation at the lowest possible taxonomic level
    # Add taxonomic level of elongation's value
    elon <- if (!is.na(elon$elon_gen)) {
      dplyr::mutate(elon,
                    elongation = elon_gen,
                    level = "genus")
    } else {
      dplyr::mutate(elon,
                    elongation = elon_fam,
                    level = "family")
    }

    # Round elongation
    elon$elongation <- round(elon$elongation, 2)

    # Tibble for each species
    dplyr::tibble(species = x,
                  elongation = elon$elongation,
                  taxon_level = elon$level)

  }, mc.cores = mc.cores)

  # List to tibble
  all_elon <- dplyr::bind_rows(all_elon)
  all_elon
}

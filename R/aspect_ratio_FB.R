#' Obtain aspect ratio of caudal fin from FishBase
#'
#' Obtain average species aspect ratio
#' from FishBase at the lowest available taxonomic level,
#' either species, genus or family level.
#'
#' @param species_list A vector of taxonomic names in the form "Genus species"
#' @param check_names Logical, check whether species names are correct
#'                    according to the accepted taxonomy in FishBase
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @return A tibble with species name, aspect ratio, and taxonomic level.
#'
#' @importFrom rfishbase load_taxa morphometrics
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup bind_rows tibble
#'
#' @examples
#' \dontrun{
#' aspect_ratio_FB("Caranx ignobilis")
#' aspect_ratio_FB(c("Chaetodon auriga", "Caesio cuning"))
#' library(rfishbase)
#' scombrids <- species_list(Family = "Scombridae")
#' aspect_ratio_FB(scombrids)
#' }
#'
#' @export
aspect_ratio_FB <- function(species_list, check_names = FALSE, mc.cores = getOption("mc.cores", 1L)) {

  # Check species names
  if (check_names) {
    check_names_FB(species_list)
  }

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]

  # All morphometrics
  morph <- rfishbase::morphometrics(taxo$Species, fields = c("Species", "AspectRatio"))

  # Get aspect ratio at the species, genus or family level for each species
  all_ar <- parallel::mclapply(species_list, function(x){

    # Get aspect ratio for all species within the family
    fam <- taxo[taxo$Species == x, "Family"]
    fam_all <- taxo[taxo$Family == fam, "Species"]
    ar  <- dplyr::filter(morph, Species %in% fam_all)

    # Join taxonomic info to the aspect ratios
    ar  <- suppressMessages(dplyr::left_join(ar, taxo[, c("Family", "Genus", "Species")]))

    # Compute species, genus and family level aspect ratio
    ar <- dplyr::group_by(ar, Species)
    ar <- dplyr::mutate(ar, ar_sp = mean(AspectRatio, na.rm = TRUE))
    ar <- dplyr::ungroup(ar)
    ar <- dplyr::group_by(ar, Genus)
    ar <- dplyr::mutate(ar, ar_gen = mean(ar_sp, na.rm = TRUE))
    ar <- dplyr::ungroup(ar)
    ar$ar_fam <- mean(ar$ar_sp, na.rm = TRUE)

    # Retain only the desired species
    ar <- ar[ar$Species == x,]

    # Retain only the first row
    ar <- ar[1,]

    # Retain the aspect ratio at the lowest possible taxonomic level
    # Add taxonomic level of aspect ratio's value
    ar <- if (!is.na(ar$ar_sp)) {
      dplyr::mutate(ar,
                    ar = ar_sp,
                    level = "species")
    } else if (is.na(ar$ar_sp) & !is.na(ar$ar_gen)) {
      dplyr::mutate(ar,
                    ar = ar_gen,
                    level = "genus")
    } else {
      dplyr::mutate(ar,
                    ar = ar_fam,
                    level = "family")
    }

    # Round aspect ratio
    ar$ar <- round(ar$ar, 2)

    # Tibble for each species
    dplyr::tibble(species = x,
                  aspect_ratio = ar$ar,
                  taxon_level = ar$level)

  }, mc.cores = mc.cores)

  # List to tibble
  all_ar <- dplyr::bind_rows(all_ar)
  all_ar
}

#' Obtain aspect ratio of caudal fin from FishBase for taxa identified to genus level
#'
#' Obtain average species aspect ratio from FishBase at the lowest
#' available taxonomic level, either genus or family level.
#'
#' @param species_list A vector of taxonomic names in the form "Genus sp." or "Genus spp."
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @return A tibble with the provided name, aspect ratio, and taxonomic level.
#'
#' @importFrom rfishbase load_taxa morphometrics
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup bind_rows tibble
#'
#' @examples
#' \dontrun{
#' aspect_ratio_FB_sp("Caranx sp.")
#' aspect_ratio_FB_sp(c("Acanthurus sp.", "Lutjanus spp."))
#' }
#'
#' @export
aspect_ratio_FB_sp <- function(species_list, mc.cores = getOption("mc.cores", 1L)) {

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]

  # All morphometrics
  morph <- rfishbase::morphometrics(taxo$Species, fields = c("Species", "AspectRatio"))

  # Get aspect ratio at the genus or family level for each species
  all_ar <- parallel::mclapply(species_list, function(x){

    # Genus
    genus <- gsub(" sp.*", "", x)

    # Get aspect ratio for all species within the family
    fam <- taxo[taxo$Genus == genus, "Family"][1]
    fam_all <- taxo[taxo$Family == fam, "Species"]
    ar  <- dplyr::filter(morph, Species %in% fam_all)

    # Join taxonomic info to the aspect ratios
    ar  <- suppressMessages(dplyr::left_join(ar, taxo[, c("Family", "Genus", "Species")]))

    # Compute species, genus and family level aspect ratio
    ar <- dplyr::group_by(ar, Species)
    ar <- dplyr::mutate(ar, ar_sp = mean(AspectRatio, na.rm = TRUE))
    ar <- dplyr::ungroup(ar)
    ar <- dplyr::group_by(ar, Genus)
    ar <- dplyr::mutate(ar, ar_gen = mean(ar_sp, na.rm = TRUE))
    ar <- dplyr::ungroup(ar)
    ar$ar_fam <- mean(ar$ar_sp, na.rm = TRUE)

    # Retain only the desired species
    ar <- ar[ar$Genus == genus,]

    # Retain only the first row
    ar <- ar[1,]

    # Retain the aspect ratio at the lowest possible taxonomic level
    # Add taxonomic level of aspect ratio's value
    ar <- if (!is.na(ar$ar_gen)) {
      dplyr::mutate(ar,
                    ar = ar_gen,
                    level = "genus")
    } else {
      dplyr::mutate(ar,
                    ar = ar_fam,
                    level = "family")
    }

    # Round aspect ratio
    ar$ar <- round(ar$ar, 2)

    # Tibble for each species
    dplyr::tibble(species = x,
                  aspect_ratio = ar$ar,
                  taxon_level = ar$level)

  }, mc.cores = mc.cores)

  # List to tibble
  all_ar <- dplyr::bind_rows(all_ar)
  all_ar

}

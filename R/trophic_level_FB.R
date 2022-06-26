#' Obtain trophic level from FishBase
#'
#' Obtain trophic level from FishBase at the lowest available taxonomic level,
#' either species, genus or family level.
#'
#' @param species_list A vector of taxonomic names in the form "Genus species"
#' @param type Character, indicating the type of estimate desired. Possible values are:
#'             "both" (default), "diet", "food items" or any abbreviation of these (see Details)
#' @param check_names Logical, check whether species names are correct
#'                    according to the accepted taxonomy in FishBase
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @details The function can retrieve three different trophic level estimates:
#' \itemize{
#' \item{\code{type="diet"}:} retrieves trophic level based on diet composition only. These estimates are available only for some species. Note that for some families there are no diet studies available;
#' \item{\code{type="food items"}:} retrieves trophic level based on food items only;
#' \item{\code{type="both"}:} computes the average of the two when both are available, otherwise return the available estimate.
#' }
#'
#' @return A tibble with species name, trophic level, and taxonomic level.
#'
#' @importFrom rfishbase load_taxa ecology
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup bind_rows tibble
#'
#' @examples
#' \dontrun{
#' trophic_level_FB("Cephalopholis argus", type = "diet")
#' trophic_level_FB(c("Acanthurus achilles", "Epibulus insidiator"), type = "food items")
#' library(rfishbase)
#' scombrids <- species_list(Family = "Scombridae")
#' trophic_level_FB(scombrids, type = "both")
#' }
#'
#' @export
trophic_level_FB <- function(species_list, type = c("both", "diet", "food items"),
                             check_names = FALSE, mc.cores = getOption("mc.cores", 1L)) {

  # Check type
  type <- match.arg(type)

  # Check species names
  if (check_names) {
    check_names_FB(species_list)
  }

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]

  # All ecology
  ecol <- rfishbase::ecology(taxo$Species, fields = c("Species", "DietTroph", "DietSeTroph", "FoodTroph", "FoodSeTroph"))
  # For a few species the FoodTroph values are shifted by one cell, correct them
  ecol <- dplyr::mutate(ecol, FoodTroph = ifelse(FoodTroph > 5, FoodSeTroph, FoodTroph))

  # Get trophic level at the species, genus or family level for each species
  all_troph <- parallel::mclapply(species_list, function(x){

    # Get trophic level for all species within the family
    fam <- taxo[taxo$Species == x, "Family"]
    fam_all <- taxo[taxo$Family == fam, "Species"]
    tl <- dplyr::filter(ecol, Species %in% fam_all)

    # Join taxonomic info to the trophic levels
    tl <- suppressMessages(dplyr::left_join(tl, taxo[, c("Family", "Genus", "Species")]))

    # Compute species, genus and family level trophic lavel
    tl$troph_sp <- if (type == "diet") {
      tl$DietTroph
    } else if (type == "food items") {
      tl$FoodTroph
    } else {
      rowMeans(tl[, c("DietTroph", "FoodTroph")], na.rm = TRUE)
    }
    tl <- dplyr::group_by(tl, Genus)
    tl <- dplyr::mutate(tl, troph_gen = mean(troph_sp, na.rm = TRUE))
    tl <- dplyr::ungroup(tl)
    tl$troph_fam <- mean(tl$troph_sp, na.rm = TRUE)

    # Retain only the desired species
    tl <- tl[tl$Species == x,]

    # Some species have different entries for different stocks of the same species
    # Retain only the first
    # Data is either identical to the first or simply missing in the additional stocks
    tl <- tl[1,]

    # Retain the trophic level at the lowest possible taxonomic level
    # Add taxonomic level of trophic level's value
    tl <- if (!is.na(tl$troph_sp)) {
      dplyr::mutate(tl,
                    troph = troph_sp,
                    level = "species")
    } else if (is.na(tl$troph_sp) & !is.na(tl$troph_gen)) {
      dplyr::mutate(tl,
                    troph = troph_gen,
                    level = "genus")
    } else {
      dplyr::mutate(tl,
                    troph = troph_fam,
                    level = "family")
    }

    # Round trophic level
    tl$troph <- round(tl$troph, 2)

    # Tibble for each species
    dplyr::tibble(species = x,
                  trophic_level = tl$troph,
                  taxon_level = tl$level)

  }, mc.cores = mc.cores)

  # List to tibble
  all_troph <- dplyr::bind_rows(all_troph)
  all_troph
}

#' Obtain trophic level from FishBase for taxa identified to genus level
#'
#' Obtain trophic level from FishBase at the lowest available taxonomic level,
#' either genus or family level.
#'
#' @param species_list A vector of taxonomic names in the form "Genus sp." or "Genus spp."
#' @param type Character, indicating the type of estimate desired. Possible values are:
#'             "both" (default), "diet", "food items" or any abbreviation of these (see Details)
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @details The function can retrieve three different trophic level estimates:
#'   \code{type="diet"} retrieves trophic level based on diet composition only. These estimates are available only for some species. Note that for some families there are no diet studies available;
#'   \code{type="food items"} retrieves trophic level based on food items only;
#'   \code{type="both"} computes the average of the two when both are available, otherwise return the available estimate.
#'
#' @return A tibble with species name, trophic level, and taxonomic level.
#'
#' @importFrom rfishbase load_taxa ecology
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup bind_rows tibble
#'
#' @examples
#' \dontrun{
#' trophic_level_FB_sp("Cephalopholis sp.", type = "diet")
#' trophic_level_FB_sp(c("Acanthurus sp.", "Lethrinus spp."), type = "food items")
#' }
#'
#' @export
trophic_level_FB_sp <- function(species_list, type = c("both", "diet", "food items"),
                                mc.cores = getOption("mc.cores", 1L)) {

  # Check type
  type <- match.arg(type)

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]

  # All ecology
  ecol <- rfishbase::ecology(taxo$Species, fields = c("Species", "DietTroph", "DietSeTroph", "FoodTroph", "FoodSeTroph"))
  # For a few species the FoodTroph values are shifted by one cell, correct them
  ecol <- dplyr::mutate(ecol, FoodTroph = ifelse(FoodTroph > 5, FoodSeTroph, FoodTroph))

  # Get trophic level at the genus or family level for each genus
  all_troph <- parallel::mclapply(species_list, function(x){

    # Genus
    genus <- gsub(" sp.*", "", x)

    # Get trophic level for all species within the family
    fam <- taxo[taxo$Genus == genus, "Family"][1]
    fam_all <- taxo[taxo$Family == fam, "Species"]
    tl <- dplyr::filter(ecol, Species %in% fam_all)

    # Join taxonomic info to the trophic levels
    tl <- suppressMessages(dplyr::left_join(tl, taxo[, c("Family", "Genus", "Species")]))

    # Compute species, genus and family level trophic lavel
    tl$troph_sp <- if (type == "diet") {
      tl$DietTroph
    } else if (type == "food items") {
      tl$FoodTroph
    } else {
      rowMeans(tl[, c("DietTroph", "FoodTroph")], na.rm = TRUE)
    }

    tl <- dplyr::group_by(tl, Genus)
    tl <- dplyr::mutate(tl, troph_gen = mean(troph_sp, na.rm = TRUE))
    tl <- dplyr::ungroup(tl)
    tl$troph_fam <- mean(tl$troph_sp, na.rm = TRUE)

    # Retain only the desired genus
    tl <- tl[tl$Genus == genus,]

    # Retain only the first row
    tl <- tl[1,]

    # Retain the trophic level at the lowest possible taxonomic level
    # Add taxonomic level of trophic level's value
    tl <- if (!is.na(tl$troph_gen)) {
      dplyr::mutate(tl,
                    troph = troph_gen,
                    level = "genus")
    } else {
      dplyr::mutate(tl,
                    troph = troph_fam,
                    level = "family")
    }

    # Round trophic level
    tl$troph <- round(tl$troph, 2)

    # Tibble for each species
    dplyr::tibble(species = x,
                  trophic_level = tl$troph,
                  taxon_level = tl$level)

  }, mc.cores = mc.cores)

  # List to tibble
  all_troph <- dplyr::bind_rows(all_troph)
  all_troph
}

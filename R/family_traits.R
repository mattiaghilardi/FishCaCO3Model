#' Obtain family-average trophic level, caudal fin aspect ratio, body elongation, and relative intestinal length
#'
#' @param family_list A vector of taxonomic families
#' @param type Character, indicating the type of trophic level desired. Possible values are:
#'             "both" (default), "diet", "food items" or any abbreviation of these (see Details)
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#' @param ... Additional arguments passed to intestinal_length()
#'
#' @details The function can retrieve three different trophic level estimates:
#' \itemize{
#' \item{\code{type="both"}:} computes the average of the trophic levels based on
#'                            diet and food items when both are available,
#'                            otherwise return the available estimate;
#' \item{\code{type="diet"}:} retrieves trophic level based on diet composition only.
#'                            These estimates are available only for some species.
#'                            Note that for some families there are no diet studies available;
#' \item{\code{type="food items"}:} retrieves trophic level based on food items only.
#' }
#'
#' @return A tibble with species name, trophic level, caudal fin aspect ratio,
#'         body elongation and relative intestinal length.
#'
#' @importFrom rfishbase load_taxa morphometrics ecology
#' @importFrom parallel mclapply
#' @importFrom dplyr filter left_join group_by mutate ungroup select rename
#' @importFrom fishtree fishtree_taxonomy
#'
family_traits <- function(family_list, type = c("both", "diet", "food items"),
                          mc.cores = getOption("mc.cores", 1L), ...) {

  # Check type
  type <- match.arg(type)

  # Retrieve taxonomic info
  taxo <- rfishbase::load_taxa()
  taxo <- as.data.frame(taxo)
  taxo <- taxo[!is.na(taxo$Species),]
  taxo <- dplyr::filter(taxo, Family %in% family_list)

  # Ecology
  ecol <- rfishbase::ecology(taxo$Species, fields = c("Species", "DietTroph", "DietSeTroph", "FoodTroph", "FoodSeTroph"))
  # For a few species the FoodTroph values are shifted by one cell, correct them
  ecol <- dplyr::mutate(ecol, FoodTroph = ifelse(FoodTroph > 5, FoodSeTroph, FoodTroph))

  # Morphometrics
  morph <- rfishbase::morphometrics(taxo$Species, fields = c("Species", "SL", "BD", "AspectRatio"))
  morph$SL <- as.numeric(morph$SL)
  morph <- dplyr::mutate(morph, elon = SL/BD)

  # Compute species and family level trophic lavel
  ecol$troph_sp <- if (type == "diet") {
    ecol$DietTroph
  } else if (type == "food items") {
    ecol$FoodTroph
  } else {
    rowMeans(ecol[, c("DietTroph", "FoodTroph")], na.rm = TRUE)
  }
  ecol <- suppressMessages(dplyr::left_join(ecol, taxo[, c("Family", "Genus", "Species")]))
  ecol <- dplyr::group_by(ecol, Family)
  ecol <- dplyr::mutate(ecol,
                        troph_fam = round(mean(troph_sp, na.rm = TRUE), 2))
  ecol <- dplyr::ungroup(ecol)
  ecol_fam <- unique(dplyr::select(ecol, Family, troph_fam))

  # Compute species and family level aspect ratio and elongation
  morph <- dplyr::group_by(morph, Species)
  morph <- dplyr::mutate(morph,
                         ar_sp = mean(AspectRatio, na.rm = TRUE),
                         elon_sp = mean(elon, na.rm = TRUE))
  morph <- dplyr::ungroup(morph)
  morph <- unique(dplyr::select(morph, Species, ar_sp, elon_sp))
  morph <- suppressMessages(dplyr::left_join(morph, taxo[, c("Family", "Genus", "Species")]))
  morph <- dplyr::group_by(morph, Family)
  morph <- dplyr::mutate(morph,
                         ar_fam = round(mean(ar_sp, na.rm = TRUE), 2),
                         elon_fam = round(mean(elon_sp, na.rm = TRUE), 2)
                         )
  morph <- dplyr::ungroup(morph)
  morph_fam <- unique(dplyr::select(morph, Family, ar_fam, elon_fam))

  # Predict species level RIL and compute family average
  sp <- fishtree::fishtree_taxonomy(ranks = "Actinopteri")[[1]]$sampled_species
  ecol <- dplyr::filter(ecol, Species %in% sp)
  ecol <- unique(dplyr::select(ecol, Family, Species, troph_fam))
  morph <- dplyr::filter(morph, Species %in% sp)
  morph <- unique(dplyr::select(morph, Family, Species, elon_fam))
  df <- suppressMessages(dplyr::left_join(ecol, morph))
  df <- dplyr::rename(df, species = Species, trophic_level = troph_fam, elongation = elon_fam)
  ril <- intestinal_length(data = dplyr::mutate(df,
                                                id = 1:nrow(df),
                                                sl = 200),
                           mc.cores = mc.cores,
                           ...)
  ril <- dplyr::mutate(ril,
                       int_length = exp(int_length),
                       ril = int_length/sl)
  ril <- dplyr::left_join(ril, dplyr::select(df, Family, species))
  ril <- dplyr::group_by(ril, Family)
  ril <- dplyr::mutate(ril,
                       ril_fam = round(mean(ril, na.rm = TRUE), 2))
  ril <- dplyr::ungroup(ril)
  ril_fam <- unique(dplyr::select(ril, Family, ril_fam))

  # Output
  traits <- dplyr::left_join(ecol_fam, morph_fam)
  traits <- dplyr::left_join(traits, ril_fam)
  traits <- dplyr::rename(traits,
                          family = Family,
                          trophic_level = troph_fam,
                          aspect_ratio = ar_fam,
                          elongation = elon_fam,
                          ril = ril_fam)
  traits
}

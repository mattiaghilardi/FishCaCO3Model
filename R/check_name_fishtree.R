#' Check species names and sampled species in the Fish Tree of Life
#'
#' @param species_list A vector of taxonomic names in the form "Genus species"
#' @param sampled Logical, if the list has to be checked for unsampled species
#'                (i.e. those without genetic data, see Rabosky et al. 2018 for detailed methodology)
#'
#' @return If "sampled=FALSE", a message and a vector of wrong species names (if any).
#'         If "sampled=TRUE", a message and a list with two components:
#' \itemize{
#' \item{"name_error":} wrong species names
#' \item{"not_sampled":} species without genetic data
#' }
#'
#' @importFrom fishtree fishtree_taxonomy
#'
#' @examples
#' \dontrun{
#' check_name_fishtree("Boops boops")
#' check_name_fishtree(c("Balistapus undulatus", "Plectropomus leopardus"),  sampled = TRUE)
#' }
#'
#' @references Rabosky D. L. et al. (2018) An inverse latitudinal gradient in
#'             speciation rate for marine fishes. Nature, 559, 392–395.
#'             https://doi.org/10.1038/s41586-018-0273-1
#'
#' @export
check_name_fishtree <- function(species_list, sampled = FALSE) {

  # Retrieve all species in the Fish tree of Life
  sp_fishtree <- fishtree::fishtree_taxonomy("Actinopteri")[[1]]

  # Check names
  sp_error <- species_list[!(species_list %in% sp_fishtree$species)]

  # If sampled=TRUE check sampled species
  if (isTRUE(sampled)) {
    sp_not_sampled <- species_list[!(species_list %in% sp_fishtree$sampled_species)]

    # Message
    if (length(sp_error) == 0) {
      message("All species names are correct")
    } else {
      message("Some species name are incorrect or not present in the Fish Tree of Life")
    }
    if (length(sp_not_sampled) == 0) {
      message("All species have genetic data")
    } else {
      message("Some species do not have genetic data")
    }

    # Output
    list(name_error = sp_error,
         not_sampled = sp_not_sampled)

  } else {
    if (length(sp_error) == 0) {
      message("All species names are correct")
    } else {
      message("These species names are incorrect or not present in the Fish Tree of Life:")
      sp_error
    }
  }
}

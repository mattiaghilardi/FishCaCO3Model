#' A modified version of "fishflux::check_name_fishbase()" that returns wrong names
#' 
#' @param species_list A vector of taxonomic names in the form "Genus species"
#' 
#' @importFrom fishflux name_errors
#' 
check_names_FB <- function(species_list) {
  sp_error <- suppressMessages(fishflux::name_errors(species_list))
  if (length(sp_error) > 0) {
    stop("Species name is incorrect: ", paste(sp_error, collapse = ", "))
  }
}


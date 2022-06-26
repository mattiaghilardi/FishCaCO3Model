#' Estimate phylogenetic effect for unobserved taxa from a Bayesian phylogenetic model
#'
#' \code{phylo_effect()} uses \code{\link[picante:phyEstimate]{picante::phyEstimate()}}
#' to estimate the phylogenetic effect for unobserved taxa through phylogenetic
#' ancestral state reconstruction. The estimation is repeated \code{ndraws} times to provide
#' a distribution comparable to that provided by the model for observed taxa.
#' The function is available for one term only (i.e. the intercept).
#'
#' @param model A Bayesian phylogenetic model of class 'brmsfit'
#' @param ... Expression in the form of `variable_name[dimension_1, ]`. (see \code{\link[tidybayes:spread_draws]{tidybayes::spread_draws()}})
#' @param ndraws The number of draws to return, or NULL to return all draws
#' @param seed A seed to use when subsampling draws (i.e. when ndraw is not NULL)
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#' @param phy An object of class "phylo" or "multiPhylo" including the new taxa
#'            for which the phylogenetic effect has to be predicted
#'
#' @return A data frame
#'
#' @importFrom picante phyEstimate
#' @importFrom tidybayes spread_draws
#' @importFrom dplyr select bind_rows
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' # Imagine you have a data frame called 'data' with these columns:
#' # 'y' - the response variable
#' # 'x1' - predictor 1
#' # 'x2' - predictor 2
#' # 'phylo' - taxonomic names
#' # Then you have a phylogenetic tree (class "phylo" or "multiphylo") for these species called 'tree1'
#' # You can create a phylogenetic correlation matrix using ape::vcv.phylo(..., corr = TRUE)
#' A <- ape::vcv.phylo(tree1, corr = TRUE)
#' # Then fit a phylogenetic model in brms
#' m1 <- brm(y ~ x1 + x2 + (1|gr(phylo, cov = A)),
#'           data = data, family = gaussian(), data2 = list(A = A), ...)
#' # With this model you can make predictions for the observed taxa
#' # (i.e., those used to train the model), but if you would like to make
#' # predictions for unobserved taxa you can only use the parameter 'allow_new_levels=TRUE'.
#' # However, this includes the uncertainty associated with the random intercept
#' # in the prediction, but does not account for the phylogenetic structure.
#' # To account for the phylogeny we could predict the phylogenetic effect for
#' # the unobserved species using 'phylo_effect()'.
#' # You need a second phylogenetic tree called 'tree2', which is 'tree1'
#' # with the additional taxa for which you would like to predict y and for which
#' # x1 and x2 are available.
#' phylo_eff_pred <- phylo_effect(model = m1, r_phylo[species,],
#'                                ndraws = 1000, seed = 123, phy = tree2)
#' Now "y" can be predicted for these new species
#' }
#'
#' @export
phylo_effect <- function(model, ..., ndraws = NULL, seed = NULL,
                         mc.cores = getOption("mc.cores", 1L), phy) {

  # Check tree
  if (!inherits(phy, "phylo") & !inherits(phy, "multiPhylo")) {
    stop("'phy' must be of class 'phylo' or 'multiPhylo'")
  }

  # Extract draws of phylogenetic effect
  phy_eff <- tidybayes::spread_draws(model = model, ... = ..., ndraws = ndraws, seed = seed)
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

# Named version of dplyr::group_split()
# List elements are identified by group's names
# From: https://github.com/tidyverse/dplyr/issues/4223
#
#' @importFrom dplyr group_by group_keys group_split
#' @importFrom rlang eval_bare expr set_names
#'
.named_group_split <- function(.tbl, ..., .keep = FALSE) {
  grouped <- dplyr::group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!dplyr::group_keys(grouped), sep = " / ")))

  list <- dplyr::group_split(grouped, .keep = .keep)
  rlang::set_names(list, names)
}

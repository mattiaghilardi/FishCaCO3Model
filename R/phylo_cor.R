#' Phylogenetic correlation matrix for class 'phylo' or 'multiPhylo'
#'
#' @param phy An object of class 'phylo' or 'multiPhylo'
#' @param mc.cores The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#'                 Must be exactly 1 on Windows (which uses the master process)
#'
#' @return For class 'phylo', a correlation matrix; for class 'multiphylo',
#'         a list with the mean and standard deviation correlation matrices
#'
#' @importFrom ape vcv.phylo
#' @importFrom parallel mclapply
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' library(fishtree)
#' # Phylo
#' tree <- fishtree_phylogeny(rank = "Balistidae")
#' phylo_cor(tree)
#' # Multiphylo
#' trees <- fishtree_complete_phylogeny(rank = "Acanthuridae")
#' phylo_cor(trees)
#' }
#'
#' @export
phylo_cor <- function(phy, mc.cores = getOption("mc.cores", 1L)) {

  # Check tree
  if (!inherits(phy, "phylo") & !inherits(phy, "multiPhylo")) {
    stop("'phy' must be of class 'phylo' or 'multiPhylo'")
  }

  # Correlation matrix for phylo object
  if (inherits(phy, "phylo")) {
    A <- ape::vcv.phylo(phy, corr = TRUE)
    return(A)
  }

  # Covariance matrix for multiphylo object
  if (inherits(phy, "multiPhylo")) {

    # Correlation matrix for each tree
    A_list <- parallel::mclapply(phy, function(x) {
      A <- ape::vcv.phylo(x, corr = TRUE)
      return(A)
    }, mc.cores = mc.cores)

    # Transform list of covariance matrices to 3d array
    A_array <- simplify2array(A_list)

    # Summarise these matrices
    A_m <- apply(A_array, 1:2, mean)
    A_sd <- apply(A_array, 1:2, stats::sd)

    # return list
    list(mean = A_m,
         sd = A_sd)
  }
}

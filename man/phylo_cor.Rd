% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo_cor.R
\name{phylo_cor}
\alias{phylo_cor}
\title{Phylogenetic correlation matrix for class 'phylo' or 'multiPhylo'}
\usage{
phylo_cor(phy, mc.cores = getOption("mc.cores", 1L))
}
\arguments{
\item{phy}{An object of class 'phylo' or 'multiPhylo'}

\item{mc.cores}{The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
Must be exactly 1 on Windows (which uses the master process)}
}
\value{
For class 'phylo', a correlation matrix; for class 'multiphylo',
a list with the mean and standard deviation correlation matrices
}
\description{
Phylogenetic correlation matrix for class 'phylo' or 'multiPhylo'
}
\examples{
\dontrun{
library(fishtree)
# Phylo
tree <- fishtree_phylogeny(rank = "Balistidae")
phylo_cor(tree)
# Multiphylo
trees <- fishtree_complete_phylogeny(rank = "Acanthuridae")
phylo_cor(trees)
}

}

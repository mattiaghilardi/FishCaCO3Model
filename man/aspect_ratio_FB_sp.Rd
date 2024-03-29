% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aspect_ratio_FB.R
\name{aspect_ratio_FB_sp}
\alias{aspect_ratio_FB_sp}
\title{Obtain aspect ratio of caudal fin from FishBase for taxa identified to genus level}
\usage{
aspect_ratio_FB_sp(species_list, mc.cores = getOption("mc.cores", 1L))
}
\arguments{
\item{species_list}{A vector of taxonomic names in the form "Genus sp." or "Genus spp."}

\item{mc.cores}{The number of cores to use for \code{\link[parallel:mclapply]{parallel::mclapply()}}.
Must be exactly 1 on Windows (which uses the master process)}
}
\value{
A tibble with the provided name, aspect ratio, and taxonomic level.
}
\description{
Obtain average species aspect ratio from FishBase at the lowest
available taxonomic level, either genus or family level.
}
\examples{
\dontrun{
aspect_ratio_FB_sp("Caranx sp.")
aspect_ratio_FB_sp(c("Acanthurus sp.", "Lutjanus spp."))
}

}

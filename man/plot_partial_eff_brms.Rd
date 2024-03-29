% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_partial_eff_brms.R
\name{plot_partial_eff_brms}
\alias{plot_partial_eff_brms}
\title{Plot partial effects of brms models}
\usage{
plot_partial_eff_brms(
  model,
  newdata,
  x,
  xlab,
  ylab,
  ndraws = 1000,
  seed = NULL,
  backtrans_y = NULL,
  backtrans_x = NULL,
  alpha = 0.2,
  color_mean = "black",
  color_draws = "grey80",
  rug = TRUE,
  data,
  color_rug = "black",
  facets = NULL,
  ...
)
}
\arguments{
\item{model}{A brmsfit object}

\item{newdata}{A data frame for which to make predictions}

\item{x}{Character, the predictor name (i.e, column name)}

\item{xlab}{A label for the x axis}

\item{ylab}{A label for the y axis}

\item{ndraws}{The number of draws to plot, or NULL to plot all draws}

\item{seed}{A seed to use when subsampling draws (i.e. when ndraw is not NULL)}

\item{backtrans_y}{Character, the type of backtransformation required for
the response variable if any. Available transformation at the moment are:
\itemize{
\item{"exp":} to backtransform from natural-log transformation
\item{"square":} to backtransform from square-root transformation
}
Default to NULL}

\item{backtrans_x}{Character, the type of backtransformation required for
the predictor if any. Available transformation at the moment are:
\itemize{
\item{"exp":} to backtransform from natural-log transformation
\item{"square":} to backtransform from square-root transformation
}
Default to NULL}

\item{alpha}{Numeric between 0 and 1 to set transparency for the 'ndraws' draws plotted}

\item{color_mean}{The colour of the line representing the average effect}

\item{color_draws}{The colour of 'ndraws' lines representing the uncertainty around the average}

\item{rug}{Logical, if a rug plot should be displayed on the x axis. Default to TRUE}

\item{data}{The data used to train the model. Only if "rug=TRUE"}

\item{color_rug}{The colour of the rug plot displaying the distribution of the predictor's values. Only if "rug=TRUE"}

\item{facets}{A set of variables or expressions quoted by \code{\link[ggplot2:vars]{ggplot2::vars()}} and
defining faceting groups on the rows or columns dimension}

\item{...}{Additional parameters passed to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}}
}
\value{
A ggplot object
}
\description{
Display partial effects of numeric predictors of a brmsfit object as spaghetti plots with the average overlaid.
}

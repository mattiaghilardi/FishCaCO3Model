#' Plot partial effects of brms models
#'
#' Display partial effects of numeric predictors of a brmsfit object as spaghetti plots with the average overlaid.
#'
#' @param model A brmsfit object
#' @param newdata A data frame for which to make predictions
#' @param x Character, the predictor name (i.e, column name)
#' @param xlab A label for the x axis
#' @param ylab A label for the y axis
#' @param ndraws The number of draws to plot, or NULL to plot all draws
#' @param seed A seed to use when subsampling draws (i.e. when ndraw is not NULL)
#' @param backtrans_y Character, the type of backtransformation required for
#'                    the response variable if any. Available transformation at the moment are:
#'                    \itemize{
#'                    \item{"exp":} to backtransform from natural-log transformation
#'                    \item{"square":} to backtransform from square-root transformation
#'                    }
#'                    Default to NULL
#' @param backtrans_x Character, the type of backtransformation required for
#'                    the predictor if any. Available transformation at the moment are:
#'                    \itemize{
#'                    \item{"exp":} to backtransform from natural-log transformation
#'                    \item{"square":} to backtransform from square-root transformation
#'                    }
#'                    Default to NULL
#' @param alpha Numeric between 0 and 1 to set transparency for the 'ndraws' draws plotted
#' @param color_mean The colour of the line representing the average effect
#' @param color_draws The colour of 'ndraws' lines representing the uncertainty around the average
#' @param rug Logical, if a rug plot should be displayed on the x axis. Default to TRUE
#' @param data The data used to train the model. Only if "rug=TRUE"
#' @param color_rug The colour of the rug plot displaying the distribution of the predictor's values. Only if "rug=TRUE"
#' @param facets A set of variables or expressions quoted by \code{\link[ggplot2:vars]{ggplot2::vars()}} and
#'               defining faceting groups on the rows or columns dimension
#' @param ... Additional parameters passed to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}
#'
#' @return A ggplot object
#'
#' @importFrom tidybayes add_epred_draws
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_line labs facet_wrap geom_rug
#'
plot_partial_eff_brms <- function(model,
                                  newdata,
                                  x,
                                  xlab, ylab,
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
                                  ...) {
  dat <- tidybayes::add_epred_draws(newdata = newdata, object = model, re_formula = NA,
                                    ndraws = ndraws, value = "pred", seed = seed)
  dat <- dplyr::mutate(dat, pred_m = mean(pred))

  if (!is.null(backtrans_y)) {
    if (backtrans_y == "square") {
      dat$pred <- dat$pred^2
      dat$pred_m <- dat$pred_m^2
    } else if (backtrans_y == "exp") {
      dat$pred <- exp(dat$pred)
      dat$pred_m <- exp(dat$pred_m)
    }
  }

  if (!is.null(backtrans_x)) {
    if (backtrans_x == "square") {
      dat[,x] <- dat[,x]^2
    } else if (backtrans_x == "exp")  {
      dat[,x] <- exp(dat[,x])
    }
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = dat[, x, drop = TRUE])) +
    ggplot2::geom_line(ggplot2::aes(y = pred, group = .draw), color = color_draws, alpha = alpha) +
    ggplot2::geom_line(ggplot2::aes(y = pred_m), color = color_mean) +
    ggplot2::labs(x = xlab, y = ylab)

  if (!is.null(facets)) {
    p <- p +
      ggplot2::facet_wrap(facets = facets, ...)
  }

  if (rug) {
    if (!is.null(backtrans_x)) {
      if (backtrans_x == "square") {
        data[,x] <- data[,x]^2
      } else if (backtrans_x == "exp")  {
        data[,x] <- exp(data[,x])
      }
    }
    p <- p +
      ggplot2::geom_rug(data = data, ggplot2::aes(x = data[, x, drop = TRUE]), sides = "b", color = color_rug)
  }

  p
}

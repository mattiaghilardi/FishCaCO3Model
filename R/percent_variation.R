#' Calculate % of variation in the response between two points along the range of a predictor "x"
#'
#' It computes the average % of variation in the response variable of a brms regression
#' model by the predictor "x"
#'
#' @param model A "brmsfit" object
#' @param x1 Data frame with a value of "x" and the mean value for the other continuous predictors (if any)
#' @param x2 Same as x1 but with a higher value of "x"
#' @param exp Logical, to back-transform when the response variable is modelled on the natural-log scale; default to TRUE
#' @param nsamples Positive integer indicating how many posterior samples should be used.
#'                 If NULL (the default) all samples are used.
#' @param summary Logical, if TRUE (the default) summary statistics are returned,
#'                if FALSE returns raw values
#' @param robust Only if "summary=TRUE". If FALSE (the default) the central tendency estimate is the mean,
#'               if TRUE is the median
#' @param probs Only if "summary=TRUE". The percentiles to be computed by the quantile function
#'
#' @return A value
#'
#' @importFrom stats quantile
#'
#' @export
percent_variation <- function(model, x1, x2, exp = TRUE, nsamples = NULL,
                              summary = TRUE, robust = FALSE, probs = c(0.025, 0.975)) {

  f1 <- fitted(model, x1, re_formula = NA, summary = FALSE, nsamples = nsamples)
  f2 <- fitted(model, x2, re_formula = NA, summary = FALSE, nsamples = nsamples)

  if (exp) {
    perc <- ((exp(f2)-exp(f1))/exp(f1))*100
  } else {
    perc <- ((f2-f1)/f1)*100
  }

  if (summary) {
    if (robust) {
      fun <- get("median", asNamespace("stats"))
    } else {
      fun <- get("mean", asNamespace("base"))
    }
    est <- fun(perc)
    q <- stats::quantile(perc, probs = probs)
    names(q) <- paste0("Q", gsub("%", "", names(q)))
    perc <- c(Estimate = est, q)
  }

  perc
}

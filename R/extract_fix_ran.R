#' Extract names of fixed and random effects from brms models
#'
#' @param model brmsfit object
#'
#' @importFrom insight get_predictors get_random
#' @importFrom dplyr case_when
#'
extract_fix_ran <- function(model) {
  pred <- insight::get_predictors(model)
  pred <- names(pred)
  if ("method" %in% pred) {
    pred <- pred[! pred == "method"]
    sigma <- "method"
  } else sigma <- "-"
  pred <- dplyr::case_when(pred == "log_weight" ~ "log(M)",
                           pred == "log_ril" ~ "log(RIL)",
                           pred == "mean_T" ~ "T",
                           pred == "sqrt_asp_ratio" ~ "sqrt(AR)",
                           pred == "mean_S" ~ "S",
                           pred == "sqrt_sampling" ~ "sqrt(ST)")
  f <- paste(pred, collapse = " + ")
  if (f == "") f <- "Intercept only"
  r <- insight::get_random(model)
  r <- names(r)
  if (is.null(r)) r <- "-"
  data.frame(fixed = f, random = r, sigma = sigma)
}

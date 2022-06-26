#' A function to predict excretion rate according to Wilson et al 2009 (p = 1 and p = 2.4) and our model
#'
#' @param w Body mass (g)
#' @param t Temperature (°C)
#' @param ril Relative intestinal length
#' @param ar Caudal fin aspect ratio
pred_C <- function(w, t, ril, ar){
  a <- 9.809e5 * w^0.75 * exp(-4727*(1/(273 + t)))
  b <- a * 2.4
  c <- exp(fitted(m_caco3_sigma,
                  newdata = data.frame(log_weight = log(w/1000),
                                       log_ril = log(ril),
                                       sqrt_asp_ratio = sqrt(ar),
                                       mean_T = t,
                                       method = "double"),
                  re_formula = NA)[, 1])
  cbind("C_wilson_p1" = a, "C_wilson_p2.4" = b, "C_ghilardi" = c)
}

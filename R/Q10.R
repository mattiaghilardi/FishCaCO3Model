#' Compute Q10 temperature coefficient
#' 
#' @param r1 rate at temperature 1
#' @param r2 rate at temperature 2
#' @param t1 temperature 1
#' @param t2 temperature 2
#' 
#' @examples 
#' \dontrun{
#' Q10(r1 = 2.2, r2 = 3.6, t1 = 15, t2 = 24)}
#' 
#' @export
Q10 <- function(r1, r2, t1, t2) {
  Q10 = (r2/r1)^(10/(t2-t1))
  Q10
}

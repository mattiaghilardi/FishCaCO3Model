#' Plot observed vs predicted values
#'
#' Display a scatterplot and linear regression of observed vs predicted values,
#' reporting the equation and R2 value.
#'
#' @param data A data frame containing observed and predicted values
#' @param obs A numeric vector, or character string if \code{data} is supplied, of observed values
#' @param pred A numeric vector, or character string if \code{data} is supplied, of predicted values
#' @param xlab A label for the x axis (predicted values)
#' @param ylab A label for the y axis (observed values)
#' @param smooth_size The size of the regression line
#' @param smooth_alpha The opacity of the confidence interval around the regression line
#' @param smooth_color The color of the regression line
#' @param smooth_fill The color of the confidence interval around the regression line
#' @param point_size The size of the points
#' @param point_alpha The opacity of the points
#' @param point_stroke The width of the points' border
#' @param point_shape The shape of the points
#' @param point_color The color of the points
#' @param point_fill The color of the inside of shapes with a border
#' @param text_size The size of the equation and R2
#' @param abline_size The size of the reference line
#' @param abline_color The color of the reference line
#' @param abline_type The type of reference line
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot aes_string geom_abline geom_smooth geom_point labs
#' @importFrom ggpubr stat_regline_equation
#'
plot_obs_vs_pred <- function(data = NULL, obs, pred, xlab, ylab,
                             smooth_size = 0.7, smooth_alpha = 0.3,
                             smooth_color = "#0D0887FF", smooth_fill = "#0D0887FF",
                             point_size = 0.8, point_alpha = 0.5, point_stroke = 0.1,
                             point_shape = 19, point_color = "#0D0887FF",
                             point_fill = "#0D0887FF",
                             text_size = 3, abline_size = 0.7,
                             abline_color = "red", abline_type = 2) {

  ggplot2::ggplot(data = data, mapping = ggplot2::aes_string(x = pred, y = obs)) +
    ggplot2::geom_abline(linetype = abline_type, color = abline_color, size = abline_size) +
    ggplot2::geom_smooth(method = "lm", size = smooth_size, alpha = smooth_alpha,
                         color = smooth_color, fill = smooth_fill) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, shape = point_shape,
                        stroke = point_stroke, color = point_color, fill = point_fill) +
    ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                                  family = "serif", size = text_size) +
    ggplot2::labs(x = xlab, y = ylab)
}

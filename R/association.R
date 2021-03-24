#' Measure of Spearman Correlation
#'
#' @param x numeric vector
#' @param y numeric vector
#' @return double
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe)
#'   anscombe_tidy <- anscombe %>%
#'   pivot_longer(cols = everything(),
#'     names_to = c(".value", "set"),
#'     names_pattern = "(.)(.)")
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_monotonic(anscombe$x1, anscombe$y1)
#'   sc_monotonic(anscombe$x2, anscombe$y2)
#'   sc_monotonic(anscombe$x3, anscombe$y3)
#'   sc_monotonic(anscombe$x4, anscombe$y4)
#' @export
sc_monotonic <- function(x, y){
  cor(x, y, method='spearman')
}

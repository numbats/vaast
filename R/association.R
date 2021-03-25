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
  stats::cor(x, y, method='spearman')
}

#' Spline based index.
#'
#' (Taken from tourr git repo)
#' Compares the variance in residuals of a fitted
#' spline model to the overall variance to find
#' functional dependence in 2D projections
#' of the data.
#'
#' @param x numeric vector
#' @param y numeric vector
#' @examples
#'
#' @export
sc_splines <- function(x,y) {
  dat <- matrix(c(x,y), nrow=length(x), byrow=FALSE)
  tourr::splines2d(dat)
}

#' Distance correlation index.
#'
#'(Taken from tourr git repo)
#' Computes the distance correlation based index on
#' 2D projections of the data.
#'
#' @keywords hplot
#' @importFrom stats na.omit
#' @param x numeric vector
#' @param y numeric vector
#' @export
#sc_dcor <- function(x,y) {
#  dat <- matrix(c(x,y), nrow=length(x), byrow=FALSE)
#  tourr::dcor2d(dat)
#}
sc_dcor <- function() {
  function(x,y) {
    xy <- na.omit(data.frame(x = x, y = y))
    measure <- with(xy, energy::dcor(x, y))
    return(measure)
  }
}

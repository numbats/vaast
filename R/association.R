#' Measure of Spearman Correlation
#'
#' @param x numeric vector
#' @param y numeric vector
#'
#' @return double
#' @export

sc_monotonic <- function(x, y){
  cor(x, y, method='spearman')
}

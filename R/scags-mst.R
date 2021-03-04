#' @export
sc_stringy <- function(x, y) UseMethod("sc_stringy")

#' @export
sc_stringy.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  mst <- x$mst(x$del, x$weights)
  vertex_counts <- igraph::degree(mst)
  sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
}

#' @export
sc_stringy <- function(x, y) UseMethod("sc_stringy")

#' @export
sc_stringy.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  mst <- gen_mst(x$del, x$weights)
  vertex_counts <- igraph::degree(mst)
  sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
}

sc_stringy.default <- function(x, y){
  sc <- scree(x, y)
  sc_stringy.scree(sc)
}

gen_mst <- function(del, weights) {
  edges <- del$mesh[, c("ind1", "ind2")]
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  igraph::mst(graph, weights =  igraph::E(graph)$weight)
}

#' Compute stirated scagnostic measure using MST
#'
#' @examples
#' @export
sc_striated <- function(x, y) UseMethod("sc_striated")

#' @export
sc_striated.scree <- function(x, y = NULL) {
  #make mst
  mst <- gen_mst(x$del, x$weights)
  #find verts with >1 connection
  vertex_counts <- igraph::degree(mst)
  weights <- igraph::E(mst)$weight
  #cos on adjacent edges in mst
  #sum of indicator func on <0.75 (=1 for ang 0 to 40 deg)
  #1/|v|scale

  }

sc_striated.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated.scree(sc)
}

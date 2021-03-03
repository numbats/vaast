#' Generate an scatter tree
#'
#' @param x,y numeric vectors
#' @param num_bins (default = 30)
#' @param ...  other args
#'
#' @return A list of graph objects with three components
#'  - minimum spanning tree
#'  - convex hull
#'  - concave hull (alpha hull)
#' @export
#' @importFrom igraph E
scree <- function(x, y, num_bins = 30, ...) {
  stopifnot(is.numeric(x),
            is.numeric(y),
            length(x) == length(y))

  # should x,y be rescaled?
  hb <- hexbin::hexbin(x, y, xbins = num_bins)

  # compute the triangulation using the hexbin centres
  del <- RTriangle::triangulate(
    RTriangle::pslg(cbind(hb@xcm, hb@ycm))
  )

  # create a graph from the edges of the triangulation
  graph <- igraph::graph_from_edgelist(del$E, directed = FALSE)
  # the vertices here are the indices of hexbinned points
  # we also add the counts
  graph <- igraph::set_vertex_attr(graph, "x", value = del$P[,1])
  graph <- igraph::set_vertex_attr(graph, "y", value = del$P[,2])
  graph <- igraph::set_vertex_attr(graph, "counts", value = hb@count)
  # add edge weights according to distance between edges of points
  graph <- igraph::set_edge_attr(graph,
                                 "weight",
                                 value = edge_weights(del$P, del$E))

  # convex hull is given by the segments of del, i.e. del$S
  chull <- igraph::graph_from_edgelist(del$S, directed = FALSE)

  # concave hull is computed via triangulation

  # spanning tree, weighted by the lengths between points
  mst<- igraph::mst(graph, weights = E(graph)$weight)

  list(mst, chull, ahull = list())

}

edge_weights <- function(points, edges) {
  sqrt(rowSums((points[edges[,1],] - points[edges[,2], ])^2))
}

# generate the edges of the alpha shape graph using Edelsbrunner's algorithm
# an alternative would be to use alphahull package for all steps,
# including the triangulation
alpha_hull <- function(del, alpha) {

}


# This is the edge filter from Wilkinson 05
psi <- function(w, q = c(0.25, 0.75)) {
  q <- quantile(w, probs = q)
  q[2] + 1.5 * diff(q)
}

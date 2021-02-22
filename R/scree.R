#' Generate an scatter tree
#'
#' @param x,y numeric vectors
#' @param num_bins (default = 30)
#' @param ...  other args
#'
#' @return A list of graph objects with three components
#'  - minimum spanning tree
#'  -

scree <- function(x, y, num_bins = 30, ...) {
  stopifnot(is.numeric(x),
            is.numeric(y),
            length(x) == length(y))

  # should x,y be rescaled?
  hb <- hexbin::hexbin(x,y, xbins = num_bins)

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
  chull <- graph_from_edgelist(del$S, directed = FALSE)

  # concave hull is computed via triangulation
  ahull <- graph_from_edgelist(
    alpha_shape(del$P, del$T, alpha = alpha_param(E(graph)$weight)),
    directed = FALSE
  )

  # spanning tree
  mst<- igraph::mst(graph, weights = E(graph)$weight)

  list(mst, chull, ahull)

}

edge_weights <- function(points, edges) {
  sqrt(rowSums((points[edges[,1],] - points[edges[,2], ])^2))
}

# generate the alpha hull, not quite right
alpha_shape <- function(points, tri, alpha) {
  # lengths of each side of the triangle
  a <- edge_weights(points, tri[, c(1,2)])
  b <- edge_weights(points, tri[, c(2,3)])
  c <- edge_weights(points, tri[, c(3,1)])

  # semi-perimeter of triangle
  s <- (a + b + c) / 2

  # area and circumference
  area <- sqrt(s * (s - a) * (s - b) * (s - c))
  circum <- (a * b * c) / (4 * area)

  # filter if below alpha
  filter_tri <- tri[circum < (1/alpha), ]
}

# from 05 Wilkinson paper
alpha_param <- function(w, q = c(0.25, 0.75)) {
  q <- quantile(w, probs = q)
  q[2] + 1.5 * diff(q)
}

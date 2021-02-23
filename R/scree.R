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
alpha_hull <- function(del, alpha) {
  stopifnot(alpha >= 0)

  # first find the distance from points on the edges of the
  # delauney triangle to the points in the voroni tesselation
  points_from <- del$P[del$E[,1], ]
  points_to <- del$P[del$E[,2], ]

  # this isn't quite right!
  voroni_from <- del$VP[del$E[,1], ]
  voroni_to <- del$VP[del$E[,2], ]

  # if voroni tesselation has infinite ray set the distance to Inf
  dm1 <- sqrt(rowSums((points_from - voroni_from)^2))
  dm2 <- sqrt(rowSums((points_from - voroni_to)^2))

  dm1[which(del$VE[,1] == -1)] <- Inf
  dm2[which(del$VE[,2] == -1)] <- Inf

  # if you're on the convex hull this will form part of the alpha complex
  alpha_extremes <- unique(as.vector(del$S))
  in_chull <- setdiff(seq_len(nrow(del$P)), alpha_extremes)

  if (length(in_chull) > 0) {
    # augment edges with distances
    aug_dist <- rbind(
      cbind(del$E, dm1),
      cbind(del$E, dm2),
      cbind(del$E[, c(2,1)], dm1),
      cbind(del$E[, c(2,1)], dm2)
    )
    # compute the max distance for each unique edge
    id <- factor(aug_dist[,1])
    alpha_maxes <- tapply(aug_dist[,3], id, max)
    # filter again by the alpha value
    alpha_inx <-  in_chull[na.omit(match(unique(aug_dist[,1])[alpha < alpha_maxes], in_chull))]
    # combine to find indices
    alpha_extremes <- c(alpha_extremes, alpha_inx)
  }

  alpha_edge_from <- match(del$E[,1], alpha_extremes)
  alpha_edge_to <- match(del$E[,2], alpha_extremes)
  valid_edge <- which(alpha_edge_from & alpha_edge_to)

  alpha_edges <- del$E[valid_edge, ]

  midpoints <- (points_from[alpha_edges[,1], ]+ points_to[alpha_edges[,2], ]) / 2

  dist_to_edge <- sqrt(rowSums((points_from[alpha_edges[,1], ]+ points_to[alpha_edges[,2], ])^2))

  check_segment <- cbind(
    voroni_to[alpha_edges, ],
    voroni_from[alpha_edges, ],
    midpoints
  )

  col_rank <- function(x) { rank(x)[length(x)] == 2}
  betw <- check_segment[,1] == check_segment[,3] &
    apply(check_segment[, c(2,4,6)], 1, col_rank) |
    apply(check_segment[, c(2,4,6)], 1, col_rank)
  betw <- as.integer(betw)
  betw[betw == 0] <- NA_integer_

  boundary_distances <- cbind(dm1[valid_edge],
                              dm2[valid_edge],
                              dist_to_edge * betw)

  boundary_min <- apply(boundary_distances, 1, min, na.rm = TRUE)
  boundary_max <- apply(boundary_distances, 1, max, na.rm = TRUE)
  valid_alpha <- (boundary_min <= alpha & alpha <= boundary_max)

  alpha_edges[valid_alpha, ]

}



# This is the edge filter from Wilkinson 05
psi <- function(w, q = c(0.25, 0.75)) {
  q <- quantile(w, probs = q)
  q[2] + 1.5 * diff(q)
}

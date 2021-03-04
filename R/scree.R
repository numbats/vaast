#' Generate an scatter tree
#'
#' @param x,y numeric vectors
#' @param binner a function that
#' @param ...  other args
#'
#' @return A list of graph objects with three components
#'  - minimum spanning tree
#'  - convex hull
#'  - concave hull (alpha hull)
#' @export
scree <- function(x, y, binner = NULL, ...) {
  # checks on x,y
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y)
  )

  if (!(is.null(binner) | is.function(binner)))
    stop("binner must be a function")

  # cast to a matrix
  xy <- cbind(unitize(x), unitize(y))

  if (is.function(binner)) {
    xy <- binner(xy)
  }

  # compute delauney triangulation
  del <- alphahull::delvor(xy)

  # edge weights from the triangulation
  weights <- gen_edge_lengths(del)

  # alpha estimator
  alpha <- psi(weights)

  structure(
    list(del = del,
         weights  = weights,
         alpha = alpha,
         mst = gen_mst,
         alpha_hull = gen_alpha_hull,
         conv_hull = gen_conv_hull
    ),
    class = "scree")
}


gen_edge_lengths <- function(del) {
  from_cols <- c("x1", "y1")
  to_cols <- c("x2", "y2")
  sqrt(rowSums((del$mesh[, from_cols] - del$mesh[, to_cols])^2))
}

gen_mst <- function(del, weights = gen_edge_lengths(del)) {
  edges <- del$mesh[, c("ind1", "ind2")]
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  igraph::mst(graph, weights =  igraph::E(graph)$weight)
}


gen_conv_hull <- function(del) {
  tripack::convex.hull(del$tri.obj)
}

gen_alpha_hull <- function(del, alpha) {
  alphahull::ahull(del, alpha = alpha)
}

# rescale input to lie in unit interval
unitize <- function(x, na.rm = TRUE) {
  rng <- range(x, na.rm = na.rm)
  (x - rng[1]) / diff(rng)
}

# This is the edge filter from Wilkinson 05
psi <- function(w, q = c(0.25, 0.75)) {
  q <- quantile(w, probs = q)
  unname(q[2] + 1.5 * diff(q))
}

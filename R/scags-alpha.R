
#' @export
sc_convex <- function(x, y) UseMethod("sc_convex")

#' @export
sc_convex.scree <- function(x,y = NULL) {
  stopifnot(is.null(y))
  chull <- x$conv_hull(x$del)
  chull_area <- splancs::areapl(cbind(chull$x, chull$y))
  ahull <- x$alpha_hull(x$del, x$alpha)
  ahull_area <- alphahull::areaahull(ahull)
  ahull_area / chull_area
}

#' @export
sc_skinny <- function(x, y) UseMethod("sc_skinny")

#' @export
sc_skinny.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  ahull <- x$alpha_hull(x$del, x$alpha)
  ahull_area <- alphahull::areaahull(ahull)
  1 - sqrt(4*pi * ahull_area) / ahull$length
}

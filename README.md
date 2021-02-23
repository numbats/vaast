
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vaast

<!-- badges: start -->

<!-- badges: end -->

A rewrite of scagnostics in pure R with minimal dependencies.

A 2d scatter plot is represented with three graph objects which can all
be computed from a `RTriangle::triangulate`:

1.  A minimum spanning tree weighted by Delauney triangle lengths
2.  Convex hull i.e. the outer segments of the triangulation
3.  Concave hull (alpha shape graph) i.e. how many triangles are
    enclosed with a ball with radius alpha.

And various measures are computed against those graphs.


<!-- README.md is generated from README.Rmd. Please edit that file -->

# -

<!-- badges: start -->
<!-- badges: end -->

**scagnostics** or scatter plot diagnostics are a technique for
automatically extracting interesting visual features from pair-wise
combinations of variables in data. This package is an implementation of
**graph theoretic scagnostics** in pure R.

## Implementation Details

### Graph based

A 2-d scatter plot can be represented by a combination of three graphs
which are computed directly from the Delauney-Voroni tesselation.

1.  A **minimum spanning tree** weighted by the lengths of the Delauney
    triangles
2.  The **convex hull** of the points i.e. the outer segments of the
    triangulation
3.  The **alpha hull ** (also called concave hull) i.e. formed by
    connect the outer edges of triangles that are enclosed within a ball
    of radius *alpha*.

All graph based scagnostic measures are computed with respect to these
three graphs.

Prior to graph construction decisions must be made about filtering
outliers ( done with respect to the distribution of edge lengths in the
triangulation ) and thinning the size of the graphs by performing
binning (for computational speed). For the moment we can forge ahead
without these but it is worth keeping in mind that the package needs to
be flexible enough to include them.  
There are also opportunities to experiment with the preprocessing here.

Two MST measures “clumpy” and “outlying” are known to cause problems

As all the graph based measures rely on the triangulation, they could be
computed lazily. More concretely, if you are only interested in
computing “skinny” you don’t need to compute the spanning tree. If you
computed “skinny” but wanted to then compute “convex” you shouldn’t need
to reconstruct the alpha-hull and so on… To begin let’s not worry about
that and focus on implementations of each measure.

### Association based

These are computed directly from the 2-d point clouds, and do not need
to be constructed from the graph.

See issue <https://github.com/numbats/vaast/issues/2> for some thoughts.

## User Interface

At the moment we have a very low level API. For pair-wise data we can
create the triangulation and other bits necessary using `scree`

``` r
x <- runif(100)
y <- runif(100)
scr <- scree(x, y)
```

Then the ‘scree’ object can be passed to various functions to compute
scagnostics for a single combination of variables:

``` r
sc_skinny(scr)
sc_stringy(scr)
```

It would also be nice to be able to do:

``` r
sc_skinny(x, y)
```

And design an interface that could handle mixing and matching various
scagnostics together.

``` r
anscombe %>% 
  sc_estimate(.vars = everything(), 
              .features = list(clumpy = sc_clumpy, skinny = sc_skinny))

#> output would look something like a tibble

#> from_var to_var clumpy skinny
#> "x1" "y1" 0.3 0.2
#> "x1", "y2", ...

anscombe %>% 
  sc_estimate(.vars = cross(x1:x4, y1:y4), 
              .features = list(clumpy = sc_clumpy, skinny = sc_skinny))

anscombe %>% 
  sc_estimate(.vars = pairwise(x1:x4, y1:y4),
              .features = list(clumpy = sc_clumpy, skinny = sc_skinny))
```


#' Compute stringy scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_stringy(anscombe$x1, anscombe$y1)
#'   sc_stringy(anscombe$x2, anscombe$y2)
#'   sc_stringy(anscombe$x3, anscombe$y3)
#'   sc_stringy(anscombe$x4, anscombe$y4)
#'
#' @export
sc_stringy <- function(x, y) UseMethod("sc_stringy")

#' @rdname sc_stringy
#' @export
sc_stringy.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  mst <- gen_mst(x$del, x$weights)
  sc_stringy.igraph(mst)
}

#' @rdname sc_stringy
#' @export
sc_stringy.default <- function(x, y){
  sc <- scree(x, y)
  sc_stringy.scree(sc)
}

#' @rdname sc_stringy
#' @export
sc_stringy.igraph <- function(mst){
  vertex_counts <- igraph::degree(mst)
  sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
}

#' Compute stirated scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_striated(anscombe$x1, anscombe$y1)
#'   sc_striated(anscombe$x2, anscombe$y2)
#'   sc_striated(anscombe$x3, anscombe$y3)
#'   sc_striated(anscombe$x4, anscombe$y4)
#'
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_striated(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_striated(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_striated(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_striated(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'
#' @export
sc_striated <- function(x, y) UseMethod("sc_striated")

#' @rdname sc_striated
#' @export
sc_striated.scree <- function(x, y = NULL) {
  mst <- gen_mst(x$del, x$weights)
  sc_striated.igraph(mst, x)

  }

#' @rdname sc_striated
#' @export
sc_striated.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated.scree(sc)
}

#' @rdname sc_striated
#' @export
sc_striated.list <- function(mst, x){
  vertex_counts <- igraph::degree(mst)
  angs <- which(vertex_counts==2)
  angles_vect <- numeric(length(angs))
  for(i in seq_len(length(angs))){
    adjs <- which(mst[angs[i]]>0)
    points <- x$del$x[adjs,]
    origin <- x$del$x[angs[i],]
    vects <- t(t(points)-origin)
    angles_vect[i] <- (vects[1,]%*%vects[2,])/(prod(mst[angs[i]][adjs]))
  }
  (sum(ifelse(angles_vect<(-0.75),1,0)))/length(vertex_counts)
}

#' Compute clumpy scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_clumpy(anscombe$x1, anscombe$y1)
#'   sc_clumpy(anscombe$x2, anscombe$y2)
#'   sc_clumpy(anscombe$x3, anscombe$y3)
#'   sc_clumpy(anscombe$x4, anscombe$y4)
#'
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_clumpy(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_clumpy(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_clumpy(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_clumpy(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'
#' @export
sc_clumpy <- function(x, y) UseMethod("sc_clumpy")

#' @rdname sc_clumpy
#' @export
sc_clumpy.scree <- function(x, y = NULL) {
  mymst <- gen_mst(x$del, x$weights)
  sc_clumpy.igraph(mymst,x)
}

#' @rdname sc_clumpy
#' @export
sc_clumpy.default <- function(x, y){
  sc <- scree(x, y)
  sc_clumpy.scree(sc)
}

#' @rdname sc_clumpy
#' @export
sc_clumpy.list <- function(mymst, x){
  mstmat <- matrix(mymst[], nrow=length(x[["del"]][["x"]][,1]))
  mstmat[upper.tri(mstmat, diag = FALSE)]=0
  #get cols and rows of each value
  edges <- which(mstmat>0)
  rows <- edges %% length(mstmat[,1])
  cols <- (edges-1) %/% length(mstmat[,1]) +1
  clumpy <- rep(0,length(edges))
  for(i in seq(length(edges))){
    #set value for edge in consideration
    jval <- mstmat[edges[i]]
    #reset inner loop values and remove target edge
    inedges <- edges[-i]
    inrows <- rows[-i]
    incols <- cols[-i]
    #make cluster with remaining edges (this second for loop is bad computationally oof)
    for(j in seq(length(inedges))){
      if(j==1){
        ind=1
      }
      checkvec <- c(inrows[ind],incols[ind]) #lower triang matrix so need to check rows and cols
      if((inrows[j]%in%checkvec | incols[j]%in%checkvec)& (j!=1)){
        ind <- c(ind,j)
      }
    }
    cluster1 <- mstmat[inedges[ind]]
    cluster2 <- mstmat[inedges[-ind]]
    kval <- ifelse(sum(cluster1)<sum(cluster2), pmax(cluster1), pmax(cluster2))
    clumpy[i] <- 1- (kval/jval)
  }
  clumpy <- clumpy[which(!is.na(clumpy))]
  max(clumpy)

}



#' Compute sparse scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_sparse(anscombe$x1, anscombe$y1)
#'   sc_sparse(anscombe$x2, anscombe$y2)
#'   sc_sparse(anscombe$x3, anscombe$y3)
#'   sc_sparse(anscombe$x4, anscombe$y4)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_sparse(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_sparse(datasaurus_dozen_wide$bullseye_x, datasaurus_dozen_wide$bullseye_y)
#'   sc_sparse(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_sparse(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_sparse(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_sparse(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$high_lines_x, datasaurus_dozen_wide$high_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$slant_down_x, datasaurus_dozen_wide$slant_up_y)
#'   sc_sparse(datasaurus_dozen_wide$star_x, datasaurus_dozen_wide$star_y)
#'   sc_sparse(datasaurus_dozen_wide$v_lines_x, datasaurus_dozen_wide$v_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$wide_lines_x, datasaurus_dozen_wide$wide_lines_y)
#'   sc_sparse(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_sparse <- function(x, y) UseMethod("sc_sparse")

#' @rdname sc_sparse
#' @export
sc_sparse.scree <- function(x, y = NULL) {
  #generate vector of MST edges
  mymst <- gen_mst(x$del, x$weights)
  sc_sparse.igraph(mymst,x)
}

#' @rdname sc_sparse
#' @export
sc_sparse.default <- function(x, y){
  sc <- scree(x, y)
  sc_sparse.scree(sc)
}

#' @rdname sc_sparse
#' @export
sc_sparse.list <- function(mymst, x){
  mstmat <- matrix(mymst[], nrow=length(x[["del"]][["x"]][,1]))
  mstmat[upper.tri(mstmat, diag = FALSE)]=0
  edges <- sort(mstmat[which(mstmat>0)])
  #calculate sparse value
  edges[floor(0.9*length(edges))]
}


#' Compute skewed scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_skewed(anscombe$x1, anscombe$y1)
#'   sc_skewed(anscombe$x2, anscombe$y2)
#'   sc_skewed(anscombe$x3, anscombe$y3)
#'   sc_skewed(anscombe$x4, anscombe$y4)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_skewed(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_skewed(datasaurus_dozen_wide$bullseye_x, datasaurus_dozen_wide$bullseye_y)
#'   sc_skewed(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_skewed(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_skewed(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_skewed(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$high_lines_x, datasaurus_dozen_wide$high_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$slant_down_x, datasaurus_dozen_wide$slant_up_y)
#'   sc_skewed(datasaurus_dozen_wide$star_x, datasaurus_dozen_wide$star_y)
#'   sc_skewed(datasaurus_dozen_wide$v_lines_x, datasaurus_dozen_wide$v_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$wide_lines_x, datasaurus_dozen_wide$wide_lines_y)
#'   sc_skewed(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_skewed <- function(x, y) UseMethod("sc_skewed")

#' @rdname sc_skewed
#' @export
sc_skewed.scree <- function(x, y = NULL) {
  #generate vector of MST edges
  mymst <- gen_mst(x$del, x$weights)
  sc_skewed.igraph(mymst, x)
}

#' @rdname sc_skewed
#' @export
sc_skewed.default <- function(x, y){
  sc <- scree(x, y)
  sc_skewed.scree(sc)
}

#' @rdname sc_skewed
#' @export
sc_skewed.list <- function(mymst, x){
  mstmat <- matrix(mymst[], nrow=length(x[["del"]][["x"]][,1]))
  mstmat[upper.tri(mstmat, diag = FALSE)]=0
  edges <- sort(mstmat[which(mstmat>0)])
  #calculate skewed value
  q90 <- edges[floor(0.9*length(edges))]
  q50 <- edges[floor(0.5*length(edges))]
  q10 <- edges[floor(0.1*length(edges))]
  (q90-q50)/(q90/q50) #no size adjustment
}


#' Compute outlying scagnostic measure using MST
#'
#' @examples
#'   require(ggplot2)
#'   require(tidyr)
#'   require(dplyr)
#'   data(anscombe_tidy)
#'   ggplot(anscombe_tidy, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~set, ncol=2, scales = "free")
#'   sc_outlying(anscombe$x1, anscombe$y1)
#'   sc_outlying(anscombe$x2, anscombe$y2)
#'   sc_outlying(anscombe$x3, anscombe$y3)
#'   sc_outlying(anscombe$x4, anscombe$y4)
#'   ggplot(datasaurus_dozen, aes(x=x, y=y)) +
#'     geom_point() +
#'     facet_wrap(~dataset, ncol=3, scales = "free")
#'   sc_outlying(datasaurus_dozen_wide$away_x, datasaurus_dozen_wide$away_y)
#'   sc_outlying(datasaurus_dozen_wide$bullseye_x, datasaurus_dozen_wide$bullseye_y)
#'   sc_outlying(datasaurus_dozen_wide$circle_x, datasaurus_dozen_wide$circle_y)
#'   sc_outlying(datasaurus_dozen_wide$dino_x, datasaurus_dozen_wide$dino_y)
#'   sc_outlying(datasaurus_dozen_wide$dots_x, datasaurus_dozen_wide$dots_y)
#'   sc_outlying(datasaurus_dozen_wide$h_lines_x, datasaurus_dozen_wide$h_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$high_lines_x, datasaurus_dozen_wide$high_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$slant_down_x, datasaurus_dozen_wide$slant_up_y)
#'   sc_outlying(datasaurus_dozen_wide$star_x, datasaurus_dozen_wide$star_y)
#'   sc_outlying(datasaurus_dozen_wide$v_lines_x, datasaurus_dozen_wide$v_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$wide_lines_x, datasaurus_dozen_wide$wide_lines_y)
#'   sc_outlying(datasaurus_dozen_wide$x_shape_x, datasaurus_dozen_wide$x_shape_y)
#'
#' @export
sc_outlying <- function(x, y) UseMethod("sc_outlying")

#' @rdname sc_outlying
#' @export
sc_outlying.scree <- function(x, y = NULL) {
  #generate vector of MST edges
  mymst <- gen_mst(x$del, x$weights)
  sc_outlying.igraph(mymst, x)
}

#' @rdname sc_outlying
#' @export
sc_outlying.default <- function(x, y){
  sc <- scree(x, y)
  sc_outlying.scree(sc)
}

#' @export
sc_outlying.list <- function(mymst, x){
  #input: original mst (mymst) and scree object (x)
  #output: outlying mst value

  #make into matrix
  mstmat <- matrix(mymst[], nrow=length(x[["del"]][["x"]][,1]))
  #identify outliers
  outliers <- outlying_identify(mymst, x)
  #calculate outlying value
  numer <- sum(mstmat[outliers,]) #sum of edges of outlying points
  denom <- sum(mstmat)
  numer/denom
}

gen_mst <- function(del, weights) {
  edges <- del$mesh[, c("ind1", "ind2")]
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  igraph::mst(graph, weights =  igraph::E(graph)$weight)
}

outlying_identify <- function(mst, scr){
  #input: takes a mst and scree
  #output: rown number of

  #make into upper tri-matrix
  mstmat <- matrix(mst[], nrow=length(scr[["del"]][["x"]][,1]))
  mstmat_diag <- mstmat
  mstmat_diag[upper.tri(mstmat, diag = FALSE)]=0

  #calculate w value
  edges <- sort(mstmat_diag[which(mstmat>0)])
  q25 <- edges[floor(0.25*length(edges))]
  q75 <- edges[floor(0.75*length(edges))]
  w <- q75 + 1.5*(q75-q25)

  #set values above w to 0 in matrix to find outlying
  mstmat_check <- mstmat
  mstmat_check[mstmat>w]=0

  #row sum of matrix, if 0 all edges are above this value
  rowsum <- mstmat_check%*%rep(1, length(mstmat_check[1,]))
  #return the outlying observation
  which(rowsum==0)
}

original_and_robust <- function(x,y){
  #input: data for 2 variables x and y
  #output: list of scree and MST objects

  #construct original scree and MST
  sc_original <- scree(x,y)
  mst_original <- gen_mst(sc_original$del, sc_original$weights)

  #identify outliers and remove them from the data
  outliers <- outlying_identify(mst_original, sc_original)
  new_x <- x[-outliers]
  new_y <- y[-outliers]

  #recalculate scree and MST
  sc_robust <- scree(new_x,new_y)
  mst_robust <- gen_mst(sc_robust$del, sc_robust$weights)

  #output 4 objects as a list
  structure(
    list(
      scree_ori = sc_original,
      mst_ori  = mst_original,
      scree_rob = sc_robust,
      mst_rob = mst_robust
    ))
}

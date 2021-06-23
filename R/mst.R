
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

#' @export
sc_stringy.scree <- function(x, y = NULL) {
  stopifnot(is.null(y))
  mst <- gen_mst(x$del, x$weights)
  vertex_counts <- igraph::degree(mst)
  sum(vertex_counts == 2) / (length(vertex_counts) - sum(vertex_counts == 1))
}

#' @export
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

#' @export
sc_striated.scree <- function(x, y = NULL) {
  mst <- gen_mst(x$del, x$weights)
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
  striated <- (sum(ifelse(angles_vect<(-0.75),1,0)))/length(vertex_counts)
  return(striated)
  }

#' @export
sc_striated.default <- function(x, y){
  sc <- scree(x, y)
  sc_striated.scree(sc)
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

#' @export
sc_clumpy.scree <- function(x, y = NULL) {
  mymst <- gen_mst(x$del, x$weights)
  mstmat <- matrix(mymst[], nrow=11)
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
    print(kval)
    print(jval)
    clumpy[i] <- 1- (kval/jval)
  }
  print(clumpy)
  return(pmax(clumpy))
}

#' @export
sc_clumpy.default <- function(x, y){
  sc <- scree(x, y)
  sc_clumpy.scree(sc)
}


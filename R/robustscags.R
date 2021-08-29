

inside_cluster <- function(mst, sc){
  #mst matrix
  matlist <- twomstmat(mst,scr)#function for matrix version of mst
  mstmat <- matlist$mat #mst matrix
  mstmat_lt <- matlist$lowertri #mst lower triangular matrix

  #calculate w
  edges <- mstmat_lt[which(mstmat>0)]
  w <- psi(edges)

  #remove outliers
  outliers <- outlying_identify(mst, sc)

  #calculate clumpy values
  value <- clumpy_value()
  #return ej, clumpy val, w
}



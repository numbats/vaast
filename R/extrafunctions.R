#' Compute scagnostics on all possible scatter plots for the given data
#'
#' @examples
#'
#' @export
sc_pairwise <- function(all_data, scags=c("outlying","stringy", "striated", "clumpy", "sparse","skewed","convex","skinny","monotonic","splines","dcor"), groups=NULL){
  all_combs <- expand.grid(colnames(all_data),colnames(all_data))%>%
    filter(!(Var1==Var2))
  all_combs <- cbind(all_combs, t(apply(all_combs, 1, intermediate_scags, data=all_data, scags=scags)))
  colnames(all_combs) <- c("Var1", "Var2", scags)
  return(all_combs)
  }

#' Intermediate function to use apply on to calculate scagnostics
#'
#'
#' @export
intermediate_scags <- function(vars, data, scags){
  #fakefunc <- function(vars, data){sc_convex(pull(data, var=vars[[1]]), pull(data, var=vars[[2]]))}
  x <- pull(data, var=vars[[1]])
  y <- pull(data, var=vars[[2]])
  return(calc_scags(x,y,scags))
#return(calc_scags(x,y, scags))
}


#' Compute selected scagnostics
#'
#' @examples
#'
#' @export
calc_scags <- function(x, y, scags=c("outlying","stringy", "striated", "clumpy", "sparse","skewed","convex","skinny","monotonic","splines","dcor")){
  #set all scagnostics to null
  outlying = NULL
  stringy = NULL
  striated = NULL
  clumpy = NULL
  sparse = NULL
  skewed = NULL
  convex = NULL
  skinny = NULL
  monotonic = NULL
  splines = NULL
  dcor = NULL
  #OUTLYING ADJUSTMENT
  #original_scree <- scree(x,y)
  #original_mst <- gen_mst(x,y)
  #should we calculate new MST?

  #CALCULATE ALPHA HULL MEASURES


  #CALCULATE ASSOCIATION MEASURES
  if("monotonic" %in% scags){
    monotonic <- sc_monotonic(x,y)
  }
  if("splines"%in%scags){
    splines <- sc_splines(x,y)
  }
  if("dcor" %in% scags){
    dcor <- sc_splines(x,y)
  }
  scagnostic_calcs <- c(outlying,stringy, striated, clumpy, sparse, skewed, convex, skinny, monotonic, splines, dcor)
  return(scagnostic_calcs)
}




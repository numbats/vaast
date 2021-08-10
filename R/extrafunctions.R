#' Compute scagnostics on all possible scatter plots for the given data
#'
#' @examples
#'
#' @export
sc_pairwise <- function(all_data, scags=c("outlying","stringy", "striated", "clumpy", "sparse","skewed","convex","skinny","monotonic","splines","dcor"), groups=NULL){
  all_combs <- expand.grid(colnames(all_data),colnames(all_data))%>%
    filter(!(Var1==Var2))
  #get rid of reversed duplicates
  all_combs <- all_combs[!duplicated(apply(all_combs,1,function(x) paste(sort(x),collapse=''))),]
  if(length(scags)==1){
    all_combs <- cbind(all_combs, apply(all_combs, 1, vaast:::intermediate_scags, data=all_data, scags=scags))
  }
  else{
    all_combs <- cbind(all_combs, t(apply(all_combs, 1, vaast:::intermediate_scags, data=all_data, scags=scags)))
  }
  scags_name <- c("outlying","stringy", "striated", "clumpy", "sparse","skewed","convex","skinny","monotonic","splines","dcor")
  scags_name <- scags_name[which(scags_name %in% scags)]
  colnames(all_combs) <- c("Var1", "Var2", scags_name)
  return(all_combs)
  }

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

  #CALC SCREE
  sc <- scree(x, y)

  #CALCULATE MST MEASURES
  mst <- gen_mst(sc$del, sc$weights)
  if("stringy" %in% scags){
    stringy <- sc_stringy.mst(mst)
  }
  if("striated" %in% scags){
    striated <- sc_striated.mst(mst, sc)
  }
  if("clumpy" %in% scags){
    clumpy <- sc_clumpy.mst(mst, sc)
  }
  if("sparse" %in% scags){
    sparse <- sc_sparse.mst(mst, sc)
  }
  if("skewed" %in% scags){
    skewed <- sc_skewed.mst(mst, sc)

  }
  if("outlying" %in% scags){
    outlying <- sc_outlying.mst(mst, sc)
  }
  #CALCULATE ALPHA HULL MEASURES
  chull <- gen_conv_hull(sc$del)
  ahull <- gen_alpha_hull(sc$del, sc$alpha)
  if("convex" %in% scags){
    convex <- sc_convex.hull(chull,ahull)
  }
  if("skinny" %in% scags){
    skinny <- sc_skinny.hull(ahull)
  }

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




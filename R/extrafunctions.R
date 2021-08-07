#' Compute scagnostics on all possible scatter plots on data
#'
#' @examples
#'
#' @export
pairwise_scags <- function(all_data, scags=NULL, groups=NULL){
  all_combs <- expand.grid(colnames(all_data),colnames(all_data))

  }

#' Compute scagnostics on all possible scatter plots on data
#'
#' @examples
#'
#' @export
calc_scags <- function(x, y scags=NULL){
  original_scree <- scree(x,y)
  original_mst <- gen_mst(x,y)
  #return row of observations
}

fakefunc <- function(gridrow,origdata){
  print(origdata[[gridrow]])
  #calcs <- tibble(total=sum(datax, datay),
  #        means_x = mean(datax),
  #        mean_y = mean(datay))
  #return(calcs)
}

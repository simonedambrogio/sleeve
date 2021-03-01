#' Create the set of Boolean Classification Matrices (BCM)
#'
#' Given a set of raters (nrt), items (ni) of a scale and their
#' investigated symptoms/attributes (na), create_BCMM_list() returns
#' a data list containg the Boolean Classification Matrix (BCM) of each
#' rater.
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#'     It can be a `as_BCM` object.
#' @return A list containing the BCM of each rater.
#' @export
#' @seealso `as_BCM`, `toyEx`
#' @examples
#' data(toyEx)
#'
#' create_BCM_list(3,12,8,data=toyEx)
create_BCM_list<-function(nrt, ni, na, data){
  first_item <- seq(1, nrt*ni, by=ni)
  last_item <- c(first_item[2:length(first_item)]-1, nrt*ni)
  BCM_list <- list(matrix(NA, ni, na))
  for(i in 1:nrt){
    if(nrt*ni!=nrow(data)){
      stop("Please, check that data's number of rows is identical to
           nrt*ni")
    } else{
      BCM_list[[i]] <- data[first_item[i]:last_item[i], 2:9]
    }
  }
  return(BCM_list)
}

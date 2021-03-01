#' Belonging Threshold (BT) estimation-weakly informative prior
#'
#' This function estimates the Belonging Threshold for each rater
#' in case of a weakly informative prior
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Data frame containing the original attribution..
#' @param NC Number of chains. Default is equal to four.
#' @param NI Number of iteration per chain. Default is equal to five.
#' @param sigma Standard Deviation Raters' BMM prior
#' @param uniform LOGICAL. If TRUE, the weakly informative prior will follow a uniform distribution.
#' @return A list containing the following elements:
#' \itemize{
#'   \item WIP - the vector containing the values of the weakly informative prior.
#'   \item BT - the vector containing the Belonging Threshold for each rater.
#'   \item BadF - the vector containing values reflecting the difference between posterior predictive and observed data
#' }
#' @export
#' @seealso `toyEx`
#' @examples
#' data(toyEx)
#'
#' ## Estimate BTs for three rater with a weakly informative prior
#' BT_WIP(3,12,8,data=toyEx,sigma=.5,uniform=TRUE)
BT_WIP<-function(nrt,ni,na,data,NC=4,NI=5,sigma=.5,uniform=FALSE){
  N<-ni*na
  df<-create_BCM_list(nrt,ni,na,data)
  mybiglist_WIP<-list()

  BCM<- as.numeric(unlist(df))
  BCM_raters<-NULL
  for (j in 1:length(df)) {
    BCM_raters<-cbind(BCM_raters,as.numeric(unlist(df[j])))

  }

if(uniform==TRUE){
  US_BMM_wi<-rep(.5,nrow(BCM_raters))
}else{
    US_BMM_wi_temp <- rowSums(BCM_raters)
  US_BMM_wi <- US_BMM_wi_temp
  US_BMM_wi_values <- seq((1/(nrt+2)), 1-(1/(nrt+2)), by= (1/(nrt+2))) #values that the prior can have
  Vs <- seq(0,nrt, by=1) # n of values
  US_BMM_wi <- US_BMM_wi_temp
  for (i in 1:(length(Vs)+1)){
    for(k in 1:length(US_BMM_wi_temp)){
      if (US_BMM_wi[k] == i-1) {
        US_BMM_wi[k] = US_BMM_wi_values[i]
      }
    }
  }

}

  weakly.informative.prior<-US_BMM_wi



  sigma <- rep(sigma, N*nrt) # NEW Standard Deviation Raters' BMM prior


  dataList <- list(N=N, nrt=nrt, BCM=BCM, PUS_BMM=weakly.informative.prior, sigma=sigma) # data for Stan model

  BT.model.wi <- rstan::sampling(stanmodels$BT_estimation,
                                 data=dataList, chains=NC, iter=NI) # run Stan model

  BT.post.wi <- rstan::extract(BT.model.wi,"BT",permuted=T)$BT   # BT posterior distribution


  BT.e.wi <- round(apply(BT.post.wi, 2, mean,trim=.25),3)  # BT estimated value apply(x, 2, function(x) mean(x, trim = .2))


  BMM.e.wi <- rstan::extract(BT.model.wi,"R_BMM")$R_BMM # BMM posterior distribution

  post.pred.wi <- NULL
  for(k in 1:nrt){
    start <- ((k-1)*N)+1
    end <- k*N
    if (k==1){
      post.pred.wi <- c(ifelse(apply(BMM.e.wi[, start:end ], 2, mean) > BMM.e.wi[k],1,0))
    }else{
      post.pred.wi <- c(post.pred.wi,c(ifelse(apply(BMM.e.wi[, start:end ], 2, mean) > BMM.e.wi[k],1,0)))
    }
  }

  badfitted.wi <- xor(post.pred.wi, BCM) # difference between posterior predictive and observed data


  tmp <- list(WIP=weakly.informative.prior,
              BT=BT.e.wi, BadF=badfitted.wi)
  mybiglist_WIP <- tmp




  return(mybiglist_WIP)

}


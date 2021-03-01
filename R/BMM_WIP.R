#' Belonging Measure Matrix (BMM) estimation-weakly informative prior
#'
#' This function estimates the Belonging Measure Matrix, starting from
#' the Belonging Threshold in case of a weakly informative prior
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#' @param NC Number of chains. Default is equal to four.
#' @param NI Number of iteration per chain. Default is equal to five
#' @param sigma Standard Deviation Raters' BMM prior
#' @param sigma.bad Standard Deviation Raters' BMM prior based on the difference between posterior predictive and observed data
#' @param uniform LOGICAL. If TRUE, the weakly informative prior will follow a uniform distribution.
#' @return A list containing the following elements:
#' \itemize{
#'   \item WIP - the vector containing the values of the weakly informative prior.
#'   \item BT - the vector containing the Belonging Threshold for each rater.
#'   \item BCM.WI - BCM estimated from BMM and BT
#'   \item BadF - the vector containing values reflecting the difference between posterior predictive and observed data
#' }
#' @export
#' @seealso `toyEx`
#' @examples
#' data(toyEx)
#'
#' ## Estimate BTs for three rater with a weakly informative prior
#' BMM_WIP(3,12,8,data=toyEx,sigma=.5,uniform=TRUE,sigma.bad=1)
BMM_WIP<-function(nrt,ni,na,data,NC=4,NI=5,sigma=.5,sigma.bad=1,uniform=FALSE){
  N<-ni*na
  df<-create_BCM_list(nrt,ni,na,data)
  sv<-BT_WIP(3,12,8,data,sigma=sigma,uniform=uniform) #BT values
  myBMM_WIP<-list()

  BCM<- as.numeric(unlist(df))
  BCM_raters<-NULL
  for (j in 1:length(df)) {
    BCM_raters<-cbind(BCM_raters,as.numeric(unlist(df[j])))

  }

  BT.e.wi=sv[[2]]
  init_values <- NULL
  for(k in 1:nrt){
    if (k==1){
      init_values <- c(ifelse(BCM_raters[,k]==0,BT.e.wi[k]/2, (BT.e.wi[k] + (1-BT.e.wi[k])/2)))
    }else{
      init_values <- c(init_values,c(ifelse(BCM_raters[,k]==0,BT.e.wi[k]/2, (BT.e.wi[k] + (1-BT.e.wi[k])/2))))
    }
  }

  init_list <- list()
  Q <- NC #per il numero di catene
  for (s in 1:Q){
    c.name <- paste("c", s, sep="")
    init_list[[c.name]] <- list(theta_r=init_values)
  }

  sigma <- rep(sigma,N*nrt)  # Standard Deviation Raters' BMM prior
  badfitted.wi=sv[[3]]
  sigma[badfitted.wi] <- sigma.bad
  dim(BT.e.wi) <- nrt
  weakly.informative.prior=sv[[1]]

  dataList <- list(N=N, nrt=nrt, BCM=BCM, PUS_BMM= weakly.informative.prior, BT=BT.e.wi, sigma=sigma) # data for Stan model


  BMM.model.wi <- rstan::sampling(stanmodels$BMM_estimation, data=dataList, chains=NC, iter=NI, init = init_list) # run Stan model


  BMM.post.wi <- rstan::extract(BMM.model.wi,"R_BMM",permuted=T)$R_BMM # BMM posterior distribution



  BMM.e.wi    <- apply(BMM.post.wi,2,mean,trim=.25)           # BMM estimated value
  BCM.e.wi.5  <- ifelse(BMM.e.wi>.5,1,0)             # BCM estimated from BMM and standard BT (.5)




  tmp <- list(WIP=weakly.informative.prior,BT=BT.e.wi,
              BCM.WI=BCM.e.wi.5,BadF=badfitted.wi)
  myBMM_WIP <- tmp



  return(myBMM_WIP)

}

#' Belonging Measure Matrix (BMM) estimation- informative prior
#'
#' This function estimates the Belonging Measure Matrix, starting from
#' the Belonging Threshold in case of an informative prior
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution..
#' @param NC Number of chains. Default is equal to four.
#' @param NI Number of iteration per chain. Default is equal to five.
#' @param data_prior A dataframe containing data representing the informative prior
#' @param sigma Standard Deviation Raters' BMM prior
#' @param sigma.bad Standard Deviation Raters' BMM prior based on the difference between posterior predictive and observed data.
#' @param correct.BT LOGICAL. If TRUE, it corrects the BT estimated value.
#' @return A list containing the following elements:
#' \itemize{
#'   \item IP - the vector containing the values of the weakly informative prior.
#'   \item BT - the vector containing the Belonging Threshold for each rater.
#'   \item BCM.I - BCM estimated from BMM and BT
#'   \item BadF - the vector containing values reflecting the difference between posterior predictive and observed data
#' }
#' @export
#' @seealso `toyEx`
#' @examples
#' data(expertBMM)
#'
#' ## Estimate BTs for three rater with a weakly informative prior
#' expertBMM<-expertBMM[,2:9]
#' BMM_IP(3,12,8,data=toyEx,data_prior=expertBMM,sigma=.5,correct.BT=FALSE,sigma.bad=1)
BMM_IP<-function(nrt,ni,na,data,NC=4,NI=5,data_prior,sigma=.5,correct.BT=FALSE,sigma.bad=1){
  N<-ni*na
  df<-create_BCM_list(nrt,ni,na,data)
  sv<-BT_IP(nrt=nrt,ni=ni,na=na,data=data,data_prior=data_prior,sigma=sigma,correct.BT=correct.BT) #BT values
  myBMM_IP<-list()

  BCM<- as.numeric(unlist(df))
  BCM_raters<-NULL
  for (j in 1:length(df)) {
    BCM_raters<-cbind(BCM_raters,as.numeric(unlist(df[j])))

  }

  BT.e.i=sv[[2]]
  init_values <- NULL
  for(k in 1:nrt){
    if (k==1){
      init_values <- c(ifelse(BCM_raters[,k]==0,BT.e.i[k]/2, (BT.e.i[k] + (1-BT.e.i[k])/2)))
    }else{
      init_values <- c(init_values,c(ifelse(BCM_raters[,k]==0,BT.e.i[k]/2, (BT.e.i[k] + (1-BT.e.i[k])/2))))
    }
  }

  init_list <- list()
  Q <- NC #per il numero di catene
  for (s in 1:Q){
    c.name <- paste("c", s, sep="")
    init_list[[c.name]] <- list(theta_r=init_values)
  }

  sigma <- rep(sigma,N*nrt)  # Standard Deviation Raters' BMM prior
  badfitted.i=sv[[3]]
  sigma[badfitted.i] <- sigma.bad
  dim(BT.e.i) <- nrt
  informative.prior=sv[[1]]


    dataList <- list(N=N, nrt=nrt, BCM=BCM, PUS_BMM= informative.prior, BT=BT.e.i, sigma=sigma) # data for Stan model


    BMM.model.i <- rstan::sampling(stanmodels$BMM_estimation, data=dataList, chains=NC, iter=NI, init = init_list) # run Stan model


    BMM.post.i <- rstan::extract(BMM.model.i,permuted=T,"R_BMM")$R_BMM # BMM posterior distribution


    BMM.e.i    <- apply(BMM.post.i,2,mean,trim=.25)           # BMM estimated value
    BCM.e.i.5  <- ifelse(BMM.e.i>.5,1,0)             # BCM estimated from BMM and standard BT (.5)


    tmp <- list(IP=informative.prior,BT=BT.e.i,
                BCM.I=BCM.e.i.5,BadF=badfitted.i)
    myBMM_IP<- tmp




  return(myBMM_IP)

}

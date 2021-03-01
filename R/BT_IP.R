#' Belonging Threshold (BT) estimation- informative prior
#'
#' This function estimates the Belonging Threshold for each rater
#' in case of an informative prior
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Data frame containing the original attribution.
#' @param NC Number of chains. Default is equal to four.
#' @param NI Number of iteration per chain. Default is equal to five.
#' @param data_prior A data frame containing data representing the informative prior
#' @param sigma Standard Deviation Raters' BMM prior
#' @param correct.BT LOGICAL. If TRUE, it corrects the BT estimated value.
#' @return A list containing the following elements:
#' \itemize{
#'   \item IP - the vector containing the values of the informative prior.
#'   \item BT - the vector containing the Belonging Threshold for each rater.
#'   \item BadF - the vector containing values reflecting the difference between posterior predictive and observed data
#' }
#' @export
#' @seealso `toyEx`
#' @examples
#' data(expertBMM)
#'
#' ## Estimate BTs for three rater with a weakly informative prior
#' expertBMM<-expertBMM[,2:9]
#' BT_IP(3,12,8,data=toyEx,data_prior=expertBMM,sigma=.5,correct.BT=FALSE)
BT_IP<-function(nrt,ni,na,data,NC=4,NI=5,data_prior,sigma=.5,correct.BT=FALSE){
  N<-ni*na
  df<-create_BCM_list(nrt,ni,na,data)
  mybiglist_IP<-list()
  BCM<- as.numeric(unlist(df))

    informative.prior <- c(as.matrix(data_prior)) # Mean Raters' BMM prior from US_BMM
    sigma <- rep(sigma, N*nrt)      # NEW Standard Deviation Raters' BMM prior

    dataList <- list(N=N, nrt=nrt, BCM=BCM, PUS_BMM=informative.prior, sigma=sigma) #  data for Stan model

    BT.model.i <- rstan::sampling(stanmodels$BT_estimation,
                       data=dataList, chains=NC, iter=NI) # run Stan model


    BT.post.i <- rstan::extract(BT.model.i,permuted=T,"BT")$BT      # BT posterior distribution



    BT.e.i <- round(apply(BT.post.i, 2, mean,trim=.25),3)  # BT estimated value

   if(correct.BT==TRUE){
         for(k in 1:length(BT.e.i)){
      if(BT.e.i[k]<=0.2) {BT.e.i [k] <- .1 #correction to possible BT
      }
      if(BT.e.i[k]> 0.2 & BT.e.i[k]<=0.4) {BT.e.i [k] <- .3
      }
      if(BT.e.i[k]> 0.4 & BT.e.i[k]<=0.6) {BT.e.i [k] <- .5
      }
      if(BT.e.i[k]> 0.6 & BT.e.i[k]<=0.8) {BT.e.i [k] <- .7
      }
      if(BT.e.i[k]> 0.8) {BT.e.i [k] <- .9
      }
    }
   }



    BMM.e.i <- rstan::extract(BT.model.i,"R_BMM")$R_BMM # BMM posterior distribution

    post.pred.i <- NULL
    for(k in 1:nrt){
      start <- ((k-1)*N)+1
      end <- k*N
      if (k==1){
        post.pred.i <- c(ifelse(apply(BMM.e.i[, start:end ], 2, mean) > BMM.e.i[k],1,0))
      }else{
        post.pred.i <- c(post.pred.i,c(ifelse(apply(BMM.e.i[, start:end ], 2, mean) > BMM.e.i[k],1,0)))
      }
    }

    badfitted.i <- xor(post.pred.i,BCM) # difference between posterior predictive and observed data

    tmp <- list(I=informative.prior,
                BT=BT.e.i, BadF=badfitted.i)
    mybiglist_IP <- tmp



  return(mybiglist_IP)

}

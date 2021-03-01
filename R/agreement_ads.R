#' Agreement's coefficients adjusted
#'
#' This function calculates different agreement's coefficients
#' between n-couples of raters
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Dataframe containing the original attribution.
#' @param NC number of chains
#' @param NI number of iterations
#' @param data_prior dataframe or vector containing the prior
#' @param prior weakly informatiove or informative prior
#' @param sigma Standard Deviation Raters' BMM prior
#' @param sigma.bad Standard Deviation Raters' BMM prior based on the difference between posterior predictive and observed data
#' @param uniform LOGICAL. If TRUE, the weakly informative prior will follow a uniform distribution.
#' @param correct.BT LOGICAL. If TRUE, it corrects the BT estimated value.
#' @param index the inter-rater agreement index to be estimated
#' @return A list containing the following elements:
#' \itemize{
#'   \item All_BCM - a data frame containing the raw data and the estimated ones, for each rater
#'   \item Results - a data frame containing the following columns:
#'   \enumerate{
#'      \item The Couple of rater compared.
#'      \item The mean BCM of Rater 1.
#'      \item The mean BCM of Rater 2.
#'      \item The BT of Rater 1.
#'      \item The BT of Rater 2.
#'      \item The agreement index before the adjustment.
#'      \item The agreement index after the adjustment.
#'
#'   }
#' }
#' @export
#' @seealso `toyEx`
#' @examples
#' data(toyEx)
#'
#' agreement_ads(nrt=3,ni=12,na=8,data=toyEx,NC=4,NI=5,prior="WI",index = "SPA",
#' sigma = .3,sigma.bad=1,uniform=TRUE)

agreement_ads<-function(nrt,ni,na,data,NC=4,NI=5,
                        data_prior=NULL,prior=NULL,index=NULL,sigma=.5,
                        sigma.bad=1,uniform=FALSE,correct.BT=FALSE){
  N<-ni*na

  df<-create_BCM_list(nrt,ni,na,data)
  nof_row<-combn(1:nrt,2)
  # namesCol<-c("Couple","BCM_R1","BCM_R2","SPA",
  #             "Kappa","BT_wi_R1","BT_wi_R2",
  #             "BT_i_R1","BT_i_R2","SPA_adj_wi",
  #             "Kappa_adj_wi","SPA_adj_i","Kappa_adj_i")

  namesCol<-c("Couple","BCM_R1","BCM_R2",
              "BT_wi_R1","BT_wi_R2", # fin qu se WI
              "BT_i_R1","BT_i_R2","index_old","index_new")   # altrimenti questi due se I

  Results<-data.frame(matrix(NA,nrow = ncol(nof_row),ncol = length(namesCol)))


  if(prior=="WI"){
    a<-BMM_WIP(nrt=nrt,ni=ni,na=na,data=data,NC=NC,NI=NI,sigma=sigma,sigma.bad=sigma.bad,uniform=uniform)

    BCM<- as.numeric(unlist(df))
    BCM_raters<-NULL
    for (m in 1:length(df)) {
      BCM_raters<-cbind(BCM_raters,as.numeric(unlist(df[m])))

    }
    colnames(BCM_raters)=paste("BCM","Rater",1:ncol(BCM_raters),sep="_")

    BCM_adj<- a[[3]]
    BCM_raters_adj<-matrix(BCM_adj,ncol=nrt)
    colnames(BCM_raters_adj)=paste("BCM","adj","Rater",1:ncol(BCM_raters_adj),sep="_")
    sleeve<-a[[2]]

    combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- BCM_raters[,couple[1]]
      BCM.r.2 <- BCM_raters[,couple[2]]
      sub_sleeve<-sleeve[c(couple[1],couple[2])]
      BCM.r.1.e.wi.5 <- BCM_raters_adj[,couple[1]]
      BCM.r.2.e.wi.5 <- BCM_raters_adj[,couple[2]]

      if(index=="SPA"){
        SPA.wi.5  <- sum(BCM.r.1.e.wi.5==BCM.r.2.e.wi.5)/N # Estimated (wi) simple Percentage Agreement whit standard threshold (.5)
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-sub_sleeve[1]
        Results[i,5]<-sub_sleeve[2]
        Results[i,6]<-NA
        Results[i,7]<-NA
        Results[i,8]<-sum(BCM.r.1==BCM.r.2)/N
        Results[i,9]<-SPA.wi.5
        namesCol[c(8,9)]<-c("SPA","SPA_adj_wi")


      } else if(index=="Kappa"){
        CK.wi.5<-irr::kappa2(cbind(BCM.r.1.e.wi.5,BCM.r.2.e.wi.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-sub_sleeve[1]
        Results[i,5]<-sub_sleeve[2]
        Results[i,6]<-NA
        Results[i,7]<-NA
        Results[i,8]<- irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,9]<-CK.wi.5
        namesCol[c(8,9)]<-c("Kappa","Kappa_adj_wi")

      } else if(index=="Maxwell"){
        mx.wi.5<-irr::maxwell(cbind(BCM.r.1.e.wi.5,BCM.r.2.e.wi.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-sub_sleeve[1]
        Results[i,5]<-sub_sleeve[2]
        Results[i,6]<-NA
        Results[i,7]<-NA
        Results[i,8]<- irr::maxwell(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,9]<-mx.wi.5
        namesCol[c(8,9)]<-c("Maxwell","Maxwell_adj_wi")

      } else {
        smx.wi.5<-irr::stuart.maxwell.mh(cbind(BCM.r.1.e.wi.5,BCM.r.2.e.wi.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-sub_sleeve[1]
        Results[i,5]<-sub_sleeve[2]
        Results[i,6]<-NA
        Results[i,7]<-NA
        Results[i,8]<- irr::stuart.maxwell.mh(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,9]<-smx.wi.5
        namesCol[c(8,9)]<-c("Stuart-Maxwell","Stuart-Maxwell_adj_wi")

      }

    }
  } else {
    a<-BMM_IP(nrt=nrt,ni=ni,na=na,data=data,NC=NC,NI=NI,data_prior = data_prior,sigma=sigma,correct.BT=correct.BT,sigma.bad=sigma.bad)

    BCM<- as.numeric(unlist(df))
    BCM_raters<-NULL
    for (m in 1:length(df)) {
      BCM_raters<-cbind(BCM_raters,as.numeric(unlist(df[m])))

    }
    colnames(BCM_raters)=paste("BCM","Rater",1:ncol(BCM_raters),sep="_")
    BCM_adj<- a[[3]]
    BCM_raters_adj<-matrix(BCM_adj,ncol=nrt)
    colnames(BCM_raters_adj)=paste("BCM","adj","Rater",1:ncol(BCM_raters_adj),sep="_")
    sleeve<-a[[2]]

    combinations<-combn(1:nrt,2)
    for (i in 1:ncol(combinations)) {
      couple<-combinations[,i]
      BCM.r.1 <- BCM_raters[,couple[1]]
      BCM.r.2 <- BCM_raters[,couple[2]]
      sub_sleeve<-sleeve[c(couple[1],couple[2])]
      BCM.r.1.e.i.5 <- BCM_raters_adj[,couple[1]]
      BCM.r.2.e.i.5 <- BCM_raters_adj[,couple[2]]


      if(index=="SPA"){
        SPA.i.5  <- sum(BCM.r.1.e.i.5==BCM.r.2.e.i.5)/N # Estimated (wi) simple Percentage Agreement whit standard threshold (.5)
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-NA
        Results[i,5]<-NA
        Results[i,6]<-sub_sleeve[1]
        Results[i,7]<-sub_sleeve[2]
        Results[i,8]<-sum(BCM.r.1==BCM.r.2)/N
        Results[i,9]<-SPA.i.5
        namesCol[c(8,9)]<-c("SPA","SPA_adj_i")

      } else if(index=="Kappa"){
        CK.i.5<-irr::kappa2(cbind(BCM.r.1.e.i.5,BCM.r.2.e.i.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-NA
        Results[i,5]<-NA
        Results[i,6]<-sub_sleeve[1]
        Results[i,7]<-sub_sleeve[2]
        Results[i,8]<- irr::kappa2(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,9]<-CK.i.5
        namesCol[c(8,9)]<-c("Kappa","Kappa_adj_i")

      } else if(index=="Maxwell"){
        mx.i.5<-irr::maxwell(cbind(BCM.r.1.e.i.5,BCM.r.2.e.i.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-NA
        Results[i,5]<-NA
        Results[i,6]<-sub_sleeve[1]
        Results[i,7]<-sub_sleeve[2]
        Results[i,8]<- irr::maxwell(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,9]<-mx.i.5
        namesCol[c(8,9)]<-c("Maxwell","Maxwell_adj_i")

      } else {
        smx.i.5<-irr::stuart.maxwell.mh(cbind(BCM.r.1.e.i.5,BCM.r.2.e.i.5))$value
        Results[i,1]<-paste(as.character(couple),collapse="-")
        Results[i,2]<-mean(BCM.r.1)
        Results[i,3]<-mean(BCM.r.2)
        Results[i,4]<-NA
        Results[i,5]<-NA
        Results[i,6]<-sub_sleeve[1]
        Results[i,7]<-sub_sleeve[2]
        Results[i,8]<- irr::stuart.maxwell.mh(cbind(BCM.r.1,BCM.r.2))$value
        Results[i,9]<-smx.i.5
        namesCol[c(8,9)]<-c("Stuart-Maxwell","Stuart-Maxwell_adj_wi")
      }
    }

  }
  Results<-data.frame(Results)
  colnames(Results)<-namesCol
  Results<-Filter(function(x) !all(is.na(x)),Results)
  Results[2:ncol(Results)]<-round(Results[2:ncol(Results)],3)

  myResList<-list(All_BCM=cbind(BCM_raters,BCM_raters_adj),
                  Results=Results)
  return(myResList)
}

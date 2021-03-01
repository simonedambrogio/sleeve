#' Agreement's coefficients adjusted
#'
#' This function calculates different agreement's coefficients
#' between n-couples of raters
#'
#' @param nrt Number of raters.
#' @param ni Number of items.
#' @param na Number of symptoms/attributes investigated by items.
#' @param data Data frame containing the original attribution.
#' @param NC number of chains
#' @param NI number of iterations
#' @param data_prior data frame or vector containing the prior
#' @param prior weakly informatiove or informative prior
#' @param sigma Standard Deviation Raters' BMM prior
#' @param sigma.bad Standard Deviation Raters' BMM prior based on the difference between posterior predictive and observed data
#' @param uniform LOGICAL. If TRUE, the weakly informative prior will follow a uniform distribution.
#' @param correct.BT LOGICAL. If TRUE, it corrects the BT estimated value.
#' @param index the inter-rater agreement index to be estimated
#' @param model model
#' @param s.levels s_level
#' @param type type
#' @param unit unit
#' @param exact logical
#' @param detail details
#' @param correct logical
#' @return A list containing the following elements:
#' \itemize{
#'   \item All_BCM - a data frame containing the raw data and the estimated ones, for each rater
#'   \item Results - a data frame containing the following columns:
#'   \enumerate{
#'      \item The mean BCM of each Rater.
#'      \item The BT of Rater CM of each Rater.
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
#' ## Compute SPA only between rater 1 and 3
#' multi_rater_adj(3,12,8,data=toyEx,prior="WI",index="kappa.light",
#' sigma=.5,sigma.bad=1,uniform=TRUE)

multi_rater_adj<-function(nrt,ni,na,data,NC=4,NI=5,sigma=.5,
                        data_prior=NULL,prior,index,model=NULL,
                        s.levels=NULL,type=NULL,unit=NULL,
                        exact=NULL,detail = NULL,correct=NULL,
                        sigma.bad=1,uniform=FALSE,correct.BT=FALSE){
  N<-ni*na
  df<-create_BCM_list(nrt,ni,na,data)
  # namesCol<-c("Couple","BCM_R1","BCM_R2","SPA",
  #             "Kappa","BT_wi_R1","BT_wi_R2",
  #             "BT_i_R1","BT_i_R2","SPA_adj_wi",
  #             "Kappa_adj_wi","SPA_adj_i","Kappa_adj_i")

  # namesCol<-c("Couple","BCM_R1","BCM_R2",
  #             "BT_wi_R1","BT_wi_R2", # fin qu se WI
  #             "BT_i_R1","BT_i_R2","index_old","index_new")   # altrimenti questi due se I
  #
  Results<-matrix(NA,nrow=1,ncol=0)

  if(prior=="WI"){
    a<-BMM_WIP(nrt=nrt,ni=ni,na=na,data=data,NC=NC,NI=NI,
               sigma=sigma,sigma.bad=sigma.bad,uniform=uniform)


    datapre<-NULL
    for (m in 1:length(df)) {
      datapre<-cbind(datapre,as.numeric(unlist(df[m])))

    }
    colnames(datapre)<-paste("BCM","Rater",1:nrt,sep="_")

    BCM_adj<- a[[3]]
    datapost<-matrix(BCM_adj,ncol=nrt)
    colnames(datapost)=paste("BCM","adj","Rater",1:ncol(datapost),sep="_")
    datapreBCM<-colMeans(datapre)
    datapostBCM<-a[[2]]
    colnames(datapre)<-paste("BT_wi","Rater",1:nrt,sep="_")



    if(index=="Finn"){

      finnPre<-irr::finn(datapre,s.levels = s.levels,model=model)$value
      finnPost<-irr::finn(datapost,s.levels = s.levels,model=model)$value
      Results<-cbind(datapreBCM,datapostBCM,finnPre,finnPost)

    } else if(index=="kappa.fleiss"){

      kappa_f_Pre<-irr::kappam.fleiss(datapre,exact=exact,detail = detail)$value
      kappa_f_Post<-irr::kappam.fleiss(datapost,exact=exact,detail = detail)$value
      Results<-cbind(datapreBCM,datapostBCM,kappa_f_Pre,kappa_f_Post)

    } else if(index=="kappa.light"){

      kappa_l_Pre<-irr::kappam.light(datapre)$value
      kappa_l_Post<-irr::kappam.light(datapost)$value
      Results<-cbind(datapreBCM,datapostBCM,kappa_l_Pre,kappa_l_Post)

    } else if(index=="kendall"){

      kendallPre<-irr::kendall(datapre,correct=correct)$value
      kendallPost<-irr::kendall(datapost,correct=correct)$value
      Results<-cbind(datapreBCM,datapostBCM,kendallPre,kendallPost)

    } else {

      ICCPre<-irr::icc(datapre,typ=type,model=model,unit=unit)$value
      ICCPost<-irr::icc(datapost,typ=type,model=model)$value
      Results<-cbind(datapreBCM,datapostBCM,ICCPre,ICCPost)
    }


  } else {
    a<-BMM_IP(nrt=nrt,ni=ni,na=na,data=data,NC=NC,NI=NI,data_prior = data_prior,
              sigma=sigma,sigma.bad=sigma.bad,correct.BT=correct.BT)


    datapre<-NULL
    for (m in 1:length(df)) {
      datapre<-cbind(datapre,as.numeric(unlist(df[m])))

    }
    colnames(datapre)<-paste("BCM","Rater",1:nrt,sep="_")

    BCM_adj<- a[[3]]
    datapost<-matrix(BCM_adj,ncol=nrt)

    datapreBCM<-colMeans(datapre)
    datapostBCM<-a[[2]]
    colnames(datapre)<-paste("BT_i","Rater",1:nrt,sep="_")




    if(index=="Finn"){

      finnPre<-irr::finn(datapre,s.levels = s.levels,model=model)$value
      finnPost<-irr::finn(datapost,s.levels = s.levels,model=model)$value
      Results<-cbind(datapreBCM,datapostBCM,finnPre,finnPost)

    } else if(index=="kappa.fleiss"){

      kappa_f_Pre<-irr::kappam.fleiss(datapre,exact=exact,detail = detail)$value
      kappa_f_Post<-irr::kappam.fleiss(datapost,exact=exact,detail = detail)$value
      Results<-cbind(datapreBCM,datapostBCM,kappa_f_Pre,kappa_f_Post)

    } else if(index=="kappa.light"){

      kappa_l_Pre<-irr::kappam.light(datapre)$value
      kappa_l_Post<-irr::kappam.light(datapost)$value
      Results<-cbind(datapreBCM,datapostBCM,kappa_l_Pre,kappa_l_Post)

    } else if(index=="kendall"){

      kendallPre<-irr::kendall(datapre,correct=correct)$value
      kendallPost<-irr::kendall(datapost,correct=correct)$value
      Results<-cbind(datapreBCM,datapostBCM,kendallPre,kendallPost)

    } else {

      ICCPre<-irr::icc(datapre,typ=type,model=model,unit=unit)$value
      ICCPost<-irr::icc(datapost,typ=type,model=model)$value
      Results<-cbind(datapreBCM,datapostBCM,ICCPre,ICCPost)
    }


  }
  myResList<-list(All_BCM=cbind(datapre,datapost),
                  Results=Results)
  return(myResList)
}

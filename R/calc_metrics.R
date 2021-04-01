#' Calculate vector of error values
#' @describeIn calculate_vector_errors absolute error
#' 
#' @param y a vector of true values
#' @param y_hat a vector of predicted values
#' 
#' @return a vector of error values
ae <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))
  
  abs(y - y_hat)
}

#' @describeIn calculate_vector_errors squared error
se <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat),
            is.numeric(y),
            is.numeric(y_hat))
  
  (y - y_hat) ^ 2
}

#' Calculate error metrics
#' @describeIn calculate_errors mean squared error
#' @inheritParams ae
#' @param na.rm boolean indicating whether NAs should be removed. 
#' Default is TRUE
#' 
#' @return one error value
mse <- function(y, y_hat, na.rm=TRUE) mean(se(y, y_hat), na.rm = na.rm)

#' @describeIn calculate_errors mean squared error
rmse <- function(y, y_hat, na.rm=TRUE) sqrt(mse(y, y_hat, na.rm=na.rm))

#' @describeIn calculate_errors mean absolute error
mae <- function(y, y_hat, na.rm=TRUE) mean(ae(y, y_hat), na.rm = na.rm)

#' @describeIn calculate_errors normalized mean absolute error 
#' @param statFUN summary statistic to use for normalization. Default is 
#' median for nmae and mean for nmse.
#' 
#' @param y_train a vector of training values
nmae <- function(y, y_hat, y_train=NULL, statFUN=stats::median, na.rm=TRUE){
  sae <- sum(ae(y, y_hat), na.rm=na.rm)
  if(!is.null(y_train)) denom <- sum(abs(y - statFUN(y_train, na.rm=na.rm)), na.rm=na.rm)
  else denom <- sum(abs(y - statFUN(y, na.rm=na.rm)), na.rm=na.rm)
  sae/denom
} 

#' @describeIn calculate_errors normalized mean squared error
nmse <- function(y, y_hat, y_train=NULL, statFUN=mean, na.rm=TRUE){
  sse <- sum(se(y, y_hat), na.rm=na.rm)
  if(!is.null(y_train)) denom <- sum((y - statFUN(y_train, na.rm=na.rm))^2, na.rm=na.rm)
  else denom <- sum((y - statFUN(y, na.rm=na.rm))^2, na.rm=na.rm)
  sse/denom
} 

#' @describeIn calculate_errors normalized root mean squared error
nrmse <- function(y, y_hat, y_train=NULL, statFUN=mean, na.rm=TRUE){
  sqrt(nmse(y, y_hat, y_train=y_train, statFUN=statFUN, na.rm=na.rm))
}

#' Calculate regression metrics
#' 
#' Calculate MAE, RMSE and utility-based regression evaluation metrics
#' @param trues a vector of true values
#' @param preds a vector of predicted values
#' @param y_train a vector of training values
#' @param norm a Boolean indicating whether to calculate normalized
#' regression metrics
#' @param aeStatFUN a function to calculate a summary of y_train for 
#' absolute error normalization. Default is median
#' @param seStatFUN a function to calculate a summary of y_train for 
#' squared error normalization. Default is mean
#' @param util a Boolean indicating whether to calculate utility-based
#' regression metrics
#' @param util.parms a named list of parameters to use for calculating utility-based
#' regression metrics. Should contain slots
#' \itemize{
#' \item \code{rel} - the relevance function. Default is "auto" which uses method="range"
#' with utilparms$cf as coefficient
#' \item \code{thr} - Relevance threshold. Default is 1
#' \item \code{beta} - Beta for F-measure. Default is 1
#' }
#' 
#' @return a named vector of calculated metrics
#' 
#' @export
regMetrics <- function(trues, preds, y_train=NULL, 
                       norm=FALSE, aeStatFUN = stats::median, seStatFUN = mean,
                       util=FALSE, util.parms=NULL){
  
  if(length(trues)==0 | all(is.na(preds))){
    metrics <- c(mae=NA, rmse=NA, bias.sq = NA, var.sq = NA)
    if(norm) metrics <- c(metrics, nmae=NA, nmse=NA)
    if(util) metrics <- c(metrics, rmse.phi = NA, 
                          prec.u=NA, rec.u=NA, Fbeta.u=NA)
  }else{
    metrics <- c(mae = mae(trues, preds),
                 rmse = rmse(trues, preds),
                 bias.sq = mean(mean(preds) - trues)^2,
                 var.sq = mean(mean(preds) - preds)^2)
    
    if(norm)
      metrics <- c(metrics, nmae_tr = nmae(trues, preds, y_train, aeStatFUN),
                   nmse_tr = nmse(trues, preds, y_train, seStatFUN),
                   nrmse_tr = nrmse(trues, preds, y_train, seStatFUN),
                   nmae = nmae(trues, preds, NULL, aeStatFUN),
                   nmse = nmse(trues, preds, NULL, seStatFUN),
                   nrmse = nrmse(trues, preds, NULL, seStatFUN))
    
    if(util){
  
      rel <- util.parms$rel
      if(is.null(rel)) rel <- "auto"
      if (is.matrix(rel)) {
        pc <- uba::phi.control(y_train, method = "range", control.pts = rel)
      }
      else if (is.list(rel)) {
        pc <- rel
      }
      else if (rel == "auto") {
        cf <- util.parms$cf
        if(is.null(cf)) cf <- 1.5
        pc <- uba::phi.control(y_train, method = "extremes", coef = cf)
      }
      else {
        stop("Argument util.params$rel should be 'auto', a list returned by phi.control or a control.pts matrix.")
      }
    
      ls <- util.parms$loss.parms
      if(is.null(ls)) ls <- uba::loss.control(y_train) 
      
      thr <- util.parms$thr
      if(is.null(thr)) thr <- 0.9
      
      beta <- util.parms$beta
      if(is.null(beta)) beta <- 1
      
      u_new <- uba::util(preds, trues, pc, ls,
                         uba::util.control(umetric="P", event.thr=thr),
                         return.uv = TRUE)
      
      phi.trues <- UBL::phi(trues, control.parms = pc)
      phi.preds <- UBL::phi(preds, control.parms = pc)
      
      pr_frm <- data.frame(Utility=u_new)
      pr_frm["phiTrues"] <- phi.trues
      pr_frm["phiPreds"] <- phi.preds
      
      prec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiPreds>thr,]$phiPreds)
      rec.u <- sum(1+pr_frm[pr_frm$phiTrues>thr & pr_frm$phiPreds>thr,]$Utility)/sum(1+pr_frm[pr_frm$phiTrues>thr,]$phiTrues)
      Fbeta.u <- (1+beta) * prec.u * rec.u / ( beta^2 * prec.u + rec.u)
      
      rmse.phi <- sqrt(mean(phi.trues[phi.trues>thr]*(trues[phi.trues>thr]-preds[phi.trues>thr])^2))
      
      #clean results
      if(is.na(rmse.phi)) prec.u=rec.u=Fbeta.u=NA
      if(!is.na(rmse.phi)){
        if(is.na(prec.u)) prec.u <- 0
        if(is.na(rec.u)) rec.u <- 0
        if(is.na(Fbeta.u)) Fbeta.u <- 0
      }
      
      metrics <- c(metrics, rmse.phi = rmse.phi, 
                   prec.u=prec.u, rec.u=rec.u, Fbeta.u=Fbeta.u) #, aucpr=aucpr, aucroc=aucroc)
    
    }
  }

  metrics
}

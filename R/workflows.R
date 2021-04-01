#' Handling NAs in train and test
#' 
#' Discard columns/rows with too many NAs and 
#' then impute with central value.
#'
#' @param train training data set
#' @param test testing data set
#' @param nORp minimum percentage of NA values in a row/column for
#' that row/column to be discarded from training set
#'
#' @return list with an entry for the training and test sets
#' (in this order), both now with no NA values
#' 
#' @export
centralImputTsNAs <- function(train, test, nORp){
  
  if(anyNA(test[ , colnames(train)])){
    # fill in test columns with central value of train column
    for (cnm in colnames(train)) if (any(idx <- is.na(test[, cnm]))){
      test[idx, cnm] <- DMwR2::centralValue(train[, cnm])          
    } 
  }
  
  test
}

#' Remove columns with single unique values from train
#' 
#'
#' @param train training data set
#' @param time the name of the column in \code{train} and
#' \code{test} containing time-stamps
#' @param site_id the name of the column in \code{train} and
#' \code{test} containing location IDs
#'
#' @param tgt name of target variable
removeSingleValCol <- function(train, time, site_id, tgt){
  # discard columns with only one value
  discCols <- setdiff(which(sapply(train, 
                                   function(y) length(unique(y)) == 1)),
                      which(colnames(train) %in% c(tgt, site_id, time)))
  
  if(length(discCols) > 0) {
    if( (ncol(train) - length(discCols) - 3) > 0 ){
      warning(paste("Dropped", length(discCols), "columns from train with a single unique value. Keeping", ncol(train) - length(discCols),"variables."), call. = FALSE)
      train <- train[, -discCols]
    }else{
      warning(paste("Should have dropped", length(discCols), "columns from train with a single unique value. Keeping all because too few would be left."), call. = FALSE)
    }
  }
  
  train
}

#' Handling NAs in train and test
#' 
#' Discard columns/rows with too many NAs and 
#' then impute with central value.
#'
#' @param train training data set
#' @param nORp minimum percentage of NA values in a row/column for
#' that row/column to be discarded from training set
#'
#' @return list with an entry for the training and test sets
#' (in this order), both now with no NA values
#' 
#' @export
centralImputTrNAs <- function(train, nORp){
  
  # discard columns with too many NAs if there would be still predictors left
  discCols <- which(sapply(train, 
                           function(y) length(which(is.na(y)))/length(y)) > nORp)
  
  if(length(discCols) > 0) {
    if( (ncol(train) - length(discCols) - 3) > 0 ){
      warning(paste("Dropped", length(discCols), "columns from train due to NAs. Keeping", ncol(train) - length(discCols),"variables."), call. = FALSE)
      train <- train[, -discCols]
    }else{
      warning(paste("Should have dropped", length(discCols), "columns from train due to NAs. Keeping all because too few would be left."), call. = FALSE)
    }
  }
    
  # discard rows with too many NAs in train
  suppressWarnings( idxs <- DMwR2::manyNAs(train, nORp = nORp) )
  if(length(idxs)){
    warning(paste0("Dropped ", length(idxs), " rows from train due to NAs. Keeping ", nrow(train) - length(idxs)," (", round(100*(nrow(train) - length(idxs))/nrow(train)),"%) rows."), call. = FALSE)
    train <- train[-idxs, ]
  } 
  
  # fill in empty value in train
  if(anyNA(train)) train <- DMwR2::centralImputation(train)
  
  train
}



#' A simple learning and prediction workflow
#' 
#' @param train a data frame for training
#' @param test a data frame for testing
#' @param time the name of the column in \code{train} and
#' \code{test} containing time-stamps
#' @param site_id the name of the column in \code{train} and
#' \code{test} containing location IDs
#' @param form a formula describing the model to learn
#' @param model the name of the algorithm to use
#' @param handleNAs string indicating how to deal with NAs.
#' If "centralImput", training observations with at least 80\%
#' of non-NA columns, will have their NAs substituted by the mean
#' value and testing observatiosn will have their NAs filled in with
#' mean value regardless.
#' @param min_train a minimum number of observations that must be
#' left to train a model. If there are not enough observations, 
#' predictions will be \code{NA}. Default is 2.
#' @param nORp a maximum number or fraction of columns with missing
#' values above which a row will be removed from train before 
#' learning the model. Only works if \code{handleNAs} was
#' set to centralImputation. Default is 0.2.
#' @param ... other parameters to feed to \code{model}
#' 
#' @return a data frame containing time-stamps, location IDs,
#' true values and predicted values
#' 
#' @export
simple_workflow <- function(train, test, form, model="lm", 
  handleNAs=NULL, min_train=2, nORp = 0.2,
  time="time", site_id="site", ...){

    #----- ARGUMENT CHECK -------#
  
  dotargs <- list(...)
  
  assertthat::assert_that(min_train>=2, 
                          msg = "Cannot train model with less than 2 observations.")
  assertthat::assert_that(time %in% colnames(train),
                          time %in% colnames(test),
                          msg = "'time' not a column in data set")
  assertthat::assert_that(site_id %in% colnames(train),
                          site_id %in% colnames(test),
                          msg = "'site_id' not a column in data set")
  
  
  #----------------------------#
  
  # get true values
  trues <- responseValues(form, test)
  tgt <- as.character(form[[2]])
  
  if(nrow(train) >= (2/3)*nrow(test)){
    
    col.inds <- which(colnames(train) %in% c(time, site_id))
    # correct default mtry if model is ranger and there is no argument given
    if(model=="ranger" & !("mtry" %in% dotargs) & is.numeric(trues))
      dotargs$mtry <- max(floor(ncol(train[,-col.inds])/3), 1)
    
    # remove training columns with only one unique value
    train <- removeSingleValCol(train, tgt = tgt, site_id=site_id, time=time)
    # pre-process NAs
    if(!is.null(handleNAs)){
      if(handleNAs=="centralImput"){
        test <- centralImputTsNAs(train, test, nORp)
        train <- centralImputTrNAs(train, nORp)
      }else{
        stop("Provide function to handle NAs.")
      }
    }
  
    if(nrow(train) >= min_train){
      # check if columns need to be removed from train
      tr.col.inds <- which(colnames(train) %in% c(time, site_id))
      ts.col.inds <- which(colnames(test) %in% c(time, site_id))
      # removing offending column indices from train and test
      if(length(tr.col.inds)) train <- as.data.frame(train)[,-tr.col.inds]
      # train model
      m <- do.call(model, c(list(form, train), dotargs))
      # make predictions
      if(model=="ranger"){
        preds <- stats::predict(m, test[,-ts.col.inds])$predictions
      } else{
        preds <- stats::predict(m, test[,-ts.col.inds])
        if (is.numeric(train[[tgt]]) && !is.null(dim(preds)))
          preds <- preds[, 1]
      } 
      # prepare result object
      res <- data.frame(trues=trues, preds=preds)
    }else{
      warning("nrow(train)<min_train", call. = FALSE)
      res <- data.frame(trues=trues, preds=as.numeric(NA))  
    }
  }else{
    warning("nrow(train)< 2/3 * nrow(test)", call. = FALSE)
    res <- data.frame(trues=trues, preds=as.numeric(NA))
  }
  
  res
}

#' Evalute the results of a predictive workflow
#' 
#' Calculate evaluation metrics from the raw results of a workflow
#' @param wfRes a data frame (or list of data frames) containing the results of
#' a predictive workflow with columns \code{trues} and \code{preds} containing
#' the real and predicted values, respectively
#' @param eval.function the function to be used to calculate error metrics from \code{wfRes}
#' @param .keptTrain a Boolean indicating whether \code{.keepTrain} was
#' set to \code{TRUE} in calls to estimation methods. Only useful
#' if evaluation metrics need training data.
#' @param ... parameters to pass to \code{eval.function}
#'
#' @return The results (or a list of results) of \code{eval.function} applied to 
#' the data frame (or list of data frames) in \code{wfRes}
#' 
#' @export
evaluate <- function(wfRes,
                     eval.function = get("regressionMetrics", asNamespace("performanceEstimation")),
                     .keptTrain = TRUE,
                     ...){
  
  if(!.keptTrain){
    if(!("results" %in% names(wfRes))) 
      fold.res <- 
        dplyr::bind_rows(lapply(wfRes, function(x) 
        eval.function(trues=x$results$trues, 
                      preds=x$results$preds, ...)), .id="fold")
    else fold.res <- as.data.frame(t(eval.function(trues=wfRes$results$trues, 
                                     preds=wfRes$results$preds, ...))) 
  }else{
    if(!("results" %in% names(wfRes))) 
      fold.res <- dplyr::bind_rows(lapply(wfRes, function(x) 
        eval.function(trues=x$results$trues, 
                      preds=x$results$preds, 
                      y_train=x$train[,3], ...)), .id="fold")
    else fold.res <- as.data.frame(t(eval.function(trues=wfRes$results$trues, 
                                     preds=wfRes$results$preds, 
                                     y_train=wfRes$train[,3], ...)) )
  }
  
  fold.res
}

#' Estimate error using a chosen method
#'
#' @param data a data frame
#' @param form a formula for learning
#' @param estimator the name of an error estimator function
#' @param est.pars a named list of arguments to feed to \code{estimator}
#' @param workflow the name of the workflow to use for making predictions
#' @param wf.pars a named list of arguments to feed to \code{workflow}
#' @param evaluator the name of the function to use to calculate evaluation results
#' @param eval.pars a named list of arguments to feed to \code{evaluator}
#' @param seed a seed to set before performing estimates
#' @param .verbose print stages
#'
#' @return The results of \code{evaluator} after applying \code{estimator} to the
#' learning task
#' 
#' @export
estimates <- function(data, form, estimator="kf_xval",
                      est.pars = list(nfolds=10, 
                                      fold.alloc.proc="Trand_SPrand"), 
                      workflow = "simple_workflow", wf.pars=NULL, 
                      evaluator = "evaluate", eval.pars=NULL,
                      seed=1234, .verbose=FALSE){
  
  if(!is.null(seed)) set.seed(seed)
  
  if(.verbose) cat("\tRunning procedure...")
  res <- do.call(estimator, c(list(data=data, form=form, 
                                   FUN=get(workflow, mode="function")), 
                                est.pars, wf.pars))
  if(.verbose) cat("\tEstimating performance...")
  est.res <- do.call(evaluator, c(list(wfRes=res), eval.pars))
  
  list(evalRes = est.res, rawRes = res)
}



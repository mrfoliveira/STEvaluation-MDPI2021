#' From spdep neighbour object to data frame usable with ggplot2
#'
#' @param nb nb object obtained with \code{spdep::dnearneigh} or similar
#' @param stations sf object containing geometry points of each location
#'
#' @return data.frame with lon, lat for each point and lon_to and lat_to
#' as xend and yend for geom_segment
#' @export
nb2ggplot <- function(nb, stations){
  nb <- spdep::nb2listw(nb, zero.policy = TRUE)
  n <- length(attributes(nb$neighbours)$region.id)
  from <- rep(1:n, sapply(nb$neighbours,length))
  to <- unlist(nb$neighbours)
  # remove lines with no neighbours
  old.to <- to
  if(length(which(to==0))){
    isolated <- from[which(to==0)]
    from <- from[-which(to==0)]
    to <- to[-which(to==0)]  
  }
  
  DA = data.frame(
    from = from,
    to = to,
    weight = unlist(nb$weights)
  )
  if(length(which(old.to==0))){
    DA <- rbind(DA, data.frame(from = isolated, 
                               to = rep(NA, NROW(isolated)), 
                               weight = rep(NA, NROW(isolated))))
  }
    
  DA = cbind(DA, sf::st_coordinates(stations)[DA$from,], 
             sf::st_coordinates(stations)[DA$to,])
  colnames(DA)[4:7] = c("lon","lat","lon_to","lat_to")
  
  DA
}


#' Summarize the results of one in-set/out-set experiment
#'
#' @param one_exp_res a list containing two slots: \code{out_estRes} containing a 
#' named vector of metrics estimated in out-set data, and \code{in_estRes} containing a
#' list of data frames where each column corresponds to a metric and each row to a 
#' repetition/iteration of an estimator used on in-set data
#'
#' @return A data frame with a first column containing a summary (e.g., the mean) of
#' metrics measured in the out-set data and further columns containing summaries of 
#' metrics estimated in the in-set data
#' 
#' @seealso \code{\link{run_one_experiment}}
#' 
#' @export
summarize_one_exp <- function(one_exp_res){
  resTab <- dplyr::left_join(tidyr::gather(one_exp_res$out_estRes$evalRes, "metric", "real"), 
                  dplyr::bind_rows(lapply(one_exp_res$in_estRes, 
                         function(y){ evalRes <- y$evalRes;
                           if("fold" %in% colnames(evalRes)){ 
                             tidyr::gather(evalRes, "metric", "estimated", -.data$fold); 
                           }else{ 
                             tidyr::gather(evalRes, "metric", "estimated");}
                           }), .id="estimator"), by="metric")
  
  args <- dplyr::bind_rows(lapply(one_exp_res$in_estRes, 
                                  function(y){ 
                                    args <- names(y$params) %in% c(
                                         "tr.perc", "ts.perc",
                                         "nfolds", "fold.alloc.proc", 
                                         "s.buffer", "t.buffer", 
                                         "window", "removeSP", "nreps");
                                    params <- as.data.frame(y$params[which(args)]);
                                    params}), .id="estimator")
  
  if(all(c("nreps", "nfolds") %in% colnames(args)))
    args <- args %>%
      dplyr::mutate(nfolds = ifelse(is.na(.data$nfolds) & !is.na(.data$nreps), 
                                    .data$nreps, .data$nfolds)) %>%
      dplyr::select(-.data$nreps);
  list(resTab = resTab, params = args)
}

#' Transform a multi-level list of summarized results into a table
#'
#' @param sumRes A multi-level list of summarized results where the first level 
#' corresponds to learning model used in the experiment, the second level 
#' contains a list for each grid size of artificial data set, the third level 
#' contains a list for each time series size, and the next level contains a 
#' data frame with the results obtained in the out-set (gold-standard or 
#' "real" error) as well as estimated errors for different estimators
#' (in wide format)
#'
#' @return A data frame containing columns identifying the learning model,
#' grid size, time series size, type of STARMA used to generate the data,
#' order of STARMA used to generate, number of iteration of the generation process
#' with those settings, lag embed order, gold standard error (that of the out-set),
#' name of error estimator and estimated error (on the in-set), in long format
#' 
#' @export
#' 
#' @import dplyr
sumRes2Tab <- function(sumRes){
  
  sumResTab <- dplyr::bind_rows(lapply(sumRes, function(d)
    dplyr::bind_rows(lapply(d, function(x) 
      dplyr::bind_rows(x,
        .id = "t_size")),
      .id = "g_size")),
    .id = "model") %>%
    tidyr::separate(.data$gen_model, c("gen_type", "gen_order"),
             sep="\\_M\\_") %>%
    tidyr::separate(.data$gen_order, c("gen_order", "gen_it"), "\\.") %>%
    dplyr::mutate(lag_order = gsub("L\\_", "", .data$lag_order)) %>% 
    dplyr::mutate_at(vars(.data$model:.data$metric), as.factor)
  
  sumRes_real <- sumResTab %>% dplyr::select(.data$model:.data$real)
  sumResTab <- sumResTab %>% dplyr::select(-.data$real)
  sumRes_others <- sumResTab %>%
    tidyr::gather(.data$estimator, .data$estimated, 9:ncol(sumResTab)) %>%
    dplyr::mutate_if(is.character, as.factor)
  sumResTab <- dplyr::left_join(sumRes_real, sumRes_others)
  
  sumResTab
}

#' Transform a multi-level list of summarized results into a table
#'
#' @param sumRes A multi-level list of summarized results where the first level 
#' corresponds to learning model used in the experiment, the second level 
#' contains results for each data set
#' 
#' @param statFUN a function to summarize the evaluation metrics. Default is \code{mean}
#' @param na.rm whether to remove NAs in function \code{statFUN}
#'
#' @return A data frame containing columns identifying the learning model,
#' data set "gold standard"/"real" error/ (that of the out-set), 
#' name of error estimator and estimated error (on the in-set), in long format
#' 
#' @export
#' 
#' @import dplyr
realSumRes2Tab <- function(sumRes, statFUN=mean,
                              na.rm = FALSE){
  
  sumResTab <- dplyr::bind_rows(lapply(sumRes, function(x) 
    dplyr::bind_rows(lapply(x, function(y) 
      summarize_one_exp(y)), 
              .id="data")), 
    .id="model") %>%
    dplyr::mutate_at(vars(.data$model:.data$metric), as.factor)
  
  sumRes_real <- sumResTab %>% 
    dplyr::select(.data$model:.data$real)
  sumResTab <- sumResTab %>% 
    dplyr::select(-.data$real)
  sumRes_others <- sumResTab %>%
    tidyr::gather(.data$estimator, .data$estimated, 4:ncol(sumResTab)) %>%
    dplyr::mutate_if(is.character, as.factor)
  sumResTab <- dplyr::left_join(sumRes_real, sumRes_others)
  
  sumResTab
}


#' Abbreviate default experiment method names for use in plots
#'
#' @param x a vector of method names
#'
#' @return A vector of abbreviated names
abbr_est_names <- function(x){
  x <- gsub("all", "a", x)
  x <- gsub("block", "b", x)
  x <- gsub("rand", "r", x)
  x <- gsub("contig", "c", x)
  x <- gsub("checker", "s", x)
  x <- gsub("sliding", "slW", x)
  x <- gsub("growing", "grW", x)
  x <- gsub("\\_SP", "S", x)
  x
}

#' Abbreviate default experiment method names for use in plots
#'
#' @param x a vector of method names
#' @export
#'
#' @return A vector of abbreviated names
abbr_fold_alloc_names <- function(x){
  x <- gsub("all", "a", x)
  x <- gsub("block", "b", x)
  x <- gsub("rand", "r", x)
  x <- gsub("contig", "c", x)
  x <- gsub("checker", "s", x)
  
  x <- gsub("\\Tr", "tR", x)
  x <- gsub("\\Ta", "tA", x)
  x <- gsub("\\Tb", "tB", x)
  x <- gsub("\\SPs", "sS", x)
  x <- gsub("\\SPc", "sC", x)
  x <- gsub("\\SPr", "sR", x)
  x <- gsub("\\SPa", "sA", x)
  
  x <- gsub("\\_t", "t", x)
  x <- gsub("\\_s", "s", x)
  x <- gsub("growing", "grW", x)
  x <- gsub("liding", "\\_slW", x)
  x <- gsub("slW", "\\_slW", x)
  x <- gsub("\\_\\_slW", "\\_slW", x)
  
  x
}

#' Compress results from all (artificial) experiments
#'
#' @param all.res a list with multiple levels (model,
#' grid size, series size and lists of full results of multiple experiments)
#' @param rmAllRaw a boolean indicating whether the whole rawRes should
#' be removed (defaults to FALSE). If TRUE, only "train" data will be
#' removed from each set of results
#'
#' @return A multi-level list containing compressed results. Either all rawRes
#' is removed, or \code{train} is substituted by a vector of the number of instances, 
#' time and location IDs in the training set, in both \code{out_estRes} and
#' \code{in_estRes}.
#' 
#' @seealso \code{\link{run_all_experiments}}
#' 
#' @import dplyr
compressAllExps <- function(all.res, rmAllRaw=F){
  
  for(model in 1:length(all.res)){
    for(g.size in 1:length(all.res[[model]])){
      for(t.size in 1:length(all.res[[model]][[g.size]])){
        for(df in 1:length(all.res[[model]][[g.size]][[t.size]])){
          for(l in 1:length(all.res[[model]][[g.size]][[t.size]][[df]])){
            res <- all.res[[model]][[g.size]][[t.size]][[df]][[l]]
            all.res[[model]][[g.size]][[t.size]][[df]][[l]] <- compressExp(res, rmAllRaw=rmAllRaw)
          }
        }
      }
    }
  }
  all.res
}

#' Compress results from one gold-standard error experiment
#'
#' @param res A list containing full results of one experiment
#' (out_estRes and in_estRes)
#' @param ytrain_keep Should training target values still be kept
#' @param rmAllRaw a boolean indicating whether the whole rawRes should
#' be removed (defaults to FALSE). If TRUE, only \code{train} data will be
#' substituted by a vector of the training values of the target variable
#'
#' @return A list containing compressed results of one experiment. 
#' Either all rawRes is removed, or \code{train} substituted by a vector 
#' of the training values of the target variable
#' in both \code{out_estRes} and \code{in_estRes}
#' 
#' @seealso \code{\link{summarize_one_exp}}, \code{\link{run_one_experiment}}
#' @export
compressExp <- function(res, ytrain_keep=T, rmAllRaw=F){
  
  res$out_estRes <- compressRes(res$out_estRes, ytrain_keep=ytrain_keep, rmAllRaw = rmAllRaw)
  
  for(in_est in 1:length(res$in_estRes)){
    res$in_estRes[[in_est]] <- compressRes(res$in_estRes[[in_est]], 
                                           ytrain_keep = ytrain_keep, 
                                           rmAllRaw = rmAllRaw)
  }
  
  res
}
  

#' Compress results from one call to estimates
#'
#' @param res A list containing full results of one call to estimates
#' @param ytrain_keep Should training target values still be kept
#' @param rmAllRaw a boolean indicating whether the whole rawRes should
#' be removed (defaults to FALSE). If TRUE, only \code{train} data will be
#' substituted by a vector of the training values of the target variable
#'
#' @return A list containing compressed results of one experiment. 
#' Either all rawRes is removed, or \code{train} substituted by a vector 
#' containing only the training target values
#' 
#' @export
compressRes <- function(res, ytrain_keep = TRUE, rmAllRaw=F){
  # substitute distance matrix in parameters by its summary
  if("t.dists" %in% names(res$params))
    res$params$t.dists <- summary(as.vector(res$params$t.dists))
  if("s.dists" %in% names(res$params))
    res$params$s.dists <- summary(as.vector(res$params$s.dists))
  
  if(rmAllRaw){
    # remove rawRes from in_est
    res$rawRes <- NULL
  }else{
    # remove train from rawRes in in_est
    if("train" %in% names(res$rawRes)){
      train <- res$rawRes$train
      if(ytrain_keep) res$rawRes$train <- list(y_train = train[,3], 
                                         stats= list(c(ntimes=length(unique(train[,1])), 
                                                       nstations=length(unique(train[,2])))))
      else res$rawRes$train <- list(stats= list(c(ntimes=length(unique(train[,1])), 
                                                  nstations=length(unique(train[,2])),
                                                  ntrain = nrow(train))))
    }else{
      if(length(res$rawRes)>0){
        for(f in 1:length(res$rawRes)){
          if("train" %in% names(res$rawRes[[f]])){
            train <- res$rawRes[[f]]$train
            if(ytrain_keep) res$rawRes[[f]]$train <- list(y_train = train[,3], 
                                                    stats= list(c(ntimes=length(unique(train[,1])), 
                                                                  nstations=length(unique(train[,2])))))
            else res$rawRes[[f]]$train <- list(stats= list(c(ntimes=length(unique(train[,1])), 
                                               nstations=length(unique(train[,2])),
                                               ntrain = nrow(train))))
          }
        }
      }
    }  
  }
  
  res
}

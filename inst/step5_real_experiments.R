cat("\n\n************************************************\nEXPERIMENTS ON REAL DATA SETS\n************************************************\n\n")

library(dplyr)
library(ranger)
library(rpart)
library(earth)
library(doParallel)

# CHANGE NUMBER OF CORES
NCORES <- 4
NUM_THREADS <- 2
# note that NCORES x NUM_THREADS will be actually used when running ranger
cat(paste("\nUsing", NCORES, "cores and up to", NCORES*NUM_THREADS, "ranger threads\n\n"))
registerDoParallel(cores=NCORES)

DATA_PATH <-  "../data/"
RES_PATH <- "../results/"

if(!dir.exists(RES_PATH)){
  dir.create(RES_PATH, recursive = TRUE)
  cat(paste0("Created folder ", RES_PATH, "\n"))
} 
cat(paste0("Saving results to ", RES_PATH, "\n\n"))

## GENERAL PARAMETERS

NORP <- 0.2
MIN_TRAIN <- 100
IN_SET_PERC <- 0.8

SEED <- 1234
# RNGkind(sample.kind = "Rounding") # used in conference version of the paper

MODELS <- c("rpart", "earth", "lm", "ranger")
IMBALANCED_DS <- c("MESApol", "NCDCPprec", "TCEQOozone",
                   "TCEQTtemp", "TCEQWwind", "RURALpm10",
                   "BEIJno", "BEIJpm10", "BEIJwind", "BEIJpm25")

ALLOC_LABELS <- c(Trand_SPrand="tRsR", Tall_SPrand="tAsR",
                  Tblock_SPall="tBsA", Trand_SPall="tRsA",
                  Tblock_SPrand="tBsR")

# CROSS-VALIDATION PARAMETERS

XVAL_ALLOCS <- c("Trand_SPrand",
                 "Tall_SPrand",  
                 "Tblock_SPall",
                 "Trand_SPall", 
                 "Tblock_SPrand")
XVAL_variants <- lapply(XVAL_ALLOCS, 
                        function(x) list(nfolds=9, fold.alloc.proc=x,
                                         .parallel = FALSE,
                                         .verbose = FALSE))
names(XVAL_variants) <- paste0("CV.9_", ALLOC_LABELS[XVAL_ALLOCS])

# PREQUENTIAL PARAMETERS

PRE_ALLOCS <- c("Tblock_SPall", "Tblock_SPrand")
PRE_variants <- list()
for(win in c("growing", "sliding")){
  for(rmSP in c(TRUE, FALSE)){
    y <- lapply(PRE_ALLOCS,
                function(x) list(nfolds=9, 
                                 fold.alloc.proc=x, 
                                 window=win, 
                                 removeSP=rmSP,
                                 .parallel = FALSE,
                                 .verbose = FALSE))
    names(y) <- paste0("PRE.9_", ALLOC_LABELS[PRE_ALLOCS],"_", ifelse(win == "growing", "grW", "slW"))
    if(rmSP) names(y) <- paste0(names(y), "_rmSP")
    PRE_variants <- c(PRE_variants, y)
  }
}

# OTHER OOS PARAMETERS

MC_variants <- list("MC.44.6"=list(tr.perc=0.44, # at 50\% of nfolds-1 size
                                   ts.perc=0.06, # at 50\% of 1 fold size
                                   nreps=9,
                                   .parallel = FALSE,
                                   .verbose = FALSE), 
                    "MC.53.7"=list(tr.perc=0.53, # at 60\% of nfolds-1 size
                                   ts.perc=0.07, # at 60\% of nfolds-1 size
                                   nreps=9,
                                   .parallel = FALSE,
                                   .verbose = FALSE)) 

HO_variants <- list("HO.80" = list(tr.perc=0.8,
                                   .verbose = FALSE),
                    "HO.89" = list(tr.perc=0.89,
                                   .verbose = FALSE))

# NON-DEPENDENT CV PARAMETERS

NDXVAL_ALLOCS <- c(rep("Trand_SPrand", 3), 
                   "Tall_SPrand",
                   "Tblock_SPall")
NDXVAL_variants <- lapply(NDXVAL_ALLOCS, 
                          function(x) list(nfolds=9, fold.alloc.proc=x,
                                           s.buffer=NA, t.buffer=NA,
                                           .parallel = FALSE,
                                           .verbose = FALSE))
names(NDXVAL_variants) <- paste0("CV.9_", ALLOC_LABELS[NDXVAL_ALLOCS])
names(NDXVAL_variants) <- paste0(names(NDXVAL_variants), c("_buffST", 
                                                           "_buffT", "_buffS", "_buffS", "_buffT"))

# set buffer of the other dimension to NULL
for(v in grep("buffS$", names(NDXVAL_variants))){
  NDXVAL_variants[[v]]$t.buffer <- NULL  
}
for(v in grep("buffT$", names(NDXVAL_variants))){
  NDXVAL_variants[[v]]$s.buffer <- NULL  
}

NDXVAL_variants <- c(NDXVAL_variants, list(CV.9_tBsR_buffSTM=list(nfolds = 9, 
                                                                  fold.alloc.proc="Tblock_SPrand",
                                                                  s.buffer = 1, t.buffer = 1,
                                                                  .parallel = FALSE,
                                                                  .verbose = FALSE)))

NDXVAL_variants <- NDXVAL_variants[c("CV.9_tAsR_buffS", "CV.9_tRsR_buffS",
                                     "CV.9_tBsR_buffSTM", "CV.9_tBsA_buffT")]

# EVALUATION METRICS OPTIONS

EVAL_FUN <- regMetrics
EVAL_PARS <- list(eval.function=EVAL_FUN, .keptTrain = TRUE, norm = TRUE, util = FALSE)

# CHECK LONGITUDE AND LATITEU TCEQ DATA
dfnms <- c('MESApol', 'NCDCPprec', 'NCDCSsol',
           'TCEQOozone', 'TCEQTtemp', 'TCEQWwind', 
           'COOKwater', 'COOKtemp', 'COOKcond', #) 
           'SACtemp', 'RURALpm10', 'BEIJno', 'BEIJpm10', #) 
           'BEIJwind', 'BEIJpm25', 'BEIJhum', 'BEIJtemp')

# RUNNING EXPERIMENTS

cat("Running experiments...\n\n")
for(m in MODELS){
  
  foreach::foreach(i=1:length(dfnms)) %do% {
    dfnm <- dfnms[i]
    
    cat(paste("\n\nExperiments with", dfnm, ":\n\n"))
    
    cat("Loading data set...\n")
    fname <- paste0(DATA_PATH, "inds_df_", dfnm, ".Rdata")
    if(file.exists(fname)) load(fname)
    
    # cat("\nCalculating spatial distance matrix...\n")
    # calculate spatial distance matrix for non-dependent X-val
    s.dist <- norm_scale(get_spatial_dist_mat(ind_df$stations, site_id = "station"))
    # cat("Calculating temporal distance matrix...\n")
    # calculate temporal distance matrix for non-dependent X-val
    t.dist <- norm_scale(get_time_dist_mat(ind_df$df$time))
    
    BETAS <- ind_df$betas
    ALPHA <- ind_df$alpha
    # fixing NDXVAL_variants
    for(v in 1:length(NDXVAL_variants)){
      if(!is.null(NDXVAL_variants[[v]]$s.buffer)){
        NDXVAL_variants[[v]][["s.dists"]] <- s.dist
        if(is.na(NDXVAL_variants[[v]]$s.buffer)) 
          NDXVAL_variants[[v]]$s.buffer <- max(BETAS)/ALPHA
      } 
    }
    for(v in 1:length(NDXVAL_variants)){
      if(!is.null(NDXVAL_variants[[v]]$t.buffer)){
        NDXVAL_variants[[v]][["t.dists"]] <- t.dist
        if(is.na(NDXVAL_variants[[v]]$t.buffer)) 
          NDXVAL_variants[[v]]$t.buffer <- max(BETAS)/(1-ALPHA)
      } 
    }
    
    ALL_variants <- c(HO_variants, MC_variants, 
                      PRE_variants, XVAL_variants,
                      NDXVAL_variants)
    IN_ESTIMATORS <- c(rep("t_oos", length(HO_variants)),
                       rep("t_oos_mc", length(MC_variants)),
                       rep("prequential_eval", length(PRE_variants)),
                       rep("kf_xval", length(XVAL_variants)),
                       rep("nd_kf_xval", length(NDXVAL_variants)))
    
    #  set workflow parameters
    WF_PARS_BASE <- list(model=m, handleNAs="centralImput", 
                         min_train = MIN_TRAIN,
                         nORp =  NORP)
    
    if(m == "ranger"){
      WF_PARS <- c(WF_PARS_BASE, list(num.trees = 500, 
                                      num.threads=NUM_THREADS, 
                                      verbose = FALSE))
    }else{
      WF_PARS <- WF_PARS_BASE
    }
    
    # set evaluation metrics parameters
    if(dfnm %in% IMBALANCED_DS){
      EVAL_PARS$util <- TRUE
      EVAL_PARS$util.parms <- list(rel="auto", thr=0.9, cf=1.5)
    }else{
      EVAL_PARS$util <- FALSE
      EVAL_PARS$util.parms <- NULL
    }
    
    # run experiment
    cat(paste0("Running experiment with ", m,"...\n"))
    one_res <- run_one_experiment(
      ind_df$df, IN_SET_PERC, as.formula("value~."), 
      in_estimators=IN_ESTIMATORS, 
      in_est.pars=ALL_variants,
      out_estimator = "t_oos",
      out_est.pars=list(tr.perc = IN_SET_PERC),
      workflow = "simple_workflow", 
      wf.pars = WF_PARS,
      evaluator = "evaluate", 
      eval.pars = EVAL_PARS,
      seed = SEED, site_id = "station", time = "time",
      .parallel = TRUE, .verbose = FALSE)
    
    # COMPRESS OR NOT?
    cat("Compressing results...\n")
    one_res <- compressExp(one_res)
    
    cat("\nSaving results for ", dfnm, " and model", m, "...\n")
    save(one_res, file=paste0(RES_PATH, "real_res_",m,"_", dfnm, ".Rdata")                         )
    cat("\nSaved results!\n")
    
    rm(ind_df)
    rm(one_res)
    gc()
    
    1
  }
}
cat("\n")

cat("\nDone!\n\n******************\n\n")
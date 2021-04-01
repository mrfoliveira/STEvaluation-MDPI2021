cat("\n\n************************************************\nCALCULATE SPATIO-TEMPORAL INDICATORS\n************************************************\n\n")

library(doParallel)

# CHANGE NUMBER OF CORES
NCORES <- 4
NUM_SPLITS <- 1
# note that NCORES x NUM_THREADS will be actually used when running ranger
cat(paste("\nUsing", NCORES, "cores and up to", NCORES, "ranger threads\n\n"))
registerDoParallel(cores=NCORES)

#--------------------------------------

# PATHS

DATA_PATH <-  "../data/"

#--------------------------------------

cat("Loading data sets...\n")
load(paste0(DATA_PATH, "dfs.Rdata"))

dfnms <- c('MESApol', 'NCDCPprec', 'NCDCSsol', 'TCEQOozone', 'TCEQTtemp', 
           'TCEQWwind', 'COOKwater', 'COOKtemp', 'COOKcond', 'SACtemp', 
           'RURALpm10', 'BEIJno', 'BEIJpm10', 'BEIJwind', 'BEIJpm25', 
           'BEIJhum', 'BEIJtemp')

IMBALANCED_DS <- c("MESApol", "NCDCPprec", "TCEQOozone",
                   "TCEQTtemp", "TCEQWwind", "RURALpm10",
                   "BEIJno", "BEIJpm10", "BEIJwind", "BEIJpm25")



inds_df <- list()
for(dfnm in dfnms){

  ALPHA <- length(unique(data_list[[dfnm]]$df$station)) / length(unique(data_list[[dfnm]]$df$time))
  ALPHA <- ifelse(ALPHA<0.5, 0.25, 0.5)
  BETAS <- c(0.0250, 0.0375, 0.0500)
  
  cat(paste0("Get spatio-temporal indicators for ", dfnm,"...\n"))
  ind_df <- get_full_indicators(data_list[[dfnm]]$df, data_list[[dfnm]]$stations,
                                k=8, var="value",
                                betas=BETAS, alpha=ALPHA,
                                stats = c("mean", "weighted.mean", "sd"), 
                                ratios2add = c(TRUE,TRUE,FALSE),
                                # .parallel=TRUE, nsplits=NUM_SPLITS,
                                .parallel=FALSE,
                                time_id="time", site_id="station") 
  
  # fix formats
  ind_df <- as.data.frame(ind_df)
  if(grepl("BEIJ", dfnm)) ind_df$time <- lubridate::ymd_hms(ind_df$time)
  else ind_df$time <- lubridate::ymd(ind_df$time)
  
  ind_df <- list(df=ind_df, alpha=ALPHA, betas=BETAS, stations=data_list[[dfnm]]$stations)
  
  cat(paste0("\nSaving indicator data for data set ", dfnm,"...\n"))

  
  save(ind_df, file=paste0(DATA_PATH, "inds_df_", dfnm,".Rdata"))
  gc()
  
  cat(paste0("Saved indicator data for data set ", dfnm,"!\n"))

  inds_df[[dfnm]] <- ind_df
 } 

save(inds_df, file=paste0(DATA_PATH, "inds_df.Rdata"))

cat("\nSpatio-temporal indicators calculated!\n\n******************\n\n")

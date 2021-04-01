cat("\n\n************************************************\nGENERATE ARTIFICIAL DATA\n************************************************\n\n")

DATA_PATH <- "../data/"

if(!dir.exists(DATA_PATH)){
  dir.create(DATA_PATH, recursive = TRUE)
  cat(paste0("Created folder ", DATA_PATH, "\n"))
} 
cat(paste0("Saving data to ", DATA_PATH, "\n\n"))


## STARMA SPECIFICATIONS
Nsites <- c((8+2)^2,
            (20+2)^2)
trash <- 100
Ntimes <- c(150,300) + 3

#Ntimes <- t.nfolds*t.fold.size*(1 + (1-in_set_perc)/in_set_perc) + LAG_avail + trash 
mtypes <- c("STARMA", "STAR", "STMA", "NL_STAR")
coef_specs <- list(M_2_10 = list(c_10=c(-2,2), c_11=c(-2,2), c_20=c(-1,1), c_21=0),
                   M_2_01 = list(c_10=c(-2,0), c_11=0, c_20=c(-1,1), c_21=c(-2,1)),
                   M_2_11.1 = list(c_10=c(-1.227,0.773), c_11=c(0.733,1.277), c_20=c(-0.227,1.773), c_21=-0.733),
                   M_2_11.2 = list(c_10=c(-1.755,0.245), c_11=c(-1.755,1.755), c_20=c(-0.755,0.755), c_21=0.245),
                   M_2_11.3 = list(c_10=c(0.227,1.773), c_11=c(-1.319,0.277), c_20=c(-0.773,0.733), c_21=-0.227),
                   M_2_11.4 = list(c_10=c(-1.378,-0.622), c_11=c(0.622,1.378), c_20=c(-0.378,0.378), c_21=0.622))
ncoefs <- c(4, 4, rep(1,4))

cat("\nGenerating simulated data...\n")
# generate the datasets
data <- generate_multiple_datasets(Nsites = Nsites, Ntimes = Ntimes, mtypes = mtypes, 
                                   coef_specs = coef_specs, ncoefs = ncoefs, trash = trash, 
                                   init_seed = 1234)
save(data, file=paste0(DATA_PATH, "data.Rdata"))


## LAGGING SPECIFICATIONS
LAG_use <- 3
SLAGS <- list(list(c(1,1,0)))
min_time <- rep(max(LAG_use) + 1, length(LAG_use))

cat("\nEmbedding data...\n")
# lag the datasets
lagged_data <- lag_multiple_datasets(data, LAG_use, SLAGS, min_time = min_time)
cat(paste("Saving lagged data to", DATA_PATH,"\n"))
save(lagged_data, file=paste0(DATA_PATH, "lagged_data.Rdata"))
cat("\nDone!\n\n******************\n\n")

source("~/send_mail.R")
send_mail(msg="Generated artificial data")

cat("\n\n************************************************\nDOWNLOAD AND READ REAL DATA\n************************************************\n\n")

library(dplyr)
library(sf)

DATA_PATH <- "../data/"
if(!dir.exists(DATA_PATH)){
  dir.create(DATA_PATH, recursive = TRUE)
  cat(paste0("Created folder ", DATA_PATH, "\n"))
} 

####
# COSTK DATA
####

download.file(url="http://www.di.uniba.it/~appice/software/COSTK/data/dataset.zip", 
              destfile = paste0(DATA_PATH, "/CostK_data"))
unzip(paste0(DATA_PATH, "/CostK_data"), exdir=DATA_PATH)

get_costk_data <- function(folder, var, base_path, df_path){
  fname <- paste0(base_path, df_path, "/", folder, "/pos.txt")
  
  print(var)
  
  stations <- read.csv(fname, sep="\t", header=T)
  if(grepl("TCEQ", folder) | grepl("SR", folder))
    stations <- df2site_sf(stations, "station", 
                           lon="y", lat="x", 4326)   
  else
    stations <- df2site_sf(stations, "station", 
                           lon="x", lat="y", 4326)   
  
  fname <- paste0(base_path, df_path, "/", folder, 
                  "/", var, ".txt")
  df <- read.csv(fname, sep="\t", header=F) %>%
    dplyr::mutate(time = as.Date(1:n(), origin=as.Date("1900-01-01"))) %>%
    tidyr::gather(station, value, -time) %>% 
    dplyr::mutate(station=as.numeric(substr(station,2,1E6))) %>%
    as.data.frame()
  
  list(stations=stations, df=df)
}

DF_PATH <- "data"

folders <- c("MESA", paste0("NCDC", c("P", "S", "T")),
             "SAC", paste0("TCEQ", c("O", "T", "W")), "SR")
vars <- c("pol", "prec", "sol", "temp", "temp", "ozone",
          "temp", "wind", "dif")
names(vars) <- folders

data_list <- list()
for(f in folders){
  print(f)
  data_list[[paste0(f, vars[f])]] <- get_costk_data(f, vars[f], DATA_PATH, DF_PATH)
}

###
# RURAL DATA
###

library(spacetime)
data(air)
rural_df <- as.data.frame(STFDF(stations, dates, 
                                data.frame(PM10 = as.vector(air)))) %>%
  select(-endTime, -timeIndex) %>%
  filter(!is.na(PM10))
colnames(rural_df)[c(3,5)] <- c("station", "value")

# read stations
rural_sf <- df2site_sf(rural_df %>% select(coords.x1, coords.x2, station) %>% distinct(), 
                       "station", lon="coords.x1", lat="coords.x2", 4326)

rural_df <- rural_df %>% select(-coords.x1, -coords.x2) %>% as.data.frame()

data_list[["RURALpm10"]] <- list(stations = rural_sf, df = rural_df)

###
# COOK FARM
###

library(GSIF)
data(cookfarm)

colnames(cookfarm$profiles)[1] <- "station"

farm_sf <- df2site_sf(cookfarm$profiles %>% 
                        dplyr::filter(UHDICM==0,
                               station %in% cookfarm$readings$SOURCEID) %>%
                        dplyr::select(station, Easting, Northing) %>% 
                        distinct() %>%
                        as.data.frame(), 
                      "station", lon="Easting", lat="Northing", crs=
                      "+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs")

farm_df <- cookfarm$readings[,1:5]
colnames(farm_df)[c(1,2)] <- c("station", "time")

vars <- c("water", "temp", "cond")
for(i in 3:5){
  x <- farm_df[,c(1,2,i)]
  colnames(x) <- c("station", "time", "value")
  
  data_list[[paste0("COOK", vars[i-2])]] <- list(stations=farm_sf, df=x)
}

####
# BEIJING DATA
####

download.file(url="https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/Air20Quality20Data.zip",
              destfile=paste0(DATA_PATH, "AirQuality"))
unzip(paste0(DATA_PATH, "AirQuality"), exdir = DATA_PATH)

DF_PATH <- "Beijing/"

beij_sf <- read.csv(paste0(DATA_PATH, DF_PATH, "Station.txt")) %>%
  dplyr::rename(station=station_id) %>%
  df2site_sf("station", "longtitude","latitude", 4326) 

beij_df <- read.csv(paste0(DATA_PATH, DF_PATH, "CrawledData.txt")) %>%
  dplyr::mutate(time = lubridate::mdy_hms(time)) %>%
  dplyr::rename(station=station_id) %>%
  dplyr::select(-weather) %>%
  as.data.frame()

vars <- paste0("BEIJ", c("pm25", "pm10", "no", "temp", "pres", "hum", "wind"))
colnames(beij_df)[c(-1,-2)] <-  vars

for(i in 1:length(vars)){
  v <- vars[i]
  x <- beij_df[,c("station", "time", v)]
  colnames(x) <- c("station", "time", "value")
  data_list[[v]] <- list(stations = beij_sf, df = x)
}

cat("Fixing to hourly...\n\n")
for(d in names(data_list)[grep("BEIJ", names(data_list))]){
  df <- data_list[[d]]$df
  df <- df %>% mutate(time = lubridate::round_date(lubridate::ymd_hms(time), "hour")) %>%
    dplyr::group_by(time, station) %>%
    dplyr::summarize(value = median(value, na.rm=T), .groups="drop")
  data_list[[d]]$df <- as.data.frame(df)
}

# Removing empty values
for(d in names(data_list)){
  df <- data_list[[d]]$df
  df <- df[!is.na(df$value), ]
  data_list[[d]]$df <- as.data.frame(df)
}

# sort by data set size
data_list <- data_list[names(sort(sapply(data_list, function(x) nrow(x$df))))]

# save data

save(data_list, file=paste0(DATA_PATH, "/dfs.Rdata"))

cat("\nData saved!\n\n******************\n\n")

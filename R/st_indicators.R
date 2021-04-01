#' Calculate temporal distance matrix
#' 
#' A function that calculates a distance matrix of time-stamps
#' @param times A vector of time-stamps
#' @param origin A date to use as origin for \code{difftime}
#' 
#' @return a matrix of distances. Row and column names
#' are a concatenation of "TIME_" and the time-stamp.
#' 
#' @export
get_time_dist_mat <- function(times, origin=min(times)){
  # unique timestamps
  times <- sort(unique(times))
  timz <- difftime(times, origin)
  dists <- as.matrix(stats::dist(timz, diag=TRUE, upper=TRUE))
  # distance to present/future will be changed to Inf later
  rownames(dists) <- paste0("TIME_", times)
  colnames(dists) <- paste0("TIME_", times)
  dists
}

#' Get spatio-temporal neighbourhood
#'
#' A function that calculates the observations that
#' are within a spatio-temporal neighbourhood of a 
#' certain radius of a time and location. 
#' 
#' @details The spatio-temporal distance is defined as
#' \deqn{D_{i,j} = d_{i,j} x \alpha + t_{i,j}
#' x (1-\alpha)}
#' where d_{i,j} is the spatial distance between
#' locations, {t_{i,j}} is the temporal distance
#' between time-stapms and \eqn{\alpha} is a weighting factor.
#' 
#' @param site a location ID
#' @param time a time-stamp
#' @param radius a radius of spatio-temporal distance
#' @param t_dist_mat a matrix of normalized temporal distances
#' between time-stamps (rownames and colnames should be
#' a concatenation of "TIME_" and the time-stamp)
#' @param s_dist_mat a matrix of normalized spatial distances
#' between locations (rownames and colnames should be
#' a concatenation of "SITE_" and the location IDs)
#' @param alpha a weighting factor for the spatio-temporal distance
#' @param time_id the name to give to the column
#' of time-stamps (Default: time)
#' @param site_id the name to give to the column
#' of location IDs (Default: site_id)
#' 
#' @details Note that \code{radius} should always be a number between
#' zero and \code{min(alpha, alpha-1)}.
#' Also note that if \code{alpha} is set to 1, then instead of a cone, 
#' the neighbourhood will have the shape of a cylinder.
#' 
#' @return A data frame where each row describes a neighbour,
#' with the first two columns containing the location ID and
#' time-stamp of the central observation, followed by two
#' columns with the neighbouring location ID and time-stamp,
#' and a final colunm containing the spatio-temporal distance
#' between the two.
#' 
#' @references \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.258.8742&rep=rep1&type=pdf}
#' 
#' @import dplyr
get_st_neighbours <- function(site, time, radius, t_dist_mat, s_dist_mat, 
                              alpha, time_id="time", site_id="site_id"){
  
  # assertions about radius given alpha
  assertthat::assert_that(alpha>0, alpha<1, radius>0, radius<min(alpha, 1-alpha))
  
  # find site and time indices in distance matrices
  site_i <- which(paste0("SITE_",site)==rownames(s_dist_mat))
  time_i <- which(paste0("TIME_",time)==rownames(t_dist_mat))
  
  neib_sites <- s_dist_mat[site_i, which(s_dist_mat[site_i,] <= radius / alpha), drop=F]
  neib_sites <- as.data.frame(t(neib_sites))
  
  neib_times <- t_dist_mat[time_i, which(t_dist_mat[time_i,] <= radius / (1 - alpha)), drop=F]
  neib_times <- as.data.frame(t(neib_times))
  
  if(nrow(neib_times) == 0 | nrow(neib_sites) == 0)
    return(NULL)
  else{
    colnames(neib_sites) <- "s_dist"
    neib_sites$site <- as.factor(rownames(neib_sites))
    
    colnames(neib_times) <- "t_dist"
    neib_times$time <- as.factor(rownames(neib_times))
    
    combns <- expand.grid(rownames(neib_sites), rownames(neib_times))
    colnames(combns) <- c("site", "time")
    combns <- combns %>%
      dplyr::left_join(neib_times, by="time") %>%
      dplyr::left_join(neib_sites, by="site") %>%
      dplyr::mutate(st_dist = .data$s_dist*alpha + .data$t_dist*(1-alpha)) %>%
      dplyr::filter(.data$st_dist <= radius) %>%
      dplyr::select(.data$site, .data$time, .data$st_dist) %>%
      as.data.frame()
    combns$site <- gsub("^SITE\\_", "", combns$site)
    combns$time <- gsub("^TIME\\_", "", combns$time)
    colnames(combns) <- c(paste0("neib_", c(site_id, time_id)), "st_dist")
    combns[[time_id]] <- as.character(time)
    combns[[site_id]] <- as.character(site)
    return(combns[,c(site_id, time_id, paste0("neib_", c(site_id, time_id)), "st_dist")])
  }
}

#' Get the spatio-temporal neighbourhoods of
#' all observations 
#' 
#' 
#' @param df a data frame of observations
#' @param time_id the name of the column containing time-stamps
#' @param site_id the name of the column containing location IDs
#' @param max_radius the maximum spatio-temporal distance
#' allowed to be included in a neighbourhood
#' @inheritParams get_st_neighbours
#' @param parallel Boolean indicating whether the code should
#' run in parallel. Default is FALSE
#' @param nsplits Number of subsets of rows to split the data 
#' frame into so they can be processed in parallel
#' @param vars Vector of character strings indicating the columns
#' whose values should be retrieved
#' 
#' @return A data frame where each row describes a neighbour,
#' with the first two columns containing the location ID and
#' time-stamp of the central observation, followed by two
#' columns with the neighbouring location ID and time-stamp,
#' a colunm containing the spatio-temporal distance between 
#' the two, and a final column containing the values of 
#' the variables in df at the neighbouring time and location.
#' 
#' @details The spatio-temporal distance is defined as
#' \deqn{D_{i,j} = d_{i,j} x \alpha + t_{i,j}
#' x (1-\alpha)}
#' where d_{i,j} is the spatial distance between
#' locations, {t_{i,j}} is the temporal distance
#' between time-stapms and \eqn{\alpha} is a weighting factor.
#' Note that the \code{radius} should always be a number between
#' zero and \code{min(alpha, alpha-1)}.
#' Also note that if \code{alpha} is set to 1, then instead of a cone, 
#' the neighbourhood will have the shape of a cylinder.
#' 
#' @references \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.258.8742&rep=rep1&type=pdf}
#' @seealso \code{\link{get_st_neighbours}}
#' 
#' @import dplyr
#' @import stringr
get_all_neib_vals <- function(df, max_radius,
                              t_dist_mat, s_dist_mat, alpha, vars, 
                              time_id, site_id, parallel=FALSE, nsplits=4){
  
  if(parallel){
    neib_df <- foreach::foreach(i=0:(nsplits-1), .combine='rbind',
                                .export=c("get_st_neighbours"),
                                .packages = c("dplyr", "stringr")) %dopar%{
                                  
                                  nrows <- nrow(df)
                                  subdf <- df[seq(i*ceiling(nrows/nsplits)+1,
                                                  min(nrows, (i+1)*ceiling(nrows/nsplits))),]
                                  dplyr::bind_rows(apply(subdf, 1, function(x) get_st_neighbours(site=x[site_id], 
                                                                                                 time=x[time_id],
                                                                                                 radius=max_radius,
                                                                                                 t_dist_mat=t_dist_mat, 
                                                                                                 s_dist_mat=s_dist_mat, 
                                                                                                 alpha=alpha,
                                                                                                 time_id=time_id, 
                                                                                                 site_id=site_id)))
                                } 
  }else{
    neib_df <- dplyr::bind_rows(apply(df, 1, function(x) get_st_neighbours(site=x[site_id], 
                                                                           time=x[time_id],
                                                                           radius=max_radius,
                                                                           t_dist_mat=t_dist_mat, 
                                                                           s_dist_mat=s_dist_mat, 
                                                                           alpha=alpha,
                                                                           time_id=time_id, 
                                                                           site_id=site_id)))    
  }
  
  
  cols <- c(site_id, time_id, vars)
  neib_df <- neib_df %>% dplyr::left_join(df %>% dplyr::select(one_of(cols)), 
                                          by=stats::setNames(c(site_id, time_id),
                                                             paste0("neib_", c(site_id, time_id)) ))
  neib_df
}

#' Get a spatio-temporal indicator from neighbourhood
#' data frame
#' 
#' Calculate a spatio-temporal indicator of a certain
#' variable within a spatio-temporal neighborhood of a 
#' certain radius.
#' 
#' @param all_neib_vals a data frame containing
#' information on observations spatio-temporal distance
#' to neighbours and variable values at the neighbouring
#' locations and times.
#' @param stat the name of a function that calculates a
#' statistic (e.g., "mean"). If the stat is "weighted.mean"
#' then the inverse of the spatio-temporal distance is used 
#' to weight the values of observations in the spatio-temporal 
#' neighbourhood
#' @param radius a value defining the maximum spatio-temporal
#' distance allowed for an observation to be considered
#' within a spatio-temporal neighbourhood
#' @param ind_name the name of the indicator column
#' @param var the name of the variable to summarize into 
#' an indicator
#' @param time_id the name of the column containing 
#' time-stamps
#' @param site_id the name of the column containing 
#' location IDs
#'
#' @return A data frame that contains a summary statistic
#' of the values found within the spatio-temporal neighbourhood
#' of a certain radius of each pair (location ID, time-stamp) in
#' \code{all_neib_vals}
#' 
#' @details The spatio-temporal distance is defined as
#' \deqn{D_{i,j} = d_{i,j} x \alpha + t_{i,j}
#' x (1-\alpha)}
#' where d_{i,j} is the spatial distance between
#' locations, {t_{i,j}} is the temporal distance
#' between time-stapms and \eqn{\alpha} is a weighting factor.
#' Note that \code{radius} should always be a number between
#' zero and \code{min(alpha, alpha-1)}, so the border conditions apply.
#' Also note that if \code{alpha} is set to 1, then instead of a cone, 
#' the neighbourhood will have the shape of a cylinder.
#' 
#' @references \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.258.8742&rep=rep1&type=pdf}
#' 
#' @import dplyr
get_st_indicator <- function(all_neib_vals, stat, radius, 
                             ind_name, var,
                             time_id="time", site_id="site"){
  
  stat_df <- all_neib_vals %>% 
    dplyr::filter(.data$st_dist <= radius) %>%
    dplyr::group_by_(time_id, site_id)
  if(stat=="weighted.mean"){
    stat_df <- stat_df %>%  
      dplyr::summarize_(value=paste0(stat, "(",var,", w=1/st_dist, na.rm=T)")) %>%
      dplyr::rename_(.dots=stats::setNames("value", ind_name))  
  }else{
    stat_df <- stat_df %>%  
      dplyr::summarize_(value=paste0(stat, "(",var,", na.rm=T)")) %>%
      dplyr::rename_(.dots=stats::setNames("value", ind_name))
  }
  stat_df
}

#' Get spatio-temporal indicators from a data frame
#' containing spatio-temporal information
#' 
#' Calculate spatio-temporal indicators of one or more
#' variables within a spatio-temporal neighborhood of one 
#' or more maximum radius (in terms of spatio-temporal 
#' distance).
#' 
#' @param df A data frame containing spatio-temporal
#' information
#' @param stations_sf An \code{sf} object containing
#' geographical information on the location of \code{df}
#' @param stats A vector containing the names of functions
#' that are to be used to calculate summarizing statistics
#' @param radiuses A vector of values defining the maximum
#' spatio-temporal distance allowed for an observation to be considered
#' within a spatio-temporal neighbourhood
#' @param vars The name of the variables to summarize into 
#' indicators
#' @param neib_type the type of neighborhood to consider.
#' Can be
#' \itemize{
#' \item \code{cone} (default) - a cone with the center of its base at 
#' the observation (spatial radius growing with time)
#' \item \code{reversed} - a cone with its peak at the observation
#' (spatial radius shrinking with time)
#' }
#' @param time_id The name of the column containing 
#' time-stamps in \code{df}
#' @param site_id The name of the column containing 
#' location IDs in \code{df}
#' @inheritParams get_all_neib_vals
#'
#' @return A data frame that contains summary statistics
#' of the values of \code{vars} found within the spatio-temporal 
#' neighbourhoods of the one or more \code{radiuses} of each pair
#'  (location ID, time-stamp) in \code{df}
#' 
#' @details The spatio-temporal distance is defined as
#' \deqn{D_{i,j} = d_{i,j} x \alpha + t_{i,j}
#' x (1-\alpha)}
#' where d_{i,j} is the spatial distance between
#' locations, {t_{i,j}} is the temporal distance
#' between time-stapms and \eqn{\alpha} is a weighting factor.
#' 
#' @references \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.258.8742&rep=rep1&type=pdf}
#' 
#' @import dplyr
get_st_indicators <- function(df, stations_sf, radiuses = c(0.1), 
                              stats = c("mean", "sd"), alpha=0.5, neib_type="cone",
                              time_id="time", site_id="site_id",
                              vars=c("value"), parallel=FALSE, nsplits=4){
  
  df[[site_id]] <- as.character(df[[site_id]])
  df[[time_id]] <- as.character(df[[time_id]])
  
  # get values within spatio-temporal neighbourhood of each observation
  inds_df <- df[,c(site_id, time_id)]
  
  # get space distance matrix
  s_dist_mat <- get_spatial_dist_mat(stations_sf, site_id)
  s_dist_mat <- norm_scale(s_dist_mat)
  
  # get time distance matrix
  t_dist_mat <- get_time_dist_mat(df[[time_id]])
  t_dist_mat <- norm_scale(t_dist_mat)
  # DISTANCE TO THE PRESENT/FUTURE IS INFINITE (distance of row to col)
  t_dist_mat[upper.tri(t_dist_mat, diag=T)] <- Inf
  
  assertthat::assert_that(all(s_dist_mat<=1), all(s_dist_mat>=0),
                          all(ifelse(t_dist_mat<=1, t_dist_mat>=0, t_dist_mat==Inf)))
  
  assertthat::assert_that(all(sapply(radiuses, function(r) any(t_dist_mat<=(r/(1-alpha))))),
      msg="Combination of radius and alpha makes it so that no timestamps are considered within 
       one or more of the neighbourhoods. Cannot calculate spatio-temporal indicators.")
  
  for(radius in radiuses){
      
    # reverse cone
    if(neib_type == "reversed"){
      new_t_dist_mat <- ifelse(t_dist_mat==Inf, Inf, radius/(1-alpha) - t_dist_mat)
      new_t_dist_mat <- ifelse(new_t_dist_mat<0, Inf, new_t_dist_mat)
    }else{
      new_t_dist_mat <- t_dist_mat
    }
    
    neib_vals <- get_all_neib_vals(df=df, max_radius=radius, 
                                 t_dist_mat=new_t_dist_mat, s_dist_mat=s_dist_mat,
                                   alpha=alpha, vars=vars, 
                                   time_id=time_id, site_id=site_id,
                                   parallel=parallel, nsplits=nsplits)
    
    for(stat in stats){
      for(var in vars){
        ind_name <- paste0(var, "_", stat, "_", radius)
        inds_df <- inds_df %>% dplyr::left_join(get_st_indicator(neib_vals, 
                                                                 stat = as.character(stat), 
                                                                 radius = as.numeric(radius),
                                                                 ind_name = ind_name,
                                                                 var=var,
                                                                 time_id=time_id, site_id=site_id),
                                                by=c(site_id, time_id))
      }
    }
  }
  
  inds_df
}


#' Add ratios between spatio-temporal neighborhood indicators
#'
#' @param df  A data frame of spatio-temporal indicators. 
#' Column names should be of type <variable name>_<indStat>_<radius>.
#' @param var A character string with the name of the variable
#' with indicators to add ratios
#' @param indStat The name of the summarizing stat that was used
#' to calculate the indicators
#'
#' @return A data frame including the original data of \code{df} 
#' and additional ratios between the indicators of subsequent radiuses.
#' 
#' @seealso \code{\link{get_st_indicators}}
#' 
add_ratios <- function(df, var, indStat="mean"){
  idxs <- grep(paste0(var,"_", indStat,"\\_[0-9]\\.[0-9]+$"), colnames(df))
  betas <- sort(sapply(idxs, function(i) as.numeric(gsub(paste0("^",var, "\\_", indStat,"_"), "", colnames(df)[i]))))
  
  nratios <- length(betas)
  for(j in seq(nratios,2,-1)){
    nm <- paste0(var,"_",indStat,"_r_", betas[j],"_", betas[j-1])
    n1 <- paste0(var,"_",indStat,"_", betas[j-1])
    n2 <- paste0(var,"_",indStat,"_", betas[j])
    df[[nm]] <- df[,n2] / df[,n1]
    df[[nm]] <- ifelse(df[,n1]==0, 0, df[[nm]])
  }
  df
}



#' Embed each time series in a spatio-temporal data set
#'
#' @param df data frame
#' @param var a character string, the name of the variable to embed
#' @param k a numeric, the embed size
#' @param time a character string, the column name identifying
#' the time of observation
#' @param station_id a character string, the column name
#' identifying the location of observation
#'
#' @return A data frame with extra columns var_Tm1, var_Tm2, ...,
#' var_Tm\(k-1\)
#' 
#' @import dplyr
embed_series <- function(df, var, k, time="time", station_id="station") {
  
  df <- as.data.frame(df)
  time_ids <- unique(df[[time]])
  assertthat::assert_that(lubridate::is.Date(time_ids) | lubridate::is.POSIXt(time_ids),
                          msg = "time should be a date or date time object")
  
  xs <- list()
  periods <- c()
  for(i in 1:length(unique(df[[station_id]]))){
    # get station ID
    s <- unique(df[[station_id]])[i]
    # get time series for this station
    x <- df[which(df[[station_id]]==s), c(time, station_id, var)]
    # get time in right format for xts
    colnames(x) <- c("time", "station", "value")
    
    # get time series periodicity
    ts <- xts::xts(x[[var]], order.by = x[["time"]])
    period <- as.numeric(xts::periodicity(ts)$difftime)
    # save period to check if it is the same across stations
    periods <- c(periods, period)
    
    # complete 
    x <- x %>% 
      tidyr::complete(time = seq(min(time), max(time), by=period)) %>%
      dplyr::arrange(time)
    
    ts.index  <- x[["time"]][-(1:(k-1))]
    ts.embed  <- stats::embed(x[ , var, drop=T], dimension = k)
    colnames(ts.embed) <- c(var, paste0(var,'_Tm', 1:(k-1)))
    ts.embed <- cbind(station=s, time=ts.index, as.data.frame(ts.embed))
    colnames(ts.embed)[1:2] <- c(station_id, time)
    ts.embed <- ts.embed[which(!is.na(ts.embed[,var])),]
    
    xs[[i]] <- ts.embed
  }
  xs <- dplyr::bind_rows(xs)
  if(length(unique(periods)) > 1) warning(paste("Warning: there are", length(unique(periods)),"different values for time series period at different stations"))
  
  df <- df %>% dplyr::left_join(xs, by=c(station_id, time, var))
  df <- df[order(df[[station_id]], df[[time]]),]
  df <- df[which(!(df[,time] %in% sort(time_ids)[1:(k-1)])),]
  
  df
}



#' Get time series embeds and spatio-temporal indicators
#'
#' @param var The name of the variable to summarize into indicators
#' @param stations An \code{sf} object containing
#' geographical information on the location of \code{df}
#' @param betas A vector of values defining the maximum
#' spatio-temporal distance allowed for an observation to be considered
#' within a spatio-temporal neighbourhood
#' @param k A numeric indicating the temporal embed size \(number\)
#' @param ratios2add A vector of Boolean values indicating, for each
#' statistic in \code{stats} whether ratios between neighorhoods of 
#' subsequent sizes should be included as extra columns
#' @param .verbose A Boolean indicating whether stages should be printed 
#' out. Defaults to TRUE
#' @inheritParams get_st_indicators
#'
#' @return A data frame that contains extra columns <var>_Tm1, <var>_Tm2, 
#' ..., <var>_Tm<k-1> with previous observations for each location,
#' summary statistics of the values of \code{var} found within the 
#' spatio-temporal neighbourhoods of the one or more \code{radiuses} 
#' of each pair (location ID, time-stamp) and ratios between them
#' 
#' 
#' @import dplyr
#' 
#' @export
get_full_indicators <- function(df, stations, k, 
                                betas, alpha=0.5, var="value",
                                stats = c("mean", "weighted.mean", "sd"), 
                                ratios2add = c(TRUE,TRUE,FALSE),
                                neib_type="cone",
                                .parallel=FALSE, nsplits=1,
                                time_id="time", site_id="station",
                                .verbose=TRUE){
  
  if(.verbose) cat("Embedding series...\t")
  emb_df <- embed_series(df, var, k=k)
  
  if(.verbose) cat("Calculating spatio-temporal indicators...\n")
  # get spatio-temporal indicators
  ind_df <- get_st_indicators(df,
                              stations_sf = stations,
                              radiuses = betas, 
                              stats = stats, 
                              alpha=alpha, 
                              neib_type=neib_type,
                              time_id=time_id, site_id=site_id,
                              vars=var, parallel=.parallel, nsplits=nsplits)
  
  if(.verbose) cat("Adding indicator ratios\n")
  # Add indicator ratios...
  # create ratio variables
  for(i in 1:length(stats)){
    if(ratios2add[i])
      ind_df <- add_ratios(ind_df, var, indStat=stats[i])
  }
  
  # Join with original data set
  df[[time_id]] <- as.character(df[[time_id]])
  df[[site_id]] <- as.character(df[[site_id]])
  
  ind_df[[time_id]] <- as.character(ind_df[[time_id]])
  ind_df[[site_id]] <- as.character(ind_df[[site_id]])
  
  emb_df[[time_id]] <- as.character(emb_df[[time_id]])
  emb_df[[site_id]] <- as.character(emb_df[[site_id]])
  
  ind_df <- dplyr::left_join(emb_df, ind_df, by=c(site_id, time_id)) %>%
    as.data.frame()
  ind_df <- dplyr::left_join(ind_df, df, by=c(site_id, time_id, var))

  if(.verbose) cat("Done!\n")
  ind_df
  
}

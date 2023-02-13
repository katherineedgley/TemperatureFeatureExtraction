library(imputeTS)
library(data.table)
library(zoo)
library(pracma)



#' Main function to extract temperature features
#' @param filename ".bin" filename of file that we want to extract temperature features from
#' @param output_dir output directory where GGIR outputs are stored in
#' Currently this is based on version 1.9-2 output. The two files that are accessed for this
#' processing are in "output_dir/meta/basic/" and "output_dir/meta/ms4.out"
#' 
#' @return A data.table containing extracted temperature features for actigraphy file
get_temp_features <- function(filename, output_dir) {
  # initialise all variables to numeric NA's
  dt <- data.table(temp.allnight.entropy = as.numeric(NA),
                   tempsd.mean = as.numeric(NA),tempmean.sd =as.numeric(NA),
                   tempsd.sd = as.numeric(NA), temp.min.hrsbeforewake = as.numeric(NA),
                   temp.max.hrsafteronset = as.numeric(NA),
                   temp.max.hrsbeforewake=as.numeric(NA),
                   temp.maxwholeday = as.numeric(NA),
                   temp.minwholeday = as.numeric(NA),temp.minmaxdiff = as.numeric(NA),
                   temp.iqr = as.numeric(NA),temp.night.entropy = as.numeric(NA),
                   temp.mrod.mean = as.numeric(NA),temp.mrod.sd = as.numeric(NA),
                   temp.mroi.mean = as.numeric(NA), temp.mroi.sd = as.numeric(NA))
  
  load(paste0(output_dir, '/meta/basic/meta_', filename, ".RData"))
  load(paste0(output_dir, '/meta/ms4.out/', filename, ".RData"))
  
  metalong<- data.table(M$metalong)
  metalong$timestamp <- as.POSIXct(metalong$timestamp,format='%Y-%m-%dT%H:%M:%OS')
  nightsummary <- data.table(nightsummary)
  nightsummary <- nightsummary[fraction_night_invalid < 0.1]
  nnights <- length(nightsummary$night)
  
  night.table <- data.table(night = 1:nnights, temp.sd = as.numeric(NA),
                            temp.mean = as.numeric(NA), temp.max.hrsafteronset = as.numeric(NA),
                            temp.max.hrsbeforewake = as.numeric(NA), temp.min.hrsbeforewake = as.numeric(NA),
                            temp.maxwholeday = as.numeric(NA),temp.minwholeday = as.numeric(NA),
                            temp.minmaxdiff = as.numeric(NA),temp.iqr = as.numeric(NA),
                            temp.day.entropy = as.numeric(NA),temp.night.entropy = as.numeric(NA),
                            temp.mrod = as.numeric(NA),temp.mroi = as.numeric(NA),
                            temp.sleepwakediff = as.numeric(NA))
  allnights <- c()
  for (n in 1:nnights) {
    print(n)
    return_ls <- extract_night_info(n, nightsummary, metalong, night.table, allnights)
    night.table <- return_ls[[1]]
    allnights <- return_ls[[2]]
  }
  
  i <- 1   # currently this is set to i=1 to extract information from a single file,
  # but this can be modified to iterate through list of files.
  dt$temp.allnight.entropy[[i]] <- (allnights %>% approx_entropy())/length(allnights)
  
  dt$tempsd.mean[[i]] <- mean(night.table$temp.sd, na.rm=TRUE)
  dt$tempmean.sd[[i]] <- sd(night.table$temp.mean, na.rm=TRUE)
  dt$tempsd.sd[[i]] <- sd(night.table$temp.sd, na.rm=TRUE)
  dt$temp.min.hrsbeforewake[[i]] <- 
    mean(night.table$temp.min.hrsbeforewake, na.rm=TRUE)
  dt$temp.max.hrsafteronset[[i]] <- 
    mean(night.table$temp.max.hrsafteronset, na.rm = TRUE)
  dt$temp.maxwholeday[[i]] <-
    mean(night.table$temp.maxwholeday, na.rm=TRUE)
  dt$temp.minwholeday[[i]] <-
    mean(night.table$temp.minwholeday, na.rm=TRUE)
  dt$temp.max.hrsbeforewake[[i]] <- 
    mean(night.table$temp.max.hrsbeforewake, na.rm=TRUE)
  dt$temp.minmaxdiff[[i]] <-
    mean(night.table$temp.minmaxdiff, na.rm = TRUE)
  dt$temp.iqr[[i]] <- mean(night.table$temp.iqr, na.rm=TRUE)
  dt$temp.day.entropy[[i]] <- 
    mean(night.table$temp.day.entropy, na.rm=TRUE)
  dt$temp.night.entropy[[i]] <- 
    mean(night.table$temp.night.entropy, na.rm=TRUE)
  dt$temp.mrod.mean[[i]] <- 
    mean(night.table$temp.mrod, na.rm=TRUE)
  dt$temp.mrod.sd[[i]] <- 
    sd(night.table$temp.mrod, na.rm=TRUE)
  dt$temp.mroi.mean[[i]] <- 
    mean(night.table$temp.mroi, na.rm=TRUE)
  dt$temp.mroi.sd[[i]] <- 
    sd(night.table$temp.mroi, na.rm=TRUE)
  dt$temp.sleepwakediff.mn[[i]] <-
    mean(night.table$temp.sleepwakediff, na.rm=TRUE)
  dt$temp.sleepwakediff.sd[[i]] <-
    sd(night.table$temp.sleepwakediff, na.rm=TRUE)
  return (dt)
}

# takes in start_time and end_time in 24 hour format - as number only, e.g., '24'
#' Function to extract a certain time period from "metalong" data processed by GGIR
#' @param metalong the metalong table extracted by GGIR, located in "/meta/basic/meta_filename" file, 
#' and extracted as M$metalong
#' @param start_date the start date of the beginning of the extract, as string
#' @param end_date the end date of the end of the extract, as string
#' @param start_time time extract begins, as string
#' @param end_time time extract ends, as string
#' 
#' @return extracted sample from full metalong
#' 
#' @example extract_from_metalong(metalong,"2018-05-12", "2018-05-13", '08:00:00','08:00:00')
#' 
extract_from_metalong <- function(metalong, start_date, end_date, 
                                  start_time, end_time) {
  start_time <- paste('T',start_time,sep='')
  end_time <- paste('T',end_time, sep='')
  extraction <-
    metalong[(timestamp >= as.POSIXct(paste(start_date, start_time, sep=''),
                                      format='%Y-%m-%dT%H:%M:%OS')) &
               (timestamp <= as.POSIXct(paste(end_date, end_time,sep=''),
                                        format='%Y-%m-%dT%H:%M:%OS'))]
  return(extraction)
}


#' Helper function to extract all temperature information for use in features from a given
#' night
#' @param n the number of the night in "nightsummary" table 
#' @param nightsummary data.table that can be loaded from GGIR output:
#' "/meta/ms4.out/filename"
#' @param metalong that can be loaded from GGIR output "/meta/basic/meta_filename"
#' @param night.table an empty table with feature names to be filled by this function
#' @param allnights a vector of the temperature across all nights to be passed through
#' each iteration of this function
#' 
#' @return list containing "night.table" filled with temperature measures, and
#' "allnights" containing the vector of temperature each valid night
extract_night_info <- function(n, nightsummary, metalong,
                               night.table, allnights) {
  night_num <- nightsummary$night[n]
  start_date <-
    unique(as.character(as.Date(metalong$timestamp)))[night_num]
  next_date <- 
    unique(as.character(as.Date(metalong$timestamp)))[night_num + 1]
  ubernext_date <-
    unique(as.character(as.Date(metalong$timestamp)))[night_num + 2]
  
  full_night <- 
    extract_from_metalong(metalong, start_date, next_date, '19:00:00','14:00:00')
  full_night[temperaturemean < 24]$temperaturemean <- NA
  
  whole_day <-
    extract_from_metalong(metalong, start_date, next_date, '08:00:00','08:00:00')
  whole_day[temperaturemean < 24]$temperaturemean <- NA
  
  full_48hrs <- 
    extract_from_metalong(metalong, start_date, ubernext_date, '00:00:00','00:00:00')
  full_48hrs[temperaturemean < 24]$temperaturemean <- NA
  
  # avoid errors when datatable is empty
  if (sum(!is.na(full_night$temperaturemean)) > 61 & # 61 is 80% of full data
      sum(!is.na(whole_day$temperaturemean)) > 77) { # 77 is 80% of full 24hrs
    print('interpolate:')
    full_night$temperaturemean <- na_interpolation(full_night$temperaturemean)
    whole_day$temperaturemean <- na_interpolation(whole_day$temperaturemean)
    full_48hrs$temperaturemean <- na_interpolation(full_48hrs$temperaturemean)
    
    night_onset <- if (nightsummary[night==night_num]$sleeplog_onset < 24) {
      as.POSIXct(paste(start_date,
                       'T', nightsummary[night == night_num]$sleeplog_onset_ts, 
                       sep=''), format='%Y-%m-%dT%H:%M:%OS')
    } else {
      as.POSIXct(paste(next_date,
                       'T', nightsummary[night == night_num]$sleeplog_onset_ts, 
                       sep=''), format='%Y-%m-%dT%H:%M:%OS')
    }
    
    night_wake <- if (nightsummary[night==night_num]$sleeplog_wake > 24) {
      as.POSIXct(paste(next_date,
                       'T', nightsummary[night == night_num]$sleeplog_wake_ts, 
                       sep=''), format='%Y-%m-%dT%H:%M:%OS')
    } else {
      as.POSIXct(paste(start_date,
                       'T', nightsummary[night == night_num]$sleeplog_wake_ts, 
                       sep=''), format='%Y-%m-%dT%H:%M:%OS')
    }
    
    night_onset <- lubridate::floor_date(night_onset, '15 minutes')
    night_wake <- lubridate::ceiling_date(night_wake, '15 minutes')
    
    full_48hrs[, sleep := 0]
    full_48hrs[(timestamp >= night_onset &
                  timestamp <= night_wake), sleep := 1]
    full_48hrs <- full_48hrs[nonwearscore == 0]
    
    sleepwakediff <- mean(full_48hrs[sleep==1]$temperaturemean, na.rm=T) -
      mean(full_48hrs[sleep==0]$temperaturemean, na.rm=T)
    
    sleeping_times <- full_night[(timestamp >= night_onset &
                                    timestamp <= night_wake),]
    # moving averages
    tempma <- sleeping_times %>% dplyr::select(timestamp, temp=temperaturemean) %>% 
      mutate(temp.ma = rollmean(temp, k = 7, fill = NA)) %>% as.data.table()
    #plot(x=tempma$timestamp, y=tempma$temp.ma, type='l')
    
    tempma.wholeday <-
      whole_day %>% dplyr::select(timestamp, temp=temperaturemean) %>%
      mutate(temp.ma = rollmean(temp, k=7, fill=NA)) %>% as.data.table()
    #plot(x=tempma.wholeday$timestamp, y=tempma.wholeday$temp.ma, type='l')
    
    
    time.series.day <- as.ts(tempma.wholeday$temp.ma) %>% na.omit()
    day.entropy <- (time.series.day %>% approx_entropy())/nrow(tempma.wholeday)
    
    maxtemp.wholeday <- 
      tempma.wholeday[temp.ma== max(temp.ma, na.rm=TRUE)][1,]
    mintemp.wholeday <- 
      tempma.wholeday[temp.ma== min(temp.ma, na.rm=TRUE)][1,]
    
    # want to calculate min/max whole day differences from 8am, as that is start of day
    # changed from 8am to midnight: 3 Dec 2021
    mintemp_wholeday_diff <- 
      difftime(mintemp.wholeday$timestamp, (as.POSIXct(paste(start_date, 'T08:00:00', sep=''),
                                                       format='%Y-%m-%dT%H:%M:%OS')), units='hours')
    # changed from 8am to midnight: 3 Dec 2021
    maxtemp_wholeday_diff <- 
      difftime(maxtemp.wholeday$timestamp,(as.POSIXct(paste(start_date, 'T08:00:00', sep=''),
                                                      format='%Y-%m-%dT%H:%M:%OS')), units='hours')
    
    night.table$temp.maxwholeday[n] <- maxtemp_wholeday_diff
    night.table$temp.minwholeday[n] <- mintemp_wholeday_diff
    night.table$temp.day.entropy[n] <- day.entropy
    night.table$temp.sleepwakediff[n] <- sleepwakediff
    
    if (nrow(tempma) > 15) {
      time.series.night <- as.ts(tempma$temp.ma) %>% na.omit()
      night.entropy <- (time.series.night %>% approx_entropy())/nrow(tempma)
      allnights <- c(allnights, tempma$temp.ma %>% na.omit())
      
      tempma[,'change.rate' := temp.ma - shift(temp.ma)]
      # maximum rate of decline  - time
      mrod.time <- tempma[change.rate == min(change.rate, na.rm=TRUE)]$timestamp
      # maximum rate of increase: time
      mroi.time <- tempma[change.rate == max(change.rate, na.rm=TRUE)]$timestamp
      
      mrod.time <- difftime(mrod.time, start_date, units='hours')
      # in rare cases there may be > 1 time with same MROD, take first
      if (length(mrod.time) > 1) {
        mrod.time <- mrod.time[1]
      }
      mroi.time <- difftime(mroi.time, start_date, units='hours')
      # in rare cases there may be > 1 time with same MROI, take first
      if (length(mroi.time) > 1) {
        mroi.time <- mroi.time[1]
      }
      
      max.temp <- tempma[temp.ma== max(temp.ma, na.rm=TRUE)][1,]
      min.temp <- tempma[temp.ma== min(temp.ma, na.rm=TRUE)][1,]
      
      
      minmax.diff <- max.temp$temp.ma - min.temp$temp.ma
      
      nightwake_diff <- difftime(night_wake, start_date, units='hours')
      nightonset_diff <- difftime(night_onset, start_date, units='hours')
      maxtemp_diff <- difftime(max.temp$timestamp, start_date, units='hours')
      mintemp_diff <- difftime(min.temp$timestamp, start_date, units='hours')
      
      time.diff.max.onset <- maxtemp_diff - nightonset_diff
      time.diff.max.wake <- nightwake_diff - maxtemp_diff
      time.diff.min <- nightwake_diff - mintemp_diff
      
      temp <- sleeping_times$temperaturemean
      night.table$temp.sd[n] <- sd(temp, na.rm=T)
      night.table$temp.mean[n] <- mean(temp, na.rm=T)
      night.table$temp.max.hrsafteronset[n] <- time.diff.max.onset
      night.table$temp.min.hrsbeforewake[n] <- time.diff.min
      night.table$temp.max.hrsbeforewake[n] <- time.diff.max.wake
      night.table$temp.minmaxdiff[n] <- minmax.diff
      night.table$temp.iqr[n] <- IQR(tempma$temp.ma, na.rm=TRUE)
      night.table$temp.night.entropy[n] <- night.entropy
      night.table$temp.mrod[n] <- mrod.time
      night.table$temp.mroi[n] <- mroi.time
    }
  }
  return(list(night.table, allnights))
}


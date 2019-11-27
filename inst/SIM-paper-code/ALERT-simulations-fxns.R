###  homemade functions for use with ALERT

###choose the largest threshold with peaks captured>85% or
## if peaks captured is <85%, call pick_best_other and success=FALSE
pick_best85 <- function (chosen_vec) {
  chosen_vec2 <- chosen_vec[chosen_vec[, "median.pct.cases.captured"] > 85.0,,drop=FALSE]
  if(nrow(chosen_vec2)<1){
    chosen_vec3 <- pick_best_other(chosen_vec)
    attr(chosen_vec3, "success") <- FALSE
  } else {
  chosen_vec3 <- chosen_vec2[chosen_vec2[, "threshold"] == max(chosen_vec2[, "threshold"], 
                                                               na.rm = TRUE),]
  attr(chosen_vec3, "success") <- TRUE
  }
  return(chosen_vec3)
}

## maximize the min.pct.cases.captured as an alternative
#to median.pct.cases.captured>85%
pick_best_other <- function(chosen_vec){
  chosen_vec2 <- chosen_vec[chosen_vec[, "min.pct.cases.captured"] == max(chosen_vec[, "min.pct.cases.captured"], 
                                                                          na.rm = TRUE),, drop=FALSE]
  chosen_vec3 <- chosen_vec2[chosen_vec2[, "threshold"] == max(chosen_vec2[, "threshold"], 
                                                               na.rm = TRUE),]
  return(chosen_vec3)
}

get_stats <- function (j, params=params, firstMonth=firstMonth_value, minWeeks=minWeeks_value) {
  result <- try(createALERT(j, firstMonth=firstMonth_value, minWeeks=minWeeks_value, target.pct=0.85))
  #print(result$out)
  if (class(result)=='try-error') {
    print("error 1: data must run in createALERT")
  } else {
    alert_stats <- result$out[,1:8]
      best_row <- pick_best85(alert_stats)
    alert_summaries <- c(best_row, params, attr(best_row, "success"))
  }
  return(alert_summaries)
}

### choose a threshold and apply across a dataset.

thresholdtestALERT <- function(data, firstMonth=firstMonth, lag=7, minWeeks=8, whichThreshold=4, k=0, target.pct=NULL, caseColumn='Cases', lastDate=NULL) {
  ## check for correct column headers
  if( !("Date" %in% colnames(data)))
    stop("data needs Date columns.")
  if( !(caseColumn %in% colnames(data)) )
    stop(paste("column named", caseColumn, "not found in data."))
    
  ## subset data if required
  if(!is.null(lastDate))
    data <- subset(data, Date<as.Date(lastDate))
  
  ## create a list where each element of the list contains the indices for rows from that season. 
  years <- unique(year(data$Date))
  idxs <- vector("list", length(years)-1) 
  for(i in 1:length(idxs)) {
    startDate2 <- as.Date(paste0(years[i], "-", firstMonth, "-01"))
    endDate2 <- as.Date(paste0(years[i]+1, "-", firstMonth, "-01"))
    idxs[[i]] <- which(data$Date >= startDate2 & data$Date < endDate2)   
  }
  
  ## threshold to test
  thresholds <- whichThreshold
  
  ## for each threshold and year, calculate metrics
  cnames <- c("threshold",
              "median.dur",
              "median.pct.cases.captured",
              "min.pct.cases.captured",
              "max.pct.cases.captured",
              "pct.peaks.captured",
              "pct.ext.peaks.captured",
              "mean.low.weeks.incl")
  if(!is.null(target.pct)) cnames <- c(cnames, "mean.duration.diff")
  out <- matrix(NA, nrow=length(thresholds), ncol=length(cnames))
  colnames(out) <- cnames
  details <- vector("list", length(thresholds))
  ## run a sample to get dim and dimnames
  samp.num <- ifelse(length(idxs[[1]])==0, 2, 1) # Used for evalALERT, if first season missing (i.e is test season) then use second season for sampleRun
  sampleRun <- special_applyALERT(data[idxs[[samp.num]],], threshold=thresholds[1], k=k, lag=lag, minWeeks=minWeeks, target.pct=target.pct, caseColumn=caseColumn)
  for(i in 1:length(thresholds)){
    tmp <- matrix(NA, nrow=length(idxs), ncol=length(sampleRun))
    colnames(tmp) <- names(sampleRun)
    for(j in 1:length(idxs)){
      if(length(idxs[[j]])==0) next # Used for evalALERT, skips missing (test) season
      tmp[j,] <- special_applyALERT(data[idxs[[j]],], threshold=thresholds[i], k=k, lag=lag, minWeeks=minWeeks, target.pct=target.pct, caseColumn=caseColumn)
    }
    details[[i]] <- tmp
    out[i,"threshold"] <- thresholds[i] ## threshold used
    out[i,"median.dur"] <- median(tmp[,"duration"], na.rm=TRUE) ## median duration
    out[i,"median.pct.cases.captured"] <- round(100*median(tmp[,"ALERT.cases.pct"], na.rm=TRUE),1) ## median % of cases captured
    out[i,"min.pct.cases.captured"] <- round(100*min(tmp[,"ALERT.cases.pct"], na.rm=TRUE),1) ## min % of cases captured
    out[i,"max.pct.cases.captured"] <- round(100*max(tmp[,"ALERT.cases.pct"], na.rm=TRUE),1) ## max % of cases captured
    out[i,"pct.peaks.captured"] <- round(100*sum(tmp[,"peak.captured"])/nrow(tmp),1) ## % of times peak captured
    out[i,"pct.ext.peaks.captured"] <- round(100*sum(tmp[,"peak.ext.captured"])/nrow(tmp),1) ## % of times peak +/- k weeks captured
    out[i,"mean.low.weeks.incl"] <- mean(tmp[,"low.weeks.incl"], na.rm=TRUE)
    if(!is.null(target.pct)) out[i,"mean.duration.diff"] <- mean(tmp[,"duration.diff"], na.rm=TRUE)
  }
  return(list(out=out, details=details))
}


### apply datebased_applyALERT across a dataset to get stats based on user-supplied dates.

datetestALERT <- function(data, CHCOstartdate, CHCOenddate, 
                          firstMonth=8, k = 0, lag = 7, 
                          minWeeks = 8, target.pct = NULL, 
                          plot = FALSE, dateColumn = "Date", caseColumn="Cases") {
  ## check for correct column headers
  if( !("Date" %in% colnames(data)))
    stop("data needs Date columns.")
  if( !(caseColumn %in% colnames(data)) )
    stop(paste("column named", caseColumn, "not found in data."))
  ## create a list where each element of the list contains the indices for rows from that season. 
  years <- unique(year(data$Date))
  idxs <- vector("list", length(years)-1) 
  for(i in 1:length(idxs)) {
    startDate2 <- as.Date(paste0(years[i], "-", firstMonth, "-01"))
    endDate2 <- as.Date(paste0(years[i]+1, "-", firstMonth, "-01"))
    idxs[[i]] <- which(data$Date >= startDate2 & data$Date < endDate2)   
  }
  ##  calculate metrics
  cnames <- c("median.dur",
              "median.pct.cases.captured",
              "min.pct.cases.captured",
              "max.pct.cases.captured",
              "pct.peaks.captured",
              "pct.ext.peaks.captured",
              "mean.low.weeks.incl")
  if(!is.null(target.pct)) cnames <- c(cnames, "mean.duration.diff")
  out <- matrix(NA, nrow=1, ncol=length(cnames))
  colnames(out) <- cnames
  ## run a sample to get dim and dimnames
  samp.num <- ifelse(length(idxs[[1]])==0, 2, 1) # Used for evalALERT, if first season missing (i.e is test season) then use second season for sampleRun
  sampleRun <- datebased_applyALERT(data[idxs[[samp.num]],], CHCOstartdate=CHCOstartdate[1], 
                                    CHCOenddate = CHCOenddate[1], k=k, lag=lag, minWeeks=minWeeks, 
                                    target.pct=target.pct, caseColumn=caseColumn,
                                    dateColumn = dateColumn)
    tmp <- matrix(NA, nrow=length(idxs), ncol=length(sampleRun))
    colnames(tmp) <- names(sampleRun)
      for(j in 1:length(idxs)){
        if(length(idxs[[j]])==0) next # Used for evalALERT, skips missing (test) season
        tmp[j,] <- datebased_applyALERT(data[idxs[[j]],], CHCOstartdate=CHCOstartdate[j], 
                                        CHCOenddate = CHCOenddate[j], k=k, lag=lag, 
                                        minWeeks=minWeeks, target.pct=target.pct, 
                                        caseColumn=caseColumn, dateColumn = dateColumn)
      }
    details <- list()
      details[[i]] <- tmp
    out[1,"median.dur"] <- median(tmp[,"duration"], na.rm=TRUE) ## median duration
    out[1,"median.pct.cases.captured"] <- round(100*median(tmp[,"ALERT.cases.pct"], na.rm=TRUE),1) ## median % of cases captured
    out[1,"min.pct.cases.captured"] <- round(100*min(tmp[,"ALERT.cases.pct"], na.rm=TRUE),1) ## min % of cases captured
    out[1,"max.pct.cases.captured"] <- round(100*max(tmp[,"ALERT.cases.pct"], na.rm=TRUE),1) ## max % of cases captured
    out[1,"pct.peaks.captured"] <- round(100*sum(tmp[,"peak.captured"])/nrow(tmp),1) ## % of times peak captured
    out[1,"pct.ext.peaks.captured"] <- round(100*sum(tmp[,"peak.ext.captured"])/nrow(tmp),1) ## % of times peak +/- k weeks captured
    out[1,"mean.low.weeks.incl"] <- mean(tmp[,"low.weeks.incl"], na.rm=TRUE)
    if(!is.null(target.pct)) out[1,"mean.duration.diff"] <- mean(tmp[,"duration.diff"], na.rm=TRUE)
  return(list(out=out, details=details))
}



#########################################
########   special applyALERT   ############
#########################################

##returns the start and end dates

special_applyALERT <- function (data, threshold, k = 0, lag = 7, minWeeks = 8, target.pct = NULL, 
          plot = FALSE, caseColumn = "Cases") 
{
  if (any(data[, caseColumn] >= threshold) == FALSE) {
    message(paste("In the season starting in", year(data$Date[1]), 
                  "the threshold of", threshold, "was not hit."))
    cnames <- c("tot.cases", "duration", "ALERT.cases", "ALERT.cases.pct", 
                "peak.captured", "peak.ext.captured", "low.weeks.incl", "start", "end",
                "duration.diff")
    out <- rep(0, length(cnames))
    names(out) <- cnames
    out[c("tot.cases")] <- sum(data[, caseColumn])
    out[c("duration.diff")] <- NA
    return(out)
  }
  idxHitDate <- min(which(data[, caseColumn] >= threshold))
  hitDate <- data[idxHitDate, "Date"]
  idxStartDate <- idxHitDate + ceiling(lag/7)
  startDate <- data[idxStartDate, "Date"]
  minEndIdx <- idxStartDate + minWeeks - 1
  if (minEndIdx > nrow(data)) {
    cnames <- c("tot.cases", "duration", "ALERT.cases", "ALERT.cases.pct", 
                "peak.captured", "peak.ext.captured", "low.weeks.incl", "start", "end",
                "duration.diff")
    out <- rep(NA, length(cnames))
    names(out) <- cnames
    out["tot.cases"] <- sum(data[, caseColumn])
    return(out)
  }
  idxEndDate <- NA
  i <- minEndIdx - 1
  while (is.na(idxEndDate)) {
    i <- i + 1
    if (is.na(data[i, caseColumn])) 
      next
    if (data[i, caseColumn] < threshold) 
      idxEndDate <- i
    if (i == nrow(data)) 
      break
  }
  endDate <- data[idxEndDate, "Date"]
  onALERT <- rep(0, nrow(data))
  if (is.na(idxEndDate)) {
    cnames <- c("tot.cases", "duration", "ALERT.cases", "ALERT.cases.pct", 
                "peak.captured", "peak.ext.captured", "low.weeks.incl", "start", "end",
                "duration.diff")
    out <- rep(NA, length(cnames))
    names(out) <- cnames
    out["tot.cases"] <- sum(data[, caseColumn])
    return(out)
  }
  onALERT[idxStartDate:idxEndDate] <- 1
  idxPeak <- min(which(data[, caseColumn] == max(data[, caseColumn], 
                                                 na.rm = TRUE)))
  if (!is.null(target.pct)) {
    postcast <- postcastALERT(data, target.pct, caseColumn = caseColumn)
  }
  cnames <- c("tot.cases", "duration", "ALERT.cases", "ALERT.cases.pct", 
              "peak.captured", "peak.ext.captured", "low.weeks.incl", 
              "start", "end")
  out <- rep(NA, length(cnames))
  names(out) <- cnames
  out["tot.cases"] <- sum(data[, caseColumn])
  #out["duration"] <- idxEndDate - idxStartDate + 1
  out["duration"] <- sum(onALERT)   #TRY this instead
  out["ALERT.cases"] <- sum(data[, caseColumn] * onALERT)
  out["ALERT.cases.pct"] <- out["ALERT.cases"]/out["tot.cases"]
  out["peak.captured"] <- idxPeak >= idxStartDate & idxPeak <= 
    idxEndDate
  out["peak.ext.captured"] <- idxPeak >= (idxStartDate + k) & 
    idxPeak <= (idxEndDate - k)
  out["low.weeks.incl"] <- sum(data[idxStartDate:idxEndDate, 
                                    caseColumn] < threshold)
  out["start"] <- startDate
  out["end"] <- endDate
  if (!is.null(target.pct)) 
    out <- c(out, duration.diff = unname(out["duration"] - 
                                           postcast["duration"]))
  else out <- c(out, duration.diff = NA)
  if (plot) 
    message("Plot option not implemented.")
  return(out)
}


############################################################################
######## use user-entered start and end dates with applyALERT   ############
############################################################################

datebased_applyALERT <- function (data, CHCOstartdate, CHCOenddate, k = 0, lag = 7, minWeeks = 8, target.pct = NULL, 
                                plot = FALSE, dateColumn = "Date", caseColumn="Cases") 
{ 
  idxStartDate <- which(abs(data[, dateColumn]-ymd(CHCOstartdate)) == min(abs(data[, dateColumn]-ymd(CHCOstartdate))))
  startDate <- data[idxStartDate, "Date"]
  minEndIdx <- idxStartDate + minWeeks - 1
  idxEndDate <- which(abs(data[, dateColumn]-ymd(CHCOenddate)) == min(abs(data[, dateColumn]-ymd(CHCOenddate))))
  endDate <- data[idxEndDate, "Date"]
  onALERT <- rep(0, nrow(data))
  if (is.na(idxEndDate)) {
    cnames <- c("tot.cases", "duration", "ALERT.cases", "ALERT.cases.pct", 
                "peak.captured", "peak.ext.captured", "low.weeks.incl", "start", "end",
                "duration.diff")
    out <- rep(NA, length(cnames))
    names(out) <- cnames
    out["tot.cases"] <- sum(data[, caseColumn])
    return(out)
  }
  onALERT[idxStartDate:idxEndDate] <- 1
  idxPeak <- min(which(data[, caseColumn] == max(data[, caseColumn], 
                                                 na.rm = TRUE)))
  if (!is.null(target.pct)) {
    postcast <- postcastALERT(data, target.pct, caseColumn = caseColumn)
  }
  cnames <- c("tot.cases", "duration", "ALERT.cases", "ALERT.cases.pct", 
              "peak.captured", "peak.ext.captured", "low.weeks.incl", 
              "start", "end")
  out <- rep(NA, length(cnames))
  names(out) <- cnames
  out["tot.cases"] <- sum(data[, caseColumn])
  #out["duration"] <- idxEndDate - idxStartDate + 1
  out["duration"] <- sum(onALERT)   #TRY this instead
  out["ALERT.cases"] <- sum(data[, caseColumn] * onALERT)
  out["ALERT.cases.pct"] <- out["ALERT.cases"]/out["tot.cases"]
  out["peak.captured"] <- idxPeak >= idxStartDate & idxPeak <= 
    idxEndDate
  out["peak.ext.captured"] <- idxPeak >= (idxStartDate + k) & 
    idxPeak <= (idxEndDate - k)
  out["low.weeks.incl"] <- sum(data[idxStartDate:idxEndDate, 
                                    caseColumn] ==0)
  out["start"] <- startDate
  out["end"] <- endDate
  if (!is.null(target.pct)) 
    out <- c(out, duration.diff = unname(out["duration"] - 
                                           postcast["duration"]))
  else out <- c(out, duration.diff = NA)
  if (plot) 
    message("Plot option not implemented.")
  return(out)
}


## make a nice figure of selected simulations
sim_compare_figs <- function(start_and_end, selected_sims, index, sims_metadata) {
  ggplot() + 
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    xlim(selected_sims[[index]]$Date[1], selected_sims[[index]]$Date[590]) +
    ylim(-15, 100) +
   # labs(y = "Simulated cases", x = "Threshold-based intervention periods") +
    geom_bar(aes(y=selected_sims[[index]]$Cases, x=selected_sims[[index]]$Date), 
             stat="identity", fill="grey60") +
    geom_vline(aes(xintercept=as.numeric(selected_sims[[index]]$Date[260])), linetype="dashed",
               alpha=0.5, show.legend = FALSE) +
    annotate("text", x = as.Date("2009-05-01"), y = 99, 
             label = paste(sims_metadata[[index]][1], "=", 
                           round(as.numeric(sims_metadata[[index]][2]), digits=4), 
                           ", firstMonth=", sims_metadata[[index]][4], ", threshold=", 
                           sims_metadata[[index]][5], ", simulation=", 
                           sims_metadata[[index]][3], sep="")) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[1], 
                  xmax = start_and_end[[index]]$xstop[1], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[2], 
                  xmax = start_and_end[[index]]$xstop[2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[3], 
                  xmax = start_and_end[[index]]$xstop[3], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[4], 
                  xmax = start_and_end[[index]]$xstop[4], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[5], 
                  xmax = start_and_end[[index]]$xstop[5], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[6], 
                  xmax = start_and_end[[index]]$xstop[6], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[7], 
                  xmax = start_and_end[[index]]$xstop[7], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[8], 
                  xmax = start_and_end[[index]]$xstop[8], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[9], 
                  xmax = start_and_end[[index]]$xstop[9], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[10], 
                  xmax = start_and_end[[index]]$xstop[10], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[[index]]$xstart[11], 
                  xmax = start_and_end[[index]]$xstop[11], 
                  ymin = -15, ymax = -5), alpha = 0.9) 
}


## make a nice figure of real data
real_compare_figs <- function(start_and_end, real_data, fig_label) {
  real_data$Cases[is.na(real_data$Cases)] <- 0
  real_data2 <- arrange(real_data, Date)
  ggplot() + 
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    xlim(real_data2$Date[1], real_data2$Date[590]) +
    ylim(-15, 100) +
    # labs(y = "Simulated cases", x = "Threshold-based intervention periods") +
    geom_bar(aes(y=real_data2$Cases, x=real_data2$Date), 
             stat="identity", fill="grey60") +
    annotate("text", x = as.Date("2009-05-01"), y = 99, 
             label = fig_label) +
    geom_rect(aes(xmin = start_and_end[1,1], 
                  xmax = start_and_end[1,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[2,1], 
                  xmax = start_and_end[2,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[3,1], 
                  xmax = start_and_end[3,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[4,1], 
                  xmax = start_and_end[4,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[5,1], 
                  xmax = start_and_end[5,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[6,1], 
                  xmax = start_and_end[6,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[7,1], 
                  xmax = start_and_end[7,2], 
                  ymin = -15, ymax = -5), alpha = 0.9) +
    geom_rect(aes(xmin = start_and_end[8,1], 
                  xmax = start_and_end[8,2], 
                  ymin = -15, ymax = -5), alpha = 0.9)#
}


  #+
   # geom_rect(aes(xmin = start_and_end[9], 
#                  xmax = start_and_end[9], 
#                  ymin = -15, ymax = -5), alpha = 0.9) +
#    geom_rect(aes(xmin = start_and_end[10], 
#                  xmax = start_and_end[10], 
#                  ymin = -15, ymax = -5), alpha = 0.9) +
#    geom_rect(aes(xmin = start_and_end[11], 
#                  xmax = start_and_end[11], 
#                  ymin = -15, ymax = -5), alpha = 0.9) 
#}

table_pretty <- function(min_val, max_val, median_val){
  char_vec <- paste("[", min_val, ", ", max_val, "]; ", median_val, sep="")
  return(char_vec)
}

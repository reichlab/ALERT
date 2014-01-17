#' Core ALERT functions.
#' Nicholas Reich
#' January 2013

require(lubridate)

#' The createALERT function calculates the ALERT thresholds table given data and other parameters.
#' @param data the historical data to use in the analysis. A data.frame with a "Date" column (must be Date objects) and a  "Cases" column.
#' @param firstMonth month number which is counted as the first month of the 'flu year' 
#' @param lag lag time in days between date of cases and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param allThresholds If TRUE, all integer threshold values between the 10th and 50th percentile are examined. If FALSE, only the 10th, 20th, 30th, 40th, and 50th percentiles are examined.
#' @param k if not NULL, the number of weeks around the peak to evaluate ALERT coverage for
#' @param target.pct optional, can specify the percentage of cases the user is targeting during the ALERT period

#' 
#' @return A matrix summarizing the performance of different ALERT thresholds on the given data. Each row represents a threshold. The columns correspond to (1) threshold used, (2) mean ALERT duration, (3/4/5/6) mean/minimum/maximum/sd of % of cases captured across seasons.
createALERT <- function(data, firstMonth=9, lag=7, minWeeks=8, allThresholds=TRUE, k, target.pct=NULL) {
        ## check for correct column headers
        if( !("Date" %in% colnames(data)) | !("Cases" %in% colnames(data)) )
                stop("data needs Date and Cases columns.")
        
        ## create a list where each element of the list contains the indices for rows from that season. 
        years <- unique(year(data$Date))
        idxs <- vector("list", length(years)-1) 
        for(i in 1:length(idxs)) {
                startDate <- as.Date(paste0(years[i], "-", firstMonth, "-01"))
                endDate <- as.Date(paste0(years[i]+1, "-", firstMonth, "-01"))
                idxs[[i]] <- which(data$Date >= startDate & data$Date < endDate)   
        }
        
        ## calculate thresholds to test
        nonZeroCaseCounts <- data$Cases[which(data$Cases>0)]
        if(allThresholds){
                tmp <- quantile(nonZeroCaseCounts, probs=c(.1, .5))
                thresholds <- unique(seq(tmp[1], tmp[2], by=1))
        } else {
                thresholds <- unique(quantile(nonZeroCaseCounts, probs=seq(.1, .5, by=.1)))
        }
        
        ## for each threshold and year, calculate metrics
        cnames <- c("threshold",
                    "mean.dur",
                    "mean.pct.cases.captured",
                    "min.pct.cases.captured",
                    "max.pct.cases.captured",
                    "sd.pct.cases.captured",
                    "pct.peaks.captured",
                    "pct.ext.peaks.captured",
                    "mean.low.weeks.incl")
        if(!is.null(target.pct)) cnames <- c(cnames, "mean.duration.diff")
        out <- matrix(NA, nrow=length(thresholds), ncol=length(cnames))
        colnames(out) <- cnames
        details <- vector("list", length(thresholds))
        ## run a sample to get dim and dimnames
        sampleRun <- applyALERT(data[idxs[[1]],], threshold=thresholds[1], k=k, lag=lag, minWeeks=minWeeks, target.pct=target.pct)
        for(i in 1:length(thresholds)){
                tmp <- matrix(NA, nrow=length(idxs), ncol=length(sampleRun))
                colnames(tmp) <- names(sampleRun)
                for(j in 1:length(idxs)){
                        tmp[j,] <- applyALERT(data[idxs[[j]],], threshold=thresholds[i], k=k, lag=lag, minWeeks=minWeeks, target.pct=target.pct)
                }
                details[[i]] <- tmp
                out[i,"threshold"] <- thresholds[i] ## threshold used
                out[i,"mean.dur"] <- mean(tmp[,"duration"], na.rm=TRUE) ## mean duration
                out[i,"mean.pct.cases.captured"] <- mean(tmp[,"ALERT.cases.pct"], na.rm=TRUE) ## mean % of cases captured
                out[i,"min.pct.cases.captured"] <- min(tmp[,"ALERT.cases.pct"], na.rm=TRUE) ## min % of cases captured
                out[i,"max.pct.cases.captured"] <- max(tmp[,"ALERT.cases.pct"], na.rm=TRUE) ## max % of cases captured
                out[i,"sd.pct.cases.captured"] <- sd(tmp[,"ALERT.cases.pct"], na.rm=TRUE) ## sd % of cases captured
                out[i,"pct.peaks.captured"] <- sum(tmp[,"peak.captured"])/nrow(tmp) ## % of times peak captured
                out[i,"pct.ext.peaks.captured"] <- sum(tmp[,"peak.ext.captured"])/nrow(tmp) ## % of times peak +/- k weeks captured
                out[i,"mean.low.weeks.incl"] <- mean(tmp[,"low.weeks.incl"], na.rm=TRUE)
                if(!is.null(target.pct)) out[i,"mean.duration.diff"] <- mean(tmp[,"duration.diff"], na.rm=TRUE)
        }
        return(list(out=out, details=details))
}

#' The applyALERT function takes one year of data and a threshold and calculates metrics.
#' @param data a single season of surveillance data
#' @param threshold the ALERT treshold to apply
#' @param k if not NULL, the number of weeks around the peak to evaluate ALERT coverage for
#' @param lag lag time in days between report date and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param target.pct optional, can specify the percentage of cases the user is targeting during the ALERT period
#' @param plot TRUE/FALSE, whether a plot should be generated
#' 
#' @return Vector with the following elements: 
#'      [1] total number of cases for the season, 
#'      [2] duration of the ALERT period
#'      [3] total number of cases in the ALERT period, 
#'      [4] fraction of cases in the ALERT period, 
#'      [5] 1 of peak was captured, 0 otherwise, 
#'      [6] 1 if peak +/- k weeks captured, 0 otherwise, 
#'      [7] 1 if counts rise above threshold for two consecutive weeks at any point after ALERT period, 0 otherwise. 
applyALERT <- function(data, threshold, k=NULL, lag, minWeeks, target.pct=NULL, plot=FALSE) {
        ## find first week where ALERT is hit
        idxHitDate <- min(which(data$Cases>=threshold))
        hitDate <- data[idxHitDate,"Date"]
        idxStartDate <- idxHitDate + ceiling(lag/7)
        startDate <- data[idxStartDate,"Date"]
        
        ## calculate end date
        minEndIdx <- idxStartDate + minWeeks - 1
        if(minEndIdx > nrow(data)) stop("start date occurred too late")
        idxEndDate <- NA
        i <- minEndIdx - 1
        while(is.na(idxEndDate)){
                i <- i + 1
                if(is.na(data$Cases[i])) next
                if(data$Cases[i] < threshold) idxEndDate <- i
                if(i==nrow(data)) break
        }
        endDate <- data[idxEndDate, "Date"]
        
        ## make 0/1 vector for ALERT period
        onALERT <- rep(0, nrow(data))
        onALERT[idxStartDate:idxEndDate] <- 1
        
        ## get peak idx
        idxPeak <- min(which(data$Cases==max(data$Cases, na.rm=TRUE))) ## could be multiple...
        
        ## run "postcasting" analysis
        ## 1: what is the maximum number of cases captured in an X-week period, where X is the same duration as the ALERT period for this year.
        ## 2: what is the shortest duration that captures X% of cases, where X is target.pct
        if(!is.null(target.pct)) {
                postcast <- postcastALERT(data, target.pct)
        }
        
        ## define output
        cnames <- c("tot.cases",
                    "duration",
                    "ALERT.cases",
                    "ALERT.cases.pct",
                    "peak.captured",
                    "peak.ext.captured",
                    "low.weeks.incl")
        out <- rep(NA, length(cnames))
        names(out) <- cnames
        out["tot.cases"] <- sum(data$Cases) ## total cases for season
        out["duration"] <- idxEndDate - idxStartDate + 1 ## duration of ALERT period
        out["ALERT.cases"] <- sum(data$Cases*onALERT) ## total number of cases in ALERT period
        out["ALERT.cases.pct"] <- out["ALERT.cases"]/out["tot.cases"] ## fraction of cases in ALERT period
        out["peak.captured"] <- idxPeak>=idxStartDate & idxPeak<=idxEndDate ## peak captured
        out["peak.ext.captured"] <- idxPeak>=(idxStartDate+k) & idxPeak<=(idxEndDate-k) ## peak +/- k weeks captured
        out["low.weeks.incl"]  <- sum(data[idxStartDate:idxEndDate, "Cases"] < threshold)
        if(!is.null(target.pct)) out <- c(out, duration.diff=unname(out["duration"]-postcast["duration"])) 
        if(plot) message("Plot option not implemented.")
        
        return(out)
}


#' The validateALERT function takes a full dataset and parameters to both create and apply the ALERT algorithm. The parameter minPercent also enalbes automatic choice of threshold. For each season in the dataset, an ALERT threshold is calculated leaving that year out, then the resulting threshold is applied to that year. The metrics are saved and summarized.
#'
#' @param data the historical data to use in the analysis
#' @param firstMonth month number which is counted as the first month of the 'flu year' 
#' @param lag lag time between report date and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param minPercent specify the minimum percent of cases to be captured by ALERT. This enables automated threshold selection.
#' @param k the number of weeks around the peak to evaluate ALERT coverage for
validateALERT <- function(data, firstMonth, lag, minWeeks, minPercent) {
        ## check for correct column headers
        if( !("Date" %in% colnames(data)) | !("Cases" %in% colnames(data)) )
                stop("Data needs Date and Cases columns.")

        ## for each season, run createALERT on other data, the applyALERT to the left out year. 
        ## get years where we have data for firstMonth
        years <- unique(year(data[month(data$Date)==firstMonth, "Date"]))        
        
        ## create a list where each element of the list contains the indices for rows from that season. 
        idxs <- vector("list", length(years)) 
        for(i in 1:length(idxs)) {
                startDate <- as.Date(paste0(years[i], "-", firstMonth, "-01"))
                endDate <- as.Date(paste0(years[i]+1, "-", firstMonth, "-01"))
                idxs[[i]] <- which(data$Date >= startDate & data$Date < endDate)   
        }
        
        ## store and report metrics. 
}

#' The postcastALERT function is meant to be called from within the applyALERT function and is used to find the optimal window of time that contains a given percentage of cases.
#' 
#' @param data a single season of surveillance data
#' @param target.pct specifies the percentage of cases the user is targeting during the ALERT period
#' 
#' @return A vector with two elements, the "pct.captured" and "duration" of the interval

postcastALERT <- function(data, target.pct) {
        totalCases <- sum(data$Cases)
        target.cases <- ceiling(totalCases*target.pct)
        minWeeksNeeded <- NA
        nWeeks <- 0
        while(is.na(minWeeksNeeded)) {
                nWeeks <- nWeeks + 1
                ## find all nWeek windows
                weekMatrix <- matrix(NA, ncol=nWeeks, nrow=nrow(data)-nWeeks+1)
                obsCases <- rep(NA, nrow(data)-nWeeks+1)
                for(i in 1:nWeeks) {
                        weekMatrix[,i] <- i:(nrow(data)-nWeeks+i)       
                }
                ## find all numbers of cases in each window
                for(i in 1:nrow(weekMatrix)){
                        obsCases[i] <- sum(data[weekMatrix[i,], "Cases"])       
                }
                ## if any over target.cases --> minWeeksNeeded <- nWeeks
                if(any(obsCases>=target.cases)) {
                        minWeeksNeeded <- nWeeks
                        pct.captured <- max(obsCases/totalCases)
                }
                
        }
        return(c(pct.captured=pct.captured, duration=minWeeksNeeded))
}


#' The plotALERT function can be called from the applyALERT function and is used to plot the details of the application of a single year of the ALERT algorithm to data.
#' 
#' @param data a single season of surveillance data
#' @param threshold the ALERT treshold to apply

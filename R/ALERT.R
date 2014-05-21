#' ALERT test data
#' 
#' Weekly influenza cases from a hospital
#' 
#' @name fluData
#' @docType data
#' @format A weekly time series.
#' @keywords fluData
#' @examples
#' 
#' data(fluData)
#' 
#' ## project the proper start date for the upcoming flu season
#' createALERT(data=fluData, firstMonth=8, lag=7, minWeeks=8, allThresholds=TRUE, k=2, target.pct=0.85)
#' 
#' ## cross validate the above projection
#' robustALERT(minPercent=c(.8, .85, .9), maxDuration=c(12, 13, 14), data=fluData, firstMonth=8, lag=7, minWeeks=8, allThresholds=TRUE, k=2, target.pct=0.85)
#' 
data(fluData)

#' Saved output from robustALERT()
#' 
#' Summarized performance of multiple ALERT rules applied to the fluData dataset. 
#' Desiged for quick loading for vignette.
#' 
#' @name alert_eval
#' @docType data
#' @format Data frame with the output from robustALERT used in the vignette.
#' @keywords alert_eval
#' @examples
#' 
#' data(alert_eval)
data(alert_eval)


#' Producing the ALERT thresholds table
#' 
#' The \code{createALERT} function calculates the ALERT thresholds table given data and other parameters.
#' 
#' @aliases createALERT
#' @param data the historical data to use in the analysis. A data frame with a 'Date' column (must be \code{Date} objects) and a 'Cases' column.
#' @param firstMonth month number which is counted as the first month of the 'flu year' 
#' @param lag lag time in days between date of cases and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param allThresholds if \code{TRUE}, all integer threshold values between the 10th and 50th percentile are examined. If \code{FALSE}, only the 10th, 20th, 30th, 40th, and 50th percentiles are examined.
#' @param k the number of weeks around the peak to evaluate ALERT coverage for
#' @param target.pct the percentage of cases the user is targeting during the ALERT period (optional)
#' @param caseColumn the name of the column with the case counts in it. Defaults to 'Cases'
#' @param lastDate a string in unambigous date format. All cases after this date will be removed prior to running the createALERT algorithm
#' @return By utilizing the \code{\link{applyALERT}} function, \code{createALERT} uses prior hospital data to prospectively determine the start and end to a period of elevated influenza incidence in a community.
#' 
#' \code{createALERT()$details} creates a list of matrices. Each matrix contains raw statistics for the performance of a threshold for each flu season in \code{data} (these statistics can be found in \code{\link{applyALERT})}.
#' 
#' \code{createALERT()$out} produces a matrix summarizing the performance of different ALERT thresholds. The columns for this matrix are as follows:
#' \item{threshold }{the minimum threshold number of cases needed to begin the ALERT period}
#' \item{median.dur }{the median ALERT period duration in weeks}
#' \item{median.pct.cases.captured }{across all seasons, the median percentage of all influenza cases contained within the ALERT period}
#' \item{min.pct.cases.captured }{the minimum percentage of annual cases captured during the ALERT period in any season}
#' \item{max.pct.cases.captured }{the maximum percentage of annual cases captured during the ALERT period in any season}
#' \item{pct.peaks.captured }{the percentage of seasons in which the ALERT period contained the peak week}
#' \item{pct.ext.peaks.captured }{the percentage of seasons in which the ALERT period contained the peak week +/- \code{k} weeks}
#' \item{median.low.weeks.incl }{the median number of weeks included in the ALERT period with counts less than \code{threshold}}
#' \item{median.duration.diff }{if \code{target.pct} specified, the median difference between the duration of the ALERT period in a season and the duration of the shortest period needed to capture \code{target.pct} of cases for that season}
#' @note %% ~~further notes~~
#' @author Nicholas G Reich and Stephen A Lauer
#' @seealso \code{\link{evalALERT}} cross-validates ALERT data over a rule
#' 
#' \code{\link{robustALERT}} cross-validates ALERT data over a series of rules
#' @references %% ~put references to the literature/web site here ~
#' @keywords createALERT
#' @examples
#' 
#' ## Find the ALERT thresholds table for fluData over all thresholds
#' data(fluData)
#' x <- createALERT(data=fluData, allThresholds=TRUE, k=2, target.pct=0.85)
#' x$out

createALERT <- function(data, firstMonth=9, lag=7, minWeeks=8, allThresholds=FALSE, k=0, target.pct=NULL, caseColumn='Cases', lastDate=NULL) {
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
                startDate <- as.Date(paste0(years[i], "-", firstMonth, "-01"))
                endDate <- as.Date(paste0(years[i]+1, "-", firstMonth, "-01"))
                idxs[[i]] <- which(data$Date >= startDate & data$Date < endDate)   
        }
        
        ## calculate thresholds to test
        nonZeroCaseCounts <- data[which(data[,caseColumn]>0), caseColumn]
        if(allThresholds){
                tmp <- quantile(nonZeroCaseCounts, probs=c(.1, .6))
                thresholds <- unique(seq(ceiling(tmp[1]), tmp[2], by=1))
        } else {
                thresholds <- unique(ceiling(quantile(nonZeroCaseCounts, probs=seq(.1, .6, by=.1))))
        }
        
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
        sampleRun <- applyALERT(data[idxs[[samp.num]],], threshold=thresholds[1], k=k, lag=lag, minWeeks=minWeeks, target.pct=target.pct, caseColumn=caseColumn)
        for(i in 1:length(thresholds)){
                tmp <- matrix(NA, nrow=length(idxs), ncol=length(sampleRun))
                colnames(tmp) <- names(sampleRun)
                for(j in 1:length(idxs)){
                        if(length(idxs[[j]])==0) next # Used for evalALERT, skips missing (test) season
                        tmp[j,] <- applyALERT(data[idxs[[j]],], threshold=thresholds[i], k=k, lag=lag, minWeeks=minWeeks, target.pct=target.pct, caseColumn=caseColumn)
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

#' Producing seasonal ALERT output
#' 
#' The \code{applyALERT} function (often called by \code{\link{createALERT}} or \code{\link{evalALERT}}) takes one year of data and a threshold and calculates metrics.
#' 
#' @aliases applyALERT
#' @param data a single season of surveillance data
#' @param threshold the ALERT threshold to apply
#' @param k the number of weeks around the peak to evaluate ALERT coverage for
#' @param lag lag time in days between report date and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param target.pct the percentage of cases the user is targeting during the ALERT period (optional)
#' @param caseColumn the name of the column with the case counts in it. Defaults to 'Cases'
#' @param plot \code{TRUE}/\code{FALSE}, whether a plot should be generated (currently unavailable)
#' 
#' @return Returns a vector with the following elements: 
#'      \item{tot.cases }{total number of cases for the season} 
#'      \item{duration }{duration of the ALERT period}
#'      \item{ALERT.cases }{total number of cases in the ALERT period}
#'      \item{ALERT.cases.pct }{fraction of cases in the ALERT period}
#'      \item{peak.captured }{1 if peak was captured, 0 otherwise}
#'      \item{peak.ext.captured }{1 if peak +/- \code{k} weeks captured, 0 otherwise}
#'      \item{low.weeks.incl }{the number of weeks included in the ALERT period with counts less than \code{threshold}}
#'      \item{duration.diff }{if \code{target.pct} specified, the difference between the duration of the ALERT period and the duration of the shortest period needed to capture \code{target.pct} using \code{\link{postcastALERT}}.}
#'      
#' @note %% ~~further notes~~
#' @author Nicholas G Reich and Stephen A Lauer
#' @seealso \code{\link{createALERT}}, \code{\link{evalALERT}}, \code{\link{robustALERT}}
#' @references %% ~put references to the literature/web site here ~
#' @keywords applyALERT
#' @examples 
#' 
#' ## Find the ALERT metrics of a season with a threshold of 3
#' data(fluData)
#' applyALERT(data=fluData, threshold=3, k=2, target.pct=0.85)

applyALERT <- function(data, threshold, k=0, lag=7, minWeeks=8, target.pct=NULL, plot=FALSE, caseColumn='Cases') {
        ## confirm that ALERT threshold is hit in this year
        if(any(data[,caseColumn]>=threshold)==FALSE) {
                message(paste("In the season starting in", year(data$Date[1]), 
                              "the threshold of", threshold, "was not hit."))
                ## if no week exceeds threshold, retun empty vector
                cnames <- c("tot.cases",
                            "duration",
                            "ALERT.cases",
                            "ALERT.cases.pct",
                            "peak.captured",
                            "peak.ext.captured",
                            "low.weeks.incl",
                            "duration.diff")
                out <- rep(0, length(cnames)) ## almost everything set to zero
                names(out) <- cnames
                out[c("tot.cases")] <- sum(data[,caseColumn])
                out[c("duration.diff")] <- NA ## leave out duration diffs when ALERT is not hit.
                return(out)
        }
        
        ## find first week where ALERT is hit
        idxHitDate <- min(which(data[,caseColumn]>=threshold))
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
                if(is.na(data[i, caseColumn])) next
                if(data[i, caseColumn] < threshold) idxEndDate <- i
                if(i==nrow(data)) break
        }
        endDate <- data[idxEndDate, "Date"]
        
        ## make 0/1 vector for ALERT period
        onALERT <- rep(0, nrow(data))
        if(is.na(idxEndDate)) stop("start date occurred too late")
        onALERT[idxStartDate:idxEndDate] <- 1
        
        ## get peak idx
        idxPeak <- min(which(data[,caseColumn]==max(data[,caseColumn], na.rm=TRUE))) ## could be multiple...
        
        ## run "postcasting" analysis
        ## 1: what is the maximum number of cases captured in an X-week period, where X is the same duration as the ALERT period for this year.
        ## 2: what is the shortest duration that captures X% of cases, where X is target.pct
        if(!is.null(target.pct)) {
                postcast <- postcastALERT(data, target.pct, caseColumn=caseColumn)
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
        out["tot.cases"] <- sum(data[,caseColumn]) ## total cases for season
        out["duration"] <- idxEndDate - idxStartDate + 1 ## duration of ALERT period
        out["ALERT.cases"] <- sum(data[,caseColumn]*onALERT) ## total number of cases in ALERT period
        out["ALERT.cases.pct"] <- out["ALERT.cases"]/out["tot.cases"] ## fraction of cases in ALERT period
        out["peak.captured"] <- idxPeak>=idxStartDate & idxPeak<=idxEndDate ## peak captured
        out["peak.ext.captured"] <- idxPeak>=(idxStartDate+k) & idxPeak<=(idxEndDate-k) ## peak +/- k weeks captured
        out["low.weeks.incl"]  <- sum(data[idxStartDate:idxEndDate, caseColumn] < threshold)
        if(!is.null(target.pct)) out <- c(out, duration.diff=unname(out["duration"]-postcast["duration"]))
        else out <- c(out, duration.diff=NA)
        if(plot) message("Plot option not implemented.")
        
        return(out)
}


#' Cross-validating ALERT output under a rule
#' 
#' The \code{evalALERT} function uses the ALERT algorithm to test a rule (either \code{minPercent} or \code{maxDuration}) on hospital influenza data. For each season in the dataset, \code{\link{createALERT}} finds an optimal ALERT \code{threshold} when leaving that year out. Then \code{\link{applyALERT}} tests that \code{threshold} in that year. The metrics are saved and summarized.
#' 
#' @aliases evalALERT
#' @param data the historical data to use in the analysis. A data frame with a 'Date' column (must be \code{Date} objects) and a 'Cases' column.
#' @param minPercent specify the minimum percent of cases to be captured at least 50\% of the time by ALERT. This enables automated threshold selection.
#' @param maxDuration specify the maximum number of weeks to be captured at least 50\% of the time by ALERT. This enables automated threshold selection.
#' @param firstMonth month number which is counted as the first month of the 'flu year' 
#' @param lag lag time between report date and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param allThresholds if \code{TRUE}, all integer threshold values between the 10th and 50th percentile are examined. If \code{FALSE}, only the 10th, 20th, 30th, 40th, and 50th percentiles are examined.
#' @param k the number of weeks around the peak to evaluate ALERT coverage for
#' @param target.pct the percentage of cases the user is targeting during the ALERT period when
#' @param caseColumn the name of the column with the case counts in it. Defaults to 'Cases'
#'  testing \code{maxDuration} (optional)
#' @return Returns a table with the following columns: 
#'      \item{season }{each flu season in the \code{data} that was able to be evaluated for the given rule, with the final row reserved for summary statistics}
#'      \item{threshold }{the minimum threshold number of cases needed to begin the ALERT period}
#'      \item{tot.cases }{total number of cases for the season} 
#'      \item{duration }{duration of the ALERT period}
#'      \item{ALERT.cases }{total number of cases in the ALERT period}
#'      \item{ALERT.cases.pct }{fraction of cases in the ALERT period}
#'      \item{peak.captured }{1 if peak was captured, 0 otherwise}
#'      \item{peak.ext.captured }{1 if peak +/- \code{k} weeks captured, 0 otherwise}
#'      \item{low.weeks.incl }{the number of weeks included in the ALERT period with counts less than \code{threshold}}
#'      \item{duration.diff }{if \code{target.pct} specified, the difference between the duration of the ALERT period and the duration of the shortest period needed to capture \code{target.pct} using \code{\link{postcastALERT}}.}
#' @return Each row in the table represents a season from the data. The final row outputs summary statistics: the median of \code{threshold}, \code{tot.cases}, \code{duration}, \code{ALERT.cases}, and \code{ALERT.cases.pct} and the mean of \code{peak.captured}, \code{peak.ext.captured}, \code{low.weeks.incl}, and \code{duration.diff}.
#'      
#' @note %% ~~further notes~~
#' @author Nicholas G Reich and Stephen A Lauer
#' @seealso \code{\link{robustALERT}} cross-validates ALERT output under a set of rules
#' @references %% ~put references to the literature/web site here ~
#' @keywords evalALERT
#' @examples
#' 
#' ## find the highest threshold captures at least 85% half of the time
#' data(fluData)
#' evalALERT(minPercent=.85, data=fluData, k=2)
#' 
#' ## find the lowest threshold that has a median duration of less than 12 weeks
#' evalALERT(maxDuration=12, data=fluData, k=2)

evalALERT <- function(data, minPercent=NULL, maxDuration=NULL, firstMonth=9, lag=7, minWeeks=8, allThresholds=FALSE, k=0, target.pct=NULL, caseColumn='Cases') {
        if(is.null(maxDuration) & is.null(minPercent))
                stop("Please choose a rule to evaluate, either with maxDuration or minPercent.")
        ## check for correct column headers
        if( !("Date" %in% colnames(data)))
                stop("data needs Date columns.")
        if( !(caseColumn %in% colnames(data)) )
                stop(paste("column named", caseColumn, "not found in data."))
        
        ## create test matrix
        eval.dat <- matrix(data=0, nrow=0, ncol=10)
        
        ## get years where we have data for firstMonth
        years <- unique(year(data[month(data$Date)==firstMonth, "Date"]))
        ## create a list where each element of the list contains the indices for rows from that season.
        idxs <- vector("list", length(years)) 
        for(i in 1:length(idxs)) {
                startDate <- as.Date(paste0(years[i], "-", firstMonth, "-01"))
                endDate <- as.Date(paste0(years[i]+1, "-", firstMonth, "-01"))
                idxs[[i]] <- which(data$Date >= startDate & data$Date < endDate)   
        }
        ## for each season, run createALERT on other data, the applyALERT to the left out year.
        if(!is.null(minPercent)){
                for(j in 1:length(years)){
                        if(length(idxs[[j]])<=3*minWeeks) next # if a season is too short to evaluate, skip it
                        ## leave one season out of data
                        data1 <- data[-idxs[[j]],]
                        ## run createALERT on other data
                        output <- as.data.frame(createALERT(data=data1, firstMonth=firstMonth, 
                                                            lag=lag, minWeeks=minWeeks, allThresholds=allThresholds, 
                                                            k=k, target.pct=minPercent,
                                                            caseColumn=caseColumn)$out)
                        output1 <- subset(output, median.pct.cases.captured>=minPercent*100)
                        ## find largest threshold that achieves target percent covered
                        if(length(output1[,1])==0){
                                print(paste0("The median captured percentage for each threshold in the year ", years[j], " fell below the minimum percentage of ", minPercent, "."))
                                next # skips any years where target percentage is not achieved
                        }
                        opt.thresh <- max(output1$threshold)
                        ## run applyALERT on left out year with threshold determined above
                        if(max(data[idxs[[j]],2])<opt.thresh) next # skips any years where threshold is not achieved (usually partial season at beginning or end of data)
                        aaa <- applyALERT(data[idxs[[j]],], threshold=opt.thresh, k=k, lag=lag, 
                                          minWeeks=minWeeks, target.pct=minPercent, 
                                          caseColumn=caseColumn, plot=FALSE)
                        ## store metrics
                        aaa <- c(season=years[j], threshold=opt.thresh, aaa)
                        eval.dat <- rbind.data.frame(eval.dat, aaa)
                }
                colnames(eval.dat) <- names(aaa)
        }
        
        if(!is.null(maxDuration)){
                for(j in 1:length(years)){
                        if(length(idxs[[j]])<=3*minWeeks) next # if a season is too short to evaluate, skip it
                        ## leave one season out of data
                        data1 <- data[-idxs[[j]],]
                        ## run createALERT on other data
                        output <- as.data.frame(createALERT(data=data1, firstMonth=firstMonth, lag=lag, minWeeks=minWeeks, allThresholds=allThresholds, k=k, caseColumn=caseColumn, target.pct=target.pct)$out)
                        output1 <- subset(output, median.dur<=maxDuration)
                        ## find smallest threshold that has a duration shorter than maxDuration
                        if(length(output1[,1])==0){
                                print(paste("The median durations for each threshold in the year", years[j], "exceeded", maxDuration, "weeks."))
                                next # skips any years where all durations are too long
                        }
                        opt.thresh <- min(output1$threshold)
                        ## run applyALERT on left out year with threshold determined above
                        if(max(data[idxs[[j]],2])<opt.thresh) next # skips any years where threshold is not achieved (usually partial season at beginning or end of data)
                        aaa <- applyALERT(data[idxs[[j]],], threshold=opt.thresh, k=k, lag=lag, 
                                          minWeeks=minWeeks, target.pct=target.pct, 
                                          caseColumn=caseColumn, plot=FALSE)
                        ## store metrics
                        aaa <- c(season=years[j], threshold=opt.thresh, aaa)
                        eval.dat <- rbind.data.frame(eval.dat, aaa)
                }
                colnames(eval.dat) <- names(aaa)
        }
        
        ## take summary statistics for final row
        bbb <- c("Summaries", round(median(eval.dat$threshold),1), 
                 round(median(eval.dat$tot.cases),1), 
                 round(median(eval.dat$duration),1), 
                 round(median(eval.dat$ALERT.cases),1), 
                 round(median(eval.dat$ALERT.cases.pct),3), 
                 round(mean(eval.dat$peak.captured),3), 
                 round(mean(eval.dat$peak.ext.captured),3), 
                 round(mean(eval.dat$low.weeks.incl),1), 
                 round(mean(eval.dat$duration.diff),1))
        eval.dat <- rbind.data.frame(eval.dat, bbb)
        return(eval.dat)
}

#' Cross-validating ALERT over a set of rules.
#' 
#' The \code{robustALERT} function finds the optimal threshold for starting an ALERT season for a vector of rules. For each rule specified, \code{\link{evalALERT}} tests the rule against every flu season for the given \code{data} and outputs summary statistics. \code{robustALERT} aggregates these statistics together for easy comparison. This function can be used to validate the results from \code{\link{createALERT}}.
#' 
#' @aliases robustALERT
#' @param data the historical data to use in the analysis. A data.frame with a
#' "Date" column (must be Date objects) and a "Cases" column.
#' @param minPercent value or vector that specifies the minimum percent of
#' cases to be captured at least 50\% of the time by ALERT. This enables automated threshold
#' selection.
#' @param maxDuration value or vector that specifies the maximum number of
#' weeks to be captured at least 50\% of the time by ALERT. This enables automated threshold
#' selection.
#' @param firstMonth firstMonth month number which is counted as the first
#' month of the 'flu year'
#' @param lag lag time in days between date of cases and action taken
#' @param minWeeks minimum number of weeks to be in ALERT
#' @param allThresholds If TRUE, all integer threshold values between the 10th
#' and 50th percentile are examined. If FALSE, only the 10th, 20th, 30th, 40th,
#' and 50th percentiles are examined.
#' @param k if not 0, the number of weeks around the peak to evaluate ALERT
#' coverage for
#' @param target.pct can specify the percentage of cases the user is targeting
#' during the ALERT period when testing maxDuration (optional)
#' @param caseColumn the name of the column with the case counts in it. Defaults to 'Cases'
#' @return A table of the median threshold and ALERT results determined by each rule with \code{\link{evalALERT}} with the following columns:
#' \item{rule }{each rule that was specified by the user, either by \code{minPercent} or \code{maxDuration}}
#'      \item{threshold }{the median threshold number of cases needed to begin an ALERT period that can satisfy the rule (as determined by cross-validation)}
#'      \item{duration }{the median duration of the ALERT period}
#'      \item{ALERT.cases }{the median number of cases in the ALERT period}
#'      \item{ALERT.cases.pct }{the median fraction of cases in the ALERT period}
#'      \item{peak.captured }{the fraction of the time the peak was captured}
#'      \item{peak.ext.captured }{the fraction of the time the peak +/- \code{k} weeks was captured}
#'      \item{low.weeks.incl }{the mean number of weeks included in the ALERT period with counts less than \code{threshold}}
#'      \item{duration.diff }{if \code{target.pct} specified, the mean difference between the duration of the ALERT period and the duration of the shortest period needed to capture \code{target.pct} using \code{\link{postcastALERT}}.}
#' @note %% ~~further notes~~
#' @author Nicholas G Reich and Stephen A Lauer
#' @seealso \code{\link{createALERT}}
#' @references %% ~put references to the literature/web site here ~
#' @keywords robustALERT
#' @examples
#' 
#' ## view the performance of ALERT over three levels of minimum case percentage and maximum duration length
#' data(fluData)
#' robustALERT(minPercent=c(.8, .85, .9), maxDuration=c(12, 13, 14), data=fluData, k=2, target.pct=0.85)

robustALERT <- function(data, minPercent=NULL, maxDuration=NULL, firstMonth=9, lag=7, minWeeks=8, allThresholds=FALSE, k=0, target.pct=NULL, caseColumn='Cases') {
        ## check for correct column headers
        if( !("Date" %in% colnames(data)))
                stop("data needs Date columns.")
        if( !(caseColumn %in% colnames(data)) )
                stop(paste("column named", caseColumn, "not found in data."))
        
        ## create test matrix
        robust.dat <- matrix(data=0, nrow=0, ncol=10)
        
        ## run evalALERT on all minPercent rules
        if(!is.null(minPercent)){
                for(i in 1:length(minPercent)){
                        one.rule <- evalALERT(data=data, firstMonth=firstMonth, lag=lag, 
                                              minWeeks=minWeeks,allThresholds=allThresholds, 
                                              k=k, caseColumn=caseColumn,
                                              minPercent=minPercent[i])
                        aaa <- one.rule[length(one.rule[,1]),]
                        aaa[1,1] <- paste("minPercent =", minPercent[i])
                        colnames(aaa)[1] <- "rule"
                        #aaa <- cbind.data.frame(aaa[,1], num.years.used=length(one.rule[,1])-1, aaa[,2:10])
                        robust.dat <- rbind.data.frame(robust.dat,aaa)
                }
        }
        
        ## run evalALERT on all maxDuration rules
        if(!is.null(maxDuration)){
                for(j in 1:length(maxDuration)){
                        one.rule <- evalALERT(data=data, firstMonth=firstMonth, lag=lag, 
                                              minWeeks=minWeeks,allThresholds=allThresholds, 
                                              k=k, caseColumn=caseColumn,
                                              target.pct=target.pct, maxDuration=maxDuration[j])
                        aaa <- one.rule[length(one.rule[,1]),]
                        aaa[1,1] <- paste("maxDuration =", maxDuration[j])
                        colnames(aaa)[1] <- "rule"
                        #aaa <- cbind.data.frame(aaa[,1], num.years.used=length(one.rule[,1])-1, aaa[,2:10])
                        robust.dat <- rbind.data.frame(robust.dat,aaa)
                }
        }
        robust.dat[,6] <- round(100*as.numeric(robust.dat[,6]),1)
        robust.dat[,7] <- round(100*as.numeric(robust.dat[,7]),1)
        robust.dat[,8] <- round(100*as.numeric(robust.dat[,8]),1)
        robust.dat <- robust.dat[,-3]
        return(robust.dat)
}

#' A component of applyALERT
#' 
#' The \code{postcastALERT} function is meant to be called from within the \code{applyALERT} function and is used to find the optimal window of time that contains a given percentage of cases.
#' 
#' @aliases postcastALERT
#' @param data a single season of surveillance data
#' @param target.pct specifies the percentage of cases the user is targeting during the ALERT period
#' @param caseColumn the name of the column with the case counts in it. Defaults to 'Cases'
#' 
#' @return A vector with two elements, the \code{pct.captured} and \code{duration} of the interval
#' 
#' @note %% ~~further notes~~
#' @author Nicholas G Reich and Stephen A Lauer
#' @seealso \code{\link{applyALERT}}, \code{\link{createALERT}}
#' @references %% ~put references to the literature/web site here ~
#' @keywords postcastALERT
#' @examples
#' 
#' data(fluData)
#' postcastALERT(fluData, target.pct=.85)

postcastALERT <- function(data, target.pct, caseColumn='Cases') {
        totalCases <- sum(data[,caseColumn])
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
                        obsCases[i] <- sum(data[weekMatrix[i,], caseColumn])       
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
#' @param threshold the ALERT threshold to apply

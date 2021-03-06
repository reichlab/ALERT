% Generated by roxygen2 (4.0.1): do not edit by hand
\name{robustALERT}
\alias{robustALERT}
\title{Cross-validating ALERT over a set of rules.}
\usage{
robustALERT(data, minPercent = NULL, maxDuration = NULL, firstMonth = 9,
  lag = 7, minWeeks = 8, allThresholds = FALSE, k = 0,
  target.pct = NULL, caseColumn = "Cases")
}
\arguments{
\item{data}{the historical data to use in the analysis. A data.frame with a
"Date" column (must be Date objects) and a "Cases" column.}

\item{minPercent}{value or vector that specifies the minimum percent of
cases to be captured at least 50\% of the time by ALERT. This enables automated threshold
selection.}

\item{maxDuration}{value or vector that specifies the maximum number of
weeks to be captured at least 50\% of the time by ALERT. This enables automated threshold
selection.}

\item{firstMonth}{firstMonth month number which is counted as the first
month of the 'flu year'}

\item{lag}{lag time in days between date of cases and action taken}

\item{minWeeks}{minimum number of weeks to be in ALERT}

\item{allThresholds}{If TRUE, all integer threshold values between the 10th
and 50th percentile are examined. If FALSE, only the 10th, 20th, 30th, 40th,
and 50th percentiles are examined.}

\item{k}{if not 0, the number of weeks around the peak to evaluate ALERT
coverage for}

\item{target.pct}{can specify the percentage of cases the user is targeting
during the ALERT period when testing maxDuration (optional)}

\item{caseColumn}{the name of the column with the case counts in it. Defaults to 'Cases'}
}
\value{
A table of the median threshold and ALERT results determined by each rule with \code{\link{evalALERT}} with the following columns:
\item{rule }{each rule that was specified by the user, either by \code{minPercent} or \code{maxDuration}}
     \item{threshold }{the median threshold number of cases needed to begin a sufficient ALERT period}
     \item{duration }{the median duration of the ALERT period}
     \item{ALERT.cases }{the median number of cases in the ALERT period}
     \item{ALERT.cases.pct }{the median fraction of cases in the ALERT period}
     \item{peak.captured }{the fraction of the time the peak was captured}
     \item{peak.ext.captured }{the fraction of the time the peak +/- \code{k} weeks was captured}
     \item{low.weeks.incl }{the mean number of weeks included in the ALERT period with counts less than \code{threshold}}
     \item{duration.diff }{if \code{target.pct} specified, the mean difference between the duration of the ALERT period and the duration of the shortest period needed to capture \code{target.pct} using \code{\link{postcastALERT}}.}
}
\description{
The \code{robustALERT} function finds the optimal threshold for starting an ALERT season for a vector of rules. For each rule specified, \code{\link{evalALERT}} tests the rule against every flu season for the given \code{data} and outputs summary statistics. \code{robustALERT} aggregates these statistics together for easy comparison. This function can be used to validate the results from \code{\link{createALERT}}.
}
\note{
%% ~~further notes~~
}
\examples{
## view the performance of ALERT over three levels of minimum case percentage and maximum duration length
data(fluData)
robustALERT(minPercent=c(.8, .85, .9), maxDuration=c(12, 13, 14), data=fluData, k=2, target.pct=0.85)
}
\author{
Nicholas G Reich and Stephen A Lauer
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
\code{\link{createALERT}}
}
\keyword{robustALERT}


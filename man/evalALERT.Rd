\name{evalALERT}
\alias{evalALERT}
\title{Cross-validating ALERT output under a rule}
\usage{
evalALERT(data, minPercent = NULL, maxDuration = NULL, firstMonth = 9,
  lag = 7, minWeeks = 8, allThresholds = FALSE, k = 0,
  target.pct = NULL)
}
\arguments{
  \item{data}{the historical data to use in the analysis. A
  data frame with a 'Date' column (must be \code{Date}
  objects) and a 'Cases' column.}

  \item{minPercent}{specify the minimum percent of cases to
  be captured at least 50\% of the time by ALERT. This
  enables automated threshold selection.}

  \item{maxDuration}{specify the maximum number of weeks to
  be captured at least 50\% of the time by ALERT. This
  enables automated threshold selection.}

  \item{firstMonth}{month number which is counted as the
  first month of the 'flu year'}

  \item{lag}{lag time between report date and action taken}

  \item{minWeeks}{minimum number of weeks to be in ALERT}

  \item{allThresholds}{if \code{TRUE}, all integer
  threshold values between the 10th and 50th percentile are
  examined. If \code{FALSE}, only the 10th, 20th, 30th,
  40th, and 50th percentiles are examined.}

  \item{k}{the number of weeks around the peak to evaluate
  ALERT coverage for}

  \item{target.pct}{the percentage of cases the user is
  targeting during the ALERT period when testing
  \code{maxDuration} (optional)}
}
\value{
Returns a table with the following columns: \item{season
}{each flu season in the \code{data} that was able to be
evaluated for the given rule, with the final row reserved
for summary statistics} \item{threshold }{the minimum
threshold number of cases needed to begin the ALERT period}
\item{tot.cases }{total number of cases for the season}
\item{duration }{duration of the ALERT period}
\item{ALERT.cases }{total number of cases in the ALERT
period} \item{ALERT.cases.pct }{fraction of cases in the
ALERT period} \item{peak.captured }{1 if peak was captured,
0 otherwise} \item{peak.ext.captured }{1 if peak +/-
\code{k} weeks captured, 0 otherwise} \item{low.weeks.incl
}{the number of weeks included in the ALERT period with
counts less than \code{threshold}} \item{duration.diff }{if
\code{target.pct} specified, the difference between the
duration of the ALERT period and the duration of the
shortest period needed to capture \code{target.pct} using
\code{\link{postcastALERT}}.}
}
\description{
The \code{evalALERT} function uses the ALERT algorithm to
test a rule (either \code{minPercent} or
\code{maxDuration}) on hospital influenza data. For each
season in the dataset, \code{\link{createALERT}} finds an
optimal ALERT \code{threshold} when leaving that year out.
Then \code{\link{applyALERT}} tests that \code{threshold}
in that year. The metrics are saved and summarized.
}
\note{
%% ~~further notes~~
}
\examples{
## find the highest threshold captures at least 85\% half of the time
data(fluData)
evalALERT(minPercent=.85, data=fluData, k=2)

## find the lowest threshold that has a median duration of less than 12 weeks
evalALERT(maxDuration=12, data=fluData, k=2)
}
\author{
Nicholas G Reich and Stephen A Lauer
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
\code{\link{robustALERT}} cross-validates ALERT output
under a set of rules
}
\keyword{evalALERT}


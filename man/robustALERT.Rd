\name{robustALERT}
\alias{robustALERT}
\title{Cross-validating ALERT over a set of rules.}
\usage{
robustALERT(minPercent = NULL, maxDuration = NULL, data, firstMonth = 9,
  lag = 7, minWeeks = 8, allThresholds = TRUE, k = 0,
  target.pct = NULL)
}
\arguments{
  \item{minPercent}{value or vector that specifies the
  minimum percent of cases to be captured on average by
  ALERT. This enables automated threshold selection.}

  \item{maxDuration}{value or vector that specifies the
  maximum number of weeks to be captured on average by
  ALERT. This enables automated threshold selection.}

  \item{data}{the historical data to use in the analysis. A
  data.frame with a "Date" column (must be Date objects)
  and a "Cases" column.}

  \item{firstMonth}{firstMonth month number which is
  counted as the first month of the 'flu year'}

  \item{lag}{lag time in days between date of cases and
  action taken}

  \item{minWeeks}{minimum number of weeks to be in ALERT}

  \item{allThresholds}{If TRUE, all integer threshold
  values between the 10th and 50th percentile are examined.
  If FALSE, only the 10th, 20th, 30th, 40th, and 50th
  percentiles are examined.}

  \item{k}{if not 0, the number of weeks around the peak to
  evaluate ALERT coverage for}

  \item{target.pct}{can specify the percentage of cases the
  user is targeting during the ALERT period when testing
  maxDuration (optional)}
}
\value{
A table of the average threshold and ALERT results
determined by each rule with \code{\link{evalALERT}} with
the following columns: \item{rule }{each rule that was
specified by the user, either by \code{minPercent} or
\code{maxDuration}} \item{threshold }{the minimum threshold
number of cases needed to begin the ALERT period}
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
The \code{robustALERT} function finds the optimal threshold
for starting an ALERT season for a vector of rules. For
each rule specified, \code{\link{evalALERT}} tests the rule
against every flu season for the given \code{data} and
outputs summary statistics. \code{robustALERT} aggregates
these statistics together for easy comparison. This
function can be used to validate the results from
\code{\link{createALERT}}.
}
\note{
%% ~~further notes~~
}
\examples{
## view the performance of ALERT over three levels of minimum case percentage and maximum duration length
data(maskData)
robustALERT(minPercent=c(.8, .85, .9), maxDuration=c(12, 13, 14), data=maskData, k=2, target.pct=0.85)
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


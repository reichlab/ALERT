\name{applyALERT}
\alias{applyALERT}
\title{Producing seasonal ALERT output}
\usage{
applyALERT(data, threshold, k = 0, lag = 7, minWeeks = 8,
  target.pct = NULL, plot = FALSE, caseColumn = "Cases")
}
\arguments{
  \item{data}{a single season of surveillance data}

  \item{threshold}{the ALERT threshold to apply}

  \item{k}{the number of weeks around the peak to evaluate
  ALERT coverage for}

  \item{lag}{lag time in days between report date and
  action taken}

  \item{minWeeks}{minimum number of weeks to be in ALERT}

  \item{target.pct}{the percentage of cases the user is
  targeting during the ALERT period (optional)}

  \item{caseColumn}{the name of the column with the case
  counts in it. Defaults to 'Cases'}

  \item{plot}{\code{TRUE}/\code{FALSE}, whether a plot
  should be generated (currently unavailable)}
}
\value{
Returns a vector with the following elements:
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
The \code{applyALERT} function (often called by
\code{\link{createALERT}} or \code{\link{evalALERT}}) takes
one year of data and a threshold and calculates metrics.
}
\note{
%% ~~further notes~~
}
\examples{
## Find the ALERT metrics of a season with a threshold of 3
data(fluData)
applyALERT(data=fluData, threshold=3, k=2, target.pct=0.85)
}
\author{
Nicholas G Reich and Stephen A Lauer
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
\code{\link{createALERT}}, \code{\link{evalALERT}},
\code{\link{robustALERT}}
}
\keyword{applyALERT}


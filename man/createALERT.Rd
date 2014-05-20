\name{createALERT}
\alias{createALERT}
\title{Producing the ALERT thresholds table}
\usage{
createALERT(data, firstMonth = 9, lag = 7, minWeeks = 8,
  allThresholds = FALSE, k = 0, target.pct = NULL, caseColumn = "Cases",
  lastDate = NULL)
}
\arguments{
  \item{data}{the historical data to use in the analysis. A
  data frame with a 'Date' column (must be \code{Date}
  objects) and a 'Cases' column.}

  \item{firstMonth}{month number which is counted as the
  first month of the 'flu year'}

  \item{lag}{lag time in days between date of cases and
  action taken}

  \item{minWeeks}{minimum number of weeks to be in ALERT}

  \item{allThresholds}{if \code{TRUE}, all integer
  threshold values between the 10th and 50th percentile are
  examined. If \code{FALSE}, only the 10th, 20th, 30th,
  40th, and 50th percentiles are examined.}

  \item{k}{the number of weeks around the peak to evaluate
  ALERT coverage for}

  \item{target.pct}{the percentage of cases the user is
  targeting during the ALERT period (optional)}

  \item{caseColumn}{the name of the column with the case
  counts in it. Defaults to 'Cases'}

  \item{lastDate}{a string in unambigous date format. All
  cases after this date will be removed prior to running
  the createALERT algorithm}
}
\value{
By utilizing the \code{\link{applyALERT}} function,
\code{createALERT} uses prior hospital data to
prospectively determine the start and end to a period of
elevated influenza incidence in a community.

\code{createALERT()$details} creates a list of matrices.
Each matrix contains raw statistics for the performance of
a threshold for each flu season in \code{data} (these
statistics can be found in \code{\link{applyALERT})}.

\code{createALERT()$out} produces a matrix summarizing the
performance of different ALERT thresholds. The columns for
this matrix are as follows: \item{threshold }{the minimum
threshold number of cases needed to begin the ALERT period}
\item{median.dur }{the median ALERT period duration in
weeks} \item{median.pct.cases.captured }{across all
seasons, the median percentage of all influenza cases
contained within the ALERT period}
\item{min.pct.cases.captured }{the minimum percentage of
annual cases captured during the ALERT period in any
season} \item{max.pct.cases.captured }{the maximum
percentage of annual cases captured during the ALERT period
in any season} \item{pct.peaks.captured }{the percentage of
seasons in which the ALERT period contained the peak week}
\item{pct.ext.peaks.captured }{the percentage of seasons in
which the ALERT period contained the peak week +/- \code{k}
weeks} \item{median.low.weeks.incl }{the median number of
weeks included in the ALERT period with counts less than
\code{threshold}} \item{median.duration.diff }{if
\code{target.pct} specified, the median difference between
the duration of the ALERT period in a season and the
duration of the shortest period needed to capture
\code{target.pct} of cases for that season}
}
\description{
The \code{createALERT} function calculates the ALERT
thresholds table given data and other parameters.
}
\note{
%% ~~further notes~~
}
\examples{
## Find the ALERT thresholds table for fluData over all thresholds
data(fluData)
x <- createALERT(data=fluData, allThresholds=TRUE, k=2, target.pct=0.85)
x$out
}
\author{
Nicholas G Reich and Stephen A Lauer
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
\code{\link{evalALERT}} cross-validates ALERT data over a
rule

\code{\link{robustALERT}} cross-validates ALERT data over a
series of rules
}
\keyword{createALERT}


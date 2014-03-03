\docType{package}
\name{ALERT-package}
\alias{ALERT}
\alias{ALERT-package}
\title{The Above Local Elevated Respiratory illness Threshold (ALERT) Algorithm}
\description{
A program that helps hospitals prospectively determine the
start and end to a period of elevated influenza incidence
in a community.
}
\details{
\tabular{ll}{ Package: \tab ALERT\cr Type: \tab Package\cr
Version: \tab 1.0\cr Date: \tab 2014-01-27\cr License: \tab
GPL-3\cr } Use createALERT to get a simple projection for
the upcoming flu season. Use robustALERT to cross-validate
createALERT under a set of rules.
}
\examples{
data(maskData)
## project the proper start date for the upcoming flu season
createALERT(data=maskData, firstMonth=8, lag=7, minWeeks=8, allThresholds=TRUE, k=2, target.pct=0.85)

## cross validate the above projection
robustALERT(minPercent=c(.8, .85, .9), maxDuration=c(12, 13, 14), data=maskData, firstMonth=8, lag=7, minWeeks=8, allThresholds=TRUE, k=2, target.pct=0.85)
}
\author{
Nicholas G Reich, Stephen A Lauer, et al.
}
\references{
TBD
}
\seealso{
createALERT
}
\keyword{ALERT}
\keyword{applyALERT}
\keyword{createALERT}
\keyword{evalALERT}
\keyword{postcastALERT}
\keyword{robustALERT}


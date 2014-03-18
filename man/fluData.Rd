\docType{data}
\name{fluData}
\alias{fluData}
\title{Masked ALERT test data}
\format{A weekly time series.}
\source{
%% ~~ reference to a publication or URL from which the data
were obtained ~~
}
\description{
Weekly influenza cases from a hospital
}
\examples{
data(fluData)

## project the proper start date for the upcoming flu season
createALERT(data=fluData, firstMonth=8, lag=7, minWeeks=8, allThresholds=TRUE, k=2, target.pct=0.85)

## cross validate the above projection
robustALERT(minPercent=c(.8, .85, .9), maxDuration=c(12, 13, 14), data=fluData, firstMonth=8, lag=7, minWeeks=8, allThresholds=TRUE, k=2, target.pct=0.85)
}
\references{
%% ~~ possibly secondary sources and usages ~~
}
\keyword{fluData}


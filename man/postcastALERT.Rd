\name{postcastALERT}
\alias{postcastALERT}
\title{A component of applyALERT}
\usage{
postcastALERT(data, target.pct)
}
\arguments{
  \item{data}{a single season of surveillance data}

  \item{target.pct}{specifies the percentage of cases the
  user is targeting during the ALERT period}
}
\value{
A vector with two elements, the \code{pct.captured} and
\code{duration} of the interval
}
\description{
The \code{postcastALERT} function is meant to be called
from within the \code{applyALERT} function and is used to
find the optimal window of time that contains a given
percentage of cases.
}
\note{
%% ~~further notes~~
}
\examples{
data(fluData)
postcastALERT(fluData, target.pct=.85)
}
\author{
Nicholas G Reich and Stephen A Lauer
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
\code{\link{applyALERT}}, \code{\link{createALERT}}
}
\keyword{postcastALERT}


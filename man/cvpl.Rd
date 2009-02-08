\name{cvpl}
\alias{cvpl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute the CVPL }
\description{
  The function computes the cross-validated partial likelihood (CVPL) for the Cox model.
}
\usage{
cvpl(x, surv.time, surv.event, strata.cox = NULL, nfold = 1, setseed = NULL, na.rm = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data matrix }
  \item{surv.time}{ vector of times to event occurrence }
  \item{surv.event}{ vector of indicators for event occurrence }
  \item{strata.cox}{ stratification variable }
  \item{nfold}{ number of folds for the cross-validation }
  \item{setseed}{ seed for the random generator }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
  \item{verbose}{ verbosity of the function }
}
\details{
  The function is not fully implemented yet
}
\value{
  -1
}
\references{ ~put references to the literature/web site here ~ }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link{logpl}}, \code{\link[survival]{coxph}} }
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }

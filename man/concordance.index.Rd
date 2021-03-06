\name{concordance.index}
\alias{concordance.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute the concordance index for survival or binary class prediction}
\description{
Function to compute the concordance index for a risk prediction, i.e. the probability that, for a pair of randomly chosen comparable samples, the sample with the higher risk prediction will experience an event before the other sample or belongs to a higher binary class.
}
\usage{
concordance.index(x, surv.time, surv.event, cl, weights, strat, alpha = 0.05, outx = TRUE, method = c("conservative", "noether", "nam"), na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of risk predictions.}
  \item{surv.time}{a vector of event times.}
  \item{surv.event}{a vector of event occurence indicators.}
  \item{cl}{a vector of binary class indicators.}
  \item{weights}{weight of each sample.}
  \item{strat}{stratification indicator.}
  \item{alpha}{apha level to compute confidence interval.}
  \item{outx}{set to \code{TRUE} to not count pairs of observations tied on \code{x} as a relevant pair. This results in a Goodman-Kruskal gamma type rank correlation.}
  \item{method}{can take the value \code{conservative}, \code{noether} or \code{name} (see paper Pencina et al. for details).}
  \item{na.rm}{\code{TRUE} if missing values should be removed.}
}
\details{
  Method \code{name} is not implemented yet.
}
\value{
  \item{c.index }{concordance index estimate.}
  \item{se }{standard error of the estimate.}
  \item{lower }{lower bound for the confidence interval.}
  \item{upper }{upper bound for the confidence interval.}
  \item{p.value }{p-value for the statistical test if the estimate if different from 0.5.}
  \item{n }{number of samples used for the estimation.}
  \item{data }{list of data used to compute the index (\code{x}, \code{surv.time} and \code{surv.event}, or \code{cl}).}
}
\references{Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) "Tutorial in biostatistics: multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing error", \emph{Statistics in Medicine}, \bold{15}, pages 361--387.

Pencina, M. J. and D'Agostino, R. B. (2004) "Overall C as a measure of discrimination in survival analysis: model specific population value and confidence interval estimation", \emph{Statistics in Medicine}, \bold{23}, pages 2109--2123, 2004.}
\author{Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link[Hmisc]{rcorr.cens}}, \code{\link[CPE]{phcpe}}, \code{\link[clinfun]{coxphCPE}}}
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
sex <- sample(0:1, 100, replace=TRUE)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
cat("survival prediction:\n")
concordance.index(x=age, surv.time=stime, surv.event=sevent, strat=strat, method="noether")
cat("binary class prediction:\n")
concordance.index(x=age, cl=sex, strat=strat, method="noether")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line

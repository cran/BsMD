\name{summary.BsProb}
\alias{summary.BsProb}
\title{Summary of Posterior Probabilities from Bayesian Screening}
\description{
    Reduced printing method for class \code{BsProb} lists. Prints
    posterior probabilities of factors and models from Bayesian screening
    procedure.}
\usage{
    summary.BsProb(object, nMod = 10, digits = 3, ...)
}
\arguments{
  \item{object}{list. \code{BsProb} class list. Output list of \code{\link{BsProb}} function.}
  \item{nMod}{integer. Number of the top ranked models to print.}
  \item{digits}{integer. Significant digits to use.}
  \item{...}{additional arguments passed to \code{summary} generic function.}
} \value{
    The function prints out the marginal factors and models posterior
    probabilities. Returns invisible list with the components:
    \item{calc}{Numeric vector with basic calculation information.}
    \item{probabilities}{Data frame with the marginal posterior factor probabilities.}
    \item{models}{Data frame with the models posterior probabilities.}
}

\references{ Box, G. E. P and R. D. Meyer (1986). "An Analysis for
Unreplicated Fractional Factorials". \emph{Technometrics.} Vol. 28. No. 1.
pp. 11--18.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors
in Fractionated Screening Experiments".
\emph{Journal of Quality Technology.} Vol. 25. No. 2. pp. 94--105.

}
\author{Ernesto Barrios.}
\seealso{
    \code{\link{BsProb}}, \code{\link{print.BsProb}}, \code{\link{plot.BsProb}}.
}
\examples{
library(BsMD)
data(BM86.data,package="BsMD")
X <- as.matrix(BM86.data[,1:15])
y <- BM86.data["y1"]
# Using prior probability of p = 0.20, and k = 10 (gamma = 2.49)
drillAdvance.BsProb <- BsProb(X = X, y = y, blk = 0, mFac = 15, mInt = 1,
            p = 0.20, g = 2.49, ng = 1, nMod = 10)
plot(drillAdvance.BsProb)
summary(drillAdvance.BsProb)

# Using prior probability of p = 0.20, and a 5 <= k <= 15 (1.22 <= gamma <= 3.74)
drillAdvance.BsProbG <- BsProb(X = X, y = y, blk = 0, mFac = 15, mInt = 1,
            p = 0.25, g = c(1.22, 3.74), ng = 3, nMod = 10)
plot(drillAdvance.BsProbG)
summary(drillAdvance.BsProbG)
}
\keyword{ design }

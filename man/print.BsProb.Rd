\name{print.BsProb}
\alias{print.BsProb}
\title{Printing Posterior Probabilities from Bayesian Screening}
\description{
    Printing method for  lists of class \code{BsProb}. Prints the posterior
    probabilities of factors and models from the Bayesian screening procedure.
}
\usage{
    print.BsProb(x, X = TRUE, resp = TRUE, factors = TRUE, models = TRUE,
            nMod = 10, digits = 3, plt = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{x}{list. Object of \code{BsProb} class, output from the \code{\link{BsProb}} function.}
  \item{X}{logical. If \code{TRUE}, the design matrix is printed.}
  \item{resp}{logical. If \code{TRUE}, the response vector is printed.}
  \item{factors}{logical. Marginal posterior probabilities are printed if \code{TRUE}.}
  \item{models}{logical. If \code{TRUE} models posterior probabilities are printed.}
  \item{nMod}{integer. Number of the top ranked models to print.}
  \item{digits}{integer. Significant digits to use for printing.}
  \item{plt}{logical. Factor marginal probabilities are plotted if \code{TRUE}.}
  \item{verbose}{logical. If \code{TRUE}, the \code{unclass}-ed list \code{x} is displayed.}
  \item{...}{additional arguments passed to \code{print} function.}
}
\value{
    The function prints out marginal factors and models posterior probabilities.
    Returns invisible list with the components:
    \item{calc}{numeric vector with general calculation information.}
    \item{probabilities}{Data frame with the marginal posterior factor probabilities.}
    \item{models}{Data frame with model the posterior probabilities.}
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
    \code{\link{BsProb}}, \code{\link{summary.BsProb}}, \code{\link{plot.BsProb}}.
}
\examples{
library(BsMD)
data(BM86.data,package="BsMD")
X <- as.matrix(BM86.data[,1:15])
y <- BM86.data["y1"]
# Using prior probability of p = 0.20, and k = 10 (gamma = 2.49)
drillAdvance.BsProb <- BsProb(X = X, y = y, blk = 0, mFac = 15, mInt = 1,
            p = 0.20, g = 2.49, ng = 1, nMod = 10)
print(drillAdvance.BsProb)
plot(drillAdvance.BsProb)

# Using prior probability of p = 0.20, and a 5 <= k <= 15 (1.22 <= gamma <= 3.74)
drillAdvance.BsProbG <- BsProb(X = X, y = y, blk = 0, mFac = 15, mInt = 1,
            p = 0.25, g = c(1.22, 3.74), ng = 3, nMod = 10)
print(drillAdvance.BsProbG, X = FALSE, resp = FALSE)
plot(drillAdvance.BsProbG)
}
\keyword{ design }

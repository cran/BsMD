\name{plot.BsProb}
\alias{plot.BsProb}
\title{Plotting of Posterior Probabilities from Bayesian Screening}
\description{
Method function for plotting marginal factor posterior probabilities for Bayesian screening.
}
\usage{
    \method{plot}{BsProb}(x, code = TRUE, prt = FALSE, cex.axis=par("cex.axis"), ...)
}
\arguments{
  \item{x}{list. List of class \code{BsProb} output from the
        \code{\link{BsProb}} function.}
  \item{code}{logical. If \code{TRUE} coded factor names are used.}
  \item{prt}{logical. If \code{TRUE}, summary of the posterior probabilities
        calculation is printed.}
  \item{cex.axis}{Magnification used for the axis annotation.
        See  \code{\link{par}}.}
  \item{\dots}{additional graphical parameters passed to \code{plot}.}
}
\details{
    A spike plot, similar to barplots, is produced with a spike for each factor.
    Marginal posterior probabilities are used for the vertical axis.
    If \code{code=TRUE}, \code{X1}, \code{X2}, ... are used to label the factors
    otherwise the original factor names are used.
    If \code{prt=TRUE}, the \code{\link{print.BsProb}} function is called
    and the posterior probabilities are displayed.
    When \code{\link{BsProb}} is called for more than one value of gamma (\code{g}),
    the spikes for each factor probability are overlapped to show the
    resulting range of each marginal probability.
}
\value{
    The function is called for its side effects. It returns an invisible
    \code{NULL}.
}

\references{
Box, G. E. P and R. D. Meyer (1986).
"An Analysis for Unreplicated Fractional Factorials".
\emph{Technometrics.} Vol. 28. No. 1. pp. 11--18.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors
in Fractionated Screening Experiments".
\emph{Journal of Quality Technology.} Vol. 25. No. 2. pp. 94--105.

}
\author{Ernesto Barrios.}
\seealso{
    \code{\link{BsProb}}, \code{\link{print.BsProb}}, \code{\link{summary.BsProb}}.
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
plot(drillAdvance.BsProbG, code = FALSE, prt = TRUE)
}
\keyword{ design }
\keyword{ hplot }

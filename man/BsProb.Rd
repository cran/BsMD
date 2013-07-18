\name{BsProb}
\alias{BsProb}
\title{Posterior Probabilities from Bayesian Screening Experiments}
\description{
Marginal factor posterior probabilities and model posterior probabilities
from designed screening experiments are calculated according to Box and Meyer's
Bayesian procedure.
}
\usage{
    BsProb(X, y, blk, mFac, mInt = 2, p = 0.25, g = 2, ng = 1, nMod = 10)
}
\arguments{
  \item{X}{Matrix. The design matrix.}
  \item{y}{vector. The response vector.}
  \item{blk}{integer. Number of blocking factors (>=0). These factors are
      accommodated in the first columns of matrix \code{X}. There are
      \code{ncol(X)-blk} design factors.}
  \item{mFac}{integer. Maximum number of factors included in the models.}
  \item{mInt}{integer <= 3. Maximum order of interactions considered in the models.}
  \item{p}{numeric. Prior probability assigned to active factors.}
  \item{g}{vector. Variance inflation factor(s) \eqn{\gamma}associated to active and interaction factors.}
  \item{ng}{integer <=20. Number of different variance inflation factors (\code{g}) used in calculations.}
  \item{nMod}{integer <=100. Number of models to keep with the highest posterior probability.}
}
\details{
Factor and model posterior probabilities are computed by Box and Meyer's Bayesian
procedure. The design factors are accommodated in the matrix \code{X} after
\code{blk} columns of the blocking factors. So, \code{ncol(X)-blk} design factors
are considered. If \code{g}, the variance inflation factor (VIF) \eqn{\gamma},
is a vector of length 1, the same VIF is used for factor main effects and interactions.
If the length of \code{g} is 2 and \code{ng} is 1, \code{g[1]} is used for factor
main effects and \code{g[2]} for the interaction effects. If \code{ng} greater than 1,
then \code{ng} values of VIFs between \code{g[1]} and \code{g[2]} are used for
calculations with the same \eqn{gamma} value for main effects and interactions.
The function calls the \acronym{FORTRAN} subroutine \file{bm} and captures summary results.
The output is a list of class \code{BsProb} for which
\code{print}, \code{plot} and \code{summary} methods are available.
}
\value{
A list with all output parameters of the \acronym{FORTRAN} subroutine \file{bm}.
The names of the list components are such that they match the original \acronym{FORTRAN}
code. Small letters used for capturing program's output.
  \item{X}{matrix. The design matrix.}
  \item{Y}{vector. The response vector.}
  \item{N}{integer. The number of runs.}
  \item{COLS}{integer. The number of design factors.}
  \item{BLKS}{integer. The number of blocking factors accommodated in the first
                columns of matrix \code{X}.}
  \item{MXFAC}{integer. Maximum number of factors considered in the models.}
  \item{MXINT}{integer. Maximum interaction order considered in the models.}
  \item{PI}{numeric. Prior probability assigned to the active factors.}
  \item{INDGAM}{integer. If 0, the same variance inflation factor (\code{GAMMA})
                is used for main and interactions effects.
                If \code{INDGAM ==1}, then \code{NGAM} different values
                of \code{GAMMA} were used.}
  \item{INDG2}{integer. If 1, the variance inflation factor \code{GAM2} was used for
                the interaction effects.}
  \item{NGAM}{integer. Number of different VIFs used for computations.}
  \item{GAMMA}{vector. Vector of variance inflation factors of length 1 or 2.}
  \item{NTOP}{integer. Number of models with the highest posterior probability}.
  \item{mdcnt}{integer. Total number of models evaluated.}
  \item{ptop}{vector. Vector of probabilities of the top \code{ntop} models.}
  \item{sigtop}{vector.  Vector of sigma-squared of the top \code{ntop} models.}
  \item{nftop}{integer.  Number of factors in each of the \code{ntop} models.}
  \item{jtop}{matrix. Matrix of the number of factors and their labels
                of the top \code{ntop} models.}
  \item{del}{numeric. Interval width of the \code{GAMMA} partition.}
  \item{sprob}{vector. Vector of posterior probabilities. If \code{ng>1} the
                probabilities are weighted averaged over \code{GAMMA}.}
  \item{pgam}{vector. Vector of values of the unscaled posterior density of \code{GAMMA}.}
  \item{prob}{matrix. Matrix of marginal factor posterior probabilities for each of
                the different values of \code{GAMMA}.}
  \item{ind}{integer. Indicator variable. \code{ind} is 1 if the \file{bm}
                subroutine exited properly. Any other number correspond to
                the format label number in the \acronym{FORTRAN} subroutine script.}

}
\references{
Box, G. E. P and R. D. Meyer (1986).
"An Analysis for Unreplicated Fractional Factorials".
\emph{Technometrics.} Vol. 28. No. 1. pp. 11--18.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors
in Fractionated Screening Experiments".
\emph{Journal of Quality Technology.} Vol. 25. No. 2. pp. 94--105.

}
\author{R. Daniel Meyer. Adapted for \R by Ernesto Barrios.}
\note{
The function is a wrapper to call the \acronym{FORTRAN} subroutine \file{bm},
modification of Daniel Meyer's original program, \file{mbcqp5.f}, for the
application of Bayesian design and analysis of fractional factorial experiments,
part of the \pkg{mdopt} bundle, available at \emph{StatLib}.
}
\seealso{
    \code{\link{print.BsProb}}, \code{\link{print.BsProb}}, \code{\link{summary.BsProb}}.
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

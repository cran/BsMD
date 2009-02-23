\name{summary.MD}
\alias{summary.MD}
\title{Summary of Best MD Follow-Up Experiments}
\description{
    Reduced printing method for lists of class \code{MD}. Displays the
    best MD criterion set of runs and their MD for follow-up experiments.
}
\usage{
summary.MD(object, digits = 3, verbose=FALSE, ...)
}
\arguments{
  \item{object}{list of \code{MD} class. Output list of \code{\link{MD}} function.}
  \item{digits}{integer. Significant digits to use in the print out.}
  \item{verbose}{logical. If \code{TRUE}, the \code{unclass}-ed \code{object} is displayed.}
  \item{...}{additional arguments passed to \code{summary} generic function.}
}
\value{
    It prints out the marginal factors and models posterior probabilities and
    the top MD follow-up experiments with their corresponding MD statistic.}
\references{
Meyer, R. D., Steinberg, D. M. and Box, G. E. P. (1996). "Follow-Up Designs
to Resolve Confounding in Multifactor Experiments (with discussion)".
\emph{Technometrics}, Vol. 38, No. 4, pp. 303--332.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors
in Fractionated Screening Experiments".
\emph{Journal of Quality Technology.} Vol. 25. No. 2. pp. 94--105.
}
\author{Ernesto Barrios.}
\seealso{\code{\link{print.MD}} and \code{\link{MD}}}
\examples{
### Reactor Experiment. Meyer et al. 1996, example 3.
library(BsMD)
data(Reactor.data,package="BsMD")

# Posterior probabilities based on first 8 runs
X <- as.matrix(cbind(blk = rep(-1,8), Reactor.data[c(25,2,19,12,13,22,7,32), 1:5]))
y <- Reactor.data[c(25,2,19,12,13,22,7,32), 6]
reactor.BsProb <- BsProb(X = X, y = y, blk = 1, mFac = 5, mInt = 3,
        p =0.25, g =0.40, ng = 1, nMod = 32)

# MD optimal 4-run design
p <- reactor.BsProb$ptop
s2 <- reactor.BsProb$sigtop
nf <- reactor.BsProb$nftop
facs <- reactor.BsProb$jtop
nFDes <- 4
Xcand <- as.matrix(cbind(blk = rep(+1,32), Reactor.data[,1:5]))
reactor.MD <- MD(X = X, y = y, nFac = 5, nBlk = 1, mInt = 3, g =0.40, nMod = 32,
        p = p,s2 = s2, nf = nf, facs = facs, nFDes = 4, Xcand = Xcand,
        mIter = 20, nStart = 25, top = 5)
print(reactor.MD)
summary(reactor.MD)
}
\keyword{ design }

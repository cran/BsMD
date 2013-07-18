\name{print.MD}
\alias{print.MD}
\title{Print Best MD Follow-Up Experiments}
\description{
    Printing method for  lists of class \code{MD}. Displays the
    best MD criterion set of runs and their MD for follow-up experiments.
}
\usage{
    \method{print}{MD}(x, X = FALSE, resp = FALSE, Xcand = TRUE, models = TRUE, nMod = x$nMod,
            digits = 3, verbose=FALSE, ...)
}
\arguments{
  \item{x}{list of class \code{MD}. Output list of the \code{\link{MD}} function.}
  \item{X}{logical. If \code{TRUE}, the initial design matrix is printed.}
  \item{resp}{logical If \code{TRUE}, the response vector of initial design is printed.}
  \item{Xcand}{logical. Prints the candidate runs if \code{TRUE}.}
  \item{models}{logical. Competing models are printed if \code{TRUE}.}
  \item{nMod}{integer. Top models to print.}
  \item{digits}{integer. Significant digits to use in the print out.}
  \item{verbose}{logical. If \code{TRUE}, the \code{unclass}-ed \code{x} is displayed.}
  \item{...}{additional arguments passed to \code{print} generic function.}
}
\value{
    The function is mainly called for its side effects. Prints out the selected
    components of the class \code{MD} objects, output of the \code{MD} function.
    For example the marginal factors and models posterior probabilities and
    the top MD follow-up experiments with their corresponding MD statistic.
    It returns invisible list with the components:
    \item{calc}{Numeric vector with basic calculation information.}
    \item{models}{Data frame with the competing models posterior probabilities.}
    \item{follow-up}{Data frame with the runs for follow-up experiments and
        their corresponding MD statistic.}
}
\references{
Meyer, R. D., Steinberg, D. M. and Box, G. E. P. (1996). "Follow-Up Designs
to Resolve Confounding in Multifactor Experiments (with discussion)".
\emph{Technometrics}, Vol. 38, No. 4, pp. 303--332.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors
in Fractionated Screening Experiments".
\emph{Journal of Quality Technology.} Vol. 25. No. 2. pp. 94--105.

}
\author{Ernesto Barrios.}
\seealso{\code{\link{MD}}, \code{\link{BsProb}}}
\examples{
# Injection Molding Experiment. Meyer et al. 1996. Example 2.
# MD for one extra experiment.
library(BsMD)
data(BM93.e3.data,package="BsMD")
X <- as.matrix(BM93.e3.data[1:16,c(1,2,4,6,9)])
y <- BM93.e3.data[1:16,10]
nBlk <- 1
nFac <- 4
mInt <- 3
g <- 2
nMod <- 5
p <- c(0.2356,0.2356,0.2356,0.2356,0.0566)
s2 <- c(0.5815,0.5815,0.5815,0.5815,0.4412)
nf <- c(3,3,3,3,4)
facs <- matrix(c(2,1,1,1,1,3,3,2,2,2,4,4,3,4,3,0,0,0,0,4),nrow=5,
    dimnames=list(1:5,c("f1","f2","f3","f4")))
nFDes <- 1
Xcand <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,
                    -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,
                    -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,
                    -1,1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1),
                    nrow=16,dimnames=list(1:16,c("blk","f1","f2","f3","f4"))
                )
mIter <- 0
startDes <- matrix(c(9,11,12,15),nrow=4)
top <- 10
injectionMolding.MD <- MD(X=X,y=y,nFac=nFac,nBlk=nBlk,mInt=mInt,g=g,
            nMod=nMod,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=nFDes,Xcand=Xcand,mIter=mIter,startDes=startDes,top=top)

print(injectionMolding.MD)
summary(injectionMolding.MD)

}
\keyword{ design }

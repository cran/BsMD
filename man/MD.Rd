\name{MD}
\alias{MD}
\title{Best Model Discrimination (MD) Follow-Up Experiments}
\description{
Best follow-up experiments based on the MD criterion
are suggested to discriminate between competing models.
}
\usage{
MD(X, y, nFac, nBlk = 0, mInt = 3, g = 2,  nMod, p, s2, nf, facs, nFDes = 4,
Xcand, mIter = 20, nStart = 5, startDes = NULL, top = 20, eps = 1e-05)
}
\arguments{
  \item{X}{matrix. Design matrix of the initial experiment.}
  \item{y}{vector. Response vector of the initial experiment.}
  \item{nFac}{integer. Number of factors in the initial experiment.}
  \item{nBlk}{integer >=1. The number of blocking factors in the initial experiment.
                They are accommodated in the first columns of matrix \code{X}.}
  \item{mInt}{integer. Maximum order of the interactions in the models.}
  \item{g}{vector. Variance inflation factor for main effects (\code{g[1]})
        and interactions effects (\code{g[2]}). If vector length
        is 1 the same inflation factor is used for main and
        interactions effects.}
  \item{nMod}{integer. Number of competing models.}
  \item{p}{vector. Posterior probabilities of the competing models.}
  \item{s2}{vector. Competing model variances.}
  \item{nf}{vector. Factors considered in each of the models.}
  \item{facs}{matrix. Matrix [\code{nMod x max(nf)}] of factor numbers in the design matrix.}
  \item{nFDes}{integer. Number of runs to consider in the follow-up experiment.}
  \item{Xcand}{matrix. Candidate runs to be chosen for the follow-up design.}
  \item{mIter}{integer. If 0, then user-entered designs \code{startDes} are
        evaluated, otherwise the maximum number of iterations for each Wynn search.}
  \item{nStart}{integer. Number of starting designs.}
  \item{startDes}{matrix. Matrix \code{[nStart x nFDes]}. Each row has the row numbers
        of the user-supplied starting design.}
  \item{top}{integer. Highest MD follow-up designs recorded.}
  \item{eps}{numeric. A small number (1e-5 by default) used for computations.}
}
\details{
The MD criterion, proposed by Meyer, Steinberg and Box is used to discriminate
among competing models. Random starting runs chosen from \code{Xcand} are used
for the Wynn search of best MD follow-up designs. \code{nStart} starting points are
tried in the search limited to \code{mIter} iterations. If \code{mIter=0} then
\code{startDes} user-provided designs are used. Posterior probabilities and
variances of the competing models are obtained from \code{\link{BsProb}}.
The function calls the \acronym{FORTRAN} subroutine \file{md} and captures
summary results. The complete output of the \acronym{FORTRAN} code is save in
the \file{MDPrint.out} file in the working directory.
}
\value{
A list with all input and output parameters of the \acronym{FORTRAN}
subroutine \code{MD}. Most of the variable names kept to match \acronym{FORTRAN} code.
    \item{NSTART}{Number of starting designs.}
    \item{NRUNS}{Number of runs used in follow-up designs.}
    \item{ITMAX}{Maximum number of iterations for each Wynn search.}
    \item{INITDES}{Number of starting points.}
    \item{NO}{Numbers of runs already completed before follow-up.}
    \item{IND}{Indicator; 0 indicates the user supplied starting designs.}
    \item{X}{Matrix for initial data (\code{nrow(X)=N0}; \code{ncol(X)=COLS+BL}).}
    \item{Y}{Response values from initial experiment (\code{length(Y)=N0}).}
    \item{GAMMA}{Variance inflation factor.}
    \item{GAM2}{If \code{IND=1}, \code{GAM2} was used for interaction factors.}
    \item{BL}{Number of blocks (>=1) accommodated in first columns of \code{X} and \code{Xcand}}.
    \item{COLS}{Number of factors.}
    \item{N}{Number of candidate runs.}
    \item{Xcand}{Matrix of candidate runs. (\code{nrow(Xcand)=N}, \code{ncol(Xcand)=ncol(X)}).}
    \item{NM}{Number of models considered.}
    \item{P}{Models posterior probability.}
    \item{SIGMA2}{Models variances.}
    \item{NF}{Number of factors per model.}
    \item{MNF}{Maximum number of factor in models. (\code{MNF=max(NF)}).}
    \item{JFAC}{Matrix with the factor numbers for each of the models.}
    \item{CUT}{Maximum interaction order considered.}
    \item{MBEST}{If \code{INITDES=0}, the first row of the \code{MBEST[1,]} matrix
        has the first user-supplied starting design. The last row the
        \code{NSTART}-th user-supplied starting design.}
    \item{NTOP}{Number of the top best designs.}
    \item{TOPD}{The D value for the best \code{NTOP} designs.}
    \item{TOPDES}{Top \code{NTOP} design factors.}
    \item{ESP}{"Small number" provided to the \file{md} \acronym{FORTRAN}
        subroutine. 1e-5 by default.}
    \item{flag}{Indicator = 1, if the \file{md} subroutine finished properly,
        -1 otherwise.}
}

\note{
The function is a wrapper to call the \acronym{FORTAN} subroutine \file{md},
modification of Daniel Meyer's original program, \file{MD.f}, part of the
\pkg{mdopt} bundle for Bayesian model discrimination of multifactor
experiments.
}
\references{
Meyer, R. D., Steinberg, D. M. and Box, G. E. P. (1996). "Follow-Up Designs
to Resolve Confounding in Multifactor Experiments (with discussion)".
\emph{Technometrics}, Vol. 38, No. 4, pp. 303--332.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors
in Fractionated Screening Experiments".
\emph{Journal of Quality Technology.} Vol. 25. No. 2. pp. 94--105.
}

\author{R. Daniel Meyer. Adapted for \R by Ernesto Barrios.}
\seealso{
    \code{\link{print.MD}}, \code{\link{BsProb}}
}
\examples{
### Injection Molding Experiment. Meyer et al. 1996, example 2.
library(BsMD)
data(BM93.e3.data,package="BsMD")
X <- as.matrix(BM93.e3.data[1:16,c(1,2,4,6,9)])
y <- BM93.e3.data[1:16,10]
p <- c(0.2356,0.2356,0.2356,0.2356,0.0566)
s2 <- c(0.5815,0.5815,0.5815,0.5815,0.4412)
nf <- c(3,3,3,3,4)
facs <- matrix(c(2,1,1,1,1,3,3,2,2,2,4,4,3,4,3,0,0,0,0,4),nrow=5,
    dimnames=list(1:5,c("f1","f2","f3","f4")))
nFDes <- 4
Xcand <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,
                    -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,
                    -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,
                    -1,1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1),
                    nrow=16,dimnames=list(1:16,c("blk","f1","f2","f3","f4"))
                )
injectionMolding.MD <- MD(X = X, y = y, nFac = 4, nBlk = 1, mInt = 3,
            g = 2, nMod = 5, p = p, s2 = s2, nf = nf, facs = facs,
            nFDes = 4, Xcand = Xcand, mIter = 20, nStart = 25, top = 10)
summary(injectionMolding.MD)



### Reactor Experiment. Meyer et al. 1996, example 3.
par(mfrow=c(1,2),pty="s")
data(Reactor.data,package="BsMD")

# Posterior probabilities based on first 8 runs
X <- as.matrix(cbind(blk = rep(-1,8), Reactor.data[c(25,2,19,12,13,22,7,32), 1:5]))
y <- Reactor.data[c(25,2,19,12,13,22,7,32), 6]
reactor8.BsProb <- BsProb(X = X, y = y, blk = 1, mFac = 5, mInt = 3,
        p =0.25, g =0.40, ng = 1, nMod = 32)
plot(reactor8.BsProb,prt=TRUE,,main="(8 runs)")

# MD optimal 4-run design
p <- reactor8.BsProb$ptop
s2 <- reactor8.BsProb$sigtop
nf <- reactor8.BsProb$nftop
facs <- reactor8.BsProb$jtop
nFDes <- 4
Xcand <- as.matrix(cbind(blk = rep(+1,32), Reactor.data[,1:5]))
reactor.MD <- MD(X = X, y = y, nFac = 5, nBlk = 1, mInt = 3, g =0.40, nMod = 32,
        p = p,s2 = s2, nf = nf, facs = facs, nFDes = 4, Xcand = Xcand,
        mIter = 20, nStart = 25, top = 5)
summary(reactor.MD)

# Posterior probabilities based on all 12 runs
X <- rbind(X, Xcand[c(4,10,11,26), ])
y <- c(y, Reactor.data[c(4,10,11,26),6])
reactor12.BsProb <- BsProb(X = X, y = y, blk = 1, mFac = 5, mInt = 3,
        p = 0.25, g =1.20,ng = 1, nMod = 5)
plot(reactor12.BsProb,prt=TRUE,main="(12 runs)")
}
\keyword{ design }

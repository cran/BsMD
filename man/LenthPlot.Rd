\name{LenthPlot}
\alias{LenthPlot}
\title{Lenth's Plot of Effects}
\description{
Plot of the factor effects with significance levels based on robust
estimation of contrast standard errors.}
}
\usage{
LenthPlot(obj, alpha = 0.05, plt = TRUE, xlab = "factors", ylab = "effects",
    faclab = NULL, cex.fac = par("cex.lab"),
    cex.axis=par("cex.axis"), adj = 1, ...)
}
\arguments{
  \item{obj}{object of class \code{lm} or vector with the factor effects.}
  \item{alpha}{numeric. Significance level used for the \emph{margin of error} (ME)
            and \emph{simultaneous margin of error} (SME). See Lenth(1989).}
  \item{plt}{logical. If \code{TRUE}, a spikes plot with the factor effects is displayed.
            Otherwise, no plot is produced.}
  \item{xlab}{character string. Used to label the x-axis. "factors" as default.}
  \item{ylab}{character string. Used to label the y-axis. "effects" as default.}
  \item{faclab}{list with components \code{fac} (numeric vector)
        and \code{lab} (character vector). The \code{fac} entries of
        effects vector (taken from \code{obj}) are labelled as \code{lab}.
        The rest of the effect names are blanked. If \code{NULL} all factors
        are labelled using the coefficients' name.}
  \item{cex.fac}{numeric. Character size used for the factor labels.}
  \item{cex.axis}{numeric. Character size used for the axis.}
  \item{adj}{numeric between 0 and 1. Determines where to place the
    "ME" (margin of error) and the "SME" (simultaneous margin of error) labels
    (character size of 0.9*\code{cex.axis}).
    0 for extreme left hand side, 1 for extreme right hand side.}
  \item{\dots}{extra parameters passed to \code{plot}.}
}
\details{
    If \code{obj} is of class \code{lm}, \code{2*coef(obj)} is used as factor
    effect with the intercept term removed. Otherwise, \code{obj} should be a
    vector with the factor effects. Robust estimate of the contrasts standard
    error is used to calculate \emph{marginal} (ME) and \emph{simultaneous margin
    of error} (SME) for the provided significance (\code{1 - alpha}) level.
    See Lenth(1989).  Spikes are used to display the factor effects.
    If \code{faclab} is \code{NULL}, factors are labelled with the effects or
    coefficient names. Otherwise, those \code{faclab\$fac} factors are labelled
    as \code{faclab\$lab}. The rest of the factors are blanked.
}
\value{
    The function is called mainly for its side effect. It returns a vector with
    the value of alpha used, ME and SME.
}
\references{Lenth, R. V. (1989).
    "Quick and Easy Analysis of Unreplicated Factorials".
    \emph{Technometrics} Vol. 31, No. 4. pp. 469--473.
}
\author{Ernesto Barrios}
\seealso{\code{\link{DanielPlot}}, \code{\link{BsProb}} and \code{\link{plot.BsProb}}}
\examples{
### Tensile Strength Experiment. Taguchi and Wu. 1980
library(BsMD)
# Data
data(BM86.data,package="BsMD")     # Design matrix and responses
print(BM86.data)    # from Box and Meyer (1986)

# Model Fitting. Box and Meyer (1986) example 2.
tensileStrength.lm <- lm(y2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    X10 + X11 + X12 + X13 + X14 + X15, data = BM86.data)
print(coef(tensileStrength.lm)) # Model coefficients

par(mfrow=c(1,2),pty="s")
DanielPlot(tensileStrength.lm, main = "Daniel Plot")
LenthPlot(tensileStrength.lm, main = "Lenth's Plot")
}
\keyword{design}

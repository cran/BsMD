\name{DanielPlot}
\alias{DanielPlot}
\title{Normal Plot of Effects}
\description{
Normal plot of effects from a two level factorial experiment.
}
\usage{
DanielPlot(fit, code = FALSE, faclab = NULL, block = FALSE,
    datax = TRUE, half = FALSE, pch = "*", cex.lab = par("cex.lab"),
    cex.pch = par("cex.axis"), ...)
}
\arguments{
  \item{fit}{object of class \code{lm}. Fitted model from \code{lm} or \code{aov}.}
  \item{code}{logical. If \code{TRUE} labels "A","B", etc are used instead of the
        names of the coefficients (factors).}
  \item{faclab}{list. If \code{NULL} points are labelled accordingly to \code{code},
    otherwise \code{faclab} should be a list with \code{idx} (integer vector) and \code{lab}
    (character vector) components. See Details.}
  \item{block}{logical. If \code{TRUE}, the first factor is labelled as "BK" (block).}
  \item{datax}{logical. If \code{TRUE}, the x-axis is used for the factor effects the
                the y-axis for the normal scores. The opposite otherwise.}
  \item{half}{logical. If \code{TRUE}, half-normal plot of effects is display.}
  \item{pch}{numeric or character. Points character.}
  \item{cex.lab}{numeric. Labels character size.}
  \item{cex.pch}{numeric. Points character size.}
  \item{\dots}{graphical parameters passed to \code{plot}.}
}
\details{
    The two levels design are assumed -1 and 1. Factor effects assumed \code{2*coef(obj)}
    ((Intercept) removed) are displayed in a \code{qqnorm} plot with the effects in
    the x-axis by default. If \code{half=TRUE} the half-normal plots of effects is
    plotted as the normal quantiles of \code{0.5*(rank(abs(effects))-0.5)/length(effects)+1}
    versus \code{abs(effects)}.
}

\value{
    The function returns invisible data frame with columns: \code{x}, \code{y}
    and \code{no}, for the coordinates and the enumeration of plotted points.
    Names of the factor effects (coefficients) are the row names of the data frame.
}
\references{
    C. Daniel (1976).
    \emph{Application of Statistics to Industrial Experimentation.}
    Wiley.

Box G. E. P, Hunter, W. C. and Hunter, J. S. (1978).
    \emph{Statistics for Experimenters}.
    Wiley.
}
\author{Ernesto Barrios.}
\seealso{\code{\link{qqnorm}}, \code{\link{LenthPlot}}}
\examples{
### Injection Molding Experiment. Box et al. 1978.
library(BsMD)
# Data
data(BM86.data,package="BsMD")     # Design matrix and response
print(BM86.data)    # from Box and Meyer (1986)

# Model Fitting. Box and Meyer (1986) example 3.
injectionMolding.lm <- lm(y3 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    X10 + X11 + X12 + X13 + X14 + X15, data = BM86.data)
print(coef(injectionMolding.lm)) # Model coefficients

# Daniel Plots
par(mfrow=c(1,3),oma=c(0,0,1,0),pty="s")
DanielPlot(injectionMolding.lm, half = TRUE, main = "Half-Normal Plot")
DanielPlot(injectionMolding.lm, main = "Normal Plot of Effects")
DanielPlot(injectionMolding.lm,
        faclab = list(idx = c(12,4,13), lab = c(" -H"," VG"," -B")),
        main = "Active Contrasts")
}
\keyword{design}

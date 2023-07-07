plot.BsProb <-
function (x, code = TRUE, prt = FALSE, cex.axis = par("cex.axis"), 
    ...) 
{
    spikes <- function(prob, lwd = 3, ...) {
        y <- prob
        n <- nrow(y)
        x <- seq(n)
        lab <- rownames(prob)
        plot(x, y[, 1], xlim = range(x), ylim = c(0, 1), type = "n", 
            xlab = "factors", ylab = "posterior probability", 
            frame = FALSE, axes = FALSE, ...)
        if (ncol(y) == 1) {
            for (i in x) segments(x[i], 0, x[i], y[i, 1], lwd = lwd, 
                col = grey(0.20000000000000001))
        }
        else {
            y[, 1] <- apply(prob, 1, min)
            y[, 2] <- apply(prob, 1, max)
            for (i in x) {
                segments(x[i], 0, x[i], y[i, 2], lwd = lwd, col = grey(0.80000000000000004), 
                  lty = 1)
                segments(x[i], 0, x[i], y[i, 1], lwd = lwd, col = grey(0.20000000000000001), 
                  lty = 1)
            }
        }
        axis(1, at = x, labels = lab, line = 0, cex.axis = cex.axis)
        axis(2, cex.axis = cex.axis)
        invisible(NULL)
    }
    if (!any(class(x) == "BsProb")) 
        return("\nArgument `x' should be class BsProb. Output of corresponding function.")
    ifelse(x$INDGAM == 0, prob <- as.matrix(x$sprob), prob <- x$prob)
    if (code) 
        rownames(prob) <- rownames(x$prob)
    else rownames(prob) <- names(x$sprob)
    spikes(prob, ...)
    if (prt) 
        summary.BsProb(x)
    invisible(NULL)
}

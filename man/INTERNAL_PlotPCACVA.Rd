\name{INTERNAL_PlotPCACVA}
\alias{INTERNAL_PlotPCACVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL_PlotPCACVA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
INTERNAL_PlotPCACVA(variables, individuals, suppIndividualsToPlot = FALSE, suppIndividuals = NULL, biplot = TRUE, variablesColors = NULL, individualsColors = NULL, suppIndividualsColors = NULL, variablesLabels = TRUE, individualsLabels = TRUE, suppIndividualsLabels = FALSE, individualsEllipses = NULL, individualsSegments = NULL, xlab = "", ylab = "", expand = NULL, mainTitle = NULL, subTitle = NULL, returnX = FALSE, returnY = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variables}{
%%     ~~Describe \code{variables} here~~
}
  \item{individuals}{
%%     ~~Describe \code{individuals} here~~
}
  \item{suppIndividualsToPlot}{
%%     ~~Describe \code{suppIndividualsToPlot} here~~
}
  \item{suppIndividuals}{
%%     ~~Describe \code{suppIndividuals} here~~
}
  \item{biplot}{
%%     ~~Describe \code{biplot} here~~
}
  \item{variablesColors}{
%%     ~~Describe \code{variablesColors} here~~
}
  \item{individualsColors}{
%%     ~~Describe \code{individualsColors} here~~
}
  \item{suppIndividualsColors}{
%%     ~~Describe \code{suppIndividualsColors} here~~
}
  \item{variablesLabels}{
%%     ~~Describe \code{variablesLabels} here~~
}
  \item{individualsLabels}{
%%     ~~Describe \code{individualsLabels} here~~
}
  \item{suppIndividualsLabels}{
%%     ~~Describe \code{suppIndividualsLabels} here~~
}
  \item{individualsEllipses}{
%%     ~~Describe \code{individualsEllipses} here~~
}
  \item{individualsSegments}{
%%     ~~Describe \code{individualsSegments} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{expand}{
%%     ~~Describe \code{expand} here~~
}
  \item{mainTitle}{
%%     ~~Describe \code{mainTitle} here~~
}
  \item{subTitle}{
%%     ~~Describe \code{subTitle} here~~
}
  \item{returnX}{
%%     ~~Describe \code{returnX} here~~
}
  \item{returnY}{
%%     ~~Describe \code{returnY} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (variables, individuals, suppIndividualsToPlot = FALSE, 
    suppIndividuals = NULL, biplot = TRUE, variablesColors = NULL, 
    individualsColors = NULL, suppIndividualsColors = NULL, variablesLabels = TRUE, 
    individualsLabels = TRUE, suppIndividualsLabels = FALSE, 
    individualsEllipses = NULL, individualsSegments = NULL, xlab = "", 
    ylab = "", expand = NULL, mainTitle = NULL, subTitle = NULL, 
    returnX = FALSE, returnY = FALSE) 
{
    par(mar = c(4, 4, 1, 1))
    par(mgp = c(2.5, 1, 0))
    par(oma = c(2, 0, 1, 0))
    if (returnX) {
        individuals[, 1] = -individuals[, 1]
        variables[, 1] = -variables[, 1]
        if (!is.null(individualsEllipses)) {
            individualsEllipses[, , 1] = -individualsEllipses[, 
                , 1]
        }
    }
    if (returnY) {
        individuals[, 2] = -individuals[, 2]
        variables[, 2] = -variables[, 2]
        if (!is.null(individualsEllipses)) {
            individualsEllipses[, , 2] = -individualsEllipses[, 
                , 2]
        }
    }
    if (is.null(expand)) {
        if (biplot) {
            max.norm.prod = max(abs(individuals[, 1]))
            max.norm.suj = max(abs(variables[, 1]))
            expand.conseil = max.norm.prod/max.norm.suj
            expand = round(expand.conseil, digits = 2)
        }
        if (!biplot) {
            expand = 1
        }
    }
    labelx = rownames(individuals)
    labely = rownames(variables)
    nb.var = dim(variables)[1]
    nb.ind = dim(individuals)[1]
    if (is.null(variablesColors)) {
        variablesColors = rep("red", nb.var)
    }
    if (is.null(individualsColors)) {
        individualsColors = rep("blue", nb.ind)
    }
    if (is.null(suppIndividualsColors)) {
        suppIndividualsColors = rainbow(nb.ind)
    }
    xMinInd = min(individuals[, 1])
    xMaxInd = max(individuals[, 1])
    yMinInd = min(individuals[, 2])
    yMaxInd = max(individuals[, 2])
    xMinVar = min(variables[, 1] * expand)
    xMaxVar = max(variables[, 1] * expand)
    yMinVar = min(variables[, 2] * expand)
    yMaxVar = max(variables[, 2] * expand)
    if (!is.null(individualsEllipses)) {
        xMinInd = min(xMinInd, min(individualsEllipses[, , 1]))
        xMaxInd = max(xMaxInd, max(individualsEllipses[, , 1]))
        yMinInd = min(yMinInd, min(individualsEllipses[, , 2]))
        yMaxInd = max(yMaxInd, max(individualsEllipses[, , 2]))
    }
    if (suppIndividualsToPlot) {
        xMinInd = min(suppIndividuals[[1]][, 1])
        xMaxInd = max(suppIndividuals[[1]][, 1])
        yMinInd = min(suppIndividuals[[1]][, 2])
        yMaxInd = max(suppIndividuals[[1]][, 2])
        if (returnX) {
            suppIndividuals[[1]][, 1] = -suppIndividuals[[1]][, 
                1]
        }
        if (returnY) {
            suppIndividuals[[1]][, 2] = -suppIndividuals[[1]][, 
                2]
        }
        for (i in 2:nb.ind) {
            xMinInd = min(xMinInd, min(suppIndividuals[[i]][, 
                1]))
            xMaxInd = max(xMaxInd, max(suppIndividuals[[i]][, 
                1]))
            yMinInd = min(yMinInd, min(suppIndividuals[[i]][, 
                2]))
            yMaxInd = max(yMaxInd, max(suppIndividuals[[i]][, 
                2]))
            if (returnX) {
                suppIndividuals[[i]][, 1] = -suppIndividuals[[i]][, 
                  1]
            }
            if (returnY) {
                suppIndividuals[[i]][, 2] = -suppIndividuals[[i]][, 
                  2]
            }
        }
    }
    if (biplot == TRUE) {
        xMinC1 = min(xMinInd, xMinVar) * 1.1
        xMaxC1 = max(xMaxInd, xMaxVar) * 1.1
        yMinC1 = min(yMinInd, yMinVar) * 1.1
        yMaxC1 = max(yMaxInd, yMaxVar) * 1.1
        offsetY1 = (yMaxC1 - yMinC1)/30
    }
    else {
        xMinC1 = xMinInd * 1.1
        xMaxC1 = xMaxInd * 1.1
        yMinC1 = yMinInd * 1.1
        yMaxC1 = yMaxInd * 1.1
        xMinC2 = xMinVar * 1.1
        xMaxC2 = xMaxVar * 1.1
        yMinC2 = yMinVar * 1.1
        yMaxC2 = yMaxVar * 1.1
        offsetY1 = (yMaxC1 - yMinC1)/30
        offsetY2 = (yMaxC2 - yMinC2)/30
        par(mfcol = c(1, 2))
    }
    plot(NULL, pch = 3, xlim = c(xMinC1, xMaxC1), ylim = c(yMinC1, 
        yMaxC1), asp = 1, xlab = xlab, ylab = ylab)
    if (!is.null(subTitle)) {
        mtext(subTitle, outer = T, side = 1, line = 1, cex = 0.8)
    }
    abline(v = 0, lty = 3)
    abline(h = 0, lty = 3)
    points(individuals, col = individualsColors, pch = 3)
    if (individualsLabels == TRUE) {
        text(individuals[, 1], individuals[, 2] + offsetY1, labels = labelx, 
            col = individualsColors, cex = 0.8)
    }
    if (!is.null(individualsEllipses)) {
        for (i in 1:nb.ind) {
            lines(individualsEllipses[i, , ], col = individualsColors[i])
        }
    }
    if (suppIndividualsToPlot) {
        for (i in 1:nb.ind) {
            points(suppIndividuals[[i]], col = suppIndividualsColors[i], 
                pch = 16)
        }
    }
    if (!is.null(individualsSegments)) {
        for (i in 1:(nb.ind - 1)) {
            ind1 = rownames(individualsSegments)[i]
            for (j in (i + 1):(nb.ind)) {
                ind2 = rownames(individualsSegments)[j]
                if (individualsSegments[i, j] == TRUE) {
                  segments(individuals[rownames(individuals) == 
                    ind1, 1], individuals[rownames(individuals) == 
                    ind1, 2], individuals[rownames(individuals) == 
                    ind2, 1], individuals[rownames(individuals) == 
                    ind2, 2])
                }
            }
        }
    }
    if (biplot) {
        if (variablesLabels == TRUE) {
            text(variables[, 1] * expand, (variables[, 2] * expand) + 
                (sign(variables[, 2]) * offsetY1), labels = labely, 
                col = variablesColors, cex = 0.8)
        }
        for (i in 1:nb.var) {
            segments(variables[i, 1] * expand, variables[i, 2] * 
                expand, 0, 0, col = variablesColors[i])
        }
    }
    else {
        plot(NULL, pch = 3, xlim = c(-1, 1), ylim = c(-1, 1), 
            asp = 1, xlab = xlab, ylab = ylab, main = mainTitle)
        x.cercle <- seq(-1, 1, by = 0.01)
        y.cercle <- sqrt(1 - x.cercle^2)
        lines(x.cercle, y = y.cercle)
        lines(x.cercle, y = -y.cercle)
        for (i in 1:nb.var) {
            segments(variables[i, 1], variables[i, 2] + (sign(variables[, 
                2]) * offsetY2), 0, 0, col = variablesColors[i])
        }
        if (variablesLabels == TRUE) {
            text(variables[, 1], variables[, 2], labels = labely, 
                col = variablesColors, cex = 0.8)
        }
        abline(v = 0, lty = 3)
        abline(h = 0, lty = 3)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

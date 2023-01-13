#'INTERNAL_PlotPCACVA
#' @importFrom graphics box abline points text lines segments arrows mtext par
#'@param variables   matrix of variables
#'@param individuals matrix of individuals
#'@param suppIndividualsToPlot TRUE or FALSE
#'@param suppIndividuals matrix of suppIndividuals to use
#'@param biplot TRUE by default, but can be FALSE
#'@param variablesColors vector of colors for variables
#'@param individualsColors vector of colors for variables
#'@param suppIndividualsColors vector of colors for supplementary individuals
#'@param individualsLabels TRUE  if to be displayed, FALSE otherwise
#'@param suppIndividualsLabels TRUE  if to be displayed, FALSE otherwise
#'@param variablesLabels TRUE  if to be displayed, FALSE otherwise
#'@param individualsEllipses array containing ellipses coordinates
#'@param individualsSegments TRUE or FALSE
#'@param xlab label for x axis
#'@param ylab label for y axis
#'@param expand useful for biplot. Number indicating how to expand the graph
#'@param mainTitle title of the graph
#'@param subTitle subtitle of the graph
#'@param returnX FALSE by default. Boolean indicating if the x axis should be reverted
#'@param returnY FALSE by default. Boolean indicating if the y axis should be reverted
#'@param listEllipsesByRep internal parameters NULL by default
#'@param cex size of text
#'@param xlim vector of two numbers indicating the limits of x axis
#'@param ylim vector of two numbers indicating the limits of y axis
#'confInt   limit for the confidence ellipses. By default 0.9
#' @export
INTERNAL_PlotPCACVA <-function (variables, individuals, suppIndividualsToPlot = FALSE,
    suppIndividuals = NULL, biplot = TRUE, variablesColors = NULL,
    individualsColors = NULL, suppIndividualsColors = NULL, variablesLabels = TRUE,
    individualsLabels = TRUE, suppIndividualsLabels = FALSE,
    individualsEllipses = NULL, individualsSegments = NULL, xlab = "",
    ylab = "", expand = NULL, mainTitle = NULL, subTitle = NULL,
    returnX = FALSE, returnY = FALSE, listEllipsesByRep = NULL,
    cex = 0.8, xlim = NULL, ylim = NULL)
{
    par(mar = c(4, 4, 1, 1))
    par(mgp = c(1, 1, 0))
    par(oma = c(2, 1, 3, 0))
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
    alphabeticalOrder = order(individualsColors)
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
    if (!is.null(listEllipsesByRep)) {
        for (x in listEllipsesByRep) {
            xMinInd = min(xMinInd, min(x[[2]][, , 1]))
            xMaxInd = max(xMaxInd, max(x[[2]][, , 1]))
            yMinInd = min(yMinInd, min(x[[2]][, , 2]))
            yMaxInd = max(yMaxInd, max(x[[2]][, , 2]))
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
    if (!(is.null(xlim)) & !(is.null(ylim))) {
        xMinC1 = xlim[1]
        xMaxC1 = xlim[2]
        yMinC1 = ylim[1]
        yMaxC1 = ylim[2]
    }
    plot(NULL, pch = 3, xlim = c(xMinC1, xMaxC1), ylim = c(yMinC1,
        yMaxC1), asp = 1, xlab = xlab, ylab = ylab, axes = FALSE)
    box()
    abline(v = 0, lty = 3)
    abline(h = 0, lty = 3)
    points(individuals, col = individualsColors, pch = 3)
    if (individualsLabels == TRUE) {
        text(individuals[, 1], individuals[, 2] + offsetY1, labels = labelx,
            col = individualsColors, cex = cex)
    }
    if (is.null(listEllipsesByRep) && !is.null(individualsEllipses)) {
        for (i in 1:nb.ind) {
            lines(individualsEllipses[i, , ], col = individualsColors[i])
        }
    }
    if (!is.null(listEllipsesByRep)) {
        rep = 1
        for (x in listEllipsesByRep) {
            for (i in 1:nb.ind) {
                lines(x[[2]][i, , ], col = individualsColors[i])
            }
            points(x[[3]], col = individualsColors, pch = 15)
            if (individualsLabels == TRUE) {
                text(x[[3]][, 1], x[[3]][, 2] + offsetY1, labels = rep,
                  col = individualsColors, cex = 0.8)
            }
            for (i in 1:nb.ind) {
                currentProductCoord = individuals[i, ]
                currentRepCoord = x[[3]][i, ]
                segments(currentProductCoord[1], currentProductCoord[2],
                  currentRepCoord[1, 1], currentRepCoord[1, 2],
                  col = individualsColors[i])
            }
            rep = rep + 1
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
                col = variablesColors, cex = cex)
        }
        for (i in 1:nb.var) {
            arrows(variables[i, 1] * expand, variables[i, 2] *
                expand, 0, 0, code = 1, length = 0.1, col = variablesColors[i])
        }
    }
    else {
        plot(NULL, pch = 3, xlim = c(-1, 1), ylim = c(-1, 1),
            asp = 1, xlab = xlab, ylab = ylab, , axes = FALSE)
        box()
        x.cercle <- seq(-1, 1, by = 0.01)
        y.cercle <- sqrt(1 - x.cercle^2)
        lines(x.cercle, y = y.cercle)
        lines(x.cercle, y = -y.cercle)
        for (i in 1:nb.var) {
            arrows(variables[i, 1], variables[i, 2], 0, 0, code = 1,
                col = variablesColors[i], length = 0.1)
        }
        if (variablesLabels == TRUE) {
            text(variables[, 1], variables[, 2] + (sign(variables[,
                2]) * offsetY2), labels = labely, col = variablesColors,
                cex = cex)
        }
        abline(v = 0, lty = 3)
        abline(h = 0, lty = 3)
    }
    if (!biplot) {
        mtext(mainTitle, 3, line = 0.5, outer = T, cex = 1.8,
            font = 2)
        if (!is.null(subTitle)) {
            mtext(subTitle, 1, line = 1, outer = T, cex = 1.2,
                font = 3)
        }
    }
    else {
        mtext(mainTitle, 3, line = 1, outer = F, cex = 1.8, font = 2)
        if (!is.null(subTitle)) {
            mtext(subTitle, 1, line = 5, outer = F, cex = 1.2,
                font = 3)
        }
    }
}

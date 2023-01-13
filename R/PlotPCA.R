#'@title  PlotPCA
#' @param res.PCA  results of the function PCA
#' @param representation  "DistanceBiplot","TwoMaps","CorrelationBiplot".
#' @param nbAxes    "Auto"
#' @param title   custumizes the title of your PCA
#' @param confInt   the confidence limit for the ellipses
#' @param ellipsesType    "barycentric" for the confidence ellipses around the product means "individual" for the representation of a proportion of confInt of the panelist scores
#' @param returnX   reverses the abscisse axis
#' @param returnY  reverses the ordinate axis
#' @param type  "R" for a R graph, "png","pdf","jpeg" for registering the graph in these formats in your working directory
#' @param fileName  name of your file if you choose type="png","pdf" or "jpeg"
#' @param ellipsesCalculation "Chi", "F" or "Sas"
#' @param expand for biplot
#' @param variablesLabels    TRUE by default.
#'@param panellists   representation of the panellists on the map: "None" if the panellist are not represented on the map, "labels" if they are represented by their names, or "points" if they are represented by points
#'@param productColors vector of colors for products
#'@param variablesColors vector of colors for variables
#'@param xlim vector of two numbers indicating the limit of x axis
#'@param  ylim vector of two numbers indicating the limit of y axis
#'@examples{
  #'  data(cheeses)
#'  PlotPCA(PCA(cheeses,productName="ProductCode",subjectName="SubjectCode",
#'  replicate="Replicate",sessionName="Session"))
#'}
#' @export
PlotPCA <- function (res.PCA, panellists = "None", representation = "DistanceBiplot",
    nbAxes = "Auto", title = "PCA", confInt = 0.95, ellipsesType = "Barycentric",
    productColors = NULL, returnX = FALSE, returnY = FALSE, type = "R",
    fileName = "PCA", ellipsesCalculation = "Chi", expand = NULL,
    variablesLabels = TRUE, variablesColors = NULL, xlim = NULL,
    ylim = NULL)
{
    res.PlotPCA = list()
    if ((representation == "DistanceBiplot" | representation ==
        "CorrelationBiplot") & res.PCA$Option == "Correlation") {
        stop("Biplots only available with covariance option.")
        scaleUnit = FALSE
    }
    biplot = FALSE
    if (representation == "DistanceBiplot") {
        matsvd = svd(as.matrix(res.PCA$B))
        U = (matsvd$u)
        D = diag(matsvd$d)
        V = matsvd$v
        ud = U %*% D
        res.PCA$IndivCoord[, 1] = ud[, 1]
        res.PCA$IndivCoord[, 2] = ud[, 2]
        res.PCA$VarCoord[, 1] = V[, 1]
        res.PCA$VarCoord[, 2] = V[, 2]
        title = paste(title, " - ","Distance Biplot",
            sep = "")
        biplot = TRUE
    }
    if (representation == "CorrelationBiplot") {
        matsvd = svd(as.matrix(res.PCA$B))
        U = (matsvd$u)
        D = diag(matsvd$d)
        V = matsvd$v
        dv = (D) %*% t(V)
        u = U
        res.PCA$IndivCoord[, 1] = U[, 1]
        res.PCA$IndivCoord[, 2] = U[, 2]
        res.PCA$VarCoord[, 1] = dv[1, ]
        res.PCA$VarCoord[, 2] = dv[2, ]
        title = paste(title, " - ", "Correlation Biplot",
            sep = "")
        biplot = TRUE
    }
    res.PlotPCA[["IndivCoord"]] = res.PCA$IndivCoord
    res.PlotPCA[["VarCoord"]] = res.PCA$VarCoord
    if (nbAxes == "Auto" || nbAxes == "auto") {
        nbAxes = max(2, res.PCA$NbDimSig)
    }
    axes = combn(nbAxes, 2)
    if (is.null(ncol(axes))) {
        axes = matrix(c(1, 2), 2, 1)
    }
    res.PlotPCA[["axes"]] = axes
    for (i in 1:ncol(axes)) {
        main = title
        axe1 = axes[1, i]
        axe2 = axes[2, i]
        inertie1 = res.PCA$EigenValues[axe1]
        inertie2 = res.PCA$EigenValues[axe2]
        inertieCumul = round(100 * (inertie1 + inertie2)/sum(res.PCA$EigenValues),
            digits = 2)
        res.PlotPCA[["cumulInertia"]] = inertieCumul
        suppIndividuals = NULL
        if (ellipsesType != "None" | panellists != "None") {
            suppIndividuals = res.PCA$IndSup
        }
        individualsEllipses = NULL
        panelistCoord = NULL
        if (ellipsesType != "None" | panellists != "None") {
            calculationsEllipse = CalculateEllipses(suppIndividualTable = suppIndividuals,
                vep = res.PCA$EigenVectors, axes = c(axe1, axe2),
                confInt = confInt, ellipsesType = ellipsesType,
                ellipsesCalculation = ellipsesCalculation)
            panelistCoord = calculationsEllipse[[1]]
            individualsEllipses = calculationsEllipse[[2]]
            res.PlotPCA[["ellipses"]] = calculationsEllipse
        }
        if (is.null(productColors)) {
            productColors = rep("blue", nrow(res.PCA$IndivCoord))
        }
        suppIndividualsLabels = FALSE
        suppIndividualsToPlot = FALSE
        if (panellists != "None") {
            suppIndividualsToPlot = TRUE
            suppIndividualsLabels = FALSE
        }
        if (representation == "CorrelationBiplot") {
            suppIndividualsLabels = FALSE
            individualsEllipses = NULL
            main = paste(main, " - ", "Correlation Biplot",
                sep = "")
        }
        main = paste(main, " (", inertieCumul, "%)", sep = "")
        xlab = paste("Dim. ", axe1, " (", round(inertie1 * 100/sum(res.PCA$EigenValues),
            digits = 2), "%)", sep = "")
        ylab = paste("Dim. ", axe2, " (", round(inertie2 * 100/sum(res.PCA$EigenValues),
            digits = 2), "%)", sep = "")
        fun = call("INTERNAL_PlotPCACVA", variables = res.PCA$VarCoord[,
            c(axe1, axe2)], individuals = res.PCA$IndivCoord[,
            c(axe1, axe2)], suppIndividuals = panelistCoord,
            suppIndividualsToPlot = suppIndividualsToPlot, biplot = biplot,
            variablesColors = variablesColors, individualsColors = productColors,
            suppIndividualsColors = NULL, variablesLabels = variablesLabels,
            individualsLabels = TRUE, suppIndividualsLabels = suppIndividualsLabels,
            individualsEllipses = individualsEllipses, returnX = returnX,
            returnY = returnY, xlab = xlab, ylab = ylab, expand = expand,
            mainTitle = main, cex = 1.5, xlim = xlim, ylim = ylim)
        if (biplot) {
            filewidth = 7
        }
        else {
            filewidth = 14
        }
        GenericPlot(type = type, fileName = paste(fileName, "Axes",
            axe1, "&", axe2), CALLFUN = fun, filewidth = filewidth,
            fileheight = 7)
    }
    return(res.PlotPCA)
}

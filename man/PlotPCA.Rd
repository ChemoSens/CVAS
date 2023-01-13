\name{PlotPCA}
\alias{PlotPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PlotPCA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotPCA(res.PCA, panellists = "None", representation = "DistanceBiplot", nbAxes = "Auto", title = "PCA", confInt = 0.95, ellipsesType = "Barycentric", returnX = FALSE, returnY = FALSE, type = "R", fileName = "PCA", ellipsesCalculation = "Chi", expand = NULL, variablesLabels = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.PCA}{
results of the function PCA
}

  \item{representation}{
"DistanceBiplot","TwoMaps","CorrelationBiplot". 
}
  \item{nbAxes}{
"Auto"
}
  \item{title}{
custumizes the title of your PCA
}
  \item{confInt}{
the confidence limit for the ellipses
}
  \item{ellipsesType}{
"barycentric" for the confidence ellipses around the product means "individual" for the representation of a proportion of confInt of the panelist scores
}

  \item{returnX}{
reverses the abscisse axis
}
  \item{returnY}{
reverses the ordinate axis
}
  \item{type}{
"R" for a R graph, "png","pdf","jpeg" for registering the graph in these formats in your working directory
}
  \item{fileName}{
name of your file if you choose type="png","pdf" or "jpeg"
}
  \item{ellipsesCalculation}{
"Chi", "F" or "Sas"
}
  \item{expand}{
%%     ~~Describe \code{expand} here~~
}
  \item{variablesLabels}{
TRUE by default. 
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
caroline.peltier@dijon.inra.fr
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(cheeses)
PlotPCA(PCA(cheeses,productName="ProductCode",subjectName="SubjectCode",replicate="Replicate",sessionName="Session"))
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (res.PCA, panellists = "None", representation = "DistanceBiplot", 
    nbAxes = "Auto", title = "PCA", confInt = 0.95, ellipsesType = "Barycentric", 
    productColors = NULL, returnX = FALSE, returnY = FALSE, type = "R", 
    fileName = "PCA", ellipsesCalculation = "Chi", expand = NULL, 
    variablesLabels = TRUE, variablesColors = NULL) 
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
        ud = U \%*\% D
        res.PCA$IndivCoord[, 1] = ud[, 1]
        res.PCA$IndivCoord[, 2] = ud[, 2]
        res.PCA$VarCoord[, 1] = V[, 1]
        res.PCA$VarCoord[, 2] = V[, 2]
        title = paste(title, " - ", TS_GetLabel("Distance Biplot"), 
            sep = "")
        biplot = TRUE
    }
    if (representation == "CorrelationBiplot") {
        matsvd = svd(as.matrix(res.PCA$B))
        U = (matsvd$u)
        D = diag(matsvd$d)
        V = matsvd$v
        dv = (D) \%*\% t(V)
        u = U
        res.PCA$IndivCoord[, 1] = U[, 1]
        res.PCA$IndivCoord[, 2] = U[, 2]
        res.PCA$VarCoord[, 1] = dv[1, ]
        res.PCA$VarCoord[, 2] = dv[2, ]
        title = paste(title, " - ", TS_GetLabel("Correlation Biplot"), 
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
            main = paste(main, " - ", TS_GetLabel("Correlation Biplot"), 
                sep = "")
        }
        main = paste(main, " (", inertieCumul, "\%)", sep = "")
        xlab = paste("Dim. ", axe1, " (", round(inertie1 * 100/sum(res.PCA$EigenValues), 
            digits = 2), "\%)", sep = "")
        ylab = paste("Dim. ", axe2, " (", round(inertie2 * 100/sum(res.PCA$EigenValues), 
            digits = 2), "\%)", sep = "")
        fun = call("PlotPCACVA", variables = res.PCA$VarCoord[, 
            c(axe1, axe2)], individuals = res.PCA$IndivCoord[, 
            c(axe1, axe2)], suppIndividuals = panelistCoord, 
            suppIndividualsToPlot = suppIndividualsToPlot, biplot = biplot, 
            variablesColors = variablesColors, individualsColors = productColors, 
            suppIndividualsColors = NULL, variablesLabels = variablesLabels, 
            individualsLabels = TRUE, suppIndividualsLabels = suppIndividualsLabels, 
            individualsEllipses = individualsEllipses, returnX = returnX, 
            returnY = returnY, xlab = xlab, ylab = ylab, expand = expand, 
            mainTitle = main)
        if (biplot) {
            filewidth = 7
        }
        else {
            filewidth = 14
        }
        GenericPlot(type = type, fileName = paste(fileName, TS_GetLabel("Axes"), 
            axe1, "&", axe2), CALLFUN = fun, filewidth = filewidth, 
            fileheight = 7)
    }
    return(res.PlotPCA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

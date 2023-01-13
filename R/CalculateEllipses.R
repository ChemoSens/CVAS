#'@import ellipse
#'@importFrom stats cov qf qchisq
CalculateEllipses <-
function (suppIndividualTable, vep, axes = c(1, 2), confInt = 0.9,
    ellipsesType = "barycentric", productName = "ProductCode",
    subjectName = "SubjectCode", ellipsesCalculation = "Chi")
{
    products = levels(as.factor(as.character(suppIndividualTable[,
        productName])))
    nb.prod = length(products)
    suj = levels(as.factor(suppIndividualTable[, subjectName]))
    nb.suj = length(suj)
    ell = array(dim = c(nb.prod, 100, 2))
    ell2 = array(dim = c(nb.prod, 100, 2))
    coord.pts = list()
    attributs = colnames(suppIndividualTable)[-c(1, 2)]
    nbAttributes = length(attributs)
    ax1 = axes[1]
    ax2 = axes[2]
    coorCentre = data.frame()
    for (i in 1:nb.prod) {
        indiv.sup = as.matrix(suppIndividualTable[suppIndividualTable[,
            productName] == products[i], 3:ncol(suppIndividualTable)])
        coord.indiv.sup = indiv.sup %*% vep
        coorCentre[i, 1] = mean(coord.indiv.sup[, ax1], na.rm = TRUE)
        coorCentre[i, 2] = mean(coord.indiv.sup[, ax2], na.rm = TRUE)
        coorCentre[i, productName] = products[i]
        coord.pts[[i]] = coord.indiv.sup[, c(ax1, ax2)]
        ddlN = length(coord.pts[[i]][, 1])
        matCov = cov(coord.pts[[i]])
        matCov2 = matCov/dim(coord.indiv.sup)[1]
        if (ellipsesType != "individual") {
            covarianceMatrix = matCov2
        }
        if (ellipsesType == "individual") {
            covarianceMatrix = matCov
        }
        quant = sqrt(2 * nb.suj * qf(confInt, 2, nb.suj - 2)/(nb.suj -
            2))
        if (ellipsesCalculation == "Chi") {
            quant = sqrt(qchisq(confInt, 2))
        }
        if (ellipsesCalculation == "F") {
            quant = sqrt(2 * nb.suj * qf(confInt, 2, nb.suj -
                2)/(nb.suj - 2))
        }
        if (ellipsesCalculation == "Sas") {
            quant = sqrt((2 * nb.suj * qf(confInt, 2, nb.suj -
                2)/(nb.suj - 2)) * (nb.suj - 1)/nb.suj)
        }
        ell[i, , 1] = ellipse::ellipse(x = covarianceMatrix,
            centre = c(coorCentre[i, 1], coorCentre[i, 2]), level = confInt,
            t = quant)[, 1]
        ell[i, , 2] = ellipse::ellipse(x = covarianceMatrix,
            centre = c(coorCentre[i, 1], coorCentre[i, 2]), level = confInt,
            t = quant)[, 2]
        ell[i, , 1] = ellipse::ellipse(covarianceMatrix, centre = c(coorCentre[i,
            1], coorCentre[i, 2]), level = confInt, t = quant)[,
            1]
        ell[i, , 2] = ellipse::ellipse(covarianceMatrix, centre = c(coorCentre[i,
            1], coorCentre[i, 2]), level = confInt, t = quant)[,
            2]
    }
    L = list(coord.pts, ell, coorCentre)
    return(L)
}

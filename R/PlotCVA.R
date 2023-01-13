PlotCVA <-
function (res.CVA, panellists = "none", confInt = 0.9, ellipsesType = "barycentric", 
    productColors = NULL, returnX = FALSE, returnY = FALSE, linkBetweenProducts = TRUE, 
    suppIndividualsToPlot = FALSE, type = "R", fileName = "CVA", 
    mainTitle = "CVA", ellipsesByRep = F, cex = 0.8, ellipsesCalculation = "Chi") 
{
    axes = combn(res.CVA$nbAxes, 2)
    eigVec = res.CVA$EigenVectors
    if (is.null(res.CVA$option)) {
        option = res.CVA$modelType
    }
    else {
        option = res.CVA$option
    }
    if (res.CVA$representation == "biplot") {
        biplotBool = TRUE
    }
    else {
        biplotBool = FALSE
    }
    if (is.null(ncol(axes))) {
        axes = matrix(c(1, 2), 2, 1)
    }
    if (type != "none") {
        for (i in 1:ncol(axes)) {
            axe1 = axes[1, i]
            axe2 = axes[2, i]
            inertie1 = res.CVA$EigenValues[axe1]
            inertie2 = res.CVA$EigenValues[axe2]
            inertieCumul = round(100 * (inertie1 + inertie2)/sum(res.CVA$EigenValues), 
                digits = 2)
            if (ellipsesType != "none" || panellists != "none") {
                suppIndividuals = res.CVA$CenteredProductSubjectTable
            }
            if (ellipsesType != "none") {
                if (option != "mam" && option != "overall" && 
                  option != "tws") {
                  calculationEllipses = CalculateEllipses(suppIndividualTable = res.CVA$CenteredProductSubjectTable, 
                    vep = res.CVA$EigenVectors, axes = c(axe1, 
                      axe2), confInt = confInt, ellipsesType = ellipsesType, 
                    productName = "ProductCode", ellipsesCalculation = ellipsesCalculation)
                }
                if (option == "mam" || option == "overall") {
                  res.CVA$decomposition[, "scoreWithoutScaling"] = res.CVA$decomposition[, 
                    "prodEffect"] + res.CVA$decomposition[, "Disag"]
                  pureData = gettingAppropriateData(res.CVA$decomposition, 
                    scoreName = "scoreWithoutScaling", subjectName = "ass", 
                    replicateName = "rep", productName = "prod", 
                    attributeName = "Attribute")
                  matrixPureData = gettingAppropriateData(res.CVA$decomposition, 
                    scoreName = "scoreWithoutScaling", subjectName = "ass", 
                    replicateName = "rep", productName = "prod", 
                    attributeName = "Attribute", asMatrix = TRUE)
                  newColnamesPureData = colnames(pureData)
                  newColnamesPureData[which(colnames(pureData) == 
                    "ass")] = "SubjectCode"
                  newColnamesPureData[which(colnames(pureData) == 
                    "prod")] = "ProductCode"
                  colnames(pureData) = newColnamesPureData
                  wDemi = res.CVA$wDemi
                  if (res.CVA$representation == "biplot") {
                    pureData[, 3:dim(pureData)[2]] = matrixPureData %*% 
                      solve(wDemi)
                    calculationEllipses = CalculateEllipses(suppIndividualTable = pureData, 
                      vep = eigVec, axes = c(axe1, axe2), confInt = confInt, 
                      ellipsesType = ellipsesType, productName = "ProductCode", 
                      subjectName = "SubjectCode", ellipsesCalculation = ellipsesCalculation)
                  }
                  else {
                    calculationEllipses = CalculateEllipses(suppIndividualTable = pureData, 
                      vep = eigVec, axes = c(axe1, axe2), confInt = confInt, 
                      ellipsesType = ellipsesType, productName = "ProductCode", 
                      subjectName = "SubjectCode", ellipsesCalculation = ellipsesCalculation)
                  }
                }
                if (option == "tws") {
                  res.CVA$decomposition[, "scoreWithoutSubject"] = res.CVA$decomposition[, 
                    "prodEffect"] + res.CVA$decomposition[, "int"]
                  matrixPureData = gettingAppropriateData(res.CVA$decomposition, 
                    scoreName = "scoreWithoutSubject", subjectName = "ass", 
                    replicateName = "rep", productName = "prod", 
                    attributeName = "Attribute", asMatrix = TRUE)
                  pureData = gettingAppropriateData(res.CVA$decomposition, 
                    scoreName = "scoreWithoutSubject", subjectName = "ass", 
                    replicateName = "rep", productName = "prod", 
                    attributeName = "Attribute")
                  newColnamesPureData = colnames(pureData)
                  newColnamesPureData[which(colnames(pureData) == 
                    "ass")] = "SubjectCode"
                  newColnamesPureData[which(colnames(pureData) == 
                    "prod")] = "ProductCode"
                  colnames(pureData) = newColnamesPureData
                  if (res.CVA$representation == "biplot") {
                    wDemi = res.CVA$wDemi
                    pureData[, 3:dim(pureData)[2]] = matrixPureData %*% 
                      solve(wDemi)
                    calculationEllipses = CalculateEllipses(suppIndividualTable = pureData, 
                      vep = eigVec, axes = c(axe1, axe2), confInt = confInt, 
                      ellipsesType = ellipsesType, productName = "ProductCode", 
                      subjectName = "SubjectCode", ellipsesCalculation = ellipsesCalculation)
                  }
                  else {
                    calculationEllipses = CalculateEllipses(suppIndividualTable = pureData, 
                      vep = eigVec, axes = c(axe1, axe2), confInt = confInt, 
                      ellipsesType = ellipsesType, productName = "ProductCode", 
                      subjectName = "SubjectCode", ellipsesCalculation = ellipsesCalculation)
                  }
                }
                panelistCoord = calculationEllipses[[1]]
                individualsEllipses = calculationEllipses[[2]]
            }
            else {
                individualsEllipses = NULL
            }
            if (linkBetweenProducts) {
                individualsSegments = res.CVA$HotellingTable > 
                  1 - confInt
            }
            else {
                individualsSegments = NULL
            }
            if (panellists == "labels" || panellists == "points") {
                suppIndividualToPlot = TRUE
                suppIndividualsLabels = FALSE
            }
            else {
                suppIndividualsLabels = FALSE
                suppIndividualToPlot = FALSE
            }
            listEllipsesByRep = NULL
            if (ellipsesByRep == T) {
                listEllipsesByRep = list()
                extendedData = res.CVA$ExtendedData
                m = apply(extendedData[, -c(1:3)], 2, "scale", 
                  center = TRUE, scale = FALSE)
                m = as.data.frame(m)
                m = cbind(as.character(extendedData$ProductCode), 
                  m)
                m = cbind(as.character(extendedData$SubjectCode), 
                  m)
                m[, -c(1, 2)] = as.matrix(m[, -c(1, 2)]) %*% 
                  solve(res.CVA$wDemi)
                colnames(m)[1] = "SubjectCode"
                colnames(m)[2] = "ProductCode"
                for (i in unique(extendedData$Replicate)) {
                  listEllipsesByRep[[i]] = CalculateEllipses(suppIndividualTable = m[extendedData$Replicate == 
                    i, ], vep = res.CVA$EigenVectors, axes = c(axe1, 
                    axe2), confInt = confInt, ellipsesType = ellipsesType, 
                    productName = "ProductCode")
                }
            }
            title = paste(mainTitle, " (", inertieCumul, "%)", 
                sep = "")
            xlab = paste("Can. ", axe1, " (", round(inertie1 * 
                100/sum(res.CVA$EigenValues), digits = 2), "%)", 
                sep = "")
            ylab = paste("Can. ", axe2, " (", round(inertie2 * 
                100/sum(res.CVA$EigenValues), digits = 2), "%)", 
                sep = "")
            subtitle = paste("NDIMSIG=", res.CVA$NbDimSig, ", Hotelling Lawley stat=", 
                round(res.CVA$Stats$Stat, digits = 3), ", F=", 
                round(res.CVA$Stats$F, digits = 3), " (p=", FriendlyPValue(res.CVA$Stats$pval), 
                ")\n", "Confidence ellipses", "=", confInt * 
                  100, "%", sep = "")
            if (biplotBool) {
                graphWidth = 7
            }
            else {
                graphWidth = 14
            }
            fun = call("INTERNAL_PlotPCACVA", variables = res.CVA$VarCoord[, 
                c(axe1, axe2)], individuals = res.CVA$IndivCoord[, 
                c(axe1, axe2)], suppIndividuals = panelistCoord, 
                suppIndividualsToPlot = suppIndividualToPlot, 
                biplot = biplotBool, variablesColors = NULL, 
                individualsColors = productColors, suppIndividualsColors = NULL, 
                variablesLabels = TRUE, individualsLabels = TRUE, 
                suppIndividualsLabels = suppIndividualsLabels, 
                individualsEllipses = individualsEllipses, returnX = returnX, 
                returnY = returnY, individualsSegments = individualsSegments, 
                cex = cex, xlab = xlab, ylab = ylab, expand = NULL, 
                subTitle = subtitle, mainTitle = title, listEllipsesByRep = listEllipsesByRep)
            GenericPlot(type = type, fileName = paste(fileName, 
                TS_GetLabel("axes"), axe1, "&", axe2), CALLFUN = fun, 
                filewidth = 7, fileheight = 7)
        }
    }
}

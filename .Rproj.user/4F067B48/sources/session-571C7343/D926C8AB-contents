CVA <-
function (data, test = "Hotelling-Lawley", option = "tw", hotellingTableBool = TRUE,
    nbAxes = "auto", alpha = 0.1, productName = "product", subjectName = "subject",
    replicateName = "rep", sessionName = "session", representation = "biplot")
{
    getNumberOfSignificantDimensionsOfCVA = function(eigVal,
        I, ddlW, p, alpha = 0.05) {
        significance = TRUE
        n = length(eigVal)
        k = 0
        K = min(p, I)
        while (significance == TRUE) {
            chi2 = ddlW * sum(eigVal[(k + 1):K], na.rm = TRUE)
            df = (p - k) * (I - k - 1)
            if (df < 1) {
                return(k)
            }
            proba = pchisq(chi2, df = df, lower.tail = FALSE)
            if (proba > alpha) {
                significance = FALSE
            }
            k = k + 1
        }
        return(k - 1)
    }
    if (!(option %in% c("ow", "tw", "mam", "overall", "tws"))) {
        stop("option does not exist. Please choose between 'tw' (Two-Way ANOVA),  'ow' (One-Way ANOVA), mam (for MAM model) or 'overall' (for the overall scaling)")
    }
    if (!productName %in% colnames(data)) {
        stop("productName is not in the colnames of data")
    }
    if (!subjectName %in% colnames(data)) {
        stop("subjectName is not in the colnames of data")
    }
    if (!replicateName %in% colnames(data)) {
        stop("replicateName is not in the colnames of data")
    }
    if (!representation %in% c("biplot", "TwoMaps")) {
        stop("representation should be 'biplot' or 'TwoMaps'")
    }
    if (!test %in% c("Hotelling-Lawley", "Wilks", "Pillai", "Roy")) {
        stop("test should be 'Hotelling-Lawley','Roy','Wilks' or 'Pillai'.")
    }
    Attributes = colnames(data)[colnames(data) != productName &
    colnames(data) != subjectName & colnames(data) != replicateName &
    colnames(data) != sessionName]
    replicates = levels(as.factor(data[, replicateName]))
    nbReplicates = length(replicates)
    nbAttributes = length(Attributes)
    products = levels(as.factor(data[, productName]))
    subjects = levels(as.factor(data[, subjectName]))
    nbSubjects = length(subjects)
    nbProducts = length(products)
    wDemi = NULL
    if (nbAttributes < 2) {
        stop("Insufficient number of attributes")
    }
    if ((nbSubjects - 1) * (nbProducts - 1) < nbAttributes) {
        stop("Insufficient number of individuals")
    }
    CenteredProductTable = data.frame(NULL)
    CenteredProductTable = aggregate(data[, Attributes[1]], by = list(data[,
        productName]), FUN = "mean", na.rm = TRUE)
    for (i in 2:length(Attributes)) {
        CenteredProductTable[, Attributes[i]] = aggregate(data[,
            Attributes[i]], by = list(data[, productName]), FUN = "mean",
            na.rm = TRUE)[, 2]
    }
    colnames(CenteredProductTable) = c(productName, Attributes)
    CenteredProductTable[, Attributes] = apply(CenteredProductTable[,
        Attributes], 2, "scale", center = TRUE, scale = FALSE)
    matrixOfCenteredProduct = as.matrix(CenteredProductTable[,
        Attributes])
    rownames(matrixOfCenteredProduct) = CenteredProductTable[,
        productName]
    CenteredProductSubjectTable = data.frame(NULL)
    CenteredProductSubjectTable = aggregate(data[, Attributes[1]],
        by = list(data[, productName], data[, subjectName]),
        FUN = "mean", na.rm = TRUE)
    for (i in 2:length(Attributes)) {
        CenteredProductSubjectTable[, Attributes[i]] = aggregate(data[,
            Attributes[i]], by = list(data[, productName], data[,
            subjectName]), FUN = "mean", na.rm = TRUE)[, 3]
    }
    colnames(CenteredProductSubjectTable) = c(productName, subjectName,
        Attributes)
    CenteredProductSubjectTable[, Attributes] = apply(CenteredProductSubjectTable[,
        Attributes], 2, "scale", center = TRUE, scale = FALSE)
    matrixOfCenteredProductSubject = as.matrix(CenteredProductSubjectTable[,
        Attributes])
    rownames(matrixOfCenteredProductSubject) = paste(CenteredProductSubjectTable[,
        productName], CenteredProductSubjectTable[, subjectName],
        sep = "")
    CenteredSubjectTable = data.frame(NULL)
    CenteredSubjectTable = aggregate(data[, Attributes[1]], by = list(data[,
        subjectName]), FUN = "mean", na.rm = TRUE)
    for (i in 2:length(Attributes)) {
        CenteredSubjectTable[, Attributes[i]] = aggregate(data[,
            Attributes[i]], by = list(data[, subjectName]), FUN = "mean",
            na.rm = TRUE)[, 2]
    }
    colnames(CenteredSubjectTable) = c(subjectName, Attributes)
    CenteredSubjectTable[, Attributes] = apply(CenteredSubjectTable[,
        Attributes], 2, "scale", center = TRUE, scale = FALSE)
    matrixOfCenteredSubject = as.matrix(CenteredSubjectTable[,
        Attributes])
    rownames(matrixOfCenteredSubject) = CenteredSubjectTable[,
        subjectName]
    decomposition = NULL
    statTest = NULL
    statF = NULL
    statPval = NULL
    if (option == "ow") {
        SSProd = t(matrixOfCenteredProduct) %*% matrixOfCenteredProduct
        SSTotal = t(matrixOfCenteredProductSubject) %*% matrixOfCenteredProductSubject
        SSres = SSTotal - SSProd
        if (det(SSres) <= 2.2e-16) {
            L = list()
            L[[1]] = det(SSres)
            names(L) = c("Det")
            return(L)
            stop("Non invertible matrix")
        }
        Y = as.matrix(CenteredProductSubjectTable[, Attributes])
        Product = as.factor(CenteredProductSubjectTable[, productName])
        Subject = as.factor(CenteredProductSubjectTable[, subjectName])
        performances = summary(manova(lm(Y ~ Product)), test = test)$stats
        statTest = performances[1, 2]
        statF = performances[1, 3]
        statPval = performances[1, 6]
    }
    else {
        if (option == "tw") {
            SSProd = nbSubjects * t(matrixOfCenteredProduct) %*%
                matrixOfCenteredProduct
            SSTotal = t(matrixOfCenteredProductSubject) %*% matrixOfCenteredProductSubject
            SSsujet = nbProducts * t(matrixOfCenteredSubject) %*%
                matrixOfCenteredSubject
            SSres = SSTotal - SSsujet - SSProd
            if (det(SSres) <= 2.2e-16) {
                L = list()
                L[[1]] = det(SSres)
                names(L) = c("Det")
                return(L)
                stop("Non invertible matrix")
            }
            dataExtended = data[, c(subjectName, productName,
                replicateName, Attributes)]
            tryCatch({
                resPerf = PanelPerformances(frame = dataExtended,
                  modelType = "classic", negativeCorrection = FALSE,
                  correctOnlyIfSignificant = FALSE, limitOfSignificance = 0.05,
                  onlySignificantDim = FALSE, manovaTest = "Hotelling")
                decomposition = resPerf$decomposition
                performances = resPerf$multiPanelPerformances
                statTest = performances["discrimination", "stat"]
                statF = performances["discrimination", "f"]
                statPval = performances["discrimination", "pvalue"]
            })
        }
        if (option == "mam") {
            SSProd = nbReplicates * nbSubjects * t(matrixOfCenteredProduct) %*%
                matrixOfCenteredProduct
            dataExtended = data[, c(subjectName, productName,
                replicateName, Attributes)]
            resMAM = PanelPerformances(frame = dataExtended,
                modelType = "mam", negativeCorrection = FALSE,
                correctOnlyIfSignificant = FALSE, limitOfSignificance = 0.05,
                onlySignificantDim = FALSE, manovaTest = "Hotelling")
            SSres = resMAM$matW
            if (det(SSres) <= 2.2e-16) {
                L = list()
                L[[1]] = det(SSres)
                names(L) = c("Det")
                return(L)
                stop("Non invertible matrix")
            }
            decomposition = resMAM$decomposition
            performances = resMAM$multiPanelPerformances
            statTest = performances["discrimination", "stat"]
            statF = performances["discrimination", "f"]
            statPval = performances["discrimination", "pvalue"]
        }
        if (option == "overall") {
            SSProd = nbReplicates * nbSubjects * t(matrixOfCenteredProduct) %*%
                matrixOfCenteredProduct
            dataExtended = data[, c(subjectName, productName,
                replicateName, Attributes)]
            resMAM = PanelPerformances(frame = dataExtended,
                modelType = "overall", negativeCorrection = FALSE,
                correctOnlyIfSignificant = FALSE, limitOfSignificance = 0.05,
                onlySignificantDim = FALSE, manovaTest = "Hotelling")
            SSres = resMAM$matW
            if (det(SSres) <= 2.2e-16) {
                L = list()
                L[[1]] = det(SSres)
                names(L) = c("Det")
                return(L)
                stop("Non invertible matrix")
            }
            decomposition = resMAM$decomposition
            performances = resMAM$multiPanelPerformances
            statTest = performances["discrimination", "stat"]
            statF = performances["discrimination", "f"]
            statPval = performances["discrimination", "pvalue"]
        }
        if (option == "tws") {
            SSProd = nbSubjects * t(matrixOfCenteredProduct) %*%
                matrixOfCenteredProduct
            SSTotal = t(matrixOfCenteredProductSubject) %*% matrixOfCenteredProductSubject
            SSsujet = nbProducts * t(matrixOfCenteredSubject) %*%
                matrixOfCenteredSubject
            SSres = SSTotal - SSsujet - SSProd
            if (det(SSres) <= 2.2e-16) {
                L = list()
                L[[1]] = det(SSres)
                names(L) = c("Det")
                return(L)
                stop("Non invertible matrix")
            }
            dataExtended = data[, c(subjectName, productName,
                replicateName, Attributes)]
            resPerf = PanelPerformances(frame = dataExtended,
                modelType = "classic", manovaTest = "Hotelling",
                onlySignificantDim = FALSE)
            decomposition = resPerf$decomposition
            performances = resPerf$multiPanelPerformances
            statTest = performances["discrimination", "stat"]
            statF = performances["discrimination", "f"]
            statPval = performances["discrimination", "pvalue"]
        }
    }
    if (is.list(option)) {
        SSProd = option[[1]]
        SSres = option[[2]]
    }
    if (det(SSres) <= 2.2e-16) {
        L = list()
        L[[1]] = det(SSres)
        names(L) = c("Det")
        return(L)
        stop("Non invertible matrix")
    }
    ConditioningOfW = sqrt(max(Re(eigen(SSres)$values))/min(Re(eigen(SSres)$values)))
    discriminationMatrix = solve(SSres) %*% SSProd
    eigVal = Re(eigen(discriminationMatrix)$values)
    vecOrtho = Re(eigen(discriminationMatrix)$vectors)
    rownames(vecOrtho) = rownames(SSProd)
    if (option == "tw" || option == "tws") {
        nbDimSig = getNumberOfSignificantDimensionsOfCVA(eigVal = eigVal,
            I = nbProducts, ddlW = (nbProducts - 1) * (nbSubjects -
                1), p = nbAttributes, alpha = 0.05)
    }
    if (option == "ow") {
        nbDimSig = getNumberOfSignificantDimensionsOfCVA(eigVal = eigVal,
            I = nbProducts, ddlW = nbSubjects * (nbProducts -
                1), p = nbAttributes, alpha = 0.05)
    }
    if (option == "overall") {
        nbDimSig = getNumberOfSignificantDimensionsOfCVA(eigVal = eigVal,
            I = nbProducts, ddlW = (nbSubjects - 1) * (nbProducts -
                2), p = nbAttributes, alpha = 0.05)
    }
    if (option == "mam") {
        nbDimSig = getNumberOfSignificantDimensionsOfCVA(eigVal = eigVal,
            I = nbProducts, ddlW = (nbSubjects - 1) * (nbProducts -
                2), p = nbAttributes, alpha = 0.05)
    }
    if (nbAxes == "auto") {
        nbAxes = max(2, nbDimSig)
    }
    else if (is.numeric(nbAxes)) {
        nbAxes = min(nbAxes, min(nbProducts - 1, nbAttributes -
            1))
    }
    normW = function(v, W) {
        return(sqrt(t(v) %*% W %*% v))
    }
    norms = apply(vecOrtho, 2, normW, SSres)
    eigVec = vecOrtho
    for (i in 1:ncol(eigVec)) {
        eigVec[, i] = vecOrtho[, i]/norms[i]
    }
    individuals = matrixOfCenteredProduct %*% eigVec
    rownames(individuals) = rownames(matrixOfCenteredProduct)
    variables = cor(matrixOfCenteredProduct, individuals)
    biplot = FALSE
    if (representation == "biplot") {
        rootAndVap = function(matSymDefPos) {
            eig = eigen(matSymDefPos)
            rac = eig$vectors %*% sqrt(diag(eig$values)) %*%
                t((eig$vectors))
            return(rac)
        }
        wDemi = rootAndVap(SSres)
        yMat = matrixOfCenteredProduct %*% solve(wDemi)
        individuals = matrix(NA, nbProducts, nbAxes)
        rownames(individuals) = rownames(matrixOfCenteredProduct)
        variables = matrix(NA, nbAttributes, nbAxes)
        rownames(variables) = colnames(matrixOfCenteredProduct)
        matsvd = svd(yMat)
        U = (matsvd$u)
        D = diag(matsvd$d)
        V = matsvd$v
        eigVec = V
        ud = U %*% D
        individuals[, 1:nbAxes] = ud[, 1:nbAxes]
        variables[, 1:nbAxes] = V[, 1:nbAxes]
        CenteredProductSubjectTable[, 3:dim(CenteredProductSubjectTable)[2]] = matrixOfCenteredProductSubject %*%
            solve(wDemi)
    }
    statResults = list()
    statResults[["Test"]] = test
    statResults[["Stat"]] = statTest
    statResults[["F"]] = statF
    statResults[["pval"]] = statPval
    tabtot = NULL
    if (hotellingTableBool) {
        if (option == "mam" || option == "overall") {
            decomposition[, "scoreWithoutScaling"] = decomposition[,
                "prodEffect"] + decomposition[, "Disag"]
            pureData = gettingAppropriateData(decomposition,
                scoreName = "scoreWithoutScaling", subjectName = "ass",
                replicateName = "rep", productName = "prod",
                attributeName = "Attribute")
            newColnamesPureData = colnames(pureData)
            newColnamesPureData[which(colnames(pureData) == "ass")] = "SubjectCode"
            newColnamesPureData[which(colnames(pureData) == "prod")] = "ProductCode"
            colnames(pureData) = newColnamesPureData
            tabtot = NULL
            tabtot = hotellingTable(matCva = pureData, vep = eigVec[,
                1:nbAxes], axes = c(1:nbAxes), colAttributes = 3:dim(CenteredProductSubjectTable)[2])
        }
        if (option == "tws") {
            decomposition[, "scoreWithoutSubject"] = decomposition[,
                "prodEffect"] + decomposition[, "int"]
            pureData = gettingAppropriateData(decomposition,
                scoreName = "scoreWithoutSubject", subjectName = "ass",
                replicateName = "rep", productName = "prod",
                attributeName = "Attribute")
            newColnamesPureData = colnames(pureData)
            newColnamesPureData[which(colnames(pureData) == "ass")] = "SubjectCode"
            newColnamesPureData[which(colnames(pureData) == "prod")] = "ProductCode"
            colnames(pureData) = newColnamesPureData
            tabtot = NULL
            tabtot = hotellingTable(matCva = pureData, vep = eigVec[,
                1:nbAxes], axes = c(1:nbAxes), colAttributes = 3:dim(CenteredProductSubjectTable)[2])
        }
        if (option != "mam" && option != "overall" && option !=
            "tws") {
            tabtot = hotellingTable(matCva = CenteredProductSubjectTable,
                vep = eigVec[, 1:nbAxes], axes = c(1:nbAxes),
                colAttributes = 3:dim(CenteredProductSubjectTable)[2],
                productName = productName)
        }
    }
    resultList = list(IndivCoord = individuals, VarCoord = variables,
        NbDimSig = nbDimSig, HotellingTable = tabtot, ConditioningOfW = ConditioningOfW,
        B = SSProd, W = SSres, EigenVectors = eigVec, EigenValues = eigVal,
        Stats = statResults, decomposition = decomposition, CenteredProductSubjectTable = CenteredProductSubjectTable,
        nbAxes = nbAxes, wDemi = wDemi, option = option, representation = representation,
        data = data)
    return(resultList)
}

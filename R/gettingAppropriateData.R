#' @importFrom stats reshape aggregate
gettingAppropriateData <-
function (selectedData, asMatrix = FALSE, scoreName = "Score",
    productName = "ProductCode", subjectName = "SubjectCode",
    attributeName = "AttributeCode", replicateName = "Replicate",
    replaceMode = "crossmean", scaleUnit = FALSE)
{
    selectedData = selectedData[, c(productName, subjectName,
        replicateName, attributeName, scoreName)]
    extendedData = reshape(selectedData, idvar = c(productName,
        subjectName, replicateName), timevar = attributeName,
        direction = "wide")
    attributeNames = substr(colnames(extendedData)[4:dim(extendedData)[2]],
        nchar(scoreName) + 2, 1000)
    newColnames = c(productName, subjectName, replicateName,
        attributeNames)
    colnames(extendedData) = newColnames
    appropriateData = data.frame(NULL)
    appropriateData = aggregate(extendedData[, attributeNames[1]],
        by = list(extendedData[, productName], extendedData[,
            subjectName]), FUN = "mean", na.rm = TRUE)
    for (i in 2:length(attributeNames)) {
        appropriateData[, attributeNames[i]] = aggregate(extendedData[,
            attributeNames[i]], by = list(extendedData[, productName],
            extendedData[, subjectName]), FUN = "mean", na.rm = TRUE)[,
            3]
    }
    colnames(appropriateData) = c(productName, subjectName, attributeNames)
    appropriateData[, attributeNames] = apply(appropriateData[,
        attributeNames], 2, "scale", center = TRUE, scale = scaleUnit)
    if (asMatrix) {
        resultingTable = as.matrix(appropriateData[, attributeNames])
        rownames(resultingTable) = paste(appropriateData[, productName],
            appropriateData[, subjectName], sep = "")
    }
    else {
        resultingTable = appropriateData
    }
    return(resultingTable)
}

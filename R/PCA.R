
#' @title  PCA
#' @description  Calculates the PCA of a dataset
#' @param data The dataset should have the following columns: subject, product, rep, attribute 1, ... attribute n
#' @param option "Covariance" or "Correlation". Covariance corresponds to the un-normalized PCA (thus, a not discriminant attribute will have a week weight in the analysis) whereas Correlation corresponds to the normalized PCA (All attributes have the same weight in the analysis)
#' @param productName name of the column containing the products. By default, "product"
#' @param subjectName name of the column containing the subjects. By default, "subject"
#' @param replicateName name of the column containing the replicates. By default, "rep"
#' @param sessionName if necessary, name of the column containing the sessions. By default "session"
#' @return{
#' #'\itemize{
#'  \item{ B }{matrix of product effect covariance}
#'  \item{ IndSup }{dataframe containing the centered scores obtained for each product*subject evaluation}
#'  \item{ EigenVectors }{eigen vectors from the diagonalization of B}
#'  \item{ EigenValues }{eigen values from the diagonalization of B}
#'  \item{ IndivCoord }{coordinates of the products (individuals) in the PCA map}
#'  \item{ VarCoord }{coordinates of the attributes (variables) in the PCA map}
#'  \item{ NbdimSig }{number of significant dimensions}
#'  \item{ Option }{choosen option}
#'}
#'}
#' @seealso{
#'    \code{\link{PlotPCA}},\code{\link{CVA}}, \code{link{PlotCVA}}
#'  }
#' @examples{
#'  data(cheeses)
#'  PCA(cheeses, productName="ProductCode",subjectName="SubjectCode",
#'  replicateName="Replicate",sessionName="Session")
#' }
#' @importFrom stats aggregate cor sd
#' @export
PCA <-
function (data, option = "Covariance", productName = "product",
    subjectName = "subject", replicateName = "rep", sessionName = "session")
{
    scaleUnit = FALSE
    if (option == "Correlation") {
        scaleUnit = TRUE
    }
    if (option == "Covariance") {
        scaleUnit = FALSE
    }
    attributs = colnames(data)[colnames(data) != productName &
        colnames(data) != subjectName & colnames(data) != replicateName &
        colnames(data) != sessionName]
    nbAttributes = length(attributs)
    CenteredProductMeansTable = data.frame(NULL)
    CenteredProductMeansTable = aggregate(data[, attributs[1]],
        by = list(data[, productName]), FUN = "mean", na.rm = TRUE)
    for (i in 2:length(attributs)) {
        CenteredProductMeansTable[, attributs[i]] = aggregate(data[,
            attributs[i]], by = list(data[, productName]), FUN = "mean",
            na.rm = TRUE)[, 2]
    }
    colnames(CenteredProductMeansTable) = c("ProductCode", attributs)
    CenteredProductMeansTable[, attributs] = apply(CenteredProductMeansTable[,
        attributs], 2, "scale", center = TRUE, scale = scaleUnit)
    if (length(which(is.na(CenteredProductMeansTable)) > 0)) {
        stop("NA in matrix.")
    }
    if (nbAttributes < 2) {
        stop("Insufficient number of attributes.")
    }
    matBmod = CenteredProductMeansTable[, attributs]
    rownames(matBmod) = rownames(CenteredProductMeansTable[,
        "ProductCode"])
    resSvd = svd(matBmod)
    eigVal = (resSvd$d)^2
    eigVec = resSvd$v
    individuals = as.matrix(CenteredProductMeansTable[, attributs]) %*%
        eigVec
    variables = cor(as.matrix(CenteredProductMeansTable[, attributs]),
        individuals)
    rownames(individuals) = CenteredProductMeansTable[, "ProductCode"]
    p = 1/nbAttributes
    contrib = eigVal/sum(eigVal)
    NbDimSig = length(contrib[contrib > p])
    suppIndividuals = NULL
    if (subjectName %in% colnames(data)) {
        CenteredProductSubjectMeanTable = data.frame(NULL)
        CenteredProductSubjectMeanTable = aggregate(data[, attributs[1]],
            by = list(data[, productName], data[, subjectName]),
            FUN = "mean", na.rm = TRUE)
        for (i in 2:length(attributs)) {
            CenteredProductSubjectMeanTable[, attributs[i]] = aggregate(data[,
                attributs[i]], by = list(data[, productName],
                data[, subjectName]), FUN = "mean", na.rm = TRUE)[,
                3]
        }
        colnames(CenteredProductSubjectMeanTable) = c("ProductCode",
            "SubjectCode", attributs)
        CenteredProductSubjectMeanTable[, attributs] = apply(CenteredProductSubjectMeanTable[,
            attributs], 2, "scale", center = TRUE, scale = FALSE)
        suppIndividuals = CenteredProductSubjectMeanTable
        if (option == "Correlation") {
            CenteredProductMeansTableNotScaled = data.frame(NULL)
            CenteredProductMeansTableNotScaled = aggregate(data[,
                attributs[1]], by = list(data[, productName]),
                FUN = "mean", na.rm = TRUE)
            for (i in 2:length(attributs)) {
                CenteredProductMeansTableNotScaled[, attributs[i]] = aggregate(data[,
                  attributs[i]], by = list(data[, productName]),
                  FUN = "mean", na.rm = TRUE)[, 2]
            }
            colnames(CenteredProductMeansTableNotScaled) = c("ProductCode",
                attributs)
            CenteredProductMeansTableNotScaled[, attributs] = apply(CenteredProductMeansTableNotScaled[,
                attributs], 2, "scale", center = TRUE, scale = FALSE)
            sdByAttributes = apply(as.matrix(CenteredProductMeansTableNotScaled[,
                attributs]), 2, sd)
            suppIndividuals[, attributs] = sweep(CenteredProductSubjectMeanTable[,
                attributs], 2, sdByAttributes, "/")
        }
    }
    resultList = list(B = matBmod, IndSup = suppIndividuals,
        EigenVectors = eigVec, EigenValues = eigVal, IndivCoord = individuals,
        VarCoord = variables, NbdimSig = NbDimSig, Option = option)
    return(resultList)
}

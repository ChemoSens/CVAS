\name{gettingAppropriateData}
\alias{gettingAppropriateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gettingAppropriateData
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gettingAppropriateData(selectedData, asMatrix = FALSE, scoreName = "Score", productName = "ProductCode", subjectName = "SubjectCode", attributeName = "AttributeCode", replicateName = "Replicate", replaceMode = "crossmean", scaleUnit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selectedData}{
%%     ~~Describe \code{selectedData} here~~
}
  \item{asMatrix}{
%%     ~~Describe \code{asMatrix} here~~
}
  \item{scoreName}{
%%     ~~Describe \code{scoreName} here~~
}
  \item{productName}{
%%     ~~Describe \code{productName} here~~
}
  \item{subjectName}{
%%     ~~Describe \code{subjectName} here~~
}
  \item{attributeName}{
%%     ~~Describe \code{attributeName} here~~
}
  \item{replicateName}{
%%     ~~Describe \code{replicateName} here~~
}
  \item{replaceMode}{
%%     ~~Describe \code{replaceMode} here~~
}
  \item{scaleUnit}{
%%     ~~Describe \code{scaleUnit} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

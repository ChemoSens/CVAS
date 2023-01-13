\name{CalculateEllipses}
\alias{CalculateEllipses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CalculateEllipses
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalculateEllipses(suppIndividualTable, vep, axes = c(1, 2), confInt = 0.9, ellipsesType = "barycentric", productName = "ProductCode", subjectName = "SubjectCode")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{suppIndividualTable}{
%%     ~~Describe \code{suppIndividualTable} here~~
}
  \item{vep}{
%%     ~~Describe \code{vep} here~~
}
  \item{axes}{
%%     ~~Describe \code{axes} here~~
}
  \item{confInt}{
%%     ~~Describe \code{confInt} here~~
}
  \item{ellipsesType}{
%%     ~~Describe \code{ellipsesType} here~~
}
  \item{productName}{
%%     ~~Describe \code{productName} here~~
}
  \item{subjectName}{
%%     ~~Describe \code{subjectName} here~~
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
function (suppIndividualTable, vep, axes = c(1, 2), confInt = 0.9, 
    ellipsesType = "barycentric", productName = "ProductCode", 
    subjectName = "SubjectCode") 
{
    LoadPackage("ellipse")
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
        coord.indiv.sup = indiv.sup \%*\% vep
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
        quant = sqrt(qchisq(confInt, 2))
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

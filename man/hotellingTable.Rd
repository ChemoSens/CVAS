\name{hotellingTable}
\alias{hotellingTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{hotellingTable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hotellingTable(matCva, vep = NULL, axes = c(1, 2), colAttributes = NULL, productName = "ProductCode")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matCva}{
%%     ~~Describe \code{matCva} here~~
}
  \item{vep}{
%%     ~~Describe \code{vep} here~~
}
  \item{axes}{
%%     ~~Describe \code{axes} here~~
}
  \item{colAttributes}{
%%     ~~Describe \code{colAttributes} here~~
}
  \item{productName}{
%%     ~~Describe \code{productName} here~~
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
function (matCva, vep = NULL, axes = c(1, 2), colAttributes = NULL, 
    productName = "ProductCode") 
{
    prod = levels(as.factor(matCva[, productName]))
    nb.prod = length(prod)
    tab = matrix(NA, nb.prod, nb.prod)
    NbAxes = length(axes)
    if (is.null(colAttributes)) {
        colAttributes = 3:dim(matCva)[2]
    }
    for (i in 1:nb.prod) {
        for (j in 1:nb.prod) {
            if (j != i) {
                if (j < i) {
                  if (!is.null(vep)) {
                    prodi = as.matrix(matCva[matCva[, productName] == 
                      prod[i], colAttributes]) \%*\% vep[, axes]
                    prodj = as.matrix(matCva[matCva[, productName] == 
                      prod[j], colAttributes]) \%*\% vep[, axes]
                  }
                  if (is.null(vep)) {
                    prodi = matCva[matCva[, productName] == prod[i], 
                      colAttributes]
                    prodj = matCva[matCva[, productName] == prod[j], 
                      colAttributes]
                  }
                  L = HotellingTest(prodi, prodj)
                  tab[i, j] = L$p.value
                  tab[j, i] = L$p.value
                }
            }
            else {
                tab[i, j] = 1
            }
        }
    }
    colnames(tab) = prod
    rownames(tab) = prod
    return(tab)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

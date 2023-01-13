\name{CVApack-package}
\alias{CVApack-package}
\alias{CVApack}
\docType{package}
\title{
CVApack-package
}
\description{
Allows plotting Canonical Variate Analysis
}
\details{
\tabular{ll}{
Package: \tab CVA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-03\cr
License: 
}
To obtain CVA results, use the \code{\link{CVA}} function on your dataset. Then, use the \code{\link{PlotCVA}} function to plot the corresponding graph.
}
\author{
Caroline Peltier

Maintainer: <caroline.peltier@dijon.inra.fr>

}

\examples{
data(cheeses)
PlotCVA(CVA(cheeses))
}

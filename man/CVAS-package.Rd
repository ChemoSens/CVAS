\name{CVAS-package}
\alias{CVAS-package}
\alias{CVAS}
\docType{package}
\title{
CVAS
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
resCVA=CVA(cheeses,option="tw",productName="ProductCode",
subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
PlotCVA(resCVA)
}

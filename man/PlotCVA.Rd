\name{PlotCVA}
\alias{PlotCVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PlotCVA
}
\description{
Plots a CVA graph 
}
\usage{
PlotCVA(res.CVA, panellists = "none", alpha = 0.1, confInt = 0.9, ellipsesType = "barycentric", productColors = NULL, returnX = FALSE, returnY = FALSE, linkBetweenProducts = TRUE, suppIndividualsToPlot = FALSE, type = "R", fileName = "CVA", mainTitle = "CVA",ellipsesCalculation = "Chi")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.CVA}{
result of the function CVA
}
  \item{panellists}{
representation of the panellists on the map: "none" if the panellist are not represented on the map, "labels" if they are represented by their names, or "points" if they are represented by points
}

  \item{confInt}{
limit for the confidence ellipses. By default 0.9
}
  \item{ellipsesType}{
"barycentric" or "individual". Barycentric if an ellipse represents a zone of confidence where the mean point is (with a probability of confInt),
 "individual" if an ellipse represents a zone where the subject scores are (with a probability of confInt)
}
  \item{productColors}{
vector of colors for products
}
  \item{returnX}{
if TRUE, the X-axis is reversed. FALSE by default
}
  \item{returnY}{
if TRUE, the Y-axis is reversed. FALSE by default
}
  \item{linkBetweenProducts}{
if TRUE, two products which were equal with Hotelling T2 test are linked by a segment. TRUE by default
}

  \item{type}{
type of the output. "R" corresponds to a classical R-graph, whereas "png","jpg" and "wmf" create a file in the working directory.
}
  \item{fileName}{
name of the file containing the CVA in the working directory
}
  \item{mainTitle}{
title of the map
}
  \item{ellipseCalculation}{"Chi", "F" or "Sas"
}
}


\references{
%% ~put references to the literature/web site here ~
}
\author{
Caroline Peltier
}

\seealso{
\code{\link{CVA}}
}
\examples{
resCVA=CVA(cheeses,representation="biplot",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate")
PlotCVA(resCVA)
resCVA=CVA(cheeses,representation="twoMaps",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate")
PlotCVA(resCVA)
}

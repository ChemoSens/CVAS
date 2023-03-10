% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotCVA.R
\name{PlotCVA}
\alias{PlotCVA}
\title{PlotCVA}
\usage{
PlotCVA(
  res.CVA,
  panellists = "none",
  confInt = 0.9,
  ellipsesType = "barycentric",
  productColors = NULL,
  returnX = FALSE,
  returnY = FALSE,
  linkBetweenProducts = TRUE,
  suppIndividualsToPlot = FALSE,
  type = "R",
  fileName = "CVA",
  mainTitle = "CVA",
  ellipsesByRep = F,
  cex = 0.8,
  ellipsesCalculation = "Chi"
)
}
\arguments{
\item{res.CVA}{result of the function CVA}

\item{panellists}{representation of the panellists on the map: "none" if the panellist are not represented on the map, "labels" if they are represented by their names, or "points" if they are represented by points}

\item{confInt}{limit for the confidence ellipses. By default 0.9}

\item{ellipsesType}{"barycentric" or "individual". Barycentric if an ellipse represents a zone of confidence where the mean point is (with a probability of confInt),
"individual" if an ellipse represents a zone where the subject scores are (with a probability of confInt)}

\item{productColors}{vector of colors for products}

\item{returnX}{if TRUE, the X-axis is reversed. FALSE by default}

\item{returnY}{if TRUE, the Y-axis is reversed. FALSE by default}

\item{linkBetweenProducts}{if TRUE, two products which were equal with Hotelling T2 test are linked by a segment. TRUE by default}

\item{suppIndividualsToPlot}{FALSE by default (for internal use)}

\item{type}{type of the output. "R" corresponds to a classical R-graph, whereas "png","jpg" and "wmf" create a file in the working directory.}

\item{fileName}{name of the file containing the CVA in the working directory}

\item{mainTitle}{title of the map}

\item{ellipsesByRep}{F or T. False by default, TRUE is an option that should not be selected (for internal use only)}

\item{cex}{size of the text, default 0.8}

\item{ellipsesCalculation}{"Chi", "F" or "Sas"}
}
\description{
Plots a CVA graph
}
\examples{
{
data(cheeses)
 resCVA=CVA(cheeses,representation="biplot",productName="ProductCode",
 subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
 PlotCVA(resCVA)
 resCVA=CVA(cheeses,representation="twoMaps",productName="ProductCode",
 subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
 PlotCVA(resCVA)}
}
\seealso{
{
 \code{\link{CVA}}
}
}

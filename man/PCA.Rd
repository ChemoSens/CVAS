\name{PCA}
\alias{PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA
}
\description{
Calculates the PCA of a dataset
}
\usage{
PCA(data, option="Covariance",productName="product",subjectName="subject",replicateName="rep",sessionName="session")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The dataset should have the following columns: subject, product, rep, attribute 1, ... attribute n
}
  \item{option}{
"Covariance" or "Correlation". Covariance corresponds to the un-normalized PCA (thus, a not discriminant attribute will have a week weight in the analysis) whereas Correlation corresponds to the normalized PCA (All attributes have the same weight in the analysis)
}
\item{productName}{name of the column containing the products. By default, "product"}
\item{subjectName}{name of the column containing the subjects. By default, "subject"}
\item{replicateName}{name of the column containing the replicates. By default, "rep"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{ B }{matrix of product effect covariance}
  \item{ IndSup }{dataframe containing the centered scores obtained for each product*subject evaluation}
  \item{ EigenVectors }{eigen vectors from the diagonalization of B}
  \item{ EigenValues }{eigen values from the diagonalization of B}
  \item{ IndivCoord }{coordinates of the products (individuals) in the PCA map}
  \item{ VarCoord }{coordinates of the attributes (variables) in the PCA map}
  \item{ NbdimSig }{number of significant dimensions}
   \item{ Option }{choosen option}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
caroline.peltier@dijon.inra.fr
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PlotPCA}},\code{\link{CVA}}, \code{link{PlotCVA}}
}
\examples{
data(cheeses)
PCA(cheeses, productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

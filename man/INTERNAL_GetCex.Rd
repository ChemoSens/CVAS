\name{INTERNAL_GetCex}
\alias{INTERNAL_GetCex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL_GetCex
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
INTERNAL_GetCex(nbPoints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbPoints}{
%%     ~~Describe \code{nbPoints} here~~
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
function (nbPoints) 
{
    cex = 1
    if (nbPoints < 5) {
        cex = 1.3
    }
    else if (nbPoints < 10) {
        cex = 1.2
    }
    else if (nbPoints < 15) {
        cex = 1.1
    }
    else if (nbPoints > 100) {
        cex = 0.5
    }
    else if (nbPoints > 60) {
        cex = 0.6
    }
    else if (nbPoints > 50) {
        cex = 0.7
    }
    else if (nbPoints > 40) {
        cex = 0.8
    }
    else if (nbPoints > 20) {
        cex = 0.9
    }
    else {
        cex = 1
    }
    return(cex)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

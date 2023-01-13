\name{FriendlyPValue}
\alias{FriendlyPValue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FriendlyPValue
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FriendlyPValue(pvalues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
%%     ~~Describe \code{pvalues} here~~
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
function (pvalues) 
{
    res = NULL
    for (p.val in pvalues) {
        if (is.na(p.val)) {
            fp = ""
        }
        else if (p.val > 0.05) {
            fp = round(p.val, 3)
        }
        else if (p.val <= 0.1 && p.val > 0.05) {
            fp = paste(round(p.val, 3), ".", sep = "")
        }
        else if (p.val <= 0.05 && p.val > 0.01) {
            fp = paste(round(p.val, 3), "*", sep = "")
        }
        else if (p.val <= 0.01 && p.val > 0.001) {
            fp = paste(round(p.val, 3), "**", sep = "")
        }
        else {
            fp = paste("<0.001***", sep = "")
        }
        res = c(res, fp)
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

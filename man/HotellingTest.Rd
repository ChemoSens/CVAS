\name{HotellingTest}
\alias{HotellingTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HotellingTest
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HotellingTest(matx, maty)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matx}{
%%     ~~Describe \code{matx} here~~
}
  \item{maty}{
%%     ~~Describe \code{maty} here~~
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
function (matx, maty) 
{
    egalite = function(x, y) {
        bool = TRUE
        for (i in 1:length(x)[1]) {
            if (x[i] != y[i]) {
                bool = FALSE
            }
        }
        return(bool)
    }
    nb.suj = dim(matx)[1]
    nb.att = dim(matx)[2]
    pval = 1
    moyx = apply(matx, 2, mean)
    moyy = apply(maty, 2, mean)
    resx = matx
    resy = maty
    for (i in 1:nb.suj) {
        resx[i, ] = matx[i, ] - moyx
    }
    for (i in 1:nb.suj) {
        resy[i, ] = maty[i, ] - moyy
    }
    mat.resx = as.matrix(resx)
    mat.resy = as.matrix(resy)
    W1 = t(mat.resx) \%*\% mat.resx
    W2 = t(mat.resy) \%*\% mat.resy
    W = (W1 + W2)/(2 * nb.suj - 2)
    if (det(W) > 1e-08) {
        Winv = solve(W)
        T2 = ((nb.suj * nb.suj)/(2 * nb.suj)) * (t(moyx - moyy) \%*\% 
            Winv \%*\% (moyx - moyy))
        t2 = T2[1, 1]
        stat = t2 * (2 * nb.suj - nb.att - 1)/((2 * nb.suj - 
            2) * nb.att)
        p.value = pf(stat, df1 = nb.att, df2 = 2 * nb.suj - 1 - 
            nb.att, lower.tail = FALSE)
    }
    else {
        print("la matrice n'est pas inversible")
        if (egalite(moyx, moyy)) {
            p.value = 1
            stat = NA
        }
        else {
            p.value = 0
            stat = NA
        }
    }
    L = list(p.value = p.value, stat = stat)
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

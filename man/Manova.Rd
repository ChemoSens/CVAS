\name{Manova}
\alias{Manova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Manova
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Manova(P, R, df2 = NULL, test = "Pillai")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
%%     ~~Describe \code{P} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{df2}{
%%     ~~Describe \code{df2} here~~
}
  \item{test}{
%%     ~~Describe \code{test} here~~
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
function (P, R, df2 = NULL, test = "Pillai") 
{
    if (det(R) > 1e-16) {
        invertible = TRUE
        p = qr(P + R)$rank
        q = qr(P)$rank
        s = min(p, q)
        m = (abs(p - q) - 1)/2
        v = df2
        n = (v - p - 1)/2
        vap = eigen(solve(R) \%*\% P)$values
        if (test == "Pillai") {
            stat = sum(vap/(1 + vap))
            fapprox = ((2 * n + s + 1)/(2 * m + s + 1)) * stat/(s - 
                stat)
            pvalue = pf(Re(fapprox), df1 = (2 * m + s + 1) * 
                s, df2 = (2 * n + s + 1) * s, lower.tail = FALSE)
        }
        if (test == "Hotelling") {
            stat = sum(vap)
            b = (p + 2 * n) * (q + 2 * n)/(2 * (2 * n + 1) * 
                (n - 1))
            c = (2 + (p * q + 2)/(b - 1))/(2 * n)
            if (n > 0) {
                fapprox = (stat/c) * ((4 + (p * q + 2)/(b - 1))/(p * 
                  q))
                pvalue = pf(Re(fapprox), df1 = p * q, df2 = 4 + 
                  (p * q + 2)/(b - 1), lower.tail = FALSE)
            }
            if (n <= 0) {
                fapprox = (2 * (s * n + 1) * stat)/(s * s * (2 * 
                  m + s + 1))
                pvalue = pf(Re(fapprox), df1 = s * (2 * m + s + 
                  1), df2 = 2 * (s * n + 1), lower.tail = FALSE)
            }
        }
        if (test == "Wilks") {
            stat = prod(1/(1 + vap))
            r = v - (p - q + 1)/2
            u = (p * q - 2)/4
            if (p * p + q * q - 5 > 0) {
                t = sqrt((p * p * q * q - 4)/(p * p + q * q - 
                  5))
            }
            else {
                t = 1
            }
            fapprox = ((r * t - 2 * u)/(p * q)) * (1 - stat^(1/t))/stat^(1/t)
            pvalue = pf(Re(fapprox), df1 = q * p, df2 = r * t - 
                2 * u, lower.tail = FALSE)
        }
        if (test == "Roy") {
            stat = max(Re(vap))
            r = max(p, q)
            fapprox = ((v - r + q)/r) * stat
            pvalue = pf(Re(fapprox), df1 = r, df2 = v - r + q, 
                lower.tail = FALSE)
        }
    }
    else {
        fapprox = NA
        pvalue = NA
        stat = NA
        invertible = FALSE
    }
    L = list(stat, fapprox, pvalue, invertible)
    names(L) = c("stat", "f", "pvalue", "invertible")
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

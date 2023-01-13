\name{GenericPlot}
\alias{GenericPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GenericPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenericPlot(type = "R", filewidth = 7, fileheight = 7, fileName = "output", CALLFUN = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{filewidth}{
%%     ~~Describe \code{filewidth} here~~
}
  \item{fileheight}{
%%     ~~Describe \code{fileheight} here~~
}
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
}
  \item{CALLFUN}{
%%     ~~Describe \code{CALLFUN} here~~
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
function (type = "R", filewidth = 7, fileheight = 7, fileName = "output", 
    CALLFUN = NULL) 
{
    CALLFUN = c(CALLFUN)
    if (!is.character(type)) 
        stop("The output format (type) must be a string!")
    if (!is.character(fileName)) 
        stop("The output name (fileName) must be a string!")
    if (type == "pdf") {
        pdf(file = paste(fileName, ".pdf", sep = ""), width = filewidth, 
            height = fileheight)
    }
    else if (type == "png") {
        png(file = paste(fileName, ".png", sep = ""), width = 100 * 
            filewidth, height = 100 * fileheight)
    }
    else if (type == "jpg") {
        jpeg(file = paste(fileName, ".jpg", sep = ""), width = 100 * 
            filewidth, height = 100 * fileheight)
    }
    else if (type == "wmf") {
        win.metafile(file = paste(fileName, ".wmf", sep = ""), 
            width = filewidth, height = fileheight)
    }
    else {
        x11()
    }
    if (!is.null(CALLFUN)) {
        lapply(CALLFUN, FUN = eval)
    }
    if (type != "R") {
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

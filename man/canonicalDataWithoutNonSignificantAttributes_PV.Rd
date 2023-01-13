\name{canonicalDataWithoutNonSignificantAttributes_PV}
\alias{canonicalDataWithoutNonSignificantAttributes_PV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
canonicalDataWithoutNonSignificantAttributes_PV
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
canonicalDataWithoutNonSignificantAttributes_PV(canonicalData, scoreName = "Score", alpha = 0.1, productName = "ProductCode", subjectName = "SubjectCode", attributeName = "AttributeCode", replicateName = "Replicate")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{canonicalData}{
%%     ~~Describe \code{canonicalData} here~~
}
  \item{scoreName}{
%%     ~~Describe \code{scoreName} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{productName}{
%%     ~~Describe \code{productName} here~~
}
  \item{subjectName}{
%%     ~~Describe \code{subjectName} here~~
}
  \item{attributeName}{
%%     ~~Describe \code{attributeName} here~~
}
  \item{replicateName}{
%%     ~~Describe \code{replicateName} here~~
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
function (canonicalData, scoreName = "Score", alpha = 0.1, productName = "ProductCode", 
    subjectName = "SubjectCode", attributeName = "AttributeCode", 
    replicateName = "Replicate") 
{
    pvalToStar = function(pvalue) {
        if (is.na(pvalue)) {
            star = ""
        }
        else {
            star = ""
            if (pvalue < 0.1) {
                star = "."
            }
            if (pvalue < 0.05) {
                star = "*"
            }
            if (pvalue < 0.01) {
                star = "**"
            }
            if (pvalue < 0.001) {
                star = "***"
            }
        }
        return(star)
    }
    if (attributeName \%in\% colnames(canonicalData) & productName \%in\% 
        colnames(canonicalData) & subjectName \%in\% colnames(canonicalData)) {
        CanonicalData = canonicalData
        CanonicalData[, attributeName] = as.factor(CanonicalData[, 
            attributeName])
        CanonicalData[, productName] = as.factor(CanonicalData[, 
            productName])
        CanonicalData[, subjectName] = as.factor(CanonicalData[, 
            subjectName])
        CanonicalData[, replicateName] = as.factor(CanonicalData[, 
            replicateName])
        Attributes = levels(CanonicalData[, attributeName])
        Replicates = levels(CanonicalData[, replicateName])
    }
    else {
        stop(paste("data requires columns whose names are", productName, 
            ",", subjectName, ",", attributeName, "and", scoreName, 
            sep = " "))
    }
    res = NULL
    statTable = data.frame()
    for (attribute in Attributes) {
        tmp = CanonicalData[CanonicalData[, attributeName] == 
            attribute, ]
        if (var(tmp[, scoreName]) != 0) {
            if (length(Replicates) > 1) {
                resAnova = anova(lm(as.formula(paste(scoreName, 
                  "~", productName, "*", subjectName, sep = "")), 
                  data = tmp))
                statTable[attribute, "F"] = resAnova[1, 3]/resAnova[3, 
                  3]
                statTable[attribute, "pvalue"] = pf(statTable[attribute, 
                  "F"], df1 = resAnova[1, 1], df2 = resAnova[3, 
                  1], lower = FALSE)
            }
            else {
                resAnova = anova(lm(as.formula(paste(scoreName, 
                  "~", productName, "+", subjectName, sep = "")), 
                  data = tmp))
                statTable[attribute, "F"] = resAnova[1, 4]
                statTable[attribute, "pvalue"] = resAnova[1, 
                  5]
            }
            statTable[attribute, "Significance"] = pvalToStar(statTable[attribute, 
                "pvalue"])
            if (statTable[attribute, "pvalue"] < alpha) {
                res = rbind(res, tmp)
            }
        }
        else {
            statTable[attribute, "F"] = NA
            statTable[attribute, "pvalue"] = 0
            statTable[attribute, "Significance"] = ""
        }
    }
    res[, attributeName] = as.factor(as.character(res[, attributeName]))
    res[, productName] = as.factor(as.character(res[, productName]))
    res[, subjectName] = as.factor(as.character(res[, subjectName]))
    resFinal = list(res, statTable)
    names(resFinal) = c("data", "statTable")
    return(resFinal)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{MultiMAM}
\alias{MultiMAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MultiMAM
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MultiMAM(frame, modelType = "total", negativeCorrection = TRUE, correctOnlyIfSignificant = FALSE, limitOfSignificance = 0.05, plotReg = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frame}{
%%     ~~Describe \code{frame} here~~
}
  \item{modelType}{
%%     ~~Describe \code{modelType} here~~
}
  \item{negativeCorrection}{
%%     ~~Describe \code{negativeCorrection} here~~
}
  \item{correctOnlyIfSignificant}{
%%     ~~Describe \code{correctOnlyIfSignificant} here~~
}
  \item{limitOfSignificance}{
%%     ~~Describe \code{limitOfSignificance} here~~
}
  \item{plotReg}{
%%     ~~Describe \code{plotReg} here~~
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
function (frame, modelType = "total", negativeCorrection = TRUE, 
    correctOnlyIfSignificant = FALSE, limitOfSignificance = 0.05, 
    plotReg = FALSE) 
{
    LoadPackage("doBy")
    attnames = labels(frame)[[2]][-1:-3]
    natt = length(attnames)
    ass = factor(as.character(frame[, 1]))
    assnames = levels(ass)
    nass = length(assnames)
    prod = factor(frame[, 2])
    prodnames = levels(prod)
    nprod = length(prodnames)
    rep = factor(frame[, 3])
    repnames = levels(rep)
    nrep = length(repnames)
    nrow = dim(frame)[1]
    ncol = dim(frame)[2]
    attnames1 = attnames2 = attnames
    asslevels = levels(unique(ass))
    frame[, 1] = ass
    names(frame)[1] = "ass"
    frame[, 2] = prod
    names(frame)[2] = "prod"
    frame[, 3] = rep
    names(frame)[3] = "rep"
    averages = rep(NA, natt)
    names(averages) = colnames(frame[, 4:ncol])
    nnegatives = rep(NA, natt)
    names(nnegatives) = colnames(frame[, 4:ncol])
    matriceDisag = matrix(NA, nprod * nass, natt)
    colnames(matriceDisag) = colnames(frame[, 4:ncol])
    matriceInter = matrix(NA, nprod * nass, natt)
    colnames(matriceInter) = colnames(frame[, 4:ncol])
    matriceProd = matrix(NA, nprod, natt)
    colnames(matriceProd) = colnames(frame[, 4:ncol])
    matriceScaling = matrix(NA, nprod * nass, natt)
    colnames(matriceScaling) = colnames(frame[, 4:ncol])
    matriceError = matrix(NA, nprod * nass * nrep, natt)
    colnames(matriceError) = colnames(frame[, 4:ncol])
    matriceXS = matrix(NA, nprod * nass, natt)
    colnames(matriceXS) = colnames(frame[, 4:ncol])
    matriceXSu = matrix(NA, nass, natt)
    colnames(matriceXSu) = colnames(frame[, 4:ncol])
    matriceXSInt = matrix(NA, nprod * nass, natt)
    colnames(matriceXSInt) = colnames(frame[, 4:ncol])
    matriceBeta = matrix(NA, nass, natt)
    colnames(matriceBeta) = colnames(frame[, 4:ncol])
    rownames(matriceBeta) = asslevels
    BetaMulti = rep(NA, nass)
    names(BetaMulti) = asslevels
    pScaling = rep(NA, natt)
    significativity = rep(NA, nass)
    significativityMatrix = matrix(NA, nass, natt)
    colnames(significativityMatrix) = colnames(frame[, 4:ncol])
    rownames(significativityMatrix) = asslevels
    decompositionDataFrame = NULL
    if (!(modelType \%in\% c("total", "classic", "multivariate"))) {
        modelType = "total"
    }
    if (modelType == "total") {
        for (p in 4:ncol) {
            X = frame[, p]
            mu = mean(X)
            frame$X = X
            xam = orderBy(~ass, data = data.frame(summaryBy(X ~ 
                ass, data = frame)))
            xpm = orderBy(~prod, data = data.frame(summaryBy(X ~ 
                prod, data = frame)))
            xapm = orderBy(~ass + prod, data = data.frame(summaryBy(X ~ 
                ass + prod, data = frame)))
            frame = frame[, which(colnames(frame) != "X")]
            xapm$int = xapm$X.mean - rep(xam$X.mean, rep(nprod, 
                nass)) - rep(xpm$X.mean, nass) + mu
            xapm$xsuj = xapm$X.mean - rep(xam$X.mean, rep(nprod, 
                nass))
            xapm$xs = rep(xpm$X.mean, nass) - mu
            matriceXS[, p - 3] = xapm[, "xs"]
            matriceXSInt[, p - 3] = xapm[, "int"]
        }
        for (i in 1:nass) {
            linearRegression = lm(as.vector(matriceXSInt[xapm[, 
                "ass"] == asslevels[i], ]) ~ as.vector(matriceXS[xapm[, 
                "ass"] == asslevels[i], ]))
            coeff = linearRegression$coefficients
            summaryResults = summary(linearRegression)
            if (plotReg) {
                plot.new()
                xpts = as.vector(matriceXS[xapm[, "ass"] == asslevels[i], 
                  ])
                ypts = as.vector(matriceXSu[xapm[, "ass"] == 
                  asslevels[i], ])
                colorpts = as.vector(rep(rainbow(natt), each = nprod))
                plot(0, 0, xlim = c(min(xpts, ypts) - 0.05, max(xpts, 
                  ypts) + 0.05), ylim = c(-0.05 + min(ypts, xpts), 
                  0.05 + max(xpts, ypts)))
                points(xpts, ypts, col = colorpts, pch = 16)
                abline(a = coeff[1], b = coeff[2] + 1, col = "red")
                abline(a = 0, b = 1, col = "black")
                abline(v = 0)
                abline(h = 0)
            }
            BetaMulti[i] = coeff[2]
            significativity[i] = summaryResults[[4]][2, 4]
        }
        names(BetaMulti) = asslevels
    }
    for (p in 4:ncol) {
        Beta = rep(NA, nass)
        BetaToUse = rep(NA, nass)
        if (modelType != "total") {
            significativity = rep(NA, nass)
        }
        X = frame[, p]
        mu = mean(X)
        averages[p - 3] = mu
        frame$X = X
        xam = orderBy(~ass, data = data.frame(summaryBy(X ~ ass, 
            data = frame)))
        xpm = orderBy(~prod, data = data.frame(summaryBy(X ~ 
            prod, data = frame)))
        xapm = orderBy(~ass + prod, data = data.frame(summaryBy(X ~ 
            ass + prod, data = frame)))
        xaprm = orderBy(~ass + prod + rep, data = data.frame(summaryBy(X ~ 
            rep + ass + prod, data = frame)))
        frame = frame[, which(colnames(frame) != "X")]
        xapm$int = xapm$X.mean - rep(xam$X.mean, rep(nprod, nass)) - 
            rep(xpm$X.mean, nass) + mu
        xapm$prodEffect = rep(xpm$X.mean, nass) - mu
        xam$sujEffect = xam$X.mean - mu
        xapm$sujEffect = rep(xam$X.mean, rep(nprod, nass)) - 
            mu
        xapm$xs = rep(xpm$X.mean, nass)
        xaprm$mu = mu
        xaprm$centeredMean = xaprm$X.mean - mu
        xaprm$sujEffect = rep(xapm$sujEffect, each = nrep)
        xaprm$prodEffect = rep(xapm$prodEffect, each = nrep)
        xaprm$int = rep(xapm$int, each = nrep)
        xaprm$err = xaprm$centeredMean - xaprm$sujEffect - xaprm$prodEffect - 
            xaprm$int
        nnegatives[p - 3] = 0
        for (i in 1:nass) {
            if (modelType == "multivariate") {
                res = lm(int ~ xs, data = subset(xapm, ass == 
                  asslevels[i]))
                aovres = anova(res)[, 2]
                Beta[i] = res$coef[2]
                BetaToUse[i] = Beta[i]
                if (aovres[2] > 1e-16) {
                  significativity[i] = summary(res)[[4]][2, 4]
                }
            }
            if (modelType == "total") {
                Beta[i] = BetaMulti[i]
                BetaToUse[i] = BetaMulti[i]
            }
            if (modelType == "classic") {
                Beta[i] = 0
                BetaToUse[i] = 0
            }
            if (BetaToUse[i] < (-1)) {
                nnegatives[p - 3] = nnegatives[p - 3] + 1
                if (negativeCorrection) {
                  BetaToUse[i] = 0
                }
            }
            xapm[xapm[, "ass"] == asslevels[i], "Beta"] = BetaToUse[i]
        }
        xapm[, "disag"] = xapm[, "int"] - xapm[, "Beta"] * (xapm[, 
            "xs"] - mu)
        xapm$scaling = xapm[, "Beta"] * (xapm[, "xs"] - mu)
        xapm[, "sujmean"] = rep(xam$X.mean, each = nprod)
        if (correctOnlyIfSignificant & modelType != "classic") {
            SSScal = nrep * t(xapm[, "scaling"]) \%*\% xapm[, "scaling"]
            SSDisag = nrep * t(xapm[, "disag"]) \%*\% xapm[, "disag"]
            SSAss = nrep * nprod * t(xam[, "sujEffect"]) \%*\% 
                xam[, "sujEffect"]
            dfScal = (nass - 1)
            dfDisag = (nass - 1) * (nprod - 2)
            FScal = (SSScal/dfScal)/(SSDisag/dfDisag)
            pScaling[p - 3] = pf(FScal, dfScal, dfDisag, lower.tail = FALSE)
            xapm[, "disag"] = xapm[, "int"]
            xapm[, "scaling"] = 0
            if (pScaling[p - 3] < limitOfSignificance) {
                xapm[, "disag"] = xapm[, "int"] - xapm[, "Beta"] * 
                  (xapm[, "xs"] - mu)
                xapm$scaling = xapm[, "Beta"] * (xapm[, "xs"] - 
                  mu)
            }
        }
        xapmass = xapm$ass
        xapmprod = xapm$prod
        xaprm$Beta = rep(xapm$Beta, each = nrep)
        xaprm$Scaling = rep(xapm[, "scaling"], each = nrep)
        xaprm$Disag = rep(xapm[, "disag"], each = nrep)
        matriceScaling[, p - 3] = xapm[, "Beta"] * (xapm[, "xs"] - 
            mu)
        matriceDisag[, p - 3] = xapm[, "disag"]
        matriceBeta[, p - 3] = Beta
        significativityMatrix[, p - 3] = significativity
        matriceInter[, p - 3] = xapm[, "int"]
        matriceProd[, p - 3] = xpm[, "X.mean"]
        matriceXS[, p - 3] = xapm[, "xs"] - mu
        matriceXSu[, p - 3] = xam[, "sujEffect"]
        matriceError[, p - 3] = xaprm$err
        xaprm2 = xaprm
        xaprm2[, "Beta"] = xaprm2[, "Beta"] + 1
        xaprm2[, "Attribute"] = attnames[p - 3]
        decompositionDataFrame = rbind(decompositionDataFrame, 
            xaprm2)
    }
    rownames(matriceProd) = xpm[, 1]
    rownames(matriceXSu) = xam[, 1]
    rownames(matriceScaling) = paste(xapm[, 1], xapm[, 2], sep = "_")
    rownames(matriceDisag) = paste(xapm[, 1], xapm[, 2], sep = "_")
    rownames(matriceInter) = paste(xapm[, 1], xapm[, 2], sep = "_")
    rownames(matriceError) = paste(xaprm[, 1], xaprm[, 2], xaprm[, 
        3], sep = "_")
    if (modelType == "multivariate") {
        resultedBeta = matriceBeta + 1
    }
    if (modelType == "total") {
        resultedBeta = BetaMulti + 1
    }
    if (modelType != "total" & modelType != "multivariate") {
        resultedBeta = NULL
    }
    L = list(decompositionDataFrame, resultedBeta, averages, 
        matriceProd, matriceXSu, matriceScaling, matriceDisag, 
        matriceError, matriceInter, nnegatives, pScaling)
    names(L) = c("decomposition", "Beta", "avgMat", "prodMat", 
        "subjMat", "scalMat", "disagMat", "errMat", "intMat", 
        "nNeg", "pvalScal")
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

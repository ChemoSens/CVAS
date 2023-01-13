\name{PanelPerformances}
\alias{PanelPerformances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PanelPerformances
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PanelPerformances(frame, modelType = "total", negativeCorrection = TRUE, correctOnlyIfSignificant = FALSE, limitOfSignificance = 0.05, onlySignificantDim = FALSE, manovaTest = "Hotelling", panelistPerf = FALSE, correlationTest = "none")
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
  \item{onlySignificantDim}{
%%     ~~Describe \code{onlySignificantDim} here~~
}
  \item{manovaTest}{
%%     ~~Describe \code{manovaTest} here~~
}
  \item{panelistPerf}{
%%     ~~Describe \code{panelistPerf} here~~
}
  \item{correlationTest}{
%%     ~~Describe \code{correlationTest} here~~
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
    onlySignificantDim = FALSE, manovaTest = "Hotelling", panelistPerf = FALSE, 
    correlationTest = "none") 
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
    if (!(modelType \%in\% c("total", "classic", "multivariate"))) {
        modelType = "total"
    }
    ListResults = MultiMAM(frame = frame, modelType = modelType, 
        negativeCorrection = negativeCorrection, correctOnlyIfSignificant = correctOnlyIfSignificant, 
        limitOfSignificance = limitOfSignificance, plotReg = FALSE)
    if (modelType == "total") {
        ListResultsMAMmultivariate = MultiMAM(frame = frame, 
            modelType = "multivariate", negativeCorrection = negativeCorrection, 
            correctOnlyIfSignificant = correctOnlyIfSignificant, 
            limitOfSignificance = limitOfSignificance, plotReg = FALSE)
        UsualBeta = ListResultsMAMmultivariate$Beta
        CorrectedBeta = apply(UsualBeta, 2, "/", ListResults$Beta)
    }
    scalingCoefficient = ListResults$Beta
    matriceScaling = ListResults$scalMat
    matriceProd = ListResults$prodMat
    matriceDisag = ListResults$disagMat
    matriceError = ListResults$errMat
    matriceInter = ListResults$intMat
    scalPval = ListResults$pvalScal
    nNeg = ListResults$nNeg
    PanelProd = nrep * nass * t(scale(matriceProd, scale = FALSE)) \%*\% 
        scale(matriceProd, scale = FALSE)
    PanelScal = nrep * t(matriceScaling) \%*\% matriceScaling
    PanelDisag = nrep * t(matriceDisag) \%*\% matriceDisag
    PanelInter = nrep * t(matriceInter) \%*\% matriceInter
    PanelError = t(matriceError) \%*\% matriceError
    if (modelType == "total") {
        disagDF = (nprod - 2) * (nass - 1)
    }
    if (modelType == "classic") {
        disagDF = (nprod - 1) * (nass - 1)
    }
    if (modelType == "multivariate") {
        if (negativeCorrection) {
            if (correctOnlyIfSignificant) {
                scalingDisc = sum(scalPval < limitOfSignificance, 
                  na.rm = TRUE)
                disagDF = (scalingDisc/natt) * ((nprod - 2) * 
                  (nass - 1) + mean(nNeg[scalingDisc], na.rm = TRUE)) + 
                  (1 - scalingDisc/natt) * (nprod - 1) * (nass - 
                    1)
                scalingDF = (scalingDisc/natt) * ((nass - 1) - 
                  mean(nNeg[scalingDisc], na.rm = TRUE))
            }
            else {
                disagDF = (nprod - 2) * (nass - 1) + mean(nNeg)
                scalingDF = (nass - 1) - mean(nNeg)
            }
        }
        else {
            if (correctOnlyIfSignificant) {
                scalingDisc = sum(scalPval < limitOfSignificance, 
                  na.rm = TRUE)
                disagDF = (scalingDisc/natt) * ((nprod - 2) * 
                  (nass - 1)) + (1 - scalingDisc/natt) * (nprod - 
                  1) * (nass - 1)
                scalingDF = (scalingDisc/natt) * (nass - 1)
            }
            else {
                disagDF = (nprod - 2) * (nass - 1)
                scalingDF = (nass - 1)
            }
        }
    }
    if (!onlySignificantDim) {
        multiPanelDiscrimination = Manova(P = PanelProd, R = PanelDisag, 
            disagDF, test = manovaTest)
        multiPanelDisagreement = Manova(PanelDisag, PanelError, 
            nass * nprod * (nrep - 1), test = manovaTest)
        multiPanelScaling = Manova(PanelScal, PanelDisag, disagDF, 
            test = manovaTest)
    }
    if (modelType == "multivariate" || modelType == "total") {
        multiPanelPerformances = matrix(NA, 4, 3)
        rownames(multiPanelPerformances) = c("discrimination", 
            "scaling", "agreement", "repeatability")
    }
    if (modelType == "classic") {
        multiPanelPerformances = matrix(NA, 3, 3)
        rownames(multiPanelPerformances) = c("discrimination", 
            "agreement", "repeatability")
    }
    colnames(multiPanelPerformances) = c("stat", "f", "pvalue")
    multiPanelPerformances["discrimination", "f"] = round(Re(multiPanelDiscrimination$f), 
        digits = 3)
    multiPanelPerformances["discrimination", "stat"] = round(Re(multiPanelDiscrimination$stat), 
        digits = 3)
    multiPanelPerformances["discrimination", "pvalue"] = round(Re(multiPanelDiscrimination$pvalue), 
        digits = 4)
    multiPanelPerformances["agreement", "stat"] = round(Re(multiPanelDisagreement$stat), 
        digits = 3)
    multiPanelPerformances["agreement", "f"] = round(Re(multiPanelDisagreement$f), 
        digits = 3)
    multiPanelPerformances["agreement", "pvalue"] = round(Re(multiPanelDisagreement$pvalue), 
        digits = 4)
    if (modelType == "multivariate" || modelType == "total") {
        multiPanelPerformances["scaling", "stat"] = round(Re(multiPanelScaling$stat), 
            digits = 3)
        multiPanelPerformances["scaling", "f"] = round(Re(multiPanelScaling$f), 
            digits = 3)
        multiPanelPerformances["scaling", "pvalue"] = round(Re(multiPanelScaling$pvalue), 
            digits = 4)
    }
    avgPanel = rep(NA, natt)
    names(avgPanel) = attnames
    FdiscrPanel = rep(NA, natt)
    names(FdiscrPanel) = attnames
    PdiscrPanel = rep(NA, natt)
    names(PdiscrPanel) = attnames
    FscalPanel = rep(NA, natt)
    names(FscalPanel) = attnames
    PscalPanel = rep(NA, natt)
    names(PscalPanel) = attnames
    FdisagPanel = rep(NA, natt)
    names(FdisagPanel) = attnames
    PdisagPanel = rep(NA, natt)
    names(PdisagPanel) = attnames
    PerrPanel = rep(NA, natt)
    names(PerrPanel) = attnames
    levelPanel = rep(NA, natt)
    names(levelPanel) = attnames
    FdiscrPanelist = matrix(NA, natt, nass)
    colnames(FdiscrPanelist) = assnames
    rownames(FdiscrPanelist) = attnames
    PdiscrPanelist = matrix(NA, natt, nass)
    colnames(PdiscrPanelist) = assnames
    rownames(PdiscrPanelist) = attnames
    FscalPanelist = matrix(NA, natt, nass)
    colnames(FscalPanelist) = assnames
    rownames(FscalPanelist) = attnames
    PscalPanelist = matrix(NA, natt, nass)
    colnames(PscalPanelist) = assnames
    rownames(PscalPanelist) = attnames
    FdisagPanelist = matrix(NA, natt, nass)
    colnames(FdisagPanelist) = assnames
    rownames(FdisagPanelist) = attnames
    PdisagPanelist = matrix(NA, natt, nass)
    colnames(PdisagPanelist) = assnames
    rownames(PdisagPanelist) = attnames
    FerrorPanelist = matrix(NA, natt, nass)
    colnames(FerrorPanelist) = assnames
    rownames(FerrorPanelist) = attnames
    PerrorPanelist = matrix(NA, natt, nass)
    colnames(PerrorPanelist) = assnames
    rownames(PerrorPanelist) = attnames
    levelPanelist = matrix(NA, natt, nass)
    colnames(levelPanelist) = assnames
    rownames(levelPanelist) = attnames
    decomposition = ListResults$decomposition
    multiLevelPanel = mean(decomposition[, "X.mean"], na.rm = TRUE)
    listAnova = list()
    for (d in 1:natt) {
        decompositionByDescr = decomposition[decomposition[, 
            "Attribute"] == attnames[d], ]
        SSprod = t(decompositionByDescr[, "prodEffect"]) \%*\% 
            decompositionByDescr[, "prodEffect"]
        DFprod = nprod - 1
        SSsuj = t(decompositionByDescr[, "sujEffect"]) \%*\% decompositionByDescr[, 
            "sujEffect"]
        DFsuj = nass - 1
        SSscal = t(decompositionByDescr[, "Scaling"]) \%*\% decompositionByDescr[, 
            "Scaling"]
        if (negativeCorrection) {
            DFscal = nass - 1 - nNeg[d]
        }
        else {
            DFscal = nass - 1
        }
        SSdisag = t(decompositionByDescr[, "Disag"]) \%*\% decompositionByDescr[, 
            "Disag"]
        if (negativeCorrection) {
            DFdisag = (nprod - 2) * (nass - 1) + nNeg[d]
        }
        else {
            DFdisag = (nprod - 2) * (nass - 1)
        }
        SSerror = t(decompositionByDescr[, "err"]) \%*\% decompositionByDescr[, 
            "err"]
        DFerror = nprod * nass * (nrep - 1)
        avgPanel[d] = mean(decompositionByDescr[, "X.mean"], 
            na.rm = TRUE)
        FdiscrPanel[d] = (SSprod/DFprod)/(SSdisag/DFdisag)
        PdiscrPanel[d] = pf(FdiscrPanel[d], DFprod, DFdisag, 
            lower.tail = FALSE)
        FscalPanel[d] = (SSscal/DFscal)/(SSdisag/DFdisag)
        PscalPanel[d] = pf(FscalPanel[d], DFscal, DFdisag, lower.tail = FALSE)
        FdisagPanel[d] = (SSdisag/DFdisag)/(SSerror/DFerror)
        PdisagPanel[d] = pf(FdisagPanel[d], DFdisag, DFerror, 
            lower.tail = FALSE)
        PerrPanel[d] = sqrt(SSerror/DFerror)
        listAnova[[attnames[d]]] = matrix(NA, 5, 5)
        rownames(listAnova[[attnames[d]]]) = c("Product", "Subject", 
            "Scaling", "Disag", "Residuals")
        colnames(listAnova[[attnames[d]]]) = c("DF", "SS", "MS", 
            "F", "P-value")
        listAnova[[attnames[d]]]["Product", "DF"] = DFprod
        listAnova[[attnames[d]]]["Product", "SS"] = SSprod
        listAnova[[attnames[d]]]["Product", "MS"] = SSprod/DFprod
        listAnova[[attnames[d]]]["Product", "F"] = FdiscrPanel[d]
        listAnova[[attnames[d]]]["Product", "P-value"] = PdiscrPanel[d]
        listAnova[[attnames[d]]]["Subject", "DF"] = DFsuj
        listAnova[[attnames[d]]]["Subject", "SS"] = SSsuj
        listAnova[[attnames[d]]]["Subject", "MS"] = SSsuj/DFsuj
        listAnova[[attnames[d]]]["Subject", "F"] = (SSsuj/DFsuj)/(SSdisag/DFdisag)
        listAnova[[attnames[d]]]["Subject", "P-value"] = pf(listAnova[[d]]["Subject", 
            "F"], df1 = DFsuj, df2 = DFdisag, lower.tail = FALSE)
        listAnova[[attnames[d]]]["Scaling", "DF"] = DFscal
        listAnova[[attnames[d]]]["Scaling", "SS"] = SSscal
        listAnova[[attnames[d]]]["Scaling", "MS"] = SSscal/DFscal
        listAnova[[attnames[d]]]["Scaling", "F"] = FscalPanel[d]
        listAnova[[attnames[d]]]["Scaling", "P-value"] = PscalPanel[d]
        listAnova[[attnames[d]]]["Disag", "DF"] = DFdisag
        listAnova[[attnames[d]]]["Disag", "SS"] = SSdisag
        listAnova[[attnames[d]]]["Disag", "MS"] = SSdisag/DFdisag
        listAnova[[attnames[d]]]["Disag", "F"] = FdisagPanel[d]
        listAnova[[attnames[d]]]["Disag", "P-value"] = PdisagPanel[d]
        listAnova[[attnames[d]]]["Residuals", "DF"] = DFerror
        listAnova[[attnames[d]]]["Residuals", "SS"] = SSerror
        listAnova[[attnames[d]]]["Residuals", "MS"] = SSerror/DFerror
        if (modelType == "total") {
            decompositionM = ListResultsMAMmultivariate$decomposition
            decompositionByDescrM = decompositionM[decompositionM[, 
                "Attribute"] == attnames[d], ]
        }
        if (panelistPerf) {
            for (i in 1:nass) {
                decompositionByDescrAndSuj = decompositionByDescr[decompositionByDescr[, 
                  "ass"] == assnames[i], ]
                effetProdInd = (decompositionByDescrAndSuj[, 
                  "prodEffect"] + decompositionByDescrAndSuj[, 
                  "int"])
                SSprodi = t(effetProdInd) \%*\% effetProdInd
                DFprodi = nprod - 1
                levelPanelist[d, i] = mean(decompositionByDescrAndSuj[, 
                  "X.mean"]) - avgPanel[d]
                if (modelType != "classic") {
                  SSscali = t(decompositionByDescrAndSuj[, "Scaling"]) \%*\% 
                    decompositionByDescrAndSuj[, "Scaling"]
                  DFscali = 1
                }
                SSdisagi = t(decompositionByDescrAndSuj[, "Disag"]) \%*\% 
                  decompositionByDescrAndSuj[, "Disag"]
                DFdisagi = nprod - 2
                if (decompositionByDescrAndSuj[1, "Beta"] < 0) {
                  DFdisagi = nprod - 1
                }
                if (correctOnlyIfSignificant & (scalPval[d] > 
                  limitOfSignificance)) {
                  DFdisagi = nprod - 1
                }
                SSerrori = t(decompositionByDescrAndSuj[, "err"]) \%*\% 
                  decompositionByDescrAndSuj[, "err"]
                DFerrori = nprod * (nrep - 1)
                FdiscrPanelist[d, i] = (SSprodi/DFprodi)/(SSerrori/DFerrori)
                PdiscrPanelist[d, i] = pf(FdiscrPanelist[d, i], 
                  DFprodi, DFerrori, lower.tail = FALSE)
                if (modelType != "classic") {
                  FscalPanelist[d, i] = (SSscali/DFscali)/(SSdisagi/DFdisagi)
                  PscalPanelist[d, i] = pf(FscalPanelist[d, i], 
                    DFscali, DFdisagi, lower.tail = FALSE)
                }
                FerrorPanelist[d, i] = (SSerrori/DFerrori)
                FerrorTMP = (SSerrori/(nprod * (nrep - 1)))/((SSerror - 
                  SSerrori)/((nass - 1) * nprod * (nrep - 1)))
                PerrorPanelist[d, i] = pf(FerrorTMP, df1 = nprod * 
                  (nrep - 1), df2 = nprod * (nrep - 1) * (nass - 
                  1), lower.tail = FALSE)
                if (modelType == "total") {
                  decompositionByDescrAndSujM = decompositionByDescrM[decompositionByDescrM[, 
                    "ass"] == assnames[i], ]
                  if (modelType != "classic") {
                    SSscali = t(decompositionByDescrAndSujM[, 
                      "Scaling"]) \%*\% decompositionByDescrAndSujM[, 
                      "Scaling"]
                    DFscali = 1
                  }
                  SSdisagi = t(decompositionByDescrAndSujM[, 
                    "Disag"]) \%*\% decompositionByDescrAndSujM[, 
                    "Disag"]
                  DFdisagi = nprod - 2
                  if (modelType != "classic") {
                    FscalPanelist[d, i] = (SSscali/DFscali)/(SSdisagi/DFdisagi)
                    PscalPanelist[d, i] = pf(FscalPanelist[d, 
                      i], DFscali, DFdisagi, lower.tail = FALSE)
                  }
                }
                if (correlationTest == "none") {
                  FdisagPanelist[d, i] = (SSdisagi/DFdisagi)/(SSerrori/DFerrori)
                  PdisagPanelist[d, i] = pf(FdisagPanelist[d, 
                    i], DFdisagi, DFerrori, lower.tail = FALSE)
                }
                if (correlationTest != "none") {
                  if (!correlationTest \%in\% c("kendall", "pearson", 
                    "spearman")) {
                    stop("correlationTest should be 'none','kendall','spearman' or 'pearson'")
                  }
                  moyennePanelist = rep(NA, nprod)
                  moyenneRestPanel = rep(NA, nprod)
                  for (z in 1:nprod) {
                    moyennePanelist[z] = mean(decompositionByDescr[decompositionByDescr[, 
                      "ass"] == assnames[i] & decompositionByDescr[, 
                      "prod"] == prodnames[z], "X.mean"], na.rm = TRUE)
                    moyenneRestPanel[z] = mean(decompositionByDescr[decompositionByDescr[, 
                      "ass"] != assnames[i] & decompositionByDescr[, 
                      "prod"] == prodnames[z], "X.mean"], na.rm = TRUE)
                  }
                  corTest = cor.test(moyennePanelist, moyenneRestPanel, 
                    method = correlationTest, alternative = "greater")
                  FdisagPanelist[d, i] = corTest$statistic
                  PdisagPanelist[d, i] = corTest$p.value
                }
            }
        }
        listPanelistPerf = list(levelPanelist, FdiscrPanelist, 
            PdiscrPanelist, FscalPanelist, PscalPanelist, FdisagPanelist, 
            PdisagPanelist, FerrorPanelist, PerrorPanelist)
        names(listPanelistPerf) = c("avg", "Fdiscr", "Pdiscr", 
            "Fscal", "Pscal", "Fdisag", "Pdisag", "Ferror", "Perror")
    }
    listPanelPerf = list(avgPanel, FdiscrPanel, PdiscrPanel, 
        FscalPanel, PscalPanel, FdisagPanel, PdisagPanel, PerrPanel)
    names(listPanelPerf) = c("avg", "Fdiscr", "Pdiscr", "Fscal", 
        "Pscal", "Fdisag", "Pdisag", "Error")
    L = list(multiPanelPerformances, ListResults$decomposition, 
        PanelDisag, listPanelPerf, listAnova)
    names(L) = c("multiPanelPerformances", "decomposition", "matW", 
        "listPanelPerf", "listAnova")
    if (modelType == "total" || modelType == "multivariate") {
        L[["Beta"]] = t(scalingCoefficient)
    }
    if (panelistPerf) {
        L[["listPanelistPerf"]] = listPanelistPerf
    }
    if (modelType == "total") {
        L[["CorrectedBeta"]] = t(CorrectedBeta)
        L[["UsualBeta"]] = t(UsualBeta)
    }
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

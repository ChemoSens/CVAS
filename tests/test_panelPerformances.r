library(CVAS)
data(cheeses)
attributes=colnames(cheeses)[5:ncol(cheeses)]
cheeses2=cheeses[,c("SubjectCode","ProductCode","Replicate",attributes)]
resPerf = PanelPerformances(frame = cheeses2,
                            modelType = "classic", manovaTest = "Hotelling",
                            onlySignificantDim = FALSE)

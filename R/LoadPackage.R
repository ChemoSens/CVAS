LoadPackage <-
function (p) 
{
    tryCatch(find.package(p), error = function(e) {
        options(repos = "http://cran.at.r-project.org/")
        install.packages(p, dependencies = TRUE)
    }, finally = library(p, character.only = T))
    tryCatch(TS_LogPackage(p), error = function(e) {
    })
}

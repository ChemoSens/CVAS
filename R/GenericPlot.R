#'@importFrom grDevices pdf png jpeg win.metafile dev.off dev.new
GenericPlot <-
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
        png(filename = paste(fileName, ".png", sep = ""), width = 100 *
            filewidth, height = 100 * fileheight)
    }
    else if (type == "jpg") {
        jpeg(filename = paste(fileName, ".jpg", sep = ""), width = 100 *
            filewidth, height = 100 * fileheight)
    }
    else if (type == "wmf") {
        win.metafile(filename = paste(fileName, ".wmf", sep = ""),
            width = filewidth, height = fileheight)
    }
    else {
       dev.new()
    }
    if (!is.null(CALLFUN)) {
        res = lapply(CALLFUN, FUN = eval)
    }
    if (type != "R") {
        dev.off()
    }
    invisible(res)
}

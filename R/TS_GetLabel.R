TS_GetLabel <-
function (label, language = NULL) 
{
    if (is.null(language)) {
        if (exists("TimeSensLanguage")) {
            language = TimeSensLanguage
        }
        else {
            language = "en"
        }
    }
    res = ""
    if (exists("TimeSensLabels")) {
        res = TimeSensLabels[TimeSensLabels$label == label, language]
    }
    if (length(res) > 0) {
        return(res[1])
    }
    return(label)
}

FriendlyPValue <-
function (pvalues) 
{
    res = NULL
    for (p.val in pvalues) {
        if (is.na(p.val)) {
            fp = ""
        }
        else if (p.val > 0.05) {
            fp = round(p.val, 3)
        }
        else if (p.val <= 0.1 && p.val > 0.05) {
            fp = paste(round(p.val, 3), ".", sep = "")
        }
        else if (p.val <= 0.05 && p.val > 0.01) {
            fp = paste(round(p.val, 3), "*", sep = "")
        }
        else if (p.val <= 0.01 && p.val > 0.001) {
            fp = paste(round(p.val, 3), "**", sep = "")
        }
        else {
            fp = paste("<0.001***", sep = "")
        }
        res = c(res, fp)
    }
    return(res)
}

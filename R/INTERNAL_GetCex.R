INTERNAL_GetCex <-
function (nbPoints) 
{
    cex = 1
    if (nbPoints < 5) {
        cex = 1.3
    }
    else if (nbPoints < 10) {
        cex = 1.2
    }
    else if (nbPoints < 15) {
        cex = 1.1
    }
    else if (nbPoints > 100) {
        cex = 0.5
    }
    else if (nbPoints > 60) {
        cex = 0.6
    }
    else if (nbPoints > 50) {
        cex = 0.7
    }
    else if (nbPoints > 40) {
        cex = 0.8
    }
    else if (nbPoints > 20) {
        cex = 0.9
    }
    else {
        cex = 1
    }
    return(cex)
}

HotellingTest <-
function (matx, maty) 
{
    egalite = function(x, y) {
        bool = TRUE
        for (i in 1:length(x)[1]) {
            if (x[i] != y[i]) {
                bool = FALSE
            }
        }
        return(bool)
    }
    nb.suj = dim(matx)[1]
    nb.att = dim(matx)[2]
    pval = 1
    moyx = apply(matx, 2, mean)
    moyy = apply(maty, 2, mean)
    resx = matx
    resy = maty
    for (i in 1:nb.suj) {
        resx[i, ] = matx[i, ] - moyx
    }
    for (i in 1:nb.suj) {
        resy[i, ] = maty[i, ] - moyy
    }
    mat.resx = as.matrix(resx)
    mat.resy = as.matrix(resy)
    W1 = t(mat.resx) %*% mat.resx
    W2 = t(mat.resy) %*% mat.resy
    W = (W1 + W2)/(2 * nb.suj - 2)
    if (det(W) > 1e-08) {
        Winv = solve(W)
        T2 = ((nb.suj * nb.suj)/(2 * nb.suj)) * (t(moyx - moyy) %*% 
            Winv %*% (moyx - moyy))
        t2 = T2[1, 1]
        stat = t2 * (2 * nb.suj - nb.att - 1)/((2 * nb.suj - 
            2) * nb.att)
        p.value = pf(stat, df1 = nb.att, df2 = 2 * nb.suj - 1 - 
            nb.att, lower.tail = FALSE)
    }
    else {
        print("la matrice n'est pas inversible")
        if (egalite(moyx, moyy)) {
            p.value = 1
            stat = NA
        }
        else {
            p.value = 0
            stat = NA
        }
    }
    L = list(p.value = p.value, stat = stat)
    return(L)
}

hotellingTable <-
function (matCva, vep = NULL, axes = c(1, 2), colAttributes = NULL, 
    productName = "ProductCode") 
{
    prod = levels(as.factor(matCva[, productName]))
    nb.prod = length(prod)
    tab = matrix(NA, nb.prod, nb.prod)
    NbAxes = length(axes)
    if (is.null(colAttributes)) {
        colAttributes = 3:dim(matCva)[2]
    }
    for (i in 1:nb.prod) {
        for (j in 1:nb.prod) {
            if (j != i) {
                if (j < i) {
                  if (!is.null(vep)) {
                    prodi = as.matrix(matCva[matCva[, productName] == 
                      prod[i], colAttributes]) %*% vep[, axes]
                    prodj = as.matrix(matCva[matCva[, productName] == 
                      prod[j], colAttributes]) %*% vep[, axes]
                  }
                  if (is.null(vep)) {
                    prodi = matCva[matCva[, productName] == prod[i], 
                      colAttributes]
                    prodj = matCva[matCva[, productName] == prod[j], 
                      colAttributes]
                  }
                  L = HotellingTest(prodi, prodj)
                  tab[i, j] = L$p.value
                  tab[j, i] = L$p.value
                }
            }
            else {
                tab[i, j] = 1
            }
        }
    }
    colnames(tab) = prod
    rownames(tab) = prod
    return(tab)
}

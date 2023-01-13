Manova <-
function (P, R, df2 = NULL, test = "Hotelling") 
{
    p = qr(P + R)$rank
    q = qr(P)$rank
    s = min(p, q)
    m = (abs(p - q) - 1)/2
    v = df2
    n = (v - p - 1)/2
    invertible = TRUE
    invertible = tryCatch({
        vap = eigen(solve(R) %*% P)$values
        s = TRUE
    }, error = function(errors) {
        errors = geterrmessage()
        print(errors)
        s = FALSE
        return(s)
    })
    if (invertible) {
        if (test == "Pillai") {
            stat = sum(vap/(1 + vap))
            fapprox = ((2 * n + s + 1)/(2 * m + s + 1)) * stat/(s - 
                stat)
            pvalue = pf(Re(fapprox), df1 = (2 * m + s + 1) * 
                s, df2 = (2 * n + s + 1) * s, lower.tail = FALSE)
        }
        if (test == "Hotelling") {
            stat = sum(vap)
            b = (p + 2 * n) * (q + 2 * n)/(2 * (2 * n + 1) * 
                (n - 1))
            c = (2 + (p * q + 2)/(b - 1))/(2 * n)
            if (n > 0) {
                fapprox = (stat/c) * ((4 + (p * q + 2)/(b - 1))/(p * 
                  q))
                pvalue = pf(Re(fapprox), df1 = p * q, df2 = 4 + 
                  (p * q + 2)/(b - 1), lower.tail = FALSE)
            }
            if (n <= 0) {
                fapprox = (2 * (s * n + 1) * stat)/(s * s * (2 * 
                  m + s + 1))
                pvalue = pf(Re(fapprox), df1 = s * (2 * m + s + 
                  1), df2 = 2 * (s * n + 1), lower.tail = FALSE)
            }
        }
        if (test == "Wilks") {
            stat = prod(1/(1 + vap))
            r = v - (p - q + 1)/2
            u = (p * q - 2)/4
            if (p * p + q * q - 5 > 0) {
                t = sqrt((p * p * q * q - 4)/(p * p + q * q - 
                  5))
            }
            else {
                t = 1
            }
            fapprox = ((r * t - 2 * u)/(p * q)) * (1 - stat^(1/t))/stat^(1/t)
            pvalue = pf(Re(fapprox), df1 = q * p, df2 = r * t - 
                2 * u, lower.tail = FALSE)
        }
        if (test == "Roy") {
            stat = max(Re(vap))
            r = max(p, q)
            fapprox = ((v - r + q)/r) * stat
            pvalue = pf(Re(fapprox), df1 = r, df2 = v - r + q, 
                lower.tail = FALSE)
        }
    }
    else {
        fapprox = NA
        pvalue = NA
        stat = NA
        invertible = FALSE
    }
    L = list(stat, fapprox, pvalue, invertible)
    names(L) = c("stat", "f", "pvalue", "invertible")
    return(L)
}

sim.com <- function(x = "list of named quadrat observations", burn = 5, relative = FALSE){
    r.names <- names(x)
    for (i in 1:burn){r.names <- sample(r.names)}
    out <- do.call(rbind, x)
    out <- split(out, f = r.names)
    out <- lapply(out, function(x) apply(x, 2, sum))
    out <- do.call(rbind, out)
    if (relative){
        out <- apply(out, 2, function(x) x / max(x))
        out <- cbind(out, rep(min(out[out != 0]) / 1000, nrow(out)))
    }
    return(out)
}

h2.reml <- function(x = "lmerTest fitted model", digits = 5){
    var.comp <- as.data.frame(VarCorr(x))[, "vcov"]
    out <- c(var.comp[1] / sum(var.comp), 
             tail(coef(summary(x))[1, ], n = 1))
    names(out) <- c("H2", "pvalue")
    return(round(out, digits))
}

test.checks <- function(x = "fitted model object"){
    print(shapiro.test(residuals(x)))
    print(car::leveneTest(x))
}


ord <- function(x){
    diag(x) <- 0
    sum(sign(x))/2
}

zeroCol <- function(x, n = 10){
    for (i in 1:ncol(x)){
        if (sum(x[, i]) < n){x[, i] <- x[, i] * 0}else{}
    }
    x
}

freqNet <- function(x, zero.diag = FALSE){
    if (class(x) == "data.frame"){x <- as.matrix(x)}
    x <- sign(x)
    x <- t(x) %*% x
    if (zero.diag){diag(x) <- 0}
    x
}

coNets <- function(x, ci.p = 95, scale = FALSE, return.signs = FALSE){
    Z.tab <- c("80" = 1.282, "85" = 1.440, "90" = 1.645, "95" = 1.960, "99" = 2.576)
    Z <- Z.tab[as.character(floor(ci.p))]
    if (class(x) == "data.frame"){x <- as.matrix(x)}
    x <- sign(x)
    N <- apply(x, 2, sum)
    A <- nrow(x)
    P <- as.matrix(N / A)
    Pab <- P %*% t(P)
    diag(Pab) <- P
    A <- array(nrow(x), dim = dim(Pab))
    Vab <- A * Pab * (1 - Pab)
    ci.u <- A * Pab + Z * Vab^(1/2)
    ci.l <- A * Pab - Z * Vab^(1/2)
    ab <- t(x) %*% x    
    net <- ab
    net[ab <= ci.u & ab >= ci.l] <- 0
    if (scale){net <- net / nrow(x)}
    if (return.signs){
        net[ab > ci.u] <- net[ab > ci.u] * 1
        net[ab < ci.l] <- net[ab < ci.l] * -1
    }
    return(net)
}

netDist <- function(x, zero.na = TRUE, method = "euclidean"){
    out <- array(0, dim = rep(length(x), 2))
    if (!is.null(names(x))){
        rownames(out) <- colnames(out) <- names(x)
    }
    if (grepl("bray", tolower(method)) | tolower(method) == "bc"){
        for (i in 1:length(x)){
            for (j in 1:length(x)){
                y <- data.frame(x[[i]][lower.tri(x[[i]])], x[[j]][lower.tri(x[[j]])])
                if (all(y == 0)){y[y == 0] <- 1}
                out[i, j] <- ecodist::bcdist(t(y))
            }
        }
    }else{
        for (i in 1:length(x)){
            for (j in 1:length(x)){
                out[i, j] <- sum((x[[i]][lower.tri(x[[i]])] - 
                                      x[[j]][lower.tri(x[[j]])])^2)^(1/2)
            }
        }
    }
    if (zero.na){out[is.na(out)] <- 0}
    dist(out)
}

netMean <- function(x, zero.na = TRUE){
    x <- Reduce("+", x) / sum(unlist(x))
    if (zero.na){x[is.na(x)] <- 0}
    x
}

rmZeros <- function(x, zero.diag = TRUE){
    if (zero.diag){diag(x) <- 0}
    x[apply(x, 1, sum) != 0, apply(x, 2, sum) != 0]
}

meanNet <- function(x, zero.diag = TRUE){
    mu <- x[[1]]
    for (i in 2:length(x)){
        mu <- mu + x[[i]]
    }
    mu <- mu / length(x)
    if (zero.diag){diag(mu) <- 0}
    return(mu)
}

varNet <- function(x, zero.diag = TRUE){
    mu <- meanNet(x, zero.diag = zero.diag)
    v <- (x[[1]] - mu)^2
    for (i in 2:length(x)){
        v <- v + (x[[i]] - mu)^2
    }
    v <- v / length(x)
    if (zero.diag){diag(v) <- 0}
    return(v)
}

distNet <- function(x, zero.diag = TRUE, dist.obj = TRUE){
    d <- matrix(0, nrow = length(x), ncol = length(x))
    for (i in 1:nrow(d)){
        for (j in 1:ncol(d)){
            xi <- x[[i]]
            xj <- x[[j]]
            if (zero.diag){diag(xi) <- diag(xj) <- 0}
            d[i, j] <- sqrt(sum((xi - xj)^2))
        }
    }
    if (dist.obj){d <- as.dist(d)}
    return(d)
}

ch.plot <- function(x = 'ordination matrix', g = 'groupings', cex = 1, plot.legend = FALSE, loc = 'topleft', mu.pch = 19){
    mu <- apply(x, 2, function(x, g) tapply(x, g, mean), g = g) 
    se <- apply(x, 2, function(x, g) tapply(x, g, function(x)
            sd(x)/sqrt(length(x))), g = g) 
    mu <- na.omit(mu) 
    se <- na.omit(se)
                                        #error bars
    cl.xu <- mu[, 1] +  se[, 1]
    cl.xl <- mu[, 1] -  se[, 1]
    cl.yu <- mu[, 2] + se[, 2]
    cl.yl <-  mu[, 2] -  se[, 2]
    if (plot.legend){
                                        #coloring
        mu.col <- rainbow(length(unique(g)))[as.numeric(unique(g))]
        plot(mu, pch = mu.pch, cex = cex, xlim = c(min(cl.xl), max(cl.xu)), ylim = c(min(cl.yl), max(cl.yu)), col = mu.col)
        for (i in 1:nrow(mu)){
            lines(x = c(cl.xl[i], cl.xu[i]), y = c(mu[i, 2], mu[i, 2]))
            lines(x = c(mu[i, 1], mu[i, 1]), y = c(cl.yl[i], cl.yu[i]))
        }
        legend(loc, legend = rownames(se), cex = cex*0.5, pch = mu.pch, col = mu.col, border = 'grey')
    }else{
                                        #coloring
        mu.col <- 'black'
        plot(mu, pch = mu.pch, cex = cex, xlim = c(min(cl.xl), max(cl.xu)), ylim = c(min(cl.yl), max(cl.yu)), col = mu.col)
        for (i in 1:nrow(mu)){
            lines(x = c(cl.xl[i], cl.xu[i]), y = c(mu[i, 2], mu[i, 2]))
            lines(x = c(mu[i, 1], mu[i, 1]), y = c(cl.yl[i], cl.yu[i]))
        }
    }
    mu
}

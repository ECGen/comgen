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

freqNets <- function(x, zero.diag = FALSE){
    if (class(x) == "data.frame"){x <- as.matrix(x)}
    x <- sign(x)
    x <- t(x) %*% x
    if (zero.diag){diag(x) <- 0}
    x
}

coNets <- function(x, zero.diag = FALSE, prob = TRUE, zero.na = TRUE, ci.p = 95){
    Z.tab <- c("80" = 1.282, "85" = 1.440, "90" = 1.645, "95" = 1.960, "99" = 2.576)
    Z <- Z.tab[as.character(floor(ci.p))]
    if (class(x) == "data.frame"){x <- as.matrix(x)}
    x <- sign(x)
    N <- apply(x, 2, sum)
    A <- sum(N)
    P <- as.matrix(N / A)
    Pab <- P %*% t(P)
    A <- array(sum(N), dim = dim(Pab))
    Vab <- A * Pab * (1 - Pab)
    ci.u <- A * Pab + Z * Vab^(1/2)
    ci.l <- A * Pab - Z * Vab^(1/2)
    ab <- t(x) %*% x    
    net <- ab
    net[ab < ci.u & ab > ci.l] <- 0
    if (zero.diag){diag(net) <- 0}
    if (prob){net <- Pab * sign(net)}
    if (zero.na){net[is.na(net)] <- 0}
    net
}

netDist <- function(x, zero.na = TRUE){
    out <- array(0, dim = rep(length(x), 2))
    if (!is.null(names(x))){
        rownames(out) <- colnames(out) <- names(x)
    }
    for (i in 1:length(x)){
        for (j in 1:length(x)){
            out[i, j] <- sum((x[[i]] - x[[j]])^2)^(1/2)
        }
    }
    if (zero.na){out[is.na(out)] <- 0}
    dist(out)
}

netMean <- function(x){
    Reduce("+", x) / sum(unlist(x))
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

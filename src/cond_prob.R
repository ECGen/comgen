## Bayes P(A|B) = (P(B|A) * P(A)) / P(B)

## P(B) = B / N
## P(A) = A / N
## P(A,B) = AB / N 
## P(B|A) = P(A,B) / P(A)
## P(A|B) = P(A,B) / P(B)

cond_prob <- function(a, b){
    n <- length(a)
    p.a <- sum(a) / n
    p.b <- sum(b) / n
    p.ab <- sum(sign(a + b == 2)) / n
    p.a_b <- p.ab / p.b
    p.b_a <- p.ab / p.a
    return(c(p.a_b, p.b_a))
}

cond_net <- function(x){
    out <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    rownames(out) <- colnames(out) <- colnames(x)
    for (i in 1:ncol(x)){
        for (j in i:ncol(x)){
            cp <- cond_prob(sign(x[, i]), sign(x[, j]))
            ##This makes P(row | col)
            out[j, i] <- cp[1]
            out[i, j] <- cp[2]
        }
    }
    out[is.na(out)] <- 0
    return(out)
}


cond.nets <- lapply(onc.q, function(x) cond_net(x) * sign(coNets(x)))

cond.dcen.onc <- unlist(lapply(cond.nets, function(x) 
    sna::centralization(x, FUN = sna::degree, normalize = FALSE)))

cond.reml <- lme4::lmer(I(log(cond.dcen.onc + 0.000001)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
cond.reml.pval <- RLRsim::exactRLRT(cond.reml)
cond.reml.pval



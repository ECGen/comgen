### Signed network metrics
as_graph_signed <- function(x){
    g = graph_from_adjacency_matrix(x, mode = "directed", weighted = TRUE)
    sv  <- numeric()
    for (i in seq(1, nrow(x))){
        for (j in seq(1, ncol(x))){
            if (x[i, j] != 0){sv  <- append(sv, sign(x[i, j]))}
        }
    }
    E(g)$sign <- sv
    return(g)
}

freeman <- function(x){
    sum(max(x) - x) 
}

a  <-  array(c(0.1, -0.1, 0, 0.5, 0.2, 0.1, 0.3, 0, 0), dim = c(3,3))
g <- as_graph_signed(a)
freeman(signnet::degree_signed(g, mode = "in", type = "neg"))
freeman(signnet::degree_signed(g, mode = "out", type = "neg"))
freeman(signnet::degree_signed(g, mode = "in", type = "pos"))
freeman(signnet::degree_signed(g, mode = "out", type = "pos"))

###
cor.test(onc.com[,"Xg"], BR)
cor.test(onc.com[,"Xg"], L)
cor.test(onc.com[,"Xg"], Cen)
cor.test(onc.com[,"Xg"], AMI)
round(cor(cbind(onc.dat[, c("BR", "CT", "pH", "CN")], onc.com[, 1:5])), 3)[1:4,5:9]
round(cor(cbind(onc.com[, 1:5], onc.dat[, c("PC", "SR", "SD", "L", "Cen", "AMI")])), 3)[1:5,6:11]

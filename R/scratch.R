
cor.test(onc.com[,"Xg"], BR)
cor.test(onc.com[,"Xg"], L)
cor.test(onc.com[,"Xg"], Cen)
cor.test(onc.com[,"Xg"], AMI)
round(cor(cbind(onc.dat[, c("BR", "CT", "pH", "CN")], onc.com[, 1:5])), 3)[1:4,5:9]
round(cor(cbind(onc.com[, 1:5], onc.dat[, c("PC", "SR", "SD", "L", "Cen", "AMI")])), 3)[1:5,6:11]

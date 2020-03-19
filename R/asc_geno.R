## SCRATCH analysis of ascendency metrics response to genotype
loadd(onc.dat)
loadd(cn.onc)
cn.onc.net <- lapply(cn.onc, as.network)
asc.l <- lapply(cn.onc.net, enaAscendency)
asc.df <- do.call(rbind, asc.l)
asc.df[is.na(asc.df)] <- 0
onc.df <- cbind(onc.dat, asc.df)

asc.reml <- lme4::lmer(I(ASC^(1/4)) ~ (1 | geno),  
                       data = onc.df, REML = TRUE)
ami.reml <- lme4::lmer(I(AMI^(1/4)) ~ (1 | geno),  
                       data = onc.df, REML = TRUE)

RLRsim::exactRLRT(asc.reml)
RLRsim::exactRLRT(ami.reml)
summary(lm(ASC~Cen, data = onc.df))
summary(lm(AMI~Cen, data = onc.df))

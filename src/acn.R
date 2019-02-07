### script for loading data for the
### ACN study
### Pit data collected 2012

library(ComGenR)
source("loadPitdata.R")

#########################
### Genotype effect
#########################

## Total abundance
abund <- apply(com.acn, 1, sum)
reml.abund.acn <- lme4::lmer(I(abund^(1/2)) ~ (1 | geno), 
                           data = acn.dat,
                           REML = TRUE)
p.reml.abund.acn <- RLRsim::exactRLRT(reml.abund.acn)
p.reml.abund.acn

RLRsim::exactRLRT(
    lme4::lmer(I(abund[acn.dat[, "leaf.type"] == "live"]^(1/2)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "live", ],
               REML = TRUE)
)
RLRsim::exactRLRT(
    lme4::lmer(I(abund[acn.dat[, "leaf.type"] == "sen"]^(1/2)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "sen", ],
               REML = TRUE)
)

## Richness
rich <- apply(com.acn, 1, function(x) sum(sign(x)))
reml.rich.acn <- lme4::lmer(I(rich^(1/2)) ~ (1 | geno), 
                           data = acn.dat,
                           REML = TRUE)
p.reml.rich.acn <- RLRsim::exactRLRT(reml.rich.acn)
p.reml.rich.acn

RLRsim::exactRLRT(
    lme4::lmer(I(rich[acn.dat[, "leaf.type"] == "live"]^(1/2)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "live", ],
               REML = TRUE)
)
RLRsim::exactRLRT(
    lme4::lmer(I(rich[acn.dat[, "leaf.type"] == "sen"]^(1/2)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "sen", ],
               REML = TRUE)
)

## Community Similarity

rel.com.acn <- rel(com.acn)
ds <- rep(min(rel.com.acn[rel.com.acn != 0]) / 100, nrow(rel.com.acn))
rel.com.acn <- cbind(rel.com.acn, ds)
set.seed(1234)
vegan::adonis2(rel.com.acn ~ geno * leaf.type, strata = acn.dat[, "tree"], 
              data = acn.dat, perm = 10000, sqrt.dist = TRUE, mrank = TRUE)

## Network similarity
set.seed(1234)
vegan::adonis2(d.cn.acn ~ geno * leaf.type, strata = acn.dat[, "tree"], 
              data = acn.dat, perm = 10000, sqrt.dist = TRUE, mrank = TRUE)

                                        # senescent versus live
vegan::adonis2(netDist(cn.acn[acn.dat[, "leaf.type"]  == "live"]) ~ geno, 
               data = acn.dat[acn.dat[, "leaf.type"] == "live", ], 
               sqrt.dist = TRUE, mrank = TRUE)
vegan::adonis2(netDist(cn.acn[acn.dat[, "leaf.type"]  == "sen"]) ~ geno, 
               data = acn.dat[acn.dat[, "leaf.type"] == "sen", ], 
               sqrt.dist = TRUE, mrank = TRUE)

## Network metrics
                                        # number of links
reml.l.acn <- lme4::lmer(I(l.cn.acn^(1/1)) ~ (1 | geno), 
                           data = acn.dat,
                           REML = TRUE)
p.reml.l.acn <- RLRsim::exactRLRT(reml.l.acn)
p.reml.l.acn

RLRsim::exactRLRT(
    lme4::lmer(I(l.cn.acn[acn.dat[, "leaf.type"] == "live"]^(1/1)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "live", ],
               REML = TRUE)
)
RLRsim::exactRLRT(
    lme4::lmer(I(l.cn.acn[acn.dat[, "leaf.type"] == "sen"]^(1/1)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "sen", ],
               REML = TRUE)
)

                                        # centralization
reml.cen.acn <- lme4::lmer(I(cen.cn.acn^(1/2)) ~ (1 | geno), 
                           data = acn.dat,
                           REML = TRUE)
p.reml.cen.acn <- RLRsim::exactRLRT(reml.cen.acn)
p.reml.cen.acn

RLRsim::exactRLRT(
    lme4::lmer(I(cen.cn.acn[acn.dat[, "leaf.type"] == "live"]^(1/2)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "live", ],
               REML = TRUE)
)
RLRsim::exactRLRT(
    lme4::lmer(I(cen.cn.acn[acn.dat[, "leaf.type"] == "sen"]^(1/2)) ~ (1 | geno), 
               data = acn.dat[acn.dat[, "leaf.type"] == "sen", ],
               REML = TRUE)
)

#########################
## Plots
#########################
                                        # abundance and richness
mdc.plot(acn.dat[, "geno"], abund, 
         xlab = "Genotype", ylab = "Standardized Score", 
         ylim = c(-1.5, 1.5), pch = 1)
mdc.plot(acn.dat[, "geno"], rich, 
         xlab = "Genotype", ylab = "Standardized Score", 
         add = TRUE, pch = 19)

                                        # Links and centralization
mdc.plot(acn.dat[, "geno"], l.cn.acn, 
         xlab = "Genotype", ylab = "Standardized Score", 
         ylim = c(-1.5, 1.5), pch = 1)
mdc.plot(acn.dat[, "geno"], cen.cn.acn, 
         xlab = "Genotype", ylab = "Standardized Score", 
         add = TRUE, pch = 19)

mdc.plot(acn.dat[, "geno"], l.cn.acn, ylim = c(-1, 1.75), 
         xlab = "Tree Genotype", ylab = "Standardized Metric",
         ord = order(tapply(l.cn.acn, acn.dat[, "geno"], mean), decreasing = TRUE))
mdc.plot(acn.dat[, "geno"], cen.cn.acn, add = TRUE, pch = 1, 
         ord = order(tapply(l.cn.acn, acn.dat[, "geno"], mean), decreasing = TRUE))
legend("topright", legend = c("Links", "Centralization"), pch = c(19, 1), bty = "none")


                                        # Network Plots
net.col <- sign(netMean(cn.acn))
net.col[net.col == -1] <- 2
net.col[net.col == 1] <- 1
coord <- gplot(abs(netMean(cn.acn)), gmode = "digraph", 
               displaylabels = TRUE, 
               edge.lwd = abs(netMean(cn.acn)) * 100, 
               edge.col = net.col,
               vertex.col = "black", 
               vertex.cex = 0.5,
               arrowhead.cex = 0.5, 
               label.cex = 1, 
               main = "All Trees")
                                        # Network ordination
coord <- ch.plot(ord.cn.acn, g = acn.dat[, "geno"], cex = 3, mu.pch = 19, pt.col = "white")
text(coord, labels = rownames(coord))
plot(vec.com.acn, col = "grey")
plot(vec.nm.acn, col = "red")
                                        # 

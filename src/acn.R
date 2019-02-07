### script for loading data for the
### ACN study
### Pit data collected 2012

library(ComGenR)
source("loadPitdata.R")

#########################
### Genotype effect
#########################

## Senescence
### Can we even use senescent leaves given the sampling design?
n.sen <- unlist(lapply(tree.arth[acn.dat[, "leaf.type"] == "sen"], nrow))
n.live <- unlist(lapply(tree.arth[acn.dat[, "leaf.type"] == "live"], nrow))
sen.trees <- tree.arth[names(tree.arth) %in% names(n.sen)[n.sen >= 20]]
sen.dat <- data.frame(do.call(rbind, strsplit(names(sen.trees), split = " ")))
colnames(sen.dat) <- c("tree", "geno", "leaf.type")
com.sen <- do.call(rbind, lapply(sen.trees, function(x) apply(x, 2 ,sum)))
cn.sen <- lapply(sen.trees, coNets, ci.p = 95, cond = TRUE)
vegan::adonis2(netDist(cn.sen) ~ geno,  data = sen.dat, 
               sqrt.dist = TRUE, mrank = TRUE)


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
set.seed(1234)
vegan::adonis2(rel.com.acn[acn.dat[, "leaf.type"] == "live", ] ~ geno, 
               data = acn.dat[acn.dat[, "leaf.type"] == "live", ], 
               perm = 10000, sqrt.dist = TRUE, mrank = TRUE)


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

##########################################
####### Species/Functional Groups ########
##########################################
## PB = P. betae 
## pb.pred = parasitic fly
## pb.abort = 
## edge.fold = sawfly
## pinch =
## mid.miner = lep
## edge.miner = lep
## tip.miner = lep
## fish.eye = ??
## tier = lep
## thrips = Thysanura
## chew.holes = ??
## chomp = large herbivore
## scrape = radula
## chew.edge = lep?


#########################
## Plots
#########################

## Main Results
## Genotypes differ in network structure on living leaves
## Response to PB is linked to senescence, more PB higher prob(senscence)

                                        # Pb frequency
plot(table(pit[, "leaf.type"]))
table(pit[pit[, "leaf.type"] == "live", "pb"])
table(pit[pit[, "leaf.type"] == "sen", "pb"])


                                        # total abundance (live vs sen)
                                        # richness (live vs sen)
mdc.plot(acn.dat[, "leaf.type"], abund, ylim = c(-1.5, 1.5), 
         xlab = "Tree Genotype", ylab = "Value", std = TRUE, 
         ord = order(tapply(abund, acn.dat[, "leaf.type"], mean), decreasing = TRUE))
mdc.plot(acn.dat[, "leaf.type"], rich, add = TRUE, pch = 1, xjit = 0,
         ord = order(tapply(abund, acn.dat[, "leaf.type"], mean), decreasing = TRUE))
legend("topright", legend = c("Abundance", "Richness"), pch = c(19, 1), bty = "none")






                                        # abundance and richness
mdc.plot(acn.dat[, "geno"], abund, ylim = c(-1.5, 1.5), 
         xlab = "Tree Genotype", ylab = "Standardized Metric",
         ord = order(tapply(abund, acn.dat[, "geno"], mean), decreasing = TRUE))
mdc.plot(acn.dat[, "geno"], rich, add = TRUE, pch = 1, xjit = 0.01,
         ord = order(tapply(abund, acn.dat[, "geno"], mean), decreasing = TRUE))
legend("topright", legend = c("Abundance", "Richness"), pch = c(19, 1), bty = "none")

                                        # Links and centralization
mdc.plot(acn.dat[, "geno"], l.cn.acn, ylim = c(-1, 1.5), 
         xlab = "Tree Genotype", ylab = "Standardized Metric",
         ord = order(tapply(abund, acn.dat[, "geno"], mean), decreasing = TRUE))
mdc.plot(acn.dat[, "geno"], cen.cn.acn, add = TRUE, pch = 1, xjit = 0.01,
         ord = order(tapply(abund, acn.dat[, "geno"], mean), decreasing = TRUE))
legend("topright", legend = c("Links", "Centralization"), pch = c(19, 1), bty = "none")


                                        # Network Plots
                                       # Live versus Sen
par(mfrow = c(1, 2))
set.seed(1234)
net.col <- sign(netMean(cn.acn))
net.col[net.col == -1] <- "red"
net.col[net.col == "1"] <- "darkgrey"
coord <- gplot(abs(netMean(cn.acn[acn.dat[, "leaf.type"] == "live"])), 
               gmode = "digraph", 
               displaylabels = TRUE, 
               edge.lwd = (abs(netMean(cn.acn[acn.dat[, "leaf.type"] == "live"]))) * 10,
               edge.col = net.col,
               vertex.col = "black", 
               vertex.cex = 0.5,
               arrowhead.cex = 0.5, 
               label.cex = 0.75, 
               main = "Live")
net.col <- sign(netMean(cn.sen))
net.col[net.col == -1] <- "red"
net.col[net.col == "1"] <- "darkgrey"
gplot(abs(netMean(cn.sen)), 
      coord = coord, 
      gmode = "digraph", 
      displaylabels = TRUE, 
      edge.lwd = (abs(netMean(cn.sen))) * 10,
      edge.col = net.col,
      vertex.col = "black", 
      vertex.cex = 0.5,
      arrowhead.cex = 0.5, 
      label.cex = 0.75, 
      main = "Senescent")

                                        # By Genotype
par(mfrow = c(1, 2))
set.seed(1234)
net.col <- sign(netMean(cn.acn))
net.col[net.col == -1] <- 2
net.col[net.col == 1] <- "darkgrey"
pdf(file = "../results/acn_live_nets.pdf", width = 9, height = 9)
par(mfrow = c(3, 3))
coord <- gplot(abs(netMean(cn.acn[acn.dat[, "leaf.type"] == "live"])), 
               gmode = "digraph", 
               displaylabels = TRUE, 
               edge.lwd = (abs(netMean(cn.acn[acn.dat[, "leaf.type"] == "live"]))) * 10,
               edge.col = net.col,
               vertex.col = "black", 
               vertex.cex = 0.5,
               arrowhead.cex = 0.5, 
               label.cex = 0.75, 
               main = "Live")
for (i in unique(acn.dat[, "geno"])){
    gplot(abs(netMean(cn.acn[acn.dat[, "geno"] == i & 
                                 acn.dat[, "leaf.type"] == "live"])), 
          coord = coord,
          gmode = "digraph", 
          displaylabels = TRUE, 
          edge.lwd = (abs(netMean(cn.acn[acn.dat[, "geno"] == i & 
                                             acn.dat[, "leaf.type"] == "live"]))) * 10, 
          edge.col = net.col,
          vertex.col = "black", 
          vertex.cex = 0.5,
          arrowhead.cex = 0.5, 
          label.cex = 0.75, 
          main = i)
}
dev.off()

                                        # Network ordination
coord <- ch.plot(ord.cn.acn, g = acn.dat[, "geno"], 
                 cex = 3, mu.pch = 19, pt.col = "white")
text(coord, labels = rownames(coord))
plot(vec.com.acn, col = "grey")
plot(vec.nm.acn, col = "red")
                                        # 

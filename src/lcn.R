###LCN: ONC Garden Analyses
###MKLau
###06Sep2018
source('lcn_load_gardens.R')
source('lcn_load_wild.R')
## Is there a manuscript associated with this project?
manuscript.dir <- "../../lcn_manuscript"

### REML

### We know from Lamit's dissertation work that lichen communities are
### heritable, largely driven by bark roughness
### Do we find similar patterns?

## Create a list to generate a results table
h2.tab <- matrix("", 8, 4)
colnames(h2.tab) <- c("Response",  "H2", "R2", "p-value")

## Total cover ~ genotype
ptc.reml <- lme4::lmer(I(ptc.onc^(1/2)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ptc.reml.pval <- RLRsim::exactRLRT(ptc.reml)
ptc.reml.pval
fligner.test(onc.dat$ptc.onc^(1/2), onc.dat$geno)
shapiro.test(residuals(ptc.reml))
h2.tab[1, "p-value"] <- ptc.reml.pval$"p.value"
h2.tab[1, "H2"] <- H2(ptc.reml, g = onc.dat$geno)
h2.tab[1, "R2"] <- R2(ptc.reml)
h2.tab[1, "Response"] <- "Percent Lichen Cover"

## Species richness ~ genotype
spr.reml <- lme4::lmer(I(spr.onc^(1/1)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
spr.reml.pval <- RLRsim::exactRLRT(spr.reml)
spr.reml.pval
shapiro.test(residuals(spr.reml))
fligner.test(onc.dat$spr.onc^(1/2), onc.dat$geno)
h2.tab[2, "p-value"] <- spr.reml.pval$"p.value"
h2.tab[2, "H2"] <- H2(spr.reml, g = onc.dat$geno)
h2.tab[2, "R2"] <- R2(spr.reml)
h2.tab[2, "Response"] <- "Lichen Species Richness"

## Bark roughness REML
prb.reml <- lme4::lmer(I(onc.rough^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
prb.reml.pval <- RLRsim::exactRLRT(prb.reml)
prb.reml.pval
fligner.test(onc.dat$onc.rough^(1/2), onc.dat$geno)
shapiro.test(residuals(prb.reml))
h2.tab[3, "p-value"] <- prb.reml.pval$"p.value"
h2.tab[3, "H2"] <- H2(prb.reml, g = onc.dat$geno)
h2.tab[3, "R2"] <- R2(prb.reml)
h2.tab[3, "Response"] <- "Percent Rough Bark"

## Is species richness correlated with percent cover?
summary(lm(spr.onc ~ ptc.onc))

## Were these correlated with bark roughness?
ptc.prb.lm <- lm(I(ptc.onc^(1/2)) ~ I(onc.rough^(1/2)))
summary(ptc.prb.lm)
fligner.test(onc.dat$ptc.onc, onc.dat$onc.rough)
shapiro.test(residuals(ptc.prb.lm))

spr.prb.lm <- lm(I(spr.onc^(1)) ~ I(onc.rough^(1/2)))
summary(spr.prb.lm)
fligner.test(onc.dat$spr.onc^(1), onc.dat$onc.rough)
shapiro.test(residuals(spr.prb.lm))

## COM ~ genotype + Bark roughness + PTC + SPR
set.seed(2)
rcom.perm <- vegan::adonis2(onc.com.rel^(1/1) ~ geno + onc.rough + ptc.onc + spr.onc, data = onc.dat, perm = 10000, mrank = TRUE)
set.seed(2)
com.perm <- vegan::adonis2(onc.com^(1/1) ~ geno + onc.rough + ptc.onc + spr.onc, data = onc.dat, perm = 10000, mrank = TRUE)
rcom.perm
com.perm

## Network ~ genotype and other varibles
set.seed(12345)
nd.perm <- vegan::adonis2(cn.d.onc ~ geno + onc.rough + spr.onc + ptc.onc, 
                          data = onc.dat, perm = 10000, mrank = TRUE)
nd.perm

## Genotype and community
set.seed(2)
rcom.g.perm <- vegan::adonis2(onc.com.rel^(1/4) ~ geno, data = onc.dat, perm = 10000, mrank = TRUE)
set.seed(2)
com.g.perm <- vegan::adonis2(onc.com^(1/4) ~ onc.geno, data = onc.dat, perm = 10000, mrank = TRUE)
rcom.g.perm
com.g.perm

h2.tab[4, "p-value"] <- unlist(rcom.g.perm)["Pr(>F)1"]
h2.tab[4, "H2"] <- H2(rcom.g.perm, g = onc.dat$geno)
h2.tab[4, "R2"] <- R2(rcom.g.perm)
h2.tab[4, "Response"] <- "Lichen Community Composition"

## Is network similarity correlated with community composition?
ecodist::mantel(cn.d.onc ~ vegdist(onc.com.rel), mrank = TRUE)
spr.d <- dist(onc.dat$spr.onc)
ptc.d <- dist(onc.dat$ptc.onc)
prb.d <- dist(onc.dat$onc.rough)
### rough -> cover -> rich -> net
ecodist::mantel(cn.d.onc ~ vegdist(onc.com.rel) + spr.d + ptc.d + prb.d, mrank = TRUE)


## Partial Mantels using RFLP distance
ecodist::mantel(cn.mu.d.onc ~ rflp.d)
ecodist::mantel(onc.com.mu.d ~ rflp.d)
ecodist::mantel(cn.mu.d.onc ~ onc.com.mu.d)


## Was lichen network similarity determined by genotype?
cn.perm <- vegan::adonis2(cn.d.onc ~ geno + onc.rough + spr.onc + ptc.onc, 
                          data = onc.dat, permutations = 10000, mrank = TRUE)
cn.perm
h2.tab[5, "p-value"] <- as.matrix(cn.perm)[1, "Pr(>F)"]
h2.tab[5, "H2"] <- H2(cn.perm, g = onc.dat[, "geno"], perm =10000)
h2.tab[5, "R2"] <- R2(cn.perm)
h2.tab[5, "Response"] <- "Lichen Network"

                                        # db rda for network similarity
dbr.cn.geno <- vegan::dbrda(cn.d.onc ~ geno, data = onc.dat, distance = "bray")
anova(dbr.cn.geno, permutations = 5000)
H2(dbr.cn.geno)

## What aspects of networks explained the similiarity?
## L = number of edges, LD = link density, C = connectivity,
## dcen = degree centrality
link.reml <- lme4::lmer(I(log(L + 0.00001)) ~ (1 | geno), 
                          data = onc.dat, REML = TRUE)
link.reml.pval <- RLRsim::exactRLRT(link.reml)
link.reml.pval
fligner.test(onc.dat$L^(1/4), onc.dat$geno)
shapiro.test(residuals(link.reml))
h2.tab[6, "p-value"] <- link.reml.pval$"p.value"
h2.tab[6, "H2"] <- H2(link.reml, g = onc.dat$geno)
h2.tab[6, "R2"] <- R2(link.reml)
h2.tab[6, "Response"] <- "Number of Network Links"

                                        # network centrality
cen.reml <- lme4::lmer(I(log(Cen + 0.000001)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
cen.reml.pval <- RLRsim::exactRLRT(cen.reml)
cen.reml.pval
fligner.test(onc.dat$L^(1/1), onc.dat$geno)
shapiro.test(residuals(cen.reml))
h2.tab[7, "p-value"] <- cen.reml.pval$"p.value"
h2.tab[7, "H2"] <- H2(cen.reml, g = onc.dat$geno)
h2.tab[7, "R2"] <- R2(cen.reml)
h2.tab[7, "Response"] <- "Network Centrality"

                                        # network modularity
mod.reml <- lme4::lmer(I(onc.ns[, "mod.lik"]^(1/4)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
mod.reml.pval <- RLRsim::exactRLRT(mod.reml)
mod.reml.pval
fligner.test(onc.ns[, "mod.lik"]^(1/4), onc.dat$geno)
shapiro.test(residuals(mod.reml))
h2.tab[8, "p-value"] <- mod.reml.pval$"p.value"
h2.tab[8, "H2"] <- H2(mod.reml, g = onc.dat$geno)
h2.tab[8, "R2"] <- R2(mod.reml)
h2.tab[8, "Response"] <- "Network Modularity"

                                        # network stats in relation to other variables
L.aov <- aov(I(log(L + 0.001)) ~ onc.rough + ptc.onc + spr.onc, data = onc.dat)
summary(L.aov)
shapiro.test(residuals(L.aov))
cen.aov <- aov(I(log(Cen + 0.00001)) ~ onc.rough + ptc.onc + spr.onc, data = onc.dat)
summary(cen.aov)
shapiro.test(residuals(cen.aov))
mod.aov <- aov(I(onc.ns[, "mod.lik"]^(1/4)) ~ onc.rough + ptc.onc + spr.onc, data = onc.dat)
summary(mod.aov)
shapiro.test(residuals((mod.aov)))
## 
summary(lm(L ~ Cen + mod.lik, data = data.frame(onc.ns)))
summary(lm(Cen ~ mod.lik + L, data = data.frame(onc.ns)))
summary(lm(mod.lik ~ Cen + L, data = data.frame(onc.ns)))
cor(onc.ns[, c("L", "Cen", "mod.lik")])
                                        # are these metrics correlated with network similarity
L.d <- dist(onc.dat$L)
cen.d <- dist(onc.dat$Cen)
mod.d <- dist(cn.mod.onc)

adonis2(cn.d.onc ~ Cen + L + cn.mod.onc, data = onc.dat, mrank = TRUE)
adonis2(cn.d.onc ~ Cen + cn.mod.onc + L, data = onc.dat, mrank = TRUE)
adonis2(cn.d.onc ~ cn.mod.onc +  Cen + L, data = onc.dat, mrank = TRUE)
adonis2(cn.d.onc ~ L + Cen + onc.ns[, "mod.lik"], data = onc.dat, mrank = TRUE)

## So, are there patterns in the centrality of individual lichen species?
sppcen.test <- apply(cen.spp[, apply(cen.spp, 2, sum) >= 2], 2, function(x)
    RLRsim::exactRLRT(lme4::lmer(I(log(x + 0.00001)) ~ (1 | geno), 
                                 data = onc.dat, REML = TRUE)))
sppcen.test <- do.call(rbind, lapply(sppcen.test, function(x)
    c(x[["statistic"]], x[["p.value"]])))
p.adjust(sppcen.test[, 2], "fdr")

### Structural Modeling
## checking variance explained by ordinations
nms.com <- nmds(vegdist(onc.com.rel), 2, 3)
nms.cn <- nmds(cn.d.onc, 1, 2)
range(nms.com$r2)
range(nms.cn$r2)
ord.cn <- nmds.min(nms.cn, 2)
ord.com <- nmds.min(nms.com, 3)
ord1.cn.reml <- lme4::lmer(I(ord.cn[, 1]^(1/1)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ord2.cn.reml <- lme4::lmer(I(ord.cn[, 2]^(1/1)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ord1.cn.reml.pval <- RLRsim::exactRLRT(ord1.cn.reml)
ord2.cn.reml.pval <- RLRsim::exactRLRT(ord2.cn.reml)
ord1.cn.reml.pval
ord2.cn.reml.pval
fligner.test(ord.cn[, 1]^(1/1), onc.dat$geno)
fligner.test(ord.cn[, 2]^(1/1), onc.dat$geno)

ord1.com.reml <- lme4::lmer(I(ord.com[, 1]^(1/1)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ord2.com.reml <- lme4::lmer(I(ord.com[, 2]^(1/1)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ord3.com.reml <- lme4::lmer(I(ord.com[, 3]^(1/1)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ord1.com.reml.pval <- RLRsim::exactRLRT(ord1.com.reml)
ord2.com.reml.pval <- RLRsim::exactRLRT(ord2.com.reml)
ord3.com.reml.pval <- RLRsim::exactRLRT(ord3.com.reml)
ord1.com.reml.pval
ord2.com.reml.pval
ord3.com.reml.pval
fligner.test(ord.com[, 1]^(1/1), onc.dat$geno)
fligner.test(ord.com[, 2]^(1/1), onc.dat$geno)
fligner.test(ord.com[, 3]^(1/1), onc.dat$geno)
summary(lm(ord.cn[, 1] ~ spr.onc + ptc.onc, data = onc.dat))
summary(lm(ord.cn[, 2] ~ spr.onc + ptc.onc, data = onc.dat))
summary(lm(ord.com[, 1] ~ spr.onc + ptc.onc, data = onc.dat))
summary(lm(ord.com[, 2] ~ spr.onc + ptc.onc, data = onc.dat))
summary(lm(ord.com[, 3] ~ spr.onc + ptc.onc, data = onc.dat))

## Network stats correlated with other lichen community aspects
## adonis2(cn.d.onc ~ ptc.onc + onc.rough + spr.onc +  L + Cen, data = onc.dat, sqrt = FALSE) 
## adonis2(cn.d.onc ~ geno * L , data = onc.dat, sqrt = FALSE) 


## Lichen size distribution
## X. gallericulata thalli are about 0.22 +/- 0.003 cm^2 on average
## with an average median size of 0.12 +/- 0.001 cm^2
## and, size does not vary significantly with genotype.
xgs.reml <- lme4::lmer(I(mean.thallus) ~ (1 | geno), 
                       data = xgs.data[xgs.data$geno %in% names(which(table(xgs.data$geno) > 2)), ],
                       REML = TRUE)
xgs.median.reml <- lme4::lmer(median.thallus ~ (1 | geno), 
                       data = xgs.data[xgs.data$geno %in% names(which(table(xgs.data$geno) > 2)), ],
                       REML = TRUE)
RLRsim::exactRLRT(xgs.reml)
RLRsim::exactRLRT(xgs.median.reml)
fligner.test(xgs.data$mean.thallus, xgs.data$geno)
fligner.test(xgs.data$median.thallus, xgs.data$geno)
mean(xgs.data$mean.thallus)
sd(xgs.data$mean.thallus) / (length(xgs.data$mean.thallus) - 1)
mean(xgs.data$median.thallus)
sd(xgs.data$median.thallus) / (length(xgs.data$median.thallus) - 1)


### What is the structure of the bipartite networks?
                                        # test for modularity 
                                        # modularity p-values
p.mod <- c(wild = length(mods.wild.sweb[mods.wild.sweb <= mod.wild]) / length(mods.wild.sweb),
           onc = length(mods.onc.sweb[mods.onc.sweb <= mod.onc]) / length(mods.onc.sweb), 
           pit = length(mods.pit.sweb[mods.pit.sweb <= mod.pit]) / length(mods.pit.sweb))
                                        # ses modularity
ses.mod <- c(wild = (mod.wild - mean(mods.wild.sweb)) / sd(mods.wild.sweb),
             onc = (mod.onc - mean(mods.onc.sweb)) / sd(mods.onc.sweb),
             pit = (mod.pit - mean(mods.pit.sweb)) / sd(mods.pit.sweb))
# nest.pit <- bipartite::nestedness(pit.com.gm.rel)
## sna::gplot(sp.up, gmode = "graph", displaylabels = TRUE, lwd = sp.up)

## Tables
h2.tab[, "H2"] <- round(as.numeric(h2.tab[, "H2"]), digits = 5)
h2.tab[, "R2"] <- round(as.numeric(h2.tab[, "R2"]), digits = 5)
h2.tab[, "p-value"] <- round(as.numeric(h2.tab[, "p-value"]), digits = 5)
h2.xtab <- xtable::xtable(h2.tab, 
                           caption = "Genotypic effects of cottonwood trees on the associated lichen community.", 
                          label = "tab:h2_table")
print(h2.xtab,
      type = "latex",
      file = "../results/h2_table.tex",
      include.rownames = FALSE,
      include.colnames = TRUE
)

## Plots
### Network diagram for explanation figure.


## Lichen size distribution
pdf("../results/xg_size.pdf", height = 5, width = 5)
plot(density(xgs.data$median.thallus),
     xlab = "Median Lichen Thallus Area (cm^2)", 
     main = "")
abline(v = median(xgs.data$median.thallus, na.rm = TRUE), lty = 2)
dev.off()


### Figure 2
pdf("../results/cn_chplot_onc.pdf", height = 5, width = 12)
par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 2.1) / 2)
gplot(netMean(cn.mu.onc), gmode = "graph", 
      displaylabels = TRUE, 
      edge.lwd = netMean(cn.mu.onc) * 20, 
      vertex.col = "darkgrey")
legend("topleft", legend = "A", bty = "n", cex = 1.5)
chp.coord <- ch.plot(cn.nms.onc, onc.geno, cex = 1.5)
plot(nv.onc, col = "darkgrey")
legend("topleft", legend = "B", bty = "n", cex = 1.5)
dev.off()

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
pdf("../results/bp_net_onc.pdf")
bipartite::plotweb(pw.onc, method = "normal", 
                   text.rot = 45, 
                   col.low = col.pal[mods.onc$tree], 
                   col.high = col.pal[mods.onc$sp],
                   bor.col.low = col.pal[mods.onc$tree], 
                   bor.col.high = col.pal[mods.onc$sp],
                   col.interaction = "grey70",
                   bor.col.interaction = "grey70", 
                   labsize = 1.5)
dev.off()

pdf("../results/chp_com_onc.pdf")
ch.plot(nms.onc, onc.geno)
##plot(cv.onc, col = "grey30")
dev.off()

## legend("topleft", legend = "A")

pdf("../results/connect_geno.pdf", width = 12, height = 7)
g.order <- tapply(ns.onc[, "C"], onc.geno, mean)
g.order <- names(g.order)[order(g.order, decreasing = TRUE)]
onc.g <- factor(onc.geno, levels = g.order)
plot(ns.onc[, "C"] ~ onc.g, xlab = "Tree Genotype", ylab = "Lichen Network Connectance (C)")
dev.off()

### Wild Stands

### Compare Stands

### Which wild uintah trees are similar to garden trees?
cn.all <- cn.wild
for (i in 1:length(cn.wild)){
    cn.all[[i]] <- cn.wild[[i]][match(rownames(cn.onc[[1]]), rownames(cn.wild[[i]])), 
                                match(colnames(cn.onc[[1]]), colnames(cn.wild[[i]]))]
}
cn.all <- append(cn.all, cn.onc)

cn.d.all <- netDist(cn.all, method = "bc")

if (!exists("cn.nms.all")){cn.nms.all <- nmds.min(nmds(cn.d.all, 2, 2))}

labs <- c(rep("wild", length(cn.wild)), og)
coords <- ch.plot(cn.nms.all, labs, mu.pch = "")
points(coords, pch = 19, col = "white", cex = 2)
text(coords, labels = rownames(coords))

### Update lichen manuscript
### Send tables and figures to manuscript directory
tabs.figs <- c("../results/h2_table.tex", 
               "../results/connect_geno.pdf",
               "../results/chp_com_onc.pdf",
               "../results/bp_net_onc.pdf",
               "../results/cn_chplot_onc.pdf",
               "../results/xg_size.pdf")
if (exists("manuscript.dir")){
    sapply(tabs.figs, file.copy, to = manuscript.dir)
}

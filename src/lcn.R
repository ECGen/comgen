###LCN: ONC Garden Analyses
###MKLau
###06Sep2018
source('lcn_load_onc.R')

### How do the unipartite networks based on the bipartite network
### compare to the unipartite networks from genotypes?
sp.up <- t(onc.com[, -ncol(onc.com)]) %*% onc.com[, -ncol(onc.com)]
sp.up <- (sp.up )^(1/5)
sna::gplot(sp.up, gmode = "graph", displaylabels = TRUE, lwd = sp.up)

### What are the heritabilities for all REML analyses?

sp.up <- onc.com[, -"ds"]

### Data objects:
## onc.com = "community" occurrences summed across all cells for each tree
## onc.q = occurrence matrices separated out for each tree
## onc.geno = genotypes
## prb.onc = percent rough bark (averaged between the upper and lower)

### Data notes:
## Trees were removed from the analysis genotype RL6 and N1.31
## No physciods
## Lecanoras merged

### REML

### We know from Lamit's dissertation work that lichen communities are
### heritable, largely driven by bark roughness
### Do we find similar patterns?

## Create a list to generate a results table
h2.tab <- matrix("", 9, 4)
colnames(h2.tab) <- c("Response", "Predictor", "p-value", "H2")

## Total cover ~ genotype
ptc.reml <- lme4::lmer(I(ptc.onc^(1/2)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ptc.reml.pval <- RLRsim::exactRLRT(ptc.reml)
fligner.test(onc.dat$ptc.onc^(1/2), onc.dat$geno)
shapiro.test(residuals(ptc.reml))
h2.tab[1, "p-value"] <- ptc.reml.pval$"p.value"
h2.tab[1, "H2"] <- H2(ptc.reml)
h2.tab[1, "Response"] <- "Percent Lichen Cover"
h2.tab[1, "Predictor"] <- "Tree Genotype"

## Species richness ~ genotype
spr.reml <- lme4::lmer(I(spr.onc^(1/2)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
spr.reml.pval <- RLRsim::exactRLRT(spr.reml)
shapiro.test(residuals(spr.reml))
fligner.test(onc.dat$spr.onc^(1/2), onc.dat$geno)
h2.tab[2, "p-value"] <- spr.reml.pval$"p.value"
h2.tab[2, "H2"] <- H2(spr.reml)
h2.tab[2, "Response"] <- "Lichen Species Richness"
h2.tab[2, "Predictor"] <- "Tree Genotype"

## Bark roughness REML
prb.reml <- lme4::lmer(I(onc.rough^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
prb.reml.pval <- RLRsim::exactRLRT(prb.reml)
fligner.test(onc.dat$onc.rough^(1/2), onc.dat$geno)
shapiro.test(residuals(prb.reml))
h2.tab[3, "p-value"] <- prb.reml.pval$"p.value"
h2.tab[3, "H2"] <- H2(prb.reml)
h2.tab[3, "Response"] <- "Percent Rough Bark"
h2.tab[3, "Predictor"] <- "Tree Genotype"

## Is species richness correlated with percent cover?
summary(lm(spr.onc ~ ptc.onc))

## Were these correlated with bark roughness?
ptc.prb.lm <- lm(I(ptc.onc^(1/2)) ~ onc.rough)
summary(ptc.prb.lm)
fligner.test(onc.dat$ptc.onc^(1/2), onc.dat$onc.rough)
shapiro.test(residuals(ptc.prb.lm))

## Is network similarity correlated with community richness?
vegan::adonis(cn.d.onc ~ spr.onc)
vegan::adonis(cn.d.onc ~ onc.rough)
vegan::adonis(cn.d.onc ~ onc.geno + onc.rough)

## Is network similarity correlated with community composition?
ecodist::mantel(cn.d.onc ~ vegdist(onc.com.rel^(1/4)))
vegan::adonis2(cen.d ~ geno, data = onc.dat)

## So, are there patterns in the centrlity of individual lichen species?
cen.d <- vegdist(cbind(cen.spp, ds = rep(1, nrow(cen.spp))))
ecodist::mantel(cen.d ~ vegdist(onc.com.rel^(1/4)))

RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Xg"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Cs"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Ls"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Ch"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
shapiro.test(residuals(lme4::lmer(I(log(cen.spp[, "Ch"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE)))
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Xm"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
RLRsim::exactRLRT(lme4::lmer(I(cen.spp[, "Xm"]) ~ (1 | geno), data = onc.dat, REML = TRUE))
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Pm"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
## Pa has centrlity of zero 
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Pu"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))
RLRsim::exactRLRT(lme4::lmer(I(log(cen.spp[, "Rs"] + 0.00001)) ~ (1 | geno), data = onc.dat, REML = TRUE))


xg.reml <- lme4::lmer(I(onc.rough^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
prb.reml.pval <- RLRsim::exactRLRT(prb.reml)
fligner.test(onc.dat$onc.rough^(1/2), onc.dat$geno)
shapiro.test(residuals(prb.reml))
h2.tab[3, "p-value"] <- prb.reml.pval$"p.value"
h2.tab[3, "H2"] <- H2(prb.reml)
h2.tab[3, "Response"] <- "Percent Rough Bark"
h2.tab[3, "Predictor"] <- "Tree Genotype"


## Was lichen network similarity determined by genotype?
cn.perm <- vegan::adonis2(cn.d.onc ~ geno, data = onc.dat, permutations = 10000)
h2.tab[4, "p-value"] <- as.matrix(cn.perm)[1, "Pr(>F)"]
h2.tab[4, "H2"] <- H2(cn.perm, g = onc.dat[, "geno"], perm =10000)
h2.tab[4, "Response"] <- "Lichen Network"
h2.tab[4, "Predictor"] <- "Genotype"

                                        # db rda for network similarity
dbr.cn.geno <- vegan::dbrda(cn.d.onc ~ geno, data = onc.dat, distance = "bray")
anova(dbr.cn.geno, permutations = 5000)
H2(dbr.cn.geno)

## What aspects of networks explained the similiarity?
## L = number of edges, LD = link density, C = connectivity,
## dcen = degree centrality
pairs(cbind(R =onc.dat$spr.onc, onc.ns[,c("L", "C", "Cen")]))

link.reml <- lme4::lmer(I(L^(1/4)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
link.reml.pval <- RLRsim::exactRLRT(link.reml)
fligner.test(onc.dat$L^(1/2), onc.dat$geno)
shapiro.test(residuals(link.reml))
h2.tab[5, "p-value"] <- link.reml.pval$"p.value"
h2.tab[5, "H2"] <- H2(link.reml)
h2.tab[5, "Response"] <- "Number of Network Links"
h2.tab[5, "Predictor"] <- "Genotype"

                                        # network centrality
cen.reml <- lme4::lmer(I(log(Cen + 0.00001)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
cen.reml.pval <- RLRsim::exactRLRT(cen.reml)
fligner.test(onc.dat$L^(1/1), onc.dat$geno)
shapiro.test(residuals(cen.reml))
h2.tab[6, "p-value"] <- cen.reml.pval$"p.value"
h2.tab[6, "H2"] <- H2(cen.reml)
h2.tab[6, "Response"] <- "Network Centrality"
h2.tab[6, "Predictor"] <- "Genotype"

## Tables
h2.tab[, "H2"] <- round(as.numeric(h2.tab[, "H2"]), digits = 2)
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

## Update figures



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


###LCN: ONC Garden Analyses
###MKLau
###06Sep2018

source('lcn_load_onc.R')

### Data objects:
## oc = "community" occurrences summed across all cells for each tree
## oq = occurrence matrices separated out for each tree
## og = genotypes

## prb.onc = percent rough bark (averaged between the upper and lower)
## prb.wild = percent rough bark (averaged between the upper and lower)
## ws = wild ses values
## wc = wild community

### Data notes:
## No physciods
## Lecanoras merged

### REML

### We know from Lamit's dissertation work that lichen communities are
### heritable, largely driven by bark roughness
## Do we find similar patterns?

## Total cover ~ genotype
ptc.reml <- lmerTest::lmer(I(ptc.onc^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
shapiro.test(residuals(ptc.reml))
hist(residuals(ptc.reml))
summary(ptc.reml)
h2.reml(ptc.reml, 10)

## Species richness ~ genotype
spr.reml <- lmerTest::lmer(I(spr.onc^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
shapiro.test(residuals(spr.reml))
hist(residuals(spr.reml))
summary(spr.reml)
h2.reml(spr.reml, 10)

## Bark roughness REML
prb.reml <- lmerTest::lmer(I(onc.rough^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
shapiro.test(residuals(prb.reml))
hist(residuals(prb.reml))
summary(prb.reml)
h2.reml(prb.reml, 10)

## Is species richness correlated with percent cover?
summary(lm(ptc.onc ~ spr.onc))

## Were these correlated with bark roughness?
ptc.prb.lm <- lm(I(ptc.onc^(1/2)) ~ onc.rough)
shapiro.test(residuals(ptc.prb.lm))
hist(residuals(ptc.prb.lm))
summary(ptc.prb.lm)

## Was composition determined by genotype?
adonis(vegdist(onc.com.rel^(1/4)) ~ geno, perm = 5000, data = onc.dat, mrank = TRUE)

## Was lichen network similarity determined by genotype?
cn.perm <- adonis(cn.d.onc ~  onc.geno, mrank = FALSE, perm = 10000)

## Is network similarity correlated with community composition?
ecodist::mantel(cn.d.onc ~ vegdist(onc.com.rel^(1/4)))

## What aspects of networks were determined by genotype?
colnames(ns.onc)

ns.C.reml <- lmerTest::lmer(I(C^(1/2)) ~ (1 | geno), data = data.frame(onc.dat, C = ns.onc[, "C"]), REML = TRUE)
shapiro.test(residuals(ns.C.reml))
hist(residuals(ns.C.reml))
summary(ns.C.reml)
h2.reml(ns.C.reml, 10)

ns.L.reml <- lmerTest::lmer(I(L^(1/2)) ~ (1 | geno), data = data.frame(onc.dat, L = ns.onc[, "L"]), REML = TRUE)
shapiro.test(residuals(ns.L.reml))
hist(residuals(ns.L.reml))
summary(ns.L.reml)
h2.reml(ns.L.reml, 10)

ns.LD.reml <- lmerTest::lmer(I(LD^(1/2)) ~ (1 | geno), data = data.frame(onc.dat, LD = ns.onc[, "LD"]), REML = TRUE)
shapiro.test(residuals(ns.LD.reml))
hist(residuals(ns.LD.reml))
summary(ns.LD.reml)
h2.reml(ns.LD.reml, 10)


## n = number of nodes, L = number of edges, C = connectivity,
## LD = link density, ppr = pathway proliferation rate,
colnames(ns.onc)


hist(unlist(nul.mod.onc), xlim = c(0.025, 0.12))
abline(v = obs.mod.onc)
bp.mod.onc


## Tables
perm.tab <- xtable::xtable(cn.perm$aov.tab, 
                           caption = "Pseudo-F Table for the perMANOVA test of genotype effect on lichen network similarity.", digits = 5, label = "tab:cn_perm")
print(perm.tab,
      type = "latex",
      file = "../results/cn_perm_onc.tex",
      include.rownames = TRUE,
      include.colnames = TRUE
)


## Plots

### Network diagram for explanation figure.

pdf("../docs/lcn_araujo_method_net.pdf")
rg10 <- sna::rgraph(n = 10, m=1, tprob=0.25, mode="graph", diag=FALSE, replace=FALSE,
         tielist=NULL, return.as.edgelist=FALSE)
gplot(rg15, gmode = "graph", vertex.col = "black")
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
                   text.rot = 90, 
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




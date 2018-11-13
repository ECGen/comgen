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
results <- list()

## Total cover ~ genotype
ptc.reml <- lme4::lmer(I(ptc.onc^(1/2)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
ptc.reml
ptc.reml.pval <- RLRsim::exactRLRT(ptc.reml)
ptc.reml.pval
fligner.test(onc.dat$ptc.onc^(1/2), onc.dat$geno)
shapiro.test(residuals(ptc.reml))
hist(residuals(ptc.reml))
H2(ptc.reml)

## Species richness ~ genotype
spr.reml <- lme4::lmer(I(spr.onc^(1/2)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
spr.reml.pval <- RLRsim::exactRLRT(spr.reml)
spr.reml.pval
shapiro.test(residuals(spr.reml))
fligner.test(onc.dat$spr.onc^(1/2), onc.dat$geno)
hist(residuals(spr.reml))
H2c(spr.reml)

## Bark roughness REML
prb.reml <- lme4::lmer(I(onc.rough^(1/2)) ~ (1 | geno), data = onc.dat, REML = TRUE)
prb.reml.pval <- RLRsim::exactRLRT(prb.reml)
prb.reml.pval
fligner.test(onc.dat$onc.rough^(1/2), onc.dat$geno)
shapiro.test(residuals(prb.reml))
hist(residuals(prb.reml))
H2c(prb.reml)

## Is species richness correlated with percent cover?
summary(lm(spr.onc ~ ptc.onc))

## Were these correlated with bark roughness?
ptc.prb.lm <- lm(I(ptc.onc^(1/2)) ~ onc.rough)
summary(ptc.prb.lm)
fligner.test(onc.dat$ptc.onc^(1/2), onc.dat$onc.rough)
shapiro.test(residuals(ptc.prb.lm))
hist(residuals(ptc.prb.lm))

## Is network similarity correlated with community richness?
vegan::adonis(cn.d.onc ~ spr.onc)
vegan::adonis(cn.d.onc ~ onc.rough)
vegan::adonis(cn.d.onc ~ onc.geno + onc.rough)

## Is network similarity correlated with community composition?
ecodist::mantel(cn.d.onc ~ vegdist(onc.com.rel^(1/4)))

## Was lichen network similarity determined by genotype?
summary(dbr.cn.geno)
vegan::adonis2(cn.d.onc ~ onc.geno)
vegan::adonis2(cn.d.onc ~ 1 | onc.geno)
anova(dbr.cn.geno, permutations = 5000)
                                        # db rda for network similarity
dbr.cn.geno <- vegan::dbrda(cn.d.onc ~ geno, data = onc.dat, distance = "bray")
ord <- vegan:::scores.cca(dbr.cn.geno)$"sites"
H2c(dbr.cn.geno)

## What aspects of networks explained the similiarity?
## L = number of edges, LD = link density, C = connectivity,
## dcen = degree centrality
pairs(cbind(R =onc.dat$spr.onc, onc.ns[,c("L", "C", "Cen")]))

link.reml <- lme4::lmer(I(L^(1/2)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
link.reml.pval <- RLRsim::exactRLRT(link.reml)
link.reml.pval
fligner.test(onc.dat$L^(1/2), onc.dat$geno)
shapiro.test(residuals(link.reml))
hist(residuals(link.reml))


cen.reml <- lme4::lmer(I(log(Cen + 0.0001)) ~ (1 | geno), 
                       data = onc.dat, REML = TRUE)
cen.reml.pval <- RLRsim::exactRLRT(cen.reml)
cen.reml.pval
fligner.test(onc.dat$L^(1/1), onc.dat$geno)
shapiro.test(residuals(cen.reml))
hist(residuals(cen.reml))


summary(lm(Cen ~ spr.onc, data = onc.dat))
summary(lm(L ~ spr.onc, data = onc.dat))

summary(lm(Cen ~ ptc.onc, data = onc.dat))
summary(lm(L ~ ptc.onc, data = onc.dat))

plot(Cen ~ spr.onc, data = onc.dat)
plot(Cen ~ ptc.onc, data = onc.dat)


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




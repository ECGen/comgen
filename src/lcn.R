###LCN: ONC Garden Analyses
###MKLau
###06Sep2018
library(enaR)
library(ComGenR)

source('lcn_load_onc.R')
source('lcn_load_wild.R')

### Data objects:
## oc = co-occurrences summed across all cells for each tree
## oq = co-occurrence matrices separated out for each tree
## og = genotypes
## osgmu = ses genotype means
## osgse = ses genotype SE
## os = onc ses values
## oco = co-occurrence counts
## och = checker counts
## prb.onc = percent rough bark (averaged between the upper and lower)
## prb.wild = percent rough bark (averaged between the upper and lower)
## ws = wild ses values
## wc = wild community

### Data notes:
## No physciods
## Lecanoras merged

### We know from Lamit's dissertation work that lichen communities are
### heritable, largely driven by bark roughness

### First test for lichen species and community heritability

adonis(cn.d.onc ~  og, mrank = FALSE, perm = 10000)
adonis(cn.d.onc ~  prb.onc, mrank = FALSE, perm = 10000)
adonis(prb.onc.d ~ og, mrank = FALSE, perm = 10000)

if (!exists("cn.nms.onc")){cn.nms.onc <- nmds.min(nmds(cn.d.onc, 2, 2))}
onc.com <- onc.com[, colnames(onc.com) != "ds"]

### Networks are correlated with community structure, once relativized
rel.com <- apply(cbind(onc.com, 
                       ds = rep(1, nrow(onc.com))), 
                 2, function(x) x / max(x))
com.d.onc <- vegdist(rel.com)
ecodist::mantel(cn.d.onc ~ com.d.onc)

### Network statistics
ns.onc <- lapply(cn.onc, pack)
ns.onc <- lapply(ns.onc, enaStructure)
ns.onc <- do.call(rbind, lapply(ns.onc, function(x) x[[2]]))

nv.onc <- envfit(cn.nms.onc, data.frame(onc.com, 
                                        R = prb.onc, 
                                        C = ns.onc[, c("C")]))



## n = number of nodes, L = number of edges, C = connectivity,
## LD = link density, ppr = pathway proliferation rate,

### Who are the players? Describe the networks.

## text(chp.coord, labels = rownames(chp.coord), bg = "white")

if (!exists("nms.onc")){
    nms.onc <- nmds.min(nmds(
        vegdist(data.frame(onc.com, ds = rep(1, nrow(onc.com)))), 2, 2))
}
cv.onc <- envfit(nms.onc, data.frame(onc.com, 
                                        R = prb.onc, 
                                        C = ns.onc[, c("C")]))

par(mfrow = c(1, 2))
ch.plot(nms.onc, og)
plot(cv.onc, col = "darkgrey")
chp.coord <- ch.plot(cn.nms.onc, og, cex = 1.5)
plot(nv.onc, col = "darkgrey")

cn.array.onc <- array(NA, dim = c(length(cn.onc), dim(cn.onc[[1]])))
for (i in 1:length(cn.onc)){cn.array.onc[i, , ] <- cn.onc[[i]]}
n.cent.onc <- sna::centralization(cn.array.onc, FUN = "degree")
sna::centralization(cn.array.onc, g=3,degree, cmode="indegree")

plot(ns.onc[, "C"] ~ factor(og))

mod.och.onc <- metaComputeModules(och, N=1)
plotModuleWeb(mod.onc)


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




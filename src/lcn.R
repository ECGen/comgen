###LCN: ONC Garden Analyses
###MKLau
###06Sep2018

rm(list=ls())
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


### We know from Lamit's dissertation work that lichen communities are heritable, largely driven by bark roughness

### First test for lichen species and community heritability

oneway.test(os ~ og) # Welch ANOVA
oneway.test(os ~ og, var.equal = TRUE) # Welch ANOVA

plot(os ~ prb.onc)

osgmu.d <- dist(osgmu[match(colnames(as.matrix(rflp.d)), names(osgmu))])
os.d <- dist(os)
prb.onc.d <- dist(prb.onc)


ecodist::mantel(osgmu.d ~ rflp.d, nperm = 10000)
ecodist::mantel(cn.mu.d.onc ~ rflp.d, nperm = 10000)

ecodist::mantel(cn.d.onc ~ os.d, nperm = 10000)
ecodist::mantel(cn.d.onc ~ prb.onc.d, nperm = 10000)

ecodist::mantel(os.d ~ prb.onc.d, nperm = 10000)

adonis(os.d ~ prb.onc)

summary(aov(prb.onc ~ og))
summary(aov(os ~ og))

w.lm <- lm(wses ~ prb.wild, 
             data = data.frame(ws, prb.wild))
summary(w.lm)

summary(w.lm)

## The SES indicated generally increased co-occurrences
## linked with bark roughness. Given that the lichen
## are functioning in a largely 2D environment and given
## that the scale at which we measured the lichen was close
## to the size of individual thalli, the lichen are likely
## interacting if they occur in the same cell, either through
## mutualism or competition.
## Because of this we chose to use each co-occurrence in a cell
## as an interaction. We do not assign a "sign" to the
## interaction as we cannot distinguish between positive or
## negative interactions (e.g., a mutualism and a parasitism
## could look identical in terms of the spatial pattern).



## Tree Genotype Influences Ecological Network Structure

## We observed significant unipartite (one-mode) network structure
## \cite{Araujo2011} in the lichen species interaction networks that
## was similar between the experimental garden and the natural stand
## (Fig. 1a and 1b; Garden: z = -6.31, p = 0.0002; Natural: z = -3.15,
## p = 0.002). 

## The two networks displayed high multivariate structural similarity
## (Mantel R = 0.51, p = 0.029).



## Node level eigen-centrality \cite{DeAngelis1989}, a measure of
## species importance that integrates indirect connections, showed
## strong correlation between the two stands (Fig. 1c; r = 0.7, t =
## 2.6135, df = 7, p = 0.035). Centrality was also highly correlated
## with total abundance in both networks (Fig. 1d; Garden: r = 0.77, t
## = 3.2427, df = 7, p = 0.014; Natural:  r = 0.86, t = 4.43, df = 7,
## p = 0.003).



## Network Response to Tree Trait Variation

## Genotype was a significant predictor of interactions on individual
## trees (Fig. 2a; F = 3.4213, num df = 12.000, denom df = 14.668,
## p-value = 0.01426). 
oneway.test(os ~ og) # Welch ANOVA
oneway.test(os ~ og, var.equal = TRUE) # Welch ANOVA

## https://mindingthebrain.blogspot.com/2014/02/three-ways-to-get-parameter-specific-p.html

nlme::lme(os ~ og, data = data.frame(os = os, og = factor(og)))
mu <- rep(mean(os), length(os))

os.adj <- os + abs(min(os)) + 1

pois <- MASS::fitdistr(os.adj, "Poisson")
car::qqp(os.adj, "pois", pois$estimate)

car::leveneTest(residuals(lme.os) ~ factor(og))

## Genetic similarity was a significant predictor or SES 
## I.e. more genetically similar trees tended to have more similar interaction nets

lcn.R
lcn_load_onc.R
lcn_onc.R

### Standardized effect size distance for genotypes is predicted by RFLP and 
mantel(oms.d ~ rflp.d, nperm = 100000)


## Network similarity is predicted by genotype
g.nd <- as.matrix(g.nd)
g.nd <- g.nd[!grepl("NP", rownames(g.nd)), !grepl("NP", colnames(g.nd))]
g.nd <- as.dist(g.nd)
g.info <- g.info[g.info[, "Garden"] == "onc", ]

adonis(g.nd ~ Geno, data = g.info)


## Individual tree genotypes with similar levels of bark roughness had
## similar levels of lichen interactions (Fig. 2a; Mantel R = 0.08, p
## = 0.013), which was similar to the correlation observed between
## bark roughness and lichen interactions in the natural stands (Fig
## 2b: r = -0.53, p = 0.050).

## And, genetic similarity was correlated with interaction network structure: 
mantel(oms.d ~ rflp.d, nperm = 5000)

## Genetic Structure Generates Foret-Scale Network Structure

## Using a bipartite (two-mode) network approach in which
## genotype-species networks were modeled using the species maximum
## relativized values of each lichen species across all
## \textit{P. angustifolia} genotypes, we found significant modularity
## in the common garden stand (Fig. 3a; z = 9.64, p < 0.001). When
## using the same analyses on individual trees in the natural stand,
## we also found significant modularity (Fig. 3b; z = 7.42, p <
## 0.001). Furthermore, nestedness of both of these networks was
## significantly lower than expected under a null model (Garden: z =
## -2.30, p < 0.001; Natural Stand: z = -2.84, p < 0.001), most likely
## as a result of module formation.


## Figures

### Fig 1. Four panel figure
## A. Unipartite network Nature Center
## B. Unipartite network Wild Stand
## C. Eigen centrality Nature Center vs Natural Stand
## D. Relative log(Total Abundance) vs Eigen Centrality

library(igraph)
wnec <- centralization.evcent(graph.adjacency(wan,weighted=TRUE))$vector
plot(graph.adjacency(wan,weighted=TRUE),vertex.size=wnec*10)
names(wnec) <- rownames(wan)
round(wnec,2)
plot(wnec[match(names(onec),names(wnec))]~onec,pch=19,xlab='Eigen Centrality (ONC)',ylab='Eigen Centrality (Uintah)',xlim=c(0,1))
abline(lm(wnec[match(names(onec),names(wnec))]~onec),lty=1)
cor.test(wnec[match(names(onec),names(wnec))],onec,method='p')

###Eigen Centrality
library(igraph)
onec <- centralization.evcent(graph.adjacency(oan,weighted=TRUE))$vector
plot(graph.adjacency(oan,weighted=TRUE),vertex.size=onec*10)
names(onec) <- rownames(oan)

par(mfrow = c(1,2))
gplot(g.mn, displaylabels=TRUE, gmode = "graph", edge.lwd = (g.mn + 1)^3)
gplot(w.mn, displaylabels=TRUE, gmode = "graph", edge.lwd = (w.mn + 1)^3)

## Fig 2. Two panel figure
## A. Mean percent Rough Bark vs Mean SES (Common Garden)
## B. Percent rough bark vs SES (Wild Stand)
ch.x <- cbind('Percent Rough Bark' = prb.onc,'SES' = os)
ch.plot(ch.x,og)
abline(lm(osgmu ~ prb.mu),lty=2)

## Fig 3. Two panel of bipartite networks
## A. Common Garden
## B. Natural Stand



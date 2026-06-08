# Analysis Notes

Make sure to run make.R first and load the necessary objects.

    library(drake)
    loadd(onc.dat)
    loadd(cn.d.onc)

## Is the network similarity driven by abundance of lichen?

We found support for the hypothesis that genotypically based variation
in lichen network structure is potentially driven by variation in bark
roughness. In an initital PerMANOVA with sequential partitioning of
variance explained, genotype has a strong, statistically significant
effect (R2 = 0.37, p-value = 0.040) when including tree traits and
lichen community indices in the model. In this ordering, condensed
tannins (R2 = 0.09, p-value = 0.486) and species richness (R2 = 0.11,
p-value = 0.017) are both significant, although with low explanatory
power. After moving genotype to the final predictor in the model,
genotype is no longer significant and has low explanatory power, while
both bark roughness (R2 = 0.14, p-value = 0.006) and species richness
(R2 = 0.20, 0.002) are now significant and have increased in the
variance in network similarity explained. Taken together, these results
suggest that there is potentially a causal relationship between
genotypically explained variance in bark roughness and lichen community
species richness that affects the similarity of lichen network
structure.

    cn.geno.trait.pc.sr <- vegan::adonis2(
                                   cn.d.onc~ geno + CT + pH + CN + BR + PC + SR,
                                   by = "term", 
                                   data = onc.dat, 
                                   mrank = TRUE,
                                   permutations = 100000)

    cn.trait.pc.sr.geno  <- vegan::adonis2(
                                       cn.d.onc~ CT + pH + CN + BR + PC + SR + geno,
                                       by = "term", 
                                       data = onc.dat, 
                                       mrank = TRUE,
                                       permutations = 100000)

# How does bark roughness affect the structure of the lichen networks? Is

it driven by lichen abudance (Percent Cover)?

Examination of the variaince in lichen network similarity explained by
species richness and the network metrics of size and centrality,
supported the hypothesis that network structure is both a function of
increasing numbers of species (i.e., nodes in the network) and other
structural differences. In analyzing two PerMANOVAs with these factors,
we found that network size (R2 = 0.77, p-value &lt; 0.001) and
centrality (R2 = 0.04, p-value &lt; 0.001) both explained vaiation in
network similarity regardless of whether species richness was the first
or last factor entered into the model. Additionally, although variation
in network similarity explained by centrality was low, it was
significant even after partitioning the variance explained by network
size. Based on this, we can conclude that network similarity is both of
function of the number of lichen species in the community determining
network size and other aspects of network strcutre, such as
centralization.

    cn.l.cen.sr  <- vegan::adonis2(
                               cn.d.onc~ L + Cen + SR,
                               by = "term", 
                               data = onc.dat, 
                               mrank = TRUE,
                               permutations = 100000)

    cn.sr.l.cen  <- vegan::adonis2(
                               cn.d.onc~ SR + L + Cen,
                               by = "term", 
                               data = onc.dat, 
                               mrank = TRUE,
                               permutations = 100000)

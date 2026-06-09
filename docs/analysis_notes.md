# Analysis Notes

    ## Make sure to run make.R first and load the necessary objects.

    library(pacman)
    p_load(drake, xtable, rmarkdown)
    loadd(onc.dat)
    loadd(cn.d.onc)

## Is the network similarity driven by abundance of lichen?

We found support for the hypothesis that genotypically based variation
in lichen network structure is potentially driven by variation in bark
roughness, not just as a function of increasing abundance or richness
of lichen species. In an initital PerMANOVA with sequential
partitioning of variance explained, genotype has a strong,
statistically significant effect (R2 = 0.37, p-value = 0.040) when
including tree traits and lichen community indices in the model. In
this ordering, condensed tannins (R2 = 0.09, p-value = 0.486) and
species richness (R2 = 0.11, p-value = 0.017) are both significant,
although with low explanatory power. After moving genotype to the
final predictor in the model, genotype is no longer significant and
has low explanatory power, while both bark roughness (R2 = 0.14,
p-value = 0.006) and species richness (R2 = 0.20, 0.002) are now
significant and have increased in the variance in network similarity
explained. Taken together, these results suggest that there is
potentially a causal relationship between genotypically explained
variance in bark roughness and lichen community species richness that
affects the similarity of lichen network structure.

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

    as.data.frame(cn.geno.trait.pc.sr)

    ##          Df  SumOfSqs         R2         F     Pr(>F)
    ## geno      9 257.29223 0.37103551 2.3523174 0.04005960
    ## CT        1  60.19753 0.08680955 4.9532522 0.04815952
    ## pH        1   7.81880 0.01127532 0.6433567 0.44353556
    ## CN        1  13.75107 0.01983011 1.1314832 0.29715703
    ## BR        1  14.71778 0.02122420 1.2110281 0.28096719
    ## PC        1  11.12976 0.01604998 0.9157931 0.35529645
    ## SR        1  73.32073 0.10573423 6.0330728 0.01687983
    ## Residual 21 255.21579 0.36804111        NA         NA
    ## Total    36 693.44369 1.00000000        NA         NA

    as.data.frame(cn.trait.pc.sr.geno)

    ##          Df   SumOfSqs         R2          F      Pr(>F)
    ## CT        1  42.567493 0.06138565  3.5025943 0.071539285
    ## pH        1   9.995145 0.01441378  0.8224336 0.375456245
    ## CN        1  19.523186 0.02815396  1.6064324 0.211677883
    ## BR        1  99.979525 0.14417829  8.2266464 0.005749943
    ## PC        1  29.877598 0.04308583  2.4584277 0.121878781
    ## SR        1 140.147628 0.20210383 11.5318108 0.001699983
    ## geno      9  96.137330 0.13863754  0.8789442 0.564134359
    ## Residual 21 255.215787 0.36804111         NA          NA
    ## Total    36 693.443693 1.00000000         NA          NA

    xtable(as.data.frame(cn.geno.trait.pc.sr))

    xtable(as.data.frame(cn.trait.pc.sr.geno))

# How does bark roughness affect the structure of the lichen networks?

Is it driven by lichen abudance (Percent Cover)?

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

    as.data.frame(cn.l.cen.sr)

    ##          Df   SumOfSqs          R2          F      Pr(>F)
    ## L         1 601.441361 0.867325448 376.372605 9.99990e-06
    ## Cen       1  36.581838 0.052753869  22.892343 2.99997e-05
    ## SR        1   2.686674 0.003874394   1.681279 1.95918e-01
    ## Residual 33  52.733819 0.076046289         NA          NA
    ## Total    36 693.443693 1.000000000         NA          NA

    as.data.frame(cn.sr.l.cen)

    ##          Df  SumOfSqs         R2         F       Pr(>F)
    ## SR        1  77.83169 0.11223938  48.70586 0.0000099999
    ## L         1 536.35579 0.77346696 335.64307 0.0000099999
    ## Cen       1  26.52240 0.03824737  16.59730 0.0001299987
    ## Residual 33  52.73382 0.07604629        NA           NA
    ## Total    36 693.44369 1.00000000        NA           NA

    xtable(as.data.frame(cn.l.cen.sr))

    xtable(as.data.frame(cn.sr.l.cen))

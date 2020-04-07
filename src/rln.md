Analysis Summary
================

-   Dead trees and non-lichen species were removed from lichen
    community analyses.
-   Lichen communities were adequately sampled, based on species
    accumulation curves, with moth resistant trees accumulating slightly
    more lichen species.
-   Lichen communities (abundance, richness, diversity, composition)
    were significantly, generally negatively, affected by
    moth susceptibility.
-   Several tree variables, including light availability, leaf litter
    abundance and rock abundance, were impacted by moth susceptibility.
-   Analysis of causal pathways supported an indirect link between moth
    susceptibility and impacts on lichen communities via decreasing
    rock (i.e. habitat) availability through increased leaf abscission
    and accumulation on rocks under trees.
-   These results support a genetically based link between intraspecific
    variation in susceptibility to an insect herbivore and community
    dynamics in an arid ecosystem.
-   Given the possible impacts of cliimate change on this system, this
    study supports the conclusion that community and ecosystem impacts
    need to be considered in an evolutionary context.

<!-- -->

    # 0. Supporting functions and libraries
    ## Support functions

    dif <- function(x){
        out=x[1]
        for (i in 2:length(x)){
            out=out-x[i]
        }
        return(out)
    }

    ## Libraries
    my.libs <- c("vegan", "ecodist")
    if (any(!(my.libs %in% installed.packages()[, 1]))){
        sapply(my.libs[!(my.libs %in% installed.packages()[, 1])], 
               install.packages)
    }else{}
    sapply(my.libs, require, character.only = TRUE)

Load Data
=========

The following are variable descriptions (Variable, Type, Range,
Definition):

-   Moth,categorical,0 or 1,Was the tree susceptible (0) or
    resistant (1) to moth attack?
-   Live/Dead,categorical,0 or 1,Was the tree dead (0) or alive (1)?
-   Litter %,continuous,0 to 100,Percent cover inside quadrat
-   Rocks &gt; 3cm? %,continuous,0 to 100,Percent cover of rocks &gt;
    3cm? inside quadrat
-   Rocks &lt; 3cm? %,continuous,0 to 100,Percent cover of rocks &lt;
    3cm? inside quadrat
-   Shrubs %,continuous,0 to 100,Percent cover of shrubs inside quadrat
-   Grass %,continuous,0 to 100,Percent cover of grass inside quadrat
-   Branches %,continuous,0 to 100,Percent cover of branches on ground
    inside quadrat
-   Distance,continuous,0 to 100,"Distance from main trunk, converted to
    percent of crown radius at that azimuth"
-   Azimuth,continuous,0 to 360,Compass direction from main trunk
-   Slope,continuous,0 to 90,Topographical steepness
-   Aspect,continuous,0 to 360,Compass direction of slope
-   Light,continuous,,Amount of light available to epiliths

<!-- -->

    ## Data are in ../data/scrl
    l.dat <- read.csv("../data/scrl/spp_env_combined.csv")

    ## Summary of data
    summary(l.dat)

    ## remove dead trees
    l.dat <- l.dat[l.dat[, "Live.Dead"] != 0, ]

    ## Lichen species list
    spp.l <- c("Acacon", "Acasup", "Acaobp", "Sterile.sp", "Brown.cr",
    "Lobalp", "Canros", "Calare", "Phydub", "Rhichr", "Xanlin", "Xanpli",
    "Xanele", "GrBr.cr", "Gray.cr")
    spp.moss <- c("Synrur", "Cerpur.Bryarg")

    ## Create a community matrix
    com <- l.dat[, colnames(l.dat) %in% spp.l]
    com.moss <- l.dat[, colnames(l.dat) %in% spp.moss]

    ## Add the tree labels to the rownames
    rownames(com) <- paste(l.dat[, "Moth"], l.dat[, "Tree.pairs"], sep = "_")
    rownames(com.moss) <- paste(l.dat[, "Moth"], l.dat[, "Tree.pairs"], sep = "_")
    rownames(l.dat) <- paste(l.dat[, "Moth"], l.dat[, "Tree.pairs"], sep = "_")

Species accumulation
====================

Are the communities on each tree type adequately sampled?

    spa.all <- specaccum(com)
    spa.res <- specaccum(com[l.dat[, "Moth"] == 0, ])
    spa.sus <- specaccum(com[l.dat[, "Moth"] == 1, ])

    plot(spa.all,
         ylim = c(0, 20),
         xlab = "Cumulative Trees Sampled",
         ylab = "Lichen Species Observed")
    lines(spa.res$sites, spa.res$richness, 
          ylim = c(0, 20), lty = 2, lwd = 3)
    lines(spa.sus$sites, spa.sus$richness, 
          ylim = c(0, 20), lty = 3, lwd = 3)
    legend("bottomright", 
           legend = c("All", "Resistant", "Susceptible"), 
           lty = c(1, 2, 3), lwd = c(1, 2, 2))

![](rln_files/figure-markdown_strict/specacum-1.png)

Moth trees have different microenvironments
===========================================

-   paired t-tests

Moth trees have different lichen communities (FIGURE ch.plot A, R, H)
=====================================================================

less abundant and diverse (paired t-tests, in text)

    abun <- apply(com, 1, sum)
    rich <- apply(com, 1, function(x) sum(sign(x)))
    shan <- apply(com, 1, diversity, index = "shannon")
    tt.a <- t.test(tapply(abun, l.dat[, "Tree.pairs"], diff))
    tt.r <- t.test(tapply(rich, l.dat[, "Tree.pairs"], diff))
    tt.h <- t.test(tapply(shan, l.dat[, "Tree.pairs"], diff))
    tt.arh <- do.call(rbind, 
                      list(a = unlist(tt.a), r = unlist(tt.r), h = unlist(tt.h)))
    data.frame(tt.arh)

    ##         statistic.t parameter.df             p.value          conf.int1
    ## a -2.35680534636893           29  0.0253991007560338  -2.89259563276878
    ## r -2.83579994251995           29 0.00824742800912123  -3.95880113980294
    ## h -2.43278934693583           29  0.0213834528180339 -0.783514237345595
    ##             conf.int2 estimate.mean.of.x null.value.mean            stderr
    ## a  -0.204737700564558  -1.54866666666667               0 0.657104189386137
    ## r  -0.641198860197056               -2.3               0 0.811058624239964
    ## h -0.0678109108683953 -0.425662574106995               0 0.174968940341312
    ##   alternative            method                                 data.name
    ## a   two.sided One Sample t-test tapply(abun, l.dat[, "Tree.pairs"], diff)
    ## r   two.sided One Sample t-test tapply(rich, l.dat[, "Tree.pairs"], diff)
    ## h   two.sided One Sample t-test tapply(shan, l.dat[, "Tree.pairs"], diff)

composition is different (PERMANOVA, in text and supplement)

    com.ds <- cbind(com, ds = rep(0.0001, nrow(com)))
    set.seed(123)
    adonis2(com.ds~ Moth, data = l.dat, 
            strata = l.dat[, "Tree.pairs"],  
            by = "margin", nperm = 100000)

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = com.ds ~ Moth, data = l.dat, by = "margin", strata = l.dat[, "Tree.pairs"], nperm = 1e+05)
    ##          Df SumOfSqs      R2      F Pr(>F)  
    ## Moth      1   0.8281 0.03894 2.3499  0.023 *
    ## Residual 58  20.4394 0.96106                
    ## Total    59  21.2676 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

three main species were reduced by moths (FDR paired t-tests, in text +
supplement)

    ind.spp <- lapply(com, function(x, p) t.test(tapply(x, p, diff)), p = l.dat[, "Tree.pairs"])
    isp <- apply(do.call(rbind, lapply(ind.spp, unlist)), 2, as.numeric)

    ## Warning in apply(do.call(rbind, lapply(ind.spp, unlist)), 2, as.numeric): NAs
    ## introduced by coercion

    ## Warning in apply(do.call(rbind, lapply(ind.spp, unlist)), 2, as.numeric): NAs
    ## introduced by coercion

    ## Warning in apply(do.call(rbind, lapply(ind.spp, unlist)), 2, as.numeric): NAs
    ## introduced by coercion

    rownames(isp) <- names(ind.spp)
    isp[, "p.value"] <- p.adjust(isp[, "p.value"], method = "fdr")
    isp <- isp[order(isp[, "p.value"]), ]
    head(isp[, 1:3])

    ##        statistic.t parameter.df    p.value
    ## Acacon   -3.377629           29 0.01390405
    ## Acasup   -3.242091           29 0.01390405
    ## Canros   -3.581884           29 0.01390405
    ## Lobalp   -2.041361           29 0.17642430
    ## Phydub   -1.922619           29 0.18031798
    ## Calare   -1.607607           29 0.22424946

Litter covering rocks was the main driver (FIGURE = ORDINATION)
===============================================================

light not litter predicted lichen composition (PERMANOVA, table 3,
Ordination)

    set.seed(123)
    adonis2(com.ds ~  Light...average + Litter.., data = l.dat, 
            strata = l.dat[, "Tree.pairs"],  
            by = "margin", nperm = 100000)

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = com.ds ~ Light...average + Litter.., data = l.dat, by = "margin", strata = l.dat[, "Tree.pairs"], nperm = 1e+05)
    ##                 Df SumOfSqs      R2      F Pr(>F)   
    ## Light...average  1   0.4115 0.01935 1.2257  0.240   
    ## Litter..         1   1.0096 0.04747 3.0072  0.007 **
    ## Residual        57  19.1362 0.89979                 
    ## Total           59  21.2676 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    nmds.out <- nmds(vegdist(com.ds), 2, 2)
    ord <- nmds.min(nmds.out, dims = 2)

    ## Minimum stress for given dimensionality:  0.2164016 
    ## r^2 for minimum stress configuration:  0.6474944

    ord.pch <- c("R", "S")[(l.dat[, "Moth"] + 1)]
    plot(X2~ X1, data = ord, pch = ord.pch)

![](rln_files/figure-markdown_strict/ord-com-plot-1.png)

litter not light was correlated with large rocks (dist cor, in text)

    cor.test(tapply(l.dat[, "Big.rocks.."], l.dat[, "Tree.pairs"], diff),
             tapply(l.dat[, "Litter.."], l.dat[, "Tree.pairs"], diff))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  tapply(l.dat[, "Big.rocks.."], l.dat[, "Tree.pairs"], diff) and tapply(l.dat[, "Litter.."], l.dat[, "Tree.pairs"], diff)
    ## t = -11.106, df = 28, p-value = 9.054e-12
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.9530598 -0.8039735
    ## sample estimates:
    ##        cor 
    ## -0.9027609

    cor.test(tapply(l.dat[, "Big.rocks.."], l.dat[, "Tree.pairs"], diff),
             tapply(l.dat[, "Light...average"], l.dat[, "Tree.pairs"], diff))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  tapply(l.dat[, "Big.rocks.."], l.dat[, "Tree.pairs"], diff) and tapply(l.dat[, "Light...average"], l.dat[, "Tree.pairs"], diff)
    ## t = 0.71624, df = 28, p-value = 0.4798
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2376184  0.4716125
    ## sample estimates:
    ##       cor 
    ## 0.1341335

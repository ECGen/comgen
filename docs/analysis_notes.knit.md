# Analysis Notes

* BR = bark roughness
* CT = bark condensed tannins
* CN = bark carbon-nitrogen ratio
* PC = percent cover (percent occupied cells in the sampling grid)
* SR = lichen species richness
* L = size of lichen network (number of nodes)
* Cen = centrality of the lichen network






``` r
## Make sure to run make.R first and load the necessary objects.

library(pacman)
p_load(drake, xtable, rmarkdown, flextable, knitr)
loadd(onc.dat)
loadd(cn.d.onc)
```


## Is the network similarity driven by abundance of lichen?

We found support for the hypothesis that genotypically based variation
in lichen network structure is potentially driven by variation in bark
roughness, not purely as a function of increasing lichen abundance or
richness. In an initital PerMANOVA with sequential partitioning of
variance explained, genotype has a strong, statistically significant
effect when including tree traits and lichen community indices in the
model. In this ordering, condensed tannins and species richness are
both significant, although with low explanatory power. After moving
genotype to the final predictor in the model, genotype is no longer
significant and has low explanatory power, while both bark roughness
and species richness are now significant and have increased in the
variance in network similarity explained. Taken together, these
results suggest that there is potentially a causal relationship
between genotypically explained variance in bark roughness and lichen
community species richness that affects the similarity of lichen
network structure.


``` r
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
```






<!-- html table generated in R 4.5.3 by xtable 1.8-8 package -->
<!-- Tue Jun 16 10:54:13 2026 -->
<table border=1>
<tr> <th> Df </th> <th> SumOfSqs </th> <th> R2 </th> <th> F </th> <th> Pr(&gt;F) </th>  </tr>
  <tr> <td align="right"> 9 </td> <td align="right"> 257.29 </td> <td align="right"> 0.37 </td> <td> 2.35 </td> <td> 0.04023 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 60.20 </td> <td align="right"> 0.09 </td> <td> 4.95 </td> <td> 0.04811 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 7.82 </td> <td align="right"> 0.01 </td> <td> 0.64 </td> <td> 0.444586 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 13.75 </td> <td align="right"> 0.02 </td> <td> 1.13 </td> <td> 0.297517 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 14.72 </td> <td align="right"> 0.02 </td> <td> 1.21 </td> <td> 0.280637 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 11.13 </td> <td align="right"> 0.02 </td> <td> 0.92 </td> <td> 0.352996 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 73.32 </td> <td align="right"> 0.11 </td> <td> 6.03 </td> <td> 0.01695 </td> </tr>
  <tr> <td align="right"> 21 </td> <td align="right"> 255.22 </td> <td align="right"> 0.37 </td> <td>  </td> <td>  </td> </tr>
  <tr> <td align="right"> 36 </td> <td align="right"> 693.44 </td> <td align="right"> 1.00 </td> <td>  </td> <td>  </td> </tr>
   </table>

<!-- html table generated in R 4.5.3 by xtable 1.8-8 package -->
<!-- Tue Jun 16 10:54:13 2026 -->
<table border=1>
<tr> <th> Df </th> <th> SumOfSqs </th> <th> R2 </th> <th> F </th> <th> Pr(&gt;F) </th>  </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 42.57 </td> <td align="right"> 0.06 </td> <td> 3.5 </td> <td> 0.070909 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 10.00 </td> <td align="right"> 0.01 </td> <td> 0.82 </td> <td> 0.378716 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 19.52 </td> <td align="right"> 0.03 </td> <td> 1.61 </td> <td> 0.208418 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 99.98 </td> <td align="right"> 0.14 </td> <td> 8.23 </td> <td> 0.00579 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 29.88 </td> <td align="right"> 0.04 </td> <td> 2.46 </td> <td> 0.120689 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 140.15 </td> <td align="right"> 0.20 </td> <td> 11.53 </td> <td> 0.00173 </td> </tr>
  <tr> <td align="right"> 9 </td> <td align="right"> 96.14 </td> <td align="right"> 0.14 </td> <td> 0.88 </td> <td> 0.563294 </td> </tr>
  <tr> <td align="right"> 21 </td> <td align="right"> 255.22 </td> <td align="right"> 0.37 </td> <td>  </td> <td>  </td> </tr>
  <tr> <td align="right"> 36 </td> <td align="right"> 693.44 </td> <td align="right"> 1.00 </td> <td>  </td> <td>  </td> </tr>
   </table>

\newpage

# How does species richness relate to network structure?

Following the above results showing the effects of bark roughness and
species richness on network structure, an examination of the variaince
in lichen network similarity explained by species richness and the
network metrics of size and centrality, supported the hypothesis that
network structure is both a function of increasing numbers of species
(i.e., nodes in the network) and other structural differences. In
analyzing two PerMANOVAs with these factors, we found that network
size and centrality both explained vaiation in network similarity
regardless of whether species richness was the first or last factor
entered into the model. Additionally, although variation in network
similarity explained by centrality was low, it was significant even
after partitioning the variance explained by network size. Based on
this, we can conclude that network similarity is both of function of
the number of lichen species in the community determining network size
and other aspects of network strcutre, such as centralization.


``` r
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


cn.geno.sr.l.cen  <- vegan::adonis2(
                                cn.d.onc~ geno + SR + L + Cen,
                                by = "term", 
                                data = onc.dat, 
                                mrank = TRUE,
                                permutations = 100000)

cn.sr.l.cen.geno  <- vegan::adonis2(
                           cn.d.onc~ SR + L + Cen + geno,
                           by = "term", 
                           data = onc.dat, 
                           mrank = TRUE,
                           permutations = 100000)
```










<!-- html table generated in R 4.5.3 by xtable 1.8-8 package -->
<!-- Tue Jun 16 10:54:13 2026 -->
<table border=1>
<tr> <th> Df </th> <th> SumOfSqs </th> <th> R2 </th> <th> F </th> <th> Pr(&gt;F) </th>  </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 601.44 </td> <td align="right"> 0.87 </td> <td> 376.37 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 36.58 </td> <td align="right"> 0.05 </td> <td> 22.89 </td> <td> 3e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 2.69 </td> <td align="right"> 0.00 </td> <td> 1.68 </td> <td> 0.195228 </td> </tr>
  <tr> <td align="right"> 33 </td> <td align="right"> 52.73 </td> <td align="right"> 0.08 </td> <td>  </td> <td>  </td> </tr>
  <tr> <td align="right"> 36 </td> <td align="right"> 693.44 </td> <td align="right"> 1.00 </td> <td>  </td> <td>  </td> </tr>
   </table>

<!-- html table generated in R 4.5.3 by xtable 1.8-8 package -->
<!-- Tue Jun 16 10:54:14 2026 -->
<table border=1>
<tr> <th> Df </th> <th> SumOfSqs </th> <th> R2 </th> <th> F </th> <th> Pr(&gt;F) </th>  </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 77.83 </td> <td align="right"> 0.11 </td> <td> 48.71 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 536.36 </td> <td align="right"> 0.77 </td> <td> 335.64 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 26.52 </td> <td align="right"> 0.04 </td> <td> 16.6 </td> <td> 0.00012 </td> </tr>
  <tr> <td align="right"> 33 </td> <td align="right"> 52.73 </td> <td align="right"> 0.08 </td> <td>  </td> <td>  </td> </tr>
  <tr> <td align="right"> 36 </td> <td align="right"> 693.44 </td> <td align="right"> 1.00 </td> <td>  </td> <td>  </td> </tr>
   </table>

<!-- html table generated in R 4.5.3 by xtable 1.8-8 package -->
<!-- Tue Jun 16 10:54:14 2026 -->
<table border=1>
<tr> <th> Df </th> <th> SumOfSqs </th> <th> R2 </th> <th> F </th> <th> Pr(&gt;F) </th>  </tr>
  <tr> <td align="right"> 9 </td> <td align="right"> 257.29 </td> <td align="right"> 0.37 </td> <td> 16.95 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 46.66 </td> <td align="right"> 0.07 </td> <td> 27.66 </td> <td> 2e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 331.43 </td> <td align="right"> 0.48 </td> <td> 196.49 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 17.57 </td> <td align="right"> 0.03 </td> <td> 10.42 </td> <td> 0.00224 </td> </tr>
  <tr> <td align="right"> 24 </td> <td align="right"> 40.48 </td> <td align="right"> 0.06 </td> <td>  </td> <td>  </td> </tr>
  <tr> <td align="right"> 36 </td> <td align="right"> 693.44 </td> <td align="right"> 1.00 </td> <td>  </td> <td>  </td> </tr>
   </table>


<!-- html table generated in R 4.5.3 by xtable 1.8-8 package -->
<!-- Tue Jun 16 10:54:14 2026 -->
<table border=1>
<tr> <th> Df </th> <th> SumOfSqs </th> <th> R2 </th> <th> F </th> <th> Pr(&gt;F) </th>  </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 77.83 </td> <td align="right"> 0.11 </td> <td> 46.14 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 536.36 </td> <td align="right"> 0.77 </td> <td> 317.97 </td> <td> 1e-05 </td> </tr>
  <tr> <td align="right"> 1 </td> <td align="right"> 26.52 </td> <td align="right"> 0.04 </td> <td> 15.72 </td> <td> 0.00029 </td> </tr>
  <tr> <td align="right"> 9 </td> <td align="right"> 12.25 </td> <td align="right"> 0.02 </td> <td> 0.81 </td> <td> 0.625554 </td> </tr>
  <tr> <td align="right"> 24 </td> <td align="right"> 40.48 </td> <td align="right"> 0.06 </td> <td>  </td> <td>  </td> </tr>
  <tr> <td align="right"> 36 </td> <td align="right"> 693.44 </td> <td align="right"> 1.00 </td> <td>  </td> <td>  </td> </tr>
   </table>



<!-- 
rmarkdown::render("./src/lcn_notebook.Rmd", 
    output_format = "pdf_document", 
    output_dir = "./results") 

rmarkdown::render("./src/lcn_notebook.Rmd", 
    output_format = "md_document", 
    output_dir = "./results") 
-->

Results
=======

Tables
------

### Heritability Values

<table>
<thead>
<tr class="header">
<th></th>
<th style="text-align: left;">Response</th>
<th style="text-align: left;">H2</th>
<th style="text-align: left;">R2</th>
<th style="text-align: left;">p-value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>prb.reml.result</td>
<td style="text-align: left;">Percent Rough Bark</td>
<td style="text-align: left;">0.385</td>
<td style="text-align: left;">0.385</td>
<td style="text-align: left;">0</td>
</tr>
<tr class="even">
<td>ph.reml.result</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">0.054</td>
<td style="text-align: left;">0.054</td>
<td style="text-align: left;">0.294</td>
</tr>
<tr class="odd">
<td>ct.reml.result</td>
<td style="text-align: left;">Condensed Tannins (CT)</td>
<td style="text-align: left;">0.28</td>
<td style="text-align: left;">0.28</td>
<td style="text-align: left;">0.014</td>
</tr>
<tr class="even">
<td>cnr.reml.result</td>
<td style="text-align: left;">Carbon-Nitrogen (CN) Ratio</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0.448</td>
</tr>
<tr class="odd">
<td>ptc.reml.result</td>
<td style="text-align: left;">Percent Lichen Cover</td>
<td style="text-align: left;">0.079</td>
<td style="text-align: left;">0.079</td>
<td style="text-align: left;">0.172</td>
</tr>
<tr class="even">
<td>spr.reml.result</td>
<td style="text-align: left;">Lichen Species Richness</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">1</td>
</tr>
<tr class="odd">
<td>spe.reml.result</td>
<td style="text-align: left;">Lichen Species Evenness</td>
<td style="text-align: left;">0.015</td>
<td style="text-align: left;">0.015</td>
<td style="text-align: left;">0.388</td>
</tr>
<tr class="even">
<td>spd.reml.result</td>
<td style="text-align: left;">Lichen Species Diversity</td>
<td style="text-align: left;">0.01</td>
<td style="text-align: left;">0.01</td>
<td style="text-align: left;">0.417</td>
</tr>
<tr class="odd">
<td>link.reml.result</td>
<td style="text-align: left;">Number of Network Links</td>
<td style="text-align: left;">0.07</td>
<td style="text-align: left;">0.07</td>
<td style="text-align: left;">0.238</td>
</tr>
<tr class="even">
<td>mod.reml.result</td>
<td style="text-align: left;">Network Modularity</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">1</td>
</tr>
<tr class="odd">
<td>cen.reml.result</td>
<td style="text-align: left;">Network Centrality</td>
<td style="text-align: left;">0.085</td>
<td style="text-align: left;">0.085</td>
<td style="text-align: left;">0.199</td>
</tr>
<tr class="even">
<td>cn.perm.h2</td>
<td style="text-align: left;">Lichen Network</td>
<td style="text-align: left;">0.16</td>
<td style="text-align: left;">0.233</td>
<td style="text-align: left;">0.025</td>
</tr>
<tr class="odd">
<td>com.perm.h2</td>
<td style="text-align: left;">Community Composition</td>
<td style="text-align: left;">0.052</td>
<td style="text-align: left;">0.173</td>
<td style="text-align: left;">0.102</td>
</tr>
</tbody>
</table>

### Predictors of Lichen Network Similarity

<table>
<thead>
<tr class="header">
<th></th>
<th style="text-align: right;">Df</th>
<th style="text-align: right;">SumOfSqs</th>
<th style="text-align: right;">R2</th>
<th style="text-align: right;">F</th>
<th style="text-align: right;">Pr(&gt;F)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>geno</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">304.927955</td>
<td style="text-align: right;">0.2334811</td>
<td style="text-align: right;">2.365069</td>
<td style="text-align: right;">0.0252975</td>
</tr>
<tr class="even">
<td>BR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">16.259420</td>
<td style="text-align: right;">0.0124497</td>
<td style="text-align: right;">1.261106</td>
<td style="text-align: right;">0.2680732</td>
</tr>
<tr class="odd">
<td>pH</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">5.037083</td>
<td style="text-align: right;">0.0038569</td>
<td style="text-align: right;">0.390684</td>
<td style="text-align: right;">0.5703430</td>
</tr>
<tr class="even">
<td>CN</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">39.666365</td>
<td style="text-align: right;">0.0303722</td>
<td style="text-align: right;">3.076585</td>
<td style="text-align: right;">0.0838916</td>
</tr>
<tr class="odd">
<td>CT</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">70.770152</td>
<td style="text-align: right;">0.0541882</td>
<td style="text-align: right;">5.489043</td>
<td style="text-align: right;">0.0333967</td>
</tr>
<tr class="even">
<td>PC</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">56.352276</td>
<td style="text-align: right;">0.0431485</td>
<td style="text-align: right;">4.370770</td>
<td style="text-align: right;">0.0376962</td>
</tr>
<tr class="odd">
<td>SR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">332.417384</td>
<td style="text-align: right;">0.2545296</td>
<td style="text-align: right;">25.782810</td>
<td style="text-align: right;">0.0001000</td>
</tr>
<tr class="even">
<td>SE</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">55.107744</td>
<td style="text-align: right;">0.0421956</td>
<td style="text-align: right;">4.274242</td>
<td style="text-align: right;">0.0422958</td>
</tr>
<tr class="odd">
<td>Residual</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">425.468500</td>
<td style="text-align: right;">0.3257781</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
<tr class="even">
<td>Total</td>
<td style="text-align: right;">50</td>
<td style="text-align: right;">1306.006880</td>
<td style="text-align: right;">1.0000000</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
</tbody>
</table>

### Predictors of Lichen Community Similarity

<table>
<thead>
<tr class="header">
<th></th>
<th style="text-align: right;">Df</th>
<th style="text-align: right;">SumOfSqs</th>
<th style="text-align: right;">R2</th>
<th style="text-align: right;">F</th>
<th style="text-align: right;">Pr(&gt;F)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>geno</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">1.8466995</td>
<td style="text-align: right;">0.1733428</td>
<td style="text-align: right;">1.3006048</td>
<td style="text-align: right;">0.1018898</td>
</tr>
<tr class="even">
<td>BR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1474919</td>
<td style="text-align: right;">0.0138445</td>
<td style="text-align: right;">1.0387653</td>
<td style="text-align: right;">0.3739626</td>
</tr>
<tr class="odd">
<td>pH</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1302223</td>
<td style="text-align: right;">0.0122235</td>
<td style="text-align: right;">0.9171375</td>
<td style="text-align: right;">0.4566543</td>
</tr>
<tr class="even">
<td>CN</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1651059</td>
<td style="text-align: right;">0.0154979</td>
<td style="text-align: right;">1.1628182</td>
<td style="text-align: right;">0.3102690</td>
</tr>
<tr class="odd">
<td>CT</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1895417</td>
<td style="text-align: right;">0.0177916</td>
<td style="text-align: right;">1.3349157</td>
<td style="text-align: right;">0.2373763</td>
</tr>
<tr class="even">
<td>PC</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2.4602836</td>
<td style="text-align: right;">0.2309377</td>
<td style="text-align: right;">17.3274361</td>
<td style="text-align: right;">0.0001000</td>
</tr>
<tr class="odd">
<td>SR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.5485856</td>
<td style="text-align: right;">0.0514937</td>
<td style="text-align: right;">3.8636125</td>
<td style="text-align: right;">0.0031997</td>
</tr>
<tr class="even">
<td>SE</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.4799261</td>
<td style="text-align: right;">0.0450489</td>
<td style="text-align: right;">3.3800531</td>
<td style="text-align: right;">0.0074993</td>
</tr>
<tr class="odd">
<td>Residual</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">4.6855957</td>
<td style="text-align: right;">0.4398195</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
<tr class="even">
<td>Total</td>
<td style="text-align: right;">50</td>
<td style="text-align: right;">10.6534523</td>
<td style="text-align: right;">1.0000000</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
</tbody>
</table>

Statistical Assumption Checks
-----------------------------

### Shapiro-Wilks Tests for Normality of Residuals

<table>
<thead>
<tr class="header">
<th style="text-align: left;">formula</th>
<th style="text-align: right;">statistic</th>
<th style="text-align: right;">p.value</th>
<th style="text-align: left;">method</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">~I(BR^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.98621</td>
<td style="text-align: right;">0.81358</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="even">
<td style="text-align: left;">~I(pH^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.91192</td>
<td style="text-align: right;">0.00108</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="odd">
<td style="text-align: left;">~I(CT^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.74548</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="even">
<td style="text-align: left;">~I(CN^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.95939</td>
<td style="text-align: right;">0.07855</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="odd">
<td style="text-align: left;">~I(PC^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.78751</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="even">
<td style="text-align: left;">~I(SR^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.71653</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="odd">
<td style="text-align: left;">~I(SE^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.65134</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="even">
<td style="text-align: left;">~I(SD^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.73027</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="odd">
<td style="text-align: left;">~I(L^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.82941</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="even">
<td style="text-align: left;">~I(mod.lik^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.42655</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
<tr class="odd">
<td style="text-align: left;">~I(Cen^(1/4))(1 | geno)</td>
<td style="text-align: right;">0.80978</td>
<td style="text-align: right;">0.00000</td>
<td style="text-align: left;">Shapiro-Wilk normality test</td>
</tr>
</tbody>
</table>

### Fligner Tests for Homogeneity of Variance

<table>
<thead>
<tr class="header">
<th style="text-align: left;">transformation</th>
<th style="text-align: left;">X1</th>
<th style="text-align: left;">X2</th>
<th style="text-align: right;">value.statistic</th>
<th style="text-align: left;">value.parameter</th>
<th style="text-align: right;">value.p.value</th>
<th style="text-align: left;">value.method</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">6.60001</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.76259</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">6.13714</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.80361</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.33549</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59610</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">6.35292</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.78479</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.82097</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54917</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">12.22552</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.27025</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">11.86070</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.29449</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.55860</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67186</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.37451</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59231</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">6.47349</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.77404</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.65884</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56476</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38987</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">6.60001</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.76259</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">6.13714</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.80361</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.33549</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59610</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">6.35292</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.78479</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.82097</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54917</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">12.22552</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.27025</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">11.86070</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.29449</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.55860</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67186</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.37451</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59231</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">6.47349</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.77404</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.65884</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56476</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38987</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">6.60001</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.76259</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">6.13714</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.80361</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.33549</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59610</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">6.35292</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.78479</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.82097</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54917</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">12.22552</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.27025</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">11.86070</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.29449</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.55860</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67186</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.37451</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59231</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">6.47349</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.77404</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.65884</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56476</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38987</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">6.60001</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.76259</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">6.13714</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.80361</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.33549</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59610</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">6.35292</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.78479</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.82097</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54917</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">12.22552</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.27025</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">11.86070</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.29449</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.55860</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67186</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.37451</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59231</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">6.47349</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.77404</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.65884</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56476</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38987</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">6.60001</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.76259</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">6.13714</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.80361</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.33549</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59610</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">6.35292</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.78479</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.82097</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54917</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">12.22552</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.27025</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">11.86070</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.29449</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.55860</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67186</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.37451</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59231</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">6.47349</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.77404</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.65884</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56476</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38987</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">4.18740</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.93850</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">3.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.96375</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.69799</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56099</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.07641</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43382</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">15.38571</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.11862</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">14.43681</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.15398</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">17.89448</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.05677</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.82940</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64550</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.85431</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45337</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.89526</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.63907</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">13.60700</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.19168</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">11.67367</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.30749</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">4.18740</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.93850</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">3.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.96375</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.69799</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56099</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.07641</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43382</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">15.38571</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.11862</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">14.43681</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.15398</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">17.89448</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.05677</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.82940</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64550</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.85431</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45337</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.89526</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.63907</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">13.60700</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.19168</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">11.67367</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.30749</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">4.18740</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.93850</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">3.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.96375</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.69799</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56099</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.07641</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43382</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">15.38571</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.11862</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">14.43681</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.15398</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">17.89448</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.05677</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.82940</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64550</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.85431</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45337</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.89526</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.63907</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">13.60700</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.19168</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">11.67367</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.30749</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">4.18740</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.93850</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">3.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.96375</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.69799</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56099</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.07641</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43382</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">15.38571</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.11862</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">14.43681</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.15398</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">17.89448</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.05677</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.82940</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64550</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.85431</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45337</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.89526</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.63907</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">13.60700</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.19168</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">11.67367</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.30749</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">4.18740</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.93850</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">3.59576</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.96375</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.69799</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56099</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.07641</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43382</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">15.38571</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.11862</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">14.43681</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.15398</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">17.89448</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.05677</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.82940</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64550</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.85431</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45337</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.89526</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.63907</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">13.60700</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.19168</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">y2</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">11.67367</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.30749</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.55162</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.48067</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.83365</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54796</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.66680</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38406</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.65092</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.84369</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.70093</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83973</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.68311</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56242</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.81387</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54985</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.53952</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67372</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.09015</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62003</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.79771</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83196</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.03537</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62538</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.25274</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41861</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.55162</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.48067</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.83365</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54796</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.66680</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38406</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.65092</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.84369</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.70093</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83973</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.68311</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56242</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.81387</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54985</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.53952</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67372</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.09015</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62003</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.79771</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83196</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.03537</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62538</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.25274</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41861</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.55162</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.48067</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.83365</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54796</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.66680</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38406</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.65092</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.84369</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.70093</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83973</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.68311</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56242</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.81387</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54985</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.53952</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67372</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.09015</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62003</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.79771</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83196</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.03537</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62538</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.25274</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41861</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.55162</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.48067</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.83365</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54796</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.66680</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38406</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.65092</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.84369</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.70093</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83973</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.68311</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56242</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.81387</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54985</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.53952</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67372</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.09015</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62003</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.79771</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83196</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.03537</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62538</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.25274</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41861</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.55162</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.48067</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.83365</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54796</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.66680</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.38406</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.65092</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.84369</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.70093</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83973</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.68311</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.56242</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.81387</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54985</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.53952</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67372</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.09015</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62003</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.79771</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.83196</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.03537</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.62538</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r2y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.25274</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41861</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">11.50535</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.31952</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.59483</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.47673</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.01958</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43878</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.08250</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.88560</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.47242</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85747</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.81243</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54999</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.27900</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.50585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.53934</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.85372</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64312</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.50477</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.00275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.72519</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.17797</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.42502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">11.50535</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.31952</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.59483</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.47673</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.01958</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43878</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.08250</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.88560</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.47242</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85747</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.81243</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54999</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.27900</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.50585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.53934</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.85372</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64312</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.50477</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.00275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.72519</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.17797</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.42502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">11.50535</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.31952</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.59483</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.47673</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.01958</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43878</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.08250</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.88560</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.47242</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85747</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.81243</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54999</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.27900</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.50585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.53934</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.85372</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64312</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.50477</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.00275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.72519</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.17797</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.42502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">11.50535</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.31952</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.59483</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.47673</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.01958</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43878</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.08250</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.88560</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.47242</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85747</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.81243</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54999</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.27900</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.50585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.53934</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.85372</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64312</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.50477</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.00275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.72519</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.17797</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.42502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">11.50535</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.31952</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.59483</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.47673</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.01958</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.43878</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.08250</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.88560</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.47242</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85747</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.81243</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.54999</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.27900</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.50585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.53934</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.85372</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.64312</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.50477</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.00275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.72519</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">r4y</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.17797</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.42502</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">12.68909</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.24158</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.93983</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.44579</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.88275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45084</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">4.97993</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.89251</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.96171</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.81847</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.37692</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49674</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.13583</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.51926</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.53939</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">8.35121</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59457</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">5.49108</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85606</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">7.46756</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.68069</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">statistic</td>
<td style="text-align: right;">10.26625</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41745</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">12.68909</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.24158</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.93983</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.44579</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.88275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45084</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">4.97993</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.89251</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.96171</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.81847</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.37692</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49674</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.13583</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.51926</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.53939</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">8.35121</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59457</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">5.49108</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85606</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">7.46756</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.68069</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">parameter</td>
<td style="text-align: right;">10.26625</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41745</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">12.68909</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.24158</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.93983</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.44579</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.88275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45084</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">4.97993</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.89251</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.96171</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.81847</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.37692</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49674</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.13583</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.51926</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.53939</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">8.35121</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59457</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">5.49108</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85606</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">7.46756</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.68069</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">p.value</td>
<td style="text-align: right;">10.26625</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41745</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">12.68909</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.24158</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.93983</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.44579</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.88275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45084</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">4.97993</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.89251</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.96171</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.81847</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.37692</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49674</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.13583</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.51926</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.53939</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">8.35121</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59457</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">5.49108</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85606</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">7.46756</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.68069</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">method</td>
<td style="text-align: right;">10.26625</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41745</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">PC</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">12.68909</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.24158</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.93983</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.44579</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SD</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.88275</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.45084</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">SE</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">4.97993</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.89251</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">BR</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.96171</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.81847</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">L</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.37692</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49674</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">Cen</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.13583</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.51926</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">mod.lik</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">9.38661</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.49585</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">C</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.53939</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.67373</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">N</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">8.35121</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.59457</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CN</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">5.49108</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.85606</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="even">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">CT</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">7.46756</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.68069</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
<tr class="odd">
<td style="text-align: left;">logy</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">data.name</td>
<td style="text-align: right;">10.26625</td>
<td style="text-align: left;">c(df = 10)</td>
<td style="text-align: right;">0.41745</td>
<td style="text-align: left;">Fligner-Killeen test of homogeneity of variances</td>
</tr>
</tbody>
</table>

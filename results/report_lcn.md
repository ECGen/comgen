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
<td style="text-align: left;">0.396</td>
<td style="text-align: left;">0.396</td>
<td style="text-align: left;">0.001</td>
</tr>
<tr class="even">
<td>ph.reml.result</td>
<td style="text-align: left;">pH</td>
<td style="text-align: left;">0.14</td>
<td style="text-align: left;">0.14</td>
<td style="text-align: left;">0.21</td>
</tr>
<tr class="odd">
<td>ct.reml.result</td>
<td style="text-align: left;">Condensed Tannins (CT)</td>
<td style="text-align: left;">0.256</td>
<td style="text-align: left;">0.256</td>
<td style="text-align: left;">0.021</td>
</tr>
<tr class="even">
<td>cnr.reml.result</td>
<td style="text-align: left;">Carbon-Nitrogen (CN) Ratio</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0.459</td>
</tr>
<tr class="odd">
<td>ptc.reml.result</td>
<td style="text-align: left;">Percent Lichen Cover</td>
<td style="text-align: left;">0.137</td>
<td style="text-align: left;">0.137</td>
<td style="text-align: left;">0.085</td>
</tr>
<tr class="even">
<td>spr.reml.result</td>
<td style="text-align: left;">Lichen Species Richness</td>
<td style="text-align: left;">0.028</td>
<td style="text-align: left;">0.028</td>
<td style="text-align: left;">0.348</td>
</tr>
<tr class="odd">
<td>spe.reml.result</td>
<td style="text-align: left;">Lichen Species Evenness</td>
<td style="text-align: left;">0.057</td>
<td style="text-align: left;">0.057</td>
<td style="text-align: left;">0.24</td>
</tr>
<tr class="even">
<td>spd.reml.result</td>
<td style="text-align: left;">Lichen Species Diversity</td>
<td style="text-align: left;">0.029</td>
<td style="text-align: left;">0.029</td>
<td style="text-align: left;">0.358</td>
</tr>
<tr class="odd">
<td>link.reml.result</td>
<td style="text-align: left;">Number of Network Links</td>
<td style="text-align: left;">0.124</td>
<td style="text-align: left;">0.124</td>
<td style="text-align: left;">0.124</td>
</tr>
<tr class="even">
<td>mod.reml.result</td>
<td style="text-align: left;">Network Modularity</td>
<td style="text-align: left;">0.005</td>
<td style="text-align: left;">0.005</td>
<td style="text-align: left;">0.444</td>
</tr>
<tr class="odd">
<td>cen.reml.result</td>
<td style="text-align: left;">Network Centrality</td>
<td style="text-align: left;">0.128</td>
<td style="text-align: left;">0.128</td>
<td style="text-align: left;">0.116</td>
</tr>
<tr class="even">
<td>com.perm.h2</td>
<td style="text-align: left;">Community Composition</td>
<td style="text-align: left;">0.161</td>
<td style="text-align: left;">0.236</td>
<td style="text-align: left;">0.011</td>
</tr>
<tr class="odd">
<td>cn.perm.h2</td>
<td style="text-align: left;">Lichen Network</td>
<td style="text-align: left;">0.069</td>
<td style="text-align: left;">0.269</td>
<td style="text-align: left;">0.03</td>
</tr>
</tbody>
</table>

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

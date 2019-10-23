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
<td style="text-align: right;">12</td>
<td style="text-align: right;">367.646117</td>
<td style="text-align: right;">0.2693747</td>
<td style="text-align: right;">2.3065203</td>
<td style="text-align: right;">0.0275972</td>
</tr>
<tr class="even">
<td>BR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">23.633515</td>
<td style="text-align: right;">0.0173163</td>
<td style="text-align: right;">1.7792496</td>
<td style="text-align: right;">0.1848815</td>
</tr>
<tr class="odd">
<td>pH</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">8.959401</td>
<td style="text-align: right;">0.0065646</td>
<td style="text-align: right;">0.6745087</td>
<td style="text-align: right;">0.4177582</td>
</tr>
<tr class="even">
<td>CN</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">37.695027</td>
<td style="text-align: right;">0.0276192</td>
<td style="text-align: right;">2.8378707</td>
<td style="text-align: right;">0.0944906</td>
</tr>
<tr class="odd">
<td>CT</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">76.220488</td>
<td style="text-align: right;">0.0558468</td>
<td style="text-align: right;">5.7382607</td>
<td style="text-align: right;">0.0286971</td>
</tr>
<tr class="even">
<td>PC</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">28.502578</td>
<td style="text-align: right;">0.0208839</td>
<td style="text-align: right;">2.1458170</td>
<td style="text-align: right;">0.1426857</td>
</tr>
<tr class="odd">
<td>SR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">332.226750</td>
<td style="text-align: right;">0.2434229</td>
<td style="text-align: right;">25.0116964</td>
<td style="text-align: right;">0.0001000</td>
</tr>
<tr class="even">
<td>SE</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">51.594856</td>
<td style="text-align: right;">0.0378036</td>
<td style="text-align: right;">3.8843196</td>
<td style="text-align: right;">0.0494951</td>
</tr>
<tr class="odd">
<td>Residual</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">438.334233</td>
<td style="text-align: right;">0.3211680</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
<tr class="even">
<td>Total</td>
<td style="text-align: right;">52</td>
<td style="text-align: right;">1364.812965</td>
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
<th style="text-align: left;">Response</th>
<th style="text-align: left;">H2</th>
<th style="text-align: left;">R2</th>
<th style="text-align: left;">p-value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>ct.reml.result</td>
<td style="text-align: left;">Condensed Tannins (CT)</td>
<td style="text-align: left;">0.25556</td>
<td style="text-align: left;">0.25556</td>
<td style="text-align: left;">0.02</td>
</tr>
<tr class="even">
<td>cen.reml.result</td>
<td style="text-align: left;">Network Centrality</td>
<td style="text-align: left;">0.20166</td>
<td style="text-align: left;">0.20166</td>
<td style="text-align: left;">0.03942</td>
</tr>
<tr class="odd">
<td>link.reml.result</td>
<td style="text-align: left;">Number of Network Links</td>
<td style="text-align: left;">0.17016</td>
<td style="text-align: left;">0.17016</td>
<td style="text-align: left;">0.06636</td>
</tr>
<tr class="even">
<td>prb.reml.result</td>
<td style="text-align: left;">Lichen Community Composition</td>
<td style="text-align: left;">0.16093</td>
<td style="text-align: left;">0.24287</td>
<td style="text-align: left;">0.0038</td>
</tr>
<tr class="odd">
<td>ptc.reml.result</td>
<td style="text-align: left;">Percent Lichen Cover</td>
<td style="text-align: left;">0.1368</td>
<td style="text-align: left;">0.1368</td>
<td style="text-align: left;">0.0859</td>
</tr>
<tr class="even">
<td>ph.reml.result</td>
<td style="text-align: left;">Lichen Network</td>
<td style="text-align: left;">0.06385</td>
<td style="text-align: left;">0.26937</td>
<td style="text-align: left;">0.0276</td>
</tr>
<tr class="odd">
<td>spe.reml.result</td>
<td style="text-align: left;">Lichen Species Evenness</td>
<td style="text-align: left;">0.05732</td>
<td style="text-align: left;">0.05732</td>
<td style="text-align: left;">0.2383</td>
</tr>
<tr class="even">
<td>mod.reml.result</td>
<td style="text-align: left;">Network Modularity</td>
<td style="text-align: left;">0.05731</td>
<td style="text-align: left;">0.05731</td>
<td style="text-align: left;">0.278</td>
</tr>
<tr class="odd">
<td>spd.reml.result</td>
<td style="text-align: left;">Lichen Species Diversity</td>
<td style="text-align: left;">0.02908</td>
<td style="text-align: left;">0.02908</td>
<td style="text-align: left;">0.35</td>
</tr>
<tr class="even">
<td>spr.reml.result</td>
<td style="text-align: left;">Lichen Species Richness</td>
<td style="text-align: left;">0.02807</td>
<td style="text-align: left;">0.02807</td>
<td style="text-align: left;">0.3447</td>
</tr>
<tr class="odd">
<td>cnr.reml.result</td>
<td style="text-align: left;">Carbon-Nitrogen (CN) Ratio</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">0.4601</td>
</tr>
</tbody>
</table>

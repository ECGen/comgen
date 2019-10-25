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
<td style="text-align: right;">12</td>
<td style="text-align: right;">367.646117</td>
<td style="text-align: right;">0.2693747</td>
<td style="text-align: right;">2.3065203</td>
<td style="text-align: right;">0.0298970</td>
</tr>
<tr class="even">
<td>BR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">23.633515</td>
<td style="text-align: right;">0.0173163</td>
<td style="text-align: right;">1.7792496</td>
<td style="text-align: right;">0.1892811</td>
</tr>
<tr class="odd">
<td>pH</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">8.959401</td>
<td style="text-align: right;">0.0065646</td>
<td style="text-align: right;">0.6745087</td>
<td style="text-align: right;">0.4121588</td>
</tr>
<tr class="even">
<td>CN</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">37.695027</td>
<td style="text-align: right;">0.0276192</td>
<td style="text-align: right;">2.8378707</td>
<td style="text-align: right;">0.0951905</td>
</tr>
<tr class="odd">
<td>CT</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">76.220488</td>
<td style="text-align: right;">0.0558468</td>
<td style="text-align: right;">5.7382607</td>
<td style="text-align: right;">0.0304970</td>
</tr>
<tr class="even">
<td>PC</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">28.502578</td>
<td style="text-align: right;">0.0208839</td>
<td style="text-align: right;">2.1458170</td>
<td style="text-align: right;">0.1455854</td>
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
<td style="text-align: right;">0.0517948</td>
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
<td style="text-align: right;">2.5029970</td>
<td style="text-align: right;">0.2359552</td>
<td style="text-align: right;">1.6694092</td>
<td style="text-align: right;">0.0105989</td>
</tr>
<tr class="even">
<td>BR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1175072</td>
<td style="text-align: right;">0.0110773</td>
<td style="text-align: right;">0.9404767</td>
<td style="text-align: right;">0.4254575</td>
</tr>
<tr class="odd">
<td>pH</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.0438228</td>
<td style="text-align: right;">0.0041311</td>
<td style="text-align: right;">0.3507391</td>
<td style="text-align: right;">0.8948105</td>
</tr>
<tr class="even">
<td>CN</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1413908</td>
<td style="text-align: right;">0.0133288</td>
<td style="text-align: right;">1.1316313</td>
<td style="text-align: right;">0.3129687</td>
</tr>
<tr class="odd">
<td>CT</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1805161</td>
<td style="text-align: right;">0.0170171</td>
<td style="text-align: right;">1.4447732</td>
<td style="text-align: right;">0.2005799</td>
</tr>
<tr class="even">
<td>PC</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2.4746960</td>
<td style="text-align: right;">0.2332873</td>
<td style="text-align: right;">19.8064012</td>
<td style="text-align: right;">0.0001000</td>
</tr>
<tr class="odd">
<td>SR</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.5885188</td>
<td style="text-align: right;">0.0554791</td>
<td style="text-align: right;">4.7102507</td>
<td style="text-align: right;">0.0013999</td>
</tr>
<tr class="even">
<td>SE</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.4353260</td>
<td style="text-align: right;">0.0410378</td>
<td style="text-align: right;">3.4841614</td>
<td style="text-align: right;">0.0089991</td>
</tr>
<tr class="odd">
<td>Residual</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">4.1231604</td>
<td style="text-align: right;">0.3886864</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
<tr class="even">
<td>Total</td>
<td style="text-align: right;">52</td>
<td style="text-align: right;">10.6079351</td>
<td style="text-align: right;">1.0000000</td>
<td style="text-align: right;">NA</td>
<td style="text-align: right;">NA</td>
</tr>
</tbody>
</table>

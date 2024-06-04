---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">vcon</span></code><a class="headerlink" href="#vcon" title="Permalink to this headline">¶</a></h1>
<p>Variable thermal conductivity information.</p>
<ol class="arabic simple">
<li>Thermal conductivity for intact salt: <span class="math notranslate nohighlight">\(\lambda_{IS} = \lambda_{IS,300}(300/T)^{\gamma_1}\)</span></li>
<li>Thermal conductivity for crushed salt: <span class="math notranslate nohighlight">\(\lambda_{CS} = \lambda_{CS,300}(300/T)^{\gamma_1}\)</span></li>
</ol>
<p>where:</p>
<p><span class="math notranslate nohighlight">\(\lambda_{CS,300} = 1.08 \left( \alpha_4 \phi^4 + \alpha_3 \phi^3 + \alpha_2 \phi^2 + \alpha_1 \phi + \alpha_0 \right)\)</span></p>
<p>Parameters related to thermal conductivity are
<span class="math notranslate nohighlight">\(\lambda_{IS,300},\gamma_1,\gamma_2,\alpha_4,\alpha_3,\alpha_2,\alpha_1,\alpha_0\)</span> and <span class="math notranslate nohighlight">\(\phi\)</span>.
An additional parameter is the specific heat of salt.</p>
<ul class="simple">
<li>Group 1 - IVCON(I), VC1F(I), VC2F(I), VC3F(I)</li>
<li>Group 2 - JA, JB, JC, IVCND</li>
</ul>
<p>The parameter (I) is incremented each time Group 1 is read. Group 2 lines will refer to this parameter. Group 1 is ended with a blank line.</p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="9%" />
<col width="75%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IVCON(i)</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Model type for ith conductivity model.</div>
<div class="line">IVCON(i) = 1, linear variation of thermal conductivity with temperature.</div>
<div class="line">IVCON(i) = 2, square root variation of thermal conductivity with liquid
saturation.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>VC1F(i)</td>
<td>real</td>
<td>Reference temperature (oC) for IVCON(i) = 1.
Conductivity <span class="math notranslate nohighlight">\(\left( \frac{W}{m \cdot K} \right)\)</span> at liquid
saturation = 1 for IVCON(i) = 2.</td>
</tr>
<tr class="row-even"><td>VC2F(i)</td>
<td>real</td>
<td>Reference conductivity ( <span class="math notranslate nohighlight">\(^oC\)</span> ) for IVCON(i) = 1. Conductivity
<span class="math notranslate nohighlight">\(\left( \frac{W}{m \cdot K} \right)\)</span> at liquid saturation = 0 for
IVCON(i) = 2.</td>
</tr>
<tr class="row-odd"><td>VC3F(i)</td>
<td>real</td>
<td>Change in conductivity with respect to temperature for IVCON(i) = 1.
Not used for IVCON(i) = 2.</td>
</tr>
<tr class="row-even"><td>IVCND</td>
<td>integer</td>
<td>Number referring to the sequence of models read in Group 1.
The default is 1.</td>
</tr>
</tbody>
</table>
<p>The following are examples of <code class="docutils literal notranslate"><span class="pre">vcon</span></code>. In the first example,
a linear conductivity model is defined and applied at each node.
The reference temperature is <span class="math notranslate nohighlight">\(20^o C\)</span>, the reference
conductivity is <span class="math notranslate nohighlight">\(1\frac{W}{m \cdot K}\)</span>, and the change
in conductivity with temperature is 0.01.</p>
<table border="1" class="docutils">
<colgroup>
<col width="23%" />
<col width="26%" />
<col width="23%" />
<col width="29%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">vcon</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td>20.0</td>
<td>1.0</td>
<td>1.e-2</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>For the second example, three models are defined for the entire domain.
Model 1 defines the constant thermal conductivity of 16.26 W/m K at 26.85 °C
(=300 °K) for a stainless steel canister (zone 1).</p>
<p>Model 2 defines all parameters for a crushed salt (zone 2).</p>
<p>Model 3 defines reference thermal conductivity of 5.4 W/m K at 26.85 °C
(=300 °K) and exponent 1.14 for the intact salt in the rest of the domain.</p>
<table border="1" class="docutils">
<colgroup>
<col width="11%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="9%" />
<col width="9%" />
<col width="11%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">vcon</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td>26.85</td>
<td>16.26</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>4</td>
<td>26.85</td>
<td>1.08</td>
<td>270.0</td>
<td>370.0</td>
<td>136.0</td>
<td>1.5</td>
<td>5.0</td>
<td>1.14</td>
</tr>
<tr class="row-even"><td>3</td>
<td>26.85</td>
<td>5.4</td>
<td>1.14</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0</td>
<td>0</td>
<td>3</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>-1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>-2</td>
<td>0</td>
<td>0</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
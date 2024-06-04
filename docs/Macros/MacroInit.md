---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">init</span></code><a class="headerlink" href="#init" title="Permalink to this headline">¶</a></h1>
<p>Set initial pressure and temperature at all nodes.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Required if macro <code class="docutils literal notranslate"><span class="pre">pres</span></code> not used.</p>
</div>
<p>Group 1 -       PEIN, TIN, TIN1, GRAD1, DEPTH, TIN2, GRAD2, QUAD</p>
<p>Note that the macro pres may overwrite some of the values that are set by macro init.</p>
<table border="1" class="docutils">
<colgroup>
<col width="4%" />
<col width="2%" />
<col width="94%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>PEIN</td>
<td>real</td>
<td>Initial value of pressure (MPa). If initial values are read from the read file (iread), then this value is ignored. If gravity is present, this is the value of the pressure at node 1, and the other nodal pressures are adjusted by applying the hydraulic head. Absolute pressures are used. Pressure as a function of depth is calculated with TIN &lt; 0.</td>
</tr>
<tr class="row-odd"><td>TIN</td>
<td>real</td>
<td>Initial value of temperature (oC). If TIN ≤ 0, then the initial temperatures are calculated using the temperature gradient formulas given below.</td>
</tr>
<tr class="row-even"><td>TIN1</td>
<td>real</td>
<td>Defined in formulas below (oC)</td>
</tr>
<tr class="row-odd"><td>GRAD1</td>
<td>real</td>
<td>Defined in formulas below (oC/m)</td>
</tr>
<tr class="row-even"><td>DEPTH</td>
<td>real</td>
<td>Defined in formulas below (m)</td>
</tr>
<tr class="row-odd"><td>TIN2</td>
<td>real</td>
<td>Defined in formulas below (oC)</td>
</tr>
<tr class="row-even"><td>GRAD2</td>
<td>real</td>
<td>Defined in formulas below (oC/m)</td>
</tr>
<tr class="row-odd"><td>QUAD</td>
<td>real</td>
<td>Defined in formulas below (oC/m2)</td>
</tr>
</tbody>
</table>
<p><span class="math notranslate nohighlight">\(T = TIN1 + GRAD1 \times Z, 0 \le Z \le DEPTH\)</span></p>
<p><span class="math notranslate nohighlight">\(T = TIN2 + GRAD2 \times Z + QUAD \times Z^2, Z \gt DEPTH\)</span></p>
<p>The following are examples of init. In the first example, the initial pressure is 3.6
MPa and the initial temperature is 240 C over the entire range of depth for the model.</p>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="12%" />
<col width="15%" />
<col width="12%" />
<col width="12%" />
<col width="15%" />
<col width="10%" />
<col width="10%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>init</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>3.6</td>
<td>0.0</td>
<td><ol class="first last arabic simple" start="240">
<li></li>
</ol>
</td>
<td>0.0</td>
<td>0.0</td>
<td><ol class="first last arabic simple" start="240">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
</tbody>
</table>
<p>In the second example, the initial pressure is 5.0 MPa and the initial temperature
field is defined using a surface temperature of 20 C and linear gradient of 0.3
C/m for depths ranging from 0 - 2500 m.</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="12%" />
<col width="12%" />
<col width="14%" />
<col width="16%" />
<col width="12%" />
<col width="12%" />
<col width="9%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>init</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>5.0</td>
<td>0.0</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>00.3</td>
<td><ol class="first last arabic simple" start="2500">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>0.3</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
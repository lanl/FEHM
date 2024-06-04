---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">rich</span></code><a class="headerlink" href="#rich" title="Permalink to this headline">¶</a></h1>
<p>Invokes Richard’s equation solution for unsaturated-saturated flow. A single
phase approach that neglects air phase flow and assumes the movement of water
is independent of air flow and pressure.</p>
<p>Uses variable switching (Pressure, Saturation).</p>
<ul class="simple">
<li>Group 1 - STRD_RICH, TOL_PHASE, PCHNG, SCHNG</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="23%" />
<col width="12%" />
<col width="65%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>STRD_RICH</td>
<td>real</td>
<td>Newton-raphson relaxation factor</td>
</tr>
<tr class="row-odd"><td>TOL_PHASE</td>
<td>real</td>
<td>Tolerance for full saturation</td>
</tr>
<tr class="row-even"><td>PCHNG</td>
<td>real</td>
<td>Pressure adjustment after variable switch</td>
</tr>
<tr class="row-odd"><td>SCHNG</td>
<td>real</td>
<td>Saturation adjustment after variable switch</td>
</tr>
</tbody>
</table>
<p>The following is an example of rich.</p>
<table border="1" class="docutils">
<colgroup>
<col width="22%" />
<col width="26%" />
<col width="26%" />
<col width="26%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>rich</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>0.95</td>
<td>1.e-5</td>
<td>1.e-3</td>
<td>1.e-3</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
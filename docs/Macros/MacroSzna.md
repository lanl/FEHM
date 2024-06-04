---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">szna</span> <span class="pre">or</span> <span class="pre">napl</span></code><a class="headerlink" href="#szna-or-napl" title="Permalink to this headline">¶</a></h1>
<ul class="simple">
<li>Group 1 - ICO2D</li>
<li>Group 2 - TREF, PREF</li>
<li>Group 3 - DENNAPL, VISCNAPL</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="7%" />
<col width="4%" />
<col width="88%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ICO2D</td>
<td>integer</td>
<td>Determines the type of air module used. ICO2D = 1, 1 degree of freedom solution to the saturated-unsaturated problem is produced. This formulation is similar to the RIchard’s Equation.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>ICO2D = 2, 1 degree of freedom solution is obtained assuming only gas flow with no liquid present.ICO2D = 3, full 2 degree of freedom solution.All other values are ignored. The default is 3.</td>
</tr>
<tr class="row-even"><td>TREF</td>
<td>real</td>
<td>Reference temperature for properties (oC).</td>
</tr>
<tr class="row-odd"><td>PREF</td>
<td>real</td>
<td>Reference pressure for properties (MPa).</td>
</tr>
<tr class="row-even"><td>DENNAPL</td>
<td>real</td>
<td>NAPL density (kg/m3).</td>
</tr>
<tr class="row-odd"><td>VISCNAPL</td>
<td>real</td>
<td>NAPL viscosity (Pa s).</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
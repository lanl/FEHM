---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">cond</span></code><a class="headerlink" href="#cond" title="Permalink to this headline">¶</a></h1>
<p>Assign thermal conductivities of the rock.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last"><code class="docutils literal notranslate"><span class="pre">cond</span></code> is required for non-isothermal problems.</p>
</div>
<p><a class="reference external" href="InputData.html#JA">Group 1 - JA, JB, JC, THXD, THYD, THZD</a></p>
<table border="1" class="docutils">
<colgroup>
<col width="22%" />
<col width="11%" />
<col width="12%" />
<col width="55%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>THXD</td>
<td>real</td>
<td>1.e-30</td>
<td>Thermal conductivity in the x-direction</td>
</tr>
<tr class="row-odd"><td>THYD</td>
<td>real</td>
<td>1.e-30</td>
<td>Thermal conductivity in the y-direction</td>
</tr>
<tr class="row-even"><td>THZD</td>
<td>real</td>
<td>1.e-30</td>
<td>Thermal conductivity in the z-direction</td>
</tr>
</tbody>
</table>
<p>The following is an example of cond. In this example all the nodes numbered 1 through 140 have thermal conductivities of 1 in the X and Y directions, and 0 in the Z direction.</p>
<table border="1" class="docutils">
<colgroup>
<col width="13%" />
<col width="11%" />
<col width="7%" />
<col width="22%" />
<col width="24%" />
<col width="22%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">cond</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td>140</td>
<td>1</td>
<td>1.00e-00</td>
<td>1e.00e-00</td>
<td>0.00e-00</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
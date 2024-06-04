---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">hyco</span></code><a class="headerlink" href="#hyco" title="Permalink to this headline">¶</a></h1>
<p>Hydraulic conductivity input.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Required if macro perm not used.</p>
</div>
<table border="1" class="docutils">
<colgroup>
<col width="19%" />
<col width="10%" />
<col width="11%" />
<col width="60%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>PNX</td>
<td>real</td>
<td>1.e-30</td>
<td>Hydraulic conductivity in the x-direction (m/s).</td>
</tr>
<tr class="row-odd"><td>PNY</td>
<td>real</td>
<td>1.e-30</td>
<td>Hydraulic conductivity in the y-direction (m/s).</td>
</tr>
<tr class="row-even"><td>PNZ</td>
<td>real</td>
<td>1.e-30</td>
<td>Hydraulic conductivity in the z-direction (m/s).</td>
</tr>
</tbody>
</table>
<p>The following is an example of the <code class="docutils literal notranslate"><span class="pre">hyco</span></code> macro. In this example, nodes 1 through 140 are specified to have hydraulic conductivities in the X, Y, and Z directions of 1.0e-5, 1.0e-5, and 0. m/s respectively.</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="11%" />
<col width="7%" />
<col width="23%" />
<col width="23%" />
<col width="23%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>hyco</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>140</td>
<td>1</td>
<td>1.00e-05</td>
<td>1.00e-05</td>
<td>0.00e-00</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
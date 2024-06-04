---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">exrl</span></code><a class="headerlink" href="#exrl" title="Permalink to this headline">¶</a></h1>
<p>Allows the user to choose linearized relative permeability. The linearized relative
permeability is formed using the nonlinear value of relative permeability at the
iteration number IEXRLP. After that iteration a relative permeability value based
on a Taylor series expansion in saturation is used.</p>
<ul class="simple">
<li>Group 1 - IEXRLP</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="9%" />
<col width="76%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IEXRLP</td>
<td>integer</td>
<td>If IEXRLP ≥ 1, then linearized relative permeability. Otherwise not enabled.</td>
</tr>
</tbody>
</table>
<p>In the following example of <code class="docutils literal notranslate"><span class="pre">exrl</span></code>, the user enables linearized relative permeability at the first iteration.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>exrl1</td>
</tr>
<tr class="row-even"><td>1</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
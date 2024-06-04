---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">rock</span></code><a class="headerlink" href="#rock" title="Permalink to this headline">¶</a></h1>
<p>Assign rock density, specific heat and porosity.</p>
<table border="1" class="docutils">
<colgroup>
<col width="13%" />
<col width="6%" />
<col width="81%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>DENRD</td>
<td>real</td>
<td>Rock density (kg/m3).</td>
</tr>
<tr class="row-odd"><td>CPRD</td>
<td>real</td>
<td>Rock specific heat <span class="math notranslate nohighlight">\(\left( \frac{MJ}{kg \cdot K} \right)\)</span>.
If CPRD &gt; 1 the code will assume the units are
<span class="math notranslate nohighlight">\(\left( \frac{J}{kg \cdot K} \right)\)</span> and multiply by <span class="math notranslate nohighlight">\(10^{-6}\)</span>.</td>
</tr>
<tr class="row-even"><td>PSD</td>
<td>real</td>
<td><div class="first line-block">
<div class="line">Porosity.</div>
</div>
<p class="last">Special note on negative porosities. If the code encounters a negative porosity,
the node at which the negative porosity occurs is effectively removed from the model.
That is, the geometric connections from that node to other nodes in the model are
removed. The volume associated with the node acts as a barrier to flow. For input purposes,
the node may still be assigned properties, though they will have no effect on the simulation results.</p>
</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">rock</span></code>. In this example the nodes numbered
1 through 140 are assigned a rock density of 2563. kg/m^3,
a rock specific heat of 1010. J/(kg K) and a porosity of 0.35.</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="15%" />
<col width="9%" />
<col width="21%" />
<col width="21%" />
<col width="18%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>rock</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>140</td>
<td>1</td>
<td><ol class="first last arabic simple" start="2563">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="1010">
<li></li>
</ol>
</td>
<td>0.35</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
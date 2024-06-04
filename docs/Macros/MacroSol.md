---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">sol</span></code><a class="headerlink" href="#sol" title="Permalink to this headline">¶</a></h1>
<ul class="simple">
<li>Group 1 - NTT, INTG</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="17%" />
<col width="10%" />
<col width="10%" />
<col width="63%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NTT</td>
<td>integer</td>
<td>1</td>
<td><div class="first last line-block">
<div class="line">Parameter that defines the type of solution required</div>
<div class="line">NTT &gt; 0 coupled solution</div>
<div class="line">NTT ≤ 0 heat transfer only solution</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>INTG</td>
<td>integer</td>
<td>-1</td>
<td><div class="first last line-block">
<div class="line">Parameter that defines element integration type</div>
<div class="line">INTG ≤ 0 Lobatto (node point) quadrature is used,
recommended for heat and mass problems without stress.</div>
<div class="line">INTG &gt; 0 Gauss quadrature is used, recommended for
problems requiring a stress solution.</div>
</div>
</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">sol</span></code>. In this example, a coupled heat-mass
solution using Lobatto quadrature is specified.</p>
<table border="1" class="docutils">
<colgroup>
<col width="50%" />
<col width="50%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>sol</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>-1</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
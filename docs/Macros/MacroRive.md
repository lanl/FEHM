---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">rive</span></code> or <code class="docutils literal notranslate"><span class="pre">well</span></code><a class="headerlink" href="#rive-or-well" title="Permalink to this headline">¶</a></h1>
<p>River or implicit well package.</p>
<ul>
<li><p class="first">Group 1 - KEYWORD</p>
<p>KEYWORD “wellmodel”</p>
<p>NRIVER, IRIVER</p>
<p>If IRIVER = 1</p>
<blockquote>
<div><p>INRIVERF, ISRIVERF, IFRIVERF,IWSP</p>
</div></blockquote>
<p>KEYWORD “<em>end macro</em>” (required)</p>
</li>
</ul>
<p>The input is terminated with keyword “<em>end rive</em>”, “<em>endrive</em>”, “<em>end well</em>” or
“<em>endwell</em>”, where the macro name should match the input macro name that was used.</p>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
<col width="9%" />
<col width="70%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KEYWORD “<em>wellmodel</em>”</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NRIVER</td>
<td>integer</td>
<td>Number of models defined for this call.</td>
</tr>
<tr class="row-odd"><td>IRIVER</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Type of surface/well flow</div>
<div class="line">IRIVER = 0 no routing, ponding, groundwater connection</div>
<div class="line">IRIVER = 1 simple fluid routing, no groundwater connection</div>
<div class="line">IRIVER = 2 simple stream definition, no groundwater connection</div>
<div class="line">IRIVER = 3 simple stream definition, with groundwater connection</div>
<div class="line">IRIVER = 4 simple stream definition, with groundwater connection (with ponding)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>INRIVERF</td>
<td>&#160;</td>
<td>Section number (id) of ith section.</td>
</tr>
<tr class="row-odd"><td>ISRIVERF</td>
<td>&#160;</td>
<td>Number of layers for the ith section.</td>
</tr>
<tr class="row-even"><td>IFRIVERF</td>
<td>&#160;</td>
<td>Number of coordinate points in the ith section.</td>
</tr>
<tr class="row-odd"><td>IWSP</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>KEYWORD “<em>endwell</em>”</td>
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
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">subm</span></code><a class="headerlink" href="#subm" title="Permalink to this headline">¶</a></h1>
<p>Create a new flow macro to represent boundary conditions on an extracted submodel.</p>
<p>Group 1 - KEYWORD, IZONE1, IZONE2</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="12%" />
<col width="74%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Keyword “flux”, “head”, or “pres” to specify type of boundary condition to output.</td>
</tr>
<tr class="row-odd"><td>IZONE1</td>
<td>integer</td>
<td>Zone defining submodel nodes.</td>
</tr>
<tr class="row-even"><td>IZONE2</td>
<td>integer</td>
<td>Zone defining nodes outside of the submodel (optional).</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
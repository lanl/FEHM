---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">bous</span></code><a class="headerlink" href="#bous" title="Permalink to this headline">¶</a></h1>
<p>Constant density and viscosity are used for the flow terms (<code class="docutils literal notranslate"><span class="pre">Boussinesq</span></code> approximation).</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Where the bous macro is used, the gravity term in the air phase is set to zero.</p>
</div>
<ul class="simple">
<li>Group 1 -     ICONS</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="17%" />
<col width="9%" />
<col width="74%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ICONS</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Parameter to enable constant density and viscosity for flow terms.</div>
<div class="line-block">
<div class="line">ICONS ≠ 0 enabled.</div>
<div class="line">ICONS = 0 disabled (default).</div>
</div>
</div>
</td>
</tr>
</tbody>
</table>
<p>The following is an example of bous. In this example the <code class="docutils literal notranslate"><span class="pre">Boussinesq</span></code> approximation is enabled.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">bous</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
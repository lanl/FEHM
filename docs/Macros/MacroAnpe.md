---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">anpe</span></code><a class="headerlink" href="#anpe" title="Permalink to this headline">¶</a></h1>
<p>Anisotropic permeability input. Adds cross terms to the perm macro.</p>
<p>The ANPE keyword implements a flux-continuous anisotropic permeability tensor with cross terms. The cross terms can either be input directly or grid rotation angles inputted and the cross terms calculated by FEHM. FEHM implements the method presented by Lee et al 2002.</p>
<p>The ANPE is incompatible with keywords GDKM, GDPM, DUAL, and DPDP.</p>
<ul>
<li><p class="first"><a href="#id1"><span class="problematic" id="id2">`</span></a>Group 1 - JA, JB, JC, ANXY, ANXZ, ANYZ</p>
<blockquote>
<div><ul class="simple">
<li>(JA, JB, JC - <a class="reference external" href="Macro20058.html">are defined here</a>)</li>
</ul>
</div></blockquote>
</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="19%" />
<col width="9%" />
<col width="11%" />
<col width="61%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ANXY</td>
<td>real</td>
<td>1.e-30</td>
<td>Anisotropic permeability in the xy-direction (m2).</td>
</tr>
<tr class="row-odd"><td>ANXZ</td>
<td>real</td>
<td>1.e-30</td>
<td>Anisotropic permeability in the xz-direction (m2).</td>
</tr>
<tr class="row-even"><td>ANYZ</td>
<td>real</td>
<td>1.e-30</td>
<td>Anisotropic permeability in the yz-direction (m2).</td>
</tr>
</tbody>
</table>
<ol class="arabic simple">
<li>Lee et al., 2002, Implementation of a Flux-Continuous Finite-Difference Method for Stratigraphic Hexahedron Grids, SPE Journal, Volume 7, Number 3, DOI 10.2118/80117-PA.</li>
</ol>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
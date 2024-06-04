---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">perm</span></code><a class="headerlink" href="#perm" title="Permalink to this headline">¶</a></h1>
<p>Assign permeabilities of the rock. Permeabilities represent average values of a volume associated with a node. Note that using rlp models 4 or 6 to describe relative permeabilities causes these values to be overwritten. Permeabilties may be entered as log values.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Required if macro <code class="docutils literal notranslate"><span class="pre">hyco</span></code> not used.</p>
</div>
<table border="1" class="docutils">
<colgroup>
<col width="22%" />
<col width="11%" />
<col width="13%" />
<col width="54%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>PNXD</td>
<td>real</td>
<td>1.e-30</td>
<td>Permeability in the x-direction (m2).</td>
</tr>
<tr class="row-odd"><td>PNYD</td>
<td>real</td>
<td>1.e-30</td>
<td>Permeability in the y-direction (m2).</td>
</tr>
<tr class="row-even"><td>PNZD</td>
<td>real</td>
<td>1.e-30</td>
<td>Permeability in the z-direction (m2).</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
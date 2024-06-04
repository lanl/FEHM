---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">pres</span></code><a class="headerlink" href="#pres" title="Permalink to this headline">¶</a></h1>
<p>The initial values defined in control statement <strong>pres</strong> supersede all others. Note that the term “saturated” refered to in IEOSD, is not the groundwater hydrology definition (volumetric fraction of pore void that is filled with water) used elsewhere in this document. Saturated here indicates that vapor and liquid phases exist simultaneously. The superheated region means that all pore space is filled with gas.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Required if macro <code class="docutils literal notranslate"><span class="pre">init</span></code> not used.</p>
</div>
<table border="1" class="docutils">
<colgroup>
<col width="6%" />
<col width="3%" />
<col width="3%" />
<col width="88%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>PHRD</td>
<td>real</td>
<td>PEIN</td>
<td>Initial pressure (MPa).</td>
</tr>
<tr class="row-odd"><td>TIND</td>
<td>real</td>
<td>&#160;</td>
<td>Initial temperature (oC) if IEOSD = 1 or 3, Initial saturation if IEOSD = 2</td>
</tr>
<tr class="row-even"><td>IEOSD</td>
<td>integer</td>
<td>1</td>
<td>Thermodynamic region parameter. IEOSD = 1, the compressed liquid regionIEOSD = 2, the saturation regionIEOSD = 3, the superheated region.If IEOSD &lt; 0 then the code uses ABS (IEOSD) and fixes the values of PHRD and TIND to the values provided above.</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
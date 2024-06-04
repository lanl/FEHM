---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">fper</span></code><a class="headerlink" href="#fper" title="Permalink to this headline">¶</a></h1>
<p>Assign permeability scaling factors.</p>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
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
<tr class="row-even"><td>SCALEX</td>
<td>real</td>
<td>1.0</td>
<td>Permeability scaling factor in the x-direction.</td>
</tr>
<tr class="row-odd"><td>SCALEY</td>
<td>real</td>
<td>1.0</td>
<td>Permeability scaling factor in the y-direction.</td>
</tr>
<tr class="row-even"><td>SCALEZ</td>
<td>real</td>
<td>1.0</td>
<td>Permeability scaling factor in the z-direction.</td>
</tr>
</tbody>
</table>
<p>The following is an example of fper. In this example, the values of the permeability (defined in a previous perm macro) are multiplied by 1.0 in the X direction, 0.5 in the Y direction, and 0.1 in the Z direction.</p>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
<col width="20%" />
<col width="10%" />
<col width="17%" />
<col width="17%" />
<col width="17%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">fper</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td>140</td>
<td>1</td>
<td>1.0</td>
<td>0.5</td>
<td>0.1</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
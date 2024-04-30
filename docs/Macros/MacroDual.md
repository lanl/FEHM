---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">dual</span></code><a class="headerlink" href="#dual" title="Permalink to this headline">¶</a></h1>
<p>Dual porosity formulation. There are three sets of parameter values at any nodal
position, for which property values must be defined. Nodes 1 to N (see macro <code class="docutils literal notranslate"><span class="pre">coor</span></code>
for definition of N) represent the fracture nodes, nodes N + 1 to 2N the first
matrix material, and nodes 2N + 1 to 3N the second matrix material.
When zones are used with the <code class="docutils literal notranslate"><span class="pre">dual</span></code> macro, additional zones are automatically
generated.</p>
<p>See instructions for the macro <code class="docutils literal notranslate"><span class="pre">zone</span></code> for a more detailed description.
The <code class="docutils literal notranslate"><span class="pre">dual</span></code> parameters are only defined for the first N nodes.</p>
<ul class="simple">
<li>Group 1 - IDUALP</li>
<li><a class="reference external" href="InputData.html#JA">Group 2 - JA, JB, JC, VOLFD1</a></li>
<li><a class="reference external" href="InputData.html#JA">Group 3 - JA, JB, JC, VOLFD2</a></li>
<li><a class="reference external" href="InputData.html#JA">Group 4 - JA, JB, JC, APUVD</a></li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="9%" />
<col width="9%" />
<col width="66%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>Input Variable</td>
<td>Format</td>
<td>Default</td>
<td>Description</td>
</tr>
<tr class="row-even"><td>IDUALP</td>
<td>integer</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">Solution descriptor for dual porosity solution.</div>
<div class="line">IDUALP = 0, information is read but not used.</div>
<div class="line">IDUALP ≠ 0, dual porosity solution is implemented</div>
<div class="line">For the special case of IDUALP = 2, the
permeabilities and conductivities are scaled
by the volume fraction, i.e., <span class="math notranslate nohighlight">\(k = vf * k\)</span>.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>VOLFD1</td>
<td>real</td>
<td>0.001</td>
<td>Volume fraction for fracture portion of the continuum.</td>
</tr>
<tr class="row-even"><td>VOLFD2</td>
<td>real</td>
<td>0.5</td>
<td>Volume fraction for the first matrix portion of the continuum.</td>
</tr>
<tr class="row-odd"><td>APUVD</td>
<td>real</td>
<td><ol class="first last arabic simple" start="5">
<li></li>
</ol>
</td>
<td>Length scale for the matrix nodes (m).</td>
</tr>
</tbody>
</table>
<p>The volume fractions VOLFD1 and VOLFD2 are related to the total volume by</p>
<p><span class="math notranslate nohighlight">\(VOLFD1 + VOLFD2 + VOLFD3 = 1.0\)</span></p>
<p>where VOLFD3 is the volume fraction of the second matrix node.
If permeability model IRLP = 4 is selected in control statement <code class="docutils literal notranslate"><span class="pre">rlp</span></code>,
VOLFD1 is calculated from RP15 (fracture porosity) in that control statement.</p>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">dual</span></code>. In this example, the dual porosity
solution is implemented for all nodes from 1 through 140. The volume fraction
for the fracture is 0.006711409, the volume fraction for the first matrix portion
is 0.335570470, and the length scale for the matrix nodes is 0.1 m.</p>
<table border="1" class="docutils">
<colgroup>
<col width="22%" />
<col width="19%" />
<col width="11%" />
<col width="48%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>dual</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>140</td>
<td>1</td>
<td>0.006711409</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>140</td>
<td>1</td>
<td>0.335570470</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>140</td>
<td>1</td>
<td>0.10</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
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
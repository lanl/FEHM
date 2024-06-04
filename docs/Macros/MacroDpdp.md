---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">dpdp</span></code><a class="headerlink" href="#dpdp" title="Permalink to this headline">¶</a></h1>
<p>Double porosity / double permeability formulation. There are two sets of parameter
values at any nodal position, for which property values must be defined.
Nodes 1 to N (see macro coor for definition of N) represent the fracture nodes
and nodes N + 1 to 2N the matrix material. When zones are used with the dpdp macro,
additional zones are automatically generated. See instructions for the macro zone
for a more detailed description. The dpdp parameters are only defined for the first N nodes.</p>
<ul class="simple">
<li>Group 1 - IDPDP</li>
<li>Group 2 - <a class="reference external" href="InputData.html#JA">JA, JB, JC, VOLFD1</a></li>
<li>Group 3 - <a class="reference external" href="InputData.html#JA">JA, JB, JC, APUV1</a></li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="5%" />
<col width="5%" />
<col width="82%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>Input Variable</td>
<td>Format</td>
<td>Default</td>
<td>Description</td>
</tr>
<tr class="row-even"><td>IDPDP</td>
<td>integer</td>
<td>&#160;</td>
<td>Solution descriptor for double porosity/double permeability solution. IDPDP = 0, information is read but not used. IDPDP ≠ 0, dpdp solution is implemented.</td>
</tr>
<tr class="row-odd"><td>VOLFD1</td>
<td>real</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>Volume fraction for fracture node.</td>
</tr>
<tr class="row-even"><td>APUV1</td>
<td>real</td>
<td><ol class="first last arabic simple" start="10">
<li></li>
</ol>
</td>
<td>Half spacing between fractures (m). See TRANSFLAG in macros MPTR and PTRK.</td>
</tr>
</tbody>
</table>
<p>The volume fraction VOLFD1 is related to the total volume by</p>
<p><span class="math notranslate nohighlight">\(VOLFD1 + VOLLFD2 = 1.0\)</span></p>
<p>where VOLFD2 is the volume fraction of the matrix node. If permeability model
<code class="docutils literal notranslate"><span class="pre">IRLP</span> <span class="pre">=</span> <span class="pre">4</span></code> is selected in control statement <strong>rlp</strong>,
VOLFD1 is calculated from RP15 (fracture porosity) in that control statement.</p>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">dpdp</span></code>. In this example,
the dual porosity/permeability solution is implemented for all nodes from
1 through 140. The fractional volume in the fractures (compared to the total volume)
is 0.005 and the length scale for matrix nodes is 0.1 m.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">dpdp</span>
<span class="mi">1</span>
<span class="mi">1</span> <span class="mi">140</span> <span class="mi">1</span> <span class="mf">0.005</span>

<span class="mi">1</span> <span class="mi">140</span> <span class="mi">1</span> <span class="mf">0.10</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
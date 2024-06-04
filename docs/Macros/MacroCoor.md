---
layout : page_macros
hero_height: is-hidden
---



<h1><code class="docutils literal notranslate"><span class="pre">coor</span></code><a class="headerlink" href="#coor" title="Permalink to this headline">¶</a></h1>
<p>Node coordinate data. These data are usually created by a mesh generation program, then cut and copied into the input file or a separate geometry data input file. The mesh must be a right handed coordinate system. Note that X, Y, and Z coordinates must be entered even if a problem is not three-dimensional. Version2.30 added the ability to provide the coordinate data in a formatted or unformatted file.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last"><code class="docutils literal notranslate"><span class="pre">coor</span></code> is required if macro <code class="docutils literal notranslate"><span class="pre">fdm</span></code> is not used.</p>
</div>
<ul class="simple">
<li>Group 1 - N</li>
<li>Group 2 - MB, CORD1, CORD2, CORD3</li>
</ul>
<p>To end the control section a blank line is entered.</p>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="4%" />
<col width="88%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>N</td>
<td>integer</td>
<td>Number of nodes in the grid</td>
</tr>
<tr class="row-odd"><td>MB</td>
<td>integer</td>
<td>Node number. If MB &lt; 0 then the difference between the absolute value of MB and the previously read absolute value of MB is used to generate intermediate values by interpolation.</td>
</tr>
<tr class="row-even"><td>CORD1</td>
<td>real</td>
<td>X-coordinate of node MB (m).</td>
</tr>
<tr class="row-odd"><td>CORD2</td>
<td>real</td>
<td>Y-coordinate of node MB (m).</td>
</tr>
<tr class="row-even"><td>CORD3</td>
<td>real</td>
<td>Z-coordinate of node MB (m).</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">coor</span></code>. In this example, there are 140 nodes in the grid. Node number 1 has X, Y, Z coordinates of 0., 200., and 0. meters respectively, node 2 has X, Y, Z coordinates of 12.5, 200., and 0. meters respectively, and so forth, with node number 140 having X, Y, Z coordinates of 300., 0., and 0. meters respectively.</p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="30%" />
<col width="30%" />
<col width="24%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>coor</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>140</td>
<td>0.00000</td>
<td>200.00000</td>
<td>0.00000</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>12.50000</td>
<td>200.00000</td>
<td>0.00000</td>
</tr>
<tr class="row-even"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-even"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>10</td>
<td>212.50000</td>
<td>200.00000</td>
<td>0.00000</td>
</tr>
<tr class="row-even"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-even"><td>.</td>
<td>.</td>
<td>.</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>140</td>
<td>300.00000</td>
<td>0.00000</td>
<td>0.00000</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
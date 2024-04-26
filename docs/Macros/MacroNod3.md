---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">nod3</span></code><a class="headerlink" href="#nod3" title="Permalink to this headline">¶</a></h1>
<p>Specify the node numbers for which detailed file output is desired and alternate nodes for terminal output.</p>
<ul class="simple">
<li>Group 1 -     M, M2, M3</li>
<li>Group 2 -     MN (1), MN (2), … , MN (M)</li>
<li>Group 3 -     X, Y, Z (as needed)</li>
<li>Group 4 -     MNI(1), MNI(2), … , MNI(M2)</li>
<li>Group 5 -     X, Y, Z (as needed)</li>
<li>Group 6 -     MNI(1), MNI(2), … , MNI(M3)</li>
<li>Group 7 -     X, Y, Z (as needed)</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="5%" />
<col width="3%" />
<col width="92%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>M</td>
<td>integer</td>
<td>Number of nodes for which information will be printed on the output file (iout). If M ≤ 0, pressure and temperature will be written on the output file for all nodes but no nodal parameter values will be printed in the history plot files. Group 2 is omitted if M ≤ 0.</td>
</tr>
<tr class="row-odd"><td>M2</td>
<td>integer</td>
<td>Number of nodes for short list (terminal printout). If M2 ≤ 0, Group 4 is omitted.</td>
</tr>
<tr class="row-even"><td>M3</td>
<td>&#160;</td>
<td>Number of nodes for short list (variable porosity model information printout). If M3 ≤ 0, Group 6 is omitted.</td>
</tr>
<tr class="row-odd"><td>MN</td>
<td>integer</td>
<td>M node numbers for which information will be printed on the output file (iout). If a MN(I) &lt; 0, then coordinates are used to define that print-out node, and the coordinate sets (X, Y, Z) for each MN(I) &lt; 0 are added after Group 2.</td>
</tr>
<tr class="row-even"><td>MNI</td>
<td>integer</td>
<td>M2 node numbers for which information will be printed on the terminal (short list). This group exists only if M2 ≠ 0. If MNI(I) &lt; 0, then coordinates are used to define the terminal output nodes, and the coordinate sets (X, Y, Z) for each MNI(I) &lt; 0 are added after Group 4.</td>
</tr>
<tr class="row-odd"><td>X</td>
<td>real</td>
<td>Coordinates of node for which information will be printed. One line for each MN or MNI &lt; 0. The code finds the node closest to the coordinate given. For 2-D problems set Z = 0. No input if no MN or MNI &lt; 0.</td>
</tr>
<tr class="row-even"><td>Y</td>
<td>real</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>Z</td>
<td>real</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>The following are examples of <code class="docutils literal notranslate"><span class="pre">nod3</span></code>. In the first example (top), detailed output to
the output file is specified for two nodes, the nodes numbered 50 and 88, and one
node is specified for terminal output, node 50. In the second example (bottom), two
nodes are specified for detailed output, the nodes numbered 50 and 88, and one node
is specified for terminal output, the node closest to the coordinates X = 100. m,
Y = 1000. m and Z = 0. m.</p>
<table border="1" class="docutils">
<colgroup>
<col width="60%" />
<col width="40%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>nod2</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>2</td>
<td>1</td>
</tr>
<tr class="row-odd"><td>50</td>
<td>88</td>
</tr>
<tr class="row-even"><td>50</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<table border="1" class="docutils">
<colgroup>
<col width="35%" />
<col width="41%" />
<col width="24%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>nod2</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>2</td>
<td>1</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>50</td>
<td>88</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>-88</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="100">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="1000">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
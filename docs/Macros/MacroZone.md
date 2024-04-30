---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">zone</span></code><a class="headerlink" href="#zone" title="Permalink to this headline">¶</a></h1>
<p>Geometric definition of grid for input parameter assignment, the default is input by nodes.</p>
<ul class="simple">
<li>Group 1 - IZONE</li>
<li>Group 2 - X1, X2, X3, X4 (for 2-D) or X1, X2, X3, X4, X5, X6, X7, X8 (for 3-D)</li>
<li>Group 3 - Y1, Y2, Y3, Y4 (for 2-D) or Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 (for 3-D)</li>
<li>Group 4 - Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 (for 3-D problems only)</li>
</ul>
<p>The following alternate form of input may be used (starting with Group 2):</p>
<ul class="simple">
<li>Group 2 - MACRO</li>
<li>Group 3 - XG, YG (for 2D) or XG, YG, ZG (for 3-D) [used with ‘list’ option]</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 3 - NIN, NODE(1), … , NODE(NIN) [used with ‘nnum’ option]</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 3 - TOL_ZONE, ZXY_MIN, ZXY_MAX [used with ‘xyli’ option]</li>
<li>Group 4 - XG, YG [used with ‘xyli’ option]</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="4%" />
<col width="3%" />
<col width="94%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IZONE</td>
<td>integer</td>
<td>Zone identification number for geometric input.</td>
</tr>
<tr class="row-odd"><td>X1-X8</td>
<td>real</td>
<td>X coordinates defining zone IZONE (m).</td>
</tr>
<tr class="row-even"><td>Y1-Y8</td>
<td>real</td>
<td>Y coordinates defining zone IZONE (m).</td>
</tr>
<tr class="row-odd"><td>Z1-Z8</td>
<td>real</td>
<td>Z coordinates defining zone IZONE (m).</td>
</tr>
<tr class="row-even"><td>MACRO</td>
<td>character*4</td>
<td>String denoting alternate input formatMACRO = “list”, read a list of X, Y, Z - coordinates, one set per line until a blank line is encountered. The nodes corresponding to these coordinates make up the zone. MACRO = “nnum”, read the number of specified nodes, followed by the node numbers. These comprise the zone.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>MACRO = “xyli”, read a column radius, followed by a list of X, Y - coordinates, one set per line until a blank line is encountered. The nodes contained in columns centered on each x, y pair and extending to the defined boundaries in the Z direction make up the zone. The column radius is necessary because there are (usually) slight variations in the Z direction of nodes above and below the prescribed X, Y coordinates.</td>
</tr>
<tr class="row-even"><td>XG</td>
<td>real</td>
<td>X coordinate of node to be included in IZONE (m).</td>
</tr>
<tr class="row-odd"><td>YG</td>
<td>real</td>
<td>Y coordinate of node to be included in IZONE (m).</td>
</tr>
<tr class="row-even"><td>ZG</td>
<td>real</td>
<td>Z coordinate of node to be included in IZONE (m).</td>
</tr>
<tr class="row-odd"><td>NIN</td>
<td>integer</td>
<td>Number of nodes in IZONE.</td>
</tr>
<tr class="row-even"><td>NODE(i)</td>
<td>integer</td>
<td>NIN node numbers of the nodes to be included in IZONE.</td>
</tr>
<tr class="row-odd"><td>TOL_ZONE</td>
<td>real</td>
<td>Column radius (m).</td>
</tr>
<tr class="row-even"><td>ZXY_MIN</td>
<td>real</td>
<td>Minimum Z coordinate for XY list (m).</td>
</tr>
<tr class="row-odd"><td>ZXY_MAX</td>
<td>real</td>
<td>Maximum Z coordinate for XY list (m).</td>
</tr>
</tbody>
</table>
<p>The geometric zone description is implemented by defining geometric regions. The coordinates given in Group 2, 3, and 4 refer to node positions. All properties defined by node (JA, JB, JC) in any control statements may be defined by <strong>zone</strong>. In the previous macro descriptions if JA &lt; 0, then the zone IZONE = ABS (JA) is referenced.</p>
<p>It is a good policy to refer to the input check file to insure that node assignments have been made as expected. When X, Y, Z coordinates are used to define zones, boundaries of those zones may be slightly different than specified. This is due to the inclusion of volume from elements adjoining included nodes.</p>
<p>When macro statements dpdp and dual are used, additional zone definitions are automatically generated. These are identified by zones 101-200 for the first set of matrix nodes and 201-300 for the second set of matrix nodes. For example, Zone 101 corresponds to the matrix material occupying the same space occupied by the fracture identified by Zone 1. Furthermore, Zone 201 refers to the second matrix layer in the dual control statement. Zones for the dpdp and dual matrix nodes may be explicitly defined and will not be superseded by automatically generated zones.</p>
<p>The macro <strong>zone</strong> must precede the usage of a ZONE reference. <strong>zone</strong> is ended with a blank line. <strong>zone</strong> can be called more than once and regions redefined. When this is done, all previous zone definitions are eliminated. A node may be included in only a single zone at a time.</p>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
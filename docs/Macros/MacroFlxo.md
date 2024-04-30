---
layout : page_macros
hero_height: is-hidden
---


<h1><code class="docutils literal notranslate"><span class="pre">flxo</span></code><a class="headerlink" href="#flxo" title="Permalink to this headline">¶</a></h1>
<p>Mass flow between two nodes is output by choosing this control statement.</p>
<ul class="simple">
<li>Group 1 -     NFLX</li>
<li>Group 2 -     IFLX1, IFLX2 (repeated NFLX times)</li>
<li>Group 3 -     X1, Y1, Z1 (as needed)</li>
<li>Group 4-      X2, Y2, Z2 (as needed)</li>
</ul>
<p>If IFLX1 &lt; 0, then after all IFLX1 and IFLX2 values are read, coordinates X1, Y1,
and Z1 are read and the node nearest to these coordinates is used. If IFLX2 &lt; 0,
coordinates for the second node are read in on another line. The code cycles
through each IFLX1 and IFLX2 in this manner, reading coordinates when needed.</p>
<p>Results are written to the screen if tty output is enabled and to the output
file <code class="docutils literal notranslate"><span class="pre">iout</span></code>.</p>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="5%" />
<col width="87%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>AIPED</td>
<td>real</td>
<td>Same as above for AIPED under keyword flow</td>
</tr>
<tr class="row-odd"><td>NFLX</td>
<td>integer</td>
<td>Number of internode mass flows to be calculated.</td>
</tr>
<tr class="row-even"><td>IFLX1</td>
<td>integer</td>
<td>First node to be used in mass flow calculation.</td>
</tr>
<tr class="row-odd"><td>IFLX2</td>
<td>integer</td>
<td>Second node to be used in mass flow calculation. If IFLX2 = 0, then the node connected to IFLX1 with the greatest internodal distance is used to calculate the mass flow.</td>
</tr>
<tr class="row-even"><td>X1</td>
<td>real</td>
<td>Coordinates of the first node to be used in mass flow calculation. Used only for those nodes where IFLX1 &lt; 0.</td>
</tr>
<tr class="row-odd"><td>Y1</td>
<td>real</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>Z1</td>
<td>real</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>X2</td>
<td>real</td>
<td>Coordinates of the second node to be used in mass flow calculation. Used only for those nodes where IFLX2 &lt; 0.</td>
</tr>
<tr class="row-even"><td>Y2</td>
<td>real</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>Z2</td>
<td>real</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>The following are examples of flxo. In these examples, one internode flux is calculated.</p>
<p>In the first case (top), from the node numbered 193 to the node numbered 195.</p>
<table border="1" class="docutils">
<colgroup>
<col width="55%" />
<col width="45%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>flxo</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>193</td>
<td>195</td>
</tr>
</tbody>
</table>
<p>In the second case, between nodes closest to coordinates 0., 0., 0. m and 20., 20., 20. m.</p>
<table border="1" class="docutils">
<colgroup>
<col width="38%" />
<col width="31%" />
<col width="31%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>flxo</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>-1</td>
<td>-7</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="20">
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
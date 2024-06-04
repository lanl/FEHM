---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">wgtu</span></code><a class="headerlink" href="#wgtu" title="Permalink to this headline">¶</a></h1>
<p>User defined destributed source boundary weight factors used in conjunction with macro boun.</p>
<ul class="simple">
<li>Group 1 - NAREA</li>
<li>Group 2 - IZONE_AREA(I), IAREAF(I), AREA01(I), AREA02(I), I = 1 to NAREA</li>
</ul>
<p>If <code class="docutils literal notranslate"><span class="pre">(IAREAF(I)</span> <span class="pre">=</span> <span class="pre">7)</span></code> a weight factor is read for each node in the current zone.</p>
<ul class="simple">
<li>Group 3 - II, (JJ, WGT_AREA(JJ), IJ = 1 to II)</li>
</ul>
<p>If <code class="docutils literal notranslate"><span class="pre">(IAREAF(I)</span> <span class="pre">=</span> <span class="pre">8)</span></code> a weight factor is read for each node in the current zone from the specified file. See Group 3 above for file format.</p>
<ul class="simple">
<li>Group 3 - WFILENAME</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="3%" />
<col width="3%" />
<col width="95%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NAREA</td>
<td>integer</td>
<td>Number of zones for which a weight factor will be calculated.</td>
</tr>
<tr class="row-odd"><td>IZONE_AREA</td>
<td>integer</td>
<td>Zone on which to calculate weight factors.</td>
</tr>
<tr class="row-even"><td>IAREAF</td>
<td>integer</td>
<td>Method for calculating weight factor for current zone.1 : Calculate area of each node in zone2 : Calculate area of each node based on total area given (area01) and portioned by volume size 3 : Calculate area of each node based on total area given (area01) and portioned by approximate area size 4 : Calculate weighting based on area (used in boun)5 : Calculate weighting based on area*perm (used in boun)6 : Calculate weighting based on vol*perm (used in boun)7 : Read a list of weights for every node in zone8 : Read a file with weights for every node in zone</td>
</tr>
<tr class="row-odd"><td>AREA01</td>
<td>real*8</td>
<td>Input area.</td>
</tr>
<tr class="row-even"><td>AREA02</td>
<td>real*8</td>
<td>Base impedance parameter.</td>
</tr>
<tr class="row-odd"><td>II</td>
<td>integer</td>
<td>Number of nodes in the current zone.</td>
</tr>
<tr class="row-even"><td>JJ</td>
<td>integer</td>
<td>Node number of the ijth node.</td>
</tr>
<tr class="row-odd"><td>WGT_AREA</td>
<td>real*8</td>
<td>Weighting factor for JJ</td>
</tr>
<tr class="row-even"><td>WFILENAME</td>
<td>character*100</td>
<td>Name of file contining weight factors for the current zone.</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
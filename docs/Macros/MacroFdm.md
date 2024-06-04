---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">fdm</span></code><a class="headerlink" href="#fdm" title="Permalink to this headline">¶</a></h1>
<p>Finite difference input.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last"><code class="docutils literal notranslate"><span class="pre">fdm</span></code> is required if macros <code class="docutils literal notranslate"><span class="pre">coor</span></code> are <code class="docutils literal notranslate"><span class="pre">elem</span></code> are not used.</p>
</div>
<ul class="simple">
<li>Group 1 - KEYWORD</li>
<li>Group 2 - MODCHAR (only if KEYWORD is “modf”)</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 2 - NX, NY, NZ (if KEYWORD is “block” or “poin”)</li>
<li>Group 3 - X0, Y0, Z0 (only if KEYWORD is ‘block)</li>
<li>Group 4 - MB, COORDINATE (X, Y or Z) (if KEYWORD is “poin”)</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 4 - MB, SPACING (DX, DY, DZ) (if KEYWORD is “bloc”)</li>
</ul>
<p>Group 4 is repeated for each NX, NY, and NZ, i.e., all X data are input, followed by Y data, followed by Z data for each division terminated by a blank line.</p>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="8%" />
<col width="84%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Keyword indicating format of finite difference input to follow (“block”, “poin”, or “modf”).</td>
</tr>
<tr class="row-odd"><td>MODCHAR</td>
<td>character*132</td>
<td>If the keyword is “modf”, the name of a modflow geometry data file is input and the finite difference input is read from that file and no other data is input.</td>
</tr>
<tr class="row-even"><td>NX</td>
<td>integer</td>
<td>Number of divisions in the x direction.</td>
</tr>
<tr class="row-odd"><td>NY</td>
<td>integer</td>
<td>Number of divisions in the y direction.</td>
</tr>
<tr class="row-even"><td>NZ</td>
<td>integer</td>
<td>Number of divisions in the z direction.</td>
</tr>
<tr class="row-odd"><td>X0</td>
<td>real</td>
<td>Coordinate of x origin point (m).</td>
</tr>
<tr class="row-even"><td>Y0</td>
<td>real</td>
<td>Coordinate of y origin point (m).</td>
</tr>
<tr class="row-odd"><td>Z0</td>
<td>real</td>
<td>Coordinate of z origin point (m).</td>
</tr>
<tr class="row-even"><td>X</td>
<td>real</td>
<td>X coordinate (m).</td>
</tr>
<tr class="row-odd"><td>Y</td>
<td>real</td>
<td>Y coordinate (m).</td>
</tr>
<tr class="row-even"><td>Z</td>
<td>real</td>
<td>Z coordinate (m).</td>
</tr>
<tr class="row-odd"><td>DX</td>
<td>real</td>
<td>Node spacing in the x direction (m).</td>
</tr>
<tr class="row-even"><td>DY</td>
<td>real</td>
<td>Node spacing in the y direction (m).</td>
</tr>
<tr class="row-odd"><td>DZ</td>
<td>real</td>
<td>Node spacing in the z direction (m).</td>
</tr>
<tr class="row-even"><td>MB</td>
<td>integer</td>
<td>Division number. If the division number is negative the code will space each divison from the previous to the current proportional to the assigned spacings.</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
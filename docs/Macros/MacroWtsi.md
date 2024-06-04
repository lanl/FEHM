---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">wtsi</span></code><a class="headerlink" href="#wtsi" title="Permalink to this headline">¶</a></h1>
<p>Water table, simplified.</p>
<p>Group 1 - NFREE, IZONE_FREE(I), I = 1 to NFREE, IFREE_IM_EX, HEAD_TOL, RLPTOL,</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="10%" />
<col width="71%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NFREE</td>
<td>integer</td>
<td>Number of zones to which water table condition are applied.</td>
</tr>
<tr class="row-odd"><td>IZONE_FREE</td>
<td>integer</td>
<td>Zone number.</td>
</tr>
<tr class="row-even"><td>IFREE_IM_EX</td>
<td>integer</td>
<td>Update parameter (0 - explicit update, 1 - implicit update).</td>
</tr>
<tr class="row-odd"><td>HEAD_TOL</td>
<td>real*8</td>
<td>Tolerance for head (m)</td>
</tr>
<tr class="row-even"><td>RLPTOL</td>
<td>real*8</td>
<td>Tolerance for saturation</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
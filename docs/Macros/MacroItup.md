---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">itup</span></code><a class="headerlink" href="#itup" title="Permalink to this headline">¶</a></h1>
<p>Controls upstream direction. The use of the itup macro is sometimes useful in
problems where the flow directions are changing rapidly. The parameter UPWGT
(in macro ctrl) must be greater than 0.5 for this macro to have any effect.</p>
<ul class="simple">
<li>Group 1 - IAD_UP</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="21%" />
<col width="12%" />
<col width="12%" />
<col width="56%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IAD_UP</td>
<td>integer</td>
<td>1000</td>
<td><div class="first last line-block">
<div class="line">Number of iterations after which the
upwind directions are held constant.</div>
<div class="line"><em>A value of 2 is suggested</em></div>
</div>
</td>
</tr>
</tbody>
</table>
<p>In the following example of itup, after 10 iterations the upwind directions are
held constant.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>itup</td>
</tr>
<tr class="row-even"><td>10</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
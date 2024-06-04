---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">impf</span></code><a class="headerlink" href="#impf" title="Permalink to this headline">¶</a></h1>
<p>Time step control based on maximum allowed variable change.</p>
<ul class="simple">
<li>Group 1 -     DELPT, DELTT, DELST, DELAT</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="7%" />
<col width="78%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>DELPT</td>
<td>real</td>
<td>Maximum allowable pressure change for which time step will be increased. (MPa)</td>
</tr>
<tr class="row-odd"><td>DELTT</td>
<td>real</td>
<td>Maximum allowable temperature change for which time step will be increased. (oC)</td>
</tr>
<tr class="row-even"><td>DELST</td>
<td>real</td>
<td>Maximum allowable saturation change for which time step will be increased.</td>
</tr>
<tr class="row-odd"><td>DELAT</td>
<td>real</td>
<td>Maximum allowable air pressure change for which time step will be increased. (MPa)</td>
</tr>
</tbody>
</table>
<p>The following is an examples of <code class="docutils literal notranslate"><span class="pre">impf</span></code>. In this example, pressure changes are
limited to 0.5 MPa, temperature changes to 20 oC, saturation changes to 0.1, and
air pressure changes to 0.05 MPa during a time step.</p>
<table border="1" class="docutils">
<colgroup>
<col width="26%" />
<col width="26%" />
<col width="22%" />
<col width="26%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>impf</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>0.5</td>
<td>20.0</td>
<td>0.1</td>
<td>0.05</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
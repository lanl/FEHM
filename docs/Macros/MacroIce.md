---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">ice</span></code> or <code class="docutils literal notranslate"><span class="pre">meth</span></code><a class="headerlink" href="#ice-or-meth" title="Permalink to this headline">¶</a></h1>
<p>Ice phase calculations, not tested.</p>
<ul class="simple">
<li>Group 1 -     ICE, SIIN, TMELT</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="11%" />
<col width="6%" />
<col width="83%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ICE</td>
<td>integer</td>
<td>Solution descriptor for ice solution.ICE = 0, information is read but not used.ICE ≠ 0, ice solution is implemented.</td>
</tr>
<tr class="row-odd"><td>SIIN</td>
<td>real</td>
<td>Default value for ice saturation (used when ice saturation SII in Group 2 is set to 0 at any node).</td>
</tr>
<tr class="row-even"><td>TMELT</td>
<td>real</td>
<td>Freezing temperature of water (oC).</td>
</tr>
<tr class="row-odd"><td>SII</td>
<td>real</td>
<td>Ice saturation. The default value is [0].</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
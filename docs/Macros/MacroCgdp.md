---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">cgdp</span></code><a class="headerlink" href="#cgdp" title="Permalink to this headline">¶</a></h1>
<p>Assign rate-limited gdpm nodes [i.e. make the connection value large for those nodes that serve only to link gdpm nodes (diffusion only) to the flow field].</p>
<p>Group 1 - wdd1</p>
<p>Group 2 - <a class="reference external" href="InputData.html#JA">JA, JB, JC, IGDPM_RATE_NODES</a></p>
<table border="1" class="docutils">
<colgroup>
<col width="24%" />
<col width="15%" />
<col width="61%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>macro</td>
<td>character</td>
<td>Process to be modified: ‘heat’, ‘tran’</td>
</tr>
<tr class="row-odd"><td>IGDPM_RATE_NODES</td>
<td>integer</td>
<td>GDPM nodes for which rate should be limited.</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
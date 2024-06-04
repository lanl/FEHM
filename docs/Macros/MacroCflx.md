---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">cflx</span></code><a class="headerlink" href="#cflx" title="Permalink to this headline">¶</a></h1>
<p>Total moles of liquid solute moving through a zone are output by choosing this control statement. Vapor solute molar flows are currently not available. When this macro is invoked, the following output is given at every solute time step:</p>
<ul class="simple">
<li>The sum of all solute source flow rates for each zone</li>
<li>The sum of all solute sink rates for each zone</li>
<li>The sum of all solute entering each zone</li>
<li>The sum of all solute leaving each zone</li>
<li>The net source/sink (boundary) solute flow for each zone</li>
</ul>
<p>The following values can be included on the macro line to specify which solute flows should be output:</p>
<ul class="simple">
<li>1 (source)</li>
<li>2 (sink)</li>
<li>3 (netin)</li>
<li>4 (netout)</li>
<li>5 (boundary)</li>
</ul>
<p>The default is to output all values.</p>
<p>Zones must be defined using macro zone prior to using this macro.</p>
<p>Group 1 - <code class="docutils literal notranslate"><span class="pre">CFLXZ</span></code></p>
<p>Group 2 - <code class="docutils literal notranslate"><span class="pre">ICFLXZ(I),</span> <span class="pre">I</span> <span class="pre">=</span> <span class="pre">1,</span> <span class="pre">CFLXZ</span></code></p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="9%" />
<col width="76%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>CFLXZ</td>
<td>integer</td>
<td>Number of zones for which output for solute flow through the zone iswritten.</td>
</tr>
<tr class="row-odd"><td>ICFLXZ</td>
<td>integer</td>
<td>Zone numbers for which solute flow output is written (NFLXZ zones)</td>
</tr>
</tbody>
</table>
<p>The following is an example of cflx. In this example solute flow through zones 1, 6 and 10 will be output.</p>
<table border="1" class="docutils">
<colgroup>
<col width="47%" />
<col width="16%" />
<col width="16%" />
<col width="21%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head" colspan="4">cflx</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>CFLXZ</td>
<td colspan="3">3</td>
</tr>
<tr class="row-odd"><td>ICFLXZ</td>
<td>1</td>
<td>6</td>
<td>10</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
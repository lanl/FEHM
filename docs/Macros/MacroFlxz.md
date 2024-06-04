---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">flxz</span></code><a class="headerlink" href="#flxz" title="Permalink to this headline">¶</a></h1>
<p>Total flow through a zone is output by choosing this control statement. When this macro is invoked, the following output is given at every heat and mass transfer time step:</p>
<ul class="simple">
<li>The sum of all source flow rates for each zone</li>
<li>The sum of all sink flow rates for each zone</li>
<li>The net quantity passing through each zone</li>
<li>The net source/sink quantity for each zone</li>
<li>The following keywords can be included on the macro line to specify which flows should be output: liquid mass, vapor mass, thermal energy. The default is to output any active quantity in the simulation.</li>
</ul>
<p>Zones must be defined using macro zone prior to using this macro.</p>
<ul class="simple">
<li>Group 1 -     NFLXZ</li>
<li>Group 2 -     IFLXZ(I), I = 1, NFLXZ</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="8%" />
<col width="77%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NFLXZ</td>
<td>integer</td>
<td>Number of zones for which output for water mass flow through the zone is written.</td>
</tr>
<tr class="row-odd"><td>IFLXZ</td>
<td>integer</td>
<td>Zone numbers for which water mass flow output is written (NFLXZ zones)</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">flxz</span></code>. In this example water mass flow through zones 1, 6 and 10 will be output.</p>
<table border="1" class="docutils">
<colgroup>
<col width="63%" />
<col width="16%" />
<col width="21%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head" colspan="3">flxz water</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>3</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
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
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">node</span></code><a class="headerlink" href="#node" title="Permalink to this headline">¶</a></h1>
<p>Specify the node numbers for which detailed output is desired. In version 2.30 macro node has been modified to allow multiple instances of the macro to be used (results are cumulative) and to allow the definition of “output zones” which are used in conjunction with macro hist. Only a single input format / keyword can be used for each instance of the node macro.</p>
<ul class="simple">
<li>Group 1 -     M</li>
<li>Group 2 -     MN (1), MN (2), … , MN (M)</li>
<li>Group 3 -     X, Y, Z (as needed)</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 1 -     KEYWORD</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="5%" />
<col width="4%" />
<col width="91%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>M</td>
<td>integer</td>
<td>Number of nodes for which information will be printed on the output (iout) and history plot (ishis, istrc) files. If M ≤ 0, pressure and temperature will be written on the output file for all nodes but no nodal parameter values will be printed in the history plot files. Group 2 is omitted if M ≤ 0.</td>
</tr>
<tr class="row-odd"><td>MN</td>
<td>integer</td>
<td>M node numbers for which information will be printed on the output file (iout). If MN(I) &lt; 0, then coordinates are used to define the print-out node, and the coordinate sets (X, Y, Z) for each MN(I) &lt; 0 are added after Group 2.</td>
</tr>
<tr class="row-even"><td>X</td>
<td>real</td>
<td>Coordinates of node for which information will be printed. One line for each MN &lt; 0. The code finds the node closest to the coordinate given. For 2-D problems set Z = 0. No input if MN &gt;0.</td>
</tr>
<tr class="row-odd"><td>Y</td>
<td>real</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>Z</td>
<td>real</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*5</td>
<td>Keyword ‘block’ to invoke node specification by JA, JB, JC format. Keyword ‘azone’ to invoke output zone specification by JA, JB, JC format. This keyword allows a single node to be included in multiple output zones.</td>
</tr>
</tbody>
</table>
<p>The following are examples of node. In the first example, 2 nodes are specified for output,
nodes 50 and 88.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">node</span>
    <span class="mi">2</span>
    <span class="mi">50</span>    <span class="mi">88</span>
</pre></div>
</div>
<p>In the second example, two nodes are specified for output, the node numbered
50 and the node closest to the coordinates X = 100. m, Y = 1000. m
and Z = 0. m.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">node</span>
    <span class="mi">2</span>
    <span class="mi">50</span>    <span class="o">-</span><span class="mi">88</span>

    <span class="mf">100.</span>  <span class="mf">1000.</span>    <span class="mf">0.</span>
</pre></div>
</div>
<p>In the third example, output is specified for the block of
nodes 1, 11, 21, 31, 41, 51, 61, 71, 81, 91 and for those nodes defined by
zone 3 (see macro zone).</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">node</span>
<span class="n">block</span>
    <span class="mi">1</span>     <span class="mi">100</span>      <span class="mi">10</span>
    <span class="o">-</span><span class="mi">3</span>      <span class="mi">0</span>       <span class="mi">0</span>
</pre></div>
</div>
<p>In the fourth example, output is specified for
two zones, the first zone contains nodes 1, 11, 21, 31, 41, 51, 61, 71, 81, 91
and the second zone is made up of the nodes defined for zone 3
(previously specified using the zone macro).</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">node</span>
<span class="n">azone</span>
   <span class="mi">1</span>      <span class="mi">100</span>      <span class="mi">10</span>
   <span class="o">-</span><span class="mi">3</span>       <span class="mi">0</span>       <span class="mi">0</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
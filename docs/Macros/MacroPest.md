---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">pest</span></code><a class="headerlink" href="#pest" title="Permalink to this headline">¶</a></h1>
<p>Output variable information for PEST parameter estimation routine.</p>
<ul class="simple">
<li>Group 1 -     MPEST</li>
<li>Group 2 -     NPEST(I), I = 1, MPEST</li>
<li>Group 3 -     X, Y, Z (as needed)</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="5%" />
<col width="87%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>MPEST</td>
<td>integer</td>
<td>Number of nodes for PEST output. At present the code outputs only pressures (heads), saturations, and temperatures.</td>
</tr>
<tr class="row-odd"><td>NPEST(I)</td>
<td>integer</td>
<td>Node numbers printed to the output file (fehmn.pest) with values of variables listed above. If NPEST(I) &lt; 0 then the node numbers are determined from the coordinates.</td>
</tr>
<tr class="row-even"><td>X, Y, Z</td>
<td>real</td>
<td>Coordinates in grid if NPEST(I) &lt; 0. The coordinates are used to find the node closest in distance to that point and that node is substituted for NPEST(I).</td>
</tr>
</tbody>
</table>
<p>The following is an example of pest. In this example pest output is specified at
5 nodes, nodes numbered 21, 23, 35, and 47, and the node closest to the coordinates
X=10. m, Y=15. m, Z=20. m.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">pest</span>
    <span class="mi">21</span>  <span class="mi">23</span>  <span class="mi">35</span> <span class="mi">47</span> <span class="o">-</span><span class="mi">50</span>
    <span class="mf">10.</span> <span class="mf">15.</span> <span class="mf">20.</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
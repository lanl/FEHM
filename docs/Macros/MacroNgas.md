---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">ngas</span></code><a class="headerlink" href="#ngas" title="Permalink to this headline">¶</a></h1>
<p>Noncondensible gas transport.</p>
<ul class="simple">
<li>Group 1 -     ICO2D</li>
<li>Group 2 - JA, JB, JC, PCO2</li>
<li>Group 3 - JA, JB, JC, CPNK</li>
<li>Group 4 - JA, JB, JC, QCD</li>
</ul>
<p>Note that all Group 2 values are entered first, followed by Group 3 values, followed by Group 4 values.</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="8%" />
<col width="8%" />
<col width="69%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ICO2D</td>
<td>integer</td>
<td>3</td>
<td><div class="first last line-block">
<div class="line">Solution descriptor for noncondensible gas transport.</div>
<div class="line">ICO2D = 1, the 3 degree of freedom solution will be
reduced to a 1 degree of freedom problem. (See macro
<code class="docutils literal notranslate"><span class="pre">iter</span></code>, the parameter ICOUPL is also set to 5 if
ICO2D = 1.)</div>
<div class="line">ICO2D = 2, the 3 degree of freedom solution will be
reduced to a 2 degree of freedom problem. (See
macro <code class="docutils literal notranslate"><span class="pre">iter</span></code>, the parameter ICOUPL is also set to
5 if ICO2D = 2.) ICO2D = 3, full 3 degree of freedom.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>PCO2</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">Initial partial pressure of noncondensible gas. If
PCO2 &lt; 0 then ABS (PCO2) is interpreted as a temperature
and the partial pressure of the noncondensible gas is
calculated according to the formula:</div>
<div class="line"><span class="math notranslate nohighlight">\(PCO2 = P_T - P_{SAT}(T)\)</span></div>
<div class="line">where <span class="math notranslate nohighlight">\(P_T\)</span> is the total pressure and <span class="math notranslate nohighlight">\(P_{SAT}(T)\)</span>
is the water saturation pressure and is a function of temperature only.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>CPNK</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">If CPNK ≤ 0, then ABS (CPNK) is the specified noncondensible
pressure and will be held at that value.</div>
<div class="line">If CPNK &gt; 0, then CPNK is the specified relative humidity
and the saturation, <span class="math notranslate nohighlight">\(S_l\)</span>, is calculated using the
vapor pressure lowering formula and the capillary pressure formula:</div>
<div class="line"><span class="math notranslate nohighlight">\(Pcap(S_l) = \mathrm{ln}(h)\rho_l RT\)</span></div>
<div class="line">where <span class="math notranslate nohighlight">\(Pcap\)</span> is the capillary function, <span class="math notranslate nohighlight">\(h\)</span> is
the humidity, <span class="math notranslate nohighlight">\(R\)</span> is the gas constant, <span class="math notranslate nohighlight">\(T\)</span>
is the temperature, and <span class="math notranslate nohighlight">\(\rho_l\)</span> is the liquid
density. Once the formula is solved, <span class="math notranslate nohighlight">\(S_l\)</span> is held
constant. The humidity condition is only enabled for the
van Genuchten capillary function model. See macro <code class="docutils literal notranslate"><span class="pre">rlp</span></code>.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>QCD</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Specified air source strength (kg/sec).</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">ngas</span></code>. In this example, a full 3 degrees of
freedom solution is specified. The initial temperature at nodes 1 to 800 is 20
oC and the code is asked to calculate the initial noncondensible gas pressure.
There is no specified noncondensible gas source.</p>
<table border="1" class="docutils">
<colgroup>
<col width="32%" />
<col width="26%" />
<col width="16%" />
<col width="26%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>ngas</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>3</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>800</td>
<td>1</td>
<td>-20</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>800</td>
<td>1</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>800</td>
<td>1</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">chea</span></code><a class="headerlink" href="#chea" title="Permalink to this headline">¶</a></h1>
<p>Convert output from pressure to head (non-head problems). The reference temperature and pressure are used to calculate density for the head calculation. For this macro the data are entered on the macro line, and if omitted the specified default values are used. (Note that all five values must be entered to override the default values.)</p>
<table border="1" class="docutils">
<colgroup>
<col width="19%" />
<col width="10%" />
<col width="11%" />
<col width="60%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>HEAD0</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Reference head (m)</td>
</tr>
<tr class="row-odd"><td>TEMP0</td>
<td>real</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>Reference temperature (oC)</td>
</tr>
<tr class="row-even"><td>PRES0</td>
<td>real</td>
<td>0.1</td>
<td>Reference pressure (MPa)</td>
</tr>
<tr class="row-odd"><td>SAT_ICH</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Saturation adjustment after variable switch</td>
</tr>
<tr class="row-even"><td>HEAD_ID</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Output head identification for small saturations</td>
</tr>
</tbody>
</table>
<p>The following is an example of chea. In this example pressures will be converted to heads for output using a reference pressure of 1 MPa and a reference temperature of 25 C.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">chea</span> <span class="mf">0.</span> <span class="mf">25.</span> <span class="mf">1.</span> <span class="mf">0.</span> <span class="mf">0.</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
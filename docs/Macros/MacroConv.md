---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">conv</span></code><a class="headerlink" href="#conv" title="Permalink to this headline">¶</a></h1>
<p>Convert input from head to pressure. Often used when converting a head-based isothermal model to a heat and mass simulation with pressure and temperature variables. The reference temperature and head are used to calculate density for the head calculation. It should be noted that this is an approximate method. Since the density is a nonlinear function of pressure and temperature this method will give slightly different answers than a calculation allowing a water column to come to thermal and mechanical equilibrium.</p>
<p>The reference head (<code class="docutils literal notranslate"><span class="pre">head0</span></code>) is converted to a pressure and added to the reference pressure (<code class="docutils literal notranslate"><span class="pre">conv1</span></code>) and this sum is used with the reference temperature to calculate a density. This density is used to convert the head to pressure for the identified zone. The option of adding a temperature gradient is provided as well.</p>
<p>The reference head (<code class="docutils literal notranslate"><span class="pre">head0</span></code>) is entered on the macro line:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">conv</span> <span class="n">HEAD0</span>
</pre></div>
</div>
<ul class="simple">
<li>Group 1 - NCONV</li>
<li>Group 2 – ZONE_CONV, ICONVF, CONV1, CONV2, CORDC, IDIRC, VARC</li>
</ul>
<p>Group 2 is entered for each zone (nconv times).</p>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
<col width="11%" />
<col width="70%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>HEAD0</td>
<td>real</td>
<td>Reference head (m)</td>
</tr>
<tr class="row-odd"><td>NCONV</td>
<td>integer</td>
<td>Number of zones for variable conversion</td>
</tr>
<tr class="row-even"><td>ZONE_CONV</td>
<td>integer</td>
<td>Zone for variable conversion</td>
</tr>
<tr class="row-odd"><td>ICONVF</td>
<td>integer</td>
<td>iconvf : 1- initial conditions, 2-boundary (fixed head)</td>
</tr>
<tr class="row-even"><td>CONV1</td>
<td>real</td>
<td>Reference pressure (MPa)</td>
</tr>
<tr class="row-odd"><td>CONV2</td>
<td>real</td>
<td>Reference temperature (oC)</td>
</tr>
<tr class="row-even"><td>CORDC</td>
<td>real</td>
<td>Reference coordinate (m)</td>
</tr>
<tr class="row-odd"><td>IDIRC</td>
<td>integer</td>
<td>Coordinate direction (1 - X, 2 - Y, 3 - Z)</td>
</tr>
<tr class="row-even"><td>VARC</td>
<td>real</td>
<td>Temperature gradient (oC/m)</td>
</tr>
</tbody>
</table>
<p>The following is an example of conv. Here the density used to convert head to pressure is calculated with a reference head of 1000 m plus 0.1 MPa and 80 °C. The nodes in zone 45 are converted from heads to pressures at 80 °C.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">conv</span>
<span class="mf">1000.</span>
<span class="mi">1</span>
<span class="mi">45</span>
<span class="mi">1</span>
<span class="mf">0.4</span>
<span class="mf">80.</span>
<span class="mf">0.</span>
<span class="mi">0</span>
<span class="mf">0.0</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
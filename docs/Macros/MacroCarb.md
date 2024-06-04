---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">carb</span></code><a class="headerlink" href="#carb" title="Permalink to this headline">¶</a></h1>
<p>Macro carb is used to set up a <span class="math notranslate nohighlight">\(CO_2\)</span> problem. Input following the problem type is grouped using sub keywords.</p>
<p>Group 1 - IPRTYPE</p>
<p>Group 1 -       KEYWORD</p>
<p>KEYWORD <code class="docutils literal notranslate"><span class="pre">co2pres</span></code></p>
<ul class="simple">
<li><a class="reference external" href="InputData.html#JA">JA, JB, JC, PHICO2, TCO2, ICES</a></li>
</ul>
<p>KEYWORD <code class="docutils literal notranslate"><span class="pre">co2flow</span></code></p>
<ul class="simple">
<li><a class="reference external" href="InputData.html#JA">JA, JB, JC, SKTMP, ESKTMP, AIPED, IFLG_FLOWMAC</a></li>
</ul>
<p>KEYWORD <code class="docutils literal notranslate"><span class="pre">co2diff</span></code></p>
<ul class="simple">
<li><a class="reference external" href="InputData.html#JA">JA, JB, JC, DIFF, TORTCO2</a></li>
</ul>
<p>KEYWORD <code class="docutils literal notranslate"><span class="pre">co2frac</span></code></p>
<ul class="simple">
<li><a class="reference external" href="InputData.html#JA">JA, JB, JC, FW, FL, YC, CSALT, INICO2FLG</a></li>
</ul>
<p>KEYWORD <code class="docutils literal notranslate"><span class="pre">userprop</span></code></p>
<ul class="simple">
<li>DENC, DENCP, DENCT, ENC, ENCP, ENCT, VISC, VISCP, VISCT</li>
<li>DENW, DENWP, DENWT, ENW, ENWP, ENWT, VISW, VISWP, VISWT</li>
</ul>
<p>KEYWORD <code class="docutils literal notranslate"><span class="pre">brine</span></code></p>
<p>Input is terminated with KEYWORD <code class="docutils literal notranslate"><span class="pre">end</span> <span class="pre">carb</span></code> or <code class="docutils literal notranslate"><span class="pre">endcarb</span></code>.</p>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="31%" />
<col width="2%" />
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
<tr class="row-even"><td>IPRTYPE</td>
<td>integer</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">1 = Water only problem (2 DOFs)</div>
<div class="line">2 = <span class="math notranslate nohighlight">\(CO_2\)</span> only problem (2 DOFs)</div>
<div class="line">3 = <span class="math notranslate nohighlight">\(CO_2\)</span>-water problem, no solubility (3 DOFs)</div>
<div class="line">4 = <span class="math notranslate nohighlight">\(CO_2\)</span>-water problem, with solubility (4 DOFs)</div>
<div class="line">5 = <span class="math notranslate nohighlight">\(CO_2\)</span>-water-air with solubility (5 DOFs)</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>&#160;</td>
<td>Remaining input is grouped using sub-macro keywords.</td>
</tr>
<tr class="row-even"><td>KEYWORD “end carb” or “endcarb”</td>
<td colspan="3">End of carb input.</td>
</tr>
<tr class="row-odd"><td>KEYWORD “co2pres”</td>
<td colspan="3">Set up the initial pressure (uses the same format as the pres macro)</td>
</tr>
<tr class="row-even"><td>PHICO2</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Initial <span class="math notranslate nohighlight">\(CO_2\)</span> pressure (MPa).</td>
</tr>
<tr class="row-odd"><td>TCO2</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Initial <span class="math notranslate nohighlight">\(CO_2\)</span> temperature (oC)</td>
</tr>
<tr class="row-even"><td>ICES</td>
<td>integer</td>
<td>0</td>
<td><div class="first last line-block">
<div class="line">Initial guess for phase state of <span class="math notranslate nohighlight">\(CO_2\)</span> (actual phase will be calculated internally), ICES settings are:</div>
<div class="line">1 = liquid</div>
<div class="line">2 =  two-phase liquid and vapor</div>
<div class="line">3 = vapor</div>
<div class="line">4 = super-critical <span class="math notranslate nohighlight">\(CO_2\)</span>.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>KEYWORD “co2flow”</td>
<td colspan="3">Set up co2 flow boundary conditions (similar to the flow macro used to set up water boundary conditions)</td>
</tr>
<tr class="row-even"><td>SKTMP</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">For IFLG_FLOWMAC = 1, 2, 3, 4, 5 or 9 <span class="math notranslate nohighlight">\(CO_2\)</span> flowing pressure (MPa).</div>
<div class="line">For SKTMP = 0 the initial value of pressure will be used for the flowing pressure.</div>
<div class="line">For IFLG_FLOWMAC = 6 or 7 water mass flow rate (kg/s).</div>
<div class="line">For IFLG_FLOWMAC = 8 <span class="math notranslate nohighlight">\(CO_2\)</span> flowing saturation.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>ESKTMP</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">For IFLG_FLOWMAC = 1, 2, 3, 6, 7 or 9 Enthalpy of fluid injected (MJ/kg). If the fluid is flowing from the rock mass, then the in-place enthalpy is used. If EFLOW&lt;0, then ABS(EFLOW) is interpreted as a temperature (C) and the enthalpy calculated accordingly.</div>
<div class="line">For IFLG_FLOWMAC = 4 or 5 <span class="math notranslate nohighlight">\(CO_2\)</span> flowing saturation.</div>
<div class="line">For IFLG_FLOWMAC = 8 mass fraction of CO2</div>
</div>
</td>
</tr>
<tr class="row-even"><td>AIPED</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">For IFLG_FLOWMAC = 1, 2, 4 or 9 <span class="math notranslate nohighlight">\(CO_2\)</span> impedance parameter.</div>
<div class="line">For IFLG_FLOWMAC = 5 or 6 value is ignored.</div>
<div class="line">For IFLG_FLOWMAC = 7 <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction in water</div>
<div class="line">For IFLG_FLOWMAC = 8 Water mass flow rate (kg/s).</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>IFLG_FLOWMAC</td>
<td>integer</td>
<td>0</td>
<td><div class="first last line-block">
<div class="line">Flag specifying boundary condition type, IFLG_FLOWMAC values:</div>
<div class="line">1 = Constant pressure boundary condition with inflow or outflow allowed. AIPED is user specified</div>
<div class="line">2 = Constant pressure boundary condition with only outflow allowed. AIPED is user specified</div>
<div class="line">3 = Constant pressure boundary condition. AIPED is calculated in the code based on block geometric parameters.</div>
<div class="line">4 = Constant pressure and constant saturation boundary condition. AIPED is user specified</div>
<div class="line">5 = Constant pressure and constant saturation boundary condition. AIPED is calculated in the code based on block geometric parameters.</div>
<div class="line">6 = Constant free phase <span class="math notranslate nohighlight">\(CO_2\)</span> mass flow rate boundary condition.</div>
<div class="line">7 = Constant source of water with specified mass fraction of CO2 (kg/s)</div>
<div class="line">8 = ????</div>
<div class="line">9 = Partial explicit update of nonlinear part of <span class="math notranslate nohighlight">\(CO_2\)</span> constant pressure</div>
</div>
</td>
</tr>
<tr class="row-even"><td>KEYWORD “co2diff”</td>
<td colspan="3">Read <span class="math notranslate nohighlight">\(CO_2\)</span> diffusivity in water</td>
</tr>
<tr class="row-odd"><td>DIFF</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Diffusion</td>
</tr>
<tr class="row-even"><td>TORTCO2</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Tortuosity for <span class="math notranslate nohighlight">\(CO_2\)</span>-water vapor diffusion.</td>
</tr>
<tr class="row-odd"><td>KEYWORD “co2frac”</td>
<td colspan="3">Read initial <span class="math notranslate nohighlight">\(CO_2\)</span>, air, water/brine saturation. FG, CO2/air-rich gas saturation (volume fraction), <span class="math notranslate nohighlight">\(FG = 1 - FW - FL.\)</span></td>
</tr>
<tr class="row-even"><td>FW</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Water-rich liquid saturation (volume fraction).</td>
</tr>
<tr class="row-odd"><td>FL</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><span class="math notranslate nohighlight">\(CO_2\)</span>-rich super-critical/liquid phase saturation (volume fraction).</td>
</tr>
<tr class="row-even"><td>YC</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Mass fraction of <span class="math notranslate nohighlight">\(CO_2\)</span> in the <span class="math notranslate nohighlight">\(CO_2\)</span>-rich phase.</td>
</tr>
<tr class="row-odd"><td>CSALT</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Initial salt concentration in water for brine (ppm) (only used if “brine” keyword is invoked.)</td>
</tr>
<tr class="row-even"><td>INICO2FLG</td>
<td>integer</td>
<td>0</td>
<td>Flag to override <span class="math notranslate nohighlight">\(CO_2\)</span> fractions read from restart file. If set to 1 the input values are used instead of those read from the restart file.</td>
</tr>
<tr class="row-odd"><td>KEYWORD “userprop”</td>
<td colspan="3">Read user defined properties for <span class="math notranslate nohighlight">\(CO_2\)</span> and brine</td>
</tr>
<tr class="row-even"><td>DENC</td>
<td>real</td>
<td>&#160;</td>
<td><span class="math notranslate nohighlight">\(CO_2\)</span> density (<span class="math notranslate nohighlight">\(kg/m^3\)</span>)</td>
</tr>
<tr class="row-odd"><td>DENCP</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of density with respect to pressure.</td>
</tr>
<tr class="row-even"><td>DENCT</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of density with respect to temperature.</td>
</tr>
<tr class="row-odd"><td>ENC</td>
<td>real</td>
<td>&#160;</td>
<td><span class="math notranslate nohighlight">\(CO_2\)</span> enthalpy (<span class="math notranslate nohighlight">\(MJ/kg\)</span>).</td>
</tr>
<tr class="row-even"><td>ENCP</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of enthalpy with respect to pressure.</td>
</tr>
<tr class="row-odd"><td>ENCT</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of enthalpy with respect to temperature.</td>
</tr>
<tr class="row-even"><td>VISC</td>
<td>real</td>
<td>&#160;</td>
<td><span class="math notranslate nohighlight">\(CO_2\)</span> viscosity (<span class="math notranslate nohighlight">\(Pa \cdot s\)</span>)</td>
</tr>
<tr class="row-odd"><td>VISCP</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of viscosity with respect to pressure.</td>
</tr>
<tr class="row-even"><td>VISCT</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of viscosity with respect to temperature.</td>
</tr>
<tr class="row-odd"><td>DENW</td>
<td>real</td>
<td>&#160;</td>
<td>Brine density (<span class="math notranslate nohighlight">\(kg/m^3\)</span>)</td>
</tr>
<tr class="row-even"><td>DENWP</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of density with respect to pressure.</td>
</tr>
<tr class="row-odd"><td>DENWT</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of density with respect to temperature.</td>
</tr>
<tr class="row-even"><td>ENW</td>
<td>real</td>
<td>&#160;</td>
<td>Brine enthalpy (<span class="math notranslate nohighlight">\(MJ/kg\)</span>).</td>
</tr>
<tr class="row-odd"><td>ENWP</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of enthalpy with respect to pressure.</td>
</tr>
<tr class="row-even"><td>ENWT</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of enthalpy with respect to temperature.</td>
</tr>
<tr class="row-odd"><td>VISW</td>
<td>real</td>
<td>&#160;</td>
<td>Brine viscosity (<span class="math notranslate nohighlight">\(Pa \cdot s\)</span>)</td>
</tr>
<tr class="row-even"><td>VISWP</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of viscosity with respect to pressure.</td>
</tr>
<tr class="row-odd"><td>VISWT</td>
<td>real</td>
<td>&#160;</td>
<td>Derivative of viscosity with respect to temperature.</td>
</tr>
<tr class="row-even"><td>KEYWORD “brine”</td>
<td colspan="3">Invoke option for brine in the simulation. (salt-concentration dependent <span class="math notranslate nohighlight">\(CO_2\)</span> solubility)</td>
</tr>
</tbody>
</table>
<p>In the following example, zone 1 is injecting <span class="math notranslate nohighlight">\(CO_2\)</span> dissolved water at 0.001 kg/s. The temperature is 20oC. The water has a dissolved <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction of 0.3. The code will check internally whether the user specified mass fraction exceeds the equilibrium mass fraction calculated using the pressure and temperature values of the injection node. In case it does exceed that value, it is fixed at the equilibrium mass fraction. The user can specify a value of “zero” and the code will automatically fix the dissolved <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction at the equilibrium value. Zone 2 is maintained at initial pressure using “aiped” calculated internally.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">carb</span>
    <span class="mi">4</span>

<span class="n">co2pres</span>
     <span class="mi">1</span>     <span class="mi">0</span>     <span class="mi">0</span>     <span class="mf">3.</span>     <span class="mf">20.</span>     <span class="mi">4</span>
    <span class="o">-</span><span class="mi">1</span>     <span class="mi">0</span>     <span class="mi">0</span>     <span class="mi">13</span>     <span class="mf">20.</span>     <span class="mi">4</span>
    <span class="o">-</span><span class="mi">2</span>     <span class="mi">0</span>     <span class="mi">0</span>     <span class="o">.</span><span class="mi">6</span>     <span class="mf">20.</span>     <span class="mi">4</span>

<span class="n">co2frac</span>
     <span class="mi">1</span>   <span class="mi">0</span>   <span class="mi">0</span>   <span class="mf">1.0</span>         <span class="mf">0.0</span>       <span class="mi">0</span>   <span class="mi">100000</span>    <span class="mi">0</span>
    <span class="o">-</span><span class="mi">1</span>   <span class="mi">0</span>   <span class="mi">0</span>    <span class="mf">0.9465</span>   <span class="o">.</span><span class="mi">0535</span>       <span class="mi">0</span>   <span class="mf">0.</span>         <span class="mf">0.</span>

<span class="n">co2flow</span>
    <span class="o">-</span><span class="mi">2</span> <span class="mi">0</span> <span class="mi">0</span>   <span class="mi">0</span>        <span class="o">-</span><span class="mf">20.</span>   <span class="o">-</span><span class="mf">1.e-1</span>    <span class="mi">3</span>
    <span class="o">-</span><span class="mi">1</span> <span class="mi">0</span> <span class="mi">0</span> <span class="o">-</span><span class="mf">0.0001</span>    <span class="o">-</span><span class="mf">20.</span>     <span class="mf">0.</span>      <span class="mi">6</span>

<span class="n">end</span> <span class="n">carb</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
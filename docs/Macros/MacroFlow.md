---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">flow</span></code><a class="headerlink" href="#flow" title="Permalink to this headline">¶</a></h1>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Either <code class="docutils literal notranslate"><span class="pre">boun</span></code> or <code class="docutils literal notranslate"><span class="pre">flow</span></code> are required for a flow problem.</p>
</div>
<p>Flow data. Source and sink parameters are input and may be used to apply boundary conditions. Note that the alternative definitions for isothermal models apply when flow is used in conjunction with control statement <a href="#id1"><span class="problematic" id="id2">**</span></a>airwater <a href="#id3"><span class="problematic" id="id4">**</span></a>(<a class="reference external" href="Macro80880.html">See Control statement airwater or air (optional)</a>).</p>
<table border="1" class="docutils">
<colgroup>
<col width="6%" />
<col width="2%" />
<col width="2%" />
<col width="90%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>Non-Isothermal model</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>SKD</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Heat and mass source strength (kg/s), heat only (MJ/s). Negative value indicates injection into the rock mass.</td>
</tr>
<tr class="row-even"><td>EFLOW</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Enthalpy of fluid injected (MJ/kg). If the fluid is flowing from the rock mass, then the in-place enthalpy is used.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>If EFLOW &lt; 0, then ABS(EFLOW) is interpreted as a temperature (oC) and the enthalpy (assuming water only) calculated accordingly. In heat only problems with EFLOW &lt; 0, the node is in contact with a large heat pipe that supplies heat to the node through an impedance AIPED so as to maintain its temperature near ABS (EFLOW). Large values (approximately 1000) of AIPED are recommended.</td>
</tr>
<tr class="row-even"><td>AIPED</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Impedance parameter. If AIPED is nonzero, the code interprets SKD as a flowing wellbore pressure (MPa) with an impedance ABS(AIPED). If AIPED &lt; 0, flow is only allowed out of the well. For heat only, AIPED is the thermal resistance. If AIPED = 0, SKD is flow rate. If AIPED ≠ 0 and SKD = 0 the initial value of pressure will be used for the flowing pressure.</td>
</tr>
</tbody>
</table>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">If the porosity of the node is zero, then there is only a temperature solution,
and the code forms a source proportional to the enthalpy difference.
The source term is given by <span class="math notranslate nohighlight">\(Q = AIPED \cdot (E - EFLOW)\)</span>,
where E is the in-place enthalpy and EFLOW is a specified enthalpy.</p>
</div>
<table border="1" class="docutils">
<colgroup>
<col width="42%" />
<col width="5%" />
<col width="5%" />
<col width="49%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td colspan="4">Isothermal model</td>
</tr>
<tr class="row-even"><td colspan="4">Case 1: AIPED = 0 (Constant Mass Rate, 1- or 2-Phase Source or Sink)</td>
</tr>
<tr class="row-odd"><td>SKD</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Mass source strength (kg/s). Negative value indicates injection into the rock mass.</td>
</tr>
<tr class="row-even"><td>EFLOW</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">a) EFLOW ≥ 0, EFLOW is the source liquid saturation,</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_w = SKD \cdot EFLOW\)</span> (kg/s)</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_a = SKD \cdot (1 - EFLOW)\)</span> (kg/s)</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">b) EFLOW &lt; 0, ABS(EFLOW) is the source air pressure (MPa)</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_w = SKD\)</span> (kg/s)</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_a = 1.0 \cdot (P_a - ABS(EEFLOW))\)</span> (kg/s)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>AIPED</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Used only as flag to indicatee Constant Mass Rate</td>
</tr>
</tbody>
</table>
<p>In the above and following relations, <span class="math notranslate nohighlight">\(Q_w\)</span> is the source term for water, <span class="math notranslate nohighlight">\(Q_a\)</span> is the source term for air,
and <span class="math notranslate nohighlight">\(P_a\)</span> is the in-place air pressure. The second case works well in situations where
inflow is specified and it is desired to hold the air pressure at a constant value.</p>
<table border="1" class="docutils">
<colgroup>
<col width="44%" />
<col width="3%" />
<col width="2%" />
<col width="50%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td colspan="4">Case 2: AIPED &gt; 0 (Constant Pressure, Constant Liquid Saturation Source or Sink)</td>
</tr>
<tr class="row-even"><td>SKD</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Specified source pressure (MPa).</td>
</tr>
<tr class="row-odd"><td>EFLOW</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><div class="first last line-block">
<div class="line">a) EFLOW &lt; 0, air only source.</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_a = AIPED \cdot (P_a - SKD)\)</span> (kg/s)</div>
<div class="line"><br /></div>
<div class="line">b) 0 &lt; EFLOW &lt; 1, EFLOW is specified source liquid saturation,</div>
<div class="line">for SKD ≥ 0, 2-phase source,</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_a = AIPED \cdot (P_a - SKD)\)</span> (kg/s)</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_w = AIPED \cdot (S_l - EEFLOW) \cdot P_0\)</span> (kg/s)</div>
<div class="line">when SKD &lt; 0, water only source <span class="math notranslate nohighlight">\(Q_a = 0\)</span>.</div>
<div class="line"><br /></div>
<div class="line">c) EFLOW = 1, water only source.</div>
<div class="line">* <span class="math notranslate nohighlight">\(Q_w = AIPED \cdot (P_l - SKD)\)</span> kg/s</div>
</div>
</td>
</tr>
<tr class="row-even"><td>AIPED</td>
<td>real</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">Impedance parameter. A large value is recommended (102 - 106) in order to create a flow</div>
<div class="line">term large enough to maintain constant pressure.</div>
</div>
</td>
</tr>
</tbody>
</table>
<table border="1" class="docutils">
<colgroup>
<col width="43%" />
<col width="5%" />
<col width="3%" />
<col width="48%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td colspan="4">Case 3: AIPED &lt; 0 (Outflow only, if <span class="math notranslate nohighlight">\(P_l &gt; SKD\)</span>)</td>
</tr>
<tr class="row-even"><td>SKD</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Pressure above which outflow occurs (MPa)</td>
</tr>
<tr class="row-odd"><td>EFLOW</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Not used.</td>
</tr>
<tr class="row-even"><td rowspan="2">AIPED</td>
<td rowspan="2">real</td>
<td rowspan="2"><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>Impedance parameter.</td>
</tr>
<tr class="row-odd"><td><span class="math notranslate nohighlight">\(Q_w = ABS(AIPED) \cdot R_l / \mu_l(P_l - SKD)\)</span> (kg/s)</td>
</tr>
</tbody>
</table>
<p>where <span class="math notranslate nohighlight">\(R_l\)</span> is the water relative peremeability and <span class="math notranslate nohighlight">\(\mu_l\)</span> is the water viscosity.</p>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="8%" />
<col width="24%" />
<col width="30%" />
<col width="13%" />
<col width="17%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>&#160;</td>
<td>PARAMETER</td>
<td colspan="2">SKD</td>
<td colspan="2">EFLOW</td>
</tr>
<tr class="row-even"><td>Physics Model</td>
<td><strong>AIPED*</strong></td>
<td><span class="math notranslate nohighlight">\(\ge 0\)</span></td>
<td><span class="math notranslate nohighlight">\(&lt; 0\)</span></td>
<td><span class="math notranslate nohighlight">\(\ge 0\)</span></td>
<td><span class="math notranslate nohighlight">\(&lt; 0\)</span></td>
</tr>
<tr class="row-odd"><td rowspan="3">Non-isothermal</td>
<td><span class="math notranslate nohighlight">\(&gt; 0\)</span></td>
<td>Flowing wellbore pressure (MPa)</td>
<td>Flowing wellbore pressure (MPa) - injection into rock mass</td>
<td rowspan="3">Enthalpy (MJ/kg)</td>
<td rowspan="3">Temperature (C)</td>
</tr>
<tr class="row-even"><td><span class="math notranslate nohighlight">\(= 0\)</span></td>
<td><div class="first last line-block">
<div class="line">Flow rate (kg/s)</div>
<div class="line">Heat only (MJ/s)</div>
</div>
</td>
<td><div class="first last line-block">
<div class="line">Flow rate (kg/s)</div>
<div class="line">Heat only (MJ/s) - injection into rock mass</div>
</div>
</td>
</tr>
<tr class="row-odd"><td><span class="math notranslate nohighlight">\(&lt; 0\)</span></td>
<td>Flowing wellbore pressure (MPa) - outflow only</td>
<td>N/A</td>
</tr>
<tr class="row-even"><td rowspan="3">Isothermal</td>
<td><span class="math notranslate nohighlight">\(&gt; 0\)</span></td>
<td>Specified source pressure (MPa)</td>
<td>Specified source pressure (MPa) - Water only source</td>
<td>Source liquid saturation</td>
<td>Air only source (used as flag)</td>
</tr>
<tr class="row-odd"><td><span class="math notranslate nohighlight">\(= 0\)</span></td>
<td>Mass source strength (kg/s)</td>
<td>Mass source strength (kg/s) - injection into rock mass</td>
<td>Source liquid saturation</td>
<td>Source air pressure (MPa)</td>
</tr>
<tr class="row-even"><td><span class="math notranslate nohighlight">\(&lt; 0\)</span></td>
<td>Pressure above which outflow occurs (MPa)</td>
<td>N/A</td>
<td>N/A</td>
<td>N/A</td>
</tr>
<tr class="row-odd"><td colspan="6">* Impedance parametere</td>
</tr>
</tbody>
</table>
<p>The following are examples of flow. In the first example, at node 88, a mass flow of 0.05 kg/s at 25 C
is being withdrawn from the model. Because fluid is being withdrawn the in-place temperature will actually be used.
For every 14th node from node 14 to 140, the pressure and temperature are being held constant at 3.6 MPa and 160 C,
respectively. This represents a constant temperature recharge along one of the problem boundaries.</p>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
<col width="17%" />
<col width="13%" />
<col width="23%" />
<col width="27%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>flow</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>88</td>
<td>88</td>
<td>1</td>
<td>0.050</td>
<td>-25.0</td>
</tr>
<tr class="row-odd"><td>14</td>
<td>140</td>
<td>14</td>
<td>3.600</td>
<td>-160.0</td>
</tr>
</tbody>
</table>
<p>In the second example, the corresponding input for airwater, is included, indicating
an isothermal air-water two-phase simulation is being run with a reference temperature
for property evaluation of 20 C and a reference pressure of 0.1 MPa. At nodes 26 and 52,
water saturation is 100% and water is being injected at 2.e-3 kg/s. At nodes 1 and 27,
there is an air only source, with a specified pressure of 0.1 MPa, and the air is
being injected at the rate of 100*(Pa - 0.1) kg/s.</p>
<table border="1" class="docutils">
<colgroup>
<col width="30%" />
<col width="15%" />
<col width="12%" />
<col width="24%" />
<col width="18%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>airwater</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>20.0</td>
<td>0.1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>flow</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>26</td>
<td>52</td>
<td>26</td>
<td>-2.e-3</td>
<td>1.0</td>
</tr>
<tr class="row-even"><td>1</td>
<td>27</td>
<td>26</td>
<td>0.1</td>
<td>-0.2</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
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
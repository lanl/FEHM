---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">eos</span></code><a class="headerlink" href="#eos" title="Permalink to this headline">¶</a></h1>
<p>Equation of State. Provide the code with alternate thermodynamic properties for the liquid and/or vapor phases. (This is one way in which the code may be instructed to simulate nonisothermal, single phase air. It may also be used to make comparisons between the code and analytical solutions that use different equations of state.)</p>
<ul class="simple">
<li>Group 1 - IIEOSD, IPSAT, ITSAT</li>
<li>Group 2 - EWI, EW2, EW3, EW4, EW5, EW6, EW7, EW8, EW9, EW10, EW11</li>
<li>Group 3 - EVI, EV2, EV3, EV4, EV5, EV6, EV7, EV8, EV9, EV10, EV11</li>
</ul>
<p>For calculation of the simplified thermodynamic equations the above data is used to generate first order equations. The exception to this is the viscosity of the liquid and use of the ideal gas law. The viscosity of the liquid uses a 1/T term. For the calculation of vapor density and its derivatives, the ideal gas law is used instead of a linear relationship. Thus, EV4 and EV5 are not used, but are included so the format is the same as that for the liquid parameters in Group 2.</p>
<table border="1" class="docutils">
<colgroup>
<col width="4%" />
<col width="2%" />
<col width="93%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IIEOSD</td>
<td>integer</td>
<td>Equation of state reference number. When IIEOSD = 1 or 2 are used, they refer to the high and low pressure data sets, respectively, in FEHM. For these values the input in Group 2 and Group 3 will be ignored after it is entered. When any value other than 1 or 2 are used, the user-defined equation of state is used with Groups 2 and 3 for input.</td>
</tr>
<tr class="row-odd"><td>IPSAT</td>
<td>integer</td>
<td>Parameter to set vapor pressure to zero. If IPSAT ≠ 0 the vapor pressure is set to zero, otherwise the vapor pressure is calculated in the code.</td>
</tr>
<tr class="row-even"><td>ITSAT</td>
<td>integer</td>
<td>Parameter to adjust the saturation temperature. If ITSAT &lt; 0, the saturation temperature is set to -1000oC. If ITSAT &gt; 0, the saturation temperature is set to 1000oC. If ITSAT = 0, the calculated value is used.</td>
</tr>
<tr class="row-odd"><td>EW1</td>
<td>real</td>
<td>Liquid reference pressure (MPa).</td>
</tr>
<tr class="row-even"><td>EW2</td>
<td>real</td>
<td>Liquid reference temperature (oC).</td>
</tr>
<tr class="row-odd"><td>EW3</td>
<td>real</td>
<td>Liquid reference density (kg/m3).</td>
</tr>
<tr class="row-even"><td>EW4</td>
<td>real</td>
<td>Derivative of liquid density with respect to pressure at reference conditions.</td>
</tr>
<tr class="row-odd"><td>EW5</td>
<td>real</td>
<td>Derivative of liquid density with respect to temperature at reference conditions.</td>
</tr>
<tr class="row-even"><td>EW6</td>
<td>real</td>
<td>Liquid reference enthalpy (MJ/kg).</td>
</tr>
<tr class="row-odd"><td>EW7</td>
<td>real</td>
<td>Derivative of liquid enthalpy with respect to pressure at reference conditions.</td>
</tr>
<tr class="row-even"><td>EW8</td>
<td>real</td>
<td>Derivative of liquid enthalpy with respect to temperature at reference conditions.</td>
</tr>
<tr class="row-odd"><td>EW9</td>
<td>real</td>
<td>Liquid reference viscosity (Pa s).</td>
</tr>
<tr class="row-even"><td>EW10</td>
<td>real</td>
<td>Derivative of liquid viscosity with respect to pressure at reference conditions.</td>
</tr>
<tr class="row-odd"><td>EW11</td>
<td>real</td>
<td>Derivative of liquid viscosity with respect to temperature at reference conditions.</td>
</tr>
<tr class="row-even"><td>EV1</td>
<td>real</td>
<td>Vapor reference pressure (MPa).</td>
</tr>
<tr class="row-odd"><td>EV2</td>
<td>real</td>
<td>Vapor reference temperature (oC).</td>
</tr>
<tr class="row-even"><td>EV3</td>
<td>real</td>
<td>Vapor reference density (kg/m3).</td>
</tr>
<tr class="row-odd"><td>EV4</td>
<td>real</td>
<td>Not used, included only to maintain a similar format to Group 2. Density variation with pressure governed by ideal gas law.</td>
</tr>
<tr class="row-even"><td>EV5</td>
<td>real</td>
<td>Not used, included only to maintain a similar format to Group 2. Density variation with temperature governed by ideal gas law.</td>
</tr>
<tr class="row-odd"><td>EV6</td>
<td>real</td>
<td>Vapor reference enthalpy (MJ/kg).</td>
</tr>
<tr class="row-even"><td>EV7</td>
<td>real</td>
<td>Derivative of vapor enthalpy with respect to pressure at reference conditions.</td>
</tr>
<tr class="row-odd"><td>EV8</td>
<td>real</td>
<td>Derivative of vapor enthalpy with respect to temperature at reference conditions.</td>
</tr>
<tr class="row-even"><td>EV9</td>
<td>real</td>
<td>Vapor reference viscosity (Pa s).</td>
</tr>
<tr class="row-odd"><td>EV10</td>
<td>real</td>
<td>Derivative of vapor viscosity with respect to pressure at reference conditions.</td>
</tr>
<tr class="row-even"><td>EV11</td>
<td>real</td>
<td>Derivative of vapor viscosity with respect to temperature at reference conditions.</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">eos</span></code>.</p>
<p>In this example, a user-defined equation of state is specified and the vapor pressure and saturation temperature are calculated in the code.
For liquid properties, the reference pressure is <span class="math notranslate nohighlight">\(0.1\:MPa\)</span>, the reference temperature is 20 °C, and the reference density is <span class="math notranslate nohighlight">\(998.\:kg/m^3\)</span>,
the derivative of density with respect to pressure is zero and with respect to temperature is <span class="math notranslate nohighlight">\(-0.2\:kg/m^3/C\)</span>.</p>
<p>The reference enthalpy is <span class="math notranslate nohighlight">\(0.88\:MJ/kg\)</span>, the derivative of enthalpy with pressure is zero, and the derivative with temperature is <span class="math notranslate nohighlight">\(4.2 \cdot 10^{-03}\:MJ/kg/C\)</span>.</p>
<p>The reference viscosity is <span class="math notranslate nohighlight">\(9 \cdot 10^{-04}\:Pa \cdot s\)</span> and the derivatives of viscosity with pressure and temperature are zero.</p>
<p>For vapor properties, the reference pressure is 0.1 MPa, the reference temperature is 20 oC, and the reference density is <span class="math notranslate nohighlight">\(1.29\:kg/m^3\)</span>.</p>
<p>The reference enthalpy is <span class="math notranslate nohighlight">\(2.5\:MJ/kg\)</span>, the derivative of enthalpy with pressure is 0, and with temperature is <span class="math notranslate nohighlight">\(0.1\:MJ/kg/C\)</span>.</p>
<p>The reference viscosity is 2.e-4 Pa⋅s and its derivatives with pressure and temperature are zero.</p>
<table border="1" class="docutils">
<colgroup>
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="8%" />
<col width="8%" />
<col width="9%" />
<col width="8%" />
<col width="13%" />
<col width="11%" />
<col width="8%" />
<col width="8%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>eos</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>3</td>
<td>0</td>
<td>0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>0.1</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>998</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ul class="first last simple">
<li></li>
</ul>
</td>
<td>0.8</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>4.2e-</td>
<td>9.e-</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-even"><td>0.1</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>1.2</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>0</td>
<td>8</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>3</td>
<td>4</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>9</td>
<td>&#160;</td>
<td>.</td>
<td>2.5</td>
<td>&#160;</td>
<td>0.1</td>
<td>2.e-</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>4</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>&#160;</td>
<td>&#160;</td>
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
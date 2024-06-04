---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">ppor</span></code><a class="headerlink" href="#ppor" title="Permalink to this headline">¶</a></h1>
<ul class="simple">
<li>Group 1 -     IPOROS</li>
<li>Group 2 - JA, JB, JC, POR1, POR2, POR3, POR4</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="26%" />
<col width="8%" />
<col width="66%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IPOROS</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Model type:</div>
<div class="line">IPOROS = 1, aquifer compressibility model</div>
<div class="line">IPOROS = -1, specific storage model (use</div>
<div class="line">only for isothermal conditions)</div>
<div class="line">IPOROS = -2, Gangi model (not available for</div>
<div class="line">air-water-heat conditions)</div>
</div>
</td>
</tr>
<tr class="row-odd"><td><strong>Model (1):</strong> IPOROS = 1,</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">Aquifer compressibility
<span class="math notranslate nohighlight">\(\phi = \phi_0 + \alpha_a(P-P_0)\)</span> where</div>
<div class="line"><span class="math notranslate nohighlight">\(\alpha_a\)</span> = aquifer compressibility (MPa-1),</div>
<div class="line"><span class="math notranslate nohighlight">\(\phi_0\)</span> = initial porosity,</div>
<div class="line"><span class="math notranslate nohighlight">\(P_0\)</span> = initial pressure (MPa)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>POR1</td>
<td>real</td>
<td>Aquifer compressibility <span class="math notranslate nohighlight">\(\alpha (MPa^{-1})\)</span></td>
</tr>
<tr class="row-odd"><td><strong>Model (-1):</strong> IPOROS = -1,</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">Specific storage
<span class="math notranslate nohighlight">\(S_s = \rho(\alpha_a + \phi \beta)\)</span> where</div>
<div class="line"><span class="math notranslate nohighlight">\(\rho\)</span> = liquid density (kg/m3),</div>
<div class="line"><span class="math notranslate nohighlight">\(g\)</span> = gravity,</div>
<div class="line"><span class="math notranslate nohighlight">\(\alpha_a\)</span> = aquifer compressibility (MPa-1),</div>
<div class="line"><span class="math notranslate nohighlight">\(\phi\)</span> = porosity,</div>
<div class="line"><span class="math notranslate nohighlight">\(\beta\)</span> = liquid compressibility (MPa-1)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>POR1</td>
<td>real</td>
<td>Specific storage <span class="math notranslate nohighlight">\(S_S (m^{-1})\)</span></td>
</tr>
<tr class="row-odd"><td><strong>Model (-2):</strong> IPOROS = -2,</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">Gangi model with calculation of initial permeability
and porosity.</div>
<div class="line"><span class="math notranslate nohighlight">\(\phi = \phi_0 \left[ 1 - \left(\frac{P_c}{P_x}\right)^m \right]\)</span>
and <span class="math notranslate nohighlight">\(P_c = \sigma - P - \alpha E(T-T_0)\)</span></div>
<div class="line">where</div>
<div class="line"><span class="math notranslate nohighlight">\(\phi_0\)</span> = initial porosity,</div>
<div class="line"><span class="math notranslate nohighlight">\(m\)</span> = Gangi exponent,</div>
<div class="line"><span class="math notranslate nohighlight">\(P_x\)</span> = fitted parameter (MPa)</div>
<div class="line"><br /></div>
<div class="line">Note: for the Gangi model the permeability is varied by
<span class="math notranslate nohighlight">\(k = k_0 \left(\frac{\phi}{\phi_0}\right)^3\)</span></div>
</div>
</td>
</tr>
<tr class="row-even"><td>POR1</td>
<td>real</td>
<td>Exponent <span class="math notranslate nohighlight">\(m\)</span> in Gangi bed of nails model.</td>
</tr>
<tr class="row-odd"><td>POR2</td>
<td>real</td>
<td><span class="math notranslate nohighlight">\(P_x\)</span> parameter (MPa) in Gangi equation.</td>
</tr>
<tr class="row-even"><td>POR3</td>
<td>real</td>
<td><span class="math notranslate nohighlight">\(\sigma\)</span> in-situ stress (MPa).</td>
</tr>
<tr class="row-odd"><td>POR4</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line"><span class="math notranslate nohighlight">\((\sigma E)\)</span> The product of the coefficient
of thermal expansion for the rock and the Young’s
modulus (MPa/C).</div>
<div class="line"><br /></div>
<div class="line">Note: For isothermal simulations the thermal term does not apply.</div>
</div>
</td>
</tr>
</tbody>
</table>
<p>In the following example of ppor, aquifer compressibility is modeled. All nodes
in the model are assigned a compressibility of 1.e-2 MPa-1.</p>
<table border="1" class="docutils">
<colgroup>
<col width="32%" />
<col width="16%" />
<col width="16%" />
<col width="37%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>ppor</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>0</td>
<td>0</td>
<td>1.e-2</td>
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
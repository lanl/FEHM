---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">rlp</span></code><a class="headerlink" href="#rlp" title="Permalink to this headline">¶</a></h1>
<p>Relative permeability and capillary pressure model. Several models are available.</p>
<ul class="simple">
<li>Group 1 - IRLP(i), RP1, RP2, RP3, RP4, RP5, RP6, RP7, RP8, RP9, RP10, RP11, RP12, RP13, RP14, RP15, RP16 (number of parameters entered depends on model selected)</li>
<li>Group 2 - JA, JB, JC, I</li>
</ul>
<p>Only those parameters defined for a given model need to be input. Group 1 is
ended when a blank line is encountered. The parameter i is incremented each
time a Group 1 line is read. Group 2 lines will refer to this parameter.
For model numbers 4, 6, and 7 (the combined van Genuchten model), the
permeability is isotropic and overwrites the input from macro <code class="docutils literal notranslate"><span class="pre">perm</span></code>.
Macro fper can be used with models 4, 6, and 7 to introduce anisotropy.</p>
<table border="1" class="docutils">
<colgroup>
<col width="44%" />
<col width="6%" />
<col width="51%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>IRLP(i)</td>
<td>integer</td>
<td>Relative permeability model type.</td>
</tr>
<tr class="row-odd"><td colspan="3"><strong>Model -1:</strong> IRLP(i) = -1, constant relative permeability, linear
capillary pressure (4 parameters required).</td>
</tr>
<tr class="row-even"><td>RP1</td>
<td>real</td>
<td>Liquid relative permeability (<span class="math notranslate nohighlight">\(m^2\)</span>).</td>
</tr>
<tr class="row-odd"><td>RP2</td>
<td>real</td>
<td>Vapor relative permeability (<span class="math notranslate nohighlight">\(m^2\)</span>).</td>
</tr>
<tr class="row-even"><td>RP3</td>
<td>real</td>
<td>Capillary pressure at zero saturation (<span class="math notranslate nohighlight">\(MPa\)</span>).</td>
</tr>
<tr class="row-odd"><td>RP4</td>
<td>real</td>
<td>Saturation at which capillary pressure goes to zero.</td>
</tr>
<tr class="row-even"><td colspan="3"><strong>Model 1:</strong>  IRLP(i) = 1, linear relative permeability, linear
capillary pressure (6 parameters required).</td>
</tr>
<tr class="row-odd"><td>RP1</td>
<td>real</td>
<td>Residual liquid saturation.</td>
</tr>
<tr class="row-even"><td>RP2</td>
<td>real</td>
<td>Residual vapor saturation.</td>
</tr>
<tr class="row-odd"><td>RP3</td>
<td>real</td>
<td>Maximum liquid saturation.</td>
</tr>
<tr class="row-even"><td>RP4</td>
<td>real</td>
<td>Maximum vapor saturation.</td>
</tr>
<tr class="row-odd"><td>RP5</td>
<td>real</td>
<td>Capillary pressure at zero saturation (MPa).</td>
</tr>
<tr class="row-even"><td>RP6</td>
<td>real</td>
<td>Saturation at which capillary pressure goes to zero.</td>
</tr>
<tr class="row-odd"><td colspan="3"><strong>Model 2</strong>:  IRLP(i) = 2, Corey relative permeability, linear
capillary pressure (4 parameters required).</td>
</tr>
<tr class="row-even"><td>RP1</td>
<td>real</td>
<td>Residual liquid saturation.</td>
</tr>
<tr class="row-odd"><td>RP2</td>
<td>real</td>
<td>Residual vapor saturation.</td>
</tr>
<tr class="row-even"><td>RP3</td>
<td>real</td>
<td>Capillary pressure at zero saturation (MPa).</td>
</tr>
<tr class="row-odd"><td>RP4</td>
<td>real</td>
<td>Saturation at which capillary pressure goes to zero.</td>
</tr>
<tr class="row-even"><td colspan="3"><strong>Model 3</strong>:  IRLP(i) = 3, van Genuchten relative permeability, van
Genuchten capillary pressure (6 parameters required). In this
model permeabilities are represented as a function of capillary
pressure [rlp(h)].</td>
</tr>
<tr class="row-odd"><td>RP1</td>
<td>real</td>
<td>Residual liquid saturation.</td>
</tr>
<tr class="row-even"><td>RP2</td>
<td>real</td>
<td>Maximum liquid saturation.</td>
</tr>
<tr class="row-odd"><td>RP3</td>
<td>real</td>
<td>Inverse of air entry head, <span class="math notranslate nohighlight">\(\alpha_G\)</span> (1/m) [note
some data is given in (1/Pa) convert using pressure = ρgΔh].</td>
</tr>
<tr class="row-even"><td>RP4</td>
<td>real</td>
<td>Power n in van Genuchten formula.</td>
</tr>
<tr class="row-odd"><td>RP5</td>
<td>real</td>
<td>Low saturation fitting parameter, multiple of
cutoff capillary pressure assigned as maximum
capillary pressure. If RP5 &lt; 0 then a linear
fit from the cutoff saturation (RP6) is used.
The slope of the cutoff saturation is used
to extend the function to saturation = 0.
If RP5 = 0, a cubic fit is used. The slope
at the cutoff saturation is matched and the
conditions <span class="math notranslate nohighlight">\(\frac{\partial}{\partial S}Pcap = 0\)</span>
and <span class="math notranslate nohighlight">\(\frac{\partial^2}{\partial S}Pcap = 0\)</span>
are forced at <span class="math notranslate nohighlight">\(S = 0\)</span>.
If RP5 &gt; 0, a multiple of the
value of the capillary pressure at
the cutoff saturation,  <span class="math notranslate nohighlight">\(RP5\cdot Pcap(S_{cutoff}\)</span>
is forced at <span class="math notranslate nohighlight">\(S = 0\)</span>.</td>
</tr>
<tr class="row-even"><td>RP6</td>
<td>real</td>
<td>Cutoff saturation used in fits described for RP5,
must be greater than RP1.</td>
</tr>
<tr class="row-odd"><td colspan="3"><strong>Model 4</strong>:  IRLP(i) = 4, van Genuchten relative permeability, van
Genuchten capillary pressure, effective continuum (15 parameters
required). In this model permeabilities are represented as a
function of capillary pressure [rlp(h)].</td>
</tr>
<tr class="row-even"><td>RP1</td>
<td>real</td>
<td>Residual liquid saturation, matrix rock material.</td>
</tr>
<tr class="row-odd"><td>RP2</td>
<td>real</td>
<td>Maximum liquid saturation, matrix rock material.</td>
</tr>
<tr class="row-even"><td>RP3</td>
<td>real</td>
<td>Inverse of air entry head, <span class="math notranslate nohighlight">\(\alpha_G\)</span> (1/m) [note
some data is given in (1/Pa) convert using
pressure = <span class="math notranslate nohighlight">\(\rho g \Delta h\)</span>], matrix rock material.</td>
</tr>
<tr class="row-odd"><td>RP4</td>
<td>real</td>
<td>Power n in van Genuchten formula, matrix rock material.</td>
</tr>
<tr class="row-even"><td>RP5</td>
<td>real</td>
<td>Low saturation fitting parameter, matrix rock
material, multiple of cutoff capillary pressure
assigned as maximum capillary pressure. If
RP5 &lt; 0 then a linear fit from the cutoff
saturation (RP6) is used. The slope of the
cutoff saturation is used to extend the
function to saturation = 0.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>If RP5 = 0, a cubic fit is used. The slope
at the cutoff saturation is matched and the
conditions  <span class="math notranslate nohighlight">\(\frac{\partial}{\partial S}Pcap = 0\)</span>
and  <span class="math notranslate nohighlight">\(\frac{\partial^2}{\partial S}Pcap = 0\)</span> are
forced at  <span class="math notranslate nohighlight">\(S = 0\)</span>. If RP5 &gt; 0, a multiple
of the value of the capillary pressure at the
cutoff saturation, <span class="math notranslate nohighlight">\(RP5\cdot Pcap(S_{cutoff})\)</span>
is forced at <span class="math notranslate nohighlight">\(S = 0\)</span>.</td>
</tr>
<tr class="row-even"><td>RP6</td>
<td>real</td>
<td>Cutoff saturation used in fits described for RP5,
must be greater than RP1.</td>
</tr>
<tr class="row-odd"><td>RP7</td>
<td>real</td>
<td>Residual liquid saturation, fracture material.</td>
</tr>
<tr class="row-even"><td>RP8</td>
<td>real</td>
<td>Maximum liquid saturation, fracture material.</td>
</tr>
<tr class="row-odd"><td>RP9</td>
<td>real</td>
<td>Inverse of air entry pressure, <span class="math notranslate nohighlight">\(\alpha_G\)</span> (1/m) [note
some data is given in(1/Pa) convert using
pressure = <span class="math notranslate nohighlight">\(\rho g \Delta h\)</span>], fracture material.</td>
</tr>
<tr class="row-even"><td>RP10</td>
<td>real</td>
<td>Power n in van Genuchten formula, fracture material.</td>
</tr>
<tr class="row-odd"><td>RP11</td>
<td>real</td>
<td>Low saturation fitting parameter, fracture material,
multiple of cutoff capillary pressure assigned as
maximum capillary pressure. If RP11 &lt; 0 then a linear
fit from the cutoff saturation (RP12) is used. The
slope of the cutoff saturation is used to extend
the function to saturation = 0.If RP11 = 0, a cubic
fit is used. The slope at the cutoff saturation is
matched and the conditions
<span class="math notranslate nohighlight">\(\frac{\partial}{\partial S}Pcap = 0\)</span>
and <span class="math notranslate nohighlight">\(\frac{\partial^2}{\partial S}Pcap = 0\)</span>
are forced at <span class="math notranslate nohighlight">\(S = 0\)</span>. If RP11 &gt; 0,
a multiple of the value of the capillary
pressure at the cutoff saturation,
<span class="math notranslate nohighlight">\(RP11\cdot Pcap(S_{cutoff})\)</span> is forced at
<span class="math notranslate nohighlight">\(S = 0\)</span>.</td>
</tr>
<tr class="row-even"><td>RP12</td>
<td>real</td>
<td>Cutoff saturation used in fits described for RP11,
must be greater than RP7.</td>
</tr>
<tr class="row-odd"><td>RP13</td>
<td>real</td>
<td>Fracture permeability (<span class="math notranslate nohighlight">\(m^2\)</span>). This is the permeability
of the individual fractures. The bulk permeability
of the fracture continuum is <span class="math notranslate nohighlight">\(RP13\times RP15\)</span>.
Can be made anisotropic with macro FPER.</td>
</tr>
<tr class="row-even"><td>RP14</td>
<td>real</td>
<td>Matrix rock saturated permeability (<span class="math notranslate nohighlight">\(m^2\)</span>). Can be made
anisotropic with macro FPER.</td>
</tr>
<tr class="row-odd"><td>RP15</td>
<td>real</td>
<td>Fracture volume fraction. Is equal to the
fracture aperture divided by the fracture
spacing (with same units). Sometimes called fracture porosity.</td>
</tr>
<tr class="row-even"><td colspan="3"><strong>Model 5</strong>:  IRLP(i) = 5, van Genuchten relative permeability, van
Genuchten capillary pressure (6 parameters required). This model
and its input are the same as for Model 3 except that
permeabilities are represented as a function of saturation
[rlp(S)] rather than capillary pressure.</td>
</tr>
<tr class="row-odd"><td colspan="3"><strong>Model 6</strong>:  IRLP(i) = 6, van Genuchten relative permeability, van
Genuchten capillary pressure, effective continuum (15 parameters
required). This model and its input are the same as for Model 4
except that permeabilities are represented as a function of
saturation [rlp(S)] rather than capillary pressure.</td>
</tr>
<tr class="row-even"><td colspan="3"><strong>Model 7</strong>:  IRLP(i) = 7, van Genuchten relative permeability, van
Genuchten capillary pressure, effective continuum with special
fracture interaction term (16 parameters required). This model
and its input are the same as for Model 6 except that the an
additional term is included which represents the fracture-matrix
interaction.</td>
</tr>
<tr class="row-odd"><td>RP16</td>
<td>real</td>
<td>Fracture-matrix interaction term. If RP16 ≤ 0,
then an additional multiplying term equal to the
relative permeability is applied to the
fracture-matrix interaction term for dual
permeability problems. If RP16 &gt; 0, then an
additional multiplying term equal to <span class="math notranslate nohighlight">\(sl**RP16\)</span>
and <span class="math notranslate nohighlight">\((1.-sl)**RP16\)</span> is applied to the
fracture-matrix interaction terms for the
liquid and vapor phases, respectively, for dual
permeability problems. Here, <span class="math notranslate nohighlight">\(sl\)</span>
is the value of saturation at the given node.</td>
</tr>
<tr class="row-even"><td colspan="3"><strong>Model 10</strong>: IRLP(i) = 10, linear relative permeability with minimum
relative permeability values, linear capillary pressure (8
parameters required).</td>
</tr>
<tr class="row-odd"><td>RP1</td>
<td>real</td>
<td>Residual liquid saturation.</td>
</tr>
<tr class="row-even"><td>RP2</td>
<td>real</td>
<td>Residual vapor saturation.</td>
</tr>
<tr class="row-odd"><td>RP3</td>
<td>real</td>
<td>Maximum liquid saturation.</td>
</tr>
<tr class="row-even"><td>RP4</td>
<td>real</td>
<td>Maximum vapor saturation.</td>
</tr>
<tr class="row-odd"><td>RP5</td>
<td>real</td>
<td>Minimum liquid permeability (<span class="math notranslate nohighlight">\(m^2\)</span>).</td>
</tr>
<tr class="row-even"><td>RP6</td>
<td>real</td>
<td>Minimum vapor permeability (<span class="math notranslate nohighlight">\(m^2\)</span>).</td>
</tr>
<tr class="row-odd"><td>RP7</td>
<td>real</td>
<td>Capillary pressure at zero saturation (<span class="math notranslate nohighlight">\(MPa\)</span>).</td>
</tr>
<tr class="row-even"><td>RP8</td>
<td>real</td>
<td>Saturation at which capillary pressure goes to zero.</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">rlp</span></code>.
In this example, Corey type relative permeability is specified,
with residual liquid saturation of 0.3, residual vapor saturation of 0.1,
a base capillary pressure of 2 MPa, and capillary pressure goes to zero at a
saturation of 1. This model is assigned to nodes numbered 1 through 140.</p>
<table border="1" class="docutils">
<colgroup>
<col width="21%" />
<col width="21%" />
<col width="21%" />
<col width="21%" />
<col width="17%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>rlp</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>2</td>
<td>0.3</td>
<td>0.1</td>
<td>2.0</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>140</td>
<td>1</td>
<td>1</td>
<td>&#160;</td>
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
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">grad</span></code><a class="headerlink" href="#grad" title="Permalink to this headline">¶</a></h1>
<p>Gradient model input.</p>
<ul class="simple">
<li>Group 1 -     NGRAD</li>
<li>Group 2 -     IZONE_GRAD, CORDG, IDIRG, IGRADF, VAR0, GRAD1</li>
</ul>
<p>Group 2 is repeated (NGRAD times) for each gradient model being defined.</p>
<table border="1" class="docutils">
<colgroup>
<col width="19%" />
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
<tr class="row-even"><td>NGRAD</td>
<td>integer</td>
<td>Number of gradient models</td>
</tr>
<tr class="row-odd"><td>IZONE_GRAD</td>
<td>integer</td>
<td>Zone associated with model</td>
</tr>
<tr class="row-even"><td>CORDG</td>
<td>real</td>
<td>Reference coordinate (m)</td>
</tr>
<tr class="row-odd"><td>IDIRG</td>
<td>integer</td>
<td>Direction of gradient (1, 2, or 3).</td>
</tr>
<tr class="row-even"><td>IGRADF</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Variable to which gradient is applied.</div>
<div class="line">1 = initial pressure</div>
<div class="line">2 = initial temperature</div>
<div class="line">3 = initial saturation</div>
<div class="line">4 = Fixed pressure</div>
<div class="line">5 = Fixed enthalpy</div>
<div class="line">-5 = Fixed inflowing temperature</div>
<div class="line">6 = Initial methane pressure</div>
<div class="line">7 = Fixed Methane pressure</div>
<div class="line">8 = depending on how <code class="docutils literal notranslate"><span class="pre">hflx</span></code> macro is configured
this is either a temperature or a heat flux</div>
<div class="line">9 = initial <span class="math notranslate nohighlight">\(CO_2\)</span> pressure</div>
<div class="line">10 = Fixed <span class="math notranslate nohighlight">\(CO_2\)</span> cell pressure</div>
<div class="line">11 = Pressure for secondary material in</div>
<div class="line-block">
<div class="line"><code class="docutils literal notranslate"><span class="pre">gdkm</span></code> or <code class="docutils literal notranslate"><span class="pre">gdpm</span></code> model</div>
</div>
<div class="line">12 = initial Temperature for matrix in <code class="docutils literal notranslate"><span class="pre">gdkm</span></code>
or <code class="docutils literal notranslate"><span class="pre">gdpm</span></code> model</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>VAR0</td>
<td>real</td>
<td>Value of variable at reference point (m).</td>
</tr>
<tr class="row-even"><td>GRAD1</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Gradient.</div>
<div class="line">Units: Pressure MPa/m, T degrees C/m,
enthalpy KJ/kg/m, heat flux MW/m</div>
</div>
</td>
</tr>
</tbody>
</table>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p>IIGRADF = 4,5,-5,7 requires that the node is previously defined as a boundary
node in a <code class="docutils literal notranslate"><span class="pre">flow</span></code> macro (or equivalent)</p>
<p>IGRAFD = 11 or 12 requires <code class="docutils literal notranslate"><span class="pre">gdkm</span></code> or <code class="docutils literal notranslate"><span class="pre">gdpm</span></code> macro</p>
<p class="last">IGRADF = 8 requires a <code class="docutils literal notranslate"><span class="pre">hflx</span></code> macro</p>
</div>
<p>The following is an example of the grad macro. A temperature gradient in the Y direction from the reference point of 0 will be applied to zone 1.</p>
<table border="1" class="docutils">
<colgroup>
<col width="21%" />
<col width="14%" />
<col width="11%" />
<col width="11%" />
<col width="18%" />
<col width="25%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>grad</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>1</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>2</td>
<td>2</td>
<td><ol class="first last arabic simple" start="10">
<li></li>
</ol>
</td>
<td>-150.</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
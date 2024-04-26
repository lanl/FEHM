---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">airwater</span> <span class="pre">or</span> <span class="pre">air</span></code><a class="headerlink" href="#airwater-or-air" title="Permalink to this headline">¶</a></h1>
<p>Isothermal air-water two-phase simulation.</p>
<ul class="simple">
<li>Group 1 -     ICO2D</li>
<li>Group 2 -     TREF, PREF</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="9%" />
<col width="76%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ICO2D</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Determines the type of air module used.</div>
<div class="line">ICO2D = 1, 1 degree of freedom solution to the saturated-unsaturated</div>
<div class="line-block">
<div class="line">problem is produced. This formulation is similar to the Richard’s</div>
<div class="line">Equation.</div>
</div>
<div class="line">ICO2D = 2, 1 degree of freedom solution is obtained assuming only gas flow</div>
<div class="line-block">
<div class="line">with no liquid present.</div>
</div>
<div class="line">ICO2D = 3, full 2 degree of freedom solution.</div>
<div class="line">All other values are ignored. The default is 3.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>TREF</td>
<td>real</td>
<td>Reference temperature for properties (oC).</td>
</tr>
<tr class="row-even"><td>PREF</td>
<td>real</td>
<td>Reference pressure for properties (MPa).</td>
</tr>
</tbody>
</table>
<p>Several macros are affected if the air module is enabled. These are:</p>
<ul class="simple">
<li><strong>pres</strong> - Because the air-water formulation is 2-phase at all times, care should be taken to insure that IEOSD is always specified to be 2. Likewise, saturations (not temperatures) are used.</li>
<li><strong>init</strong> - This macro should not be used because the saturation values cannot be specified.</li>
<li><strong>flow</strong> - A variety of different flow and boundary values are input with this macro when the macro <strong>airwater</strong> is also used. See description of control statement <strong>flow</strong>.</li>
</ul>
<p>The following is an example of <strong>airwater</strong>. In this example, a full 2-degrees-of-freedom solution is specified with a reference temperature for property evaluation of 20 oC and a reference pressure of 0.1 MPa.</p>
<table border="1" class="docutils">
<colgroup>
<col width="43%" />
<col width="57%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head" colspan="2">airwater</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>3</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>0.1</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
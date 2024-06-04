---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">cden</span></code><a class="headerlink" href="#cden" title="Permalink to this headline">¶</a></h1>
<p>Use concentration-dependent density for flow.</p>
<p>The following restrictions apply to the use of this macro:</p>
<ol class="arabic simple">
<li>It cannot be used with the macro “head”, which assumes constant fluid density</li>
<li>The updating of density is explicit, based on the concentration values at the previous time step. Therefore, accuracy of the solution must be tested by using a smaller time step and ensuring that the results have converged</li>
<li>The fluid flow time steps should be small enough that only one or two solute time steps are carried out before the next fluid time step, because relatively small changes in the concentration field are required for accuracy; and</li>
<li>The heat and mass transfer solution must be kept on during the entire simulation for the results to be meaningful (see macro trac).</li>
</ol>
<p>Group 1 - ISPCDEN</p>
<p>Group 2 - FACTCDEN</p>
<table border="1" class="docutils">
<colgroup>
<col width="6%" />
<col width="4%" />
<col width="90%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ISPCDEN</td>
<td>integer</td>
<td>The number of the chemical component in trac that is used for applying the concentration-dependent density.</td>
</tr>
<tr class="row-odd"><td>FACTCDEN</td>
<td>real</td>
<td>The factor used in the following relationship for fluid density (kg/m3):density = density_water + FACTCDEN*C where density_water = the density of pure water (kg/m3), and C is the concentration of chemical component ISPCDEN</td>
</tr>
</tbody>
</table>
<p>The following is an example of cden. In this example, component number 1 in trac is used. For concentrations of order 1, the density correction would be 100, of order 10% of the nominal value of water density of 1000 kg/m3.</p>
<table border="1" class="docutils">
<colgroup>
<col width="63%" />
<col width="38%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head" colspan="2">cden</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ISPCDEN</td>
<td>1</td>
</tr>
<tr class="row-odd"><td>FACTCDEN</td>
<td><ol class="first last arabic simple" start="100">
<li></li>
</ol>
</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
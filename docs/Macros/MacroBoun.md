---
layout : page_macros
hero_height: is-hidden
---


<h1><code class="docutils literal notranslate"><span class="pre">boun</span></code><a class="headerlink" href="#boun" title="Permalink to this headline">¶</a></h1>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Either boun or flow is required for a flow problem.</p>
</div>
<p>Implement boundary conditions and sources or sinks. Input may be time dependent and cyclic. Time step sizes may also be adjusted.</p>
<ul class="simple">
<li>Group 1 - <code class="docutils literal notranslate"><span class="pre">KEYWORD</span></code></li>
</ul>
<p>The Group 1 KEYWORD ‘model’, which starts each model sequence, is followed immediately by a Group 2 KEYWORD of ‘ti’, ‘ti_linear’, ‘cy’ or ‘cy_linear’.</p>
<ul class="simple">
<li>Group 2 - <code class="docutils literal notranslate"><span class="pre">KEYWORD</span></code></li>
<li>Group 3 - <code class="docutils literal notranslate"><span class="pre">NTIMES,</span> <span class="pre">TIME(I),</span> <span class="pre">I=1,NTIMES</span></code></li>
</ul>
<p>The Group 4 KEYWORDs define the various boundary condition parameters being entered. These KEYWORDs and associated data, Group 5, are repeated as needed for each model. Note that some keywords do not have associated variables.</p>
<ul class="simple">
<li>Group 4 - <code class="docutils literal notranslate"><span class="pre">KEYWORD</span></code></li>
<li>Group 5 - <code class="docutils literal notranslate"><span class="pre">VARIABLE(I),</span> <span class="pre">I=1,NTIMES</span></code></li>
</ul>
<p>Additional models are entered by beginning again with Group 1. The <code class="docutils literal notranslate"><span class="pre">MODEL_NUMBER</span></code> is incremented each time a new model is read, and is used to assign boundary conditions to specified nodes or zones in Group 6. After all models have been entered, the section is terminated with KEYWORD ‘end’ or a blank line.</p>
<table border="1" class="docutils">
<colgroup>
<col width="12%" />
<col width="9%" />
<col width="79%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td><div class="first last line-block">
<div class="line">Keyword specifying a model designation, time for boundary condition</div>
<div class="line">or source/sink changes, or actual variable or source change.</div>
<div class="line">Keywords, which must be entered starting in the first column, are:</div>
<div class="line-block">
<div class="line"><code class="docutils literal notranslate"><span class="pre">model</span></code> - new model definition to follow</div>
</div>
<div class="line"><em>Note</em>: Descriptive text, such as the model number, may be appended after</div>
<div class="line-block">
<div class="line">the ‘model’ keyword as long as it is contained on a single line, and</div>
<div class="line">begins after column four.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">ti</span></code> - time sequence for changes to follow (days). The ‘ti’ keyword results</div>
<div class="line-block">
<div class="line">in step function changes in boundary condition</div>
<div class="line">VARIABLE(i) at each TIME(i).</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">ti_linear</span></code> - time sequence for changes to follow (days). The ‘ti_linear’</div>
<div class="line-block">
<div class="line">keyword will apply a boundary condition that changes linearly</div>
<div class="line">with time. This option does not impose any control on time</div>
<div class="line">step size, so it is possible that a single time step can span</div>
<div class="line">an entire time interval and the linear change will not be seen.</div>
<div class="line">If time step size control is important it should be imposed</div>
<div class="line">in the time or ctrl macros.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">cy</span></code> - cyclic time sequence for changes to follow (days). As with the ‘ti’</div>
<div class="line-block">
<div class="line">keyword, boundary condition changes are step functions.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">cy_linear</span></code> - cyclic time sequence for changes to follow (days). As with the</div>
<div class="line-block">
<div class="line"><code class="docutils literal notranslate"><span class="pre">ti_linear</span></code> keyword, boundary condition changes linearly with time.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">sec</span></code> - Time sequence input is in seconds.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">min</span></code> - Time sequence input is in minutes.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">day</span></code> - Time sequence input is in days. (Default)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">year</span></code> - Time sequence input is in years.</div>
<div class="line"><em>Note</em>: The keywords, ‘ti’, ‘ti_linear’, ‘cy’ and ‘cy_linear’, require the</div>
<div class="line-block">
<div class="line">time to start at 0.0. This provides the initial boundary and</div>
<div class="line">source/sink information. If the input for ‘ti’ or ‘cy’ does not</div>
<div class="line">start at 0.0 time the code assumes boundary conditions and</div>
<div class="line">source/sinks are 0.0 at time 0.0. The ‘cy’ keyword involves a</div>
<div class="line">cyclic changing of conditions. In our procedure the cycle ends at</div>
<div class="line">the last specified time. Thus the code reverts to the first</div>
<div class="line">specified time values. Because of this, the boundary conditions</div>
<div class="line">and source/sinks for the last time change are always set to the first</div>
<div class="line">time values.The default units for time input is days. Input in seconds,</div>
<div class="line">minutes or years is converted to days. Time input units have no</div>
<div class="line">associated variable input.               |</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">tran</span></code> - Keyword to indicate boundary conditions will not be invoked until</div>
<div class="line-block">
<div class="line">the steady state portion of the simulation is completed and a</div>
<div class="line">transient run is initiated. See macro stea for more details.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">sa</span></code> - air source sequence for changes to follow (kg/s)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">sw</span></code> - water source sequence for changes to follow (kg/s)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">swf</span></code> - source water factor sequence for changes to follow (multiplier</div>
<div class="line-block">
<div class="line">for existing mass flow rates)</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">se</span></code> - enthalpy source sequence for changes to follow (MW)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">sf</span></code> - water seepage face sequence with pressures for changes to follow (MPa)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">sfh</span></code>  - water seepage face sequence with heads for changes to follow (m)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">fd</span></code>  - water drainage area sequence for changes to follow (m2)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">dsa</span></code> - distributed air source sequence for changes to follow (kg/s)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">dsw</span></code> - distributed water source sequence for changes to follow(kg/s)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">dse</span></code> - distributed enthalpy source sequence for changes to follow (MW)</div>
<div class="line"><em>Note</em>: A distributed source (keywords ‘dsa’, ‘dsw’, and ‘dse’) is a source</div>
<div class="line-block">
<div class="line">term divided over a group of nodes or a zone proportional to the nodal volume.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgt</span></code> - Distributed source is weighted using the nodal control volume.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtx</span></code> - Distributed source is weighted using nodal area = control volume / x length scale.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgty</span></code> - Distributed source is weighted using nodal area = control volume / y length scale.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtz</span></code> - Distributed source is weighted using nodal area = control volume / z length scale.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtp</span></code> - Distributed source is weighted using nodal control volume * permeability.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtpx</span></code> - Distributed source is weighted using nodal control volume * permeability / x length scale.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtpy</span></code> - Distributed source is weighted using nodal control volume * permeability / y length scale.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtpz</span></code> - Distributed source is weighted using nodal volume * permeability / z length scale.</div>
<div class="line"><em>Note</em>: The length scale term is the dimension of the control volume bounding box, xmax-xmin,</div>
<div class="line-block">
<div class="line">ymax-ymin, zmax-zmin, depending upon the suffix, x,y,z. This option is useful when one</div>
<div class="line">wants to apply a distributed source on a mesh with variable size mesh cells and would</div>
<div class="line">like the source percentage to be allocated based on surface area of each node.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtr</span></code> - Distributed source is weighted using nodal volume * permeability * relative permeability</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgtu</span></code> - Distributed source is weighted with user specified values (See macro wgtu).</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">wgww</span></code> - Distributed source is weighted using nodal volume * permeability * relative</div>
<div class="line-block">
<div class="line">permeability * exponentially weigthed distance from pump</div>
</div>
<div class="line"><em>Note</em>: The distributed source weighting options have no associated variable input.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">s</span></code> - fixed saturation sequence for changes to follow</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">hd</span></code> - fixed hydraulic head sequence for changes to follow (m)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">pw</span></code> - fixed water pressure sequence for changes to follow (MPa)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">pa</span></code> - fixed air pressure sequence for changes to follow (MPa)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">hdo</span></code> - fixed hydraulic head sequence for changes to follow (m) (constrained to outflow only)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">pwo</span></code> - fixed water pressure sequence for changes to follow (MPa) (constrained to outflow only)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">pao</span></code> - fixed air pressure sequence for changes to follow (MPa) (constrained to outflow only)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">en</span></code> - fixed enthalpy sequence for changes to follow (MW)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">t</span></code> - fixed temperature sequence for changes to follow (oC)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">h</span></code> - fixed humidity sequence for changes to follow (must be used with van Genuchten relative</div>
<div class="line-block">
<div class="line">permeability model)</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">ft</span></code> - fixed flowing temperature sequence for change to follow (oC). By flowing temperature</div>
<div class="line-block">
<div class="line">we mean the temperature of the inflow stream for a specified source. If no source</div>
<div class="line">inflow occurs where this condition is applied, it will be ignored.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">kx</span></code> - fixed X permeability sequence for changes to follow (m2)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">ky</span></code> - fixed Y permeability sequence for changes to follow (m2)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">kz</span></code> - fixed Z permeability sequence for changes to follow (m2)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">if</span></code> - impedance factor for use with fixed water pressure boundary condition.</div>
<div class="line-block">
<div class="line">If left out the impedance factor will be set to the volume of the grid cell.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">si</span></code>  - initial value saturation sequence for changes to follow</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">pai</span></code>  - initial value air pressure sequence for changes to follow (MPa)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">pwi</span></code>  - initial value water pressure sequence for changes to follow (MPa)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">tmi</span></code>  - initial value temperature sequence for changes to follow (oC)</div>
<div class="line"><em>Note</em>: The keywords ‘si’, ‘pai’, ‘pwi’, and ‘tmi’ refer to changes for a variable</div>
<div class="line-block">
<div class="line">that is NOT fixed. They are similar to specifying initial conditions in that</div>
<div class="line">regard but may be changed according to a time sequence. At present these 4</div>
<div class="line">keywords only work with isothermal air-water calculations.</div>
</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">chmo</span></code> - model number sequence for changes to follow</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">ts</span></code> - timestep sequence for changes to follow (days)</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">end</span></code> - signifies end of keyword input, a blank line will also work.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>NTIMES</td>
<td>integer</td>
<td>Number of time changes for boundary condition or source/sink specification.</td>
</tr>
<tr class="row-even"><td>TIME</td>
<td>real</td>
<td>NTIMES times for changes in boundary conditions or source/sinks.</td>
</tr>
<tr class="row-odd"><td>VARIABLE</td>
<td>real</td>
<td>NTIMES new values for boundary conditions or source/sinks.</td>
</tr>
<tr class="row-even"><td>MODEL_NUMBER</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Boundary condition model to be associated with designated nodes or zones</div>
<div class="line">(the number corresponds to the numerical order in which the models were</div>
<div class="line">input, i.e., beginning with KEYWORD <code class="docutils literal notranslate"><span class="pre">model</span></code>)</div>
</div>
</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">boun</span></code>. In this example two models are defined.
The first model cyclically changes the water source in a 1.e05 day cycle, i.e.,
the <code class="docutils literal notranslate"><span class="pre">cy</span></code> keyword entry shows that the time cycle ends at 1.e05 days and at this
time the cycle reverts to 0.0 days. Note that the water source at 1.e05 days
equals that at 0.0 days. Also in model 1 the flowing temperature was alternated
between 20oC and 50oC. The second model uses a time change that occurs at 1.e20 days.
This effectively removes any time variance from model 2. Model 2 has a fixed
water pressure and flowing temperature condition. The models are applied at
nodes 26 and 27 in the last two lines. It should be noted that the model numbers
included in the example (following KEYWORD <code class="docutils literal notranslate"><span class="pre">model</span></code>) are not part of the required
input but are descriptive text used to enhance readability of the macro.</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="16%" />
<col width="16%" />
<col width="16%" />
<col width="16%" />
<col width="16%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head" colspan="6">boun</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>model 1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>cy</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>4</td>
<td>0.0</td>
<td>1.e1</td>
<td>1.e2</td>
<td>1.e5</td>
</tr>
<tr class="row-odd"><td>sw</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>-1.e-4</td>
<td>-1.e-5</td>
<td>-1.e-3</td>
<td>-1.e-4</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>ft</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>20.0</td>
<td>50.0</td>
<td>50.0</td>
<td>20.0</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>model 2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>ti</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>2</td>
<td>0.0</td>
<td>1.e20</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>pw</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>0.1</td>
<td>0.1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>ft</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>20.0</td>
<td>20.0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>26</td>
<td>26</td>
<td>1</td>
<td>1</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>27</td>
<td>27</td>
<td>1</td>
<td>2</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>In the second example, a distributed water source is used to model a zone where
production is turned on and off. Keyword <code class="docutils literal notranslate"><span class="pre">kz</span></code> is used to specify a higher
permeability when production is occurring. The <code class="docutils literal notranslate"><span class="pre">ts</span></code> keyword is used to reset the
time step to 1.0 days, at the beginning of each time interval. The model is
applied to zone 100 in the last line.</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="16%" />
<col width="16%" />
<col width="16%" />
<col width="16%" />
<col width="16%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head" colspan="6">boun</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>model 1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>ti</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>4</td>
<td>0.0</td>
<td>91.325</td>
<td>182.62</td>
<td>273.93</td>
</tr>
<tr class="row-odd"><td>ts</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>1.0</td>
<td>1.0</td>
<td>1.0</td>
<td>1.0</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>dsw</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>29.248</td>
<td>0.0</td>
<td>29.248</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>kz</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>8e-12</td>
<td>2e-12</td>
<td>8e-12</td>
<td>2e-12</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>-100</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
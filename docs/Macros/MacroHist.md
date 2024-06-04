---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">hist</span></code><a class="headerlink" href="#hist" title="Permalink to this headline">¶</a></h1>
<p>History data output selection, output timestep intervals, and time intervals. Parameters will be output for the nodes specified in the node or nod2 macro in individual files for each parameter selected. If output zones are defined (node macro) the output will be a volume weighted average for the zone. Currently zone averaged values can be output for pressure, head, temperature, and enthalpy. History files will be named using the root portion of the history file name (<code class="docutils literal notranslate"><span class="pre">root_name.his</span></code>) provided as input, e.g., pressure output would be in a file named: <code class="docutils literal notranslate"><span class="pre">root_name_pres.his</span></code>. The named history output file will contain run information and the list of selected output nodes (with their coordinates) and zones (with node list).</p>
<ul class="simple">
<li>Group 1 -     CHDUM</li>
</ul>
<p>or using optional input or keywords</p>
<ul class="simple">
<li>Group 1 -     CHDUM, NHIST, HISTIME</li>
</ul>
<p>where CHDUM is ‘years’, ‘days’, ‘hrs’, or ‘seconds’, or</p>
<ul class="simple">
<li>Group 1 -     CHDUM, CHDUM1, … , CHDUMn</li>
</ul>
<p>where CHDUM is mpa, pressure, density, viscosity, or global</p>
<table border="1" class="docutils">
<colgroup>
<col width="5%" />
<col width="3%" />
<col width="92%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>CHDUM</td>
<td>character*80</td>
<td>Keyword specifying type of history plot data files to be created. Keywords are entered one per line and terminated with ‘end hist’ or a blank line. Keywords must be entered starting in the 1st column. Valid keywords (case insensitive) are:</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘tecplot’ - data will be output using tecplot style headers and format</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘csv’ or ‘surfer’ - data and parameter headers will be output as comma separated variables (‘csv’ format)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘years’ - output time in years</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘days’ - output time in days</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘hrs’ - output time in hours</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘seconds’ - output time in seconds</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘mpa’ or ‘pressure’ - output pressure in MPa</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘deg’ or ‘temperature’ - output temperature in oC</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘head’ or ‘meters’ - output head in meters</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘feet’ - output head in feet</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘saturation’ - output saturation</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘wco’ - water content</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘flow’ or ‘kgs’ - output flow in kg/s</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘enthalpy’ - output enthalpy in MJ/kg</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘efl’ or ‘mjs’ - output enthalpy flow (MJ/s)</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘density’ - output density (<span class="math notranslate nohighlight">\(kg/m^3\)</span>)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘humidity’ - output relative humidity</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘viscosity’ - output viscosity (Pa-s)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘zflux’ - output zone fluxes</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘concentration’ - output species concentration (concentrations for each specie will be output in a separate file)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘wt’ - output water table elevation</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2s’ - output <span class="math notranslate nohighlight">\(CO_2\)</span> saturations (volume fractions) An ‘l’ or ‘g’ may be appended to co2s to specify that only liquid or gas saturation should be output.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2m’ - output total <span class="math notranslate nohighlight">\(CO_2\)</span> mass (kg), free <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction, and dissolved <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction. Other wise use the form listed below to output specified quantity:</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2mt’ - output total <span class="math notranslate nohighlight">\(CO_2\)</span> mass (kg)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2mf’ - output free <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2md’ - output dissolved <span class="math notranslate nohighlight">\(CO_2\)</span> mass fraction</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘cfluxz’ - output <span class="math notranslate nohighlight">\(CO_2\)</span> zone fluxes</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘displacements’ - output displacements (m), ‘disx’, ‘disy’ or ‘disz’ may be used to select only the x, y, or z displacement.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘stress’ - output stresses, ‘strsx’, strsy’, ‘strsz’, ‘strsxy’, ‘strsxz’ or ‘strsyz’ may be used to select specific stress components.</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘strain’ - output strain</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘rel’ - output a table of relative permeability values for each input model.</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘global’ - output global parameters</td>
</tr>
<tr class="row-odd"><td>NHIST</td>
<td>integer</td>
<td>Optional: <em>Time step</em> interval for history plots (number of timesteps). Output history information each NHIST timesteps. If not entered NHIST = 1.</td>
</tr>
<tr class="row-even"><td>HISTIME</td>
<td>real</td>
<td>Optional: <em>Time</em> interval for history plots. In addition to output each NHIST timesteps, output history information each HISTIME. Units correspond to units specified by selected time output keyword (years, days, hours, or seconds). If not entered HISTIME = 1.e30.</td>
</tr>
<tr class="row-odd"><td>CHDUM1 … CHDUMn</td>
<td>character*80</td>
<td>Optional keywords specifying selections for history plot data files to be created. Optional keywords are entered on the same line as primary keywords. If no optional keywords are used the code will determine what will be output based on problem input. Up to 3 optional keywords may be entered.</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>Valid keywords (case insensitive) used with keyword ‘pressure’ or ‘mpa’ are:</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘total’ or ‘water’- output total or water pressure</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘air’ - output air pressure</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘capillary’ - output capillary pressure</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2’ - output <span class="math notranslate nohighlight">\(CO_2\)</span> pressure</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>Valid keywords (case insensitive) used with keyword ‘density’ or ‘viscosity’ are:</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘water’ - output water density or viscosity</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘air’ - output air/vapor density or viscosity</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘co2’ - output <span class="math notranslate nohighlight">\(CO_2\)</span> liquid and gas density or viscosity. An ‘l’ or ‘g’ may be appended to co2 to specify that only liquid or gas density should be output.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>Valid keywords (case insensitive) used with keyword ‘global’ are:</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘mass’ - output mass balances only for problem (excluding steam) (used with keyword ‘global’)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘water’ - output water balances only for problem (excluding steam) (used with keyword ‘global’)</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘steam’ - output mass/water balance only for problem including steam (used with keyword ‘global’)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>‘air’ = output air / vapor balances only for problem (used with keyword ‘global’)</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>‘energy’ - output energy balances only for problem (used with keyword ‘global’)</td>
</tr>
</tbody>
</table>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">If a time keyword (years, days, hrs, or seconds) is not entered, output time will be in days and data will be output for each timestep. The time output keywords may be used with optional input NHIST and HISTIME.</p>
</div>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">If a file format keyword is being used, it must immediately follow the macro name. Alternatively, it may be entered on the macro line. The default is for the headers and data to be output using plain text and spaces.</p>
</div>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">If no optional keywords are used with the ‘global’ keyword the code will determine which balances will be output based on problem input (mass/energy or water/air). Currently only 1 optional keyword may be used with global to specify a single balance type. Balance output includes: Total (mass, water, air in kg, energy in MJ) in system, total discharge, total input, current discharge, current input, and net discharge.</p>
</div>
<p>The following are examples of <code class="docutils literal notranslate"><span class="pre">hist</span></code>. For this first example, time will be output in years and temperatures in oC. Data will be output each 100000 timesteps or at time intervals of 50 years.</p>
<table border="1" class="docutils">
<colgroup>
<col width="35%" />
<col width="40%" />
<col width="25%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">hist</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>years</td>
<td>100000</td>
<td><ol class="first last arabic simple" start="50">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td>deg</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>In this second example, pressures in MPa (water and air) and temperatures in oC
will be written each timestep and time will be output in days. The global
mass balance for water will also be output at each time step.</p>
<table border="1" class="docutils">
<colgroup>
<col width="40%" />
<col width="35%" />
<col width="25%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">hist</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>mpa</td>
<td>total</td>
<td>air</td>
</tr>
<tr class="row-odd"><td>deg</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>global</td>
<td>mass</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>end</td>
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
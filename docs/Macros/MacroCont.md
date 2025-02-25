---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">cont</span></code><a class="headerlink" href="#cont" title="Permalink to this headline">¶</a></h1>
<p>Contour data output format, output timestep intervals, and time intervals.</p>
<ul class="simple">
<li>Group 1 - <code class="docutils literal notranslate"><span class="pre">NCNTR</span></code>, <code class="docutils literal notranslate"><span class="pre">CONTIM</span></code></li>
</ul>
<p>An alternative form of input for macro cont is possible. This is</p>
<ul class="simple">
<li>Group 1 - <code class="docutils literal notranslate"><span class="pre">ALTC</span></code>, <code class="docutils literal notranslate"><span class="pre">NCNTR</span></code>, <code class="docutils literal notranslate"><span class="pre">CONTIM</span></code>, <code class="docutils literal notranslate"><span class="pre">KEYWORD</span></code></li>
<li>Group 2 - CHDUM (only input if ALTC is ‘avs’, ‘avsx’, ‘surf’, or ‘tec’)</li>
</ul>
<p>If CHDUM = <code class="docutils literal notranslate"><span class="pre">'zone'</span></code> that line is followed by</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">NSURF</span>
<span class="n">IZONE_ISURF</span><span class="p">(</span><span class="n">I</span><span class="p">),</span> <span class="n">I</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">NSURF</span>
</pre></div>
</div>
<table border="1" class="docutils">
<colgroup>
<col width="1%" />
<col width="1%" />
<col width="98%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ALTC</td>
<td>character*4</td>
<td>Keyword specifying the type of contour output wanted (avs, avsx, fehm, free, ment, ptrn, surf, tec, vtk): 'avs' produces contour plot files compatible with the AVS postprocessor. ‘avsx’ produces contour plot files compatible with the AVS Express postprocessor.‘fehm’ produces a binary output file. The same contour plot file is produced using the first form of Group1 input. ‘free’ produces a free format contour plot file. ‘surf’ produces a contour plot file compatible with the SURFER postprocessor. ‘tec’ produces a contour plot file compatible with the TECPLOT postprocessor. ‘vtk’ produces a VTK contour plot file compatible with Paraview. </td>
</tr>
<tr class="row-odd"><td>NCNTR</td>
<td>integer</td>
<td>Time step interval for contour plots (number of timesteps). Output contour information each NCNTR timesteps.</td>
</tr>
<tr class="row-even"><td>CONTIM</td>
<td>real</td>
<td>Time interval for contour plots (days). In addition to output each NCNTR timesteps, output contour information each CONTIM days.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘time”, use time (days) in file name instead of number (used when altc is “avs” or “avsx” or “sur” or “tec”)</td>
</tr>
<tr class="row-even"><td>CHDUM</td>
<td>character*72</td>
<td>Keyword specifying type of contour plot data files to be created in AVS UCD, AVS Express, SURFER or TECPLOT format. Keywords are entered one per line and terminated with ‘endcont’ or ‘end cont’. Valid keywords (case insensitive) are:</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>(f)ormatted - output data in ASCII format.(m)aterial - output contour values for material properties.(l)iquid - output contour values for liquid phase.(va)por - output contour values for vapor phase.(dp)dp - output contour values for dual permeability nodes.(g)eo - output geometry values (coordinates and connectivity, avs style or old tecplot style).(gr)id - output grid geometry and connectivity in a tecplot grid file format. Parameter files will be output using tecplot variable format.(n)odit - do not output a contour file at each dit (see time macro).</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>(c)oncentration - output solute concentration values.(ca)pillary - output capillary pressure values.(co2) - output saturation values (liquid/supercritical liquid and gas).(de)nsity - output density values.(di)splacement - output x, y, and z displacements for stress problem. When ‘reldisp’ is specified in the strs macro relative displacements are output.(fh)ydrate - output htdrate fraction.(fl)ux - output node flux (additional keywords ‘net’ ‘volume’ ‘vwg’)(fw)water - output water fraction.(h)ead - output head values.(hy)drate - output hydrate values.(pe)rmeability - output permeability values.(po)rosity - output porosity values.(p)ressure - output pressure values.(s)aturation - output saturation values.(so)urce - output source values.(stra)in - output strain for a stress problem.(stre)ss - output defined stresses (x, y, z, xy, xz, yz).(t)emperature - output temperature values.(ve)locity - output velocity values.(wt) - output water table elevation.(x)yz - output node coordinates(zi)d - output number of zone containing this node (as defined at end of input file)(z)one - output values for specified zones (entered on following lines)(e)ndcont - last keyword entered.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>If a format keyword is not entered, the default is ‘formatted’. In the current version of the code this is the only format option supported. The default for data keywords is “off” (no output). The letters given in ( ) are sufficient to identify the keyword. The ‘zone’ and ‘geo’ keywords can not be used together. Geometry data will not be output if both keywords are specified.</td>
</tr>
<tr class="row-even"><td>NSURF</td>
<td>integer</td>
<td>Number of output zones (entered following ‘zone’ keyword).</td>
</tr>
<tr class="row-odd"><td>IZONE_SURF</td>
<td>integer</td>
<td>List of nsurf zone numbers (entered following ‘zone’ keyword).</td>
</tr>
</tbody>
</table>

<p>FEHM will automatically distinguish between the alternative input formats. When keywords are used they must be entered starting in the first column. The contour data will be output whenever either of the interval criteria are satisfied.</p>
<p>For keyword output, if the material keyword is selected, the following material property values (at the initial time) will be written for each node: permeability in the x, y, and z directions, thermal conductivity in the x, y, and z directions, porosity, rock specific heat, capillary pressure, relative permeability model being used, and capillary pressure model being used. If vapor and/or liquid are selected, pressure, velocity, or density must also be defined (otherwise, no data for these values will be written). velocity will result in vector values, other values will be scalar. If concentration is selected, values will be output only if nspeci is defined for tracer solutions. See the control statement trac for a description of nspeci for solutes.</p>
<p>The following are examples of cont. For the first example, FEHM binary format contour output files will be written every 100 timesteps and for each 1.e20 days. The second example invokes AVS contour output. AVS UCD formatted files will be written for every 100 time steps and 1.e20 days. The resulting files will include a log file, geometry file, plus header and data files for the following: material properties, solute concentrations, liquid velocities, pressures and temperatures.
The third example write a vtk contour file. The first number is for Time step interval, and the second is for Time interval.
</p>

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">cont</span>
<span class="mi">100</span>
<span class="mf">1.e20</span>
<span class="n">contavsmatconliquidvelocitypressuretempformattedendavs</span>
<span class="mi">100</span>
<span class="mf">1.e20</span>
</pre></div>
</div>

<div>
<pre>
cont
vtk   10000    10
</pre>
</div>

  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>

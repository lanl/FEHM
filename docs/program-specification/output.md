---
title : Output
layout : page_getting-started
permalink: /program-specification/output
hero_height: is-small
---

 Output is found in the code generated files (output file, write file, history plot file, solute plot file, contour plot file, contour plot file for dual or dpdp, stiffness matrix data file, input check file, submodel output file, PEST output files, streamline particle tracking output files, and AVS output files) described in <a href="/getting-started/data-files"> Data Files </a>. 

Macro commands (input options) dealing with output control are <a class="reference external" href="../Macros/MacroCont ">cont</a>, <a class="reference external"  href="../Macros/MacroCtrl ">ctrl</a>, <a class="reference external"  href="../Macros/MacroNod2 ">nod2</a>, <a class="reference external"  href="../Macros/MacroNode ">node</a>, <a class="reference external"  href="../Macros/MacroMptr ">mptr</a>, <a class="reference external"  href="../Macros/MacroPest ">pest</a>, <a class="reference external"  href="../Macros/MacroPtrk ">ptrk</a>, <a class="reference external"  href="../Macros/MacroSptr ">sptr</a>, <a class="reference external"  href="../Macros/MacroSubm ">subm</a>, <a class="reference external"  href="../Macros/MacroWflo ">wflo</a>, and <a class="reference external"  href="../Macros/MacroTime ">time</a>: cont is used to specify output format and time intervals for contour data output ( ```fehmn.con```,  ```fehmn.dp```); ctrl is used to specify if element coefficients calculated in the code should be saved ( ```fehmn.stor```); node and nod2 are used to provide nodal or coordinate positions for which general information and history data will be output ( ```fehmn.out```,  ```fehmn.his```,  ```fehmn.trc```, and terminal output); mptr has an option to specify whether or not particle tracking information is written to the restart file ( ```fehmn.fin```); pest is used to specify PEST parameter estimation routine output format ( ```fehmn.pest```,```fehmn.pest1```); ptrk has an option to specify whether or not particle tracking information is written to the restart file ( ```fehmn.fin```) and what information to output; sptr has options to specify what streamline particle tracking information will be output ( ```fehmn.sptr1```,```fehmn.sptr2```,```fehmn.sptr3```); submand wflo are used to specify nodes and boundary conditions should be output for an extracted submodel region; and time provides input on the time printout interval for nodal information ( ```fehmn.out``` and terminal output).

 The code itself provides no graphical capabilities. History plots of the energy source, source strength, temperature, pressure, capillary pressure, and saturation are made from the filen.his FEHM output files. Data from the filen.trc files is used to make history tracer plots of the 10 species concentrations. Contour plots can be made from the  ```filen.con```,  ```filen.dp```, and AVS FEHM output files. 

 AVS provides tools for visualization and analysis of volumetric scalar and vector data Contour plots using 2-d quad grids and 3-d hex grids for material properties, temperature, saturation, pressure, velocities, and solute concentrations can be made. The plots can be rotated, zoomed, scaled, translated, and printed. Axis values and the color bar can be customized. AVS FEHM output files are available for the following node data: material properties, liquid and vapor phase values, velocity and pressure values, temperature, saturation, concentrationand, and for the dua and dpdp models. The AVS output files from FEHM are written in an ASCII format that can be imported into AVS UCD graphics routines for viewing. 

 Additional information on the data found in the output files is given below. 

## Output file (```filen.out```)
 Information contained in the general output file is mostly self explanatory. The file starts with the code version, date, and time followed by the user input problem title. A summary of the I/O files used, macro control statements read, and array storage follow. Variable information for user specified nodes at user selected time intervals is written. The file ends with a summary of simulation time, number of time steps in the problem, the number of iterations taken, and total cpu time. 

## Write file (```filen.fin```)

 The write file contains the final values of pressure, temperature, saturation, and simulation time for the run. The final version of the file is generally written at the end of the simulation. This information is also written if the minimum user supplied time step has been realized and execution is terminated. The primary use of the write file is as a restart file. The write file contains the following: 

  
* Code version number, date, time 
* Problem title 
* Simulation time (days) 
* Gas flag [h20 (default), ngas, air] 
* Tracer flag [trac, ptrk, ntra (default - no output)] 
* Stress flag [strs, nstr (default - no output)] 
* Dpdp flag [dpdp, ndpd (default - no output)] 
* Dual flag [dual, ndua (default - no output)] 
* If the gas flag is ‘h20’ (neither air or ngas are set), followed by 
* Final temperature (oC) at each node 
* Final saturation (dimensionless) at each node 
* Final pressure (MPa) at each node 
* Or if ‘ngas’ flag is set, followed by 
* Final temperature (oC) at each node 
* Final saturation (dimensionless) at each node 
* Final pressure (MPa) at each node 
* Final capillary pressure (MPa) at each node 
* Or if ‘air’ flag is set, followed by 
* Final saturation (dimensionless) at each node 
* Final pressure (MPa) at each node 
* If fluxes have been specified in the rest macro (or for compatibility with older versions of the code if (ABS (PRNT_RST) = 2 in ptrk macro) 
* Label: ‘all fluxes’ or ‘liquid flux’ or ‘vapor flux’ 
* Number of mass flux values 
* Mass flux value (kg/s) for each connection of each node, starting with node 1. Note: mass flux values for the fracture domain are listed first followed by the mass flux values in the matrix domain. The mass flux between fracture and matrix elements are listed last. If both liquid and vapor fluxes as written liquid flux will be output first. 
 
 Otherwise: 
  
* Label: ‘no fluxes’ 
* If ‘trac’ flag is set followed by 
* Number of species 
* Species concentration (vapor or liquid, dimensionless) for each node for each species 
* Or if ‘ptrk’ flag is set followed by 
* Number of particles, final random number seed for transport calculations, final random number seed for particle release position (use by GoldSim) 
* Final node position for each particle. If the value is negative, the particle left the model domain at a fluid sink at that node. 
* Fractional time remaining at current node for each particle. 
* Multiplier to the plug flow residence time for each particle at the current node position, accounting for dispersion, sorption, and matrix diffusion effects. 
* Age for each particle, i.e. the time since the particle entered the system. However, if the particle has left the system, this value is the time that the particle left the system. 
* If the random number seed for transport calculations in the file is negative, the arrays for the fractional time remaining and the multiplier to the plug flow time have been omitted using the PRNT_RST = -1 or PRNT_RST = -2 option (see PRNT_RST description in the PTRK macro). A restart simulation using this input file will only approximate the behavior of particles since each particle will be assumed to have just entered the node. It is preferable to restart a particle tracking simulation using a file that contains the full restart information. 
* If strs (not implemented in this version) 
* If ‘dpdp’ or ‘dual’ flag is set 
 
 The above information includes dual porosity/dual permeability nodes. 

## Example of FEHM restart (```fin```) file using original format:

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V3</span><span class="o">.</span><span class="mi">1</span><span class="n">gf</span> <span class="mi">12</span><span class="o">-</span><span class="mi">02</span><span class="o">-</span><span class="mi">02</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">02</span><span class="o">/</span><span class="mi">14</span><span class="o">/</span><span class="mi">2012</span> <span class="mi">10</span><span class="p">:</span><span class="mi">33</span><span class="p">:</span><span class="mi">49</span>
<span class="n">Unsaturated</span> <span class="n">Diffusion</span> <span class="n">tests</span>
<span class="mf">5000.0000000000000</span>
<span class="n">ngas</span>
<span class="n">trac</span>
<span class="n">nstr</span>
<span class="n">ndpd</span>
<span class="n">ndua</span>
<span class="mf">34.99999999987494</span>      <span class="mf">34.99999999987494</span>      <span class="mf">29.99740954219060</span>      <span class="mf">29.99740954219060</span>
<span class="mf">24.99481908388880</span>      <span class="mf">24.99481908388880</span>      <span class="mf">19.99222863160355</span>      <span class="mf">19.99222863160355</span>
<span class="mf">14.99935303204482</span>      <span class="mf">14.99935303204482</span>      <span class="mf">10.00000000012507</span>      <span class="mf">10.00000000012507</span>
<span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span>
<span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1727371363921276</span>     <span class="mf">0.1727371363921281</span>
<span class="mf">0.4344871249926068</span>     <span class="mf">0.4344871249926068</span>     <span class="mf">0.7817833455822488</span>     <span class="mf">0.7817833455822516</span>
<span class="mf">0.1001154694602094</span>     <span class="mf">0.1001154694602094</span>     <span class="mf">0.1001154694628803</span>     <span class="mf">0.1001154694628803</span>
<span class="mf">0.1001154694707533</span>     <span class="mf">0.1001154694707533</span>     <span class="mf">0.1001154694901246</span>     <span class="mf">0.1001154694901246</span>
<span class="mf">0.1001154722096991</span>     <span class="mf">0.1001154722096991</span>     <span class="mf">0.1001154822144740</span>     <span class="mf">0.1001154822144740</span>
<span class="mf">0.9766094917448133E-01</span> <span class="mf">0.9766094917448133E-01</span> <span class="mf">0.9770095207050181E-01</span> <span class="mf">0.9770095207050181E-01</span>
<span class="mf">0.9774097492577727E-01</span> <span class="mf">0.9774097492577727E-01</span> <span class="mf">0.9778102503811041E-01</span> <span class="mf">0.9778102503811041E-01</span>
<span class="mf">0.9841762151617499E-01</span> <span class="mf">0.9841762151617499E-01</span> <span class="mf">0.9888735221221216E-01</span> <span class="mf">0.9888735221221216E-01</span>
<span class="n">no</span> <span class="n">fluxes</span>
<span class="mi">1</span>
  <span class="mf">1.040511468</span>  <span class="mf">1.040511468</span>  <span class="mf">1.023454397</span>  <span class="mf">1.023454397</span>
  <span class="mf">1.006402317</span>  <span class="mf">1.006402317</span>  <span class="mf">0.9893551455</span> <span class="mf">0.9893551455</span>
  <span class="mf">0.9701457197</span> <span class="mf">0.9701457197</span> <span class="mf">0.9516070487</span> <span class="mf">0.9516070487</span>
</pre></div>
</div>
 
 In FEHM Version 3.00 or newer, the format of the file was simplified to allow user selection of output parameters. The modified header contains: 
  
* Code version number, date, time 
* Problem title 
* Simulation time (days) 
* Number of nodes (n), type keyword (‘dual’, ‘dpdp’ or ‘nddp) to designate if dual porosity or double permeability were invoked or not for the simulation. 
 
 The remainder of the data is output following a ‘keyword’ header specifying the type of data to follow. n values will be output for each specified output parameter. Flux data will be output using the original format following the parameter output and preceding any particle tracking ‘ptrk’ or transport ‘trac’ output. Particle tracking output or transport output will be preceded by a ‘keyword’ header and use the same format as before. 

## Example of FEHM restart (```fin```) file using new format:

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V3</span><span class="o">.</span><span class="mi">1</span><span class="n">gf</span> <span class="mi">12</span><span class="o">-</span><span class="mi">02</span><span class="o">-</span><span class="mi">09</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">02</span><span class="o">/</span><span class="mi">09</span><span class="o">/</span><span class="mi">2012</span> <span class="mi">11</span><span class="p">:</span><span class="mi">48</span><span class="p">:</span><span class="mi">27</span>
<span class="n">Unsaturated</span> <span class="n">Diffusion</span> <span class="n">tests</span>
<span class="mf">5000.0000000000000</span>
<span class="mi">12</span> <span class="n">nddp</span>
<span class="n">temperature</span>
<span class="mf">34.99999999987494</span> <span class="mf">34.99999999987494</span> <span class="mf">29.99740954219060</span> <span class="mf">29.99740954219060</span>
<span class="mf">24.99481908388880</span> <span class="mf">24.99481908388880</span> <span class="mf">19.99222863160355</span> <span class="mf">19.99222863160355</span>
<span class="mf">14.99935303204482</span> <span class="mf">14.99935303204482</span> <span class="mf">10.00000000012507</span> <span class="mf">10.00000000012507</span>
<span class="n">saturation</span>
<span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span>
<span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1000000000000000E-98</span> <span class="mf">0.1727371363921276</span> <span class="mf">0.1727371363921281</span>
<span class="mf">0.4344871249926068</span> <span class="mf">0.4344871249926068</span> <span class="mf">0.7817833455822488</span> <span class="mf">0.7817833455822516</span>
<span class="n">pressure</span>
<span class="mf">0.1001154694602094</span> <span class="mf">0.1001154694602094</span> <span class="mf">0.1001154694628803</span> <span class="mf">0.1001154694628803</span>
<span class="mf">0.1001154694707533</span> <span class="mf">0.1001154694707533</span> <span class="mf">0.1001154694901246</span> <span class="mf">0.1001154694901246</span>
<span class="mf">0.1001154722096991</span> <span class="mf">0.1001154722096991</span> <span class="mf">0.1001154822144740</span> <span class="mf">0.1001154822144740</span>
<span class="n">gaspressure</span>
<span class="mf">0.9766094917448133E-01</span> <span class="mf">0.9766094917448133E-01</span> <span class="mf">0.9770095207050181E-01</span> <span class="mf">0.9770095207050181E-01</span>
<span class="mf">0.9774097492577727E-01</span> <span class="mf">0.9774097492577727E-01</span> <span class="mf">0.9778102503811041E-01</span> <span class="mf">0.9778102503811041E-01</span>
<span class="mf">0.9841762151617499E-01</span> <span class="mf">0.9841762151617499E-01</span> <span class="mf">0.9888735221221216E-01</span> <span class="mf">0.9888735221221216E-01</span>
<span class="n">no</span> <span class="n">fluxes</span>
<span class="n">trac</span>
<span class="mi">1</span>
<span class="mf">1.040511468</span> <span class="mf">1.040511468</span> <span class="mf">1.023454397</span> <span class="mf">1.023454397</span>
<span class="mf">1.006402317</span> <span class="mf">1.006402317</span> <span class="mf">0.9893551455</span> <span class="mf">0.9893551455</span>
<span class="mf">0.9701457197</span> <span class="mf">0.9701457197</span> <span class="mf">0.9516070487</span> <span class="mf">0.9516070487</span>
</pre></div>
</div>


## History plot file (```filen.his```)
 The history plot file contains the following: 
  
* Code version number, date, time 
* Problem title 
* Gas flag (‘ngas’, ‘airw’, or blank) 
* Tracer flag (‘trac’ or blank) 
* Stress flag (‘strs’ or blank) 
* Number of nodes for which data are output 
* Node number and X, Y, and Z coordinate (m) of each node for which data are output 
* ‘headings’ 
* Depending on problem flags, 2 lines with field descriptors, Case 1 (default): 
* “node flow enthalpy(Mj/kg) flow(kg/s) temperature (deg C) total pressure (Mpa)” 
* “capillary pressure (Mpa) saturation (kg/kg)” 
 
 or Case 2 (hydraulic head): 
  
* “node flow enthalpy(Mj/kg) flow(kg/s) temperature (deg C) hydraulic head (m)” 
* “total pressure (Mpa) saturation (kg/kg)” 
 
 or Case 3 (ngas): 
  
* “node flow enthalpy(Mj/kg) flow(kg/s) temperature (deg C) total pressure (Mpa)” 
* “air pressure (Mpa) capillary pressure (Mpa) saturation (kg/kg) relative humidity” 
 
 And for each time step 
  
* Time (days) followed by 
 
 For Case 1: 
  
* Node number, Energy source (MJ/s), Source strength (kg/s), Temperature (oC), Total pressure (MPa), Capillary pressure (MPa), and Saturation (dimensionless) for each specified output node. 
 
 For Case 2 (hydraulic head): 
  
* Node number, Energy source (MJ/s), Source strength (kg/s), Temperature (oC), Hydraulic head (m), Total pressure (MPa), and Saturation (dimensionless) for each specified output node. 
 
 For Case 3 (ngas): 
  
* Node number, Energy source (MJ/s), Source strength (kg/s), Temperature (oC), Total pressure (MPa), Air pressure (MPa), Capillary pressure (MPa), Saturation (dimensionless), and Relative humidity for each specified output node 
 
## Example of history output file, ```filen.his```

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V2</span><span class="o">.</span><span class="mi">30</span><span class="n">sun</span> <span class="mi">05</span><span class="o">-</span><span class="mi">04</span><span class="o">-</span><span class="mi">19</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">04</span><span class="o">/</span><span class="mi">20</span><span class="o">/</span><span class="mi">2005</span> <span class="mi">08</span><span class="p">:</span><span class="mi">45</span><span class="p">:</span><span class="mi">05</span>
<span class="o">*****</span> <span class="mi">2</span><span class="o">-</span><span class="n">D</span> <span class="n">Heat</span> <span class="n">Conduction</span> <span class="n">Model</span> <span class="o">*****</span>
<span class="mi">1</span>
<span class="mi">111</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="n">headings</span>
<span class="n">node</span> <span class="n">flow</span> <span class="n">enthalpy</span><span class="p">(</span><span class="n">Mj</span><span class="o">/</span><span class="n">kg</span><span class="p">)</span> <span class="n">flow</span><span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">s</span><span class="p">)</span> <span class="n">temperature</span><span class="p">(</span><span class="n">deg</span> <span class="n">C</span><span class="p">)</span> <span class="n">total</span> <span class="n">pressure</span><span class="p">(</span><span class="n">Mpa</span><span class="p">)</span>
<span class="n">capillary</span> <span class="n">pressure</span><span class="p">(</span><span class="n">Mpa</span><span class="p">)</span> <span class="n">saturation</span><span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">kg</span><span class="p">)</span>
<span class="mf">0.0E+0</span>
<span class="mi">111</span> <span class="mf">0.100000000E-19</span> <span class="mf">0.00000000</span> <span class="mf">200.000000</span> <span class="mf">10.0000000</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span>
<span class="mf">5.0E-3</span>
<span class="mi">111</span> <span class="mf">0.100000000E-19</span> <span class="mf">0.00000000</span> <span class="mf">199.999999</span> <span class="mf">10.0000000</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mf">4.000049999999937</span>
<span class="mi">111</span> <span class="mf">0.100000000E-19</span> <span class="mf">0.00000000</span> <span class="mf">100.183607</span> <span class="mf">10.0000000</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span>
<span class="o">-</span><span class="mf">4.000049999999937</span>
<span class="mi">111</span> <span class="mf">0.100000000E-19</span> <span class="mf">0.00000000</span> <span class="mf">100.183607</span> <span class="mf">10.0000000</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span>
</pre></div>
</div>


## Alternate History plot files (```filen.his```,```filen_param.his```)
 The history plot file (```filen.his```) contains the following: 
  
* Code version number, date, time 
* Problem title 
* Gas flag (‘ngas’, ‘airw’, or blank) 
* Tracer flag (‘trac’ or blank) 
* Stress flag (‘strs’ or blank) 
* List of parameters written to individual history files (possible paramenters are: pressure, temperature, head, saturation, flow, enthalpy, humidity, zone flux, and global) 
* Number of nodes for which data are output 
* Node number and X, Y, and Z coordinate (m) of each node for which data are output 
 
## Example of alternate history output file, ```filen.his```

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V2</span><span class="o">.</span><span class="mi">30</span><span class="n">sun</span> <span class="mi">05</span><span class="o">-</span><span class="mi">04</span><span class="o">-</span><span class="mi">20</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">04</span><span class="o">/</span><span class="mi">20</span><span class="o">/</span><span class="mi">2005</span> <span class="mi">08</span><span class="p">:</span><span class="mi">50</span><span class="p">:</span><span class="mi">56</span>
<span class="o">*****</span> <span class="mi">2</span><span class="o">-</span><span class="n">D</span> <span class="n">Heat</span> <span class="n">Conduction</span> <span class="n">Model</span> <span class="o">*****</span>
<span class="n">Parameters</span> <span class="n">written</span> <span class="n">to</span> <span class="n">individual</span> <span class="n">history</span> <span class="n">files</span><span class="p">:</span>
<span class="n">temperature</span>
<span class="k">for</span> <span class="n">the</span> <span class="n">following</span> <span class="n">nodes</span><span class="p">:</span>
<span class="mi">1</span>
<span class="mi">111</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span>
</pre></div>
</div>

 If zones for output are specified: 
  
* Number of zones over which output is averaged 
 
 And for each zone: 
  
* Number of nodes in the zone and zone number followed by a list of nodes in the zone. 
 
## Example of alternate history output file including zones, ```filen.his```

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V2</span><span class="o">.</span><span class="mi">30</span><span class="n">sun</span> <span class="mi">05</span><span class="o">-</span><span class="mi">04</span><span class="o">-</span><span class="mi">20</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">04</span><span class="o">/</span><span class="mi">20</span><span class="o">/</span><span class="mi">2005</span> <span class="mi">09</span><span class="p">:</span><span class="mi">10</span><span class="p">:</span><span class="mi">21</span>
<span class="o">*****</span> <span class="mi">2</span><span class="o">-</span><span class="n">D</span> <span class="n">Heat</span> <span class="n">Conduction</span> <span class="n">Model</span> <span class="o">*****</span>
<span class="n">Parameters</span> <span class="n">written</span> <span class="n">to</span> <span class="n">individual</span> <span class="n">history</span> <span class="n">files</span><span class="p">:</span>
<span class="n">temperature</span>
<span class="k">for</span> <span class="n">the</span> <span class="n">following</span> <span class="n">nodes</span> <span class="ow">and</span> <span class="n">zones</span><span class="p">:</span>
<span class="mi">1</span>
<span class="mi">111</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span> <span class="mf">0.00000000</span>
<span class="n">Number</span> <span class="n">of</span> <span class="n">averaged</span> <span class="n">output</span> <span class="n">zones</span><span class="p">:</span> <span class="mi">2</span>
<span class="mi">121</span> <span class="n">Nodes</span> <span class="n">averaged</span> <span class="ow">in</span> <span class="n">Zone</span> <span class="o">-</span><span class="mi">1</span>
<span class="mi">0000001</span> <span class="mi">0000002</span> <span class="mi">0000003</span> <span class="mi">0000004</span> <span class="mi">0000005</span> <span class="mi">0000006</span> <span class="mi">0000007</span> <span class="mi">0000008</span> <span class="mi">0000009</span> <span class="mi">0000010</span>
<span class="mi">0000011</span> <span class="mi">0000012</span> <span class="mi">0000013</span> <span class="mi">0000014</span> <span class="mi">0000015</span> <span class="mi">0000016</span> <span class="mi">0000017</span> <span class="mi">0000018</span> <span class="mi">0000019</span> <span class="mi">0000020</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mi">0000121</span>
<span class="mi">11</span> <span class="n">Nodes</span> <span class="n">averaged</span> <span class="ow">in</span> <span class="n">Zone</span> <span class="o">-</span><span class="mi">2</span>
<span class="mi">0000001</span> <span class="mi">0000012</span> <span class="mi">0000023</span> <span class="mi">0000034</span> <span class="mi">0000045</span> <span class="mi">0000056</span> <span class="mi">0000067</span> <span class="mi">0000078</span> <span class="mi">0000089</span> <span class="mi">0000100</span>
<span class="mi">0000111</span>
</pre></div>
</div>
 The history plot parameter files (filen_param.his) contain the following: 
  
* Code version number, date, time 
* Problem title 
* Blank line 
* Output parameter title and units 
* Heading: Time (units) Nodes: Node number 1 … Node number n 
 
 -or- if zones are specified 
  
* Heading: Time (units) Nodes: Node number 1 … Node number n Zones: Zone number 1 … Zode number n 
 
 And for each time step (time units may be seconds, days, or years as specified in hist macro) 
  
* Time (units) followed by parameter value for each specified node and zone. 
 
## Example of alternate history output file, ```filen_temp.his```

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V2</span><span class="o">.</span><span class="mi">30</span><span class="n">sun</span> <span class="mi">05</span><span class="o">-</span><span class="mi">04</span><span class="o">-</span><span class="mi">20</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">04</span><span class="o">/</span><span class="mi">20</span><span class="o">/</span><span class="mi">20</span> <span class="mi">09</span><span class="p">:</span><span class="mi">39</span><span class="p">:</span><span class="mi">26</span>
<span class="o">*****</span> <span class="mi">2</span><span class="o">-</span><span class="n">D</span> <span class="n">Heat</span> <span class="n">Conduction</span> <span class="n">Model</span> <span class="o">*****</span>
<span class="n">Temperature</span> <span class="p">(</span><span class="n">C</span><span class="p">)</span>
<span class="n">Time</span> <span class="p">(</span><span class="n">seconds</span><span class="p">)</span> <span class="n">Nodes</span><span class="p">:</span> <span class="mi">111</span>
<span class="mf">0.0E+0</span> <span class="mf">200.0</span>
<span class="mf">432.0</span> <span class="mf">199.99999943712558</span>
<span class="mf">864.0</span> <span class="mf">199.99999502157842</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mf">345599.9999999945</span> <span class="mf">100.18362295565089</span>
<span class="mf">345604.31999999453</span> <span class="mf">100.18360733099975</span>
</pre></div>
</div>

## Example of alternate history output file including zones, filen_temp.his

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V2</span><span class="o">.</span><span class="mi">30</span><span class="n">sun</span> <span class="mi">05</span><span class="o">-</span><span class="mi">04</span><span class="o">-</span><span class="mi">20</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">04</span><span class="o">/</span><span class="mi">20</span><span class="o">/</span><span class="mi">20</span> <span class="mi">09</span><span class="p">:</span><span class="mi">10</span><span class="p">:</span><span class="mi">21</span>
<span class="o">*****</span> <span class="mi">2</span><span class="o">-</span><span class="n">D</span> <span class="n">Heat</span> <span class="n">Conduction</span> <span class="n">Model</span> <span class="o">*****</span>
<span class="n">Temperature</span> <span class="p">(</span><span class="n">C</span><span class="p">)</span>
<span class="n">Time</span> <span class="p">(</span><span class="n">seconds</span><span class="p">)</span> <span class="n">Nodes</span><span class="p">:</span> <span class="mi">111</span> <span class="n">Zones</span><span class="p">:</span> <span class="o">-</span><span class="mi">1</span> <span class="o">-</span><span class="mi">2</span>
<span class="mf">0.0E+0</span> <span class="mf">200.0</span> <span class="mf">199.99999999999985</span> <span class="mf">200.0</span>
<span class="mf">432.0</span> <span class="mf">199.99999943712558</span> <span class="mf">187.43603443819043</span> <span class="mf">193.49769265774958</span>
<span class="mf">864.0</span> <span class="mf">199.99999502157842</span> <span class="mf">184.97995686212098</span> <span class="mf">192.16893175137713</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mf">345599.9999999945</span> <span class="mf">100.18362295565089</span> <span class="mf">100.07411373139031</span> <span class="mf">100.11665754330934</span>
<span class="mf">345604.31999999453</span> <span class="mf">100.18360733099975</span> <span class="mf">100.07410742498257</span> <span class="mf">100.11664761680842</span>
</pre></div>
</div>

 The history plot global parameter file (filen_param.his) contains the following: 
  
* Code version number, date, time 
* Problem title 
* Blank line 
 
 Depending on output selected the following possible headings: 
  
* mass / energy: Time (days) Total mass in system (kg) Total mass of steam in system (kg) Water discharge (kg) Water input (kg) Total water discharge (kg) Total water input (kg) Net (kg) water discharge Total enthalpy in system (MJ) Enthalpy discharge (MJ) Enthalpy input (MJ) Total enthalpy discharge (MJ) Total enthalpy input (MJ) Net (MJ) enthalpy discharge 
* water / air: Time (days) Total water in system (kg) Total mass of steam in system (kg) Water discharge (kg) Water input (kg) Net (kg) water discharge Total water discharge (kg) Total water input (kg) Total air in system (kg) Air discharge (kg) Air input (kg) Total air discharge (kg) Total air input kg (kg/s)Net (kg) air discharge 
* mass / water only (no steam): Time (days) Total water in system (kg) Water discharge (kg) Water input (kg) Total water discharge (kg) Total water input (kg) Net (kg) water discharge 
* mass / water only (steam): Time (days) Total mass in system (kg) Total mass of steam in system (kg) Water discharge (kg) Water input (kg) Total water discharge (kg) Total water input (kg) Net (kg) water discharge 
* air only: Time (days) Total air in system (kg) Air discharge (kg) Air input (kg) Total air discharge (kg) Total air input kg (kg/s) Net (kg) air discharge 
* energy only: Time (days) Total enthalpy in system (MJ) Enthalpy discharge (MJ) Enthalpy input (MJ) Total enthalpy discharge MJ Total enthalpy input (MJ)Net (MJ) enthalpy discharge 
 
 And for each time step 
  
* Time (days) followed by selected global parameter values. 
 
## Solute plot file (```filen.trc```)
 Solute data is output for the same nodes used for the history plot file. The solute plot file contains: 
  
* Code version number, date, time 
* Problem title 
* Number of nodes for which data are output 
* Node number and X, Y, and Z coordinate (m) of each node for which data are output 
* Number of different species/components for tracer solution, Number of liquid components, Number of immobile components, Number of vapor components, and Number of aqueous complexes 
 
 and for each time step and each species 
  
* Time (days), species number followed by 
* Species concentration (dimensionless) for each specified output node. 
 
 When particle tracking is used, the concentration can be output in several different forms (number of particles, number per fluid mass, or number per total volume). The choice of which form to use is input in the ptrk macro. 

## Example of solute data history plot file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V2</span><span class="o">.</span><span class="mi">30</span><span class="n">sun</span> <span class="mi">05</span><span class="o">-</span><span class="mi">04</span><span class="o">-</span><span class="mi">06</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">04</span><span class="o">/</span><span class="mi">06</span><span class="o">/</span><span class="mi">2005</span> <span class="mi">19</span><span class="p">:</span><span class="mi">58</span><span class="p">:</span><span class="mi">25</span>
<span class="n">Check</span> <span class="n">of</span> <span class="n">FEHMN</span> <span class="n">against</span> <span class="n">SORBEQ</span><span class="p">,</span> <span class="n">All</span> <span class="n">isotherms</span>
<span class="mi">1</span>
<span class="mi">201</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.100000000E+01</span>
<span class="mi">5</span> <span class="mi">5</span> <span class="mi">0</span> <span class="mi">0</span> <span class="mi">0</span>
<span class="mf">1.550709E-4</span> <span class="mi">1</span> <span class="n">species</span> <span class="c1">#001</span>
<span class="mf">4.855185258201169E-29</span>
<span class="mf">1.550709E-4</span> <span class="mi">2</span> <span class="n">species</span> <span class="c1">#002</span>
<span class="mf">3.7325204274237394E-53</span>
<span class="mf">1.550709E-4</span> <span class="mi">3</span> <span class="n">species</span> <span class="c1">#003</span>
<span class="mf">3.7325204274237394E-53</span>
<span class="mf">1.550709E-4</span> <span class="mi">4</span> <span class="n">species</span> <span class="c1">#004</span>
<span class="mf">1.0E-90</span>
<span class="mf">1.550709E-4</span> <span class="mi">5</span> <span class="n">species</span> <span class="c1">#005</span>
<span class="mf">1.0E-90</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mf">1.3E-3</span> <span class="mi">1</span> <span class="n">species</span> <span class="c1">#001</span>
<span class="mf">0.999994799250939</span>
<span class="mf">1.3E-3</span> <span class="mi">2</span> <span class="n">species</span> <span class="c1">#002</span>
<span class="mf">0.9999246650715027</span>
<span class="mf">1.3E-3</span> <span class="mi">3</span> <span class="n">species</span> <span class="c1">#003</span>
<span class="mf">0.9999947753362215</span>
<span class="mf">1.3E-3</span> <span class="mi">4</span> <span class="n">species</span> <span class="c1">#004</span>
<span class="mf">0.9999947668674314</span>
<span class="mf">1.3E-3</span> <span class="mi">5</span> <span class="n">species</span> <span class="c1">#005</span>
<span class="mf">0.9999947643072644</span>
</pre></div>
</div>

## Contour plot file (```filen.con```)
 The contour plot file contains: 
  
* Code version number, date, time 
* Problem title 
* Tracer (‘trac’) solution or blank 
* Stress (‘strs’) solution or blank 
* Number of nodes for which data are output 
* X, Y, and Z coordinate (m) of each node for which data are output 
* Number of nodes per element, total number of elements 
* Nodal connectivity information for each node of each element 
* X, Y, Z permeability (<span class="math notranslate nohighlight">\(m^2\)</span>) for each node 
* X, Y, Z thermal conductivity <span class="math notranslate nohighlight">\(\left(\frac{W}{m \cdot K}\right)\)</span> for each node 
* Porosity, Rock specific heat <span class="math notranslate nohighlight">\(\left(\frac{MJ}{kg \cdot K}\right)\)</span>, Capillary pressure (<span class="math notranslate nohighlight">\(MPa\)</span>) for each node 
* Number of degrees of freedom per node for the current problem, Direction of gravity in problem, Value of gravity 
 
 If tracer solution is present: 
  
* Number of species 
 
 and for each specified time: 
  
* Time (days), injection phase (Š 0 liquid, &lt; 0 vapor) followed by 
 
 If injection phase is liquid: 
  
* Liquid transmissibility / density, Liquid density (kg/m3), Pressure - Capillary Pressure (MPa), Temperature (oC) 
 
 and if tracer solution is present: 
  
* Species concentration of liquid phase 
 
 Or if injection phase is vapor: 
  
* Vapor transmissibility / density, Vapor density (kg/m3), Pressure (MPa), Temperature (oC) 
 
 and if tracer solution is present: 
  
* Species concentration of vapor phase. 
 
## Contour plot file for dual or dpdp (```filen.dp```)
 The contour plot file for dual or dpdp contains the same information as the regular contour plot file only the parameter values are for the dual porosity / dual permeability nodes. 

## Stiffness matrix data (```filen.stor```)
 The stiffness matrix data file is used to store the finite element coefficients for each node. It eliminates the need for the code to recompute the coefficients for subsequent runs. It contains the following: 
  
* Code version number, date, time 
* Problem title 
* Number of storage locations needed to store geometric input types, Number of nodes, Size of connectivity array 
* Volume associated with each node 
* Nodal connectivity information for each connection 
* Position of geometric coefficient for each connection 
* Diagonal position in connectivity array for each node 
* Finite element geometric coefficient for each storage location 
* If stress solution is enabled 
* Finite element geometric coefficient for each storage location for the stress module. 
 
## Input check file (```filen.chk```) 
 This file contains a summary of input information that may be of use in debugging or memory management of the code. The positions of maximum and minimum values of input parameters and derived quantities are given. Also provided is an analysis of array storage requirements. 

## Submodel output file (```filen.subbc```)
 The submodel output file contains “flow” macro data that represents boundary conditions for an extracted submodel (i.e., the output will use the format of the “flow” input macro). The file contains: 
  
* Heading: “flow Boundary Conditions Output:”, code version number, date, time for each submodel node, if boundary type is head or pressure, 
* Node number, Node number, 1, Head (m) or Pressure (MPa), Impedance parameter, #, X coordinate, Y coordinate, Z cooordinate 
 
 or if boundary type is flux: 
  
* Node number, Node number, 1, Flux (kg/s), 0.0d00, #, X coordinate, Y coordinate, Z cooordinate 
 
 A blank line to signal end of flow macro input followed by the file termination macro ```stop```. 
 An example is provided with <a class="reference external" href="MacroSubm.html">subm</a> input on. 

## Error output file (```fehmn.err```)
 This file contains the code version number, date, and time followed by any error or warning messages issued by the code during a run. 

## Multiple simulations script files (```fehmn.pre```,```fehmn.post```)
 The multiple simulations script file fehmn.pre contains UNIX shell script style instructions for pre-processing input data, while the script file fehmn.post contains UNIX shell script style instructions for post-processing data. 

## PEST output files (```filen.pest```,```filen.pest1```)
 The PEST output file is used to output data in a format suitable for use by the Parameter Estimation Program (PEST) (Watermark Computing, 1994). The first file (filen.pest) contains: 
  
* Heading: “PEST Output:”, code version number, date, time 
* First parameter label: “pressures” or “heads” 
* node number and pressure (MPa) or head (ft) for each specified output node 
* Second parameter label: “saturations” 
* node number and saturation for each specified output node 
* Third parameter label: “temperatures” 
* node number and temperature (oC) for each specified output node 
* Fourth parameter label: “permeabilities” 
* node number and x, y, and z permeability (m2) for each specified output node 
* Heading:      “Total Flux (kg/s) Leaving Zone (flxz macro option)” 
* “Zone Number  Source  Sink    Net     Boundary” 
* zone number, source flux, sink flux, net flux, boundary flux 
 
 The second file (filen.pest1) contains: 
  
* Heading: “PEST Output:”, code version number, date, time 
* Parameter label: “pressures” or “heads” 
* node number, relative permeability model used, pressure (MPa) or head (ft), saturation, and temperature (oC) for each specified output node 
 
## Particle statistics file (```filen.ptrk```)
 The data found in the “.ptrk” file was previously reported in the general output file. From version 2.25 of the code and forward that data will only be reported in the optional, “.ptrk” file unless a coupled GoldSim-FEHM simulation is being run (see the PRNT_RST flag, see Control Statement <a class="reference external" href="MacroMptr.html">mptr</a>). In addition, the user has the option of selecting which statistics parameters are reported. The default is to report all statistics parameters. 
 For the default case, the ```.ptrk``` file contains: 
  
* Title line: <strong>TITLE=”V1=Number Having Entered System, V2=Number Currently In System, V3=Number Having Left System, V4=Number Having Decayed, V5=Number Having Been Filtered, V6=Number That Left This Time”</strong> 
* Header line: <strong>VARIABLES=”Time (days)” “Sp001 V1” “Sp001 V2” “Sp001 V3” “Sp001 V4” “Sp001 V5” “Sp001 V6” … “Spnnn V1” “Spnnn V2” “Spnnn V3” “Spnnn V4” “Spnnn V5” “Spnnn V6”</strong> 
 
 A heading is output for each variable (V1 to V6) for each species in the model and nnn is the total number of species in the simulation. The header line is followed by (for each output time step) the simulation time (days), and the six output variables for each species modeled (See the first example below). When the user selects a subset of the statistics parameters the header line and data will only contain those variables that have been selected for output (See the second example below). 

## Example of default “.ptk” particle statistics file.

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">TITLE</span><span class="o">=</span><span class="s2">&quot;V1=Number Having Entered System, V2=Number Currently In System, V3=Number Having Left</span>
<span class="n">System</span><span class="p">,</span> <span class="n">V4</span><span class="o">=</span><span class="n">Number</span> <span class="n">Having</span> <span class="n">Decayed</span><span class="p">,</span> <span class="n">V5</span><span class="o">=</span><span class="n">Number</span> <span class="n">Having</span> <span class="n">Been</span> <span class="n">Filtered</span><span class="p">,</span> <span class="n">V6</span><span class="o">=</span><span class="n">Number</span> <span class="n">That</span> <span class="n">Left</span> <span class="n">This</span> <span class="n">Time</span><span class="s2">&quot;</span>
<span class="n">VARIABLES</span><span class="o">=</span><span class="s2">&quot;Time (days)&quot;</span> <span class="s2">&quot;Sp001 V1&quot;</span> <span class="s2">&quot;Sp001 V2&quot;</span> <span class="s2">&quot;Sp001 V3&quot;</span> <span class="s2">&quot;Sp001 V4&quot;</span> <span class="s2">&quot;Sp001 V5&quot;</span> <span class="s2">&quot;Sp001 V6&quot;</span>
<span class="s2">&quot;Sp002 V1&quot;</span> <span class="s2">&quot;Sp002 V2&quot;</span> <span class="s2">&quot;Sp002 V3&quot;</span> <span class="s2">&quot;Sp002 V4&quot;</span> <span class="s2">&quot;Sp002 V5&quot;</span> <span class="s2">&quot;Sp002 V6&quot;</span>
<span class="mf">365.25000000000</span> <span class="mi">18760</span> <span class="mi">18686</span> <span class="mi">74</span> <span class="mi">0</span> <span class="mi">0</span> <span class="mi">74</span> <span class="mi">18760</span> <span class="mi">18670</span> <span class="mi">90</span> <span class="mi">0</span> <span class="mi">0</span> <span class="mi">90</span>
<span class="mf">1278.3750000000</span> <span class="mi">18760</span> <span class="mi">17860</span> <span class="mi">900</span> <span class="mi">0</span> <span class="mi">0</span> <span class="mi">826</span> <span class="mi">18760</span> <span class="mi">17895</span> <span class="mi">865</span> <span class="mi">0</span> <span class="mi">0</span> <span class="mi">775</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mf">7270752.1054688</span> <span class="mi">18760</span> <span class="mi">932</span> <span class="mi">16710</span> <span class="mi">1118</span> <span class="mi">0</span> <span class="mi">6</span> <span class="mi">18760</span> <span class="mi">912</span> <span class="mi">16692</span> <span class="mi">1156</span> <span class="mi">0</span> <span class="mi">14</span>
<span class="mf">7305000.0000000</span> <span class="mi">18760</span> <span class="mi">929</span> <span class="mi">16712</span> <span class="mi">1119</span> <span class="mi">0</span> <span class="mi">2</span> <span class="mi">18760</span> <span class="mi">905</span> <span class="mi">16695</span> <span class="mi">1160</span> <span class="mi">0</span> <span class="mi">3</span>
</pre></div>
</div>

## Example of “.ptk” particle statistics file with five output variables selected.

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">TITLE</span><span class="o">=</span><span class="s2">&quot;V1=Number Having Entered System, V2=Number Currently In System, V3=Number Having Left</span>
<span class="n">System</span><span class="p">,</span> <span class="n">V4</span><span class="o">=</span><span class="n">Number</span> <span class="n">Having</span> <span class="n">Decayed</span><span class="p">,</span> <span class="n">V5</span><span class="o">=</span><span class="n">Number</span> <span class="n">Having</span> <span class="n">Been</span> <span class="n">Filtered</span><span class="p">,</span> <span class="n">V6</span><span class="o">=</span><span class="n">Number</span> <span class="n">That</span> <span class="n">Left</span> <span class="n">This</span> <span class="n">Time</span><span class="s2">&quot;</span>
<span class="n">VARIABLES</span><span class="o">=</span><span class="s2">&quot;Time (days)&quot;</span> <span class="s2">&quot;Sp001 V1&quot;</span> <span class="s2">&quot;Sp001 V2&quot;</span> <span class="s2">&quot;Sp001 V3&quot;</span> <span class="s2">&quot;Sp001 V4&quot;</span> <span class="s2">&quot;Sp001 V6&quot;</span> <span class="s2">&quot;Sp002 V1&quot;</span>
<span class="s2">&quot;Sp002 V2&quot;</span> <span class="s2">&quot;Sp002 V3&quot;</span> <span class="s2">&quot;Sp002 V4&quot;</span> <span class="s2">&quot;Sp002 V6&quot;</span>
<span class="mf">365.25000000000</span> <span class="mi">18760</span> <span class="mi">18686</span> <span class="mi">74</span> <span class="mi">0</span> <span class="mi">74</span> <span class="mi">18760</span> <span class="mi">18670</span> <span class="mi">90</span> <span class="mi">0</span> <span class="mi">90</span>
<span class="mf">1278.3750000000</span> <span class="mi">18760</span> <span class="mi">17860</span> <span class="mi">900</span> <span class="mi">0</span> <span class="mi">826</span> <span class="mi">18760</span> <span class="mi">17895</span> <span class="mi">865</span> <span class="mi">0</span> <span class="mi">775</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mf">7270752.1054688</span> <span class="mi">18760</span> <span class="mi">932</span> <span class="mi">16710</span> <span class="mi">1118</span> <span class="mi">6</span> <span class="mi">18760</span> <span class="mi">912</span> <span class="mi">16692</span> <span class="mi">1156</span> <span class="mi">14</span>
<span class="mf">7305000.0000000</span> <span class="mi">18760</span> <span class="mi">929</span> <span class="mi">16712</span> <span class="mi">1119</span> <span class="mi">2</span> <span class="mi">18760</span> <span class="mi">905</span> <span class="mi">16695</span> <span class="mi">1160</span> <span class="mi">3</span>
</pre></div>
</div>

## Mass Output from GoldSim Particle Tracking Simulation
 To provide a simplified method for tracking solute mass from a FEHM/GoldSim coupled simulation an optional output file may be written that contains cumulative mass output (mg/l) (see the PRNT_RST flag!,, <a class="reference external" href="MacroMptr.html">mptr</a>). The data found in ```FEHM_GSM_Mass_balance.txt``` has a format similar to that used for the particle statistics output (see above) and contains: 
  
* Title line: <strong>TITLE=”V1=Mass Having Entered System, V2=Mass Currently In System, V3=Mass Having Left System, V4=Mass Having Decayed, V5=Mass Having Been Filtered, V6=Mass Having Decayed Outside The UZ, V7=Filtered Mass Having Decayed”</strong> 
* Header line: <strong>VARIABLES=”Time (years)” “Sp001 V1” “Sp001 V2” “Sp001 V3” “Sp001 V4” “Sp001 V5” “Sp001 V6” “Sp001 V7”… “Spnnn V1” “Spnnn V2” “Spnnn V3” “Spnnn V4” “Spnnn V5” “Spnnn V6” “Spnnn V7”</strong> 
 
 A heading is output for each variable (V1 to V7) for each species in the model and nnn is the total number of species in the simulation. The header line is followed by (for each output time step) the simulation time (days), and the seven output variables for each species modeled. 

## Particle Exit Locations and Count Output
 To facilitate use of particle tracking simulation statistics an optional output file may be written that contains particle exit locations and count (see the PRNT_RST flag, section, see <a class="reference external" href="MacroMptr.html">mptr</a>). The data found in the “.ptrk_fin” file can also be extracted from the “.fin” files through use of post-processors and/or file editors. 
 The “.ptrk_fin” file output file contains: 
  
* Header line: <strong>VARIABLES = “Node” “X” “Y” “Z” “Num_exited” “Zone”</strong> 
 
 Followed by (for each node where particles have exited the system) the node number, the X, Y, and Z coordinates of the node (m), the number of particles that exited at that node, and the number of the zone (if defined) that contains that node. 

## Streamline particle tracking output files (filen.sptr1, filen.sptr2, filen.sptr3)
 The streamline particle tracking output files contain information generated during a streamline particle tracking simulation. Depending on output options selected (macro sptr) zero, one, two or three output files are generated. 
 When option IPRT ≥ 1, the first file (filen.sptr1) contains: 
  
* Code version number, date, time 
* Problem title 
* “days=”, Current time of streamline particle tracking iteration (days), streamline particle tracking timestep number for current iteration, and 
 
 For each particle: 
  
* Particle number, X coordinate of particle, Y coordinate of particle, Z coordinate of particle, Element or node number where the particle is located 
 
 When option IPRTO ≥ 1, the second file (filen.sptr2) contains: 
  
* Code version number, date, time 
* Problem title 
* Heading: “particle_number x y z time porosity saturation permeability rock_density pressure temperature zone old_node new_node” (Note that the heading only includes property titles for the default properties or those properties specified by keyword.) 
 
 For each particle: 
  
* Particle number, X coordinate of particle, Y coordinate of particle, Z coordinate of particle, Current time that the particle has reached (days), Property value of the unit the particle is residing in for each specified keyword [in the following order, if specified: porosity, saturation, permeability (m2), density (kg/m3), pressure (MPa), temperature (oC), zone number], Element or node number where the particle is located, Previous node number 
 
 If option “zbtc” is invoked, the third file (filen.sptr3) contains: 
  
* Code version number, date, time 
* Problem title 
* Heading: ” Time (days) Zone1 Particles …” 
* Time (days), Cumulative number of particles that have arrived at each specified zone for breakthrough curve output 
 
 or when option ‘alt’ is specified: 
  
* Code version number, date, time 
* Problem title 
* Heading: Time (days) Particle# ID Zone Node 
* Followed by breakthrough time, particle number, particle ID, breakthrough zone, break through node for each particle that reaches the breakthrough zone. 
 
 Or when option ‘alt xyz’ is specified: 
  
* Heading: x(m) y(m) z(m) Time (days) Particle # ID Zone Node 
* Followed by x, y, z coordinate where particle entered breakthrough zone, breakthrough time, particle number, breakthrough zone, break through node for each particle that reaches the breakthrough zone. 
 
## Example of default “.sptr3” file.

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V3</span><span class="o">.</span><span class="mi">1</span><span class="n">win32</span> <span class="mi">12</span><span class="o">-</span><span class="mi">02</span><span class="o">-</span><span class="mi">02</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">02</span><span class="o">/</span><span class="mi">02</span><span class="o">/</span><span class="mi">2012</span> <span class="mi">10</span><span class="p">:</span><span class="mi">22</span><span class="p">:</span><span class="mi">34</span>
<span class="o">***</span> <span class="n">Validation1</span> <span class="n">Test</span> <span class="n">Problem</span><span class="p">:</span> <span class="mi">3</span><span class="o">-</span><span class="n">D</span> <span class="n">Homogeneous</span> <span class="n">Flow</span> <span class="ow">and</span> <span class="n">Transport</span> <span class="o">***</span>
<span class="n">Time</span> <span class="p">(</span><span class="n">days</span><span class="p">)</span> <span class="n">Zone1</span> <span class="n">Particles</span> <span class="o">.</span> <span class="o">.</span> <span class="o">.</span>
<span class="mf">208000.0</span> <span class="mi">0</span>
<span class="mf">210000.0</span> <span class="mi">40</span>
</pre></div>
</div>

## Example of “.sptr3” file generated using option “alt”.

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V3</span><span class="o">.</span><span class="mi">1</span><span class="n">win32</span> <span class="mi">12</span><span class="o">-</span><span class="mi">02</span><span class="o">-</span><span class="mi">02</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">02</span><span class="o">/</span><span class="mi">02</span><span class="o">/</span><span class="mi">2012</span> <span class="mi">10</span><span class="p">:</span><span class="mi">22</span><span class="p">:</span><span class="mi">41</span>
<span class="o">***</span> <span class="n">Validation1</span> <span class="n">Test</span> <span class="n">Problem</span><span class="p">:</span> <span class="mi">3</span><span class="o">-</span><span class="n">D</span> <span class="n">Homogeneous</span> <span class="n">Flow</span> <span class="ow">and</span> <span class="n">Transport</span> <span class="o">***</span>
<span class="n">Time</span> <span class="p">(</span><span class="n">days</span><span class="p">)</span> <span class="n">Particle</span><span class="c1"># ID Zone Node</span>
<span class="mf">209391.24384033</span> <span class="mi">1</span> <span class="mi">1888</span> <span class="mi">5</span> <span class="mi">1938</span>
<span class="mf">209391.24384033</span> <span class="mi">2</span> <span class="mi">1888</span> <span class="mi">5</span> <span class="mi">1938</span>
<span class="o">...</span>
<span class="mf">209391.24381900</span> <span class="mi">39</span> <span class="mi">51256</span> <span class="mi">5</span> <span class="mi">51306</span>
<span class="mf">209391.24381900</span> <span class="mi">40</span> <span class="mi">51256</span> <span class="mi">5</span> <span class="mi">51306</span>
</pre></div>
</div>

## Example of “.sptr3” file generated using option “alt xyz”.

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">FEHM</span> <span class="n">V3</span><span class="o">.</span><span class="mi">1</span><span class="n">win32</span> <span class="mi">12</span><span class="o">-</span><span class="mi">02</span><span class="o">-</span><span class="mi">02</span> <span class="n">QA</span><span class="p">:</span><span class="n">NA</span> <span class="mi">02</span><span class="o">/</span><span class="mi">02</span><span class="o">/</span><span class="mi">2012</span> <span class="mi">10</span><span class="p">:</span><span class="mi">22</span><span class="p">:</span><span class="mi">48</span>
<span class="o">***</span> <span class="n">Validation1</span> <span class="n">Test</span> <span class="n">Problem</span><span class="p">:</span> <span class="mi">3</span><span class="o">-</span><span class="n">D</span> <span class="n">Homogeneous</span> <span class="n">Flow</span> <span class="ow">and</span> <span class="n">Transport</span> <span class="o">***</span>
<span class="n">x</span><span class="p">(</span><span class="n">m</span><span class="p">)</span> <span class="n">y</span><span class="p">(</span><span class="n">m</span><span class="p">)</span> <span class="n">z</span><span class="p">(</span><span class="n">m</span><span class="p">)</span> <span class="n">Time</span> <span class="p">(</span><span class="n">days</span><span class="p">)</span> <span class="n">Particle</span><span class="c1"># ID Zone Node</span>
<span class="mf">19800.0000</span> <span class="mf">0.100000000</span> <span class="o">-</span><span class="mf">5.50000000</span> <span class="mf">209391.24384033</span> <span class="mi">1</span> <span class="mi">1888</span> <span class="mi">5</span> <span class="mi">1938</span>
<span class="mf">19800.0000</span> <span class="mf">0.200000000</span> <span class="o">-</span><span class="mf">5.50000000</span> <span class="mf">209391.24384033</span> <span class="mi">2</span> <span class="mi">1888</span> <span class="mi">5</span> <span class="mi">1938</span>
<span class="o">...</span>
<span class="mf">19800.0000</span> <span class="o">-</span><span class="mf">2800.90000</span> <span class="o">-</span><span class="mf">200.000000</span> <span class="mf">209391.24381900</span> <span class="mi">39</span> <span class="mi">51256</span> <span class="mi">5</span> <span class="mi">51306</span>
<span class="mf">19800.0000</span> <span class="o">-</span><span class="mf">2801.00000</span> <span class="o">-</span><span class="mf">200.000000</span> <span class="mf">209391.24381900</span> <span class="mi">40</span> <span class="mi">51256</span> <span class="mi">5</span> <span class="mi">51306</span>
</pre></div>
</div>

## Log output file (filen.avs_log)
 The log output file is identical for AVS, AVS Express, Surfer and Tecplot. It contains: 
  
* Code version number, date 
* AVS log identifier 
* Problem title 
 
 and for each specified time: 
  
* Output file prefix, Call number, and Time. The output time units will correspond to those selected in the cont macro. 
 
## Example of contour log output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># FEHM V3.1gf 12-02-09 QA:NA 02/09/2012</span>
<span class="c1"># LOG AVS OUTPUT</span>
<span class="c1"># Unsaturated Diffusion tests</span>
<span class="c1"># Root filename Output Time (days)</span>
<span class="n">output</span><span class="o">/</span><span class="n">box</span><span class="o">.</span><span class="mi">00001</span> <span class="mf">0.00000000</span>
<span class="n">output</span><span class="o">/</span><span class="n">box</span><span class="o">.</span><span class="mi">00002</span> <span class="mf">1001.68908</span>
<span class="n">output</span><span class="o">/</span><span class="n">box</span><span class="o">.</span><span class="mi">00003</span> <span class="mf">2002.77258</span>
<span class="n">output</span><span class="o">/</span><span class="n">box</span><span class="o">.</span><span class="mi">00004</span> <span class="mf">3002.77258</span>
<span class="n">output</span><span class="o">/</span><span class="n">box</span><span class="o">.</span><span class="mi">00005</span> <span class="mf">4003.12657</span>
<span class="n">output</span><span class="o">/</span><span class="n">box</span><span class="o">.</span><span class="mi">00006</span> <span class="mf">5000.00000</span>
</pre></div>
</div>

## AVS header output files (filen.type_head)
 The AVS ASCII (formatted) header files are identical for AVS and [AVS Express output](#AVS-header). The data types, “mat”, “sca”, “vec” or “con”, are described below. The header files contain: 
 20 lines of text with information about the FEHM AVS output files. The text is followed by a one line AVS UCD file header containing: 
  
* number of nodes 
* number of cells 
* number of data components for the nodes 
* number of data components for the cells (currently 0) 
* number of data components for the model (currently 0) 
 
## <a id="AVS-header"></a>Example of AVS header output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># FEHM V3.1gf 12-02-02 QA:NA 02/14/2012</span>
<span class="c1"># AVS UNSTRUCTURED CELL DATA (UCD) FROM FEHM</span>
<span class="c1"># Unsaturated Diffusion tests</span>
<span class="c1"># *****************************************************</span>
<span class="c1"># To prepare files for input to avs one must</span>
<span class="c1"># concatinate header/geometry/node_value files.</span>
<span class="c1"># For example, if your FEHM input file was fe.dat,</span>
<span class="c1"># headers are fe10001_sca_head fe10001_vec_head, ...,</span>
<span class="c1"># mesh geometry will be in fe10001_geo,</span>
<span class="c1"># field output will be in fe10001_sca_node,</span>
<span class="c1"># fe10001_vec_node, fe10001_con_dual_node</span>
<span class="c1">#</span>
<span class="c1"># A UCD input file can be produced using</span>
<span class="c1"># cat fe10001_sca_head fe10001_geo fe10001_sca_node &gt;</span>
<span class="c1"># fe10001_sca_node.inp</span>
<span class="c1">#</span>
<span class="c1"># The UNIX foreach command is useful for processing</span>
<span class="c1"># multiple files. Also use the shell script fehm2avs</span>
<span class="c1"># to perform automatic processing of all output.</span>
<span class="c1"># *****************************************************</span>
<span class="mi">0000000012</span> <span class="mi">5</span> <span class="mi">5</span> <span class="mi">0</span> <span class="mi">0</span>
</pre></div>
</div>

## Geometry output files (filen.geo, filen_grid.dat)
 Geometry data will be output when keyword “geom” or “grid” are included in the output variable list. For AVS the geometry data is output in a separate file when keyword ‘geom’ is used. For Tecplot, output of the geometry data depends on the input keyword, ‘geom’ or ‘grid’. If keyword “geom’ is used the geometry data is contained in the first contour file for each type of data requested. If keyword ‘grid’ is used the geometry data is output to a separate Tecplot “grid” file. 
 The ASCII (formatted) geometry file for AVS contains the following: 
  
* Node id and X, Y, Z coordinates for each node 
* Cell id, Material id, Cell type, and the list of Cell vertices 
* The ASCII (formatted) geometry file for AVS Express contains one additional line of data at the beginning of the file, followed by the data specified above: 
* Number of nodes, Number of elements, 0, 0, 0 
 
## Example of AVS geometry output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="mi">0000000001</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000002</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000003</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.200000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000004</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.200000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000005</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.400000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000006</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.400000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000007</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.600000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000008</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.600000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000009</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.800000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000010</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.800000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000011</span> <span class="mf">0.000000000E+00</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000012</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.100000000E+01</span> <span class="mf">0.000000000E+00</span>
<span class="mi">0000000001</span> <span class="mi">1</span> <span class="n">quad</span> <span class="mi">1</span> <span class="mi">2</span> <span class="mi">4</span> <span class="mi">3</span>
<span class="mi">0000000002</span> <span class="mi">1</span> <span class="n">quad</span> <span class="mi">3</span> <span class="mi">4</span> <span class="mi">6</span> <span class="mi">5</span>
<span class="mi">0000000003</span> <span class="mi">1</span> <span class="n">quad</span> <span class="mi">5</span> <span class="mi">6</span> <span class="mi">8</span> <span class="mi">7</span>
<span class="mi">0000000004</span> <span class="mi">1</span> <span class="n">quad</span> <span class="mi">7</span> <span class="mi">8</span> <span class="mi">10</span> <span class="mi">9</span>
<span class="mi">0000000005</span> <span class="mi">1</span> <span class="n">quad</span> <span class="mi">9</span> <span class="mi">10</span> <span class="mi">12</span> <span class="mi">11</span>
</pre></div>
</div>
 The Tecplot “grid” file contains the following (note that only the grid coordinates used are output): 
  
* Header line with problem title 
* Filetype header 
* Variable header (coordinates used) 
* Zone header with grid specification and type (N = number of nodes, E = number of elements) 
 
 Followed by: 
  
* “N” node coordinate sets 
* Element connectivity for the grid for “E” elements. 
 
## Example of Tecplot grid output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">TITLE</span> <span class="o">=</span> <span class="s2">&quot;Unsaturated Diffusion tests&quot;</span>
<span class="n">FILETYPE</span> <span class="o">=</span> <span class="s2">&quot;GRID&quot;</span>
<span class="n">VARIABLES</span> <span class="o">=</span> <span class="s2">&quot;X coordinate (m)&quot;</span> <span class="s2">&quot;Y coordinate (m)&quot;</span>
<span class="n">ZONE</span> <span class="n">T</span> <span class="o">=</span> <span class="s2">&quot;GRID&quot;</span><span class="p">,</span> <span class="n">N</span> <span class="o">=</span> <span class="mi">12</span><span class="p">,</span> <span class="n">E</span> <span class="o">=</span> <span class="mi">5</span><span class="p">,</span> <span class="n">DATAPACKING</span> <span class="o">=</span> <span class="n">POINT</span><span class="p">,</span> <span class="n">ZONETYPE</span>
<span class="o">********************************************</span>
<span class="n">FEQUADRILATERAL</span><span class="p">,</span> <span class="n">STRANDID</span> <span class="mi">0</span><span class="p">,</span> <span class="n">SOLUTIONTIME</span> <span class="mf">0.</span>
<span class="o">********************************************</span>
<span class="mf">0.000000000E+00</span> <span class="mf">0.000000000E+00</span>
<span class="mf">0.100000000E+01</span> <span class="mf">0.000000000E+00</span>
<span class="mf">0.000000000E+00</span> <span class="mf">0.200000000E+00</span>
<span class="mf">0.100000000E+01</span> <span class="mf">0.200000000E+00</span>
<span class="mf">0.000000000E+00</span> <span class="mf">0.400000000E+00</span>
<span class="mf">0.100000000E+01</span> <span class="mf">0.400000000E+00</span>
<span class="mf">0.000000000E+00</span> <span class="mf">0.600000000E+00</span>
<span class="mf">0.100000000E+01</span> <span class="mf">0.600000000E+00</span>
<span class="mf">0.000000000E+00</span> <span class="mf">0.800000000E+00</span>
<span class="mf">0.100000000E+01</span> <span class="mf">0.800000000E+00</span>
<span class="mf">0.000000000E+00</span> <span class="mf">0.100000000E+01</span>
<span class="mf">0.100000000E+01</span> <span class="mf">0.100000000E+01</span>
<span class="mi">1</span> <span class="mi">2</span> <span class="mi">4</span> <span class="mi">3</span>
<span class="mi">3</span> <span class="mi">4</span> <span class="mi">6</span> <span class="mi">5</span>
<span class="mi">5</span> <span class="mi">6</span> <span class="mi">8</span> <span class="mi">7</span>
<span class="mi">7</span> <span class="mi">8</span> <span class="mi">10</span> <span class="mi">9</span>
<span class="mi">9</span> <span class="mi">10</span> <span class="mi">12</span> <span class="mi">11</span>
</pre></div>
</div>
 Geometry data when contained in the normal tecplot data files uses shared variables. Coordinates and connectivity are output only in the first file. The data files contain the following: 
  
* The Tecplot “grid” file contains the following: 
* Header line with code version number, date, time and problem title 
* Variable header 
* Zone header with time, grid specification and type 
 
 Followed by: 
  
* “N” node coordinate, node number and datasets 
* Element connectivity for the grid for “E” elements 
 
## Example of Tecplot data output file with geometry data included

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">TITLE</span> <span class="o">=</span> <span class="s2">&quot;FEHM V3.1gf 12-02-09 QA:NA 02/09/2012 11:48:26 Unsaturated</span>
<span class="n">Diffusion</span> <span class="n">tests</span><span class="s2">&quot;</span>
<span class="n">VARIABLES</span> <span class="o">=</span> <span class="s2">&quot;X (m)&quot;</span> <span class="s2">&quot;Y (m)&quot;</span> <span class="s2">&quot;Node&quot;</span> <span class="s2">&quot;Vapor_Species_001&quot;</span>
<span class="n">ZONE</span> <span class="n">T</span> <span class="o">=</span> <span class="s2">&quot;Simulation time 0.00000000 days&quot;</span><span class="p">,</span> <span class="n">N</span> <span class="o">=</span> <span class="mi">12</span><span class="p">,</span> <span class="n">E</span> <span class="o">=</span> <span class="mi">5</span><span class="p">,</span>
<span class="n">DATAPACKING</span> <span class="o">=</span> <span class="n">POINT</span><span class="p">,</span> <span class="n">ZONETYPE</span> <span class="o">=</span> <span class="n">FEQUADRILATERAL</span><span class="p">,</span> <span class="n">N</span> <span class="o">=</span> <span class="mi">12</span><span class="p">,</span> <span class="n">E</span> <span class="o">=</span> <span class="mi">5</span><span class="p">,</span>
<span class="n">DATAPACKING</span> <span class="o">=</span> <span class="n">POINT</span><span class="p">,</span> <span class="n">ZONETYPE</span> <span class="o">=</span> <span class="n">FEQUADRILATERAL</span>
<span class="mf">0.00000000</span> <span class="mf">0.00000000</span> <span class="mi">0000000001</span> <span class="mf">1.00000000</span>
<span class="mf">1.00000000</span> <span class="mf">0.00000000</span> <span class="mi">0000000002</span> <span class="mf">1.00000000</span>
<span class="mf">0.00000000</span> <span class="mf">0.200000000</span> <span class="mi">0000000003</span> <span class="mf">1.00000000</span>
<span class="mf">1.00000000</span> <span class="mf">0.200000000</span> <span class="mi">0000000004</span> <span class="mf">1.00000000</span>
<span class="mf">0.00000000</span> <span class="mf">0.400000000</span> <span class="mi">0000000005</span> <span class="mf">1.00000000</span>
<span class="mf">1.00000000</span> <span class="mf">0.400000000</span> <span class="mi">0000000006</span> <span class="mf">1.00000000</span>
<span class="mf">0.00000000</span> <span class="mf">0.600000000</span> <span class="mi">0000000007</span> <span class="mf">1.00000000</span>
<span class="mf">1.00000000</span> <span class="mf">0.600000000</span> <span class="mi">0000000008</span> <span class="mf">1.00000000</span>
<span class="mf">0.00000000</span> <span class="mf">0.800000000</span> <span class="mi">0000000009</span> <span class="mf">1.00000000</span>
<span class="mf">1.00000000</span> <span class="mf">0.800000000</span> <span class="mi">0000000010</span> <span class="mf">1.00000000</span>
<span class="mf">0.00000000</span> <span class="mf">1.00000000</span> <span class="mi">0000000011</span> <span class="mf">1.00000000</span>
<span class="mf">1.00000000</span> <span class="mf">1.00000000</span> <span class="mi">0000000012</span> <span class="mf">1.00000000</span>
<span class="mi">1</span> <span class="mi">2</span> <span class="mi">4</span> <span class="mi">3</span>
<span class="mi">3</span> <span class="mi">4</span> <span class="mi">6</span> <span class="mi">5</span>
<span class="mi">5</span> <span class="mi">6</span> <span class="mi">8</span> <span class="mi">7</span>
<span class="mi">7</span> <span class="mi">8</span> <span class="mi">10</span> <span class="mi">9</span>
<span class="mi">9</span> <span class="mi">10</span> <span class="mi">12</span> <span class="mi">11</span>
</pre></div>
</div>

## Contour data output files (filen.number_type_node.suffix)
 All the ASCII (formatted) node data files for AVS (suffix ‘avs’) contain the following headers: 
  
* Number of data components and size of each component 
* A label/unit string for each data component 
 
 followed by, for each node: 
  
* the associated node data (described by data type below). 
* All the ASCII (formatted) node data files for AVS Express (suffix ‘avsx’) contain the following headers, on a single line delimited by ” : ” : 
* Current simulation time (with format “nodes at time days”) 
* A label/unit string for each data component 
 
 followed by, for each node: 
  
* the associated node data (described by data type below), delimited by ” : “. 
 
 All of the node data files for Surfer (suffix ‘csv’) contain a single header line containing: 
  
* A label/unit string for each data component separated by “,” 
 
 followed by for each node: 
  
* the associated node data (described by data type below), delimited by ” , “. 
 
 All the node data files for Tecplot (suffix ‘dat’) contain the following headers: 
  
* Header line with code version number, date, time and problem title 
* Filetype header when keyword “grid” is used 
* Variable header (the variable header is only output to the first data file unless keyword ‘grid’ is used) 
* Zone header with time (not output for material property files if keyword ‘grid’ is used) 
 
 or 
  
* Zone header with time and grid specification and type if keyword geom is used 
 
 followed by for each node: 
  
* the associated node data (described by data type below). 
 
 The dual or dpdp values for each of these fields will be written to a file with “dual” in the file name. 

| Contour File Content |   Data Type Designation  | Output Parameters* |
|:---------------|:-----------------------|:--------------------|
| Material properties | mat, mat_dual | Permeability in each active direction (m^2) |
|   |   | Thermal conductivity in each active direction <span class="math notranslate nohighlight">\(\left( \frac{W}{m \cdot K} \right)\)</span> |
|   |   | Porosity |
|   |   | Rock bulk density (kg/m^3) |
|   |   | Rock specific heat <span class="math notranslate nohighlight">\(\left( \frac{MJ}{kg \cdot K} \right)\)</span> |
|   |   | Capillary pressure (MPa) |
|   |   | Relative permeability model |
|   |   | Capillary pressure model. |
|   |   | Note: Output to the material properities file is dependent on the simulation being performed. |
| Scalar parameters | sca, sca_dual | Zone number |
|   |   | Pressure (MPa) - liquid, vapor, capillary |
|   |   | Temperature (oC) |
|   |   | Saturation |
|   |   | CO2 (Water volume fraction, Liquid CO2 fraction, Gaseous CO2 fraction, Dissolved CO2 mass fraction, Phase state of CO2 |
|   |   | Head (m) |
|   |   | Porosity |
|   |   | Density (kg/m^3) - liquid, vapor |
|   |   | Permeability in each active direction (m^2) |
|   |   | Source (kg/s) |
|   |   | Mass Flux (kg/s) |
|   |   | Volume weighted Mass Flux (kg/s/m^3) |
|   |   | Displacement (m) for each specified direction |
|   |   | Stress (MPa) for each specified direction |
|   |   | Volume Strain |
| Vector parameters | vec, vec_dual | Volume Flux (m^3/m^2/s) - liquid, vapor |
| Solute concentrations | con, con_dual | Species concentration (moles/kg) |
| * Output parameters are dependent upon the simulation being performed and keywords specified in the cont macro. | | |


## Example of AVS material properties output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="mi">11</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="n">Permeability</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span><span class="p">,</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span>
<span class="n">Permeability</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span><span class="p">,</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span>
<span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span><span class="p">,</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span>
<span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span><span class="p">,</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span>
<span class="n">Porosity</span><span class="p">,</span> <span class="p">(</span><span class="n">non</span> <span class="n">dim</span><span class="p">)</span>
<span class="n">Rock</span> <span class="n">bulk</span> <span class="n">density</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">),</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">)</span>
<span class="n">Rock</span> <span class="n">specific</span> <span class="n">heat</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">),</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">)</span>
<span class="n">Capillary</span> <span class="n">pressure</span> <span class="p">(</span><span class="n">MPa</span><span class="p">),</span> <span class="p">(</span><span class="n">MPa</span><span class="p">)</span>
<span class="n">Relative</span> <span class="n">permeability</span> <span class="n">model</span><span class="p">,</span> <span class="p">(</span><span class="n">flag</span><span class="p">)</span>
<span class="n">Capillary</span> <span class="n">pressure</span> <span class="n">model</span><span class="p">,</span> <span class="p">(</span><span class="n">flag</span><span class="p">)</span>
<span class="mi">0000000001</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000002</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000003</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000004</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000005</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000006</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000007</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000008</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000009</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000010</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000011</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000012</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
</pre></div>
</div>

## Example of AVS Express material properties output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">node</span> <span class="p">:</span> <span class="n">Permeability</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span> <span class="p">:</span> <span class="n">Permeability</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span> <span class="p">:</span> <span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span>
<span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span> <span class="p">:</span> <span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span> <span class="p">:</span> <span class="n">Porosity</span> <span class="p">:</span> <span class="n">Rock</span> <span class="n">bulk</span> <span class="n">density</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">)</span>
<span class="p">:</span> <span class="n">Rock</span> <span class="n">specific</span> <span class="n">heat</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="p">:</span> <span class="n">Capillary</span> <span class="n">pressure</span> <span class="p">(</span><span class="n">MPa</span><span class="p">)</span> <span class="p">:</span> <span class="n">Relative</span> <span class="n">permeability</span> <span class="n">model</span>
<span class="p">:</span> <span class="n">Capillary</span> <span class="n">pressure</span> <span class="n">model</span>
<span class="mi">0000000001</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000002</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000003</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000004</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000005</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000006</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000007</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000008</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000009</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000010</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000011</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
<span class="mi">0000000012</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-12</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">1.000000E-05</span> <span class="p">:</span> <span class="mf">0.800000</span>
<span class="p">:</span> <span class="mf">1000.00</span> <span class="p">:</span> <span class="mf">0.00000</span> <span class="p">:</span> <span class="mi">1</span> <span class="p">:</span> <span class="mi">1</span>
</pre></div>
</div>

## Example of Surfer material properties output file

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">node</span> <span class="p">,</span> <span class="n">Permeability</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span> <span class="p">,</span> <span class="n">Permeability</span> <span class="p">(</span><span class="n">m</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span> <span class="p">,</span> <span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span>
<span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span> <span class="p">,</span> <span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span> <span class="p">,</span> <span class="n">Porosity</span> <span class="p">,</span> <span class="n">Rock</span> <span class="n">bulk</span> <span class="n">density</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">)</span>
<span class="p">,</span> <span class="n">Rock</span> <span class="n">specific</span> <span class="n">heat</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="p">,</span> <span class="n">Capillary</span> <span class="n">pressure</span> <span class="p">(</span><span class="n">MPa</span><span class="p">)</span> <span class="p">,</span> <span class="n">Relative</span> <span class="n">permeability</span> <span class="n">model</span>
<span class="p">,</span> <span class="n">Capillary</span> <span class="n">pressure</span> <span class="n">model</span>
<span class="mi">0000000001</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000002</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000003</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000004</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000005</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000006</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000007</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000008</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000009</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000010</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000011</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
<span class="mi">0000000012</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-12</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">1.000000E-05</span><span class="p">,</span> <span class="mf">0.800000</span>
<span class="p">,</span> <span class="mf">1000.00</span><span class="p">,</span> <span class="mf">0.00000</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span>
</pre></div>
</div>

## Example of Tecplot material properties output file without geometry keyword

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">TITLE</span> <span class="o">=</span> <span class="s2">&quot;FEHM V3.1gf 12-02-15 QA:NA 02/15/2012 11:29:32 Unsaturated Diffusion tests&quot;</span>
<span class="n">VARIABLES</span> <span class="o">=</span> <span class="s2">&quot;node&quot;</span> <span class="s2">&quot;Permeability (m**2) in X&quot;</span> <span class="s2">&quot;Permeability (m**2) in Y&quot;</span> <span class="s2">&quot;Thermal</span>
<span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span><span class="s2">&quot; &quot;</span><span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span><span class="s2">&quot; &quot;</span><span class="n">Porosity</span><span class="s2">&quot; &quot;</span><span class="n">Rock</span> <span class="n">bulk</span>
<span class="n">density</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Rock</span> <span class="n">specific</span> <span class="n">heat</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Capillary</span> <span class="n">pressure</span> <span class="p">(</span><span class="n">MPa</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Relative</span>
<span class="n">permeability</span> <span class="n">model</span><span class="s2">&quot; &quot;</span><span class="n">Capillary</span> <span class="n">pressure</span> <span class="n">model</span><span class="s2">&quot;</span>
<span class="mi">0000000001</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000002</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000003</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000004</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000005</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000006</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000007</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000008</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000009</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000010</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000011</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000012</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
</pre></div>
</div>

## Example of Tecplot material properties output file with ‘grid’ keyword

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">TITLE</span> <span class="o">=</span> <span class="s2">&quot;FEHM V3.1gf 12-02-15 QA:NA 02/15/2012 07:58:21 Unsaturated Diffusion tests&quot;</span>
<span class="n">FILETYPE</span> <span class="o">=</span> <span class="s2">&quot;SOLUTION&quot;</span>
<span class="n">VARIABLES</span> <span class="o">=</span> <span class="s2">&quot;node&quot;</span> <span class="s2">&quot;Permeability (m**2) in X&quot;</span> <span class="s2">&quot;Permeability (m**2) in Y&quot;</span> <span class="s2">&quot;Thermal</span>
<span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">X</span><span class="s2">&quot; &quot;</span><span class="n">Thermal</span> <span class="n">Conductivity</span> <span class="p">(</span><span class="n">W</span><span class="o">/</span><span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span><span class="s2">&quot; &quot;</span><span class="n">Porosity</span><span class="s2">&quot; &quot;</span><span class="n">Rock</span> <span class="n">bulk</span>
<span class="n">density</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Rock</span> <span class="n">specific</span> <span class="n">heat</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Capillary</span> <span class="n">pressure</span> <span class="p">(</span><span class="n">MPa</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Relative</span>
<span class="n">permeability</span> <span class="n">model</span><span class="s2">&quot; &quot;</span><span class="n">Capillary</span> <span class="n">pressure</span> <span class="n">model</span><span class="s2">&quot;</span>
<span class="mi">0000000001</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000002</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000003</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000004</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000005</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000006</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000007</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000008</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000009</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000010</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000011</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">0000000012</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span>
<span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
</pre></div>
</div>

## Example of Tecplot material properties output file with ‘geom’ keyword

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">VARIABLES</span> <span class="o">=</span> <span class="s2">&quot;X coordinate (m)&quot;</span> <span class="s2">&quot;Y coordinate (m)&quot;</span> <span class="s2">&quot;node&quot;</span> <span class="s2">&quot;Permeability (m**2) in X&quot;</span>
<span class="s2">&quot;Permeability (m**2) in Y&quot;</span> <span class="s2">&quot;Thermal Conductivity (W/m*K) in X&quot;</span> <span class="s2">&quot;Thermal Conductivity (W/</span>
<span class="n">m</span><span class="o">*</span><span class="n">K</span><span class="p">)</span> <span class="ow">in</span> <span class="n">Y</span><span class="s2">&quot; &quot;</span><span class="n">Porosity</span><span class="s2">&quot; &quot;</span><span class="n">Rock</span> <span class="n">bulk</span> <span class="n">density</span> <span class="p">(</span><span class="n">kg</span><span class="o">/</span><span class="n">m</span><span class="o">**</span><span class="mi">3</span><span class="p">)</span><span class="s2">&quot; &quot;</span><span class="n">Rock</span> <span class="n">specific</span> <span class="n">heat</span> <span class="p">(</span><span class="n">MJ</span><span class="o">/</span><span class="n">kg</span><span class="o">*</span><span class="n">K</span><span class="p">)</span><span class="s2">&quot;</span>
<span class="s2">&quot;Capillary pressure (MPa)&quot;</span> <span class="s2">&quot;Relative permeability model&quot;</span> <span class="s2">&quot;Capillary pressure model&quot;</span>
<span class="n">ZONE</span> <span class="n">T</span> <span class="o">=</span> <span class="s2">&quot;Material properties&quot;</span><span class="p">,</span> <span class="n">N</span> <span class="o">=</span> <span class="mi">12</span><span class="p">,</span> <span class="n">E</span> <span class="o">=</span> <span class="mi">5</span><span class="p">,</span> <span class="n">DATAPACKING</span> <span class="o">=</span> <span class="n">POINT</span><span class="p">,</span> <span class="n">ZONETYPE</span>
<span class="o">***************</span>
<span class="n">FEQUADRILATERAL</span>
<span class="o">***************</span>
<span class="mf">0.00000000</span> <span class="mf">0.00000000</span> <span class="mi">0000000001</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">1.00000000</span> <span class="mf">0.00000000</span> <span class="mi">0000000002</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">0.00000000</span> <span class="mf">0.200000000</span> <span class="mi">0000000003</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">1.00000000</span> <span class="mf">0.200000000</span> <span class="mi">0000000004</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">0.00000000</span> <span class="mf">0.400000000</span> <span class="mi">0000000005</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">1.00000000</span> <span class="mf">0.400000000</span> <span class="mi">0000000006</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">0.00000000</span> <span class="mf">0.600000000</span> <span class="mi">0000000007</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">1.00000000</span> <span class="mf">0.600000000</span> <span class="mi">0000000008</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">0.00000000</span> <span class="mf">0.800000000</span> <span class="mi">0000000009</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">1.00000000</span> <span class="mf">0.800000000</span> <span class="mi">0000000010</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">0.00000000</span> <span class="mf">1.00000000</span> <span class="mi">0000000011</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mf">1.00000000</span> <span class="mf">1.00000000</span> <span class="mi">0000000012</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000E-12</span> <span class="mf">1.000000</span><span class="n">E</span><span class="o">-</span>
<span class="mi">05</span> <span class="mf">1.000000E-05</span> <span class="mf">0.800000</span> <span class="mf">1000.00</span> <span class="mf">1.000000E-03</span> <span class="mf">0.00000</span> <span class="mi">1</span> <span class="mi">1</span>
<span class="mi">1</span> <span class="mi">2</span> <span class="mi">4</span> <span class="mi">3</span>
<span class="mi">3</span> <span class="mi">4</span> <span class="mi">6</span> <span class="mi">5</span>
<span class="mi">5</span> <span class="mi">6</span> <span class="mi">8</span> <span class="mi">7</span>
<span class="mi">7</span> <span class="mi">8</span> <span class="mi">10</span> <span class="mi">9</span>
<span class="mi">9</span> <span class="mi">10</span> <span class="mi">12</span> <span class="mi">11</span>
</pre></div>
</div>

## SURFER and TECPLOT contour output files with specified ‘zone’
 The content of the contour files generated when the ‘zone’ keyword is used in the cont macro is the same as that for the regular output with the following exceptions: 

1. Geometry keywords, ‘geom’ and ‘grid’, are ignored; 
1. Data is output only for the nodes in the specified zones; 
1. For “surfer”, a separate file is written for each output zone and the file names generated include the output zone number (using 4 digits, e.g., 0001); 
1. For tecplot, the simulation time is written into a text string, and the zone headers include only the zone number, and output is separated by zone. 

  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">mptr</span></code><a class="headerlink" href="#mptr" title="Permalink to this headline">¶</a></h1>
<p>Multiple species ingrowth particle tracking. Note that data for each numbered group must be input. The other input is optional.</p>
<ul class="simple">
<li>Group 1 -     NSPECI, MAXLAYERS, MAX_PARTICLES, RIPFEHM, MAX1D</li>
<li>Group 2 -     POUT, PRNT_RST<ul>
<li>or when PRNT_RST ≥ 20, selected output parameters</li>
</ul>
</li>
<li>Group 2 -     POUT, PRNT_RST, PRNT_VARNUM ( 1 … 6)<ul>
<li>Optional keyword “tcurve” is input to indicate that transfer function curves should be input to model matrix diffusion. It is followed by NUMPARAMS and TFILENAME.<ul>
<li>KEYWORD</li>
<li>NUMPARAMS, FFMAX</li>
<li>TFILENAME</li>
</ul>
</li>
<li>Optional keyword “zptr” designates zones for breakthrough curves will be defined. It is followed by IPZONE and IDZONE.<ul>
<li>KEYWORD ‘zptr’</li>
<li>IPZONE</li>
<li>IDZONE(I) I = 1 to IPZONE</li>
</ul>
</li>
</ul>
</li>
<li>Group 3 -     RSEED, RSEED_RELEASE<ul>
<li>Optional keyword “wtri” is input to indicate a water table rise calculation should be performed. It is followed by WATER_TABLE. For GoldSim the water table rise calculation is controlled by passing a new water table elevation to the code during the simulation and the keyword is not required.<ul>
<li>KEYWORD ‘wtri’</li>
<li>WATER_TABLE</li>
</ul>
</li>
</ul>
</li>
<li>Group 4 -     DAYCS, DAYCF, DAYHF, DAYHS<ul>
<li>An optional, flexible input structure involving the assignment of transport parameters is implemented in the particle tracking input to allow multiple realization simulations to use different parameters for each realization. The user invokes this option using the keyword “file” before Group 5, followed by the name of the file that the transport parameters reside in. The applicable transport parameters are defined in Group 5 and Group 9.<ul>
<li>KEYWORD ‘file’</li>
<li>PFILENAME</li>
</ul>
</li>
<li>The structure of the alternate parameter file is:<ul>
<li>NINPUTS</li>
<li>PTRPARAM(I) I=1 to NINPUTS</li>
<li>.</li>
<li>.</li>
<li>.</li>
</ul>
</li>
<li>[a line of parameters is present for each realization]</li>
</ul>
</li>
</ul>
<p>The method for assigning a given value of the particle tracking parameter (PTRPARAM) to a specific transport parameter, defined in Group 5 or Group 9, is discussed below. There are an arbitrary number of input lines each representing a given realization of parameter values. In a multiple-realization scenario, the code enters the input file for each realization, and for this input, reads down the corresponding number of lines to obtain the parameters for that realization. For example, for realization number 10, the code reads down to the 10th line of data (line 11 in the file) and uses those parameter values.</p>
<p>Once these parameters are read in for a given realization, they must be assigned to specific transport parameters. This is done in the following way in Group 5 or Group 9. If any of the inputs other than TRANSFLAG are negative, the code takes the absolute value of the number and interprets it as the column number from which to assign the transport parameter. For example, if DIFFMFL = -5, then the diffusion coefficient is the fifth input number in the PTRPARAM array. In this way, any of the transport parameter inputs can be assigned through the alternate input file rather than the input line in mptr. It should be noted that for the colloid diversity model, only K_REV need be negative to indicate values should be read from the parameter file, if K_REV is negative then all five parameters are read from the file, otherwise the equation parameters will be read from the mptr macro. This is to accommodate the fact that the SLOPE_KF may have negative values.</p>
<p>Group 5 is used to define models in which identical transport parameters are assumed to apply. Group 5 data are read until a blank line is encountered. The model number ITRC is incremented by 1 each time a line is read. Model parameters defined in Group 5 are assigned to nodes or zones using Group6.</p>
<p>Optional keyword “afm” indicates the Active Fracture Model input for determining fracture spacing should be used. Optional keyword “dfree” is input to indicate that a free water diffusion coefficient and tortuosity will be entered instead of the molecular diffusion coefficient.</p>
<blockquote>
<div><p>KEYWORD ‘afm’</p>
<p>KEYWORD ‘dfre’</p>
</div></blockquote>
<ul class="simple">
<li>Group 5 -     TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), APERTUR(ITRC), MATRIX_POR(ITRC)<ul>
<li>or when ‘afm’ is implemented:</li>
<li>Group 5 -   TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), APERTUR(ITRC), MATRIX_POR(ITRC), SRESIDUAL(ITRC), GAMMA_AFM(ITRC)</li>
</ul>
</li>
<li>Group 6 -     JA, JB, JC, ITRC (JA, JB, JC - defined on <a class="reference external" href="Macro20058.html">See JA, JB, JC, PROP1, PROP2, …</a>)</li>
</ul>
<p>The following groups (Group 7 - 12) are repeated for each species.</p>
<ul class="simple">
<li>Group 7 -     ITH_SPECI, TRAK_TYPE, HALF_LIFE, IDAUGHTER, CONFACTOR, NEWCONFACTOR, CONFTIME, GMOL, P_FRACTION, ASTEP, CFRACTION</li>
</ul>
<p>Optional keyword “size” is input to indicate that the colloid size distribution model option is being used. It is followed by PART_SIZE and PROBSIZE.</p>
<blockquote>
<div><p>KEYWORD ‘size’</p>
<p>PART_SIZE(I), PROBSIZE(I) - an arbitrary numbers of lines of input, terminated by a blank line.</p>
</div></blockquote>
<p>Optional keyword “dive” is input to indicate that the colloid diversity model is being used. It is followed by FLAG_COL_DAUGHTER, optional keyword “file” and the name of the file containing the CDF table or equation data (a description of the format for the file is provided with the second mptr example below), the TPRPFLAG with optional SIMNUM, optional CDF equation parameters (when “file” is not used), and keyword “irreversible” or “reversible” with FLAG_LOG.</p>
<table border="1" class="docutils">
<colgroup>
<col width="41%" />
<col width="5%" />
<col width="54%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>KEYWORD ‘dive’</td>
<td>or</td>
<td>KEYWORD ‘dive</td>
</tr>
<tr class="row-even"><td>FLAG_COL_DAUGHTER</td>
<td>&#160;</td>
<td>FLAG_COL_DAUGHTER</td>
</tr>
<tr class="row-odd"><td>KEYWORD ‘file’</td>
<td>&#160;</td>
<td>TPRPFLAG</td>
</tr>
<tr class="row-even"><td>CDFFILENAME</td>
<td>&#160;</td>
<td>or</td>
</tr>
<tr class="row-odd"><td>TPRPFLAG</td>
<td>&#160;</td>
<td>TPRPFLAG, SIMNUM</td>
</tr>
<tr class="row-even"><td>or</td>
<td>&#160;</td>
<td>K_REV, R_MIN,R_MAX, SLOPE_KF, CINT_KF</td>
</tr>
<tr class="row-odd"><td>TPRPFLAG, SIMNUM</td>
<td>&#160;</td>
<td>KEYWORD ‘irreversible’</td>
</tr>
<tr class="row-even"><td>KEYWORD ‘irreversible’</td>
<td>&#160;</td>
<td>or</td>
</tr>
<tr class="row-odd"><td>or</td>
<td>&#160;</td>
<td>KEYWORD ‘reversible’, FLAG_LOG</td>
</tr>
<tr class="row-even"><td>KEYWORD ‘reversible’, FLAG_LOG</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>Note that optional KEYWORDs “size” and “dive” are only used when colloid transport is enabled.</p>
<ul class="simple">
<li>Group 8 - LAYERS</li>
<li>Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, DIFFMFL</li>
</ul>
<p>or for simulations using “dfree”:</p>
<ul class="simple">
<li>Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, H2O_DIFF, TORT_DIFF</li>
</ul>
<p>or for simulations with colloid (<code class="docutils literal notranslate"><span class="pre">TRANSFLAG</span> <span class="pre">&lt;</span> <span class="pre">0</span></code>):</p>
<ul class="simple">
<li>Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, DIFFMFL, KCOLL, RCOLL, FCOLL</li>
</ul>
<p>or for simulations with colloid using “dfree”:</p>
<ul class="simple">
<li>Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, H2O_DIFF, TORT_DIFF, KCOLL, RCOLL, FCOLL</li>
<li>Group 10 - NS</li>
<li>Group 11 - JA, JB, JC, TMPCNSK</li>
</ul>
<p>Note that because the number of source terms is controlled by the value entered for NS, Group 11 input is not terminated with a blank line.</p>
<ul class="simple">
<li>Group 12 - PINMASS, T1SK, T2SK</li>
</ul>
<p>For transient source terms, Group 12 is repeated for each time interval and terminated with a blank line. Groups 11 and 12 are repeated for each source term (from 1 to NS).</p>
<p>For decay-ingrowth calculations, when the particle injection period is too small (for example, 1.E-4 days) compared to the half-life of the radionuclides and the half-life is large (for example 1.E+9 days), numerical errors in the decay-ingrowth calculation may arise due to truncation error. To get better accuracy, the user should try to increase the length of the injection period.</p>
<p>For particle tracking simulations using the transfer function method (see <a class="reference external" href="Macro49660.html">See Transfer function curve data input file</a> for input file format), it is sometimes desirable to identify the parameter ranges over which the two- and three-parameter type curves are accessed, so that an assessment can be made regarding the density of transfer function curves in a given part of the parameter space. If the flag output_flag in the transfer function file is set to “out”, the code writes the real*8 array param_density to the <a href="#id1"><span class="problematic" id="id2">*</span></a>.out file in the following format:</p>
<p>For regular parameter spacings, the output is:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">i</span> <span class="o">=</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nump1</span>
   <span class="n">j</span> <span class="o">=</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nump2</span>
       <span class="n">k</span> <span class="o">=</span> <span class="n">nump3</span>

           <span class="n">write</span><span class="p">(</span><span class="n">iout</span><span class="o">.*</span><span class="p">)</span> <span class="n">param_density</span><span class="p">(</span><span class="n">i</span><span class="p">,</span><span class="n">j</span><span class="p">,</span><span class="n">k</span><span class="p">)</span>

       <span class="n">end</span> <span class="n">do</span>
   <span class="n">end</span> <span class="n">do</span>
<span class="n">end</span> <span class="n">do</span>
</pre></div>
</div>
<p>For two-parameter models, only the i and j loops are used. The value of param_density is the number of times any particle passes through any node at those values of the parameters. This allows the user to identify regions in which a greater density of transfer functions may be required. For the option ‘free’ in which there is no structure to the parameter grid used for the transfer function curves, nump1 is the total number of curves, and nump2 and nump3 are equal to 1.</p>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="11%" />
<col width="73%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NSPECI</td>
<td>integer</td>
<td>Number of species in the simulation.</td>
</tr>
<tr class="row-odd"><td>MAXLAYERS</td>
<td>integer</td>
<td>Maximum number of property layers in the model.
The actual number of layers used in the model must be ≤ MAXLAYERS.</td>
</tr>
<tr class="row-even"><td>MAX_PARTICLES</td>
<td>integer</td>
<td>Maximum number of particles used for individual species.</td>
</tr>
<tr class="row-odd"><td>RIPFEHM</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to indicate if simulation is coupled with GoldSim.</div>
<div class="line-block">
<div class="line">RIPFEHM = 0, FEHM standalone simulation</div>
<div class="line">RIPFEHM = 1, GoldSim-FEHM coupling simulation</div>
</div>
</div>
</td>
</tr>
<tr class="row-even"><td>MAX1D</td>
<td>integer</td>
<td>Maximum 1-D array size for holding particle tracking information for
all simulated species. The value of MAX1D depends on number of species,
number of time steps, number of radionuclide release bins, number of
species involved in ingrowth, and the length of the decay-ingrowth chain.</td>
</tr>
<tr class="row-odd"><td>POUT</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to specify the concentration output format:</div>
<div class="line">1 -  Concentrations computed as number of particles per unit total volume
(rock and fluid)</div>
<div class="line">2 -  Concentrations computed as number of particles per unit fluid volume
(the fluid is liquid for TRAK_TYPE = 1 and gas for TRAK_TYPE = 2).</div>
<div class="line">3 -  Concentrations computed as number of particles at a given node point.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>PRNT_RST</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to specify whether particle information is written to the “.fin”,
“.ptrk_fin”, or “.ptrk” files:</div>
<div class="line">If PRNT_RST = 0, Particle information is not written to the output files.</div>
<div class="line">If PRNT_RST = 1, 11, 21, 31, 41 All particle information necessary for a
restart is written to the “.fin” file.</div>
<div class="line">If PRNT_RST = -1, -11, -21, -31, -41 Only particle positions and ages are
written to the “.fin” file.</div>
<div class="line">If ABS (PRNT_RST) = 2, 12, 22, 32, 42 Mass flux values are written to the
“.fin” file followed by particle information.</div>
<div class="line">If 10 ≤ ABS(PRNT_RST) &lt; 30 Particle exit locations and count are written
to the “.ptrk_fin” file.</div>
<div class="line">If ABS(PRNT_RST) ≥ 20 Cumulative particle counts versus time are written
to the “.ptrk” file, for variables specified by PRNT_VARNUM (the default
is to output all variables).</div>
<div class="line">If ABS(PRNT_RST) ≥ 40, Cumulative mass output from a FEHM/GoldSim coupled
simulation will be written to file <code class="docutils literal notranslate"><span class="pre">FEHM_GSM_Mass_balance.txt</span></code>. Note that to
track cumulative mass an additional array of size <code class="docutils literal notranslate"><span class="pre">maxparticles*nspeci</span></code> must
be allocated so caution should be used when specifying this option to ensure
sufficient system memory is available.</div>
<div class="line"><br /></div>
<div class="line">When particle tracking data or mass fluxes are written to the <code class="docutils literal notranslate"><span class="pre">.fin</span></code> file,
the arrays are written after all of the heat and mass simulation information.
The mass fluxes can be read into the code in a subsequent ptrk or mptr simulation
and the code can simulate transport on this steady state flow field (see macro
<code class="docutils literal notranslate"><span class="pre">rflo</span></code>).The particle information written is sufficient to perform a restart of the
particle tracking simulation and to post-process the data to compile statistics
on the particle tracking run. However, for a large number of particles, this
file can become quite large, so particle tracking information should only be
written when necessary. Thus, 0 should be used for <code class="docutils literal notranslate"><span class="pre">PRNT_RST</span></code> unless restarting
or post-processing to obtain particle statistics is required. Selecting the
“-” options allows a subset of the full set of information needed for a
restart (particle positions and ages) to be written. Restart runs that use
this file as input will only be approximate, since the particle is assumed
to have just entered its current cell. For restart runs, <code class="docutils literal notranslate"><span class="pre">PRNT_RST</span> <span class="pre">=</span> <span class="pre">1</span></code> is
preferred, while <code class="docutils literal notranslate"><span class="pre">PRNT_RST</span> <span class="pre">=</span> <span class="pre">-1</span></code> is appropriate for output of particle
statistics for post- processing.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>PRNT_VARNUM</td>
<td>integer</td>
<td><div class="first line-block">
<div class="line">A list of integers specifying which particle counts should be output. For each
value entered <code class="docutils literal notranslate"><span class="pre">PRNT_VAR(PRNT_VARNUM)</span></code> is set to true. If no values are entered
the default is to print all variables.</div>
<div class="line">1 – Number of particles that have entered the system</div>
<div class="line">2 – Number of particles currently in the system</div>
<div class="line">3 – Number of particles that have left the system</div>
<div class="line">4 – Number of particles that have decayed</div>
<div class="line">5 – Number of particles that have been filtered</div>
<div class="line">6 – Number of particles that left this time interval</div>
</div>
<div class="last line-block">
<div class="line">Note: The data found in the “.ptrk” file was previously reported in the
general output file. From version 2.25 of the code and forward that data
will be reported in the optional, “.ptrk” file unless a coupled GoldSim-FEHM
simulation is being run. In addition, the user has the option of selecting
which statistics parameters are reported. The default is to report all
statistics parameters.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “tcurve” indicating transfer function curve data should be
input to model matrix diffusion. If the keyword is found then NUMPARAMS and
FILENAME are entered, otherwise they are omitted.</td>
</tr>
<tr class="row-odd"><td>NUMPARAMS</td>
<td>integer</td>
<td>Number of parameters that define the transfer function curves being used.</td>
</tr>
<tr class="row-even"><td>FFMAX</td>
<td>real</td>
<td>The maximum fracture flow fraction used in the transfer function curve
data. Default value: 0.99.</td>
</tr>
<tr class="row-odd"><td>TFILENAME</td>
<td>character</td>
<td>Name of input file containing the transfer function curve data.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘zptr’ designating zones for breakthrough curves will be
defined. If no keyword is input, IPZONE and IDZONE are also omitted.</td>
</tr>
<tr class="row-odd"><td>IPZONE</td>
<td>integer</td>
<td>Number of zones for which breakthrough curves are to be output</td>
</tr>
<tr class="row-even"><td>IDZONE</td>
<td>integer</td>
<td>A list of zones for which particle breakthrough data are required. The code
outputs the number of particles that leave the system at each zone IDZONE
at the current time step. This information is written to the “.out” file
at each heat and mass transfer time step.</td>
</tr>
<tr class="row-odd"><td>RSEED</td>
<td>integer</td>
<td>6-digit integer random number seed.</td>
</tr>
<tr class="row-even"><td>RSEED_RELEASE</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">6-digit integer random number seed for particle release location calculation.
If a value is not entered for RSEED_RELEASE it will be set equal to RSEED.</div>
<div class="line"><br /></div>
<div class="line">Note that for GoldSim-FEHM coupled simulations the random seeds are controlled
by GoldSIM and the values input in the mptr macro are not used.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘wtri” indicatiing a water table rise calculation should
be performed.</td>
</tr>
<tr class="row-even"><td>WATER_TABLE</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Water table elevation to be used for water table rise calculation.</div>
<div class="line">Note that for GoldSim-FEHM coupled simulations the water table
rise calculations are controlled by GoldSIM and the values input
in the mptr macro are not used and may be omitted.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>DAYCS</td>
<td>real</td>
<td>Time which the particle tracking solution is enabled (days).</td>
</tr>
<tr class="row-even"><td>DAYCF</td>
<td>real</td>
<td>Time which the particle tracking solution is disabled (days).</td>
</tr>
<tr class="row-odd"><td>DAYHF</td>
<td>real</td>
<td>Time which the flow solution is disabled (days).</td>
</tr>
<tr class="row-even"><td>DAYHS</td>
<td>real</td>
<td>Time which the flow solution is enabled (days).</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating alternate transport
parameter file input for multiple simulation realizations.</td>
</tr>
<tr class="row-even"><td>PFILENAME</td>
<td>character*80</td>
<td>Name of file from which to read transport parameters.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘afm’ designating the Active Fracture
Model input for determining fracture spacing should be used.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*5</td>
<td>Optional keyword ‘dfree’ designates that the free water
diffusion coefficient and tortuosity will be input instead
of the molecular diffusion coefficient.</td>
</tr>
<tr class="row-odd"><td>TCLX</td>
<td>real</td>
<td>Dispersivity in the x-direction (m). The input value is
ignored when dispersion is turned off.</td>
</tr>
<tr class="row-even"><td>TCLY</td>
<td>real</td>
<td>Dispersivity in the y-direction (m). The input value is
ignored when dispersion is turned off.</td>
</tr>
<tr class="row-odd"><td>TCLZ</td>
<td>real</td>
<td>Dispersivity in the z-direction (m). The input value is
ignored when dispersion is turned off.</td>
</tr>
<tr class="row-even"><td>APERTUR</td>
<td>real</td>
<td>Mean fracture aperture (m). The input value is ignored
when matrix diffusion is turned off.</td>
</tr>
<tr class="row-odd"><td>MATRIX_POR</td>
<td>real</td>
<td>Porosity of the rock matrix. Used to simulate diffusion
and sorption in the rock matrix when matrix diffusion
is invoked, otherwise the input value of MATRIX_POR is ignored.</td>
</tr>
<tr class="row-even"><td>SRESIDUAL</td>
<td>real</td>
<td>Residual saturation in the Active Fracture Model used for
determining the spacing between active fractures.
This parameter is only needed when the keyword ‘afm’ is
included, in which case the input must be entered.
However, the model is only used in dual permeability
simulations at locations where the finite spacing matrix
diffusion model is invoked.</td>
</tr>
<tr class="row-odd"><td>GAMMA_AFM</td>
<td>real</td>
<td>Exponent in the Active Fracture Model used for determining
the spacing between active fractures. See comments for
SRESIDUAL above.</td>
</tr>
<tr class="row-even"><td>ITRC</td>
<td>integer</td>
<td>Model number for parameters defined in group 5.</td>
</tr>
<tr class="row-odd"><td>ITH_SPECI</td>
<td>integer</td>
<td>Number index of the ith species.</td>
</tr>
<tr class="row-even"><td>TRAK_TYPE</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to denote the fluid phase of the particles:</div>
<div class="line">1 - liquid phase particles</div>
<div class="line">2 - vapor phase particles</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>HALF_LIFE</td>
<td>real</td>
<td>Half-life for irreversible first order decay reaction(s)
(days). Set HALF_LIFE = 0 for no decay.</td>
</tr>
<tr class="row-even"><td>IDAUGHTER</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Index of the daughter species (i.e., the index number of
the species to which the current species decays)</div>
<div class="line">If IDAUGHTER = 0, there is no decay and no ingrowth,</div>
<div class="line">If IDAUGHTER = -1, there is decay but no ingrowth.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>CONFACTOR</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Initial conversion factor for GoldSim-FEHM coupling and
FEHM standalone simulations (# of particles/mole).</div>
<div class="line"><br /></div>
<div class="line">For FEHM stand alone simulations:</div>
<div class="line">If CONFACTOR = 0, no conversion is necessary. The input
value of PINMASS is the number of particles.</div>
<div class="line"><br /></div>
<div class="line">For GoldSim-FEHM coupling:</div>
<div class="line">If CONFACTOR = 0, at each time step, the code selects a
conversion factor based on the available memory and the
remaining simulation time (end time - current time). The
code then uses the selected conversion factor to calculate
the number of particles to be injected at the current time step.</div>
<div class="line"><br /></div>
<div class="line">For both stand alone and GoldSim-FEHM coupling cases:</div>
<div class="line">If CONFACTOR &gt; 0, the code assumes the input mass is in moles
and uses the product of the CONFACTOR and the input mass to calculate
the input number of particles at each time step.</div>
<div class="line">When CONFACTOR &gt;0, FEHM may use an updated conversion factor from
previous time step(s) as the input for the current time step instead
of using the original input CONFACTOR for improved results.</div>
<div class="line"><br /></div>
<div class="line">If CONFACTOR &lt; 0, the code uses the product of the absolute value
of CONFACTOR and the input mass (in moles) to calculate the input
number of particles at each time step. A CONFACTOR updated from a
previous time step will not be used.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>NEWCONFACTOR</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Replace the initial value of CONFACTOR with that specified by
NEWCONFACTOR.</div>
<div class="line">If NEWCONFACTOR = 0, use automatic conversion factors.</div>
<div class="line">If NEWCONFACTOR &gt; 0, then use the product of the CONFACTOR
and the input mass (in moles) to calculate the input number of
particles at each time step starting from CONFTIME. In this case,
FEHM may use an updated conversion factor from previous time step(s)
as a modification to CONFACTOR.</div>
<div class="line">If NEWCONFACTOR &lt; 0, then FEHM uses the product of the absolute value
of NEWCONFACTOR and the input mass (in moles) to calculate the input
number of particles at each time step (<code class="docutils literal notranslate"><span class="pre">CONFACTOR</span> <span class="pre">=</span> <span class="pre">-NEWCONFACTOR</span></code>).</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>CONFTIME</td>
<td>real</td>
<td>The time at which to change the CONFACTOR value to that specified
by NEWCONFACTOR.</td>
</tr>
<tr class="row-even"><td>GMOL</td>
<td>real</td>
<td>The molecular weight of the ith species. The code uses GMOL and
CONFACTOR to convert the mass from number of particles to grams
in the final output for GoldSim-FEHM coupling.</td>
</tr>
<tr class="row-odd"><td>P_FRACTION</td>
<td>real</td>
<td>The decay-ingrowth particle release factor (percentage of the
maximum number of particles released for the current species).
For decay-ingrowth simulations, P_FRACTION is used to reduce
the number of particles released by parent or daughter species,
thus, avoiding memory overflow in the daughter species due to
parent decay-ingrowth where multiple parents decay to the same
daughter species. The normal range of P_FRACTION is from 0 to 1.
The default value is 0.25. A user should select an appropriate
value based on the mass input of parent and daughter species,
half-lives, importance of each species to the transport results,
and simulation time period.</td>
</tr>
<tr class="row-even"><td>ASTEP</td>
<td>integer</td>
<td>Maximum length of array used to hold particle tracking information
for the ith species. Its value depends on number of time steps,
number of release bins, and number of parent species.
The sum of ASTEP for all species should be equal to or smaller
than MAX1D.</td>
</tr>
<tr class="row-odd"><td>CFRACTION</td>
<td>real</td>
<td>The fraction of the user determined maximum number of particles
(MAX_PARTICLES) to be assigned by mass, (1 – cfraction) will
then be the fraction of particles assigned by time step.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘size’ designating that the colloid size
distribution model option is being used (combined with
the interface filtration option in the itfc macro). If the
keyword is not input, PART_SIZE and PROBSIZE are also omitted.
Colloid size is only sampled once for each realization.</td>
</tr>
<tr class="row-odd"><td>PART_SIZE</td>
<td>real</td>
<td>Colloid particle size for this entry of the particle size
distribution table (paired with a value of PROBSIZE).
An arbitrary number of entries can be input, terminated with
a blank line. The code assigns each particle a size based on
this distribution of particle sizes, and decides if particles
are irreversibly filtered based on the pore size distribution
assigned in the itfc macro.</td>
</tr>
<tr class="row-even"><td>PROBSIZE</td>
<td>real</td>
<td>Colloid cumulative probability for the distribution of sizes
(paired with a value of PART_SIZE). See description of
PART_SIZE above for details. The final entry of the table
must have PROBSIZE = 1, since the distribution is assumed
to be normalized to unity.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword “dive” signifying that the specie being
specified is either a colloid species using the colloid
diversity model or a non-colloid daughter species of a
colloid species.</td>
</tr>
<tr class="row-even"><td>FLAG_COL_DAUGHTER</td>
<td>integer</td>
<td>When FLAG_COL_DAUGHTER = 1 signals that the species being
specified is a non-colloid species that can result as a
daughter product of a colloid parent species. If the species
is not a daughter product or the daughter product is a
colloid, FLAG_COL_DAUGHTER = 0.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating the cumulative probability
distribution function (CDF) retardation parameters for the
colloid diversity model should be read from an external file.</td>
</tr>
<tr class="row-even"><td>CDF_FILENAME</td>
<td>character*80</td>
<td><div class="first last line-block">
<div class="line">Name of the file containing the cumulative probability
distribution function (CDF) (entered if optional keyword
‘file’ follows keyword ‘dive’). See below for file formats.</div>
<div class="line">If TPRPFLAG = 11 or 12, Table option</div>
<div class="line">If TPRPFLAG = 13 or 14, Equation option</div>
<div class="line"><br /></div>
<div class="line">The following equations are used for <span class="math notranslate nohighlight">\(R_{min} \le R \le R_{max}\)</span>,</div>
<div class="line"><span class="math notranslate nohighlight">\(R = 1 + K_f / K_{rev}\)</span>,
<span class="math notranslate nohighlight">\(\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>TPRP_FLAG</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Values of TPRPFLAG between 11 and 14 signify that the colloid
diversity model with equal weight sampling will be used:</div>
<div class="line">TPRPFLAG = 11: CDF vs retardation factor specified in a table</div>
<div class="line">TPRPFLAG = 12: similar to 11, but the SQRT(CDF) is used instead
of CDF for sampling</div>
<div class="line">TPRPFLAG = 13: CDF vs <span class="math notranslate nohighlight">\(K_f\)</span> (Attachment rate constant)
specified as a straight line equation in the log-log space</div>
<div class="line">TPRPFLAG = 14: similar to 13, but the SQRT(CDF) is used
instead of CDF for sampling</div>
</div>
</td>
</tr>
<tr class="row-even"><td>SIMNUM</td>
<td>integer</td>
<td>Simulation number, used for selecting the table/equation from
the colloid diversity file. For GoldSim-FEHM coupled simulations
or FEHM runs using the ‘msim’ option this parameter is passed
to the code. For non-coupled simulations it is an optional
input. (Default value = 1)</td>
</tr>
<tr class="row-odd"><td>K_REV</td>
<td>real</td>
<td>Detachment rate constant for reversible filtration of
irreversible colloids.</td>
</tr>
<tr class="row-even"><td>R_MIN</td>
<td>real</td>
<td>Minimum value of the retardation factor for reversible
filtration of irreversible colloids.</td>
</tr>
<tr class="row-odd"><td>R_MAX</td>
<td>real</td>
<td>Maximum value of the retardation factor for reversible
filtration of irreversible colloids</td>
</tr>
<tr class="row-even"><td>SLOPE_KF</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Value of the slope (<span class="math notranslate nohighlight">\(m\)</span>) in the log-log space
for the equation:</div>
<div class="line"><span class="math notranslate nohighlight">\(\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>CINT_KF</td>
<td>real</td>
<td>Value of the intercept (<span class="math notranslate nohighlight">\(b\)</span>) in the log-log space
for the above equation</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Keyword specifying whether the colloid species is
‘irreversible’ or ‘reversible’.</td>
</tr>
<tr class="row-odd"><td>FLAG_LOG</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">For reversible colloids an average retardation factor is used:</div>
<div class="line"><br /></div>
<div class="line">If FLAG_LOG = 0: a linear average of the distribution is used</div>
<div class="line">If FLAG_LOG = 1: a log-linear average of the distribution
is used</div>
</div>
</td>
</tr>
<tr class="row-even"><td>LAYERS</td>
<td>integer</td>
<td>Number of layers in which the transport properties of the
ith species are to be modified. If no property is altered,
then set layers=0.</td>
</tr>
<tr class="row-odd"><td>LAYER_I</td>
<td>integer</td>
<td>The index number of the ith layer defined in group 5</td>
</tr>
<tr class="row-even"><td>TRANSFLAG</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to specify which transport mechanisms apply [abs(TRANSFLAG)]:</div>
<div class="line">1 - advection only (no dispersion or matrix diffusion)</div>
<div class="line">2 - advection and dispersion (no matrix diffusion)</div>
<div class="line">3 - advection and matrix diffusion, infinite fracture spacing
solution (no dispersion)</div>
<div class="line">4 - advection, dispersion, and matrix diffusion, infinite fracture
spacing solution</div>
<div class="line">5 - advection and matrix diffusion, finite fracture spacing
solution (no dispersion)</div>
<div class="line">6 - advection, dispersion, and matrix diffusion, finite fracture
spacing solution</div>
<div class="line">8 - use the the transfer function approach with 3 dimensionless
parameters and type curves for handling fracture-matrix interactions.</div>
<div class="line"><br /></div>
<div class="line">For TRANSFLAG &lt; 0, transport simulations include colloids.</div>
<div class="line">For equivalent continuum solutions, the fracture spacing in the
finite spacing model is determined using</div>
<div class="line"><br /></div>
<div class="line"><span class="math notranslate nohighlight">\(SPACING = APERTURE / POROSITY\)</span></div>
<div class="line"><br /></div>
<div class="line">For dual permeability models, the fracture spacing input parameter
APUV1 in the <code class="docutils literal notranslate"><span class="pre">dpdp</span></code> macro is used as the half-spacing between fractures.
If the Active Fracture Model (see keyword ‘afm’) is used, APUV1 is the
geometric fracture half-spacing, and the additional terms SRESIDUAL
and GAMMA_AFM are used to determine the spacing between active
fractures (see below).</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>KD</td>
<td>real</td>
<td>Sorption coefficient (linear, reversible, equilibrium sorption).
Units are kg-fluid / kg-rock (these units are equivalent to the
conventional units of cc/g when the carrier fluid is water at
standard conditions). This value applies to the medium as a whole
when matrix diffusion is turned off, whereas for simulations invoking
matrix diffusion, the value applies to the rock matrix.
For the latter case, sorption in the flowing system (fractures)
is modeled using the RD_FRAC variable.</td>
</tr>
<tr class="row-even"><td>RD_FRAC</td>
<td>real</td>
<td>Retardation factor within the primary porosity (fractures) for a
matrix diffusion particle tracking simulation (use 1 for no
sorption on fracture faces). The input value is ignored unless
matrix diffusion is invoked.</td>
</tr>
<tr class="row-odd"><td>DIFFMFL</td>
<td>real</td>
<td>Molecular diffusion coefficient in the rock matrix (m2/s).
The input value is ignored unless matrix diffusion is invoked.</td>
</tr>
<tr class="row-even"><td>H2O_DIFF</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Free water diffusion coefficient. The molecular diffusion
coefficient is calculated as
<span class="math notranslate nohighlight">\(H2O\_DIFF \times TORT\_DIFF\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>TORT_DIFF</td>
<td>real</td>
<td>Tortuosity</td>
</tr>
<tr class="row-even"><td>KCOLL</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Colloid distribution parameter, the ratio of contaminant mass
residing on colloids to the mass present in aqueous form.
It is used to compute an effective aperture via the following:</div>
<div class="line"><span class="math notranslate nohighlight">\(APWID = APERERTURE \cdot (1 + KCOLL)\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>RCOLL</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Colloid retardation factor. Used, in conjunction with kcoll,
to adjust colloid retardation in fractures using the following
formula:</div>
<div class="line"><span class="math notranslate nohighlight">\(FRACRD = \frac{RD\_FRAC + KCOLL \cdot RCOLL}{1+KCOLL}\)</span></div>
</div>
</td>
</tr>
<tr class="row-even"><td>FCOLL</td>
<td>real</td>
<td>Colloid filtration parameter. Used to compute the probability a colloid
will be irreversibly filtered along the path between two nodes using
the following:
<span class="math notranslate nohighlight">\(PROBFILT = 1 - \exp(DISTANCE/FCOLL)\)</span> where
<span class="math notranslate nohighlight">\(DISTANCE\)</span> iis the leength of the path between nodes.</td>
</tr>
<tr class="row-odd"><td>NS</td>
<td>&#160;</td>
<td>Number of spatial source terms for the ith species</td>
</tr>
<tr class="row-even"><td>TMPCNSK</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Particle injection parameter assigned for nodes defined by JA,
JB, and JC. Two options are available:</div>
<div class="line"><br /></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td><div class="first last line-block">
<div class="line">TMPCNSK &gt; 0. - particles are injected at each node in
proportion to the source mass flow rate at the node. This
boundary condition is equivalent to injecting a solute of
a given concentration into the system. Note: the source
flow rates used to assign the number and timing of particle
injections are those at the beginning of the particle
tracking simulation (time DAYCS). Transient changes in this
source flow rate during the particle tracking simulation do
not change the number of particles input to the system.</div>
<div class="line">TMPCNSK &lt; 0. - particles are introduced at the node(s),
regardless of whether there is a fluid source at the node.</div>
<div class="line"><br /></div>
<div class="line">Default is 0. for all unassigned nodes, meaning that no
particles are injected at that node.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>PINMASS</td>
<td>real</td>
<td>Input mass. If CONFACTOR = 0, PINMASS is the number of particles
to be injected at locations defined by TMPCNSK.
If CONFACTOR &gt; 0, PINMASS is the input mass expressed in moles.
The code uses CONFACTOR to convert PINMASS into number of particles.</td>
</tr>
<tr class="row-odd"><td>T1SK</td>
<td>real</td>
<td>Time (days) when particle injection begins. Default is 0.</td>
</tr>
<tr class="row-even"><td>T2SK</td>
<td>real</td>
<td>Time (days) when particle injection ends. Default is 0.</td>
</tr>
</tbody>
</table>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last"><strong>Notes on Restarting:</strong> As with all restart runs for FEHM, a “.ini” file is
specified to be read to set the initial conditions upon restarting. However,
there are two possibilities for restart calculations with particle tracking
(mptr or ptrk): 1) the heat and mass transfer solution is being restarted, but
the particle tracking simulation is initiated during the restart run (it was not
carried out in the simulation that generated the “.ini” file); or 2) the heat
and mass transfer solution and the particle tracking simulation are both being
restarted. If the code does not find the “ptrk” key word at the top of the “.ini”
file, then the original run did not employ particle tracking, and Case 1 is assumed.
A common example is a preliminary calculation that establishes a fluid flow steady
state, followed by a restart simulation of transport.</p>
</div>
<p>If “ptrk” was written into the “.ini” file in the original run, the particle data in the “.ini” file are read and used to initialize the particle tracking simulation (Case 2). In this instance, the number of particles (NPART) must be set the same for the restart run as in the original run or the results will be unpredictable.When restarting a particle tracking simulation, certain input data are overwritten by information in the “.ini” file. These parameters include RSEED, RSEED_RELEASE, PCNSK, T1SK, and T2SK. Other input parameters can be set to different values in the restart run than they were in the original run, but of course care must be taken to avoid physically unrealistic assumptions, such as an abrupt change in transport properties (Group 4 input) part way through a simulation.</p>
<p>A final note on restart calculations is in order. A common technique in FEHM restart calculations is to reset the starting time at the top of the “.ini” file to 0 or in the time macro so that the starting time of the restart simulation is arbitrarily 0, rather than the ending time of the original simulation. This is useful for the example of the steady state flow calculation, followed by a restart solute transport calculation. Although this technique is acceptable for particle tracking runs that are initiated only upon restart (Case 1), it is invalid when a particle tracking run is being resumed (Case 2). The reason is that all particle times read from the “.ini” file are based on the starting time of the original simulation during which the particle tracking simulation was initiated.</p>
<p>The following is an example of mptr. A multiple-species decay-chain
(<span class="math notranslate nohighlight">\(\rightarrow 2 \rightarrow 3\)</span>) is simulated, with decay half lives of the
species equaling 10,000, 3,000, 10,000, and 4,000 years, respectively. In this
simulation a maximum of 3 property layers are specified although only 1 layer is used,
the maximum number of particles is specified to be 1100100, and FEHM is run in stand-alone
mode. Concentrations will be computed as number of particles per unit fluid volume and
no output will be written to the “.fin” file. Use of the ‘zptr’ keyword indicates that
a single zone will be defined for breakthrough curve output which will be written to
the “.out” file. The random number seed is defined to be 244562. The particle tracking
solution is enabled at 0.1 days, and disabled at 3.65e8 days, while the flow solution
is disabled at 38 days and re-enabled at 3.65e8 days. Dispersivity in the X-, Y-, and
Z- directions are defined to be 0.005 m, the mean fracture aperture is 0.0001 m, and
the matrix porosity is 0.3. Particles for species 1 are injected at a constant rate
from 0 to 5,000 years, and species 2, 3, and 4 are formed through the decay reactions,
with no input at the inlet. Advection and dispersion (without matrix diffusion) is being
modeled. The retardation factors for the four species are 1, 1, 1.9, and 1, respectively
(i.e. only species 3 sorbs).</p>
<table border="1" class="docutils">
<colgroup>
<col width="10%" />
<col width="15%" />
<col width="15%" />
<col width="14%" />
<col width="10%" />
<col width="5%" />
<col width="5%" />
<col width="5%" />
<col width="6%" />
<col width="13%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>mptr</td>
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
<tr class="row-even"><td>4</td>
<td>3</td>
<td>1100100</td>
<td>0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 1</td>
</tr>
<tr class="row-odd"><td>2</td>
<td>0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 2</td>
</tr>
<tr class="row-even"><td>zptr</td>
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
<tr class="row-odd"><td>1</td>
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
<tr class="row-even"><td>1</td>
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
<tr class="row-odd"><td>244562</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 3</td>
</tr>
<tr class="row-even"><td>0.1</td>
<td>3.65e8</td>
<td>38</td>
<td>3.65e8</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 4</td>
</tr>
<tr class="row-odd"><td>0.005</td>
<td>0.005</td>
<td>0.005</td>
<td>1.e-4</td>
<td>0.3</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 5</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
<tr class="row-odd"><td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 6</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
<tr class="row-odd"><td>1</td>
<td>1</td>
<td>3.652485E6</td>
<td>2</td>
<td>1</td>
<td>-1</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>0.5</td>
<td>Group 7</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 8</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>2</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>1.e-14</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 9</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 10</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>202</td>
<td>201</td>
<td>-1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 11</td>
</tr>
<tr class="row-even"><td><ol class="first last arabic simple" start="10000">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>&#160;</td>
<td>365.25E2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="10000">
<li></li>
</ol>
</td>
<td>365.25E2</td>
<td>&#160;</td>
<td>730.5E2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-even"><td><ol class="first last arabic simple" start="10000">
<li></li>
</ol>
</td>
<td>730.5E2</td>
<td>&#160;</td>
<td>1.09575E5</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="10000">
<li></li>
</ol>
</td>
<td>1.09575E5</td>
<td>&#160;</td>
<td>1.461E5</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-even"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>.</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>.</td>
</tr>
<tr class="row-even"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>.</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="10000">
<li></li>
</ol>
</td>
<td>1.7532E6</td>
<td>1.789725E6</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-even"><td><ol class="first last arabic simple" start="10000">
<li></li>
</ol>
</td>
<td>1.789725E6</td>
<td>1.82625E6</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
<tr class="row-even"><td>2</td>
<td>1</td>
<td>1.095745E6</td>
<td>3</td>
<td>0</td>
<td>-1</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>0.5</td>
<td>Group 7</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 8</td>
</tr>
<tr class="row-even"><td>1</td>
<td>2</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>1.e-14</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 9</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 10</td>
</tr>
<tr class="row-even"><td>1</td>
<td>202</td>
<td>201</td>
<td>-1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 11</td>
</tr>
<tr class="row-odd"><td>0</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>1.825E6</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
<tr class="row-odd"><td>3</td>
<td>1</td>
<td>3.652485E6</td>
<td>4</td>
<td>0</td>
<td>-1</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>0.5</td>
<td>Group 7</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 8</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>2</td>
<td>0.108</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>1.e-14</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 9</td>
</tr>
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 10</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>202</td>
<td>201</td>
<td>-1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 11</td>
</tr>
<tr class="row-even"><td>0</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>1.825E6</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
<tr class="row-even"><td>4</td>
<td>1</td>
<td>1.460972E6</td>
<td>-1</td>
<td>0</td>
<td>-1</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>0.5</td>
<td>Group 7</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 8</td>
</tr>
<tr class="row-even"><td>1</td>
<td>2</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>1.e-14</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 9</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 10</td>
</tr>
<tr class="row-even"><td>1</td>
<td>202</td>
<td>201</td>
<td>-1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 11</td>
</tr>
<tr class="row-odd"><td>0</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>1.825E6</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
</tbody>
</table>
<p>In the second example, transfer function data is used along with the active fracture and colloid diversity models. The the cumulative probability distribution function (CDF) retardation parameters for the colloid diversity model are entered in an external file using the table format. The format for the new input files associated with the colloid diversity model are:</p>
<p>For ptrk/sptr/mptr simulations with TPRP_FLAG = 11 or 12:</p>
<ul class="simple">
<li>Header line indicating the species number (always 1 for ptrk/sptr simulations)</li>
<li>Multiple tables, each with the following format:<ul>
<li>One line specifying the realization number</li>
<li>Multiple lines with two columns of data (real) representing rcdiv and probdiv. Note that probdiv should start at 0 and end with 1.</li>
<li>A blank line is used to specify the end of the distribution table.</li>
</ul>
</li>
<li>A blank line is used to specify the end of the file</li>
</ul>
<p>For mptr, a distribution will need to be entered for each colloid species. Therefore,
the header line and tables are repeated for each colloid species.</p>
<p>For ptrk/mptr simulations with TPRP_FLAG = 13 or 14:</p>
<ul class="simple">
<li>Header line containing comments</li>
<li>Multiple lines, each line containing realization number, <span class="math notranslate nohighlight">\(b, m, k_r, R_{min}, R_{max}\)</span></li>
<li>A blank line is used to specify the end of the file</li>
</ul>
<p>For sptr simulations with TPRP_FLAG = 13 or 14:</p>
<ul class="simple">
<li>Header line containing comments</li>
<li>Multiple lines, each line containing realization number, <span class="math notranslate nohighlight">\(b, m, k_r, R_{min}, R_{max}, \alpha_L\)</span></li>
<li>A blank line is used to specify the end of the file</li>
</ul>
<p>Particle statistics data for the cumulative number of particles that have left the sytem and the number of particles that left during the current timestep are output.</p>
<table border="1" class="docutils">
<colgroup>
<col width="31%" />
<col width="5%" />
<col width="5%" />
<col width="7%" />
<col width="7%" />
<col width="5%" />
<col width="9%" />
<col width="17%" />
<col width="7%" />
<col width="7%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>mptr</td>
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
<td>100</td>
<td>500000</td>
<td>0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 1</td>
</tr>
<tr class="row-odd"><td>0</td>
<td>30</td>
<td>3</td>
<td>6</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 2</td>
</tr>
<tr class="row-even"><td>tcurve</td>
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
<tr class="row-odd"><td>3</td>
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
<tr class="row-even"><td>../colloid_cell/input/uz_tfcurves_nn_3960.in</td>
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
<tr class="row-odd"><td>244562</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 3</td>
</tr>
<tr class="row-even"><td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>1.e20</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>1.e20</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 4</td>
</tr>
<tr class="row-odd"><td>afm</td>
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
<tr class="row-even"><td>1.e-03</td>
<td>1.e-03</td>
<td>1.e-03</td>
<td>1.e-03</td>
<td>0.2</td>
<td>0.01</td>
<td>0.6</td>
<td>// layer 1 tcwm1 zone 1</td>
<td>&#160;</td>
<td>Group 5</td>
</tr>
<tr class="row-odd"><td>1.e-03</td>
<td>1.e-03</td>
<td>1.e-03</td>
<td>0.00e+00</td>
<td>0.2</td>
<td>0.00</td>
<td>0.0</td>
<td>// layer 1 tcwm1 zone 2</td>
<td>&#160;</td>
<td>Group 5</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
<tr class="row-odd"><td>1</td>
<td>10</td>
<td>1</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 6</td>
</tr>
<tr class="row-even"><td>11</td>
<td>20</td>
<td>1</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 6</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
<tr class="row-even"><td>1</td>
<td>1</td>
<td>0</td>
<td>-1</td>
<td>1</td>
<td>0</td>
<td>1.00E+15243</td>
<td>1.0</td>
<td>speci1</td>
<td>Group 7</td>
</tr>
<tr class="row-odd"><td>diveresity</td>
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
<tr class="row-even"><td>0</td>
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
<tr class="row-odd"><td>file</td>
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
<tr class="row-even"><td>../colloid_cell/input/rcoll_data.dat</td>
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
<tr class="row-odd"><td>11 1</td>
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
<tr class="row-even"><td>irreversible</td>
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
<tr class="row-odd"><td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 8</td>
</tr>
<tr class="row-even"><td>1</td>
<td>-2</td>
<td>0.0</td>
<td>1.0e+00</td>
<td>1.00e-30</td>
<td>1e+20</td>
<td>1</td>
<td>1.00</td>
<td>#1 tcwM1</td>
<td>Group 9</td>
</tr>
<tr class="row-odd"><td>2</td>
<td>-2</td>
<td>0.0</td>
<td>1.0e+00</td>
<td>1.00e-30</td>
<td>1e+20</td>
<td>1</td>
<td>1.00</td>
<td>#2 tcwM2</td>
<td>Group 9</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
<tr class="row-odd"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 10</td>
</tr>
<tr class="row-even"><td>1</td>
<td>1</td>
<td>1</td>
<td>-1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 11</td>
</tr>
<tr class="row-odd"><td>100000</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>0.01</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-even"><td>&#160;</td>
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
<tr class="row-odd"><td>3</td>
<td>1</td>
<td>0</td>
<td>-1</td>
<td>1</td>
<td>0</td>
<td>1.00EE+15243</td>
<td>1.0</td>
<td>speci3</td>
<td>Group 7</td>
</tr>
<tr class="row-even"><td>diversity</td>
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
<tr class="row-odd"><td>0</td>
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
<tr class="row-even"><td>file</td>
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
<tr class="row-odd"><td>../colloid_cell/input/rcoll_data.dat</td>
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
<tr class="row-even"><td>11 3</td>
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
<tr class="row-odd"><td>irreversible</td>
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
<tr class="row-even"><td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 8</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>-2</td>
<td>0.0</td>
<td>1.0e+00</td>
<td>1.00e-30</td>
<td>1e+20</td>
<td>1</td>
<td>1.00</td>
<td>#1 tcwM1</td>
<td>Group 9</td>
</tr>
<tr class="row-even"><td>2</td>
<td>-2</td>
<td>0.0</td>
<td>1.0e+00</td>
<td>1.00e-30</td>
<td>1e+20</td>
<td>1</td>
<td>1.00</td>
<td>#2 tcwM2</td>
<td>Group 9</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
<tr class="row-even"><td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 10</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>1</td>
<td>1</td>
<td>-1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 11</td>
</tr>
<tr class="row-even"><td>0</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>0.01</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
</tbody>
</table>
<p>With the file <code class="docutils literal notranslate"><span class="pre">rcoll_data.dat</span></code> as:</p>
<table border="1" class="docutils">
<colgroup>
<col width="26%" />
<col width="29%" />
<col width="23%" />
<col width="23%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td colspan="4">Test file for 1D importance sampling</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0.0</td>
<td>0.0</td>
<td>0.0</td>
</tr>
<tr class="row-odd"><td>1.0</td>
<td>0.0</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>2.0</td>
<td>0.125</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>3.0</td>
<td>0.25</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>4.0</td>
<td>0.375</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>5.0</td>
<td>0.5</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>6.0</td>
<td>0.625</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>7.0</td>
<td>0.75</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>8.0</td>
<td>0.875</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>9.0</td>
<td>1.0</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>3</td>
<td>0.0</td>
<td>0.0</td>
<td>0.0</td>
</tr>
<tr class="row-even"><td>1.0</td>
<td>0.0</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>2.0</td>
<td>0.125</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>9.0</td>
<td>1.0</td>
<td>&#160;</td>
<td>&#160;</td>
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
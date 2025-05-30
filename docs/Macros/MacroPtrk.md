---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">ptrk</span></code><a class="headerlink" href="#ptrk" title="Permalink to this headline">¶</a></h1>
<p>Particle tracking simulation input. Note that data for each numbered group must be input. The other input is optional.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Macro <code class="docutils literal notranslate"><span class="pre">ptrk</span></code> <strong>cannot</strong> be used with <code class="docutils literal notranslate"><span class="pre">trac</span></code></p>
</div>
<ul class="simple">
<li>Group 1 - NPART, RSEED</li>
</ul>
<p>Optional keyword “wtri” is input to indicate a water table rise calculation should be performed. It is followed by WATER_TABLE. For GoldSim the water table rise calculation is controlled by passing a new water table elevation to the code during the simulation and the keyword is not required.</p>
<blockquote>
<div><p>KEYWORD ‘wtri’</p>
<p>WATER_TABLE</p>
<p>KEYWORD ‘rip’</p>
</div></blockquote>
<p>Group 2 - DAYCS, DAYCF, DAYHF, DAYHS</p>
<p>or if the optional keyword ‘rip’ follows Group 1</p>
<p>Group 2 - DAYCS, DAYCF, DAYHF, DAYHS, RIPFEHM, CONFREAD, GMOL, P_FRACTION</p>
<p>Group 3 - TRAK_TYPE, HALF_LIFE, POUT, PRNT_RST</p>
<blockquote>
<div>or when PRNT_RST ≥ 20, selected output parameters</div></blockquote>
<p>Group 3 - TRAK_TYPE, HALF_LIFE, POUT, PRNT_RST PRNT_VARNUM (V1 … V6)</p>
<p>Optional keyword “tcurve” is input to indicate that transfer function curves should be input to model matrix diffusion. It is followed by NUMPARAMS and TFILENAME.</p>
<blockquote>
<div><p>KEYWORD</p>
<p>NUMPARAMS, FFMAX</p>
<p>TFILENAME</p>
</div></blockquote>
<p>Optional keyword “zptr” designates zones for breakthrough curves will be defined. It is followed by IPZONE and IDZONE.</p>
<blockquote>
<div><p>KEYWORD ‘zptr’</p>
<p>IPZONE</p>
<p>IDZONE(I) I = 1 to IPZONE</p>
</div></blockquote>
<p>Group 4 is used to define models in which identical sorption and transport parameters are assumed to apply. Group 4 data are read until a blank line is encountered. The model number ITRC is incremented by 1 each time a line is read. Model parameters defined in Group 4 are assigned to nodes or zones using Group 5.</p>
<p>An optional, flexible input structure involving the assignment of transport parameters is implemented in the particle tracking input to allow multiple realization simulations to use different parameters for each realization. The user invokes this option using the keyword ‘file’ before Group 4, followed by the name of the file that the transport parameters reside in.</p>
<blockquote>
<div><p>KEYWORD ‘file’</p>
<p>PFILENAME</p>
</div></blockquote>
<p>The structure of the alternate parameter file is:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">NINPUTS</span>
<span class="n">PTRPARAM</span><span class="p">(</span><span class="n">I</span><span class="p">)</span> <span class="n">I</span><span class="o">=</span><span class="mi">1</span> <span class="n">to</span> <span class="n">NINPUTS</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="p">[</span><span class="n">a</span> <span class="n">line</span> <span class="n">of</span> <span class="n">parameters</span> <span class="ow">is</span> <span class="n">present</span> <span class="k">for</span> <span class="n">each</span> <span class="n">realization</span><span class="p">]</span>
</pre></div>
</div>
<p>The method for assigning a given value of PTRPARAM to a specific transport parameter, defined in Group 4, is discussed below. There are an arbitrary number of input lines each representing a given realization of parameter values. In a multiple-realization scenario, the code enters the input file for each realization, and for this input, reads down the corresponding number of lines to obtain the parameters for that realization. For example, for realization number 10, the code reads down to the 10th line of data (line 11 in the file) and uses those parameter values.</p>
<p>Once these parameters are read in for a given realization, they must be assigned to specific transport parameters. This is done in the following way in Group 4. If any of the inputs other than TRANSFLAG are negative, the code takes the absolute value of the number and interprets it as the column number from which to assign the transport parameter. For example, if DIFFMFL = -5, then the diffusion coefficient is the fifth input number in the PTRPARAM array. In this way, any of the transport parameter inputs can be assigned through the alternate input file rather than the input line in ptrk. It should be noted that for the colloid diversity model, only K_REV need be negative to indicate values should be read from the parameter file, if K_REV is negative then all five parameters are read from the file, otherwise the equation parameters will be read from the ptrk macro. This is to accommodate the fact that the SLOPE_KF may have negative values.</p>
<blockquote>
<div><p>KEYWORD ‘afm’</p>
<p>KEYWORD ‘dfree’</p>
</div></blockquote>
<p>Optional keyword “size” is input to indicate that the colloid size distribution model option is being used. It is followed by PART_SIZE and PROBSIZE.</p>
<blockquote>
<div><p>KEYWORD ‘size’</p>
<p>PART_SIZE(I), PROBSIZE(I) - an arbitrary numbers of lines of input, terminated by a blank line.</p>
</div></blockquote>
<p>Optional keyword “dive” is input to indicate that the colloid diversity model is being used. It is followed by optional keyword “file” and the name of the file containing the CDF table or equation data (a description of the format for the file is provided with the second mptr example), the TPRPFLAG with optional SIMNUM, optional CDF equation parameters (when “file” is not used), and keyword “irreversible” or “reversible” with FLAG_LOG.</p>
<table border="1" class="docutils">
<colgroup>
<col width="41%" />
<col width="6%" />
<col width="53%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><div class="first last line-block">
<div class="line">KEYWORD ‘dive’</div>
<div class="line">KEYWORD ‘file”</div>
<div class="line">CDFFILENAME TPRPFLAG</div>
<div class="line">or</div>
<div class="line">TPRPFLAG, SIMNUM</div>
<div class="line">KEYWORD ‘irreversible’</div>
<div class="line">or</div>
<div class="line">KEYWORD ‘reversible’,</div>
<div class="line">FLAG_LOG</div>
</div>
</td>
<td>or</td>
<td><div class="first last line-block">
<div class="line">KEYWORD ‘dive’</div>
<div class="line">TPRPFLAG</div>
<div class="line">or</div>
<div class="line">TPRPFLAG, SIMNUM</div>
<div class="line">K_REV, R_MIN,R_MAX, SLOPE_KF,</div>
<div class="line">CINT_KF</div>
<div class="line">KEYWORD ‘irreversible’</div>
<div class="line">or</div>
<div class="line">KEYWORD ‘reversible’, FLAG_LOG</div>
</div>
</td>
</tr>
</tbody>
</table>
<p>Note that optional KEYWORDs “size” and “dive” are only used when colloid transport is enabled.</p>
<p>There are eight possible forms of input for Group 4, which depend on whether or not the Active Fracture Model is implemented (optional KEYWORD “afm”), the free water diffusion coefficient is used (optional KEYWORD “dfree”), and colloid transport is enabled (TRANSFLAG &lt; 0). If colloid transport is enabled, then ABS(TRANSFLAG) is used to determine the transport mechanism.</p>
<p>If there is no colloid transport, TRANSFLAG &gt; 0</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), DIFFMAT(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC)</p>
<p>or with ‘dfree”’</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), H2O_DIFF(ITRC), TORT_DIFF(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC)</p>
<p>or when ‘afm’ is implemented</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), DIFFMAT(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC), SRESIDUAL(ITRC), GAMMA_AFM(ITRC)</p>
<p>or when ‘afm’ is implemented with ‘dfree’</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), H2O_DIFF(ITRC), TORT_DIFF(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC), SRESIDUAL(ITRC), GAMMA_AFM(ITRC)</p>
<p>Or when colloid transport is enabled, TRANSFLAG &lt; 0</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), DIFFMAT(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC), KCOLL(ITRC), RCOLL(ITRC), FCOLL(ITRC)</p>
<p>or with ‘dfree’</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), H2O_DIFF(ITRC), TORT_DIFF(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC), KCOLL(ITRC), RCOLL(ITRC), FCOLL(ITRC)</p>
<p>or when ‘afm’ is implemented</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), DIFFMAT(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC), KCOLL(ITRC), RCOLL(ITRC), FCOLL(ITRC), SRESIDUAL(ITRC), GAMMA_AFM(ITRC)</p>
<p>or when ‘afm’ is implemented with ‘dfree’</p>
<p>Group 4 - TRANSFLAG(ITRC), KD(ITRC), TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), H2O_DIFF(ITRC), TORT_DIFF(ITRC), RD_FRAC(ITRC), MATRIX_POR(ITRC), APERTURE(ITRC), KCOLL(ITRC), RCOLL(ITRC), FCOLL(ITRC), SRESIDUAL(ITRC), GAMMA_AFM(ITRC)</p>
<p>Group 5 - JA, JB, JC, ITRC</p>
<p>Group 6 - JA, JB, JC, PCNSK, T1SK, T2SK</p>
<p>If POUT = 5, an additional Group is included at the end of the ptrk input</p>
<p>Group 7 - JA, JB, JC, NODEPCONC (JA, JB, JC - defined on <a class="reference external" href="Macro20058.html">See JA, JB, JC, PROP1, PROP2, …</a>)</p>
<p>The concentration output is written to the “.trc”, “.out”, and AVS concentration output files. The “.fin” file is used only when specified (a non-zero value is input for PRNT_RST).</p>
<table border="1" class="docutils">
<colgroup>
<col width="2%" />
<col width="2%" />
<col width="96%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NPART</td>
<td>integer</td>
<td>Number of particles in the simulation. Note: the actual number may be slightly less than the number specified by the user because when the code divides the particles among the starting nodes as specified in Group 7, the code must input an integer number of particles at each node.</td>
</tr>
<tr class="row-odd"><td>RSEED</td>
<td>integer</td>
<td>6-digit integer random number seed.</td>
</tr>
<tr class="row-even"><td>RSEED_RELEASE</td>
<td>integer</td>
<td>6-digit integer random number seed for particle release location calculation. If a value is not entered for RSEED_RELEASE it will be set equal to RSEED.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>Note that for GoldSim-FEHM coupled simulations the random seeds are controlled by GoldSIM and the values input in the ptrk macro are not used.</td>
</tr>
<tr class="row-even"><td>MAX1D</td>
<td>integer</td>
<td>Maximum 1-D array size for holding particle tracking information for all simulated species. The value of MAX1D depends on number of species, number of time steps, number of radionuclide release bins, number of species involved in ingrowth, and the length of the decay-ingrowth chain.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘wtri” indicatiing a water table rise calculation should be performed.</td>
</tr>
<tr class="row-even"><td>WATER_TABLE</td>
<td>real</td>
<td>Water table elevation to be used for water table rise calculation.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>Note that for GoldSim-FEHM coupled simulations the water table rise calculations are controlled by GoldSIM and the values input in the ptrk macro are not used and may be omitted.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘rip’ designating this is a GoldSim coupled simulation and the alternate Group 2 data format should be used.</td>
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
<tr class="row-odd"><td>RIPFEHM</td>
<td>integer</td>
<td>Parameter for assigning the particle starting locations based on radionuclide flux input from the computer code GoldSim. Used when GoldSim is the driver program and FEHM is dynamically linked to perform particle tracking transport. If RIPFEHM = 1 Use input from GoldSim to assign particle starting locations.If RIPFEHM ≠ 1 Assign starting locations in the normal way.</td>
</tr>
<tr class="row-even"><td>CONFREAD</td>
<td>real</td>
<td>Initial conversion factor for GoldSim-FEHM coupling (# of particles/mole). If CONFREAD=0, at each time step, the code selects a conversion factor based on the available memory and the remaining simulation time (end time - current time). The code then uses the selected conversion factor to calculate the number of particles to be injected at the current time step. If CONFREAD &gt;0, the code uses the product of the CONFREAD and the input mass (in moles) to calculate the input number of particles at each time step.</td>
</tr>
<tr class="row-odd"><td>GMOL</td>
<td>real</td>
<td>The molecular weight of the ith species. The code uses GMOL and CONFREAD to convert the mass from number of particles to grams in the final output.</td>
</tr>
<tr class="row-even"><td>P_FRACTION</td>
<td>real</td>
<td>The decay-ingrowth particle release factor (percentage of the maximum number of particles released for the parent species). Values are in the range 0. - 1. If the value is omitted it will default to 0.5.</td>
</tr>
<tr class="row-odd"><td>TRAK_TYPE</td>
<td>integer</td>
<td>Flag to denote the fluid phase of the particles:1 - liquid phase particles2 - vapor phase particles</td>
</tr>
<tr class="row-even"><td>HALF_LIFE</td>
<td>real</td>
<td>Half-life for irreversible first order decay reaction (days). Set HALF_LIFE = 0 for no decay.</td>
</tr>
<tr class="row-odd"><td>POUT</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to specify the concentration output:0 - Concentration output is a running total of the number of particles which have left each node, divided by the fluid or vapor mass at that node depending on trak_type.</div>
<div class="line">1 - Concentrations computed as number of particles per unit total volume (rock and fluid)</div>
<div class="line">2 - Concentrations computed as number of particles per unit fluid volume (the fluid is liquid for TRAK_TYPE = 1 and gas for TRAK_TYPE = 2).</div>
<div class="line">3 - Concentrations computed as number of particles at a given node point.</div>
<div class="line">5 - If this option is invoked, the particles injected at a particular node are assigned concentrations according to the input concentration defined in the NODEPCONC array (Group 7). The code then outputs the mixed mean concentration at each node in the model based on the assumption of steady state flow.</div>
<div class="line">6 - Used for C-14 radioactive decay particle mixing model (only liquid tracer). For meaningful results the particles must all be injected simultaneously in a pulse (give a very short duration of injection starting at time 0). The code contains data describing the function f(t) vs. time where f(t) is given as: <span class="math notranslate nohighlight">\(\int_0^t \exp(-kt)dt\)</span> where <span class="math notranslate nohighlight">\(t\)</span> is the time the particle enters the system and <span class="math notranslate nohighlight">\(k\)</span> is the radioactive decay constant for C-14. The output is the final concentration after all the particles have left the system.</div>
<div class="line">-1, -2, -3, or -6 - Concentrations computed as specified above for abs(pout). The “.trc” file contains breakthrough output for the first node specified in the node macro.</div>
<div class="line">-7 - Output is written every time a particle leaves a cell. This output is particle number, cell number that the particle is leaving, zone number of the cell, and time the particle leaves the cell.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>PRNT_RST</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag to specify whether particle information is written to the “.fin”, “.ptrk_fin”, or “.ptrk” files:</div>
<div class="line">If PRNT_RST = 0, Particle information is not written to the output files.</div>
<div class="line">If PRNT_RST = 1, 11, 21, 31, 41, All particle information necessary for a restart is written to the “.fin” file.</div>
<div class="line">If PRNT_RST = -1, -11, -21, -31, -41, Only particle positions and ages are written to the “.fin” file.</div>
<div class="line">If ABS (PRNT_RST) = 2, 12, 22, 32 or 42, Mass flux values are written to the “.fin” file followed by particle information.</div>
<div class="line">If 10 ≤ ABS(PRNT_RST) &lt; 30, Particle exit locations and count are written to the “.ptrk_fin” file.</div>
<div class="line">If ABS(PRNT_RST) ≥ 20, Cumulative particle counts are written to the “.ptrk” file, for variables specified by PRNT_VARNUM (the default is to output all variables).</div>
<div class="line">If ABS(PRNT_RST) ≥ 40, Cumulative mass output from a FEHM/GoldSim coupled simulation will be written to file FEHM_GSM_Mass_balance.txt. Note that to track cumulative mass an additional array of size maxparticles*nspeci must be allocated so caution should be used when specifying this option to ensure sufficient system memory is available.</div>
<div class="line">When particle tracking data or mass fluxes are written to the “.fin” file, the arrays are written after all of the heat and mass simulation information. The mass fluxes can be read into the code in a subsequent ptrk or mptr simulation and the code can simulate transport on this steady state flow field (see macro rflo).The particle information written is sufficient to perform a restart of the particle tracking simulation and to post-process the data to compile statistics on the particle tracking run.</div>
<div class="line">However, for a large number of particles, this file can become quite large, so particle tracking information should only be written when necessary. Thus, 0 should be used for PRNT_RST unless restarting or post-processing to obtain particle statistics is required. Selecting the -1 option allows a subset of the full set of information needed for a restart (particle positions and ages) to be written. Restart runs that use this file as input will only be approximate, since the particle is assumed to have just entered its current cell. For restart runs, PRNT_RST = 1 is preferred, while PRNT_RST = -1 is appropriate for output of particle statistics for post-processing.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>PRNT_VARNUM</td>
<td>integer</td>
<td>A list of integers specifying which particle counts should be output. For each value entered PRNT_VAR(PRNT_VARNUM) is set to .true. If no values are entered the default is to print all variables. 1 – Number of particles that have entered the system 2 – Number of particles currently in the system 3 – Number of particles that have left the system 4 – Number of particles that have decayed 5 – Number of particles that have been filtered 6 – Number of particles that left this time interval</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>Note: The data found in the “.ptrk” file was previously reported in the general output file. From version 2.25 of the code and forward that data will be reported in the optional, “.ptrk” file unless a coupled GoldSim-FEHM simulation is being run. In addition, the user has the option of selecting which statistics parameters are reported. The default is to report all statistics parameters.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “tcurve” indicating transfer function curve data should be input to model matrix diffusion. If the keyword is found then NUMPARAMS and FILENAME are entered, otherwise they are omitted.</td>
</tr>
<tr class="row-even"><td>NUMPARAMS</td>
<td>integer</td>
<td>Number of parameters that define the transfer function curves being used.</td>
</tr>
<tr class="row-odd"><td>TFILENAME</td>
<td>character</td>
<td>Name of input file containing the transfer function curve data.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘zptr’ designating zones for breakthrough curves will be defined. If no keyword is input, IPZONE and IDZONE are also omitted.</td>
</tr>
<tr class="row-odd"><td>IPZONE</td>
<td>integer</td>
<td>Number of zones for which breakthrough curves are to be output</td>
</tr>
<tr class="row-even"><td>IDZONE</td>
<td>integer</td>
<td>A list of zones for which particle breakthrough data are required.The code outputs the number of particles that leave the system at each zone IDZONE at the current time step. This information is written to the “.out” file at each heat and mass transfer time step.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating alternate transport parameter file input for multiple simulation realizations.</td>
</tr>
<tr class="row-even"><td>PFILENAME</td>
<td>character*80</td>
<td>Name of file from which to read transport parameters.</td>
</tr>
<tr class="row-odd"><td>NINPUTS</td>
<td>integer</td>
<td>Number of inputs in each row of data in the alternate transport parameter input file (PFILENAME).</td>
</tr>
<tr class="row-even"><td>PTRPARAM</td>
<td>real</td>
<td>Array of transport parameter values in the alternate transport parameter input file (PFILENAME). The parameters that may be input are those entered for Group 4. Only those parameters being changed at each realization need be entered.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘afm’ designating the Active Fracture Model input for determining fracture spacing should be used.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*5</td>
<td>Optional keyword ‘dfree’ designates that the free water diffusion coefficient and tortuosity will be input instead of the molecular diffusion coefficient.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘size’ designating that the colloid size distribution model option is being used (combined with the interface filtration option in the itfc macro). If the keyword is not input, PART_SIZE and PROBSIZE are also omitted.</td>
</tr>
<tr class="row-even"><td>PART_SIZE</td>
<td>real</td>
<td>Colloid particle size for this entry of the particle size distribution table (paired with a value of PROBSIZE). An arbitrary number of entries can be input, terminated with a blank line. The code assigns each particle a size based on this distribution of particle sizes, and decides if particles are irreversibly filtered based on the pore size distribution assigned in the itfc macro.</td>
</tr>
<tr class="row-odd"><td>PROBSIZE</td>
<td>real</td>
<td>Colloid cumulative probability for the distribution of sizes (paired with a value of PART_SIZE). See description of PART_SIZE above for details. The final entry of the table must have PROBSIZE = 1, since the distribution is assumed to be normalized to unity.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword “dive” signifying that the specie being specified is either a colloid species using the colloid diversity model or a non-colloid daughter species of a colloid species.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating the cumulative probability distribution function (CDF) retardation parameters for the colloid diversity model should be read from an external file.</td>
</tr>
<tr class="row-even"><td>CDF_FILENAME</td>
<td>character*80</td>
<td><div class="first last line-block">
<div class="line">Name of the file containing the cumulative probability distribution function (CDF) (entered if optional keyword ‘file’ follows keyword ‘dive’).</div>
<div class="line">If TPRPFLAG = 11 or 12, Table option</div>
<div class="line">If TPRPFLAG = 13 or 14, Equation option</div>
<div class="line">The following equations are used for <span class="math notranslate nohighlight">\(R_{min} \le R \le R_{max}\)</span>, <span class="math notranslate nohighlight">\(R = 1 + K_f / K_{rev}\)</span>, <span class="math notranslate nohighlight">\(\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>TPRP_FLAG</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Values of TPRPFLAG between 11 and 14 signify that the colloid diversity model with equal weight sampling will be used:</div>
<div class="line">TPRPFLAG = 11: CDF vs retardation factor specified in a table</div>
<div class="line">TPRPFLAG = 12: similar to 11, but the SQRT(CDF) is used instead of CDF for sampling</div>
<div class="line">TPRPFLAG = 13: CDF vs <span class="math notranslate nohighlight">\(K_f\)</span> (Attachment rate constant) specified as a straight line equation in the log-log space</div>
<div class="line">TPRPFLAG = 14: similar to 13, but the SQRT(CDF) is used instead of CDF for sampling</div>
</div>
</td>
</tr>
<tr class="row-even"><td>SIMNUM</td>
<td>integer</td>
<td>Simulation number, used for selecting the table/equation from the colloid diversity file. For GoldSim-FEHM coupled simulations or FEHM runs using the ‘msim’ option this parameter is passed to the code. For non-coupled simulations it is an optional input. (Default value = 1)</td>
</tr>
<tr class="row-odd"><td>K_REV</td>
<td>real</td>
<td>Detachment rate constant for reversible filtration of irreversible colloids.</td>
</tr>
<tr class="row-even"><td>R_MIN</td>
<td>real</td>
<td>Minimum value of the retardation factor for reversible filtration of irreversible colloids.</td>
</tr>
<tr class="row-odd"><td>R_MAX</td>
<td>real</td>
<td>Maximum value of the retardation factor for reversible filtration of irreversible colloids</td>
</tr>
<tr class="row-even"><td>SLOPE_KF</td>
<td>real</td>
<td>Value of the slope (<span class="math notranslate nohighlight">\(m\)</span>) in the log-log space for the equation: <span class="math notranslate nohighlight">\(\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)\)</span></td>
</tr>
<tr class="row-odd"><td>CINT_KF</td>
<td>real</td>
<td>Value of the intercept (//b//) in the log-log space for the above equation</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Keyword specifying whether the colloid species is ‘irreversible’ or ‘reversible’.</td>
</tr>
<tr class="row-odd"><td>FLAG_LOG</td>
<td>integer</td>
<td>For reversible colloids an average retardation factor is used: If FLAG_LOG = 0: a linear average of the distribution is usedIf FLAG_LOG = 1: a log-linear average of the distribution is used</td>
</tr>
<tr class="row-even"><td>TRANSFLAG</td>
<td>integer</td>
<td>Flag to specify which transport mechanisms apply [abs(TRANSFLAG)]:1 - advection only (no dispersion or matrix diffusion)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>2 - advection and dispersion (no matrix diffusion)</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>3 - advection and matrix diffusion, infinite fracture spacing solution (no dispersion)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>4 - advection, dispersion, and matrix diffusion, infinite fracture spacing solution</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>5 - advection and matrix diffusion, finite fracture spacing solution (no dispersion)</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>6 - advection, dispersion, and matrix diffusion, finite fracture spacing solution</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>8 - use the the transfer function approach with 3 dimensionless parameters and type curves for handling fracture-matrix interactions.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>For TRANSFLAG &lt; 0, transport simulations include colloids.</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>For equivalent continuum solutions, the fracture spacing in the finite spacing model is determined using <span class="math notranslate nohighlight">\(SPACING = APERTURE / POROSITY\)</span>.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>For dual permeability models, the fracture spacing input parameter APUV1 in the dpdp macro is used as the half-spacing between fractures. If the Active Fracture Model (see keyword ‘afm’) is used, APUV1 is the geometric fracture half-spacing, and the additional terms SRESIDUAL and GAMMA_AFM are used to determine the spacing between active fractures (see below).</td>
</tr>
<tr class="row-even"><td>KD</td>
<td>real</td>
<td>Sorption coefficient (linear, reversible, equilibrium sorption). Units are kg-fluid / kg-rock (these units are equivalent to the conventional units of cc/g when the carrier fluid is water at standard conditions). This value applies to the medium as a whole when matrix diffusion is turned off, whereas for simulations invoking matrix diffusion, the value applies to the rock matrix. For the latter case, sorption in the flowing system (fractures) is modeled using the RD_FRAC variable.</td>
</tr>
<tr class="row-odd"><td>TCLX</td>
<td>real</td>
<td>Dispersivity in the x-direction (m). The input value is ignored when dispersion is turned off.</td>
</tr>
<tr class="row-even"><td>TCLY</td>
<td>real</td>
<td>Dispersivity in the y-direction (m). The input value is ignored when dispersion is turned off.</td>
</tr>
<tr class="row-odd"><td>TCLZ</td>
<td>real</td>
<td>Dispersivity in the z-direction (m). The input value is ignored when dispersion is turned off.</td>
</tr>
<tr class="row-even"><td>DIFFMAT</td>
<td>real</td>
<td>Molecular diffusion coefficient in the rock matrix (m2/s). The input value is ignored unless matrix diffusion is invoked.</td>
</tr>
<tr class="row-odd"><td>RD_FRAC</td>
<td>real</td>
<td>Retardation factor within the primary porosity (fractures) for a matrix diffusion particle tracking simulation (use 1 for no sorption on fracture faces). The input value is ignored unless matrix diffusion is invoked.</td>
</tr>
<tr class="row-even"><td>MATRIX_POR</td>
<td>real</td>
<td>Porosity of the rock matrix. Used to simulate diffusion and sorption in the rock matrix when matrix diffusion is invoked, otherwise the input value of MATRIX_POR is ignored.</td>
</tr>
<tr class="row-odd"><td>APERTURE</td>
<td>real</td>
<td>Mean fracture aperture (m). The input value is ignored when matrix diffusion is turned off.</td>
</tr>
<tr class="row-even"><td>KCOLL</td>
<td>real</td>
<td>Colloid distribution parameter, the ratio of contaminant mass residing on colloids to the mass present in aqueous form. It is used to compute an effective aperture via the following: <span class="math notranslate nohighlight">\(APWID = APERRTURE \cdot (1 + KCOLL)\)</span></td>
</tr>
<tr class="row-odd"><td>RCOLL</td>
<td>real</td>
<td>Colloid retardation factor. Used, in conjunction with kcoll, to adjust colloid retardation in fractures using the following formula: <span class="math notranslate nohighlight">\(FRACRD = \frac{RD_FRAC + KCOLL \cdot ROLL}{1 + KCOLL}\)</span></td>
</tr>
<tr class="row-even"><td>FCOLL</td>
<td>real</td>
<td>Colloid filtration parameter. Used to compute the probability a colloid will be irreversibly filtered along the path between two nodes using the following: <span class="math notranslate nohighlight">\(PROBFILT = 1 - \exp(DISTANCE / FCOLL)\)</span>&nbsp;where <span class="math notranslate nohighlight">\(DISTANCE\)</span> is the length of the path between nodes.</td>
</tr>
<tr class="row-odd"><td>SRESIDUAL</td>
<td>real</td>
<td>Residual saturation in the Active Fracture Model used for determining the spacing between active fractures. This parameter is only needed when the keyword ‘afm’ is included, in which case the input must be entered. However, the model is only used in dual permeability simulations at locations where the finite spacing matrix diffusion model is invoked [abs(TRANSFLAG) = 5 or 6].</td>
</tr>
<tr class="row-even"><td>GAMMA_AFM</td>
<td>real</td>
<td>Exponent in the Active Fracture Model used for determining the spacing between active fractures. See comments for SRESIDUAL above.</td>
</tr>
<tr class="row-odd"><td>ITRC</td>
<td>integer</td>
<td>Model number for parameters defined in group 4. Default is 1.</td>
</tr>
<tr class="row-even"><td>PCNSK</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Particle injection parameter assigned for nodes defined by JA, JB, and JC. When multiple lines of input are given for Group 6, all PCNSK values must have the same sign (i.e. the two options, described below, cannot be invoked in the same simulation).</div>
<div class="line">PCNSK &gt; 0 - particles are injected at each node in proportion to the source mass flow rate at the node. When multiple lines of input are given for Group 6, PCNSK is proportional to the particle injection concentration. This boundary condition is equivalent to injecting a solute of a given concentration into the system. Note: the source flow rates used to assign the number and timing of particle injections are those at the beginning of the particle tracking simulation (time DAYCS). Transient changes in this source flow rate during the particle tracking simulation do not change the input of particles to the system.</div>
<div class="line">PCNSK &lt; 0 - particles are introduced at the node(s), regardless of whether there is a fluid source at the node. When multiple lines of input are given for Group 6, abs(PCNSK) is proportional to the number of particles introduced at the node(s).Default is 0 for all unassigned nodes, meaning that no particles are injected at that node.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>T1SK</td>
<td>real</td>
<td>Time (days) when particle injection begins. Default is 0.</td>
</tr>
<tr class="row-even"><td>T2SK</td>
<td>real</td>
<td>Time (days) when particle injection ends. Default is 0.</td>
</tr>
<tr class="row-odd"><td>NODEPCONC</td>
<td>real</td>
<td>Input particle concentrations. The concentration associated with a particle entering the system at a specified node.</td>
</tr>
</tbody>
</table>
<p>The following is an example of ptrk:</p>
<p>In this example, 100,000 nondecaying, liquid-borne particles are introduced as a
sharp pulse (from time 10 to 10.0001 days) with the injection fluid in zone 3
(an injection well defined in the zone macro preceding ptrk). The particle
tracking simulation starts as the heat and mass transfer simulation is turned
off at day 10, after having established a fluid flow steady state. Two models
are defined for assigning transport properties of the particles. All nodes are
assigned to model 1, after which model 2 properties are assigned for zone 2.
A combined advection, dispersion, and matrix diffusion model is used for all nodes.
However, sorption in the matrix occurs only for model 2 (which is zone 2 in this
simulation), and the matrix transport properties (porosity, fracture spacing,
diffusion coefficient) differ for this model as well.</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="14%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="16%" />
<col width="7%" />
<col width="11%" />
<col width="12%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>ptrk</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>100000</td>
<td>122945</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td><ol class="first last arabic simple" start="10">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="10">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="20">
<li></li>
</ol>
</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0</td>
<td>20</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>4</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td>5.e-11</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>0.1</td>
<td>0.333</td>
</tr>
<tr class="row-even"><td>4</td>
<td><ol class="first last arabic simple" start="3">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
<td>1.e-10</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>0.28</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
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
</tr>
<tr class="row-even"><td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>-2</td>
<td>0</td>
<td>0</td>
<td>2</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
</tr>
<tr class="row-odd"><td>-3</td>
<td>0</td>
<td>0</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="10">
<li></li>
</ol>
</td>
<td>10.0001</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
</tr>
</tbody>
</table>
<p>In the second example, transfer function data is used along with the active
fracture and colloid diversity models. Particle statistics data for the cumulative
number of particles that have left the sytem and the number of particles that
left during the current timestep are output.</p>
<table border="1" class="docutils">
<colgroup>
<col width="32%" />
<col width="6%" />
<col width="6%" />
<col width="6%" />
<col width="7%" />
<col width="6%" />
<col width="4%" />
<col width="7%" />
<col width="5%" />
<col width="2%" />
<col width="4%" />
<col width="4%" />
<col width="4%" />
<col width="8%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>ptrk</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<tr class="row-even"><td>100000</td>
<td>244562</td>
<td>244562</td>
<td>&#160;</td>
<td>&#160;</td>
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
<tr class="row-odd"><td><ol class="first last arabic simple" start="0">
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0.230</td>
<td>3</td>
<td>6</td>
<td>&#160;</td>
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
<tr class="row-odd"><td>tcurve</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<tr class="row-odd"><td>../colloid_cell/input/uz_tfcurves_nn_3960.in</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<tr class="row-even"><td>afm</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<tr class="row-odd"><td>diversity</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>../colloid_cell/input/rcoll_equation.dat</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<tr class="row-even"><td>13 5</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>-2</td>
<td>0.0</td>
<td>1.e-03</td>
<td>1.e-03</td>
<td>1.00e-30</td>
<td>1.0e+00</td>
<td>0.2</td>
<td>1.e-03</td>
<td>1e+20</td>
<td>1</td>
<td>1.00</td>
<td>0.01</td>
<td>0.6</td>
<td>//layer 1</td>
</tr>
<tr class="row-odd"><td>-2</td>
<td>0.0</td>
<td>1.e-03</td>
<td>1.e-03</td>
<td>1.00e-03</td>
<td>1.0e+00</td>
<td>0.2</td>
<td>0.00e+00</td>
<td>1e+20</td>
<td>1</td>
<td>1.00</td>
<td>0.00</td>
<td>0.0</td>
<td>//layer 2</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>1</td>
<td>1</td>
<td>-1. 0.</td>
<td>0.01</td>
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
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>With the file <code class="docutils literal notranslate"><span class="pre">rcoll_equation.dat</span></code> as:</p>
<table border="1" class="docutils">
<colgroup>
<col width="10%" />
<col width="13%" />
<col width="13%" />
<col width="23%" />
<col width="23%" />
<col width="18%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td colspan="6">test file for 3Dmptr equal weight sampling</td>
</tr>
<tr class="row-even"><td>1</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td>3</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="5">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="5">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-even"><td>4</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="5">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="5">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td>5</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="9">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>-0.9</td>
</tr>
<tr class="row-even"><td>6</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="5">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple" start="85">
<li></li>
</ol>
</td>
<td>1.16667</td>
<td>-2.25</td>
</tr>
<tr class="row-odd"><td>8</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>100000</td>
<td><ol class="first last arabic simple">
<li></li>
</ol>
</td>
<td>-5.</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
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
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">sptr</span></code><a class="headerlink" href="#sptr" title="Permalink to this headline">¶</a></h1>
<p>Streamline particle tracking is invoked with this control statement.</p>
<p>Group 1 -       DTMX, IPRT, IPRTO, RSEED, TPLIM, DISTLIM, LINELIM, DIVD_WEIGHT</p>
<p>Group 2 -       COURANT_FACTOR, IPRTR, ITENSOR, IREVERS, FREEZ_TIME, MAX_JUMP</p>
<p>Keywords and their associated input are described below. These keywords may be entered in any order following Group 2, but must be directly followed by any associated data. All keywords are optional except “tprp” which flags input of the particle transport properties.</p>
<p>Optional keyword “tcurve” is input to indicate that transfer function curves should be input to model matrix diffusion. It is followed by NUMPARAMS and TFILENAME.</p>
<ul class="simple">
<li>KEYWORD “tcurve”</li>
<li>NUMPARAMS</li>
<li>TFILENAME</li>
</ul>
<p>Optional keyword “omr” is input to indicate the grid has octree mesh refinement. It is followed by the number of refined (OMR) nodes in the grid. An optional keyword, “file”, which if input is followed by the name of a file from which to read or to write omr initialization calculation arrays, and the file format. The “file” option should only be used for particle tracking runs using steady state flow. It should also be noted that the file should be re-generated if parameter changes that affect velocity are made since velocities will not be recalculated if the file exists.</p>
<ul class="simple">
<li>KEYWORD “omr”</li>
<li>OMR_NODES</li>
<li>KEYWORD “file”</li>
<li>OMRFILENAME</li>
<li>OMR_FORM</li>
</ul>
<p>Optional keyword “tpor” is input to indicate tracer porosity will be input and is followed by tracer porosity values input using JA, JB, JC format or the keyword “file” and the name of a file containing the tracer porosity values.</p>
<ul class="simple">
<li>KEYWORD “tpor”</li>
<li>JA, JB, JC, PS_TRAC (JA, JB, JC - defined on page 33)</li>
</ul>
<p>-or-</p>
<ul class="simple">
<li>KEYWORD “tpor”</li>
<li>KEYWORD “file”</li>
<li>TPORFILENAME</li>
</ul>
<p>Optional keyword “wtdt” to indicate initial particle locations should be moved below the water table by distance DELTAWT.</p>
<ul class="simple">
<li>KEYWORD “wtdt”, DELTAWT</li>
</ul>
<p>Optional keyword “volum” to indicate control volumes associated with computation of sptr velocities should be written to an output file. These volumes are used with PLUMECALC to account for the approximate control volumes used for the velocity calculations on an OMR grid. The “volume” keyword may be followed by the file format.</p>
<ul class="simple">
<li>KEYWORD “volum”</li>
</ul>
<p>-or-</p>
<ul class="simple">
<li>KEYWORD “volum”, SPTRX_FORMAT</li>
</ul>
<p>Optional keywords “po”, “sa”, “pe”, “de”, “pr”, “te”, “zo”, “id” define parameters to be output. No data is associated with the parameter output flags. These parameters are output in the ‘<a href="#id1"><span class="problematic" id="id2">*</span></a>.sptr2” file in the order entered.</p>
<p>Optional keyword “xyz” indicates that coordinate data should be included in the abbreviated ‘<a href="#id3"><span class="problematic" id="id4">*</span></a>.sptr2’ output file (see IPRTO below).</p>
<p>The optional keyword “zbtc” indicates that breakthrough curves will be computed for specified zones. It is followed by NZBTC … and ZBTC. If keyword “alt” follows the “zbtc” keyword an alternate output format will be used where,</p>
<ul class="simple">
<li>KEYWORD “zbtc”</li>
</ul>
<p>-or-</p>
<ul class="simple">
<li>KEYWORD “zbtc” “alt”</li>
<li>NZBTC, DTMN, PART_MULT, PART_FRAC, DELTA_PART, PART_STEPS, TIME_BTC</li>
<li>ZBTC</li>
</ul>
<p>Optional keyword “cliff”</p>
<p>Optional keyword “corner”</p>
<p>Optional keyword “capture”</p>
<p>Optional keyword “spring”</p>
<p>Keyword (‘tprp’) specifies that the particle transport properties will be entered on subsequent lines. The transport properties input (Group 4) depends on the value of TPRP_FLAG, the first value input for that group. For TPRP_FLAG = 2 or 4, format of Group 4 input also depends on the form of the dispersion coefficient tensor, as selected using the flag ITENSOR (Group 2).</p>
<ul class="simple">
<li>KEYWORD “tprp”</li>
</ul>
<p>TPRP_FLAG = 1:</p>
<p>Group 3 -       TPRP_FLAG, KD</p>
<p>TPRP_FLAG = 2:</p>
<p>ITENSOR = 1</p>
<p>Group 3 -       TPRP_FLAG, KD, DM, A1, A2, A3, A4, ASX, ASY, VRATIO</p>
<p>ITENSOR = 2</p>
<p>Group 3 -       TPRP_FLAG, KD, DM, AL, ATH, ATV, VRATIO</p>
<p>ITENSOR = 3</p>
<p>Group 3 -       TPRP_FLAG, KD, DM, ALH, ALV, ATH, ATV, VRATIO</p>
<p>ITENSOR = 4</p>
<p>Group 3 -       TPRP_FLAG, KD, DM, AL, AT, VRATIO</p>
<p>ITENSOR = 5</p>
<p>Group 3 -       TPRP_FLAG, KD, AL, ATH, ATV, DM, VRATIO</p>
<p>TPRP_FLAG = 3:</p>
<p>Group 3 -       TPRP_FLAG, KD, DIFM, RD_FRAC, POR_MATRIX, APERTURE</p>
<p>TPRP_FLAG = 4:</p>
<p>ITENSOR = 1</p>
<p>Group 3 -       TPRP_FLAG, KD, DIFM, RD_FRAC, POR_MATRIX, APERTURE, DM, A1, A2, A3, A4, ASX, ASY, VRATIO</p>
<p>ITENSOR = 2</p>
<p>Group 3 -       TPRP_FLAG, KD, DIFM, RD_FRAC, POR_MATRIX, APERTURE, DM, AL, ATH, ATV, VRATIO</p>
<p>ITENSOR = 3</p>
<p>Group 3 -       TPRP_FLAG, KD, DIFM, RD_FRAC, POR_MATRIX, APERTURE, DM, ALH, ALV, ATH, ATV, VRATIO</p>
<p>ITENSOR = 4</p>
<p>Group 3 -       TPRP_FLAG, KD, DIFM, RD_FRAC, POR_MATRIX, APERTURE, DM, AL, AT, VRATIO</p>
<p>ITENSOR = 5 (not recommended, included for compatibility with older versions of the code)</p>
<p>Group 3 -       TPRP_FLAG, KD, DIFM, RD_FRAC, POR_MATRIX, APERTURE, AL, ATH, ATV, DM, VRATIO</p>
<p>TPRP_FLAG =11 or TPRP_FLAG = 12:</p>
<p>Group 3-        TPRP_FLAG, SIMNUM</p>
<ul class="simple">
<li>KEYWORD ‘file’</li>
<li>CDF_FILENAME</li>
</ul>
<p>TPRP_FLAG =13 or TPRP_FLAG = 14:</p>
<p>Group 3 -       TPRP_FLAG, SIMNUM</p>
<ul class="simple">
<li>KEYWORD ‘file’</li>
<li>CDF_FILENAME</li>
</ul>
<p>or</p>
<p>Group 3 -       TPRP_FLAG, K_REV, R_MIN, R_MAX, SLOPE_KF, CINT_KF, AL, ATH, ATV, DM, VRATIO</p>
<p>Group 4 -       JA, JB, JC, MODEL_NUMBER</p>
<p>Group 3 is ended when a blank line is encountered. The MODEL_NUMBER is incremented each time a Group 3 line is read, and Group 4 lines refer to this parameter.</p>
<p>Group 5 -       ITM, IST</p>
<p>-or-</p>
<p>Group 5 -       ITM, IST, COUNT_STEPS_MAX, SPTR_FLAG</p>
<p>Group 6 -       NX, NY, NZ</p>
<p>Group 7 -       X10, Y10, Z10</p>
<p>Group 8 -       XDIM, YDIM, ZDIM</p>
<p>Group 9 -       IJKV(I), X1(I), Y1(I), Z1(I) for I = 1 to NUMPART</p>
<p>Group 9 input is terminated with a blank line.</p>
<p>-or-</p>
<p>Group 9 -       KEYWORD “file”</p>
<ul class="simple">
<li>SPTR_FILENAME</li>
</ul>
<p>Restart runs for particle tracking may be accomplished by reading particle starting locations from a file. A particle restart file is generated by adding the optional SPTR_FLAG, keyword “save” to group 5.</p>
<p>Note when IST = 0 or 1, Group 10 is used and place holders are inserted for Groups 7-9 and NUMPART is equal to the number of particle starting locations that are entered; however, when IST = 2, Group 10 is not implemented and Groups 7-9 are used followed by a blank line. NUMPART equals NX*NY*NZ in this case.</p>
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
<tr class="row-even"><td>DTMX</td>
<td>real</td>
<td>Time step control (seconds). FEHM will account for all particles every abs(dtmx) seconds and write information to the “.sptr3” file if the “zbtc” keyword is present. This controls the output density for breakthrough curve information only. If you are not using/creating breakthrough curves, set DTMX very large (e.g. 1e20). If DTMX is negative, the time step for streamline calculations is forced to be abs(DTMX) seconds.</td>
</tr>
<tr class="row-odd"><td>IPRT</td>
<td>integer</td>
<td>Flag to denote whether individual particle positions are written at specified intervals to the “.sptr1” file. The particle coordinate positions are used to get a snapshot of the particle plume at various times during the simulation.IPRT = 0, No outputIPRT &gt; 0, Output is written to the “.sptr1” file every IPRT time steps.</td>
</tr>
<tr class="row-even"><td>IPRTO</td>
<td>integer</td>
<td>Flag to denote if particle streamline information is written to the “.sptr2” file. The information is used to draw complete particle streamlines (for a relatively small number of particles).IPRTO = 0, No outputIPRTO &gt; 0, Extended output is written to the “.sptr2” file.IPRTO &lt; 0, Abbreviated output is written to the “.sptr2” file.If abs(IPRTO) = 1,output is formatted, abs(IPRTO) = 2 output is unformatted, and abs(IPRTO) = 3 output is in binary format.</td>
</tr>
<tr class="row-odd"><td>RSEED</td>
<td>integer</td>
<td>Random number seed for the random number generator. For compatibility with earlier versions of FEHM in which this input did not exist, if no value of RSEED is input, the code assigns a value of 466201.</td>
</tr>
<tr class="row-even"><td>TPLIM</td>
<td>real</td>
<td>Minimum amount of time (days) a particle should move before location is output. Default is 0 days.</td>
</tr>
<tr class="row-odd"><td>DISTLIM</td>
<td>real</td>
<td>Minimum distance (m) a particle should move before location is output. Default is 0 meters.</td>
</tr>
<tr class="row-even"><td>LINELIM</td>
<td>integer</td>
<td>Maximum number of lines that will be written to sptr2 file.</td>
</tr>
<tr class="row-odd"><td>DIVD_WEIGHT</td>
<td>real</td>
<td>Weight factor for the derivative of the dispersion tensor term. Default is 1. If a value of zero is entered, the derivative term is not used.</td>
</tr>
<tr class="row-even"><td>COURANT_FACTOR</td>
<td>integer</td>
<td>Fraction of the distance through a cell that a particle moves in a single time step. This is used to ensure that the particle, on average, traverses less than one cell before a random-walk dispersion step is performed. For example, a factor of 0.25 indicates that the particle should take at least 4 time steps to move through a cell.</td>
</tr>
<tr class="row-odd"><td>IPRTR</td>
<td>integer</td>
<td>Flag for choosing the method for computing concentrations in cells based on the particle tracking information that will be written to the “.trc” or AVS output files.IPRTR Š 0, particle concentrations are computed as number of particles residing in the cell divided by the fluid mass in the cell.IPRTR &lt; 0, an integral of the particle concentration specified above is made and reported. This integral is the normalized cumulative concentration, which for a steady state flow field is equivalent to the response to a step change in particle concentration (note that the particles are input as a pulse).</td>
</tr>
<tr class="row-even"><td>ITENSOR</td>
<td>integer</td>
<td>Flag indicating the mathematical form of the dispersion coefficient tensor to be selected.ITENSOR = 0, No dispersion.ITENSOR = 1, Generalized form of the axisymmetric tensor, from Lichtner et al. (2002)ITENSOR = 2, Axisymmetric form of the dispersion coefficient tensor of Burnett and Frind (1987)ITENSOR = 3, Modified form of the dispersion coefficient tensor of Burnett and Frind (1987). See Lichtner et al. (2002) for detailsITENSOR = 4, Isotropic form of the dispersion coefficient tensor of Tompson et al. (1987)ITENSOR = 5, Original form of the Burnett and Frind (1987) tensor as implemented in FEHM Version 2.10.</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>Note: for Version 2.10 and earlier, the variable ITENSOR did not exist. For compatibility with these earlier versions, when ITENSOR is omitted from the input file, the code uses the ITENSOR = 5 formulation and the pre-existing input format. It is recommended that new simulations use one of the other tensor formulations (ITENSOR = 1 to 4).</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>In addition, the sign of ITENSOR is used as a switch as follows: if ITENSOR &lt; 0, abs(ITENSOR) is the flag, but the <span class="math notranslate nohighlight">\(\nabla \cdot D\)</span> term is not included in the computation of particle displacements. Under normal circumstances, an approximation of the term <span class="math notranslate nohighlight">\(\nabla \cdot D\)</span> is used in the particle tracking algorithm to obtain accurate solutions in cases where there are gradients in <span class="math notranslate nohighlight">\(D\)</span>.</td>
</tr>
<tr class="row-odd"><td>IREVERS</td>
<td>integer</td>
<td>Flag indicating if reverse particle tracking should be performed. If omitted, forward tracking is performed.IREVERS = 0,  Standard forward trackingIREVERS =  -1, Forward tracking only after exiting the time loop (this is needed for comparing results with reverse tracking)IREVERS = +1, Reverse tracking.Note: When using reverse tracking, turn off the dispersion, ITENSOR = 0, as it does not make sense to try to reverse the random part of the displacement. The value for ITENSOR must be entered to use this option.</td>
</tr>
<tr class="row-even"><td>FREEZ_TIME</td>
<td>real</td>
<td>If greater than zero, time (days) at which flow solution is frozen and only particle transport is computed after that time. If omitted, the flow solution continues for the entire simulation. Values for ITENSOR and IREVERS must be entered to use this option.</td>
</tr>
<tr class="row-odd"><td>MAX_JUMP</td>
<td>integer</td>
<td>When using random walk, the maximum number of cells a particle is allowed to jump in a single step. (Default is 10).</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “tcurve” indicating transfer function curve data should be input to model matrix diffusion. If the keyword is found then NUMPARAMS and FILENAME are entered, otherwise they are omitted.</td>
</tr>
<tr class="row-odd"><td>NUMPARAMS</td>
<td>integer</td>
<td>Number of parameters that define the transfer function curves being used.</td>
</tr>
<tr class="row-even"><td>TFILENAME</td>
<td>character</td>
<td>Name of input file containing the transfer function curve data.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “omr” to indicate the grid has octree mesh refinement. If the keyword is found then OMR_NODES is entered, and optionally keyword “file” with OMRFILENAME and OMR_FORM, otherwise they are omitted.</td>
</tr>
<tr class="row-even"><td>OMR_NODES</td>
<td>integer</td>
<td>Number of refined (omr) nodes in the grid.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “file” indicating that the omr initialization calculation arrays should be written to or read from a file. This option should only be used with steday-state flow.</td>
</tr>
<tr class="row-even"><td>OMRFILENAME</td>
<td>character</td>
<td>Name of file from which to read or to write omr arrays.</td>
</tr>
<tr class="row-odd"><td>OMR_FORM</td>
<td>character</td>
<td>Format of the omr file, formatted or unformatted.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “tpor” to indicate tracer porosities should be read.</td>
</tr>
<tr class="row-odd"><td>PS_TRAC</td>
<td>real</td>
<td>Tracer porosity</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “file” indicating that the tracer porosities should be read from a file.</td>
</tr>
<tr class="row-odd"><td>TPORFILENAME</td>
<td>character</td>
<td>Name of file from which to read tracer porosity.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “wtdt”</td>
</tr>
<tr class="row-odd"><td>DELTAWT</td>
<td>real</td>
<td>Distance below the water table that particles should be started.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “volum” to indicate control volumes should be output. Output is written to the “.sptrx” file.</td>
</tr>
<tr class="row-odd"><td>SPTRX_FORMAT</td>
<td>character</td>
<td>File format for control volume output file “formatted” or “unformatted”. Default is formatted.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Optional keywords “po” (porosity), “sa” (saturation), “pe” (permeability), “de” (density), “pr” (pressure), “te” (temperature), “zo” (zone number), “id” (particle identifier) 1 per line, indicating which parameters will be output along the particle path (written to “.sptr2” file). If no keywords are present no parameter ddata will be output. Note that in older versions of FEHM porosity and saturation were output if no keywords were entered.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Keyword ‘tprp’ specifying that transport properties are to follow on subsequent lines.</td>
</tr>
<tr class="row-even"><td>TPRP_FLAG</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag indicating what type of transport property information is to follow on the line</div>
<div class="line">1 - KD only</div>
<div class="line">2 -  KD and 5 terms of dispersivity tensor</div>
<div class="line">3 - (Dual porosity) - Matrix KD, diffusion coefficient, retardation factor in fracture, and fracture aperture. No dispersion</div>
<div class="line">4 - (Dual porosity) - Matrix KD, diffusion coefficient, retardation factor in fracture, fracture aperture, and 5 terms of dispersivity tensor</div>
<div class="line">11 -  Colloid diversity model with importance sampling, CDF vs Retardation Factor specified as a table in the optional file specified by CDF_FILENAME</div>
<div class="line">12 -  Similar to case 11 above , except the SQRT(CDF) is used instead of CDF for importance sampling</div>
<div class="line">13 - Colloid diversity model with importance sampling, CDF vs <span class="math notranslate nohighlight">\(K_f\)</span> (attachment rate constant) specified as a straight line equation in the log-log space either on this line or in the optional file specified by CDF_FILENAME</div>
<div class="line">14 -  Similar to case 13 above, except the SQRT(CDF) is used instead of CDF for importance sampling</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>SIMNUM</td>
<td>integer</td>
<td>Simulation number, used for selecting the table/equation from the colloid diversity file.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating the cumulative probability distribution function (CDF) retardation parameters for the colloid diversity model should be read from an external file</td>
</tr>
<tr class="row-odd"><td>CDF_FILENAME</td>
<td>character*80</td>
<td><div class="first last line-block">
<div class="line">Name of the file containing the cumulative probability distribution function (CDF) (entered if optional keyword ‘file’ follows keyword ‘dive’).</div>
<div class="line">If TPRPFLAG = 11 or 12, Table option.</div>
<div class="line">If TPRPFLAG = 13 or 14, Equation option.</div>
<div class="line">The following equations are used for <span class="math notranslate nohighlight">\(R_{min} \le R \le R_{max}\)</span>, <span class="math notranslate nohighlight">\(R = 1 + K_f / K_{rev}\)</span>, <span class="math notranslate nohighlight">\(\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)\)</span></div>
</div>
</td>
</tr>
<tr class="row-even"><td>KD</td>
<td>real</td>
<td>Matrix sorption coefficient</td>
</tr>
<tr class="row-odd"><td>DIFM</td>
<td>real</td>
<td>Diffusion coefficient applying to matrix diffusion submodel (m2/s)</td>
</tr>
<tr class="row-even"><td>RD_FRAC</td>
<td>real</td>
<td>Retardation factor in fracture media</td>
</tr>
<tr class="row-odd"><td>POR_MATRIX</td>
<td>real</td>
<td>Matrix porosity (fracture volume fraction is specified in rock macro)</td>
</tr>
<tr class="row-even"><td>APERTURE</td>
<td>real</td>
<td>Fracture aperture (m)</td>
</tr>
<tr class="row-odd"><td>AL</td>
<td>real</td>
<td>Longitudinal dispersivity, αL (m). ITENSOR = 2, 4, or 5</td>
</tr>
<tr class="row-even"><td>ALH</td>
<td>real</td>
<td>Horizontal longitudinal dispersivity, αLH (m). ITENSOR = 3</td>
</tr>
<tr class="row-odd"><td>ALV</td>
<td>real</td>
<td>Vertical longitudinal dispersivity, αLT (m). ITENSOR = 3</td>
</tr>
<tr class="row-even"><td>AT</td>
<td>real</td>
<td>Transverse dispersivity, αT (m). ITENSOR = 4</td>
</tr>
<tr class="row-odd"><td>ATH</td>
<td>real</td>
<td>Transverse horizontal dispersivity, αTH (m). ITENSOR = 2, 3, or 5</td>
</tr>
<tr class="row-even"><td>ATV</td>
<td>real</td>
<td>Transverse vertical dispersivity, αTV (m). ITENSOR = 2, 3, or 5</td>
</tr>
<tr class="row-odd"><td>A1</td>
<td>real</td>
<td>Generalized dispersivity term α1 (m) from Lichtner et al. (2002)</td>
</tr>
<tr class="row-even"><td>A2</td>
<td>real</td>
<td>Generalized dispersivity term α2 (m) from Lichtner et al. (2002)</td>
</tr>
<tr class="row-odd"><td>A3</td>
<td>real</td>
<td>Generalized dispersivity term α3 (m) from Lichtner et al. (2002)</td>
</tr>
<tr class="row-even"><td>A4</td>
<td>real</td>
<td>Generalized dispersivity term α4 (m) from Lichtner et al. (2002)</td>
</tr>
<tr class="row-odd"><td>ASX</td>
<td>real</td>
<td>Direction cosine of the axis of symmetry from Lichtner et al. (2002)</td>
</tr>
<tr class="row-even"><td>ASY</td>
<td>real</td>
<td>Direction cosine of the axis of symmetry from Lichtner et al. (2002)</td>
</tr>
<tr class="row-odd"><td>DM</td>
<td>real</td>
<td>Molecular diffusion coefficient (m2/s)</td>
</tr>
<tr class="row-even"><td>VRATIO</td>
<td>real</td>
<td>Parameter to control the movement of particles into low velocity cells via random walk. Used to restrict the artificial migration of particles into low permeability zones due to dispersion. The value of VRATIO is used as a ratio for determining if random walk into a new cell is allowed. If the ratio of the average velocity in the new cell divided by the velocity in the previous cell is less than VRATIO, then the particle is not allowed to migrate into the new cell. It is returned to its previous location, and a new random walk is computed and applied. Up to 10 attempts at a random walk are allowed, after which the particle location is left at the current location for the next advective step.</td>
</tr>
<tr class="row-odd"><td>MODEL_NUMBER</td>
<td>integer</td>
<td>Number of model (referring to the sequence of models read) to be assigned to the designated nodes or zone.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>keyword</td>
<td>Optional keyword ‘zbtc’ specifying that zone breakthrough curves will be computed. Output will be written to the “.sptr3” file. If ‘zbtc’ is omitted, so are NZBTC and ZBTC. Note that the zones must be specified in a zone macro preceding the sptr macro in the input file before they are invoked using the keyword ‘zbtc’.</td>
</tr>
<tr class="row-odd"><td>NZBTC</td>
<td>integer</td>
<td>Number of zones for which breakthrough curves will be computed.</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>Note that DTMN, PART_MULT, PART_FRAC, DELTA_PART, PART_STEPS, and TIME_BTC are optional input. They must be entered in the order given. When not entered the default values will be used.</td>
</tr>
<tr class="row-odd"><td>DTMN</td>
<td>real</td>
<td>Time step control (seconds). FEHM will account for all particles at time step intervals starting with dtmn seconds and write information to the “.sptr3” file if the “zbtc” keyword is present. This controls the output density for breakthrough curve information only. (Default DTMX)</td>
</tr>
<tr class="row-even"><td>PART_MULT</td>
<td>real</td>
<td>Time step multiplication factor. (Default 2.)</td>
</tr>
<tr class="row-odd"><td>PART_FRAC</td>
<td>real</td>
<td>Fraction of particles that should break through before checking if time step should be increased. (Default 0.1*NUMPART)</td>
</tr>
<tr class="row-even"><td>DELTA_PART</td>
<td>real</td>
<td>Fraction of particles that should break through during a time step so that the time step is not increased</td>
</tr>
<tr class="row-odd"><td>PART_STEPS</td>
<td>integer</td>
<td>Number of time steps that should be checked for DELTA_PART before increasing the time step</td>
</tr>
<tr class="row-even"><td>TIME_BTC</td>
<td>real</td>
<td>Time to start using small breakthrough time steps (DTMN) for late initial breakthrough (days).</td>
</tr>
<tr class="row-odd"><td>ZBTC</td>
<td>integer</td>
<td>NZBTC zone numbers of the zone(s) for which breakthrough curves will be computed.</td>
</tr>
<tr class="row-even"><td>ITM</td>
<td>integer</td>
<td>Maximum number of time steps to accomplish the FEHM time step ‘day’</td>
</tr>
<tr class="row-odd"><td>IST</td>
<td>integer</td>
<td>Flag to specify type of input for particlesIST = 0, local position and corresponding element number (Group 10)IST = 1, global position (Group 10)IST = 2, specify a zone of particles (Groups 7-9)</td>
</tr>
<tr class="row-even"><td>COUNT_STEPS_MAX</td>
<td>integer</td>
<td>Maximim number of steps a particle is allowed to take in a sptr run. (Default 1000000) Input of this value is optional. If this value is omitted, the default will be used. The value must precede SPTR_FLAG if being used.</td>
</tr>
<tr class="row-odd"><td>SPTR_FLAG</td>
<td>character</td>
<td>Optional keyword “save” to signal that final particle locations and times should be written to a file, <a href="#id5"><span class="problematic" id="id6">*</span></a>.sptrs, for a particle restart run.</td>
</tr>
<tr class="row-even"><td>NX</td>
<td>integer</td>
<td>Number of divisions in the x-direction</td>
</tr>
<tr class="row-odd"><td>NY</td>
<td>integer</td>
<td>Number of divisions in the y-direction</td>
</tr>
<tr class="row-even"><td>NZ</td>
<td>integer</td>
<td>Number of divisions in the z-direction</td>
</tr>
<tr class="row-odd"><td>X10</td>
<td>real</td>
<td>X-coordinate of the origin (xmin)</td>
</tr>
<tr class="row-even"><td>Y10</td>
<td>real</td>
<td>Y-coordinate of the origin (ymin)</td>
</tr>
<tr class="row-odd"><td>Z10</td>
<td>real</td>
<td>Z-coordinate of the origin (zmin)</td>
</tr>
<tr class="row-even"><td>XDIM</td>
<td>real</td>
<td>Length of X-direction</td>
</tr>
<tr class="row-odd"><td>YDIM</td>
<td>real</td>
<td>Length of Y-direction</td>
</tr>
<tr class="row-even"><td>ZDIM</td>
<td>real</td>
<td>Length of Z-direction</td>
</tr>
<tr class="row-odd"><td>IJKV(I)</td>
<td>integer</td>
<td>Node or element number</td>
</tr>
<tr class="row-even"><td>X1(I)</td>
<td>real</td>
<td>Starting X-coordinate for a particle</td>
</tr>
<tr class="row-odd"><td>Y1(I)</td>
<td>real</td>
<td>Starting Y-coordinate for a particle</td>
</tr>
<tr class="row-even"><td>Z1(I)</td>
<td>real</td>
<td>Starting Z-coordinate for a particle</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Optional keyword “file” indicating that the particle starting locations should be read from a file. If this file has been generated by the code using the “save” keyword in Group 6, particle starting times will also be read.</td>
</tr>
<tr class="row-even"><td>SPTRFILENAME</td>
<td>character</td>
<td>Name of file from which to read initial particle locations.</td>
</tr>
</tbody>
</table>
<p>The following are examples of sptr. In the first example 10000 particles are
inserted at the inlet within a single cell, and the breakthrough curve at a
downstream location (defined in a call to zone) is recorded for the case of
longitudinal dispersion with a dispersivity of 100 m and sorption with a KD of
0.0223715. Breakthrough concentration is output every 1.728e8 seconds to the
“.sptr3” file.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">sptr</span>
    <span class="mf">1.728e8</span>   <span class="mi">0</span>    <span class="mi">0</span>   <span class="mi">0</span>                    <span class="n">Group</span> <span class="mi">1</span>
    <span class="mf">0.25</span>      <span class="mi">0</span>    <span class="mi">5</span>   <span class="mi">0</span>                    <span class="n">Group</span> <span class="mi">2</span>
<span class="n">tprp</span>
    <span class="mi">2</span>         <span class="mf">0.0223715</span> <span class="mf">100.</span> <span class="mf">0.</span> <span class="mf">0.</span> <span class="mf">0.</span> <span class="mf">0.</span>    <span class="n">Group</span> <span class="mi">3</span>

    <span class="mi">1</span>         <span class="mi">0</span>         <span class="mi">0</span>    <span class="mi">1</span>              <span class="n">Group</span> <span class="mi">4</span>

<span class="n">zbtc</span>
    <span class="mi">1</span>
    <span class="mi">5</span>
 <span class="mi">1000</span>         <span class="mi">2</span>                             <span class="n">Group</span> <span class="mi">5</span>
    <span class="mi">1</span>       <span class="mi">100</span>      <span class="mi">1000</span>                   <span class="n">Group</span> <span class="mi">6</span>
    <span class="mf">0.</span>    <span class="o">-</span><span class="mf">1500.</span>        <span class="mf">0.</span>                  <span class="n">Group</span> <span class="mi">7</span>
   <span class="mf">10.</span>     <span class="mf">3000.</span>       <span class="mf">12.5</span>                 <span class="n">Group</span> <span class="mi">8</span>
                                            <span class="n">Group</span> <span class="mi">9</span>
</pre></div>
</div>
<p>In the second example, both longitudinal and transverse dispersion are invoked,
but no sorption. The solute is input as a patch on the inlet face of the model.
The dimensions of the patch will be 3,000 m in the y-direction and 12.5 m in the
vertical direction, starting at the surface, and 100000 particles are injected.
Data to generate a steady state concentration plume is output in the “.trc” file.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">sptr</span>
    <span class="mf">2.88e7</span>      <span class="mi">0</span>           <span class="mi">0</span>       <span class="mi">0</span>                                   <span class="n">Group</span> <span class="mi">1</span>
    <span class="mf">0.25</span>        <span class="mi">0</span>           <span class="mi">5</span>       <span class="mi">0</span>                                   <span class="n">Group</span> <span class="mi">2</span>

<span class="n">tprp</span>
    <span class="mi">2</span>           <span class="mf">0.</span>       <span class="mf">100.</span>     <span class="mf">0.1</span>       <span class="mf">0.1</span>     <span class="mf">0.</span>      <span class="o">-</span><span class="mf">1.e-10</span>     <span class="n">Group</span> <span class="mi">3</span>

    <span class="mi">1</span>           <span class="mi">0</span>          <span class="mi">0</span>        <span class="mi">1</span>                                   <span class="n">Group</span> <span class="mi">4</span>

    <span class="mi">1000</span>        <span class="mi">2</span>                                                       <span class="n">Group</span> <span class="mi">5</span>
       <span class="mi">1</span>      <span class="mi">100</span>       <span class="mi">1000</span>                                            <span class="n">Group</span> <span class="mi">6</span>
       <span class="mf">0.</span>   <span class="o">-</span><span class="mf">1500.</span>         <span class="mf">0.</span>                                           <span class="n">Group</span> <span class="mi">7</span>
      <span class="mf">10.</span>    <span class="mf">3000.</span>       <span class="o">-</span><span class="mf">12.5</span>                                          <span class="n">Group</span> <span class="mi">8</span>
</pre></div>
</div>
<p>The third example uses the colloid diversity model with importance sampling specified
as an equation using an external input file, using the third set of parameters in
the rcoll_eqn.dat, with the file rcoll_eqn.dat as:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">Colloid</span> <span class="n">diversity</span> <span class="n">model</span> <span class="n">equation</span> <span class="n">parameters</span>
<span class="mi">1</span> <span class="mf">1.5641426E-5</span> <span class="mf">1.0</span> <span class="mf">63933.785</span> <span class="mf">0.7081742</span> <span class="mf">0.0E+0</span> <span class="mf">100.</span> <span class="mf">10.</span> <span class="mf">0.1</span> <span class="mf">5.e-12</span> <span class="mf">0.1</span>
<span class="mi">2</span> <span class="mf">1.1755084E-3</span> <span class="mf">1.0</span> <span class="mf">851.69573</span> <span class="mf">0.7676392</span> <span class="mf">0.0E+0</span> <span class="mf">100.</span> <span class="mf">10.</span> <span class="mf">0.1</span> <span class="mf">5.e-12</span> <span class="mf">0.1</span>
<span class="mi">3</span> <span class="mf">1.0417102E-5</span> <span class="mf">1.0</span> <span class="mf">95996.984</span> <span class="mf">0.7438557</span> <span class="mf">0.0E+0</span> <span class="mf">100.</span> <span class="mf">10.</span> <span class="mf">0.1</span> <span class="mf">5.e-12</span> <span class="mf">0.1</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="o">.</span>
<span class="mi">100</span> <span class="mf">2.0808208E-4</span> <span class="mf">1.0</span> <span class="mf">4806.7954</span> <span class="mf">0.62846046</span> <span class="mf">0.0E+0</span> <span class="mf">100.</span> <span class="mf">10.</span> <span class="mf">0.1</span> <span class="mf">5.e-12</span> <span class="mf">0.1</span>
</pre></div>
</div>
<p>The fourth example uses the colloid diversity model with importance sampling
specified as an equation with parameters input in the sptr macro.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">sptr</span>
<span class="mf">1.728e8</span>         <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">0</span>       <span class="n">Group</span> <span class="mi">1</span>
<span class="mf">0.25</span>            <span class="mi">0</span>       <span class="mi">5</span>       <span class="mi">0</span>       <span class="n">Group</span> <span class="mi">2</span>
<span class="n">tprp</span>
<span class="mi">13</span>              <span class="mi">3</span>                       <span class="n">Group</span> <span class="mi">3</span>
<span class="n">rcoll_eqn</span><span class="o">.</span><span class="n">dat</span>

<span class="mi">1</span>               <span class="mi">0</span>       <span class="mi">0</span>       <span class="mi">1</span>       <span class="n">Group</span> <span class="mi">4</span>

<span class="n">zbtc</span>
<span class="mi">1</span>
<span class="mi">5</span>
<span class="mi">1000</span>            <span class="mi">2</span>                       <span class="n">Group</span> <span class="mi">5</span>
<span class="mi">1</span>               <span class="mi">100</span>     <span class="mi">1000</span>            <span class="n">Group</span> <span class="mi">6</span>
<span class="mf">0.</span>              <span class="o">-</span><span class="mf">1500.</span>  <span class="mf">0.</span>              <span class="n">Group</span> <span class="mi">7</span>
<span class="mf">10.</span>             <span class="mf">3000.</span>   <span class="o">-</span><span class="mf">12.5</span>           <span class="n">Group</span> <span class="mi">8</span>
                                        <span class="n">Group</span> <span class="mi">9</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
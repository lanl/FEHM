---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">trac</span></code><a class="headerlink" href="#trac" title="Permalink to this headline">¶</a></h1>
<ul class="simple">
<li>Group 1 - KEYWORD ‘userc’, ANO, AWC, EPC, UPWGTA<ul>
<li>Optional keyword “file” is used to specify the name of the data file that will contain input for the userc subroutine.</li>
<li>KEYWORD ‘file’</li>
<li>USERC_FILENAME</li>
</ul>
</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 1 - ANO, AWC, EPC, UPWGTA</li>
<li>Group 2 - DAYCS, DAYCF, DAYHF, DAYHS</li>
<li>Group 3 - IACCMX, DAYCM, DAYCMM, DAYCMX, NPRTTRC</li>
<li>Group 4 - KEYWORD ‘tpor’</li>
<li>Group 5 - JA, JB, JC, PS_TRAC</li>
</ul>
<p>Tracer porosity is entered only if the Group 4 keyword (‘tpor’), which specifies tracer porosity input, is present, otherwise Groups 4 and 5 are omitted.</p>
<ul class="simple">
<li>Group 6 - NSPECI</li>
<li>Group 7 - KEYWORD ‘ldsp’</li>
</ul>
<p>The Group 7 keyword (‘ldsp’) specifies longitudinal / transverse dispersion should be used. If X, Y, Z dispersion is desired Group 7 is omitted, and dispersivities are input in X, Y, Z order in Group 9 or Group 12. When longitudinal / transverse dispersion is invoked the Z-components of dispersivity are omitted from the Group 9 or Group 12 input, and X and Y represent longitudinal and transverse dispersion respectively. Note that an “L” or “V” added to the Group 9 or Group 12 variable names (MFLAG, SEHDIFF, TCX, TCY, TCZ, IADSF, A1ADSF, A2ADSF, BETADF, DIFFM) indicates the value is for the liquid or vapor phase, respectively.</p>
<ul class="simple">
<li>Group 8 - KEYWORD ‘dspl’ or ‘dspv’ or ‘dspb’</li>
</ul>
<p>The Group 8 keyword specifies that the same diffusion coefficient and dispersivities are to be used for all species of the same type (liquid and/or vapor). This will make the calculations more efficient and thus should be used if applicable. If Group 8 is omitted, Groups 9 and 10 are also omitted, and input resumes with Group 11.</p>
<p>If only liquid species are present (keyword ‘dspl’) or only vapor species are present (keyword ‘dspv’) with no longitudinal / transverse dispersion, Group 9 is defined as follows:</p>
<ul class="simple">
<li>Group 9 - MFLAG, SEHDIFF, TCX, TCY, TCZ</li>
</ul>
<p>Otherwise if both liquid and vapor are present (keyword ‘dspb’), parameters for both must be entered.</p>
<ul class="simple">
<li>Group 9 - MFLAGL, SEHDIFFL, TCLX, TCLY, TCLZ, MFLAGV, SEHDIFFV, TCVX, TCVY, TCVZ</li>
<li>Groups 9 is used to define transport models for which diffusion and dispersion parameters are identical. Group 9 is read in until a blank line is encountered. The model number is incremented by 1 each time a line is read.</li>
<li>Group 10 -  JA, JB, JC, ITRCDSP</li>
<li>Group 11 - ICNS [SPNAM]</li>
</ul>
<p>There are two options for group twelve. If the same diffusion coefficient and dispersivities are to be used for all species of the same type (liquid and/ or vapor - keyword ‘dspl’, ‘dspv’, or ‘dspb’) only sorption parameters are input:</p>
<ul class="simple">
<li>Group 12 - IADSF, A1ADSF, A2ADSF, BETADF</li>
</ul>
<p>or for a Henry’s Law Species (both liquid and vapor)</p>
<ul class="simple">
<li>Group 12 - IADSFL, A1ADSFL, A2ADSFL, BETADFL, IADSFV, A1ADSFV, A2ADSFV, BETADFV</li>
</ul>
<p>In the absence of a Group 8 keyword (“dspl’, ‘dspv’, or ‘dspb’) the following input (for liquid or vapor) which includes the sorption and dispersion parameters is used:</p>
<ul class="simple">
<li>Group 12 - IADSF, A1ADSF, A2ADSF, BETADF, MFLAG, DIFFM, TCX, TCY, TCZ</li>
</ul>
<p>For a Henry’s Law Species (both liquid and vapor) if DIFFML ≥ 0</p>
<ul class="simple">
<li>Group 12 - IADSFL, A1ADSFL, A2ADSFl, BETADFL, MFLAGL, DIFFML, TCLX, TCLY, TCLZ, IADSFV, A1ADSFV, A2ADSFV, BETADFV, MFLAGV, DIFFMV, TCVX, TCVY, TCVZ</li>
<li>Group 13 - JA, JB, JC, ITRCD</li>
<li>Group 14 - HENRY_MODEL, HAWWA(1), HAWWA(2), HAWWA(3), HAWWA(4), HAWWA(5) (only input for a Henry’s Law species, otherwise omitted)</li>
<li>Group 15 - JA, JB, JC, ANQO</li>
<li>Group 16 - JA, JB, JC, CNSK, T1SK, T2SK</li>
</ul>
<p>Groups 11, 12, 13, 14, 15, and 16 are entered as a unit for each solute. However, for a solid species, only groups 11, 15, and 16 are entered (groups 12, 13, and 14 are not applicable for a solid species). Groups 12 and 13 are used to define transport models for which sorption, diffusion and dispersion parameters are identical. For a liquid or vapor species, only one set of Group 12 parameters should be entered per region. However, for a Henry’s Law species, two sets of parameters per region must be entered. For this case, the liquid sorption parameters should be entered on the first line and the vapor sorption parameters on a second line or as a continuation of the first line. Group 12 is read in until a blank line is encountered. The model number is incremented by 1 each time a line is read. Group 13 then assigns a transport model number to every node.</p>
<p>Injection nodes must be specified in control statement <strong>flow</strong>.</p>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="13%" />
<col width="72%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character*5</td>
<td>Keyword for invoking a solute transport user subroutine.
If the word <code class="docutils literal notranslate"><span class="pre">userc</span></code> is placed in this position, then the
code invokes a solute transport user subroutine at each time
step. Omit this key word if there is no solute user subroutine
for the simulation.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating the name for the user
subroutine input transport parameter file will be input.
If this keyword and the following line are omitted, the
default name will be <code class="docutils literal notranslate"><span class="pre">userc_data.dat</span></code>.</td>
</tr>
<tr class="row-even"><td>USERC_FILENAME</td>
<td>character*80</td>
<td>Name of file from which to read transport parameters for
optional user subroutine.</td>
</tr>
<tr class="row-odd"><td>ANO</td>
<td>real</td>
<td>Initial solute concentration, set at all nodes for all
species unless overwritten by a restart file input or
values in group 14 below (moles/kg fluid).</td>
</tr>
<tr class="row-even"><td>AWC</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Implicitness factor for solute solution.</div>
<div class="line">AWC &gt; 1.0 gives 2nd order solution</div>
<div class="line">AWC ≤ 1.0 gives 1st order solution</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>EPC</td>
<td>real</td>
<td>Equation tolerance for solute solution. When the square
root of the sum of the squared residuals is lower than
EPC, the solution is assumed to be converged.</td>
</tr>
<tr class="row-even"><td>UPWGTA</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Upstream weighting term for the solute solution.</div>
<div class="line">UPWGTA &lt; 0.5 UPWGTA is set to 0.5</div>
<div class="line">UPWGTA &gt; 1.0 UPWGTA is set to 1.0</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>DAYCS</td>
<td>real</td>
<td>Time which the solute solution is enabled (days).</td>
</tr>
<tr class="row-even"><td>DAYCF</td>
<td>real</td>
<td>Time which the solute solution is disabled (days).</td>
</tr>
<tr class="row-odd"><td>DAYHF</td>
<td>real</td>
<td>Time which the flow solution is disabled (days).</td>
</tr>
<tr class="row-even"><td>DAYHS</td>
<td>real</td>
<td>Time which the flow solution is enabled (days).</td>
</tr>
<tr class="row-odd"><td>IACCMX</td>
<td>integer</td>
<td>Maximum number of iterations allowed in solute
solution if time step multiplier is enabled</td>
</tr>
<tr class="row-even"><td>DAYCM</td>
<td>real</td>
<td>Time step multiplier for solute solution</td>
</tr>
<tr class="row-odd"><td>DAYCMM</td>
<td>real</td>
<td>Initial time step for solute solution (days)</td>
</tr>
<tr class="row-even"><td>DAYCMX</td>
<td>real</td>
<td>Maximum time step for solute solution (days)</td>
</tr>
<tr class="row-odd"><td>NPRTTRC</td>
<td>integer</td>
<td>Print-out interval for solute information. Data for
every NPRTTRC solute time step will be written to the
“.trc” file. If this parameter is omitted (for
compatibility with old input files) the default
value is 1. Note that the first and last solute
time step within a heat and mass transfer step
automatically get printed.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td>Keyword ‘tpor’ specifying optional tracer porosity
should be input. If group 4 is omitted, porosities
assigned in macro rock are used.</td>
</tr>
<tr class="row-odd"><td>PS_TRAC</td>
<td>real</td>
<td>Tracer porosity</td>
</tr>
<tr class="row-even"><td>NSPECI</td>
<td>integer</td>
<td>Number of solutes simulated.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Keyword ‘ldsp’ specifying longitudinal / transverse dispersion.
If x, y, z dispersion is desired group 7 is omitted, and
dispersivities are input in x, y, and then z order
(group 9 or group 12). Otherwise, if longitudinal
/ transverse dispersion is desired the keyword
‘ldsp’ is entered and dispersivities are instead
input in longitudinal and then transverse order
with values for the third dimension omitted.</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character*4</td>
<td><div class="first last line-block">
<div class="line">Keyword specifying the same diffusion coefficient and
dispersivities are to be used for all species of the same
type (liquid and/or vapor).</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">dspl</span></code> indicates that only liquid species exist.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">dspv</span></code> indicates that only vapor species exist.</div>
<div class="line"><code class="docutils literal notranslate"><span class="pre">dspb</span></code> indicates that both liquid and vapor species exist.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>ICNS</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Phase designation for the ith solute</div>
<div class="line">* -2 - Henry’s Law species (input and output concentration
values are gas concentrations).</div>
<div class="line">* -1 - Vapor species.</div>
<div class="line">* 0 - Solid species</div>
<div class="line">* 1 - Liquid species</div>
<div class="line">* 2 - Henry’s Law species (input and output concentration
values are liquid concentrations)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>SPNAM</td>
<td>character*20</td>
<td><dl class="first last docutils">
<dt>For each species, the name of the species (e.g. Sulfate).</dt>
<dd>This is an optional identifier that may be input when
macro <code class="docutils literal notranslate"><span class="pre">rxn</span></code> is not being used.</dd>
</dl>
</td>
</tr>
<tr class="row-odd"><td>MFLAG</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag denoting type of diffusion model to be used</div>
<div class="line">0 - the molecular diffusion coefficient is a constant.</div>
<div class="line">1 - Millington Quirk diffusion model for liquid or vapor.</div>
<div class="line">2 - Conca and Wright diffusion model for liquid, alternate
Millington Quirk diffusion model for vapor.</div>
<div class="line">3 - vapor diffusion coefficient is calculated as a function
of pressure and temperature using tortuosity from adif
macro, of the “Models and Methods Summary” of the FEHM Application
(Zyvoloski et al. 1999).</div>
<div class="line"><br /></div>
<div class="line">FEHM calculates liquid contaminant flux as
J = (Water Content)x(D*)x(GradC) and vapor contaminant flux as
J = (Air Content)x(D*)x(GradC) where D* is the diffusion
coefficient input in this macro. Water content is defined as
porosity x saturation and air content is defined as
porosity x (1 - saturation). For more explanation on the Millington
Quirk and Conca/Wright models see Stauffer, PH, JA Vrugt, HJ Turin,
CW Gable, and WE Soll (2009) <strong>Untangling diffusion from advection
in unsaturated porous media: Experimental data, modeling, and
parameter uncertainty assessment</strong>. Vadose Zone Journal,
8:510-522, doi:10.2136/vzj2008.0055.</div>
<div class="line"><br /></div>
<div class="line">SEHDIFF real Molecular diffusion coefficient (m2/s)</div>
<div class="line"><br /></div>
<div class="line">When MFLAG = 0, the input diffusion coefficient is used
directly in the contaminant flux equations presented above.
However, MFLAG = 1 or 2, the free air or free water diffusion
coefficient is input and the correct porous diffusion is
calculated within FEHM. For MFLAG = 3, the code assumes a
free air diffusion coefficient of 2.33e-5 m2/s for water
vapor in air as described in the Models and Methods
Summary Eq. 21.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>SEHDIFF</td>
<td>real</td>
<td>Molecular diffusion coefficient (m^2^/s)</td>
</tr>
<tr class="row-odd"><td>TCX</td>
<td>real</td>
<td>Dispersivity in x-direction (m)</td>
</tr>
<tr class="row-even"><td>TCY</td>
<td>real</td>
<td>Dispersivity in y-direction (m)</td>
</tr>
<tr class="row-odd"><td>TCZ</td>
<td>real</td>
<td>Dispersivity in z-direction (m)</td>
</tr>
<tr class="row-even"><td>ITRCDSP</td>
<td>integer</td>
<td>Region number for dispersion parameters given in group 9
(keyword dspl, dspv, or dspv). Default is 1.</td>
</tr>
<tr class="row-odd"><td>IADSF</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Adsorption model type for the ith species, ith region</div>
<div class="line">0 - conservative solute</div>
<div class="line">1 - linear sorption isotherm</div>
<div class="line">2 - Freundlich sorption isotherm</div>
<div class="line">3 - Modified Freundlich sorption isotherm</div>
<div class="line">4 - Langmuir sorption isotherm</div>
</div>
</td>
</tr>
<tr class="row-even"><td>A1ADSF</td>
<td>real</td>
<td>α1 parameter in adsorption model</td>
</tr>
<tr class="row-odd"><td>A2ADSF</td>
<td>real</td>
<td>α2 parameter in adsorption model</td>
</tr>
<tr class="row-even"><td>BETADF</td>
<td>real</td>
<td>β parameter in adsorption model</td>
</tr>
<tr class="row-odd"><td>DIFFM</td>
<td>real</td>
<td>Molecular diffusion coefficient (m^2^/s) See discussion
for SEHDIFF.</td>
</tr>
<tr class="row-even"><td>ITRCD</td>
<td>integer</td>
<td>Region number for group 12 sorption parameters or
for sorption and dispersion parameters (no keyword).
Default is 1.</td>
</tr>
<tr class="row-odd"><td>HENRY_MODEL</td>
<td>integer</td>
<td><div class="first last line-block">
<div class="line">Flag denoting which model is to be used for
defining the temperature dependence of the Henry’s law constant</div>
<div class="line">1 - van’t Hoff model</div>
<div class="line">2 - Multi-parameter fit to experimental data (used for carbonate system)</div>
<div class="line">3 - Henry’s model uses water vapor pressure (<span class="math notranslate nohighlight">\(H = P_{wv}\)</span>)</div>
</div>
</td>
</tr>
<tr class="row-even"><td>HAWWA(1)</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Term in Henry’s Law temperature dependence model:</div>
<div class="line">For model 1 or 3 - parameter value is <span class="math notranslate nohighlight">\(A_H\)</span></div>
<div class="line">For model 2 - parameter value is <span class="math notranslate nohighlight">\(A_{H,1}\)</span></div>
<div class="line">For model 3 - not used</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>HAWWA(2)</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Term in Henry’s Law temperature dependence model:</div>
<div class="line">For model 1 - parameter value is <span class="math notranslate nohighlight">\(\Delta H_H\)</span></div>
<div class="line">For model 2 - parameter value is <span class="math notranslate nohighlight">\(A_{H,2}\)</span></div>
<div class="line">For model 3 - Henry’s constant modifier,
<span class="math notranslate nohighlight">\(H = P_{wv} \cdot \Delta H_H\)</span></div>
</div>
</td>
</tr>
<tr class="row-even"><td>HAWWA(3)</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Term in Henry’s Law temperature dependence model:</div>
<div class="line">For model 1 - not used</div>
<div class="line">For model 2 - parameter value is <span class="math notranslate nohighlight">\(A_{H,3}\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>HAWWA(4)</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Term in Henry’s Law temperature dependence model:</div>
<div class="line">For model 1 - not used</div>
<div class="line">For model 2 - parameter value is <span class="math notranslate nohighlight">\(A_{H,4}\)</span></div>
</div>
</td>
</tr>
<tr class="row-even"><td>HAWWA(5)</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Term in Henry’s Law temperature dependence model:</div>
<div class="line">For model 1 - not used</div>
<div class="line">For model 2 - parameter value is <span class="math notranslate nohighlight">\(A_{H,5}\)</span></div>
</div>
</td>
</tr>
<tr class="row-odd"><td>ANQO</td>
<td>real</td>
<td>Initial concentration of tracer, which will supersede
the value given in group 1. Note that if initial
values are read from a restart file, these values
will be overwritten. Units are moles per kg vapor or
liquid for a liquid, vapor, or Henry’s law species,
and moles per kg of solid for a solid species.
Default is 0.</td>
</tr>
<tr class="row-even"><td>CNSK</td>
<td>real</td>
<td>Injection concentration at inlet node (moles per kg
liquid or vapor). If fluid is exiting at a node, then
the in-place concentration is used. If CNSK &lt; 0, then
the concentration at that particular node will be held
at a concentration of abs(cnsk) (default is 0 for all
unassigned nodes).</td>
</tr>
<tr class="row-odd"><td>T1SK</td>
<td>real</td>
<td>Time (days) when tracer injection begins. Default is 0.</td>
</tr>
<tr class="row-even"><td>T2SK</td>
<td>real</td>
<td>Time (days) when tracer injection ends. Default is 0.
If T2SK &lt; 0, the absolute value of T2SK is used for this
parameter, and the code interprets the negative value as a
flag to treat the node as a zero-solute-flux node for cases
in which a fluid sink is defined for that node. For this
case, the solute will stay within the model at the node
despite the removal of fluid at that location.
If a fluid source is present at the node, CNSK is
the concentration entering with that fluid, as in the
normal implementation of a solute source. Note that the
code cannot handle the case of T2SK &lt; 0 and CNSK &lt; 0
(fixed concentration), as these are incompatible inputs.
Therefore, the code prints an error message and stops
for this condition.</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
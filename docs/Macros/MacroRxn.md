---
layout : page_macros
hero_height: is-hidden
---


<h1><code class="docutils literal notranslate"><span class="pre">rxn</span></code><a class="headerlink" href="#rxn" title="Permalink to this headline">¶</a></h1>
<p>Chemical reactions between components are invoked with this control statement.
It is used in conjunction with control statement <code class="docutils literal notranslate"><span class="pre">trac</span></code>.
For facilitating the construction of the rxn input, a header describing the
input is required before each group whether data is entered for that group or
not (see examples) unless otherwise noted.</p>
<p>The header is an arbitrary text string that must be contained on a single line.
Note that for components that do not react with other components, <code class="docutils literal notranslate"><span class="pre">rxn</span></code> is unnecessary.
Specifically, conservative tracers or tracers that follow the equilibrium sorption
isotherms can be modeled with just the <code class="docutils literal notranslate"><span class="pre">trac</span></code> macro.</p>
<p>Note the parameters NCPNT, NIMM, NVAP (used as indices for input) are determined
by the code using information input for the <code class="docutils literal notranslate"><span class="pre">trac</span></code> macro.
NCPNT is equal to the number of liquid components, NIMM is equal to the number
of immobile components, and NVAP is equal to the number of vapor components
specified in the <code class="docutils literal notranslate"><span class="pre">trac</span></code> macro.</p>
<ul class="simple">
<li>Group 1 - NCPLX, NUMRXN</li>
<li>Group 2 - NGROUPS<ul>
<li>GROUP (ICPNT), ICPNT = 1, NCPNT (repeated NGROUPS times, once for each group).</li>
</ul>
</li>
<li>Group 3 - IDCPNT, CPNTNAM, IFXCONC, CPNTPRT, CPNTGS (repeated NCPNT times)</li>
<li>Group 4 - IDCPLX, CPLXNAM, CPLXPRT (repeated NCPLX times)</li>
<li>Group 5 - IDIMM, IMMNAM, IMMPRT (repeated NIMM times;)</li>
<li>Group 6 - IDVAP, VAPNAM, VAPPRT (repeated NVAP times)</li>
<li>Group 7 - ISKIP</li>
<li>Group 8 - RSDMAX<ul>
<li>HEADING</li>
</ul>
</li>
</ul>
<p>Note that this is an additional heading line that precedes the LOGKEQ heading.</p>
<ul class="simple">
<li>Group 9 - LOGKEQ</li>
<li>Group 10 - CKEQ, HEQ (NCPLX times)</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 10 - KEYWORD, NUM_TEMPS<ul>
<li>EQTEMP(I), I = 1, NUM_TEMPS</li>
<li>LKEQ(I), I = 1, NUM_TEMPS</li>
</ul>
</li>
</ul>
<p>For group 10, the keyword ‘lookup’ can be used in place of CKEQ and HEQ. LOOKUP allows a lookup table to be used to describe the equilibrium constant (K) as a function of temperature. After the keyword, the user must specify the number of values of temperature and K that will be used to describe K as a function of temperature. On the next lines the temperatures, and then the K values are entered. FEHM performs a piecewise linear interpolation between the values given.</p>
<ul class="simple">
<li>Group 11 - STOIC(ICPNT), ICPNT = 1, NCPNT (repeated NCPLX times, once for each aqueous complex)</li>
</ul>
<p>Input for groups 9, 10 and 11 is omitted if NCPLX, the number of aqueous complexes, is zero.</p>
<p>The remaining groups are entered as a unit for each kinetic reaction. If there are no kinetic reactions specified, none of the following groups (or their headers) are input. The input for groups 14 and on, depend on the kinetic reaction type specified in group 12. The input for each kinetic reaction type is described below.</p>
<ul class="simple">
<li>Group 12 - IDRXN</li>
<li>Group 13 - JA, JB, JC (JA, JB, JC - defined on page 22)</li>
<li><strong>IDRXN = 1: Linear kinetic reaction</strong></li>
<li>Group 14 - IAQUEOUS, IIMMOBILE</li>
<li>Group 15 - KD</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 15 - KEYWORD, NUM_TEMPS<ul>
<li>EQTEMP(I), I = 1, NUM_TEMPS</li>
<li>TCOEFF(I), I = 1, NUM_TEMPS</li>
</ul>
</li>
</ul>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">The distribution coefficient can be replaced with the keyword ‘lookup’ for temperature dependent coefficients. The keyword is described in Group 10.</p>
</div>
<ul class="simple">
<li>Group 16 - RATE</li>
<li><strong>IDRXN = 2: Langmuir kinetic reaction</strong></li>
<li>Group 14 - IAQUEOUS, IIMMOBILE</li>
<li>Group 15 - DISTCOEEF</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 15 - KEYWORD, NUM_TEMPS<ul>
<li>EQTEMP(I), I = 1, NUM_TEMPS</li>
<li>TCOEFF(I), I = 1, NUM_TEMPS</li>
</ul>
</li>
</ul>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">The distribution coefficient can be replaced with the keyword ‘lookup’ for temperature dependent coefficients. The keyword is described in Group 10.</p>
</div>
<ul class="simple">
<li>Group 16 - RATE</li>
<li>Group 17 - MAXCONC</li>
<li><strong>IDRXN = 3: General reaction</strong></li>
<li>Group 14 - NIMMOBILE, NAQUEOUS, NVAPOR</li>
<li>Group 15 - KFOR, KREV</li>
<li>Group 16 - IIMMOBILE (I = 1, NIMMOBILE)</li>
<li>Group 17 - IMSTOIC (I = 1, NIMMOBILE)</li>
</ul>
<p>Omit groups 16 and 17 (including headers) if NIMMOBILE is zero.</p>
<ul class="simple">
<li>Group 18 - IAQUEOUS (I = 1, NAQUEOUS)</li>
<li>Group 19 - AQSTOIC (I = 1, NAQUEOUS)</li>
</ul>
<p>Omit groups 18 and 19 (including headers) if NAQUEOUS is zero.</p>
<ul class="simple">
<li>Group 20 - IVAPOR (I = 1, NVAPOR)</li>
<li>Group 21 - IVSTOIC (I = 1, NVAPOR)</li>
</ul>
<p>Omit groups 20 and 21 (including headers) if NVAPOR is zero.</p>
<ul class="simple">
<li><strong>IDRXN = 4: Dual Monod kinetics biodegradation reaction</strong></li>
<li>Group 14 - NAQUEOUS (must be &gt;2 and &lt; 5), NIMMOBILE</li>
<li>Group 15 - SUBSTRATE, ELECACC, COMP3, COMP4, COMP5</li>
</ul>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">The users choice of NAQUEOUS determines whether COMP3, COMP4, and COMP5 need
to be entered. NIMMOBILE must be 1. NAQUEOUS also determines whether the
COMP3, COMP4, and COMP5 stoichiometries are to be used by the code.
However values must always be entered for groups 21-23 regardless of
the value of NAQUEOUS.</p>
</div>
<ul class="simple">
<li>Group 16 - BIOMASS</li>
<li>Group 17 - KS</li>
<li>Group 18 - KA</li>
<li>Group 19 - DECAY</li>
<li>Group 20 - ELECACCSTOIC</li>
<li>Group 21 - COMP3STOIC</li>
<li>Group 22 - COMP4STOIC</li>
<li>Group 23 - COMP5STOIC</li>
<li>Group 24 - PHTHRESH</li>
<li>Group 25 - QM</li>
<li>Group 26 - YIELD</li>
<li>Group 27 - XMINIT</li>
<li>Group 28 - NBIOFRM</li>
</ul>
<p>Omit Group 29 if NBIOFRM = 0</p>
<ul class="simple">
<li>Group 29 - ICBIO (I = 1, NBIOFRM)</li>
<li><strong>IDRXN=5: Radioactive Decay Reaction</strong></li>
<li>Group 14 - HALFLIFE</li>
<li>Group 15 - RXNTYPE</li>
<li>Group 16 - PARENT, DAUGHTER</li>
<li><strong>IDRXN=6: Kinetic Henry’s Law reaction</strong></li>
<li>Group 14 - IAQUEOUS, IVAPOR</li>
<li>Group 15 - KH</li>
<li>Group 16 - RATE</li>
<li><strong>IDRXN = 7: Kinetic precipitation/dissolution reaction for total component concentrations</strong></li>
<li>Group 14 - IIMMOBILE</li>
<li>Group 15 - NAQUEOUS</li>
<li>Group 16 - IAQUEOUS (I = 1, NAQUEOUS)</li>
<li>Group 17 - IMSTOIC</li>
<li>Group 18 - AQSTOIC (I = 1, NAQUEOUS)</li>
<li>Group 19 - SOLUBILITY</li>
</ul>
<p>or</p>
<ul class="simple">
<li>Group 19 - KEYWORD, NUM_TEMPS<ul>
<li>EQTEMP(I), I = 1, NUM_TEMPS</li>
<li>TCOEFF(I), I = 1, NUM_TEMPS</li>
</ul>
</li>
</ul>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">Solubility can be replaced with the keyword ‘lookup’ for temperature dependent solubilities. The keyword is described in Group 10.</p>
</div>
<ul class="simple">
<li>Group 20 - RATE</li>
<li>Group 21 - SAREA</li>
<li><strong>IDRXN = 8: Kinetic precipitation/dissolution reaction for total component concentrations with rates based on free-ion concentrations</strong></li>
</ul>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">The input is identical to that for reaction model 7.</p>
</div>
<table border="1" class="docutils">
<colgroup>
<col width="10%" />
<col width="1%" />
<col width="46%" />
<col width="43%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NCPLX</td>
<td>integer</td>
<td>Number of aqueous complexes (equal to number of equilibrium reactions)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>NUMRXN</td>
<td>integer</td>
<td>Number of kinetic reactions</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NGROUPS</td>
<td>integer</td>
<td>Number of groups. See GROUP to determine how to set this parameter.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>GROUP</td>
<td>integer</td>
<td>This variable controls the selective coupling solution method employed by FEHM. NCPNT values are entered for each line of input, and NGROUPS lines of input are required, one for each group. If a value is non-zero, then that aqueous component is present in the group. A value of zero denotes that the species is not present in the group. Grouping of aqueous components that take part in rapid kinetic reactions is required for convergence. However, memory requirements increase as the square of the maximum number of aqueous components in a group.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IDCPNT</td>
<td>integer</td>
<td>For each total aqueous component, the number identifying each total aqueous component (e.g. 1, 2, etc.)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>CPNTNAM</td>
<td>character*20</td>
<td>For each total aqueous component, the name of the total aqueous component (e.g. Sulfate)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IFXCONC</td>
<td>integer</td>
<td>For each total aqueous component, the Flag denoting the type of total aqueous component1 - total aqueous concentration is specified in TRAC macro2 - specify log of free ion concentration in TRAC macro (use for pH). For example, if H+ is the component, IFXCONC of 2 allows for pH to be directly input in place of concentration values in the TRAC macro.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>CPNTPRT</td>
<td>integer</td>
<td>For each total aqueous component, the Flag denoting which total aqueous component concentrations are printed to the “.trc” file and the AVS files.0 - Print to file1 - Do not print to file</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>CPNTGS</td>
<td>real</td>
<td>Guess for the initial uncomplexed component concentration used in speciation reactions. We recommend 1.0e-9. On rare occasions, the chemical speciation solver may have trouble converging. Choosing more representative values for CPNTGS will help convergence.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDCPLX</td>
<td>integer</td>
<td>For each aqueous complex, the number identifying each aqueous complex. By convention, the first complex should be given the number 101, the second, 102, etc.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>CPLXNAM</td>
<td>character*20</td>
<td>For each aqueous complex, the name of the aqueous complex (e.g. H2SO4)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>CPLXPRT</td>
<td>integer</td>
<td>For each aqueous complex, the Flag denoting which aqueous complex concentrations are printed to the “.trc” file and the AVS files.0 - Print to file1 - Do not print to file</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IDIMM</td>
<td>integer</td>
<td>For each immobile component, the number identifying each immobile component (e.g. 1, 2, etc.)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IMMNAM</td>
<td>character*20</td>
<td>For each immobile component, the name of the immobile component (e.g. Calcite, Co[adsorbed] )</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IMMPRT</td>
<td>integer</td>
<td>For each immobile component, the Flag denoting which immobile component concentrations are printed to the “.trc” file and the AVS files.0 - Print to file1 - Do not print to file</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDVAP</td>
<td>integer</td>
<td>Fore each vapor component, the number identifying each vapor component (e.g. 1, 2, etc.)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>VAPNAM</td>
<td>character*20</td>
<td>For each vapor component, the name of the vapor component (e.g. CO2[gas])</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>VAPPRT</td>
<td>integer</td>
<td>For each vapor component, the Flag denoting which vapor component concentrations are printed to the “.trc” file and the AVS files.0 - Print to file1 - Do not print to file</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>ISKIP</td>
<td>integer</td>
<td>Flag denoting whether chemical speciation calculations should be done at nodes which have already converged in a previous transport iteration0 -  Do chemical speciation calculations at each node for every iteration (recommended option)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>1 - To save computational time, this options tells FEHM to do equilibrium speciation calculations only at nodes which have not converged during the previous transport iteration. Sometimes, this option can lead to mass balance errors (mass balances can be checked in the “.out” file to see if results are satisfactory). This option is only recommended for very large problems.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>RSDMAX</td>
<td>real</td>
<td>The tolerance for the equilibrium speciation calculations. We recommend 1x10-9 for most problems.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>HEADING</td>
<td>character</td>
<td>One line descriptive comment which precedes the LOGKEQ heading</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>LOGKEQ</td>
<td>integer</td>
<td>Flag denoting the whether K or log K is entered by the user0 - constants are given as K1 - constants are given as log K</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>CKEQ</td>
<td>real</td>
<td>For each aqueous complex, the equilibrium constant</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>HEQ</td>
<td>real</td>
<td>For each aqueous complex, the enthalpy of the equilibrium reaction. The Van Hoff equation is used to determine the value of the equilibrium constant as a function of temperature. Note the keyword ‘lookup’ (see [[wiki:Macro14037</td>
<td>See For group 10, the keyword ‘lookup’ can be used in place of CKEQ and HEQ. LOOKUP allows a lookup table to be used to describe the equilibrium constant (K) as a function of temperature. After the keyword, the user must specify the number of values of temperature and K that will be used to describe K as a function of temperature. On the next lines the temperatures, and then the K values are entered. FEHM performs a piecewise linear interpolation between the values given.]]) can be used in the place of CKEQ and HEQ.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Keyword ‘lookup’ designating a lookup table will be used to describe the equilibrium constant as function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NUM_TEMPS</td>
<td>integer</td>
<td>Number of values of temperature and K that will be used to describe K as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>EQTEMP</td>
<td>real</td>
<td>Temperatures for K function (oC)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>LKEQ</td>
<td>real</td>
<td>Equilibrium constants for K function</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>STOIC</td>
<td>real</td>
<td>For each aqueous complex, the stoichiometry describing how to “make” the complex from the total aqueous components must be entered. If positive, the solute is a reactant; if negative, the solute is a product; and if 0, the solute is not present in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IDRXN</td>
<td>integer</td>
<td>For each kinetic reaction, this parameter specifies kinetic reaction model. Currently, the following reaction models are available. Additional kinetic formulations can be added without significant code development.1 - linear kinetic reaction2 - langmuir kinetic reaction3 - general kinetic reaction4 - dual Monod biodegradation reaction5 - radioactive decay reaction6 - kinetic Henry’s law reaction7 - precipitation/dissolution reaction8 - precipitation/dissolution reaction with rates based on free-ion concentration</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>JA, JB, JC</td>
<td>integer</td>
<td>JA, JB, JC are described on page 22. Here these parameters are used to specify the nodes at which the current kinetic reaction takes place. If the reaction takes place throughout the problem domain simply enter 1 0 0.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IDRXN = 1: Linear Kinetic Reaction</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IAQUEOUS</td>
<td>integer</td>
<td>The aqueous component number (e.g. 1, 2, etc.) or the aqueous complex number (e.g. 101, 102, etc.) which corresponds to the sorbing component</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IIMMOBILE</td>
<td>integer</td>
<td>The immobile component number (e.g. 1, 2, etc.) which corresponds to the sorbed component</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KD</td>
<td>real</td>
<td>Distribution coefficient (kg water / kg rock)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td>Keyword ‘lookup’ designating a lookup table will be used to describe the distribution coefficient as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>NUM_TEMPS</td>
<td>integer</td>
<td>Number of values (temperatures and distribution coefficients) that will be used to describe KD as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>EQTEMP</td>
<td>real</td>
<td>Temperatures for distribution coefficient function (oC)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>TCOEFF</td>
<td>real</td>
<td>Distribution coefficients corresponding to temperatures.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>RATE</td>
<td>real</td>
<td>Reaction rate parameter (1/hr)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDRXN = 2: Langmuir Kinetic Reaction</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IAQUEOUS</td>
<td>integer</td>
<td>The aqueous component number (e.g. 1, 2, etc.) or the aqueous complex number (e.g. 101, 102, etc.) which corresponds to the sorbing component</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IIMMOBILE</td>
<td>integer</td>
<td>The immobile component number (e.g. 1, 2, etc.) which corresponds to the sorbed component</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>DISTCOEFF</td>
<td>real</td>
<td>Distribution coefficient (kg water/ moles)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Keyword ‘lookup’ designating a lookup table will be used to describe the distribution coefficient as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NUM_TEMPS</td>
<td>integer</td>
<td>Number of values (temperatures and distribution coefficients) that will be used to describe DISTCOEFF as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>EQTEMP</td>
<td>real</td>
<td>Temperatures for distribution coefficient function (oC)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>TCOEFF</td>
<td>real</td>
<td>Distribution coefficients for DISTCOEFF as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>RATE</td>
<td>real</td>
<td>Reaction rate parameter (1/hr)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>MAXCONC</td>
<td>real</td>
<td>Maximum concentration (moles/kg rock)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDRXN = 3: General Kinetic Reaction</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NIMMOBILE</td>
<td>integer</td>
<td>The number of immobile components which participate in the reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>NAQUEOUS</td>
<td>integer</td>
<td>The number of aqueous components and complexes which participate in the reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NVAPOR</td>
<td>integer</td>
<td>The number of vapor species which participate in the reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KFOR</td>
<td>real</td>
<td>The forward reaction rate parameter. [(concentration units)p x s]-1, where p is the sum of the exponents of all concentrations in the forward reaction minus 1. Thus the units of the reaction rate are (concentration units)/hr.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>KREV</td>
<td>real</td>
<td>The reverse reaction rate parameter. [(concentration units)p x s]-1, where p is the sum of the exponents of all concentrations in the reverse reaction minus 1. Thus the units of the reaction rate are (concentration units)/hr.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IIMMOBILE</td>
<td>integer</td>
<td>The immobile component numbers which correspond to the immobile reactants and products in the reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IMSTOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to each immobile component participating in the reaction. If positive, the solute is a reactant; if negative, the solute is a product; and if 0, the solute is not present in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IAQUEOUS</td>
<td>integer</td>
<td>The aqueous component or aqueous complex numbers which correspond to the aqueous reactants and products in the reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>AQSTOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to each aqueous component or aqueous complex participating in the reaction. If positive, the solute is a reactant; if negative, the solute is a product; and if 0, the solute is not present in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IVAPOR</td>
<td>integer</td>
<td>The vapor component numbers which correspond to the vapor reactants and products in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IVSTOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to each vapor component participating in the reaction. If positive, the solute is a reactant; if negative, the solute is a product; and if 0, the solute is not present in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDRXN = 4: Dual Monod Biodegradation Reaction</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NAQUEOUS</td>
<td>integer</td>
<td>The number of aqueous species which participate in the reaction. At least 2 aqueous species must participate, the substrate (e.g. organic carbon) and the electron acceptor (e.g. Oxygen). Up to 5 aqueous species can participate. The third, fourth and fifth aqueous components are either reactants or products of the biodegradation reaction. The value entered for NAQUEOUS determines whether COMP3, COMP4, and COMP5 stoichiometries are to be used by the code.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>NIMMOBILE</td>
<td>integer</td>
<td>The number of immobile components which participate in the reaction. For the biodegradation reaction this value is 1.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>SUBSTRATE</td>
<td>integer</td>
<td>The aqueous component number which corresponds to the substrate (a.k.a the electron donor) for the biodegradation reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>ELECACC</td>
<td>integer</td>
<td>The aqueous component number which corresponds to the electron acceptor for the biodegradation reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>COMP3</td>
<td>integer</td>
<td>The aqueous component number which corresponds to a reactant or product in the biodegradation reaction (e.g. CO2, NH3, H+, etc.). Note that this parameter is optional. COMP3 should only be entered if NAQUEOUS&gt;2.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>COMP4</td>
<td>integer</td>
<td>The aqueous component number which corresponds to a reactant or product in the biodegradation reaction (e.g. CO2, NH3, H+, etc.).Note that this parameter is optional. COMP4 should only be entered if NAQUEOUS&gt;3.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>COMP5</td>
<td>integer</td>
<td>The aqueous component number which corresponds to a reactant or product in the biodegradation reaction (e.g. CO2, NH3, H+, etc.)Note that this parameter is optional. COMP5 should only be entered if NAQUEOUS&gt;4.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>BIOMASS</td>
<td>real</td>
<td>The solid component number which corresponds to the biomass (this is the immobile component).</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>KS</td>
<td>real</td>
<td>The Monod half-maximum-rate concentration for the substrate (moles/ kg water)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KA</td>
<td>real</td>
<td>The Monod half-maximum-rate concentration for the electron acceptor (moles/ kg water)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>DECAY</td>
<td>real</td>
<td>First order microbial decay coefficient (hr-1)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>ELECACCSTOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to the electron acceptor. Note that the stoichiometry of the substrate is 1 by definition.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>COMP3STOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to COMP3. A value is always entered whether or not it is used.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>COMP4STOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to COMP4. A value is always entered whether or not it is used.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>COMP5STOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to COMP5. A value is always entered whether or not it is used.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>PHTHRESH</td>
<td>real</td>
<td>In many systems, the biodegradation reaction will stop as the pH becomes either too acidic or basic. This parameter can be used to stop the biodegradation reaction once the simulated pH approaches a certain value. For example, if PHTHRESH = 10, the biodegradation reaction will cease if the pH is above 10 in the simulation. Note that PHTHRESH is an upper threshold for pH.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>QM</td>
<td>real</td>
<td>The maximum specific rate of substrate utilization (moles/kg biomass/hr)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>YIELD</td>
<td>real</td>
<td>The microbial yield coefficient (kg biomass/mole substrate)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>XMINIT</td>
<td>real</td>
<td>In many systems, biomass does not decay below a certain concentration. The biomass concentration is not allowed to fall below XMINIT (moles/ kg rock).</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>NBIOFRM</td>
<td>integer</td>
<td>Depending on the problem setup, many aqueous complexes can be formed from a total aqueous component. Only some of these complexes may be biodegradable. This parameter is used to specify which forms (the total aqueous concentration, the free ion concentration, or various complex concentrations) of the component degrade. If NBIOFRM = 0, the total aqueous concentration of the component will degrade and the next parameter ICBIO should be omitted.If NBIOFRM = 2, then two forms of this component degrade. These forms are specified using ICBIO.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>ICBIO</td>
<td>integer</td>
<td>Specify the aqueous component numbers (e.g. 1) or aqueous complex numbers (e.g. 101, 102, etc.) corresponding to the biodegradable form of the substrate.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDRXN = 5: Radioactive Decay Reaction</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>HALFLIFE</td>
<td>real</td>
<td>Half life (years)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>RXNTYPE</td>
<td>integer</td>
<td>Flag denoting the type of component participating in the reaction:  0 - Solid 1 - Liquid  -1  - Vapor</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>PARENT</td>
<td>integer</td>
<td>The number of the component which corresponds to the parent in the radioactive decay reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>DAUGHTER</td>
<td>integer</td>
<td>The number of the component which corresponds to the daughter in the radioactive decay reaction. If the simulation does not model the daughter product set daughter = 0.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IDRXN = 7: Kinetic Precipitation/Dissolution Reaction</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IIMMOBILE</td>
<td>integer</td>
<td>The immobile component number (e.g. 1, 2, etc.) which corresponds to the dissolving mineral</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NAQUEOUS</td>
<td>integer</td>
<td>The number of aqueous species which participate in the reaction</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IAQUEOUS</td>
<td>integer</td>
<td>The aqueous component numbers which correspond to the aqueous components which enter into the solubility product expression. Note that the total aqueous concentration of the component will dissolve.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>IMSTOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to the immobile component participating in the reaction. If positive, the solute is a reactant; if negative, the solute is a product; and if 0, the solute is not present in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>AQSTOIC</td>
<td>real</td>
<td>The stoichiometry corresponding to the aqueous components participating in the reaction. If positive, the solute is a reactant; if negative, the solute is a product; and if 0, the solute is not present in the reaction.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>SOLUBILITY</td>
<td>real</td>
<td>The solubility product. The units of the solubility product depend on the number of aqueous components participating in the reaction. For example, if there are two aqueous components participating the units would be (moles2/[kg water] 2)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character</td>
<td>Keyword ‘lookup’ designating a lookup table will be used to describe the solubility as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>NUM_TEMPS</td>
<td>integer</td>
<td>Number of values of temperature and solubility that will be used to describe solubility as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>EQTEMP</td>
<td>real</td>
<td>Temperatures for solubility function (oC)</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>TCOEFF</td>
<td>real</td>
<td>Solubilities for solubility as a function of temperature.</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>RATE</td>
<td>real</td>
<td>Reaction rate parameter (moles/[m2s])</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>SAREA</td>
<td>real</td>
<td>Surface area of the mineral (m2/m3 rock)</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>IDRXN = 8: Kinetic Precipitation/Dissolution Reaction (rates based on free-ion concentrations)
This model and its input are the same as for IDRXN = 7 except that the rates are based on the free-ion concentration
instead of the total concentration.</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>In general, the following examples illustrate only portions of the rxn macro and
putting all of these example inputs together will not result in a working FEHM <code class="docutils literal notranslate"><span class="pre">rxn</span></code>
macro. However, the dissolution example (the last example in this section)
provides an example of a complete rxn macro and corresponds to the first example
given for the trac macro.
In addition, the Reactive Transport Example and the Validation Test Plan for the
FEHM Application Version 2.30 (10086-VTP-2.30-00) include full example problems
with input files which demonstrate the use of rxn. These input files can
be used to see how rxn fits in with the other macros.</p>
<p>Specifically, the information in rxn must be consistent with the trac macro.
For example, if a linear kinetic sorption reaction is invoked by rxn,
a liquid component and solid component must be specified in the trac macro.</p>
<p><strong>General Reaction Parameters.</strong> In the following example two aqueous complexes and four kinetic reactions are specified. Three liquid components, Co, Fe, and EDTA, two complexes, CoEDTA and FeEDTA, and three immobile species, Co, CoEDTA and FeEDTA, are identified. Aqueous components 1 and 3 are coupled during solution while 2 is solved independently. Note that group 6 data has been omitted since there are no vapor species in this example.</p>
<p>[example]</p>
<p><strong>Equilibrium Reaction Parameters.</strong> Equilibrium speciation reactions modeled by
FEHM can be written in the following general form:</p>
<p><span class="math notranslate nohighlight">\(\sum_{j=1}^{N_c} a_{ij}\hat{C_j} \Leftrightarrow \hat{X_i} \;\;\;\; i = 1,...,N_x\)</span></p>
<p>where <span class="math notranslate nohighlight">\(\hat{C_j}\)</span> is the chemical formula for the aqueous component j, and
<span class="math notranslate nohighlight">\(\hat{X_i}\)</span> is the cchemical formula for the aqueous complex i, <span class="math notranslate nohighlight">\(a_{ij}\)</span>
is a stoichiometric coefficient <strong>{STOIC}</strong> giving the number of moles of component
j in complex i, and <span class="math notranslate nohighlight">\(N_x\)</span> is the number of aaqueous components.</p>
<p>Here is a simple example of an equilibrium speciation reaction:</p>
<p><span class="math notranslate nohighlight">\(aA + bB \Leftrightarrow A_aB_b\)</span></p>
<p>where a and b are the stoichiometric coefficients of components A and B, respectively.
At equilibrium, the concentrations of A, B, and <span class="math notranslate nohighlight">\(A_aB_b\)</span> must satisfy the
law of mass action for this reaction:</p>
<p><span class="math notranslate nohighlight">\(K = \frac{[A_aB_b]}{[A]^a[B]^b}\)</span></p>
<p>where K is the equilibrium constant <strong>{CKEQ}</strong> for the reaction and [X] is the
molar concentration of species X. The total aqueous concentrations of components
A and B are given by:</p>
<p><span class="math notranslate nohighlight">\(U_A = [A] + a[A_aB_b]\)</span></p>
<p><span class="math notranslate nohighlight">\(U_b = [B] + b[A_aB_b]\)</span></p>
<p>Therefore, the total aqueous concentrations, U, are functions of the uncomplexed component concentrations and the complex concentrations. Note that the kinetic reactions discussed in the next section are functions of the uncomplexed component concentrations and complex concentrations with the exception of the radioactive decay and precipitation/dissolution reaction in which the total aqueous concentration is used in the reaction rate law.</p>
<p>The following input describes the equilibrium speciation relations between the total aqueous components and aqueous complexes. Recall that the names of the components and aqueous complexes are given in Groups 3-6. The skip node option, ISKIP, is turned off and the tolerance for the speciation solver, RSDMAX, is set to 1e-9 which is the recommended value. The stoichiometry, STOIC, is specified so that the components Co and EDTA make up the first complex CoEDTA, and the components Fe and EDTA make up the second complex FeEDTA. The equilibrium constants, CKEQ, for CoEDTA and FeEDTA are 1.e18 and 6.31e27, respectively. The enthalpy change is 0.</p>
<table border="1" class="docutils">
<colgroup>
<col width="44%" />
<col width="19%" />
<col width="10%" />
<col width="2%" />
<col width="2%" />
<col width="2%" />
<col width="23%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>** ISKIP <strong>0</strong> RSDMAX1e-9** Chemical ** LOGKEQ0** CKEQ 1.0e186.31e27** STOIC <a href="#id1"><span class="problematic" id="id2">**</span></a>1.00.0</td>
<td>**&nbsp;Information **&nbsp;HEQ <a href="#id3"><span class="problematic" id="id4">**</span></a>00&nbsp;0.01.0</td>
<td>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.01.0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 7&nbsp;Group 8&nbsp;&nbsp;Group 9&nbsp;Group 10&nbsp;&nbsp;Group 11</td>
</tr>
</tbody>
</table>
<p>Below is an example of using the ‘lookup’ keyword in place of CKEQ and HEQ. This example describes the temperature dependence of the equilibrium constant of complex HCO3-. Note that LOGKEQ = 1 in this example, so log K is specified in the input rather than K.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">** LOGKEQ1** KEYWORDlookup</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>&#160;</td>
</tr>
<tr class="row-odd"><td>6.5</td>
</tr>
<tr class="row-even"><td>**&nbsp;NUM_TEMPS</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
</tr>
<tr class="row-even"><td>&#160;</td>
</tr>
<tr class="row-odd"><td>6.344</td>
</tr>
<tr class="row-even"><td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
</tr>
<tr class="row-even"><td>6.268</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
</tr>
<tr class="row-even"><td>0</td>
</tr>
<tr class="row-odd"><td>6.388</td>
</tr>
<tr class="row-even"><td>&#160;</td>
</tr>
<tr class="row-odd"><td>5</td>
</tr>
<tr class="row-even"><td>6.723</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
</tr>
<tr class="row-even"><td>0</td>
</tr>
<tr class="row-odd"><td>7.2</td>
</tr>
<tr class="row-even"><td>Group 9&nbsp;Group 10</td>
</tr>
</tbody>
</table>
<p>General Kinetic Parameters. In FEHM, eight kinetic reaction models are supported. Additional kinetic subroutines can be added without significant code development. The following is an example of input for the general kinetic parameters. A linear kinetic reaction, IDRXN = 1, is specified as occurring at each node in the problem (JA = 1, JB and JC = 0).</p>
<table border="1" class="docutils">
<colgroup>
<col width="31%" />
<col width="11%" />
<col width="15%" />
<col width="5%" />
<col width="5%" />
<col width="5%" />
<col width="29%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>** IDRXN <strong>1</strong> JA1</td>
<td>JB0</td>
<td>JC <a href="#id5"><span class="problematic" id="id6">**</span></a>0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12&nbsp;Group 13</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 1: Linear Kinetic Reaction.</strong> The retardation of contaminants due to
adsorption/desorption can be modeled with a linear kinetic sorption/desorption
expression. The rate of adsorption/desorption of component j is given by:</p>
<p><span class="math notranslate nohighlight">\(R_j = -k_m \left( c_j - \frac{m_j}{K_D} \right)\)</span></p>
<p>where <span class="math notranslate nohighlight">\(c_j\)</span> denotes the uncomplexed aqueous concentration of j <strong>{IAQUEOUS}</strong>,
<span class="math notranslate nohighlight">\(m_j\)</span> denotes the adsorbed concentration of species j <strong>{IIMMOBILE}</strong>,
<span class="math notranslate nohighlight">\(k_m\)</span> is the mass transfer coefficient <strong>{RATE}</strong>,
and <span class="math notranslate nohighlight">\(K_D\)</span> is the distribution coefficient <strong>{KD}</strong>.</p>
<p>As <span class="math notranslate nohighlight">\(k_m \rightarrow \infty\)</span>, this expression reduces to the linear equilibrium
isotherm. The following example illustrates input for a linear kinetic reaction.</p>
<p>In this kinetic reaction, aqueous component 1 adsorbs to form the immobile component 3.
The <span class="math notranslate nohighlight">\(K_D\)</span> for the reaction is 5.07 and the mass transfer coefficient is 1.0.</p>
<table border="1" class="docutils">
<colgroup>
<col width="39%" />
<col width="12%" />
<col width="9%" />
<col width="2%" />
<col width="2%" />
<col width="2%" />
<col width="34%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>** IDRXN <strong>1</strong> JA1** IAQ1** KD <strong>5.07</strong> RATE <a href="#id7"><span class="problematic" id="id8">**</span></a>1.0</td>
<td>JB0IMMOBILE3</td>
<td>JC <strong>0</strong></td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12&nbsp;Group 13&nbsp;Group 14&nbsp;Group 15&nbsp;Group 16</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 2: Langmuir Kinetic Reaction.</strong> The Langmuir kinetic reaction rate law is given by:</p>
<p><span class="math notranslate nohighlight">\(R_j = -k_m \frac{\rho}{\theta} \left( K_L c_j \left( m_j^{MAX} - m_j \right) - m_j \right)\)</span></p>
<p>where <span class="math notranslate nohighlight">\(k_m\)</span> is the rate constant for desorption <strong>{RATE}</strong>, <span class="math notranslate nohighlight">\(\rho\)</span> is the bulk rock density
<strong>{DENR}</strong>, <span class="math notranslate nohighlight">\(\theta\)</span> is the porosity <strong>{POR}</strong>. <span class="math notranslate nohighlight">\(K_L\)</span> is the distribution coefficient
<strong>{DISTCOEFF}</strong>, and <span class="math notranslate nohighlight">\(m_j^{MAX}\)</span> is the maximum concentration that can adsorb onto the j
solid <strong>{MAXCONC}</strong>. Ask <span class="math notranslate nohighlight">\(k_m \rightarrow \infty\)</span>, this expression reduces to the Langmuir equilibrium
isotherm. Example input for a Langmuir kinetic reaction follows.
In this kinetic reaction, aqueous complex 101 sorbs to form immobile component 1.
The distribution coefficient for the reaction is 2.e5 and the mass transfer coefficient
is 0.05. The maximum sorbed concentration is 1.69e-5.</p>
<table border="1" class="docutils">
<colgroup>
<col width="46%" />
<col width="9%" />
<col width="7%" />
<col width="2%" />
<col width="2%" />
<col width="2%" />
<col width="32%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>** IDRXN <strong>2</strong> JA1** IAQ101** DISTCO <strong>2.0e5</strong> RATE <strong>0.05</strong> MAXCO <a href="#id9"><span class="problematic" id="id10">**</span></a>1.69e-5</td>
<td>JB0IMMOBILE1</td>
<td>JC <strong>0</strong></td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12&nbsp;Group 13&nbsp;Group 14&nbsp;Group 15&nbsp;Group 16&nbsp;Group 17</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 3: General Kinetic Reaction.</strong> Many reactions fall under the category of the general kinetic reaction. The reaction is described by a forward rate constant {KFOR}, a reverse rate constant {KREV}, and a set of stoichiometric coefficients. The form of the general reversible reaction is given by:</p>
<p><span class="math notranslate nohighlight">\(\sum_{j=1}^{N_c+N_x+N_m+N_v} \eta_j \hat{Z}_j \Leftrightarrow \sum_{k=1}^{N_c+N_x+N_m+N_v} \eta_k' \hat{z}_k\)</span></p>
<p>where <span class="math notranslate nohighlight">\(N_c\)</span> is the number of aqueous components, <span class="math notranslate nohighlight">\(N_x\)</span> is the number
of aqueous complexes { <strong>NAQUEOUS</strong> = Nc + Nx }, <span class="math notranslate nohighlight">\(N_m\)</span> is the number of
immobile components <strong>{NIMMOBILE}</strong>, <span class="math notranslate nohighlight">\(N_v\)</span> is the number of vapor components
<strong>{NVAPOR}</strong>, <span class="math notranslate nohighlight">\(\eta_j\)</span> are reactant stoichiometric
coefficients, <span class="math notranslate nohighlight">\(\eta_k'\)</span> are product stoichiometric coefficients, and
<span class="math notranslate nohighlight">\(\hat{z}_i\)</span> is the chemical formula for species i, which may be an uncomplexed aqueous
component, aqueous complex, immobile component or vapor component.</p>
<p>The rate law for a general reversible reaction is given by the following expression:</p>
<p><span class="math notranslate nohighlight">\(R(z_i) = (\eta_i - \eta_i') \left[ k_f \prod_{j=1}^{N_c+N_x+N_m+N_v} z_j^{\eta_j} - k_r \prod_{k=1}^{N_c+N_x+N_m+N_v} z_k^{\eta_k'} \right]\)</span></p>
<p>where <span class="math notranslate nohighlight">\(z_i\)</span> is the concentration of species i.</p>
<p>The following is example input for a
general kinetic reaction. Solid component 1 reacts to form solid components 2 and 3.
The forward reaction rate is 1.26e-2 and the reverse reaction rate is 0. Therefore, this
is an irreversible reaction. Note also that only solid components are reacting so groups 18-21 have been omitted.</p>
<table border="1" class="docutils">
<colgroup>
<col width="34%" />
<col width="16%" />
<col width="13%" />
<col width="4%" />
<col width="2%" />
<col width="2%" />
<col width="29%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>** IDRXN <strong>3</strong> JA1** NIMM3** KFOR1.26e-2** IIMM <strong>1</strong> IMSTOIC1.0</td>
<td>JB0NAQSP0KREV <strong>0.0&nbsp;2</strong>-1.0</td>
<td>JC <a href="#id11"><span class="problematic" id="id12">**</span></a>0NVAPOR0&nbsp;&nbsp;&nbsp;3&nbsp;-1.0</td>
<td><a href="#id13"><span class="problematic" id="id14">**</span></a></td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12&nbsp;Group 13&nbsp;Group 14&nbsp;Group 15&nbsp;Group 16&nbsp;Group 17</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 4: Dual Monod Biodegradation Reaction.</strong> Biodegradation is an irreversible
process in which bacteria oxidize an organic substrate to produce energy and biomass.
In addition to biomass, the biodegradation process requires the presence of an electron
acceptor (e.g. oxygen, nitrate, etc.) and nutrients (e.g. nitrogen and phosphorous).
An example of a biodegradation reaction is given by the following reaction:</p>
<p><span class="math notranslate nohighlight">\(\mathrm{Substrate} + \mathrm{Electron Acceptor} + H^+ + \mathrm{Nutrients} \rightarrow \mathrm{cells} + CO_2 + H_2O + NH_3\)</span></p>
<p>FEHM models the rate of biodegradation of a substrate with a multiplicative Monod model, which is given by:</p>
<p><span class="math notranslate nohighlight">\(R_s = -q_m m_b \frac{[S]}{K_s + [S]K_A} \frac{[A]}{K_A + [A]}\)</span></p>
<p>where [S] is the aqueous concentration of substrate (a.k.a the electron donor) <strong>{SUBSTRATE}</strong>,
[A] is the aqueous concentration of the electron acceptor <strong>{ELECACC}</strong>,
and <span class="math notranslate nohighlight">\(m_b\)</span> is the concentration of the immobile biomass <strong>{BIOMASS}</strong>.</p>
<p>The parameter qm is the maximum specific rate of substrate utilization <strong>{QM}</strong>,
which represents the maximum amount of substrate that can be consumed per unit
mass of bacteria per unit time. The parameters <span class="math notranslate nohighlight">\(K_S\)</span> <strong>{KS}</strong> and <span class="math notranslate nohighlight">\(K_A\)</span>
<strong>{KA}</strong> are the Monod half-maximum-rate concentrations for the electron
donor and electron acceptor, respectively. The rate of microbial growth is
given by the synthesis rate (which is proportional to the rate of substrate
degradation) minus a first-order decay rate.</p>
<p><span class="math notranslate nohighlight">\(R_{cells} = -YR_s - b(m_b - m_{b,init})\)</span></p>
<p>where Y is the microbial yield coefficient <strong>{YIELD}</strong> and b is the first-order
microbial decay coefficient <strong>{DECAY}</strong>. In the above equation,
the assumption is made that the background conditions are sufficient to
sustain a microbial population of a given size; therefore, the biomass
concentration is not allowed to fall below its initial background concentration
(<span class="math notranslate nohighlight">\(m_{b,init}\)</span>) <strong>{XMINIT}</strong>.
In the following example of input for a dual monod biodegradation reaction,
aqueous component 3 is the substrate and aqueous component 4 is the electron acceptor.
Note that there are only two entries in Group 15 and Groups 21-23 are omitted since
NAQUEOUS = 2. In addition Group 29 data is left out since NBIOFRM = 0.</p>
<table border="1" class="docutils">
<colgroup>
<col width="41%" />
<col width="13%" />
<col width="4%" />
<col width="1%" />
<col width="1%" />
<col width="13%" />
<col width="27%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">** IDRXN <strong>4</strong> JA1** NAQSP 2** SUBSTRA3** BIOMASS4** KS <strong>0.201e-3</strong> KA ** 0.00625e-3** DECAY <strong>0.0020833</strong> EASTOIC3.10345</th>
<th class="head">JB0NIMMOBILE1ELECACC <strong>4</strong>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**</th>
<th class="head">JC <strong>0</strong></th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">Group 12&nbsp;Group 13&nbsp;Group 14&nbsp;Group 15&nbsp;Group 16&nbsp;Group 17&nbsp;Group 18&nbsp;Group 19&nbsp;Group 20</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>** COMP3STOIC <strong>0</strong> COMP4STOIC <strong>0</strong> COMP5STOIC <strong>0</strong> PHTHRSH <a href="#id15"><span class="problematic" id="id16">**</span></a>8</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 21&nbsp;Group 22&nbsp;Group 23&nbsp;Group 24</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>** QM <a href="#id17"><span class="problematic" id="id18">**</span></a>8.0226 e-4 ** YIELD <a href="#id19"><span class="problematic" id="id20">**</span></a>44.8732</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 25&nbsp;Group 26</td>
</tr>
<tr class="row-even"><td>** XMINIT <strong>0.0</strong> NBIOFRM0** ICBIO **</td>
<td><a href="#id21"><span class="problematic" id="id22">**</span></a></td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 27&nbsp;Group 28&nbsp;Group 29</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 5: Radioactive Decay Reaction.</strong> Radioactive decay is a simple first order decay process given by:</p>
<p><span class="math notranslate nohighlight">\(A \rightarrow B\)</span></p>
<p>where A is the parent <strong>{PARENT}</strong> and B is the daughter product <strong>{DAUGHTER}</strong>.
The half life of the reaction is defined as the time it takes for the concentration
of A to decrease by a factor of 2. In the following example of input for a
radioactive decay reaction, aqueous component 2, the parent, reacts to form aqueous
component 3, the daughter product. The half life for the reaction is 432.0 years.</p>
<table border="1" class="docutils">
<colgroup>
<col width="38%" />
<col width="15%" />
<col width="11%" />
<col width="2%" />
<col width="2%" />
<col width="2%" />
<col width="31%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>** IDRXN <strong>5</strong> JA1** HALFLIFE432.0** RXNTYPE1** PARENT2</td>
<td>JB0**&nbsp;**&nbsp;DAUGHTER3</td>
<td>JC <a href="#id23"><span class="problematic" id="id24">**</span></a>0&nbsp;&nbsp;&nbsp;&nbsp;**</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12&nbsp;Group 13&nbsp;Group 14&nbsp;Group 15&nbsp;Group 16</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 7: Kinetic Precipitation/Dissolution Reaction.</strong> A general reaction describing
the precipitation/dissolution of a mineral p <strong>{IIMMOBILE}</strong> can be written in the
following form:</p>
<p><span class="math notranslate nohighlight">\(m_p \Leftrightarrow \mu_{p1}c_1 + \mu_{p1}c_2 + ... + \mu_{p1}c_{N_c}\)</span></p>
<p>where the <span class="math notranslate nohighlight">\(c_j\)</span> are the aqueous concentrations <strong>{IAQUEOUS}</strong> and the <span class="math notranslate nohighlight">\(\mu_{pj}\)</span>
are stoichiometric coefficients <strong>{AQSTOIC}</strong>.
The equilibrium constant for this reaction is known as the solubility product.
Since the activity of a pure solid is equal to one, the reaction quotient <span class="math notranslate nohighlight">\(Q_p\)</span>
is defined as follows:</p>
<p><span class="math notranslate nohighlight">\(Q_p = \prod_{j=1}^{N_c} c_j^{\mu_{pj}}\)</span></p>
<p>At equilibrium, <span class="math notranslate nohighlight">\(Q_p\)</span> is equal to the solubility product. The surface-controlled
rate of precipitation/dissolution of a mineral is given by:</p>
<p><span class="math notranslate nohighlight">\(R(m_p) = sign\left(\log\frac{Q_p}{K_{sp}} \right) A_p k_p \left| \frac{Q_p}{K_{sp}} -1 \right|\)</span></p>
<p>where Ap is reactive surface area of the mineral <strong>{AREA}</strong>, <span class="math notranslate nohighlight">\(k_p\)</span> is the precipitation rate
constant <strong>{RATE}</strong>, and <span class="math notranslate nohighlight">\(K_{sp}\)</span> is the solubility product <strong>{SOLUBILITY}</strong>.
Currently, this precipitation/dissolution subroutine only allows for the total
aqueous concentration of a component to dissolve.
The dissolution of uncomplexed aqueous concentration and complex concentrations
is not currently supported.</p>
<p>The following is example input for a kinetic precipitation/dissolution reaction (calcite dissolution). This example corresponds to the example used for the trac macro (<a class="reference external" href="Macro37141.html">See In the following example of trac, calcite dissolution is simulated. The input groups are given to the right of the table to facilitate review of the example. The initial solute concentration is set to 0 (but is later overwritten by group 14 input), the implicitness factor is 1 resulting in a 1st order solution, the equation tolerance is 1.e-7, and the upstream weighting is set to 0.5. The solute transport solution is turned on as the heat and mass solution is turned off at day 1. The heat and mass solution resumes on day 1000. Two solutes are simulated in this example. Solute 1 is a nonsorbing (conservative) liquid species (a1 = a2 = 0., b = 1.) with a molecular diffusion coefficient of 1.e-9 m2/s, and dispersivity of 0.0067 m in the X-direction. This transport model applies over the entire problem domain. The initial solute concentration for solute 1 is 6.26e-5 mol/kg-water. Solute 2 is a solid species with an initial solute concentration of 2.e-5 mol/kg-solid. There is no solute source for either solute. The corresponding data for the rxn macro, that would complete this example, is given on page 159.</a>). A single kinetic reaction is specified with no aqueous complexes. The liquid and immobile components are identified. Aqueous component 1 participates in the precipitation/dissolution reaction with immobile component 1. No data is input for groups 4, 6, 9, 10, and 11. The solubility product for the reaction is 6.26e-5 and the reaction rate constant is 100.</p>
<table border="1" class="docutils">
<colgroup>
<col width="27%" />
<col width="19%" />
<col width="15%" />
<col width="7%" />
<col width="8%" />
<col width="4%" />
<col width="20%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>rxn** NCPLX0** GROUP <strong>11</strong> IDCPNT1** IDCPLX** IDIMM1** IDVAP</td>
<td>NUMRXN <a href="#id25"><span class="problematic" id="id26">**</span></a>1&nbsp;&nbsp;&nbsp;CPNTNAMCa[aq]CPLXNAMIMMNAMCa[s]VAPNAM</td>
<td>IFXCONC0CPLXPRT IMMPRT <a href="#id27"><span class="problematic" id="id28">**</span></a>0VAPPRT</td>
<td>CPNTPRT 0**</td>
<td>CPNTGS 1.0e-9</td>
<td><a href="#id29"><span class="problematic" id="id30">**</span></a></td>
<td>Group 1&nbsp;Group 2&nbsp;&nbsp;Group 3&nbsp;Group 4Group 5&nbsp;Group 6</td>
</tr>
<tr class="row-even"><td>** ISKIP <strong>0</strong> RSDMAX1e-13** Chemical ** LOGKEQ** CKEQ ** STOIC **</td>
<td>**&nbsp; reaction**HEQ **</td>
<td>information</td>
<td><a href="#id31"><span class="problematic" id="id32">**</span></a></td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 7&nbsp;Group 8&nbsp;&nbsp;Group 9Group 10Group 11</td>
</tr>
<tr class="row-odd"><td>** IDRXN <strong>7</strong> JA1</td>
<td>JB0</td>
<td>JC <a href="#id33"><span class="problematic" id="id34">**</span></a>0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 12&nbsp;Group 13</td>
</tr>
<tr class="row-even"><td>** IIMMOBLE1** NAQSP <strong>1</strong> IAQ <strong>1</strong> IMSTOIC1** AQSTOIC1** SOLPROD6.26e-5</td>
<td>**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**&nbsp;**&nbsp;**</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 14&nbsp;Group 15&nbsp;Group 16&nbsp;Group 17&nbsp;Group 18&nbsp;Group 19</td>
</tr>
<tr class="row-odd"><td>** RATE <strong>100</strong> AREA <a href="#id35"><span class="problematic" id="id36">**</span></a>1.0</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>Group 20&nbsp;Group 21</td>
</tr>
</tbody>
</table>
<p><strong>IDRXN = 8: Kinetic Precipitation/Dissolution Reaction (rates based on free-ion concentrations).</strong>
The reaction modeled is analogous to that for IDRXN =7, except that rates are based
on the uncomplexed (free-ion) concentration of the species. The total concentration
is equal to the free ion concentration + all of the complex concentrations.
For a more detailed discussion of the differences between total aqueous and free-ion
concentration, see the “Models and Methods Summary” of the
FEHM Application [Zyvoloski et al. 1999, [[FEHM-MMS.htm#14829|]], [[FEHM-MMS.htm#14829|]]].</p>
<p>For example, for a species such as Cobalt from the multisolute problem
(see [[FEHM-UM.9.5.htm#76956|See Reactive Transport Example]]),
the free ion concentration is simply the concentration of Cobalt in it’s
uncomplexed state. The total Cobalt would be Free Ion Cobalt + all Cobalt Complexes
(e.g. CoEDTA from the multisolute verification problem). Using IDRXN = 7
allows the total Cobalt to dissolve, while IDRXN= 8 allows only the
free ion Cobalt to dissolve.</p>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
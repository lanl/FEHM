---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">iter</span></code><a class="headerlink" href="#iter" title="Permalink to this headline">¶</a></h1>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">This control statement is optional but recommended. However, if the user is not
familiar with the linear equation solver routines in FEHM (Zyvoloski and Robinson, 1995)
control statement <strong>iter</strong> should not be used.</p>
</div>
<ul class="simple">
<li>Group 1 -     G1, G2, G3, TMCH, OVERF</li>
<li>Group 2 -     IRDOF, ISLORD, IBACK, ICOUPL, RNMAX</li>
</ul>
<p>The parameters G1, G2, and G3 are used to calculate the completion criteria for
the linear equation solver. The equation for the stopping criteria is:</p>
<p><span class="math notranslate nohighlight">\(EPE = G3 * \mathrm{max}(TMCH, \mathrm{MAX}(F0, \mathrm{MIN(G1*\mathrm{SQRT(R^2), G2*R^2})}))\)</span></p>
<p>where <span class="math notranslate nohighlight">\(R^2\)</span> is the sum-squared of the equation residuals, and F0 is the
<span class="math notranslate nohighlight">\(SQRT(R0^2)*EPM\)</span> for the first iterration (see macro <code class="docutils literal notranslate"><span class="pre">ctrl</span></code> for a
definition of EPM). The other parameters are defined below:</p>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="11%" />
<col width="9%" />
<col width="65%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>G1</td>
<td>real</td>
<td>1.e-6</td>
<td>Multiplier for the linear convergence region of the
Newton-Raphson iteration.</td>
</tr>
<tr class="row-odd"><td>G2</td>
<td>real</td>
<td>1.e-6</td>
<td>Multiplier for the quadratic convergence region of
the Newton-Raphson iteration.</td>
</tr>
<tr class="row-even"><td>G3</td>
<td>real</td>
<td>1.e-3</td>
<td>Multiplier relating Newton Raphson residuals to
stopping criteria for linear solver</td>
</tr>
<tr class="row-odd"><td>TMCH</td>
<td>real</td>
<td>1.e-9</td>
<td>Machine tolerance if TMCH &gt; 0. If satisfied by
the residual norm, the Newton iteration is assumed to
be complete. Newton-Raphson stopping criteria if TMCH &lt; 0
(recommended). If TMCH &lt; 0 then the ABS(TMCH) is used as a
tolerance for each equation at each node. Convergence is
achieved if the residual of every equation at every
node is &lt; ABS(TMCH).</td>
</tr>
<tr class="row-even"><td>OVERF</td>
<td>real</td>
<td>1.1</td>
<td>Over relaxation factor for passive nodes in adaptive
implicit method.</td>
</tr>
<tr class="row-odd"><td>IRDOF</td>
<td>integer</td>
<td>0</td>
<td><div class="first last line-block">
<div class="line">Enables the reduced degree of freedom method. If
IRDOF = 0, reduced degrees of freedom are not required.
When IRDOF = 1, a reduced degree of freedom from 3 to 2
or 3 to 1 is used. When IRDOF = 2, a reduced degree of
freedom from 3 to 2 is used. If IRDOF =11, then an air
only solution is found for the isothermal air-water
process model. If IRDOF = -11, then the residual for
the air equation with the airwater macro is ignored.
If IRDOF = 13, then a liquid only solution for the
airwater macro is assumed. {0}</div>
<div class="line"><br /></div>
<div class="line">Examples of 1, 2, 3, 4 and 6 degrees of freedom models are:</div>
<div class="line">1 - heat only or mass only.</div>
<div class="line">2 - heat and mass, or air-water (isothermal)</div>
<div class="line">3 - air-water with heat (non-isothermal)</div>
<div class="line">4 - heat and mass, double permeability or air-water
(isothermal), double permeability</div>
<div class="line">6 - air-water with heat, double permeability</div>
<div class="line"><br /></div>
<div class="line">See Tseng and Zyvoloski (2000) for more information
on the reduced degree of freedom method.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>ISLORD</td>
<td>integer</td>
<td>0</td>
<td>Reordering parameter. The value of ISLORD and the
corresponding equation order is given below. The ordering
has an effect on the speed of convergence of several
solution algorithms, but will not affect most users.
For problems of order 2 or greater, the ordering can be
understood by labeling each equation. For example for a
3-degree of freedom problem with mass, heat, and noncondensible
gas, label the mass equation as 1, the heat equation as 2,
and the noncondensible gas equation as 3. In general mass
(water), heat or air, air. For double permeability problems
fracture equations precede matrix equations, i.e., for an
air-water problem - mass water fracture, mass air fracture,
mass water matrix, mass air matrix. {0}</td>
</tr>
<tr class="row-odd"><td>IBACK</td>
<td>integer</td>
<td>0</td>
<td>IRDOF parameter. If IBACK = 0, SOR iterations are not
performed before call to solver. If IBACK = 1,
SOR iterations are performed before call to solver.
If IBACK = 2, SOR iterations are performed before
call to SOLVER, and SOLVER is called twice. {0}</td>
</tr>
<tr class="row-even"><td>ICOUPL</td>
<td>integer</td>
<td>0</td>
<td>Number of SOR iterations used in reduced degree of
freedom methods. {0}</td>
</tr>
<tr class="row-odd"><td>RNMAX</td>
<td>real</td>
<td>1.0e+11</td>
<td>Maximum running time for problem before the solution
is stopped (cpu minutes).</td>
</tr>
</tbody>
</table>
<div class="section" id="id1">
<h2>{0}<a class="headerlink" href="#id1" title="Permalink to this headline">¶</a></h2>
<table border="1" class="docutils">
<colgroup>
<col width="8%" />
<col width="23%" />
<col width="23%" />
<col width="23%" />
<col width="23%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>ISLORD</td>
<td>2 Degrees of Freedom</td>
<td>3 Degrees of Freedom</td>
<td>4 Degrees of Freedom</td>
<td>6 Degrees of Freedom</td>
</tr>
<tr class="row-even"><td>0</td>
<td>1,2</td>
<td>1,2,3</td>
<td>1,2,3,4</td>
<td>1,2,3,4,5,6</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>2,1</td>
<td>1,3,2</td>
<td>1,3,2,4</td>
<td>1,4,2,5,3,6</td>
</tr>
<tr class="row-even"><td>2</td>
<td>&#160;</td>
<td>2,1,3</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>3</td>
<td>&#160;</td>
<td>2,3,1</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">iter</span></code>.
In this example, the tolerances for the linear and quadratic convergence regions
for the Newton-Raphson method are specified to be 1.e-5 times the initial residual,
tolerance for the adaptive-implicit method is 1.e-5, machine tolerance is 1.e-9,
and over-relaxation factor is 1.2. The reduced degree of freedom method is enabled,
reordering is not done, SOR iterations are not performed before calling the solver,
two SOR iterations are used in the reduced degree of freedom method, and the solution
procedure is terminated if not completed within 200 CPU minutes.</p>
<table border="1" class="docutils">
<colgroup>
<col width="20%" />
<col width="20%" />
<col width="20%" />
<col width="20%" />
<col width="20%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td colspan="5">iter</td>
</tr>
<tr class="row-even"><td>1.e-5</td>
<td>1.e-5</td>
<td>1.e-5</td>
<td>1.e-9</td>
<td>1.2</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>0</td>
<td>0</td>
<td>2</td>
<td>200.0</td>
</tr>
</tbody>
</table>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
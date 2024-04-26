---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">ctrl</span></code><a class="headerlink" href="#ctrl" title="Permalink to this headline">¶</a></h1>
<p>Assign various control parameters needed for equation solvers and matrix solver routines. Suggested values for the control parameters are shown in <strong>“{ }”</strong> in the table. For older input files where MAXSOLVE and ACCM were not input, the default is <code class="docutils literal notranslate"><span class="pre">ACCM</span> <span class="pre">=</span> <span class="pre">gmre</span></code> and <code class="docutils literal notranslate"><span class="pre">MAXSOLVE</span> <span class="pre">=</span> <span class="pre">3*NORTH</span></code>.</p>
<ul class="simple">
<li>Group 1 - MAXIT, EPM, NORTH, MAXSOLVE, ACCM</li>
<li>Group 2 - <a class="reference external" href="InputData.html#JA">JA, JB, JC, NAR</a></li>
<li>Group 3 - AAW, AGRAV, UPWGT</li>
<li>Group 4 - IAMM, AIAA, DAYMIN, DAYMAX</li>
<li>Group 5 - ICNL, LDA</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="1%" />
<col width="1%" />
<col width="1%" />
<col width="98%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>Input Variable</td>
<td>Format</td>
<td>Default</td>
<td>Description</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">MAXIT</span></code></td>
<td>integer</td>
<td>&#160;</td>
<td>Maximum number of iterations allowed in either the overall Newton cycle or the inner cycle to solve for the corrections at each iteration. If <code class="docutils literal notranslate"><span class="pre">MAXIT</span> <span class="pre">&lt;</span> <span class="pre">0</span></code> then the maximum number of iterations is <code class="docutils literal notranslate"><span class="pre">ABS(MAXIT)</span></code> but the minimum number of iterations is set to 2. {10}</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">EPM</span></code></td>
<td>real</td>
<td>&#160;</td>
<td>Tolerance for Newton cycle (nonlinear equation tolerance). Note - EPM gets overwritten by TMCH in ITER macro if that variable is defined. {1.e-5}</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">NORTH</span></code></td>
<td>integer</td>
<td>&#160;</td>
<td>Number of orthogonalizations in the linear equation solver. Note - for more complicated problems, increase <code class="docutils literal notranslate"><span class="pre">NORTH</span></code>. For example, for fully coupled stress problems recommend using a value of 80 and gmre. {8 for gmre, 1 for bcgs}</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">MAXSOLVE</span></code></td>
<td>integer</td>
<td>&#160;</td>
<td>Maximum number of solver iterations per Newton iteration allowed. {100}</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">ACCM</span></code></td>
<td>character*4</td>
<td>&#160;</td>
<td>Acceleration method for solver bcgs - Biconjugate gradient stabilized acceleration. Recommended for isothermal steady-state saturated flow problems. gmre - Generalized minimum residual acceleration. Recommended for all other types of problems.</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">NAR</span></code></td>
<td>integer</td>
<td>1</td>
<td>The order of partial Gauss elimination {1 or 2 is recommended}. Larger values increase memory utilization but may be necessary for convergence.</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">AAW</span></code></td>
<td>real</td>
<td>&#160;</td>
<td>Implicitness factor. {1}AAW ‚â§ 1, use standard pure implicit formulation.AAW &gt; 1, use second-order implicit method.</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">AGRAV</span></code></td>
<td>integer</td>
<td>&#160;</td>
<td>Direction of gravity AGRAV = 0, no gravity is used.AGRAV = 1, X-direction.AGRAV = 2, Y-direction.AGRAV = 3, Z-direction.A value for gravity of 9.81 m/s2 is used in the code when AGRAV ‚â  0. If AGRAV &gt; 3, AGRAV is set equal to 3.</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">UPWGT</span></code></td>
<td>real</td>
<td>&#160;</td>
<td>Value of upstream weighting {1.0}.If UPWGT &lt; 0.5, UPWGT is set to 0.5If UPWGT &gt; 1.0, UPWGT is set to 1.0</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">IAMM</span></code></td>
<td>integer</td>
<td>&#160;</td>
<td>Maximum number of iterations for which the code will multiply the time step size. If this number of time steps is exceeded at any time, the time step will not be increased for the next time. Set IAMM &lt; MAXIT {7-10}.</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">AIAA</span></code></td>
<td>real</td>
<td>1</td>
<td>Time step multiplier {1.2-2.0}</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">DAYMIN</span></code></td>
<td>real</td>
<td>1.0e-05</td>
<td>Minimum time step size (days)</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">DAYMAX</span></code></td>
<td>real</td>
<td>30.0</td>
<td>Maximum time step size (days)</td>
</tr>
<tr class="row-odd"><td><code class="docutils literal notranslate"><span class="pre">ICNL</span></code></td>
<td>integer</td>
<td>&#160;</td>
<td>Parameter that specifies the geometry: <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">0</span></code>, three-dimensional. <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">1</span></code>, X - Y plane. <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">2</span></code>, X - Z plane. <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">3</span></code>, Y - Z plane. <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">4</span></code>, X - Y radial plane, (radius is X). <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">5</span></code>, X - Z radial plane, (radius is X). <code class="docutils literal notranslate"><span class="pre">ICNL</span> <span class="pre">=</span> <span class="pre">6</span></code>, Y - Z radial plane, (radius is Y)</td>
</tr>
<tr class="row-even"><td><code class="docutils literal notranslate"><span class="pre">LDA</span></code></td>
<td>integer</td>
<td>0</td>
<td>Parameter that specifies the external storage of geometric coefficients: <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">-2</span></code>, element coefficients are calculated in the code and saved, unformatted, on file filen.stor. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">-1</span></code>, element coefficients are calculated in the code and saved on file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code>. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">0</span></code>, element coefficients are calculated in the code and not saved. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">+1</span></code>, element coefficients are read from file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code> and no coefficients are calculated in the code. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">+2</span></code>, element coefficients are read, unformatted, from file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code> and no coefficients are calculated in the code. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">+5</span></code>, element coefficients are read from file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code> and no coefficients are calculated in the code. Coefficients are re-saved unformatted, on <code class="docutils literal notranslate"><span class="pre">filen_UNF.stor</span></code>. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">+6</span></code>, element coefficients are read from file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code> and no coefficients are calculated in the code. Coefficients are re-saved formatted, on <code class="docutils literal notranslate"><span class="pre">filen_FOR.stor</span></code>. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">+7</span></code>, element coefficients are read, unformatted, from file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code> and no coefficients are calculated in the code. Coefficients are re-saved unformatted, on <code class="docutils literal notranslate"><span class="pre">filen_UNF.stor</span></code>. <code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">=</span> <span class="pre">+8</span></code>, element coefficients are read, unformatted, from file <code class="docutils literal notranslate"><span class="pre">filen.stor</span></code> and no coefficients are calculated in the code. Coefficients are re-saved formatted, on <code class="docutils literal notranslate"><span class="pre">filen_FOR.stor</span></code>. It should be noted that if the coefficients are read from a file (<code class="docutils literal notranslate"><span class="pre">LDA</span> <span class="pre">&gt;</span> <span class="pre">0</span></code>) then the macro nfinv is ignored as well as information read from macros elem and coor since the coefficients are not being calculated.</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">ctrl</span></code>. In this example, the maximum number of iterations allowed is 40, tolerance for solving the nonlinear equations using Newton iterations is 1.e-7, and the number of orthogonalizations in the linear equation solver is 8. The order of partial Gauss elimination for all nodes 1 through 140 is 1. A forward implicit formulation is used for the time derivative, there is no gravity, and full upstream weighting is used. The number of iterations for which the time step is increased is 40, the time step is increased by a factor of 1.2 at each iteration, the minimum time step size is 0.1 days, and the maximum time step size is 60 days. The geometry of the problem is 2-dimensional in the X-Y plane and the finite element coefficients are calculated during the run and not saved.</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">ctrl</span>
<span class="mi">40</span>
<span class="mi">1</span>
<span class="mf">1.e-7</span>
<span class="mi">140</span>
<span class="mi">8</span>
<span class="mi">1</span>
<span class="mi">24</span>
<span class="mi">1</span>
<span class="n">gmre</span>
<span class="mf">1.0</span>
<span class="mi">40</span>
<span class="mi">1</span>
<span class="mf">0.0</span>
<span class="mf">1.2</span>
<span class="mi">00</span>
<span class="mf">1.00</span>
<span class="mf">0.1</span>
<span class="mf">60.0</span>
</pre></div>
</div>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">stea</span></code><a class="headerlink" href="#stea" title="Permalink to this headline">¶</a></h1>
<p>The macro <code class="docutils literal notranslate"><span class="pre">stea</span></code> is used to manage steady state simulations and is available with all
physics modules in FEHM. The macro directs FEHM to monitor changes in variables
from timestep to timestep and stop the steady state run when the changes are less
than some prescribed tolerance. Alternatively the steady state run is directed
to finish when the global “in” and “out” fluxes are less than a prescribed
tolerance or the simulated time exceeds the input value.</p>
<p>After the steady state portion of the simulation is completed, a transient run may
be performed. This is accomplished with the boun macro and the key word <code class="docutils literal notranslate"><span class="pre">tran</span></code>.
See the description of the <code class="docutils literal notranslate"><span class="pre">boun</span></code> macro for details.</p>
<p>The user should be aware that when the <code class="docutils literal notranslate"><span class="pre">stea</span></code> macro is used, the parameters
associated with the <code class="docutils literal notranslate"><span class="pre">time</span></code> macro pertain to the transient portion of the
simulation if a transient part exists. Values for these parameters may be
input using a keyword but if not entered will default to the values specified
for the <code class="docutils literal notranslate"><span class="pre">time</span></code> macro.</p>
<ul class="simple">
<li>Group 1 - KEYWORD, VALUE</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="15%" />
<col width="10%" />
<col width="75%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KEYWORD</td>
<td>character</td>
<td><div class="first last line-block">
<div class="line">The following keywords are used with steady to specify
the variables to be checked for steady state:</div>
<div class="line">shea - Head (m)</div>
<div class="line">spre - Pressure (MPa)</div>
<div class="line">stem - Temperature (oC)</div>
<div class="line">ssat - Saturation</div>
<div class="line">sair - Partial pressure of air/gas (MPa)</div>
<div class="line">sflu - Mass flux (kg/s)</div>
<div class="line">sent - Enthalpy (MJ/s)</div>
<div class="line">stim - Maximum time for steady state simulation (days)</div>
<div class="line">sday - Initial time step size for steady state simulation (days)</div>
<div class="line">smul - Time step multiplication factorsmst - Minimum
number of time steps to be used for steady state simulation</div>
<div class="line">snst - Maximum number of time steps to be used for steady
state simulation</div>
<div class="line">shtl - Option to reduce the head_tol factor as the
solution approaches steady-state</div>
<div class="line">stmc - Option to reduce the machine tolerancs factor (tmch)
factor as the solution approaches steady-state</div>
<div class="line">sacc - Maximum change allowed in the accumulation term when
flux is being checked</div>
<div class="line">sper - The tolerance is interpreted as a fractional change
in the variable being checked [i.e.,
(new_value - old_value)/old_value]. Without this keyword
it is an absolute change in the variable value.</div>
<div class="line">endstea - Signifies end of keyword input, a blank line will also work.</div>
</div>
</td>
</tr>
<tr class="row-odd"><td>VALUE</td>
<td>real</td>
<td>Variable tolerance or time control parameter value.</td>
</tr>
</tbody>
</table>
<p>In the following example a steady state solution is specified. The tolerance for
<code class="docutils literal notranslate"><span class="pre">head</span></code> is specified to be 0.1 m and for <code class="docutils literal notranslate"><span class="pre">flux</span></code> 0.00001kg/s. The steady state solution
will be allowed to run for a maximum of 1.e12 days and the time step multiplier
is set to 2.</p>
<table border="1" class="docutils">
<colgroup>
<col width="53%" />
<col width="47%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">stea</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>shead</td>
<td>1.d-1</td>
</tr>
<tr class="row-odd"><td>stime</td>
<td>1.e12</td>
</tr>
<tr class="row-even"><td>smult</td>
<td><ol class="first last arabic simple" start="2">
<li></li>
</ol>
</td>
</tr>
<tr class="row-odd"><td>sflux</td>
<td>1.d-5</td>
</tr>
<tr class="row-even"><td>end</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">time</span></code><a class="headerlink" href="#time" title="Permalink to this headline">¶</a></h1>
<p>Time step and time of simulation data.</p>
<ul class="simple">
<li>Group 1 - DAY, TIMS, NSTEP, IPRTOUT, YEAR, MONTH, INITTIME</li>
<li>Group 2 - DIT1, DIT2, DIT3, ITC, DIT4 (as needed)</li>
</ul>
<p>DAY should be larger than DAYMIN defined in control statement <strong>ctrl</strong>.
The code proceeds to the next control statement when a blank line is encountered for Group 2.
Group 2 can be used to generate output at specific times (with multiple Group 2s).
Contour plot output will be written at each DIT1 regardless of the input in control statement <strong>cont</strong>.
The restart file will be written (or rewritten if one already exists) at each DIT1. If DIT4 is omitted
(for compatibility with older input files where DIT4 was not input) the maximum
time step defined in the control statement ctrl will be used.</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="10%" />
<col width="72%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>DAY</td>
<td>real</td>
<td>Initial time step size (days).</td>
</tr>
<tr class="row-odd"><td>TIMS</td>
<td>real</td>
<td>Final simulation time (days).</td>
</tr>
<tr class="row-even"><td>NSTEP</td>
<td>integer</td>
<td>Maximum number of time steps allowed.</td>
</tr>
<tr class="row-odd"><td>IPRTOUT</td>
<td>integer</td>
<td>Print-out interval for nodal information (pressure,
enthalpy etc.), as set up under control statement node.
(i.e., number of time steps).</td>
</tr>
<tr class="row-even"><td>YEAR</td>
<td>integer</td>
<td>Year that simulation starts.</td>
</tr>
<tr class="row-odd"><td>MONTH</td>
<td>integer</td>
<td>Month that simulation starts.</td>
</tr>
<tr class="row-even"><td>INITTIME</td>
<td>real</td>
<td>Initial time of simulation (days). For compatibility
with older versions, if this parameter is absent the
initial time of simulation will be 0 if no restart
file is used, or the time in the restart file if one is used.</td>
</tr>
<tr class="row-odd"><td>DIT1</td>
<td>real</td>
<td>Time (days) for time step change.</td>
</tr>
<tr class="row-even"><td>DIT2</td>
<td>real</td>
<td>New time step size (days). If DIT2 &lt; 0 then ABS (DIT2)
is the new time step multiplier.</td>
</tr>
<tr class="row-odd"><td>DIT3</td>
<td>real</td>
<td><div class="first last line-block">
<div class="line">Implicitness factor for new time step.</div>
<div class="line">DIT3 ≤ 1.0 backward Euler.</div>
<div class="line">DIT3 &gt; 1.0 for second-order implicit scheme.</div>
</div>
</td>
</tr>
<tr class="row-even"><td>ITC</td>
<td>integer</td>
<td>New print-out interval.</td>
</tr>
<tr class="row-odd"><td>DIT4</td>
<td>real</td>
<td>Maximum time step size for next time interval (days).</td>
</tr>
</tbody>
</table>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">time</span></code>. In this example, the initial time step size is
30 days, the final simulation time is 3650 days, the number of time steps allowed is
20, nodal information is printed out for every 5th time step, the simulation starts
in the 10th month of 1989, and the initial time of simulation is assigned a value of
0.</p>
<p>The time step multiplier is changed after 1 day, and the new time step multiplier
is 1.2, backward Euler is used from this time on and the printout interval is every
10th time step. The maximum time step size for the next interval is omitted so the
default value entered in the <code class="docutils literal notranslate"><span class="pre">ctrl</span></code> macro will be used.</p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="21%" />
<col width="13%" />
<col width="11%" />
<col width="16%" />
<col width="11%" />
<col width="13%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>time</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>30.0</td>
<td>3650.0</td>
<td>20</td>
<td>5</td>
<td>1989</td>
<td>10</td>
<td>0.0</td>
</tr>
<tr class="row-odd"><td>1.0</td>
<td>-1.2</td>
<td>1.0</td>
<td>10</td>
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
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
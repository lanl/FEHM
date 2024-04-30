---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">head</span></code><a class="headerlink" href="#head" title="Permalink to this headline">¶</a></h1>
<p>Hydraulic head values are used for input and output instead of pressures. Use of this macro enables the Boussinesq approximation (bous macro) and isothermal air-water two-phase simulation (airwater macro) automatically. It affects the pres and flow macros by requiring head information where pressure values were previously required. The default is to have no input associated with this macro. However, an optional head increment can be given after the head macro keyword. This value will be added to all input head values to ensure a single phase fluid state. Note that the value will be subtracted before output is written.</p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="8%" />
<col width="9%" />
<col width="66%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>HEAD0</td>
<td>real</td>
<td><ol class="first last arabic simple" start="0">
<li></li>
</ol>
</td>
<td>An incremental value that will be added to all input heads (m).</td>
</tr>
</tbody>
</table>
<p>The following is an example of head. In this example the optional head increment is included and a value of 1000. m is added to all input head values.</p>
<table border="1" class="docutils">
<colgroup>
<col width="46%" />
<col width="54%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>head</td>
<td><ol class="first last arabic simple" start="1000">
<li></li>
</ol>
</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
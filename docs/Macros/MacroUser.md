---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">user</span></code><a class="headerlink" href="#user" title="Permalink to this headline">¶</a></h1>
<p>Group 1 - KK</p>
<table border="1" class="docutils">
<colgroup>
<col width="4%" />
<col width="3%" />
<col width="93%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>KK</td>
<td>integer</td>
<td>Integer number passed to subroutine user for user defined input parameters. This user subroutine call differs from the one invoked in the control file in that whereas that subroutine is called at every time step, this one is called only at the beginning of the simulation to set parameters that do not change later in the simulation.</td>
</tr>
</tbody>
</table>
<p>The following is an example of user. In this example, the number 5 is passed to the user subroutine.</p>
<table border="1" class="docutils">
<colgroup>
<col width="100%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">user</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>5</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
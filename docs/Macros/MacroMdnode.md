---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">mdnode</span></code><a class="headerlink" href="#mdnode" title="Permalink to this headline">¶</a></h1>
<p>Enables extra connections to be made to nodes. This is useful for simulating wellbore connections, faults, and flow across internal boundaries.</p>
<ul class="simple">
<li>Group 1 -     NUM_MD, MAX_CON, IELIM, SX_MULT</li>
<li>Group 2 -     NODE, IPAR, NPAR (repeated NUM_MD times)</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="10%" />
<col width="10%" />
<col width="61%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Default</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>NUM_MD</td>
<td>integer</td>
<td>0</td>
<td>Number of new connections to be entered.</td>
</tr>
<tr class="row-odd"><td>MDMAX</td>
<td>integer</td>
<td>0</td>
<td>Maximum number of new connections to a given node.
This does not include old connections. Thus,
if a node was already connected to 5 neighboring
nodes and two new connections were added to this
node in this macro statement and this was the
maximum number of connections added in this
macro statement, then MDMAX = 2.</td>
</tr>
<tr class="row-even"><td>I_ELIM</td>
<td>integer</td>
<td>0</td>
<td>IF I_ELIM Š 0, then no action. IF I_ELIM &lt; 0,
then nodal connections are eliminated as needed
if redundant.</td>
</tr>
<tr class="row-odd"><td>SX_MULT</td>
<td>real*8</td>
<td>1.0</td>
<td>Multiplier for equilibrium conditions.</td>
</tr>
<tr class="row-even"><td>NODE</td>
<td>integer</td>
<td>0</td>
<td>Node to which new connection is established.</td>
</tr>
<tr class="row-odd"><td>IPAR</td>
<td>integer</td>
<td>0</td>
<td>IPAR is not used at present. Its value is ignored.
However the entered number must be an integer.</td>
</tr>
<tr class="row-even"><td>NPAR</td>
<td>integer</td>
<td>0</td>
<td>NPAR is the new connected node. If NPAR = NODE,
no new connection is established.</td>
</tr>
</tbody>
</table>
<p>The following are examples of mdnode. In the first example (top), 3 new connections
are specified, node 10 is connected to node 15, node 100 is connected to node
106, and node 10 is connected to node 320. A maximum of 2 new connections are
allowed per node. The multiplier for equilibrium conditions is set to 10. In the
second example (bottom), 4 new connections are specified, node 1 is connected to
node 16, node 2 is connected to node 1, node 4 is connected to node 1 and node
10 is connected to node 203. A maximum of 3 new connections are allowed per node.
The multiplier for equilibrium conditions is set to 100.</p>
<table border="1" class="docutils">
<colgroup>
<col width="40%" />
<col width="15%" />
<col width="25%" />
<col width="20%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>mdnode</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>3</td>
<td>2</td>
<td>0</td>
<td>10</td>
</tr>
<tr class="row-odd"><td>10</td>
<td>0</td>
<td>15</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>100</td>
<td>0</td>
<td>106</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>10</td>
<td>0</td>
<td>320</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
<table border="1" class="docutils">
<colgroup>
<col width="38%" />
<col width="14%" />
<col width="24%" />
<col width="24%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>mdnode</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>4</td>
<td>3</td>
<td>0</td>
<td>100</td>
</tr>
<tr class="row-odd"><td>1</td>
<td>0</td>
<td>16</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>2</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>4</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>10</td>
<td>0</td>
<td>203</td>
<td>&#160;</td>
</tr>
</tbody>
</table>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
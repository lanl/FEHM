---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">itfc</span></code><a class="headerlink" href="#itfc" title="Permalink to this headline">¶</a></h1>
<p>Data to define flow and transport parameters at interfaces between pairs of zones.</p>
<ul class="simple">
<li>Group 1 -     ZONE_PAIR(I,1), ZONE_PAIR(I,2), RED_FACTOR(I)- an arbitrary number of lines of input, terminated by a blank line.</li>
<li>Group 2 -     (FILTER_FLAG(J), J= 1,NSPECI)</li>
<li>Group 3 -     ZONEC_PAIR(K,1), ZONEC_PAIR(K,2), FTN_FACTOR(K)- an arbitrary number of lines of input, terminated by a blank line.<ul>
<li>KEYWORD ‘file’</li>
<li>SFILENAME</li>
<li>ITFCPORSIZE(I), ITFCPROBSIZE(I)- an arbitrary number of lines of input, terminated by a blank line.</li>
</ul>
</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="16%" />
<col width="66%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>ZONE_PAIR</td>
<td>integer</td>
<td>Zone number for the zones for which the code
identifies the interface connections when
applying the permeability reduction factor.</td>
</tr>
<tr class="row-odd"><td>RED_FACTOR</td>
<td>real</td>
<td>Reduction factor multiplying the harmonically
weighted saturated permeability for all connections
at the interface identified by ZONE_PAIR</td>
</tr>
<tr class="row-even"><td>FILTER_FLAG</td>
<td>integer</td>
<td>FEHM has a provision to apply transport mechanisms
for size exclusion or filtration at interfaces
defined in the itfc macro. These provisions can
be used to simulate conditions in which, for
example, abrupt changes in properties occur at
interfaces, or hydrologic conditions not
explicitly incorporated in a model (a thin clay
layer, for example) are thought to be present
that affect transport across the interface.
The means for specifying these interface transport
conditions is the itfc macro. Thus, this parameter
is a flag used to distinguish whether the size
exclusion or filtration is to be implemented
(a value 1) or not (a value 0) for each species
identified in the trac, ptrk, or mptr macros.
The default value is 0. See the definition of
FTN_FACTOR below for details on how to invoke
the size exclusion or filtration model.</td>
</tr>
<tr class="row-odd"><td>ZONEC_PAIR</td>
<td>integer</td>
<td>Zone number for the zones for which the code
identifies the interface connections when applying
the transport filtration or size exclusion factors.</td>
</tr>
<tr class="row-even"><td>FTN_FACTOR</td>
<td>real</td>
<td>Filtration or size exclusion factor applied for
all connections at the interface identified by
ZONEC_PAIR. For the trac macro, a size exclusion
model is implemented, where FTN_FACTOR = 0 (size
exclusion) or 1 (no exclusion) are options. For
ptrk or mptr, a filtration model is implemented,
where the parameter is the probability of the
particle passing through the interface (if 0,
filtration is guaranteed; if 1, there is no filtration).
For the particle tracking model, FTN_FACTOR &lt; 0
denotes that the pore size distribution is being
used. This option is used with the particle size
distribution option in ptrk and mptr, so that
each particle is assigned a size. The cumulative
pore size distribution is then used as a
probability distribution function, and when a
particle encounters the interface, a pore size is
randomly selected from the distribution. If the
particle is larger than the pore, it is filtered.
Note that filtered particles remain at that location
in the model and are no longer transported.</td>
</tr>
<tr class="row-odd"><td>KEYWORD</td>
<td>character*4</td>
<td>Optional keyword ‘file’ designating that the pore
size distribution information is being input in
a separate file. This input is entered only for
interfaces in which FTN_FACTOR &lt; 0 is used.</td>
</tr>
<tr class="row-even"><td>SFILENAME</td>
<td>character*80</td>
<td>Optional file name containing the pore size distribution
table. This input is entered only for interfaces
in which FTN_FACTOR &lt; 0 is used.</td>
</tr>
<tr class="row-odd"><td>ITFCPORSIZE</td>
<td>real</td>
<td>Pore size for this entry of the pore size distribution
table (paired with a value of ITFCPROBSIZE).
An arbitrary number of entries can be input,
terminated with a blank line. These entries are
located in the file SFILENAME if specified, or
in the itfc input file if the alternate input
file is not used. The code decides if particles
are irreversibly filtered by comparing the
particle size to the randomly selected pore size.
This input is entered only for interfaces in which
FTN_FACTOR &lt; 0 is used.</td>
</tr>
<tr class="row-even"><td>ITFCPROBSIZE</td>
<td>real</td>
<td>Cumulative probability for the distribution of
pore sizes (paired with a value of ITFCPORSIZE).
See description of ITFCPORSIZE above for details.
The final entry of the table must have
ITFCPROBSIZE = 1, since the distribution is assumed
to be normalized to unity. This input is entered
only for interfaces in which FTN_FACTOR &lt; 0 is used.</td>
</tr>
</tbody>
</table>
<p>Note that data for each numbered group must be input. The other input is optional.
If filtration is not implemented for any species, a single blank line is input
for Groups 2 and 3, signaling the end of itfc input.</p>
<p>The following is an example of itfc. In this example, the permeability reduction
factor of 0.1 is applied to all node connections at the interface between zones
6 and 10, or 6 and 11.</p>
<table border="1" class="docutils">
<colgroup>
<col width="40%" />
<col width="27%" />
<col width="33%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td>itfc</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>6</td>
<td>10</td>
<td>0.1</td>
</tr>
<tr class="row-odd"><td>6</td>
<td>11</td>
<td>0.1</td>
</tr>
<tr class="row-even"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
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
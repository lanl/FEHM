---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">gdkm</span></code><a class="headerlink" href="#gdkm" title="Permalink to this headline">¶</a></h1>
<p>Generalized dual permeability model.</p>
<p>The input structure for the gdkm module is the same as the gdpm input in FEHM.</p>
<ul class="simple">
<li>Group 1 -     GDKM_FLAG, NGDKMNODES</li>
<li>Group 2 -     NGDKM_LAYERS(I), VFRAC_PRIMARY(I), (GDKM_X(I,J), J=1,NGDKM_LAYERS(I))</li>
</ul>
<p>An arbitrary numbers of lines of input, terminated by a blank line.</p>
<p>Group 3 -       JA, JB, JC, IGDKM (JA, JB, JC - defined on <a class="reference external" href="Macro20058.html">JA, JB, JC, PROP1, PROP2, …</a>)</p>
<table border="1" class="docutils">
<colgroup>
<col width="5%" />
<col width="3%" />
<col width="91%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>GDKM_FLAG</td>
<td>integer</td>
<td>Flag to denote that the GDKM model option is being invoked. The default is 0 if GDKM is not being used. At present the only model allowed is model 11. This is a parallel fracture type model.</td>
</tr>
<tr class="row-odd"><td>NGDKMNODES</td>
<td>integer</td>
<td>Total number of matrix nodes present in the simulation. The GDKM gridblocks, in contrast to GDPM gridblocks, are restricted to a single secondary node in a GDKM gridblock. Thus NGDPMNODES is equal to the number of gridblocks that have a GDKM model applied to them.</td>
</tr>
<tr class="row-even"><td>NGDKM_LAYERS</td>
<td>integer</td>
<td>The number of matrix nodes specified for this model number. This is always 1 for GDKM grid blocks. All primary nodes assigned to this model number (using the IGDPM input below) will have 1 matrix node.</td>
</tr>
<tr class="row-odd"><td>VFRAC_PRIMARY</td>
<td>real</td>
<td>The fraction of the total control volume that is assigned to the primary porosity. Then, 1-VFRAC_PRIMARY is the fraction of the control volume that is assigned to the secondary porosity node.</td>
</tr>
<tr class="row-even"><td>GDKM_X</td>
<td>real</td>
<td>The matrix discretization distance for the matrix node associated with this model (units of meters). For the one secondary node allowed in the GDKM formulation, the average distance to the secondary node from the primary node.</td>
</tr>
<tr class="row-odd"><td>IGDKM</td>
<td>integer</td>
<td>Model number for parameters defined in group 2. These values are assigned only for the primary nodes. The default is 0, which denotes that there are no dual permeability nodes at that primary nodes.</td>
</tr>
</tbody>
</table>
<p>Input Variable Format Description</p>
  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
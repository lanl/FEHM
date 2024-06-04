---
layout : page_macros
hero_height: is-hidden
---

<h1><code class="docutils literal notranslate"><span class="pre">gdpm</span></code><a class="headerlink" href="#gdpm" title="Permalink to this headline">¶</a></h1>
<p>Data to define the parameters in the Generalized Dual Porosity model formulation.</p>
<ul class="simple">
<li>Group 1 -     GDPM_FLAG, NGDPMNODES</li>
<li>Group 2 -     NGDPM_LAYERS(I), VFRAC_PRIMARY(I), (GDPM_X(I,J), J=1,NGDPM_LAYERS(I))- an arbitrary numbers of lines of input, terminated by a blank line.</li>
<li>Group 3 -     JA, JB, JC, IGDPM (JA, JB, JC - defined on page 27)</li>
</ul>
<table border="1" class="docutils">
<colgroup>
<col width="4%" />
<col width="2%" />
<col width="94%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Input Variable</th>
<th class="head">Format</th>
<th class="head">Description</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>GDPM_FLAG</td>
<td>integer</td>
<td>Flag to denote that the GDPM model option is being invoked. The default is 0 if GDPM is not being used. If 1, matrix node geometry is parallel fractures; if 2, matrix node geometry is spherical, with the fractured medium residing at the exterior of an idealized spherical block, and transport occurs into the block.</td>
</tr>
<tr class="row-odd"><td>NGDPMNODES</td>
<td>integer</td>
<td>Total number of matrix nodes present in the simulation. Since this number may not be known at runtime, the code may be run once with a placeholder value for NGDPMNODES. If the number is incorrect, the code will report the appropriate value and stop. This value can then be entered and the simulation will proceed when the code is rerun.</td>
</tr>
<tr class="row-even"><td>NGDPM_LAYERS</td>
<td>integer</td>
<td>The number of matrix nodes specified for this model number. All primary nodes assigned to this model number (using the IGDPM input below) will have NGDPM_LAYERS matrix nodes.</td>
</tr>
<tr class="row-odd"><td>VFRAC_PRIMARY</td>
<td>real</td>
<td>The fraction of the total control volume that is assigned to the primary porosity. Then, 1-VFRAC_PRIMARY is the fraction of the control volume that is divided among the dual porosity nodes.</td>
</tr>
<tr class="row-even"><td>GDPM_X</td>
<td>real</td>
<td>The matrix discretization distances for the matrix nodes associated with this model (units of meters). Grid points are placed at these values to discretize each matrix block. There must be NGDPM_LAYERS values, entered in ascending order. For the parallel plate geometry, the final value is the distance to the centerline between the fractures, and for the spherical geometry, the final value is the radius of the sphere.</td>
</tr>
<tr class="row-odd"><td>IGDPM</td>
<td>integer</td>
<td>Model number for parameters defined in group 2. These values are assigned only for the primary nodes. The default is 0, which denotes that there are no dual porosity nodes at that primary node.</td>
</tr>
</tbody>
</table>
<p>Based on the input in this macro, the code internally assigns node numbers, finite element coefficients,
and reconstructs the connectivity array for the grid. The original nodes in the grid (the primary nodes) retain the node numbers 1 to
NEQ_PRIMARY, where NEQ_PRIMARY is the number of nodes assigned in the macro coor. The matrix nodes are assigned numbers from NEQ_PRIMARY + 1 to NEQ_PRIMARY + NGDPMNODES.
To assign matrix node numbers, the code loops through each primary node, and if GDPM nodes are specified, assigns the node numbers in increasing order within the matrix block of that primary node.</p>
<div class="admonition note">
<p class="first admonition-title">Note</p>
<p class="last">The user is responsible for assigning rock, hydrologic, and transport properties
for the matrix nodes as well as the primary nodes. For input using zones, this process is facilitated with
the convention that zone numbers associated with matrix nodes are set to ZONE_DPADD + the zone number
for the corresponding fracture node. This convention is overwritten for any matrix node for which the user assigns a zone number
using the ‘nnum’ option in the macro zone.</p>
</div>
<p>For output, the code can report time-varying values in the “.out”, “.his”, and “.trc” files for both primary and matrix nodes, but fields written for the entire grid (for example, in the AVS output using the macro cont) are output only for the primary nodes.</p>
<p>The following is an example of <code class="docutils literal notranslate"><span class="pre">gdpm</span></code>. In this example the matrix node geometry is parallel to the fractures and there are 1479 matrix nodes distributed in 29 layers. A single model is defined which is applied to the entire problem domain.</p>
<table border="1" class="docutils">
<colgroup>
<col width="14%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">gdpm</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
<th class="head">&#160;</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>1</td>
<td>1479</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>29</td>
<td>.0001</td>
<td>.001</td>
<td>.002</td>
<td>.003</td>
<td>.004</td>
<td>.006</td>
<td>.009</td>
</tr>
<tr class="row-even"><td>.019</td>
<td>.02901</td>
<td>.03901</td>
<td>.04901</td>
<td>.05901</td>
<td>.09902</td>
<td>.19904</td>
<td>.29906</td>
</tr>
<tr class="row-odd"><td>.39908</td>
<td>.49910</td>
<td>.59912</td>
<td>.69914</td>
<td>.79916</td>
<td>.89918</td>
<td>.99920</td>
<td>1.4993</td>
</tr>
<tr class="row-even"><td>1.9994</td>
<td>2.4995</td>
<td>2.9996</td>
<td>3.4997</td>
<td>3.9998</td>
<td>4.4999</td>
<td>5.0000</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-even"><td>1</td>
<td>0</td>
<td>0</td>
<td>1</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
<td>&#160;</td>
</tr>
<tr class="row-odd"><td>&#160;</td>
<td>&#160;</td>
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
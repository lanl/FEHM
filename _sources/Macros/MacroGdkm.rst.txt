========
``gdkm``
========

Generalized dual permeability model.

The input structure for the gdkm module is the same as the gdpm input in FEHM.

* Group 1 -	GDKM_FLAG, NGDKMNODES

* Group 2 -	NGDKM_LAYERS(I), VFRAC_PRIMARY(I), (GDKM_X(I,J), J=1,NGDKM_LAYERS(I))

An arbitrary numbers of lines of input, terminated by a blank line.

Group 3 -	JA, JB, JC, IGDKM (JA, JB, JC - defined on `JA, JB, JC, PROP1, PROP2, . . . <Macro20058.html>`_)

+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                                                                                                                                                                                                              |
+================+=========+==========================================================================================================================================================================================================================================================================+
| GDKM_FLAG      | integer | Flag to denote that the GDKM model option is being invoked. The default is 0 if GDKM is not being used. At present the only model allowed is model 11. This is a parallel fracture type model.                                                                           |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NGDKMNODES     | integer | Total number of matrix nodes present in the simulation. The GDKM gridblocks, in contrast to GDPM gridblocks, are restricted to a single secondary node in a GDKM gridblock. Thus NGDPMNODES is equal to the number of gridblocks that have a GDKM model applied to them. |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NGDKM_LAYERS   | integer | The number of matrix nodes specified for this model number. This is always 1 for GDKM grid blocks. All primary nodes assigned to this model number (using the IGDPM input below) will have 1 matrix node.                                                                |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| VFRAC_PRIMARY  | real    | The fraction of the total control volume that is assigned to the primary porosity. Then, 1-VFRAC_PRIMARY is the fraction of the control volume that is assigned to the secondary porosity node.                                                                          |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| GDKM_X         | real    | The matrix discretization distance for the matrix node associated with this model (units of meters). For the one secondary node allowed in the GDKM formulation, the average distance to the secondary node from the primary node.                                       |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IGDKM          | integer | Model number for parameters defined in group 2. These values are assigned only for the primary nodes. The default is 0, which denotes that there are no dual permeability nodes at that primary nodes.                                                                   |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Input Variable Format Description


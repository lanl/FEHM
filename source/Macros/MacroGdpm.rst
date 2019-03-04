========
``gdpm``
========

Data to define the parameters in the Generalized Dual Porosity model formulation.

* Group 1 -	GDPM_FLAG, NGDPMNODES

* Group 2 -	NGDPM_LAYERS(I), VFRAC_PRIMARY(I), (GDPM_X(I,J), J=1,NGDPM_LAYERS(I))- an arbitrary numbers of lines of input, terminated by a blank line.

* Group 3 -	JA, JB, JC, IGDPM (JA, JB, JC - defined on page 27)

+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                                                                                                                                                                                                                                                                                                                                                                          |
+================+=========+======================================================================================================================================================================================================================================================================================================================================================================================================================================+
| GDPM_FLAG      | integer | Flag to denote that the GDPM model option is being invoked. The default is 0 if GDPM is not being used. If 1, matrix node geometry is parallel fractures; if 2, matrix node geometry is spherical, with the fractured medium residing at the exterior of an idealized spherical block, and transport occurs into the block.                                                                                                          |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NGDPMNODES     | integer | Total number of matrix nodes present in the simulation. Since this number may not be known at runtime, the code may be run once with a placeholder value for NGDPMNODES. If the number is incorrect, the code will report the appropriate value and stop. This value can then be entered and the simulation will proceed when the code is rerun.                                                                                     |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NGDPM_LAYERS   | integer | The number of matrix nodes specified for this model number. All primary nodes assigned to this model number (using the IGDPM input below) will have NGDPM_LAYERS matrix nodes.                                                                                                                                                                                                                                                       |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| VFRAC_PRIMARY  | real    | The fraction of the total control volume that is assigned to the primary porosity. Then, 1-VFRAC_PRIMARY is the fraction of the control volume that is divided among the dual porosity nodes.                                                                                                                                                                                                                                        |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| GDPM_X         | real    | The matrix discretization distances for the matrix nodes associated with this model (units of meters). Grid points are placed at these values to discretize each matrix block. There must be NGDPM_LAYERS values, entered in ascending order. For the parallel plate geometry, the final value is the distance to the centerline between the fractures, and for the spherical geometry, the final value is the radius of the sphere. |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IGDPM          | integer | Model number for parameters defined in group 2. These values are assigned only for the primary nodes. The default is 0, which denotes that there are no dual porosity nodes at that primary node.                                                                                                                                                                                                                                    |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Based on the input in this macro, the code internally assigns node numbers, finite element coefficients,
and reconstructs the connectivity array for the grid. The original nodes in the grid (the primary nodes) retain the node numbers 1 to
NEQ_PRIMARY, where NEQ_PRIMARY is the number of nodes assigned in the macro coor. The matrix nodes are assigned numbers from NEQ_PRIMARY + 1 to NEQ_PRIMARY + NGDPMNODES.
To assign matrix node numbers, the code loops through each primary node, and if GDPM nodes are specified, assigns the node numbers in increasing order within the matrix block of that primary node.

.. note::

  The user is responsible for assigning rock, hydrologic, and transport properties 
  for the matrix nodes as well as the primary nodes. For input using zones, this process is facilitated with
  the convention that zone numbers associated with matrix nodes are set to ZONE_DPADD + the zone number 
  for the corresponding fracture node. This convention is overwritten for any matrix node for which the user assigns a zone number 
  using the ‘nnum' option in the macro zone.

For output, the code can report time-varying values in the ".out", ".his", and ".trc" files for both primary and matrix nodes, but fields written for the entire grid (for example, in the AVS output using the macro cont) are output only for the primary nodes.

The following is an example of ``gdpm``. In this example the matrix node geometry is parallel to the fractures and there are 1479 matrix nodes distributed in 29 layers. A single model is defined which is applied to the entire problem domain.

+---------+--------+--------+--------+--------+--------+--------+--------+
| gdpm    |        |        |        |        |        |        |        |
+=========+========+========+========+========+========+========+========+
| 1       | 1479   |        |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+--------+--------+
| 29      | .0001  | .001   | .002   | .003   | .004   | .006   | .009   |
+---------+--------+--------+--------+--------+--------+--------+--------+
| .019    | .02901 | .03901 | .04901 | .05901 | .09902 | .19904 | .29906 |
+---------+--------+--------+--------+--------+--------+--------+--------+
| .39908  | .49910 | .59912 | .69914 | .79916 | .89918 | .99920 | 1.4993 |
+---------+--------+--------+--------+--------+--------+--------+--------+
| 1.9994  | 2.4995 | 2.9996 | 3.4997 | 3.9998 | 4.4999 | 5.0000 |        |
+---------+--------+--------+--------+--------+--------+--------+--------+
|         |        |        |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+--------+--------+
| 1       | 0      | 0      | 1      |        |        |        |        |
+---------+--------+--------+--------+--------+--------+--------+--------+
|         |        |        |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+--------+--------+

========
``dual``
========

Dual porosity formulation. There are three sets of parameter values at any nodal
position, for which property values must be defined. Nodes 1 to N (see macro ``coor``
for definition of N) represent the fracture nodes, nodes N + 1 to 2N the first
matrix material, and nodes 2N + 1 to 3N the second matrix material.
When zones are used with the ``dual`` macro, additional zones are automatically
generated.

See instructions for the macro ``zone`` for a more detailed description.
The ``dual`` parameters are only defined for the first N nodes.

* Group 1 - IDUALP
* `Group 2 - JA, JB, JC, VOLFD1 <InputData.html#JA>`_
* `Group 3 - JA, JB, JC, VOLFD2 <InputData.html#JA>`_
* `Group 4 - JA, JB, JC, APUVD <InputData.html#JA>`_

+----------------+---------+---------+-------------------------------------------------------------------+
| Input Variable | Format  | Default | Description                                                       |
+----------------+---------+---------+-------------------------------------------------------------------+
| IDUALP         | integer |         | | Solution descriptor for dual porosity solution.                 |
|                |         |         | | IDUALP = 0, information is read but not used.                   |
|                |         |         | | IDUALP ≠ 0, dual porosity solution is implemented               |
|                |         |         | | For the special case of IDUALP = 2, the                         |
|                |         |         |   permeabilities and conductivities are scaled                    |
|                |         |         |   by the volume fraction, i.e., :math:`k = vf * k`.               |
+----------------+---------+---------+-------------------------------------------------------------------+
| VOLFD1         | real    | 0.001   | Volume fraction for fracture portion of the continuum.            |
+----------------+---------+---------+-------------------------------------------------------------------+
| VOLFD2         | real    | 0.5     | Volume fraction for the first matrix portion of the continuum.    |
+----------------+---------+---------+-------------------------------------------------------------------+
| APUVD          | real    | 5.      | Length scale for the matrix nodes (m).                            |
+----------------+---------+---------+-------------------------------------------------------------------+


The volume fractions VOLFD1 and VOLFD2 are related to the total volume by

:math:`VOLFD1 + VOLFD2 + VOLFD3 = 1.0`

where VOLFD3 is the volume fraction of the second matrix node.
If permeability model IRLP = 4 is selected in control statement ``rlp``,
VOLFD1 is calculated from RP15 (fracture porosity) in that control statement.

The following is an example of ``dual``. In this example, the dual porosity
solution is implemented for all nodes from 1 through 140. The volume fraction
for the fracture is 0.006711409, the volume fraction for the first matrix portion
is 0.335570470, and the length scale for the matrix nodes is 0.1 m. 

+------+-----+---+-------------+
| dual |     |   |             |
+------+-----+---+-------------+
| 1    |     |   |             |
+------+-----+---+-------------+
| 1    | 140 | 1 | 0.006711409 |
+------+-----+---+-------------+
|      |     |   |             |
+------+-----+---+-------------+
| 1    | 140 | 1 | 0.335570470 |
+------+-----+---+-------------+
|      |     |   |             |
+------+-----+---+-------------+
| 1    | 140 | 1 | 0.10        |
+------+-----+---+-------------+
|      |     |   |             |
+------+-----+---+-------------+

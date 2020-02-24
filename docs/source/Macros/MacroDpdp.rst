========
``dpdp``
========

Double porosity / double permeability formulation. There are two sets of parameter
values at any nodal position, for which property values must be defined.
Nodes 1 to N (see macro coor for definition of N) represent the fracture nodes
and nodes N + 1 to 2N the matrix material. When zones are used with the dpdp macro,
additional zones are automatically generated. See instructions for the macro zone
for a more detailed description. The dpdp parameters are only defined for the first N nodes.

* Group 1 - IDPDP
* Group 2 - `JA, JB, JC, VOLFD1  <InputData.html#JA>`_
* Group 3 - `JA, JB, JC, APUV1 <InputData.html#JA>`_ 

+----------------+---------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Default | Description                                                                                                                                                 |
+----------------+---------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IDPDP          | integer |         | Solution descriptor for double porosity/double permeability solution. IDPDP = 0, information is read but not used. IDPDP ≠ 0, dpdp solution is implemented. |
+----------------+---------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| VOLFD1         | real    | 1.      | Volume fraction for fracture node.                                                                                                                          |
+----------------+---------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
| APUV1          | real    | 10.     | Half spacing between fractures (m). See TRANSFLAG in macros MPTR and PTRK.                                                                                  |
+----------------+---------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+

The volume fraction VOLFD1 is related to the total volume by

:math:`VOLFD1 + VOLLFD2 = 1.0`

where VOLFD2 is the volume fraction of the matrix node. If permeability model
``IRLP = 4`` is selected in control statement **rlp**,
VOLFD1 is calculated from RP15 (fracture porosity) in that control statement.

The following is an example of ``dpdp``. In this example,
the dual porosity/permeability solution is implemented for all nodes from
1 through 140. The fractional volume in the fractures (compared to the total volume)
is 0.005 and the length scale for matrix nodes is 0.1 m.

.. code::
   
   dpdp
   1
   1 140 1 0.005

   1 140 1 0.10

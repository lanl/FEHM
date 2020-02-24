========
``rock``
========

Assign rock density, specific heat and porosity.

+----------------+--------+-------------------------------------------------------------------------------------------------------+
| Input Variable | Format | Description                                                                                           |
+================+========+=======================================================================================================+
| DENRD          | real   | Rock density (kg/m3).                                                                                 |
+----------------+--------+-------------------------------------------------------------------------------------------------------+
| CPRD           | real   | Rock specific heat :math:`\left( \frac{MJ}{kg \cdot K} \right)`.                                      |
|                |        | If CPRD > 1 the code will assume the units are                                                        |
|                |        | :math:`\left( \frac{J}{kg \cdot K} \right)` and multiply by :math:`10^{-6}`.                          |
+----------------+--------+-------------------------------------------------------------------------------------------------------+
| PSD            | real   | | Porosity.                                                                                           |
|                |        | Special note on negative porosities. If the code encounters a negative porosity,                      |
|                |        | the node at which the negative porosity occurs is effectively removed from the model.                 |
|                |        | That is, the geometric connections from that node to other nodes in the model are                     |
|                |        | removed. The volume associated with the node acts as a barrier to flow. For input purposes,           |
|                |        | the node may still be assigned properties, though they will have no effect on the simulation results. |
+----------------+--------+-------------------------------------------------------------------------------------------------------+

The following is an example of ``rock``. In this example the nodes numbered
1 through 140 are assigned a rock density of 2563. kg/m^3,
a rock specific heat of 1010. J/(kg K) and a porosity of 0.35.

+------+-----+---+-------+-------+------+
| rock |     |   |       |       |      |
+------+-----+---+-------+-------+------+
| 1    | 140 | 1 | 2563. | 1010. | 0.35 |
+------+-----+---+-------+-------+------+


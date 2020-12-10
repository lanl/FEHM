========
``dvel``
========

Velocity between two nodes is output by choosing this control statement.
The input for this macro is identical to macro ``flxo``, except that
velocities instead of fluxes are calculated (see
`flxo (optional) <MacroFlxo.html>`_).

In the following example of ``dvel``, a single internode velocity is
calculated between nodes 101 and 102.

+------+-----+
| dvel |     |
+------+-----+
| 1    |     |
+------+-----+
| 101  | 102 |
+------+-----+
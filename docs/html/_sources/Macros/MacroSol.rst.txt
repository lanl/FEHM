=======
``sol``
=======

* Group 1 - NTT, INTG

+----------------+---------+---------+-----------------------------------------------------------+
| Input Variable | Format  | Default | Description                                               |
+================+=========+=========+===========================================================+
| NTT            | integer | 1       | | Parameter that defines the type of solution required    |
|                |         |         | | NTT > 0 coupled solution                                |
|                |         |         | | NTT ≤ 0 heat transfer only solution                     |
+----------------+---------+---------+-----------------------------------------------------------+
| INTG           | integer | -1      | | Parameter that defines element integration type         |
|                |         |         | | INTG ≤ 0 Lobatto (node point) quadrature is used,       |
|                |         |         |   recommended for heat and mass problems without stress.  |
|                |         |         | | INTG > 0 Gauss quadrature is used, recommended for      |
|                |         |         |   problems requiring a stress solution.                   |
+----------------+---------+---------+-----------------------------------------------------------+

The following is an example of ``sol``. In this example, a coupled heat-mass
solution using Lobatto quadrature is specified.

+-----+-----+
| sol |     |
+-----+-----+
|  1  |  -1 |
+-----+-----+
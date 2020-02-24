========
``rich``
========

Invokes Richard's equation solution for unsaturated-saturated flow. A single
phase approach that neglects air phase flow and assumes the movement of water
is independent of air flow and pressure.

Uses variable switching (Pressure, Saturation).

* Group 1 - STRD_RICH, TOL_PHASE, PCHNG, SCHNG

+----------------+--------+---------------------------------------------+
| Input Variable | Format | Description                                 |
+================+========+=============================================+
| STRD_RICH      | real   | Newton-raphson relaxation factor            |
+----------------+--------+---------------------------------------------+
| TOL_PHASE      | real   | Tolerance for full saturation               |
+----------------+--------+---------------------------------------------+
| PCHNG          | real   | Pressure adjustment after variable switch   |
+----------------+--------+---------------------------------------------+
| SCHNG          | real   | Saturation adjustment after variable switch |
+----------------+--------+---------------------------------------------+

The following is an example of rich. 

+------+-------+-------+-------+
| rich |       |       |       |
+------+-------+-------+-------+
| 0.95 | 1.e-5 | 1.e-3 | 1.e-3 |
+------+-------+-------+-------+


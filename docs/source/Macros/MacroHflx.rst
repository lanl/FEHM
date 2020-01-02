========
``hflx``
========

Heat flow input.

A negative heat flow indicates heat flow into the reservoir. 

+----------------+--------+---------+----------------------------------------------------------------------------+
| Input Variable | Format | Default | Description                                                                |
+================+========+=========+============================================================================+
| QFLUX          | real   | 0.      | | If QFLXM = 0, then QFLUX is the heat flow (MW). If QFLXM ≠ 0,            |
|                |        |         |   then QFLUX is a temperature (oC) and the heat flow is                    |
|                |        |         |   calculated according to the formula:                                     |
|                |        |         | | :math:`Q_H = QFLXM \cdot (T-QFLUX)` (MW).                                |
+----------------+--------+---------+----------------------------------------------------------------------------+
| QFLXM          | real   | 0.      | | If QFLXM > 0, multiplier for heat flow equation given in QFLUX           |
|                |        |         |   description (MW/oC). This must be large for large temperature gradients, |
|                |        |         |   or when a constant temperature must be maintained.                       |
|                |        |         | | If QFLXM < 0, then QFLUX is interpreted as a fixed saturation and        |
|                |        |         | | :math:`Q_H = ABS(QFLXM) \cdot (S_l - QFLUX)` (MW).                       |
+----------------+--------+---------+----------------------------------------------------------------------------+

The following is an example of ``hflx``. In this example, at each node from 401 to 410,
a heat flow of 0.001 MW is being injected into the model.

+------+-----+---+--------+-----+
| hflx |     |   |        |     |
+------+-----+---+--------+-----+
| 401  | 410 | 1 | -0.001 | 0.0 |
+------+-----+---+--------+-----+
|      |     |   |        |     |
+------+-----+---+--------+-----+
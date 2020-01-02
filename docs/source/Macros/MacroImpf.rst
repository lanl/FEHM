========
``impf``
========

Time step control based on maximum allowed variable change.

* Group 1 -	DELPT, DELTT, DELST, DELAT

+----------------+--------+------------------------------------------------------------------------------------+
| Input Variable | Format | Description                                                                        |
+================+========+====================================================================================+
| DELPT          | real   | Maximum allowable pressure change for which time step will be increased. (MPa)     |
+----------------+--------+------------------------------------------------------------------------------------+
| DELTT          | real   | Maximum allowable temperature change for which time step will be increased. (oC)   |
+----------------+--------+------------------------------------------------------------------------------------+
| DELST          | real   | Maximum allowable saturation change for which time step will be increased.         |
+----------------+--------+------------------------------------------------------------------------------------+
| DELAT          | real   | Maximum allowable air pressure change for which time step will be increased. (MPa) |
+----------------+--------+------------------------------------------------------------------------------------+

The following is an examples of ``impf``. In this example, pressure changes are
limited to 0.5 MPa, temperature changes to 20 oC, saturation changes to 0.1, and
air pressure changes to 0.05 MPa during a time step.

+------+------+-----+------+
| impf |      |     |      |
+------+------+-----+------+
| 0.5  | 20.0 | 0.1 | 0.05 |
+------+------+-----+------+

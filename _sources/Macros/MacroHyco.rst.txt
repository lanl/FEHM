``hyco``
=========

Hydraulic conductivity input.

.. note::

  Required if macro perm not used.

+----------------+--------+---------+--------------------------------------------------+
| Input Variable | Format | Default | Description                                      |
+================+========+=========+==================================================+
| PNX            | real   | 1.e-30  | Hydraulic conductivity in the x-direction (m/s). |
+----------------+--------+---------+--------------------------------------------------+
| PNY            | real   | 1.e-30  | Hydraulic conductivity in the y-direction (m/s). |
+----------------+--------+---------+--------------------------------------------------+
| PNZ            | real   | 1.e-30  | Hydraulic conductivity in the z-direction (m/s). |
+----------------+--------+---------+--------------------------------------------------+

The following is an example of the ``hyco`` macro. In this example, nodes 1 through 140 are specified to have hydraulic conductivities in the X, Y, and Z directions of 1.0e-5, 1.0e-5, and 0. m/s respectively. 

+------+-----+---+----------+----------+----------+
| hyco |     |   |          |          |          |
+------+-----+---+----------+----------+----------+
| 1    | 140 | 1 | 1.00e-05 | 1.00e-05 | 0.00e-00 |
+------+-----+---+----------+----------+----------+
|      |     |   |          |          |          |
+------+-----+---+----------+----------+----------+
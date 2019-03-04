
========
``fper``
========

Assign permeability scaling factors.

+----------------+--------+---------+-------------------------------------------------+
| Input Variable | Format | Default | Description                                     |
+================+========+=========+=================================================+
| SCALEX         | real   | 1.0     | Permeability scaling factor in the x-direction. |
+----------------+--------+---------+-------------------------------------------------+
| SCALEY         | real   | 1.0     | Permeability scaling factor in the y-direction. |
+----------------+--------+---------+-------------------------------------------------+
| SCALEZ         | real   | 1.0     | Permeability scaling factor in the z-direction. |
+----------------+--------+---------+-------------------------------------------------+

The following is an example of fper. In this example, the values of the permeability (defined in a previous perm macro) are multiplied by 1.0 in the X direction, 0.5 in the Y direction, and 0.1 in the Z direction. 

+------+------+---+-----+-----+-----+
| fper |      |   |     |     |     |
+======+======+===+=====+=====+=====+
|   1  |  140 | 1 | 1.0 | 0.5 | 0.1 |
+------+------+---+-----+-----+-----+

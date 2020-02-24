``cond``
========

Assign thermal conductivities of the rock.

.. note::

  ``cond`` is required for non-isothermal problems.

`Group 1 - JA, JB, JC, THXD, THYD, THZD <InputData.html#JA>`_

+----------------+--------+---------+-----------------------------------------+
| Input Variable | Format | Default | Description                             |
+================+========+=========+=========================================+
| THXD           | real   | 1.e-30  | Thermal conductivity in the x-direction |
+----------------+--------+---------+-----------------------------------------+
| THYD           | real   | 1.e-30  | Thermal conductivity in the y-direction |
+----------------+--------+---------+-----------------------------------------+
| THZD           | real   | 1.e-30  | Thermal conductivity in the z-direction |
+----------------+--------+---------+-----------------------------------------+

The following is an example of cond. In this example all the nodes numbered 1 through 140 have thermal conductivities of 1 in the X and Y directions, and 0 in the Z direction.

+------+-----+---+----------+-----------+----------+
| cond |     |   |          |           |          |
+======+=====+===+==========+===========+==========+ 
|  1   | 140 | 1 | 1.00e-00 | 1e.00e-00 | 0.00e-00 |
+------+-----+---+----------+-----------+----------+
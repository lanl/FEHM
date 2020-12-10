========
``perm``
========

Assign permeabilities of the rock. Permeabilities represent average values of a volume associated with a node. Note that using rlp models 4 or 6 to describe relative permeabilities causes these values to be overwritten. Permeabilties may be entered as log values.

.. note::

  Required if macro ``hyco`` not used.

+----------------+--------+---------+---------------------------------------+
| Input Variable | Format | Default | Description                           |
+================+========+=========+=======================================+
| PNXD           | real   | 1.e-30  | Permeability in the x-direction (m2). |
+----------------+--------+---------+---------------------------------------+
| PNYD           | real   | 1.e-30  | Permeability in the y-direction (m2). |
+----------------+--------+---------+---------------------------------------+
| PNZD           | real   | 1.e-30  | Permeability in the z-direction (m2). |
+----------------+--------+---------+---------------------------------------+



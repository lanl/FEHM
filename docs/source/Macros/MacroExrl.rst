========
``exrl``
========

Allows the user to choose linearized relative permeability. The linearized relative
permeability is formed using the nonlinear value of relative permeability at the
iteration number IEXRLP. After that iteration a relative permeability value based
on a Taylor series expansion in saturation is used.

* Group 1 - IEXRLP

+----------------+---------+--------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                    |
+================+=========+================================================================================+
| IEXRLP         | integer | If IEXRLP ≥ 1, then linearized relative permeability. Otherwise not enabled.   |
+----------------+---------+--------------------------------------------------------------------------------+

In the following example of ``exrl``, the user enables linearized relative permeability at the first iteration. 

+-------+
| exrl1 |
+-------+
| 1     |
+-------+

========
``bous``
========

Constant density and viscosity are used for the flow terms (``Boussinesq`` approximation).

.. note::

  Where the bous macro is used, the gravity term in the air phase is set to zero.

* Group 1 -	ICONS

+----------------+---------+----------------------------------------------------------------------+
| Input Variable | Format  | Description                                                          |
+================+=========+======================================================================+
| ICONS          | integer | | Parameter to enable constant density and viscosity for flow terms. |
|                |         | |  ICONS ≠ 0 enabled.                                                |
|                |         | |  ICONS = 0 disabled (default).                                     |
+----------------+---------+----------------------------------------------------------------------+

The following is an example of bous. In this example the ``Boussinesq`` approximation is enabled. 

+-------+
| bous  |
+=======+
|   1   |
+-------+
========
``pres``
========

The initial values defined in control statement **pres** supersede all others. Note that the term "saturated" refered to in IEOSD, is not the groundwater hydrology definition (volumetric fraction of pore void that is filled with water) used elsewhere in this document. Saturated here indicates that vapor and liquid phases exist simultaneously. The superheated region means that all pore space is filled with gas.

.. note::

  Required if macro ``init`` not used.

+----------------+---------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Default | Description                                                                                                                                                                                                                                              |
+================+=========+=========+==========================================================================================================================================================================================================================================================+
| PHRD           | real    | PEIN    | Initial pressure (MPa).                                                                                                                                                                                                                                  |
+----------------+---------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| TIND           | real    |         | Initial temperature (oC) if IEOSD = 1 or 3, Initial saturation if IEOSD = 2                                                                                                                                                                              |
+----------------+---------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IEOSD          | integer | 1       | Thermodynamic region parameter. IEOSD = 1, the compressed liquid regionIEOSD = 2, the saturation regionIEOSD = 3, the superheated region.If IEOSD < 0 then the code uses ABS (IEOSD) and fixes the values of PHRD and TIND to the values provided above. |
+----------------+---------+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

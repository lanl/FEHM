========
``cflx``
========

Total moles of liquid solute moving through a zone are output by choosing this control statement. Vapor solute molar flows are currently not available. When this macro is invoked, the following output is given at every solute time step:

* The sum of all solute source flow rates for each zone
* The sum of all solute sink rates for each zone
* The sum of all solute entering each zone
* The sum of all solute leaving each zone
* The net source/sink (boundary) solute flow for each zone

The following values can be included on the macro line to specify which solute flows should be output:

* 1 (source)
* 2 (sink)
* 3 (netin)
* 4 (netout)
* 5 (boundary)

The default is to output all values.

Zones must be defined using macro zone prior to using this macro.

Group 1 - ``CFLXZ``

Group 2 - ``ICFLXZ(I), I = 1, CFLXZ``

+----------------+---------+------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                  |
+================+=========+==============================================================================+
| CFLXZ          | integer | Number of zones for which output for solute flow through the zone iswritten. |
+----------------+---------+------------------------------------------------------------------------------+
| ICFLXZ         | integer | Zone numbers for which solute flow output is written (NFLXZ zones)           |
+----------------+---------+------------------------------------------------------------------------------+

The following is an example of cflx. In this example solute flow through zones 1, 6 and 10 will be output. 

+---------+---+---+----+
|       cflx           |
+=========+===+===+====+
| CFLXZ   | 3          |
+---------+---+---+----+
| ICFLXZ  | 1 | 6 | 10 |
+---------+---+---+----+
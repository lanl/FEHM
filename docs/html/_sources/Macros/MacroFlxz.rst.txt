========
``flxz``
========

Total flow through a zone is output by choosing this control statement. When this macro is invoked, the following output is given at every heat and mass transfer time step:

* The sum of all source flow rates for each zone

* The sum of all sink flow rates for each zone

* The net quantity passing through each zone

* The net source/sink quantity for each zone

* The following keywords can be included on the macro line to specify which flows should be output: liquid mass, vapor mass, thermal energy. The default is to output any active quantity in the simulation.

Zones must be defined using macro zone prior to using this macro.

* Group 1 -	NFLXZ

* Group 2 -	IFLXZ(I), I = 1, NFLXZ

+----------------+---------+-----------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                       |
+================+=========+===================================================================================+
| NFLXZ          | integer | Number of zones for which output for water mass flow through the zone is written. |
+----------------+---------+-----------------------------------------------------------------------------------+
| IFLXZ          | integer | Zone numbers for which water mass flow output is written (NFLXZ zones)            |
+----------------+---------+-----------------------------------------------------------------------------------+

The following is an example of ``flxz``. In this example water mass flow through zones 1, 6 and 10 will be output. 

+------------+---+----+
| flxz water          |
+============+===+====+
| 3          |   |    |
+------------+---+----+
| 1          | 6 | 10 |
+------------+---+----+
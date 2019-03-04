========
``ngas``
========

Noncondensible gas transport. 

* Group 1 -	ICO2D

* Group 2 - JA, JB, JC, PCO2

* Group 3 - JA, JB, JC, CPNK

* Group 4 - JA, JB, JC, QCD

Note that all Group 2 values are entered first, followed by Group 3 values, followed by Group 4 values.

+----------------+---------+---------+-----------------------------------------------------------------------------+
| Input Variable | Format  | Default | Description                                                                 |
+================+=========+=========+=============================================================================+
| ICO2D          | integer | 3       | | Solution descriptor for noncondensible gas transport.                     |
|                |         |         | | ICO2D = 1, the 3 degree of freedom solution will be                       |
|                |         |         |   reduced to a 1 degree of freedom problem. (See macro                      |
|                |         |         |   ``iter``, the parameter ICOUPL is also set to 5 if                        |
|                |         |         |   ICO2D = 1.)                                                               |
|                |         |         | | ICO2D = 2, the 3 degree of freedom solution will be                       |
|                |         |         |   reduced to a 2 degree of freedom problem. (See                            |
|                |         |         |   macro ``iter``, the parameter ICOUPL is also set to                       |
|                |         |         |   5 if ICO2D = 2.) ICO2D = 3, full 3 degree of freedom.                     |
+----------------+---------+---------+-----------------------------------------------------------------------------+
| PCO2           | real    | 0.      | | Initial partial pressure of noncondensible gas. If                        |
|                |         |         |   PCO2 < 0 then ABS (PCO2) is interpreted as a temperature                  |
|                |         |         |   and the partial pressure of the noncondensible gas is                     |
|                |         |         |   calculated according to the formula:                                      |
|                |         |         | | :math:`PCO2 = P_T - P_{SAT}(T)`                                           |
|                |         |         | | where :math:`P_T` is the total pressure and :math:`P_{SAT}(T)`            |
|                |         |         |   is the water saturation pressure and is a function of temperature only.   |
+----------------+---------+---------+-----------------------------------------------------------------------------+
| CPNK           | real    | 0.      | | If CPNK ≤ 0, then ABS (CPNK) is the specified noncondensible              |
|                |         |         |   pressure and will be held at that value.                                  |
|                |         |         | | If CPNK > 0, then CPNK is the specified relative humidity                 |
|                |         |         |   and the saturation, :math:`S_l`, is calculated using the                  |
|                |         |         |   vapor pressure lowering formula and the capillary pressure formula:       |
|                |         |         | | :math:`Pcap(S_l) = \mathrm{ln}(h)\rho_l RT`                               |
|                |         |         | | where :math:`Pcap` is the capillary function, :math:`h` is                |
|                |         |         |   the humidity, :math:`R` is the gas constant, :math:`T`                    |
|                |         |         |   is the temperature, and :math:`\rho_l` is the liquid                      |
|                |         |         |   density. Once the formula is solved, :math:`S_l` is held                  |
|                |         |         |   constant. The humidity condition is only enabled for the                  |
|                |         |         |   van Genuchten capillary function model. See macro ``rlp``.                |
+----------------+---------+---------+-----------------------------------------------------------------------------+
| QCD            | real    | 0.      | Specified air source strength (kg/sec).                                     |
+----------------+---------+---------+-----------------------------------------------------------------------------+


The following is an example of ``ngas``. In this example, a full 3 degrees of
freedom solution is specified. The initial temperature at nodes 1 to 800 is 20
oC and the code is asked to calculate the initial noncondensible gas pressure.
There is no specified noncondensible gas source.


+------+-----+---+-----+
| ngas |     |   |     |
+------+-----+---+-----+
| 3    |     |   |     |
+------+-----+---+-----+
| 1    | 800 | 1 | -20 |
+------+-----+---+-----+
|      |     |   |     |
+------+-----+---+-----+
| 1    | 800 | 1 | 0.  |
+------+-----+---+-----+
|      |     |   |     |
+------+-----+---+-----+
| 1    | 800 | 1 | 0.  |
+------+-----+---+-----+
|      |     |   |     |
+------+-----+---+-----+
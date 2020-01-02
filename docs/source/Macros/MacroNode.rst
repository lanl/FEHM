========
``node``
========

Specify the node numbers for which detailed output is desired. In version 2.30 macro node has been modified to allow multiple instances of the macro to be used (results are cumulative) and to allow the definition of "output zones" which are used in conjunction with macro hist. Only a single input format / keyword can be used for each instance of the node macro.

* Group 1 -	M

* Group 2 -	MN (1), MN (2), . . . , MN (M)

* Group 3 - 	X, Y, Z (as needed)

or

* Group 1 - 	KEYWORD

+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format      | Description                                                                                                                                                                                                                                                                                                 |
+================+=============+=============================================================================================================================================================================================================================================================================================================+
|  M             | integer     | Number of nodes for which information will be printed on the output (iout) and history plot (ishis, istrc) files. If M ≤ 0, pressure and temperature will be written on the output file for all nodes but no nodal parameter values will be printed in the history plot files. Group 2 is omitted if M ≤ 0. |
+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| MN             | integer     | M node numbers for which information will be printed on the output file (iout). If MN(I) < 0, then coordinates are used to define the print-out node, and the coordinate sets (X, Y, Z) for each MN(I) < 0 are added after Group 2.                                                                         |
+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| X              | real        | Coordinates of node for which information will be printed. One line for each MN < 0. The code finds the node closest to the coordinate given. For 2-D problems set Z = 0. No input if MN >0.                                                                                                                |
+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Y              | real        |                                                                                                                                                                                                                                                                                                             |
+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Z              | real        |                                                                                                                                                                                                                                                                                                             |
+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| KEYWORD        | character*5 | Keyword ‘block' to invoke node specification by JA, JB, JC format. Keyword ‘azone' to invoke output zone specification by JA, JB, JC format. This keyword allows a single node to be included in multiple output zones.                                                                                     |
+----------------+-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The following are examples of node. In the first example, 2 nodes are specified for output,
nodes 50 and 88.

.. code::

  node
      2
      50    88


In the second example, two nodes are specified for output, the node numbered
50 and the node closest to the coordinates X = 100. m, Y = 1000. m
and Z = 0. m.

.. code::

  node
      2
      50    -88
  
      100.  1000.    0.


In the third example, output is specified for the block of
nodes 1, 11, 21, 31, 41, 51, 61, 71, 81, 91 and for those nodes defined by
zone 3 (see macro zone).

.. code::

  node
  block
      1     100      10
      -3      0       0

In the fourth example, output is specified for
two zones, the first zone contains nodes 1, 11, 21, 31, 41, 51, 61, 71, 81, 91
and the second zone is made up of the nodes defined for zone 3
(previously specified using the zone macro).

.. code::

  node
  azone
     1      100      10
     -3       0       0  
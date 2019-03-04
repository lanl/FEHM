========
``pest``
========

Output variable information for PEST parameter estimation routine.

* Group 1 -	MPEST

* Group 2 -	NPEST(I), I = 1, MPEST

* Group 3 -	X, Y, Z (as needed)

+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                                                                                                            |
+================+=========+========================================================================================================================================================================+
| MPEST          | integer | Number of nodes for PEST output. At present the code outputs only pressures (heads), saturations, and temperatures.                                                    |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NPEST(I)       | integer | Node numbers printed to the output file (fehmn.pest) with values of variables listed above. If NPEST(I) < 0 then the node numbers are determined from the coordinates. |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| X, Y, Z        | real    | Coordinates in grid if NPEST(I) < 0. The coordinates are used to find the node closest in distance to that point and that node is substituted for NPEST(I).            |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The following is an example of pest. In this example pest output is specified at
5 nodes, nodes numbered 21, 23, 35, and 47, and the node closest to the coordinates
X=10. m, Y=15. m, Z=20. m.

.. code::

  pest
      21  23  35 47 -50
      10. 15. 20.
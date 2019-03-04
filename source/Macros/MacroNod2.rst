========
``nod2``
========

Specify the node numbers for which detailed file output is desired and alternate nodes for terminal output.

* Group 1 -	M, M2

* Group 2 -	MN (1), MN (2), . . . , MN (M)

* Group 3 -	X, Y, Z (as needed)

* Group 4 -	MNI(1), MNI(2), . . . , MNI(M2)

* Group 5 -	X, Y, Z (as needed) 

+----------------+---------+--------------------------------------------------------------+
| Input Variable | Format  | Description                                                  |
+================+=========+==============================================================+
|  M             | integer | Number of nodes for which information will be printed        |
|                |         | on the output file (iout). If M ≤ 0, pressure and            |
|                |         | temperature will be written on the output file for           |
|                |         | all nodes but no nodal parameter values will be              |
|                |         | printed in the history plot files. Group 2 is omitted        |
|                |         | if M ≤ 0.          |                                         |
+----------------+---------+--------------------------------------------------------------+
| M2             | integer | Number of nodes for short list (terminal printout).          |
|                |         | If M2 ≤ 0, Group 4 is omitted.                               |
|                |         |                                                              |
|                |         |                                                              |
|                |         |                                                              |
+----------------+---------+--------------------------------------------------------------+
| MN             | integer | M node numbers for which information will be printed         | 
|                |         | on the output file (iout). If a MN(I) < 0, then              | 
|                |         | coordinates are used to define that print-out node,          | 
|                |         | and the coordinate sets (X, Y, Z) for each MN(I) < 0         | 
|                |         | are added after Group 2.                                     | 
+----------------+---------+--------------------------------------------------------------+
| MNI            | integer | M2 node numbers for which information will be printed        |
|                |         | on the terminal (short list). This group exists only         |
|                |         | if M2 ≠ 0. If MNI(I) < 0, then coordinates are used          |
|                |         | to define the terminal output nodes, and the coordinate      |
|                |         | sets (X, Y, Z) for each MNI(I) < 0 are added after Group 4.  |
+----------------+---------+--------------------------------------------------------------+
| X              | real    | Coordinates of node for which information will be printed.   |
|                |         | One line for each MN or MNI < 0. The code finds the node     |
|                |         | closest to the coordinate given. For 2-D problems set        |
|                |         | Z = 0. No input if no MN or MNI < 0.                         |
+----------------+---------+--------------------------------------------------------------+
| Y              | real    |                                                              |
+----------------+---------+--------------------------------------------------------------+
| Z              | real    |                                                              |
+----------------+---------+--------------------------------------------------------------+

The following are examples of ``nod2``. In the first example (top), detailed output to the
output file is specified for two nodes, the nodes numbered 50 and 88, and one node is
specified for terminal output, node 50. In the second example (bottom), two nodes are
specified for detailed output, the nodes numbered 50 and 88, and one node is specified
for terminal output, the node closest to the coordinates X = 100. m, Y = 1000. m and
Z = 0. m. 

+------+----+
| nod2 |    |
+------+----+
| 2    | 1  |
+------+----+
| 50   | 88 |
+------+----+
| 50   |    |
+------+----+


+------+-------+----+
| nod2 |       |    |
+------+-------+----+
| 2    | 1     |    |
+------+-------+----+
| 50   | 88    |    |
+------+-------+----+
| -88  |       |    |
+------+-------+----+
| 100. | 1000. | 0. |
+------+-------+----+


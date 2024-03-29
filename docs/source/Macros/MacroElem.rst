========
``elem``
========

Element connectivity data. These data are usually created by a mesh generation program, then cut and copied into the input file or a separate geometry data input file. 

.. note::

  ``elem`` is required if macro ``fdm`` is not used.

* Group 1 - ``NS``, ``NEI``
* Group 2 - ``MB``, ``NELM (1)``, ``NELM (2)``, . . ., ``NELM (NS)``

If :math:`NS < 0` then :math:`ABS(NS)` is interpreted as the number of nodes per element.
:math:`NS < 0` signals the code to make rectangles (or bricks in three dimensions)
a sum of triangles (or tetrahedrals). This provides more stability in nonlinear problems
with a distorted mesh. `See Elements available with FEHM in 2-D and 3-D problems showing
nodal numbering convention. <Macro39151.html>`_ shows available element types
and the nodal numbering convention. To end the control section a blank line is entered. 

+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                                                                                                                                    |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``NS``         | integer | Number of nodes per element.                                                                                                                                                                   |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``NEI``        | integer | Number of elements                                                                                                                                                                             |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``MB``         | integer | Element number. ``If MB < 0`` then the difference between the absolute value of MB and the previous absolute value of MB is used to generate intermediate values by interpolation in the code. |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``NELM (1)``   | integer | First node of element MB                                                                                                                                                                       |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``NELM (2)``   | integer | Second node of element MB                                                                                                                                                                      |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ...                                                                                                                                                                                                                       |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``NELM (NS)``  | integer | Last node of element MB                                                                                                                                                                        |
+----------------+---------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The following is an example of **elem**. In this example there are 4 nodes per element,
i.e., the elements are 2-dimensional quadrilaterals. There are a total of
117 elements in the model, element number 1 is defined by nodes 15, 16, 2, and 1,
element number 2 is defined by nodes 16, 17, 3 and2, and so on.

+------+-----+-----+-----+-----+
| elem |     |     |     |     |
+------+-----+-----+-----+-----+
| 4    | 117 |     |     |     |
+------+-----+-----+-----+-----+
| 1    | 15  | 16  | 2   | 1   |
+------+-----+-----+-----+-----+
| 2    | 16  | 17  | 3   | 2   |
+------+-----+-----+-----+-----+
| .    | .   | .   | .   | .   |
+------+-----+-----+-----+-----+
| .    | .   | .   | .   | .   |
+------+-----+-----+-----+-----+
| .    | .   | .   | .   | .   |
+------+-----+-----+-----+-----+
| 10   | 24  | 25  | 11  | 10  |
+------+-----+-----+-----+-----+
| 11   | 25  | 26  | 12  | 11  |
+------+-----+-----+-----+-----+
| 12   | 26  | 27  | 13  | 12  |
+------+-----+-----+-----+-----+
| .    | .   | .   | .   | .   |
+------+-----+-----+-----+-----+
| .    | .   | .   | .   | .   |
+------+-----+-----+-----+-----+
| .    | .   | .   | .   | .   |
+------+-----+-----+-----+-----+
| 116  | 138 | 139 | 125 | 124 |
+------+-----+-----+-----+-----+
| 117  | 139 | 140 | 126 | 125 |
+------+-----+-----+-----+-----+
|      |     |     |     |     |
+------+-----+-----+-----+-----+


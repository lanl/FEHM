========
``grad``
========

Gradient model input.

* Group 1 - 	NGRAD

* Group 2 - 	IZONE_GRAD, CORDG, IDIRG, IGRADF, VAR0, GRAD1

Group 2 is repeated (NGRAD times) for each gradient model being defined.

+----------------+---------+----------------------------------------------------------+
| Input Variable | Format  | Description                                              |
+================+=========+==========================================================+
| NGRAD          | integer | Number of gradient models                                |
+----------------+---------+----------------------------------------------------------+
| IZONE_GRAD     | integer | Zone associated with model                               |
+----------------+---------+----------------------------------------------------------+
| CORDG          | real    | Reference coordinate (m)                                 |
+----------------+---------+----------------------------------------------------------+
| IDIRG          | integer | Direction of gradient (1, 2, or 3).                      |
+----------------+---------+----------------------------------------------------------+
| IGRADF         | integer | | Variable to which gradient is applied.                 |
|                |         | | 1 = initial pressure                                   |
|                |         | | 2 = initial temperature                                |
|                |         | | 3 = initial saturation                                 |
|                |         | | 4 = Fixed pressure                                     |
|                |         | | 5 = Fixed enthalpy                                     |
|                |         | | -5 = Fixed inflowing temperature                       |
|                |         | | 6 = Initial methane pressure                           |
|                |         | | 7 = Fixed Methane pressure                             |
|                |         | | 8 = depending on how ``hflx`` macro is configured      |
|                |         |   this is either a temperature or a heat flux            |
|                |         | | 9 = initial :math:`CO_2` pressure                      |
|                |         | | 10 = Fixed :math:`CO_2` cell pressure                  |
|                |         | | 11 = Pressure for secondary material in                |
|                |         | |   ``gdkm`` or ``gdpm`` model                           |
|                |         | | 12 = initial Temperature for matrix in ``gdkm``        |
|                |         |   or ``gdpm`` model                                      |
+----------------+---------+----------------------------------------------------------+
| VAR0           | real    | Value of variable at reference point (m).                |
+----------------+---------+----------------------------------------------------------+
| GRAD1          | real    | | Gradient.                                              |
|                |         | | Units: Pressure MPa/m, T degrees C/m,                  |
|                |         |   enthalpy KJ/kg/m, heat flux MW/m                       |
+----------------+---------+----------------------------------------------------------+

.. note::

  IIGRADF = 4,5,-5,7 requires that the node is previously defined as a boundary
  node in a ``flow`` macro (or equivalent)

  IGRAFD = 11 or 12 requires ``gdkm`` or ``gdpm`` macro
  
  IGRADF = 8 requires a ``hflx`` macro


The following is an example of the grad macro. A temperature gradient in the Y direction from the reference point of 0 will be applied to zone 1. 

+------+----+---+---+-----+-------+
| grad |    |   |   |     |       |
+------+----+---+---+-----+-------+
| 1    |    |   |   |     |       |
+------+----+---+---+-----+-------+
| 1    | 0. | 2 | 2 | 10. | -150. |
+------+----+---+---+-----+-------+
   

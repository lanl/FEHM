========
``subm``
========

Create a new flow macro to represent boundary conditions on an extracted submodel.

Group 1 - KEYWORD, IZONE1, IZONE2 

+----------------+-------------+------------------------------------------------------------------------------------+
| Input Variable | Format      | Description                                                                        |
+================+=============+====================================================================================+
| KEYWORD        | character*4 | Keyword "flux", "head", or "pres" to specify type of boundary condition to output. |
+----------------+-------------+------------------------------------------------------------------------------------+
| IZONE1         | integer     | Zone defining submodel nodes.                                                      |
+----------------+-------------+------------------------------------------------------------------------------------+
| IZONE2         | integer     | Zone defining nodes outside of the submodel (optional).                            |
+----------------+-------------+------------------------------------------------------------------------------------+

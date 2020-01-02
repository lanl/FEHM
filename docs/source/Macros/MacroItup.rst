========
``itup``
========

Controls upstream direction. The use of the itup macro is sometimes useful in
problems where the flow directions are changing rapidly. The parameter UPWGT
(in macro ctrl) must be greater than 0.5 for this macro to have any effect.

* Group 1 - IAD_UP 

+----------------+---------+---------+-------------------------------------------+
| Input Variable | Format  | Default | Description                               |
+================+=========+=========+===========================================+
| IAD_UP         | integer | 1000    | | Number of iterations after which the    |
|                |         |         |   upwind directions are held constant.    |
|                |         |         | | *A value of 2 is suggested*             |
+----------------+---------+---------+-------------------------------------------+

In the following example of itup, after 10 iterations the upwind directions are
held constant. 

+------+
| itup |
+------+
| 10   |
+------+
========
``wtsi``
========

Water table, simplified.

Group 1 - NFREE, IZONE_FREE(I), I = 1 to NFREE, IFREE_IM_EX, HEAD_TOL, RLPTOL,

+----------------+---------+--------------------------------------------------------------+
| Input Variable | Format  | Description                                                  |
+================+=========+==============================================================+
| NFREE          | integer | Number of zones to which water table condition are applied.  |
+----------------+---------+--------------------------------------------------------------+
| IZONE_FREE     | integer | Zone number.                                                 |
+----------------+---------+--------------------------------------------------------------+
| IFREE_IM_EX    | integer | Update parameter (0 - explicit update, 1 - implicit update). |
+----------------+---------+--------------------------------------------------------------+
| HEAD_TOL       | real*8  | Tolerance for head (m)                                       |
+----------------+---------+--------------------------------------------------------------+
| RLPTOL         | real*8  | Tolerance for saturation                                     |
+----------------+---------+--------------------------------------------------------------+

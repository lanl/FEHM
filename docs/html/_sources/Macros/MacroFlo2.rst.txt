========
``flo2``
========

* Group 1 - JA, JB, JC, JD, SKD, EFLOW, AIPED

Multiple lines of input may be used, terminated by a blank line.

+----------------+---------+--------------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                          |
+================+=========+======================================================================================+
| JA             | integer | Indices used to define planes in a 3-D simulation with a regular numbering pattern.  |
+----------------+---------+                                                                                      +
| JB             | integer |                                                                                      |
+----------------+---------+                                                                                      +
| JC             | integer |                                                                                      |
+----------------+---------+                                                                                      +
| JD             | integer |                                                                                      |
+----------------+---------+--------------------------------------------------------------------------------------+


The flow rates are defined within the inner loop of the do loops:

.. code::

  DO JK = JA, JB
    KL = JK - JA
    DO IJ = JA + KL, JC + KL, JD
        ⋅ ⋅ ⋅
    ENDDO
  ENDDO
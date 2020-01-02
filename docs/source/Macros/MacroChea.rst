========
``chea``
========

Convert output from pressure to head (non-head problems). The reference temperature and pressure are used to calculate density for the head calculation. For this macro the data are entered on the macro line, and if omitted the specified default values are used. (Note that all five values must be entered to override the default values.)

+----------------+--------+---------+--------------------------------------------------+
| Input Variable | Format | Default | Description                                      |
+================+========+=========+==================================================+
| HEAD0          | real   | 0.      | Reference head (m)                               |
+----------------+--------+---------+--------------------------------------------------+
| TEMP0          | real   | 20.     | Reference temperature (oC)                       |
+----------------+--------+---------+--------------------------------------------------+
| PRES0          | real   | 0.1     | Reference pressure (MPa)                         |
+----------------+--------+---------+--------------------------------------------------+
| SAT_ICH        | real   | 0.      | Saturation adjustment after variable switch      |
+----------------+--------+---------+--------------------------------------------------+
| HEAD_ID        | real   | 0.      | Output head identification for small saturations |
+----------------+--------+---------+--------------------------------------------------+

The following is an example of chea. In this example pressures will be converted to heads for output using a reference pressure of 1 MPa and a reference temperature of 25 C.

.. code::
  
  chea 0. 25. 1. 0. 0.
========
``conv``
========

Convert input from head to pressure. Often used when converting a head-based isothermal model to a heat and mass simulation with pressure and temperature variables. The reference temperature and head are used to calculate density for the head calculation. It should be noted that this is an approximate method. Since the density is a nonlinear function of pressure and temperature this method will give slightly different answers than a calculation allowing a water column to come to thermal and mechanical equilibrium.

The reference head (``head0``) is converted to a pressure and added to the reference pressure (``conv1``) and this sum is used with the reference temperature to calculate a density. This density is used to convert the head to pressure for the identified zone. The option of adding a temperature gradient is provided as well.

The reference head (``head0``) is entered on the macro line:

.. code::

  conv HEAD0

* Group 1 - NCONV

* Group 2 – ZONE_CONV, ICONVF, CONV1, CONV2, CORDC, IDIRC, VARC

Group 2 is entered for each zone (nconv times).

+----------------+---------+---------------------------------------------------------+
| Input Variable | Format  | Description                                             |
+================+=========+=========================================================+
| HEAD0          | real    | Reference head (m)                                      |
+----------------+---------+---------------------------------------------------------+
| NCONV          | integer | Number of zones for variable conversion                 |
+----------------+---------+---------------------------------------------------------+
| ZONE_CONV      | integer | Zone for variable conversion                            |
+----------------+---------+---------------------------------------------------------+
| ICONVF         | integer | iconvf : 1- initial conditions, 2-boundary (fixed head) |
+----------------+---------+---------------------------------------------------------+
| CONV1          | real    | Reference pressure (MPa)                                |
+----------------+---------+---------------------------------------------------------+
| CONV2          | real    | Reference temperature (oC)                              |
+----------------+---------+---------------------------------------------------------+
| CORDC          | real    | Reference coordinate (m)                                |
+----------------+---------+---------------------------------------------------------+
| IDIRC          | integer | Coordinate direction (1 - X, 2 - Y, 3 - Z)              |
+----------------+---------+---------------------------------------------------------+
| VARC           | real    | Temperature gradient (oC/m)                             |
+----------------+---------+---------------------------------------------------------+

The following is an example of conv. Here the density used to convert head to pressure is calculated with a reference head of 1000 m plus 0.1 MPa and 80 °C. The nodes in zone 45 are converted from heads to pressures at 80 °C. 

.. code::

   conv
   1000.
   1
   45
   1
   0.4
   80.
   0.
   0
   0.0


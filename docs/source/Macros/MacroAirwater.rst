===================
``airwater or air``
===================

Isothermal air-water two-phase simulation. 

* Group 1 -	ICO2D

* Group 2 -	TREF, PREF

+----------------+---------+------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                  |
+================+=========+==============================================================================+
| ICO2D          | integer | | Determines the type of air module used.                                    |
|                |         | | ICO2D = 1, 1 degree of freedom solution to the saturated-unsaturated       |
|                |         | |   problem is produced. This formulation is similar to the Richard's        |
|                |         | |   Equation.                                                                |
|                |         | | ICO2D = 2, 1 degree of freedom solution is obtained assuming only gas flow |
|                |         | |   with no liquid present.                                                  |
|                |         | | ICO2D = 3, full 2 degree of freedom solution.                              |
|                |         | | All other values are ignored. The default is 3.                            |
+----------------+---------+------------------------------------------------------------------------------+
| TREF           | real    | Reference temperature for properties (oC).                                   |
+----------------+---------+------------------------------------------------------------------------------+
| PREF           | real    | Reference pressure for properties (MPa).                                     |
+----------------+---------+------------------------------------------------------------------------------+

Several macros are affected if the air module is enabled. These are:

* **pres** - Because the air-water formulation is 2-phase at all times, care should be taken to insure that IEOSD is always specified to be 2. Likewise, saturations (not temperatures) are used.

* **init** - This macro should not be used because the saturation values cannot be specified.

* **flow** - A variety of different flow and boundary values are input with this macro when the macro **airwater** is also used. See description of control statement **flow**.

The following is an example of **airwater**. In this example, a full 2-degrees-of-freedom solution is specified with a reference temperature for property evaluation of 20 oC and a reference pressure of 0.1 MPa.

+------------+----------------+
|        airwater             |
+============+================+
| 3          |                |
+------------+----------------+
| 20.        | 0.1            |
+------------+----------------+

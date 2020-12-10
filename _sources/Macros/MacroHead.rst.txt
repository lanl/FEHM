========
``head``
========

Hydraulic head values are used for input and output instead of pressures. Use of this macro enables the Boussinesq approximation (bous macro) and isothermal air-water two-phase simulation (airwater macro) automatically. It affects the pres and flow macros by requiring head information where pressure values were previously required. The default is to have no input associated with this macro. However, an optional head increment can be given after the head macro keyword. This value will be added to all input head values to ensure a single phase fluid state. Note that the value will be subtracted before output is written.

+----------------+--------+---------+-----------------------------------------------------------------+
| Input Variable | Format | Default | Description                                                     |
+================+========+=========+=================================================================+
| HEAD0          | real   | 0.      | An incremental value that will be added to all input heads (m). |
+----------------+--------+---------+-----------------------------------------------------------------+

The following is an example of head. In this example the optional head increment is included and a value of 1000. m is added to all input head values.

+------+-------+
| head | 1000. |
+------+-------+

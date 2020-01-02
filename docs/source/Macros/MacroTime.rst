========
``time``
========

Time step and time of simulation data. 

* Group 1 - DAY, TIMS, NSTEP, IPRTOUT, YEAR, MONTH, INITTIME

* Group 2 - DIT1, DIT2, DIT3, ITC, DIT4 (as needed)

DAY should be larger than DAYMIN defined in control statement **ctrl**.
The code proceeds to the next control statement when a blank line is encountered for Group 2.
Group 2 can be used to generate output at specific times (with multiple Group 2s).
Contour plot output will be written at each DIT1 regardless of the input in control statement **cont**.
The restart file will be written (or rewritten if one already exists) at each DIT1. If DIT4 is omitted
(for compatibility with older input files where DIT4 was not input) the maximum
time step defined in the control statement ctrl will be used. 

+----------------+---------+-----------------------------------------------------------------+
| Input Variable | Format  | Description                                                     |
+================+=========+=================================================================+
| DAY            | real    | Initial time step size (days).                                  |
+----------------+---------+-----------------------------------------------------------------+
| TIMS           | real    | Final simulation time (days).                                   |
+----------------+---------+-----------------------------------------------------------------+
| NSTEP          | integer | Maximum number of time steps allowed.                           |
+----------------+---------+-----------------------------------------------------------------+
| IPRTOUT        | integer | Print-out interval for nodal information (pressure,             |
|                |         | enthalpy etc.), as set up under control statement node.         |
|                |         | (i.e., number of time steps).                                   |
+----------------+---------+-----------------------------------------------------------------+
| YEAR           | integer | Year that simulation starts.                                    |
+----------------+---------+-----------------------------------------------------------------+
| MONTH          | integer | Month that simulation starts.                                   |
+----------------+---------+-----------------------------------------------------------------+
| INITTIME       | real    | Initial time of simulation (days). For compatibility            |
|                |         | with older versions, if this parameter is absent the            |
|                |         | initial time of simulation will be 0 if no restart              |
|                |         | file is used, or the time in the restart file if one is used.   |
+----------------+---------+-----------------------------------------------------------------+
| DIT1           | real    | Time (days) for time step change.                               |
+----------------+---------+-----------------------------------------------------------------+
| DIT2           | real    | New time step size (days). If DIT2 < 0 then ABS (DIT2)          |
|                |         | is the new time step multiplier.                                |
+----------------+---------+-----------------------------------------------------------------+
| DIT3           | real    | | Implicitness factor for new time step.                        |
|                |         | | DIT3 ≤ 1.0 backward Euler.                                    |
|                |         | | DIT3 > 1.0 for second-order implicit scheme.                  |
+----------------+---------+-----------------------------------------------------------------+
| ITC            | integer | New print-out interval.                                         |
+----------------+---------+-----------------------------------------------------------------+
| DIT4           | real    | Maximum time step size for next time interval (days).           |
+----------------+---------+-----------------------------------------------------------------+

The following is an example of ``time``. In this example, the initial time step size is
30 days, the final simulation time is 3650 days, the number of time steps allowed is
20, nodal information is printed out for every 5th time step, the simulation starts
in the 10th month of 1989, and the initial time of simulation is assigned a value of
0.

The time step multiplier is changed after 1 day, and the new time step multiplier
is 1.2, backward Euler is used from this time on and the printout interval is every
10th time step. The maximum time step size for the next interval is omitted so the
default value entered in the ``ctrl`` macro will be used.

+------+--------+-----+----+------+----+-----+
| time |        |     |    |      |    |     |
+------+--------+-----+----+------+----+-----+
| 30.0 | 3650.0 | 20  | 5  | 1989 | 10 | 0.0 |
+------+--------+-----+----+------+----+-----+
| 1.0  | -1.2   | 1.0 | 10 |      |    |     |
+------+--------+-----+----+------+----+-----+
|      |        |     |    |      |    |     |
+------+--------+-----+----+------+----+-----+

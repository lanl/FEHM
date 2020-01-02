========
``stea``
========

The macro ``stea`` is used to manage steady state simulations and is available with all
physics modules in FEHM. The macro directs FEHM to monitor changes in variables
from timestep to timestep and stop the steady state run when the changes are less
than some prescribed tolerance. Alternatively the steady state run is directed
to finish when the global "in" and "out" fluxes are less than a prescribed
tolerance or the simulated time exceeds the input value. 

After the steady state portion of the simulation is completed, a transient run may
be performed. This is accomplished with the boun macro and the key word ``tran``.
See the description of the ``boun`` macro for details.

The user should be aware that when the ``stea`` macro is used, the parameters
associated with the ``time`` macro pertain to the transient portion of the
simulation if a transient part exists. Values for these parameters may be
input using a keyword but if not entered will default to the values specified
for the ``time`` macro.


* Group 1 - KEYWORD, VALUE

+----------------+-----------+-------------------------------------------------------------------------------+
| Input Variable | Format    | Description                                                                   |
+================+===========+===============================================================================+
| KEYWORD        | character | | The following keywords are used with steady to specify                      |
|                |           |   the variables to be checked for steady state:                               |
|                |           | | shea - Head (m)                                                             |
|                |           | | spre - Pressure (MPa)                                                       |
|                |           | | stem - Temperature (oC)                                                     |
|                |           | | ssat - Saturation                                                           |
|                |           | | sair - Partial pressure of air/gas (MPa)                                    |
|                |           | | sflu - Mass flux (kg/s)                                                     |
|                |           | | sent - Enthalpy (MJ/s)                                                      |
|                |           | | stim - Maximum time for steady state simulation (days)                      |
|                |           | | sday - Initial time step size for steady state simulation (days)            |
|                |           | | smul - Time step multiplication factorsmst - Minimum                        |
|                |           |   number of time steps to be used for steady state simulation                 |
|                |           | | snst - Maximum number of time steps to be used for steady                   |
|                |           |   state simulation                                                            |
|                |           | | shtl - Option to reduce the head_tol factor as the                          |
|                |           |   solution approaches steady-state                                            |
|                |           | | stmc - Option to reduce the machine tolerancs factor (tmch)                 |
|                |           |   factor as the solution approaches steady-state                              |
|                |           | | sacc - Maximum change allowed in the accumulation term when                 |
|                |           |   flux is being checked                                                       |
|                |           | | sper - The tolerance is interpreted as a fractional change                  |
|                |           |   in the variable being checked [i.e.,                                        |
|                |           |   (new_value - old_value)/old_value]. Without this keyword                    |
|                |           |   it is an absolute change in the variable value.                             |
|                |           | | endstea - Signifies end of keyword input, a blank line will also work.      |
+----------------+-----------+-------------------------------------------------------------------------------+
| VALUE          | real      | Variable tolerance or time control parameter value.                           |
+----------------+-----------+-------------------------------------------------------------------------------+

In the following example a steady state solution is specified. The tolerance for
``head`` is specified to be 0.1 m and for ``flux`` 0.00001kg/s. The steady state solution
will be allowed to run for a maximum of 1.e12 days and the time step multiplier
is set to 2. 

+--------+-------+
| stea   |       |
+========+=======+
| shead  | 1.d-1 |
+--------+-------+
| stime  | 1.e12 |
+--------+-------+
| smult  | 2.    |
+--------+-------+
| sflux  | 1.d-5 |
+--------+-------+
| end    |       |
+--------+-------+

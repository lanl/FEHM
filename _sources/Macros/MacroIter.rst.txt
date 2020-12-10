========
``iter``
========

.. note::

  This control statement is optional but recommended. However, if the user is not
  familiar with the linear equation solver routines in FEHM (Zyvoloski and Robinson, 1995)
  control statement **iter** should not be used. 

* Group 1 -	G1, G2, G3, TMCH, OVERF

* Group 2 -	IRDOF, ISLORD, IBACK, ICOUPL, RNMAX

The parameters G1, G2, and G3 are used to calculate the completion criteria for
the linear equation solver. The equation for the stopping criteria is:

:math:`EPE = G3 * \mathrm{max}(TMCH, \mathrm{MAX}(F0, \mathrm{MIN(G1*\mathrm{SQRT(R^2), G2*R^2})}))`

where :math:`R^2` is the sum-squared of the equation residuals, and F0 is the
:math:`SQRT(R0^2)*EPM` for the first iterration (see macro ``ctrl`` for a 
definition of EPM). The other parameters are defined below:

+----------------+-----------+---------+--------------------------------------------------------------------+
| Input Variable | Format    | Default | Description                                                        |
+================+===========+=========+====================================================================+
| G1             | real      | 1.e-6   | Multiplier for the linear convergence region of the                |
|                |           |         | Newton-Raphson iteration.                                          |
+----------------+-----------+---------+--------------------------------------------------------------------+
| G2             | real      | 1.e-6   | Multiplier for the quadratic convergence region of                 |
|                |           |         | the Newton-Raphson iteration.                                      |
+----------------+-----------+---------+--------------------------------------------------------------------+
| G3             | real      | 1.e-3   | Multiplier relating Newton Raphson residuals to                    |
|                |           |         | stopping criteria for linear solver                                |
+----------------+-----------+---------+--------------------------------------------------------------------+
| TMCH           | real      | 1.e-9   | Machine tolerance if TMCH > 0. If satisfied by                     |
|                |           |         | the residual norm, the Newton iteration is assumed to              |
|                |           |         | be complete. Newton-Raphson stopping criteria if TMCH < 0          |
|                |           |         | (recommended). If TMCH < 0 then the ABS(TMCH) is used as a         |
|                |           |         | tolerance for each equation at each node. Convergence is           |
|                |           |         | achieved if the residual of every equation at every                |
|                |           |         | node is < ABS(TMCH).                                               |
+----------------+-----------+---------+--------------------------------------------------------------------+
| OVERF          | real      | 1.1     | Over relaxation factor for passive nodes in adaptive               |
|                |           |         | implicit method.                                                   |
+----------------+-----------+---------+--------------------------------------------------------------------+
| IRDOF          | integer   | 0       | | Enables the reduced degree of freedom method. If                 |
|                |           |         |   IRDOF = 0, reduced degrees of freedom are not required.          |
|                |           |         |   When IRDOF = 1, a reduced degree of freedom from 3 to 2          |
|                |           |         |   or 3 to 1 is used. When IRDOF = 2, a reduced degree of           |
|                |           |         |   freedom from 3 to 2 is used. If IRDOF =11, then an air           |
|                |           |         |   only solution is found for the isothermal air-water              |
|                |           |         |   process model. If IRDOF = -11, then the residual for             |
|                |           |         |   the air equation with the airwater macro is ignored.             |
|                |           |         |   If IRDOF = 13, then a liquid only solution for the               |
|                |           |         |   airwater macro is assumed. {0}                                   |
|                |           |         | |                                                                  |
|                |           |         | | Examples of 1, 2, 3, 4 and 6 degrees of freedom models are:      |
|                |           |         | | 1 - heat only or mass only.                                      |
|                |           |         | | 2 - heat and mass, or air-water (isothermal)                     |
|                |           |         | | 3 - air-water with heat (non-isothermal)                         |
|                |           |         | | 4 - heat and mass, double permeability or air-water              |
|                |           |         |   (isothermal), double permeability                                |
|                |           |         | | 6 - air-water with heat, double permeability                     |
|                |           |         | |                                                                  |
|                |           |         | | See Tseng and Zyvoloski (2000) for more information              |
|                |           |         |   on the reduced degree of freedom method.                         |
+----------------+-----------+---------+--------------------------------------------------------------------+
| ISLORD         | integer   | 0       | Reordering parameter. The value of ISLORD and the                  |
|                |           |         | corresponding equation order is given below. The ordering          |
|                |           |         | has an effect on the speed of convergence of several               |
|                |           |         | solution algorithms, but will not affect most users.               |
|                |           |         | For problems of order 2 or greater, the ordering can be            |
|                |           |         | understood by labeling each equation. For example for a            |
|                |           |         | 3-degree of freedom problem with mass, heat, and noncondensible    |
|                |           |         | gas, label the mass equation as 1, the heat equation as 2,         |
|                |           |         | and the noncondensible gas equation as 3. In general mass          |
|                |           |         | (water), heat or air, air. For double permeability problems        |
|                |           |         | fracture equations precede matrix equations, i.e., for an          |
|                |           |         | air-water problem - mass water fracture, mass air fracture,        |
|                |           |         | mass water matrix, mass air matrix. {0}                            |
+----------------+-----------+---------+--------------------------------------------------------------------+
| IBACK          | integer   | 0       | IRDOF parameter. If IBACK = 0, SOR iterations are not              |
|                |           |         | performed before call to solver. If IBACK = 1,                     |
|                |           |         | SOR iterations are performed before call to solver.                |
|                |           |         | If IBACK = 2, SOR iterations are performed before                  |
|                |           |         | call to SOLVER, and SOLVER is called twice. {0}                    |
+----------------+-----------+---------+--------------------------------------------------------------------+
| ICOUPL         | integer   | 0       | Number of SOR iterations used in reduced degree of                 |
|                |           |         | freedom methods. {0}                                               |
+----------------+-----------+---------+--------------------------------------------------------------------+
| RNMAX          | real      | 1.0e+11 | Maximum running time for problem before the solution               |
|                |           |         | is stopped (cpu minutes).                                          |
+----------------+-----------+---------+--------------------------------------------------------------------+


{0}
---

+--------+----------------------+----------------------+----------------------+----------------------+
| ISLORD | 2 Degrees of Freedom | 3 Degrees of Freedom | 4 Degrees of Freedom | 6 Degrees of Freedom |
+--------+----------------------+----------------------+----------------------+----------------------+
| 0      | 1,2                  | 1,2,3                | 1,2,3,4              | 1,2,3,4,5,6          |
+--------+----------------------+----------------------+----------------------+----------------------+
| 1      | 2,1                  | 1,3,2                | 1,3,2,4              | 1,4,2,5,3,6          |
+--------+----------------------+----------------------+----------------------+----------------------+
| 2      |                      | 2,1,3                |                      |                      |
+--------+----------------------+----------------------+----------------------+----------------------+
| 3      |                      | 2,3,1                |                      |                      |
+--------+----------------------+----------------------+----------------------+----------------------+


The following is an example of ``iter``.
In this example, the tolerances for the linear and quadratic convergence regions
for the Newton-Raphson method are specified to be 1.e-5 times the initial residual,
tolerance for the adaptive-implicit method is 1.e-5, machine tolerance is 1.e-9,
and over-relaxation factor is 1.2. The reduced degree of freedom method is enabled,
reordering is not done, SOR iterations are not performed before calling the solver,
two SOR iterations are used in the reduced degree of freedom method, and the solution
procedure is terminated if not completed within 200 CPU minutes.

+-------+-------+-------+-------+-------+
| iter                                  |
+-------+-------+-------+-------+-------+
| 1.e-5 | 1.e-5 | 1.e-5 | 1.e-9 | 1.2   |
+-------+-------+-------+-------+-------+
| 1     | 0     | 0     | 2     | 200.0 |
+-------+-------+-------+-------+-------+




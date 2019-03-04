========
``cden``
========

Use concentration-dependent density for flow.

The following restrictions apply to the use of this macro:

1) It cannot be used with the macro "head", which assumes constant fluid density
2) The updating of density is explicit, based on the concentration values at the previous time step. Therefore, accuracy of the solution must be tested by using a smaller time step and ensuring that the results have converged
3) The fluid flow time steps should be small enough that only one or two solute time steps are carried out before the next fluid time step, because relatively small changes in the concentration field are required for accuracy; and
4) The heat and mass transfer solution must be kept on during the entire simulation for the results to be meaningful (see macro trac).

Group 1 - ISPCDEN

Group 2 - FACTCDEN

+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                                                                                                                                                                    |
+================+=========+================================================================================================================================================================================================================================+
| ISPCDEN        | integer | The number of the chemical component in trac that is used for applying the concentration-dependent density.                                                                                                                    |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| FACTCDEN       | real    | The factor used in the following relationship for fluid density (kg/m3):density = density_water + FACTCDEN*C where density_water = the density of pure water (kg/m3), and C is the concentration of chemical component ISPCDEN |
+----------------+---------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The following is an example of cden. In this example, component number 1 in trac is used. For concentrations of order 1, the density correction would be 100, of order 10% of the nominal value of water density of 1000 kg/m3. 

+----------+------+
|      cden       |
+==========+======+
| ISPCDEN  | 1    |
+----------+------+
| FACTCDEN | 100. |
+----------+------+


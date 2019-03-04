========
``adif``
========

Air-water vapor diffusion. The air-water diffusion equation is given as Equation (21)
of the "Models and Methods Summary" of the FEHM Application (Zyvoloski et al. 1999). 

* Group 1 - TORT

+----------------+--------+----------------------------------------------------------------------------+
| Input Variable | Format | Description                                                                |
+================+========+============================================================================+
| TORT           | real   | | Tortuosity for air-water vapor diffusion.                                |
|                |        | | If TORT > 0, :math:`\tau` of eqn 21, otherwise                           |
|                |        | | if TORT < 0, :math:`abs(\tau \phi S_v)` of the same equation.            |
|                |        | | If TORT > 1, water-vapor diffusion coefficient is set equal to the value |
|                |        | | specified for the first vapor species defined in the trac macro.         |
+----------------+--------+----------------------------------------------------------------------------+

The following is an example of ``adif``. In this example the tortuosity (:math:`\tau`)
for vapor diffusion is specified to be 0.8. 

+------+
| adif |
+------+
| 0.8  |
+------+


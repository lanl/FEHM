========
``ppor``
========

* Group 1 -	IPOROS
* Group 2 - JA, JB, JC, POR1, POR2, POR3, POR4 

+-------------------------------+---------+-------------------------------------------------------------------------------+
| Input Variable                | Format  | Description                                                                   |
+===============================+=========+===============================================================================+
| IPOROS                        | integer | | Model type:                                                                 |
|                               |         | | IPOROS = 1, aquifer compressibility model                                   |
|                               |         | | IPOROS = -1, specific storage model (use                                    |
|                               |         | | only for isothermal conditions)                                             |
|                               |         | | IPOROS = -2, Gangi model (not available for                                 |
|                               |         | | air-water-heat conditions)                                                  |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| **Model (1):** IPOROS = 1,    |         | | Aquifer compressibility                                                     |
|                               |         |   :math:`\phi = \phi_0 + \alpha_a(P-P_0)` where                               |
|                               |         | | :math:`\alpha_a` = aquifer compressibility (MPa-1),                         |
|                               |         | | :math:`\phi_0` = initial porosity,                                          |
|                               |         | | :math:`P_0` = initial pressure (MPa)                                        |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| POR1                          | real    | Aquifer compressibility :math:`\alpha (MPa^{-1})`                             |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| **Model (-1):** IPOROS = -1,  |         | | Specific storage                                                            |
|                               |         |   :math:`S_s = \rho(\alpha_a + \phi \beta)` where                             |
|                               |         | | :math:`\rho` = liquid density (kg/m3),                                      |
|                               |         | | :math:`g` = gravity,                                                        |
|                               |         | | :math:`\alpha_a` = aquifer compressibility (MPa-1),                         |
|                               |         | | :math:`\phi` = porosity,                                                    |
|                               |         | | :math:`\beta` = liquid compressibility (MPa-1)                              |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| POR1                          | real    | Specific storage :math:`S_S (m^{-1})`                                         |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| **Model (-2):** IPOROS = -2,  |         | | Gangi model with calculation of initial permeability                        |
|                               |         |   and porosity.                                                               |
|                               |         | | :math:`\phi = \phi_0 \left[ 1 - \left(\frac{P_c}{P_x}\right)^m \right]`     |
|                               |         |   and :math:`P_c = \sigma - P - \alpha E(T-T_0)`                              |
|                               |         | | where                                                                       |
|                               |         | | :math:`\phi_0` = initial porosity,                                          |
|                               |         | | :math:`m` = Gangi exponent,                                                 |
|                               |         | | :math:`P_x` = fitted parameter (MPa)                                        |
|                               |         | |                                                                             |
|                               |         | | Note: for the Gangi model the permeability is varied by                     |
|                               |         |   :math:`k = k_0 \left(\frac{\phi}{\phi_0}\right)^3`                          |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| POR1                          | real    | Exponent :math:`m` in Gangi bed of nails model.                               |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| POR2                          | real    | :math:`P_x` parameter (MPa) in Gangi equation.                                |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| POR3                          | real    | :math:`\sigma` in-situ stress (MPa).                                          |
+-------------------------------+---------+-------------------------------------------------------------------------------+
| POR4                          | real    | | :math:`(\sigma E)` The product of the coefficient                           |
|                               |         |   of thermal expansion for the rock and the Young's                           |
|                               |         |   modulus (MPa/C).                                                            |
|                               |         | |                                                                             |
|                               |         | | Note: For isothermal simulations the thermal term does not apply.           |
+-------------------------------+---------+-------------------------------------------------------------------------------+


In the following example of ppor, aquifer compressibility is modeled. All nodes
in the model are assigned a compressibility of 1.e-2 MPa-1. 

+------+---+---+-------+
| ppor |   |   |       |
+------+---+---+-------+
| 1    |   |   |       |
+------+---+---+-------+
| 1    | 0 | 0 | 1.e-2 |
+------+---+---+-------+
|      |   |   |       |
+------+---+---+-------+



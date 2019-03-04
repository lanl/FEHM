=======
``rlp``
=======

Relative permeability and capillary pressure model. Several models are available. 

* Group 1 - IRLP(i), RP1, RP2, RP3, RP4, RP5, RP6, RP7, RP8, RP9, RP10, RP11, RP12, RP13, RP14, RP15, RP16 (number of parameters entered depends on model selected)

* Group 2 - JA, JB, JC, I

Only those parameters defined for a given model need to be input. Group 1 is
ended when a blank line is encountered. The parameter i is incremented each
time a Group 1 line is read. Group 2 lines will refer to this parameter.
For model numbers 4, 6, and 7 (the combined van Genuchten model), the
permeability is isotropic and overwrites the input from macro ``perm``.
Macro fper can be used with models 4, 6, and 7 to introduce anisotropy.

 

+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| Input Variable                                                        | Format  | Description                                                                       |
+=======================================================================+=========+===================================================================================+
| IRLP(i)                                                               | integer | Relative permeability model type.                                                 |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model -1:** IRLP(i) = -1, constant relative permeability, linear                                                                                                  |
| capillary pressure (4 parameters required).                                                                                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP1                                                                   | real    | Liquid relative permeability (:math:`m^2`).                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP2                                                                   | real    | Vapor relative permeability (:math:`m^2`).                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP3                                                                   | real    | Capillary pressure at zero saturation (:math:`MPa`).                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP4                                                                   | real    | Saturation at which capillary pressure goes to zero.                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 1:**  IRLP(i) = 1, linear relative permeability, linear                                                                                                     |
| capillary pressure (6 parameters required).                                                                                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP1                                                                   | real    | Residual liquid saturation.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP2                                                                   | real    | Residual vapor saturation.                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP3                                                                   | real    | Maximum liquid saturation.                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP4                                                                   | real    | Maximum vapor saturation.                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP5                                                                   | real    | Capillary pressure at zero saturation (MPa).                                      |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP6                                                                   | real    | Saturation at which capillary pressure goes to zero.                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 2**:  IRLP(i) = 2, Corey relative permeability, linear                                                                                                      |
| capillary pressure (4 parameters required).                                                                                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP1                                                                   | real    | Residual liquid saturation.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP2                                                                   | real    | Residual vapor saturation.                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP3                                                                   | real    | Capillary pressure at zero saturation (MPa).                                      |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP4                                                                   | real    | Saturation at which capillary pressure goes to zero.                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 3**:  IRLP(i) = 3, van Genuchten relative permeability, van                                                                                                 |
| Genuchten capillary pressure (6 parameters required). In this                                                                                                       |
| model permeabilities are represented as a function of capillary                                                                                                     |
| pressure [rlp(h)].                                                                                                                                                  |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP1                                                                   | real    | Residual liquid saturation.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP2                                                                   | real    | Maximum liquid saturation.                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP3                                                                   | real    | Inverse of air entry head, :math:`\alpha_G` (1/m) [note                           |
|                                                                       |         | some data is given in (1/Pa) convert using pressure = ρgΔh].                      |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP4                                                                   | real    | Power n in van Genuchten formula.                                                 |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP5                                                                   | real    | Low saturation fitting parameter, multiple of                                     |
|                                                                       |         | cutoff capillary pressure assigned as maximum                                     |
|                                                                       |         | capillary pressure. If RP5 < 0 then a linear                                      |
|                                                                       |         | fit from the cutoff saturation (RP6) is used.                                     |
|                                                                       |         | The slope of the cutoff saturation is used                                        |
|                                                                       |         | to extend the function to saturation = 0.                                         |
|                                                                       |         | If RP5 = 0, a cubic fit is used. The slope                                        |
|                                                                       |         | at the cutoff saturation is matched and the                                       |
|                                                                       |         | conditions :math:`\frac{\partial}{\partial S}Pcap = 0`                            |
|                                                                       |         | and :math:`\frac{\partial^2}{\partial S}Pcap = 0`                                 |
|                                                                       |         | are forced at :math:`S = 0`.                                                      |
|                                                                       |         | If RP5 > 0, a multiple of the                                                     |
|                                                                       |         | value of the capillary pressure at                                                |
|                                                                       |         | the cutoff saturation,  :math:`RP5\cdot Pcap(S_{cutoff}`                          |
|                                                                       |         | is forced at :math:`S = 0`.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP6                                                                   | real    | Cutoff saturation used in fits described for RP5,                                 |
|                                                                       |         | must be greater than RP1.                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 4**:  IRLP(i) = 4, van Genuchten relative permeability, van                                                                                                 |
| Genuchten capillary pressure, effective continuum (15 parameters                                                                                                    |
| required). In this model permeabilities are represented as a                                                                                                        |
| function of capillary pressure [rlp(h)].                                                                                                                            |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP1                                                                   | real    | Residual liquid saturation, matrix rock material.                                 |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP2                                                                   | real    | Maximum liquid saturation, matrix rock material.                                  |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP3                                                                   | real    | Inverse of air entry head, :math:`\alpha_G` (1/m) [note                           |
|                                                                       |         | some data is given in (1/Pa) convert using                                        |
|                                                                       |         | pressure = :math:`\rho g \Delta h`], matrix rock material.                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP4                                                                   | real    | Power n in van Genuchten formula, matrix rock material.                           |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP5                                                                   | real    | Low saturation fitting parameter, matrix rock                                     |
|                                                                       |         | material, multiple of cutoff capillary pressure                                   |
|                                                                       |         | assigned as maximum capillary pressure. If                                        |
|                                                                       |         | RP5 < 0 then a linear fit from the cutoff                                         |
|                                                                       |         | saturation (RP6) is used. The slope of the                                        |
|                                                                       |         | cutoff saturation is used to extend the                                           |
|                                                                       |         | function to saturation = 0.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
|                                                                       |         | If RP5 = 0, a cubic fit is used. The slope                                        |
|                                                                       |         | at the cutoff saturation is matched and the                                       |
|                                                                       |         | conditions  :math:`\frac{\partial}{\partial S}Pcap = 0`                           |
|                                                                       |         | and  :math:`\frac{\partial^2}{\partial S}Pcap = 0` are                            |
|                                                                       |         | forced at  :math:`S = 0`. If RP5 > 0, a multiple                                  |
|                                                                       |         | of the value of the capillary pressure at the                                     |
|                                                                       |         | cutoff saturation, :math:`RP5\cdot Pcap(S_{cutoff})`                              |
|                                                                       |         | is forced at :math:`S = 0`.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP6                                                                   | real    | Cutoff saturation used in fits described for RP5,                                 |
|                                                                       |         | must be greater than RP1.                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP7                                                                   | real    | Residual liquid saturation, fracture material.                                    |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP8                                                                   | real    | Maximum liquid saturation, fracture material.                                     |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP9                                                                   | real    | Inverse of air entry pressure, :math:`\alpha_G` (1/m) [note                       |
|                                                                       |         | some data is given in(1/Pa) convert using                                         |
|                                                                       |         | pressure = :math:`\rho g \Delta h`], fracture material.                           |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP10                                                                  | real    | Power n in van Genuchten formula, fracture material.                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP11                                                                  | real    | Low saturation fitting parameter, fracture material,                              |
|                                                                       |         | multiple of cutoff capillary pressure assigned as                                 |
|                                                                       |         | maximum capillary pressure. If RP11 < 0 then a linear                             |
|                                                                       |         | fit from the cutoff saturation (RP12) is used. The                                |
|                                                                       |         | slope of the cutoff saturation is used to extend                                  |
|                                                                       |         | the function to saturation = 0.If RP11 = 0, a cubic                               |
|                                                                       |         | fit is used. The slope at the cutoff saturation is                                |
|                                                                       |         | matched and the conditions                                                        |
|                                                                       |         | :math:`\frac{\partial}{\partial S}Pcap = 0`                                       |
|                                                                       |         | and :math:`\frac{\partial^2}{\partial S}Pcap = 0`                                 |
|                                                                       |         | are forced at :math:`S = 0`. If RP11 > 0,                                         |
|                                                                       |         | a multiple of the value of the capillary                                          |
|                                                                       |         | pressure at the cutoff saturation,                                                |
|                                                                       |         | :math:`RP11\cdot Pcap(S_{cutoff})` is forced at                                   |
|                                                                       |         | :math:`S = 0`.                                                                    |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP12                                                                  | real    | Cutoff saturation used in fits described for RP11,                                |
|                                                                       |         | must be greater than RP7.                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP13                                                                  | real    | Fracture permeability (:math:`m^2`). This is the permeability                     |
|                                                                       |         | of the individual fractures. The bulk permeability                                |
|                                                                       |         | of the fracture continuum is :math:`RP13\times RP15`.                             |
|                                                                       |         | Can be made anisotropic with macro FPER.                                          |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP14                                                                  | real    | Matrix rock saturated permeability (:math:`m^2`). Can be made                     |
|                                                                       |         | anisotropic with macro FPER.                                                      |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP15                                                                  | real    | Fracture volume fraction. Is equal to the                                         |
|                                                                       |         | fracture aperture divided by the fracture                                         |
|                                                                       |         | spacing (with same units). Sometimes called fracture porosity.                    |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 5**:  IRLP(i) = 5, van Genuchten relative permeability, van                                                                                                 |
| Genuchten capillary pressure (6 parameters required). This model                                                                                                    |
| and its input are the same as for Model 3 except that                                                                                                               |
| permeabilities are represented as a function of saturation                                                                                                          |
| [rlp(S)] rather than capillary pressure.                                                                                                                            |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 6**:  IRLP(i) = 6, van Genuchten relative permeability, van                                                                                                 |
| Genuchten capillary pressure, effective continuum (15 parameters                                                                                                    |
| required). This model and its input are the same as for Model 4                                                                                                     |
| except that permeabilities are represented as a function of                                                                                                         |
| saturation [rlp(S)] rather than capillary pressure.                                                                                                                 |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 7**:  IRLP(i) = 7, van Genuchten relative permeability, van                                                                                                 |
| Genuchten capillary pressure, effective continuum with special                                                                                                      |
| fracture interaction term (16 parameters required). This model                                                                                                      |
| and its input are the same as for Model 6 except that the an                                                                                                        |
| additional term is included which represents the fracture-matrix                                                                                                    |
| interaction.                                                                                                                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP16                                                                  | real    | Fracture-matrix interaction term. If RP16 ≤ 0,                                    |
|                                                                       |         | then an additional multiplying term equal to the                                  |
|                                                                       |         | relative permeability is applied to the                                           |
|                                                                       |         | fracture-matrix interaction term for dual                                         |
|                                                                       |         | permeability problems. If RP16 > 0, then an                                       |
|                                                                       |         | additional multiplying term equal to :math:`sl**RP16`                             |
|                                                                       |         | and :math:`(1.-sl)**RP16` is applied to the                                       |
|                                                                       |         | fracture-matrix interaction terms for the                                         |
|                                                                       |         | liquid and vapor phases, respectively, for dual                                   |
|                                                                       |         | permeability problems. Here, :math:`sl`                                           |
|                                                                       |         | is the value of saturation at the given node.                                     |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| **Model 10**: IRLP(i) = 10, linear relative permeability with minimum                                                                                               |
| relative permeability values, linear capillary pressure (8                                                                                                          |
| parameters required).                                                                                                                                               |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP1                                                                   | real    | Residual liquid saturation.                                                       |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP2                                                                   | real    | Residual vapor saturation.                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP3                                                                   | real    | Maximum liquid saturation.                                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP4                                                                   | real    | Maximum vapor saturation.                                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP5                                                                   | real    | Minimum liquid permeability (:math:`m^2`).                                        |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP6                                                                   | real    | Minimum vapor permeability (:math:`m^2`).                                         |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP7                                                                   | real    | Capillary pressure at zero saturation (:math:`MPa`).                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+
| RP8                                                                   | real    | Saturation at which capillary pressure goes to zero.                              |
+-----------------------------------------------------------------------+---------+-----------------------------------------------------------------------------------+


The following is an example of ``rlp``.
In this example, Corey type relative permeability is specified,
with residual liquid saturation of 0.3, residual vapor saturation of 0.1,
a base capillary pressure of 2 MPa, and capillary pressure goes to zero at a
saturation of 1. This model is assigned to nodes numbered 1 through 140.


+-----+-----+-----+-----+----+
| rlp |     |     |     |    |
+-----+-----+-----+-----+----+
| 2   | 0.3 | 0.1 | 2.0 | 1. |
+-----+-----+-----+-----+----+
|     |     |     |     |    |
+-----+-----+-----+-----+----+
| 1   | 140 | 1   | 1   |    |
+-----+-----+-----+-----+----+
|     |     |     |     |    |
+-----+-----+-----+-----+----+


========
``vcon``
========

Variable thermal conductivity information.

1. Thermal conductivity for intact salt: :math:`\lambda_{IS} = \lambda_{IS,300}(300/T)^{\gamma_1}`

2. Thermal conductivity for crushed salt: :math:`\lambda_{CS} = \lambda_{CS,300}(300/T)^{\gamma_1}`

where:

:math:`\lambda_{CS,300} = 1.08 \left( \alpha_4 \phi^4 + \alpha_3 \phi^3 + \alpha_2 \phi^2 + \alpha_1 \phi + \alpha_0 \right)`

Parameters related to thermal conductivity are
:math:`\lambda_{IS,300},\gamma_1,\gamma_2,\alpha_4,\alpha_3,\alpha_2,\alpha_1,\alpha_0` and :math:`\phi`.
An additional parameter is the specific heat of salt.

* Group 1 - IVCON(I), VC1F(I), VC2F(I), VC3F(I)

* Group 2 - JA, JB, JC, IVCND

The parameter (I) is incremented each time Group 1 is read. Group 2 lines will refer to this parameter. Group 1 is ended with a blank line. 

+----------------+---------+-----------------------------------------------------------------------------+
| Input Variable | Format  | Description                                                                 |
+================+=========+=============================================================================+
| IVCON(i)       | integer | | Model type for ith conductivity model.                                    |
|                |         | | IVCON(i) = 1, linear variation of thermal conductivity with temperature.  |
|                |         | | IVCON(i) = 2, square root variation of thermal conductivity with liquid   |
|                |         |   saturation.                                                               |
+----------------+---------+-----------------------------------------------------------------------------+
| VC1F(i)        | real    | Reference temperature (oC) for IVCON(i) = 1.                                |
|                |         | Conductivity :math:`\left( \frac{W}{m \cdot K} \right)` at liquid           |
|                |         | saturation = 1 for IVCON(i) = 2.                                            |
+----------------+---------+-----------------------------------------------------------------------------+
| VC2F(i)        | real    | Reference conductivity ( :math:`^oC` ) for IVCON(i) = 1. Conductivity       |
|                |         | :math:`\left( \frac{W}{m \cdot K} \right)` at liquid saturation = 0 for     |
|                |         | IVCON(i) = 2.                                                               |
+----------------+---------+-----------------------------------------------------------------------------+
| VC3F(i)        | real    | Change in conductivity with respect to temperature for IVCON(i) = 1.        |
|                |         | Not used for IVCON(i) = 2.                                                  |
+----------------+---------+-----------------------------------------------------------------------------+
| IVCND          | integer | Number referring to the sequence of models read in Group 1.                 |
|                |         | The default is 1.                                                           |
+----------------+---------+-----------------------------------------------------------------------------+


The following are examples of ``vcon``. In the first example,
a linear conductivity model is defined and applied at each node.
The reference temperature is :math:`20^o C`, the reference
conductivity is :math:`1\frac{W}{m \cdot K}`, and the change
in conductivity with temperature is 0.01.

+-------+--------+-------+---------+
| vcon  |        |       |         |
+=======+========+=======+=========+
| 1     |  20.0  |  1.0  |  1.e-2  |
+-------+--------+-------+---------+
|       |        |       |         |
+-------+--------+-------+---------+
| 1     | 0      | 0     | 1       |
+-------+--------+-------+---------+
|       |        |       |         |
+-------+--------+-------+---------+

For the second example, three models are defined for the entire domain.
Model 1 defines the constant thermal conductivity of 16.26 W/m K at 26.85 °C
(=300 °K) for a stainless steel canister (zone 1).

Model 2 defines all parameters for a crushed salt (zone 2).

Model 3 defines reference thermal conductivity of 5.4 W/m K at 26.85 °C
(=300 °K) and exponent 1.14 for the intact salt in the rest of the domain. 

+------+-------+-------+-------+-------+-------+-----+-----+------+
| vcon |       |       |       |       |       |     |     |      |
+======+=======+=======+=======+=======+=======+=====+=====+======+
| 1    | 26.85 | 16.26 | 0.    |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+
| 4    | 26.85 | 1.08  | 270.0 | 370.0 | 136.0 | 1.5 | 5.0 | 1.14 |
+------+-------+-------+-------+-------+-------+-----+-----+------+
| 3    | 26.85 | 5.4   | 1.14  |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+
|      |       |       |       |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+
| 1    | 0     | 0     | 3     |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+
| -1   | 0     | 0     | 1     |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+
| -2   | 0     | 0     | 2     |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+
|      |       |       |       |       |       |     |     |      |
+------+-------+-------+-------+-------+-------+-----+-----+------+



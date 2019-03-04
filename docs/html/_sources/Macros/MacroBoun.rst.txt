========
``boun`` 
========

.. note::

  Either boun or flow is required for a flow problem.

Implement boundary conditions and sources or sinks. Input may be time dependent and cyclic. Time step sizes may also be adjusted. 

* Group 1 - ``KEYWORD``

The Group 1 KEYWORD 'model', which starts each model sequence, is followed immediately by a Group 2 KEYWORD of 'ti', 'ti_linear', 'cy' or 'cy_linear'. 

* Group 2 - ``KEYWORD``

* Group 3 - ``NTIMES, TIME(I), I=1,NTIMES``

The Group 4 KEYWORDs define the various boundary condition parameters being entered. These KEYWORDs and associated data, Group 5, are repeated as needed for each model. Note that some keywords do not have associated variables.

* Group 4 - ``KEYWORD``

* Group 5 - ``VARIABLE(I), I=1,NTIMES``

Additional models are entered by beginning again with Group 1. The ``MODEL_NUMBER`` is incremented each time a new model is read, and is used to assign boundary conditions to specified nodes or zones in Group 6. After all models have been entered, the section is terminated with KEYWORD 'end' or a blank line. 

+----------------+-------------+------------------------------------------------------------------------------------------------------------+
| Input Variable | Format      | Description                                                                                                |
+================+=============+============================================================================================================+
| KEYWORD        | character*4 | | Keyword specifying a model designation, time for boundary condition                                      |
|                |             | | or source/sink changes, or actual variable or source change.                                             |
|                |             | | Keywords, which must be entered starting in the first column, are:                                       |
|                |             | |   ``model`` - new model definition to follow                                                             |
|                |             | | *Note*: Descriptive text, such as the model number, may be appended after                                |
|                |             | |       the 'model' keyword as long as it is contained on a single line, and                               |
|                |             | |       begins after column four.                                                                          |
|                |             | | ``ti`` - time sequence for changes to follow (days). The 'ti' keyword results                            |
|                |             | |          in step function changes in boundary condition                                                  |
|                |             | |          VARIABLE(i) at each TIME(i).                                                                    |
|                |             | | ``ti_linear`` - time sequence for changes to follow (days). The 'ti_linear'                              |
|                |             | |                 keyword will apply a boundary condition that changes linearly                            |
|                |             | |                 with time. This option does not impose any control on time                               |
|                |             | |                 step size, so it is possible that a single time step can span                            |
|                |             | |                 an entire time interval and the linear change will not be seen.                          |
|                |             | |                 If time step size control is important it should be imposed                              |
|                |             | |                 in the time or ctrl macros.                                                              |
|                |             | | ``cy`` - cyclic time sequence for changes to follow (days). As with the 'ti'                             |
|                |             | |          keyword, boundary condition changes are step functions.                                         |
|                |             | | ``cy_linear`` - cyclic time sequence for changes to follow (days). As with the                           |
|                |             | |               ``ti_linear`` keyword, boundary condition changes linearly with time.                      |
|                |             | | ``sec`` - Time sequence input is in seconds.                                                             |
|                |             | | ``min`` - Time sequence input is in minutes.                                                             |
|                |             | | ``day`` - Time sequence input is in days. (Default)                                                      |
|                |             | | ``year`` - Time sequence input is in years.                                                              |
|                |             | | *Note*: The keywords, 'ti', 'ti_linear', 'cy' and 'cy_linear', require the                               |
|                |             | |       time to start at 0.0. This provides the initial boundary and                                       |
|                |             | |       source/sink information. If the input for 'ti' or 'cy' does not                                    |
|                |             | |       start at 0.0 time the code assumes boundary conditions and                                         |
|                |             | |       source/sinks are 0.0 at time 0.0. The 'cy' keyword involves a                                      |
|                |             | |       cyclic changing of conditions. In our procedure the cycle ends at                                  |
|                |             | |       the last specified time. Thus the code reverts to the first                                        |
|                |             | |       specified time values. Because of this, the boundary conditions                                    |
|                |             | |       and source/sinks for the last time change are always set to the first                              |
|                |             | |       time values.The default units for time input is days. Input in seconds,                            |
|                |             | |       minutes or years is converted to days. Time input units have no                                    |
|                |             | |       associated variable input.               |                                                         |
|                |             | | ``tran`` - Keyword to indicate boundary conditions will not be invoked until                             |
|                |             | |            the steady state portion of the simulation is completed and a                                 |
|                |             | |            transient run is initiated. See macro stea for more details.                                  |
|                |             | | ``sa`` - air source sequence for changes to follow (kg/s)                                                |
|                |             | | ``sw`` - water source sequence for changes to follow (kg/s)                                              |
|                |             | | ``swf`` - source water factor sequence for changes to follow (multiplier                                 |
|                |             | |           for existing mass flow rates)                                                                  |
|                |             | | ``se`` - enthalpy source sequence for changes to follow (MW)                                             |
|                |             | | ``sf`` - water seepage face sequence with pressures for changes to follow (MPa)                          |
|                |             | | ``sfh``  - water seepage face sequence with heads for changes to follow (m)                              |
|                |             | | ``fd``  - water drainage area sequence for changes to follow (m2)                                        |
|                |             | | ``dsa`` - distributed air source sequence for changes to follow (kg/s)                                   |
|                |             | | ``dsw`` - distributed water source sequence for changes to follow(kg/s)                                  |
|                |             | | ``dse`` - distributed enthalpy source sequence for changes to follow (MW)                                |
|                |             | | *Note*: A distributed source (keywords 'dsa', 'dsw', and 'dse') is a source                              |
|                |             | |         term divided over a group of nodes or a zone proportional to the nodal volume.                   |
|                |             | | ``wgt`` - Distributed source is weighted using the nodal control volume.                                 |
|                |             | | ``wgtx`` - Distributed source is weighted using nodal area = control volume / x length scale.            |
|                |             | | ``wgty`` - Distributed source is weighted using nodal area = control volume / y length scale.            |
|                |             | | ``wgtz`` - Distributed source is weighted using nodal area = control volume / z length scale.            |
|                |             | | ``wgtp`` - Distributed source is weighted using nodal control volume * permeability.                     |
|                |             | | ``wgtpx`` - Distributed source is weighted using nodal control volume * permeability / x length scale.   |
|                |             | | ``wgtpy`` - Distributed source is weighted using nodal control volume * permeability / y length scale.   |
|                |             | | ``wgtpz`` - Distributed source is weighted using nodal volume * permeability / z length scale.           |
|                |             | | *Note*: The length scale term is the dimension of the control volume bounding box, xmax-xmin,            |
|                |             | |         ymax-ymin, zmax-zmin, depending upon the suffix, x,y,z. This option is useful when one           |
|                |             | |         wants to apply a distributed source on a mesh with variable size mesh cells and would            |
|                |             | |         like the source percentage to be allocated based on surface area of each node.                   |
|                |             | | ``wgtr`` - Distributed source is weighted using nodal volume * permeability * relative permeability      |
|                |             | | ``wgtu`` - Distributed source is weighted with user specified values (See macro wgtu).                   |
|                |             | | ``wgww`` - Distributed source is weighted using nodal volume * permeability * relative                   |
|                |             | |            permeability * exponentially weigthed distance from pump                                      |
|                |             | | *Note*: The distributed source weighting options have no associated variable input.                      |
|                |             | | ``s`` - fixed saturation sequence for changes to follow                                                  |
|                |             | | ``hd`` - fixed hydraulic head sequence for changes to follow (m)                                         |
|                |             | | ``pw`` - fixed water pressure sequence for changes to follow (MPa)                                       |
|                |             | | ``pa`` - fixed air pressure sequence for changes to follow (MPa)                                         |
|                |             | | ``hdo`` - fixed hydraulic head sequence for changes to follow (m) (constrained to outflow only)          |
|                |             | | ``pwo`` - fixed water pressure sequence for changes to follow (MPa) (constrained to outflow only)        |
|                |             | | ``pao`` - fixed air pressure sequence for changes to follow (MPa) (constrained to outflow only)          |
|                |             | | ``en`` - fixed enthalpy sequence for changes to follow (MW)                                              |
|                |             | | ``t`` - fixed temperature sequence for changes to follow (oC)                                            |
|                |             | | ``h`` - fixed humidity sequence for changes to follow (must be used with van Genuchten relative          |
|                |             | |         permeability model)                                                                              |
|                |             | | ``ft`` - fixed flowing temperature sequence for change to follow (oC). By flowing temperature            |
|                |             | |          we mean the temperature of the inflow stream for a specified source. If no source               |
|                |             | |          inflow occurs where this condition is applied, it will be ignored.                              |
|                |             | | ``kx`` - fixed X permeability sequence for changes to follow (m2)                                        |
|                |             | | ``ky`` - fixed Y permeability sequence for changes to follow (m2)                                        |
|                |             | | ``kz`` - fixed Z permeability sequence for changes to follow (m2)                                        |
|                |             | | ``if`` - impedance factor for use with fixed water pressure boundary condition.                          |
|                |             | |          If left out the impedance factor will be set to the volume of the grid cell.                    |
|                |             | | ``si``  - initial value saturation sequence for changes to follow                                        |
|                |             | | ``pai``  - initial value air pressure sequence for changes to follow (MPa)                               |
|                |             | | ``pwi``  - initial value water pressure sequence for changes to follow (MPa)                             |
|                |             | | ``tmi``  - initial value temperature sequence for changes to follow (oC)                                 |
|                |             | | *Note*: The keywords 'si', 'pai', 'pwi', and 'tmi' refer to changes for a variable                       |
|                |             | |         that is NOT fixed. They are similar to specifying initial conditions in that                     |
|                |             | |         regard but may be changed according to a time sequence. At present these 4                       |
|                |             | |         keywords only work with isothermal air-water calculations.                                       |
|                |             | | ``chmo`` - model number sequence for changes to follow                                                   |
|                |             | | ``ts`` - timestep sequence for changes to follow (days)                                                  |
|                |             | | ``end`` - signifies end of keyword input, a blank line will also work.                                   |
+----------------+-------------+------------------------------------------------------------------------------------------------------------+
| NTIMES         | integer     | Number of time changes for boundary condition or source/sink specification.                                |
+----------------+-------------+------------------------------------------------------------------------------------------------------------+
| TIME           | real        | NTIMES times for changes in boundary conditions or source/sinks.                                           |
+----------------+-------------+------------------------------------------------------------------------------------------------------------+
| VARIABLE       | real        | NTIMES new values for boundary conditions or source/sinks.                                                 |
+----------------+-------------+------------------------------------------------------------------------------------------------------------+
| MODEL_NUMBER   | integer     | | Boundary condition model to be associated with designated nodes or zones                                 |
|                |             | | (the number corresponds to the numerical order in which the models were                                  |
|                |             | | input, i.e., beginning with KEYWORD ``model``)                                                           |
+----------------+-------------+------------------------------------------------------------------------------------------------------------+


The following is an example of ``boun``. In this example two models are defined.
The first model cyclically changes the water source in a 1.e05 day cycle, i.e.,
the ``cy`` keyword entry shows that the time cycle ends at 1.e05 days and at this
time the cycle reverts to 0.0 days. Note that the water source at 1.e05 days
equals that at 0.0 days. Also in model 1 the flowing temperature was alternated
between 20oC and 50oC. The second model uses a time change that occurs at 1.e20 days.
This effectively removes any time variance from model 2. Model 2 has a fixed
water pressure and flowing temperature condition. The models are applied at
nodes 26 and 27 in the last two lines. It should be noted that the model numbers
included in the example (following KEYWORD ``model``) are not part of the required
input but are descriptive text used to enhance readability of the macro.

+---------+--------+--------+--------+--------+--------+
|                         boun                         |
+=========+========+========+========+========+========+
| model 1 |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
| cy      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |   4    |   0.0  |   1.e1 |   1.e2 |  1.e5  |
+---------+--------+--------+--------+--------+--------+
| sw      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         | -1.e-4 | -1.e-5 | -1.e-3 | -1.e-4 |        |
+---------+--------+--------+--------+--------+--------+
| ft      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |  20.0  |  50.0  | 50.0   |  20.0  |        |
+---------+--------+--------+--------+--------+--------+
| model 2 |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
| ti      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |    2   |  0.0   | 1.e20  |        |        |
+---------+--------+--------+--------+--------+--------+
| pw      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |   0.1  |  0.1   |        |        |        |
+---------+--------+--------+--------+--------+--------+
| ft      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |  20.0  | 20.0   |        |        |        |
+---------+--------+--------+--------+--------+--------+
| end     |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |   26   |   26   |   1    |   1    |        |
+---------+--------+--------+--------+--------+--------+
|         |   27   |   27   |   1    |   2    |        |
+---------+--------+--------+--------+--------+--------+


In the second example, a distributed water source is used to model a zone where
production is turned on and off. Keyword ``kz`` is used to specify a higher
permeability when production is occurring. The ``ts`` keyword is used to reset the
time step to 1.0 days, at the beginning of each time interval. The model is
applied to zone 100 in the last line.


+---------+--------+--------+--------+--------+--------+
|                         boun                         |
+=========+========+========+========+========+========+
| model 1 |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
| ti      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         |   4    |   0.0  | 91.325 | 182.62 | 273.93 |
+---------+--------+--------+--------+--------+--------+
| ts      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         | 1.0    | 1.0    | 1.0    | 1.0    |        |
+---------+--------+--------+--------+--------+--------+
| dsw     |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         | 29.248 | 0.0    | 29.248 | 0.     |        |
+---------+--------+--------+--------+--------+--------+
| kz      |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         | 8e-12  | 2e-12  | 8e-12  | 2e-12  |        |
+---------+--------+--------+--------+--------+--------+
| end     |        |        |        |        |        |
+---------+--------+--------+--------+--------+--------+
|         | -100   | 0      | 0      | 1      |        |
+---------+--------+--------+--------+--------+--------+


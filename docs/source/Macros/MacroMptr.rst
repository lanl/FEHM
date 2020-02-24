========
``mptr``
========

Multiple species ingrowth particle tracking. Note that data for each numbered group must be input. The other input is optional.

* Group 1 -	NSPECI, MAXLAYERS, MAX_PARTICLES, RIPFEHM, MAX1D

* Group 2 -	POUT, PRNT_RST

  - or when PRNT_RST ≥ 20, selected output parameters

* Group 2 -	POUT, PRNT_RST, PRNT_VARNUM ( 1 . . . 6)

  - Optional keyword "tcurve" is input to indicate that transfer function curves should be input to model matrix diffusion. It is followed by NUMPARAMS and TFILENAME.

    + KEYWORD

    + NUMPARAMS, FFMAX

    + TFILENAME

  - Optional keyword "zptr" designates zones for breakthrough curves will be defined. It is followed by IPZONE and IDZONE.

    + KEYWORD ‘zptr'

    + IPZONE

    + IDZONE(I) I = 1 to IPZONE

* Group 3 -	RSEED, RSEED_RELEASE

  - Optional keyword "wtri" is input to indicate a water table rise calculation should be performed. It is followed by WATER_TABLE. For GoldSim the water table rise calculation is controlled by passing a new water table elevation to the code during the simulation and the keyword is not required.

    + KEYWORD ‘wtri'

    + WATER_TABLE

* Group 4 -	DAYCS, DAYCF, DAYHF, DAYHS

  - An optional, flexible input structure involving the assignment of transport parameters is implemented in the particle tracking input to allow multiple realization simulations to use different parameters for each realization. The user invokes this option using the keyword "file" before Group 5, followed by the name of the file that the transport parameters reside in. The applicable transport parameters are defined in Group 5 and Group 9.

    + KEYWORD ‘file'

    + PFILENAME

  - The structure of the alternate parameter file is:

    + NINPUTS	

    + PTRPARAM(I) I=1 to NINPUTS

    + .

    + .

    + .

  - [a line of parameters is present for each realization]

The method for assigning a given value of the particle tracking parameter (PTRPARAM) to a specific transport parameter, defined in Group 5 or Group 9, is discussed below. There are an arbitrary number of input lines each representing a given realization of parameter values. In a multiple-realization scenario, the code enters the input file for each realization, and for this input, reads down the corresponding number of lines to obtain the parameters for that realization. For example, for realization number 10, the code reads down to the 10th line of data (line 11 in the file) and uses those parameter values.

Once these parameters are read in for a given realization, they must be assigned to specific transport parameters. This is done in the following way in Group 5 or Group 9. If any of the inputs other than TRANSFLAG are negative, the code takes the absolute value of the number and interprets it as the column number from which to assign the transport parameter. For example, if DIFFMFL = -5, then the diffusion coefficient is the fifth input number in the PTRPARAM array. In this way, any of the transport parameter inputs can be assigned through the alternate input file rather than the input line in mptr. It should be noted that for the colloid diversity model, only K_REV need be negative to indicate values should be read from the parameter file, if K_REV is negative then all five parameters are read from the file, otherwise the equation parameters will be read from the mptr macro. This is to accommodate the fact that the SLOPE_KF may have negative values.

Group 5 is used to define models in which identical transport parameters are assumed to apply. Group 5 data are read until a blank line is encountered. The model number ITRC is incremented by 1 each time a line is read. Model parameters defined in Group 5 are assigned to nodes or zones using Group6. 

Optional keyword "afm" indicates the Active Fracture Model input for determining fracture spacing should be used. Optional keyword "dfree" is input to indicate that a free water diffusion coefficient and tortuosity will be entered instead of the molecular diffusion coefficient.

	KEYWORD ‘afm'

	KEYWORD ‘dfre'

* Group 5 -	TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), APERTUR(ITRC), MATRIX_POR(ITRC)

  - or when ‘afm' is implemented:

  - Group 5 -	TCLX(ITRC), TCLY(ITRC), TCLZ(ITRC), APERTUR(ITRC), MATRIX_POR(ITRC), SRESIDUAL(ITRC), GAMMA_AFM(ITRC)

* Group 6 -	JA, JB, JC, ITRC (JA, JB, JC - defined on `See JA, JB, JC, PROP1, PROP2, . . . <Macro20058.html>`_)

The following groups (Group 7 - 12) are repeated for each species.

* Group 7 -	ITH_SPECI, TRAK_TYPE, HALF_LIFE, IDAUGHTER, CONFACTOR, NEWCONFACTOR, CONFTIME, GMOL, P_FRACTION, ASTEP, CFRACTION

Optional keyword "size" is input to indicate that the colloid size distribution model option is being used. It is followed by PART_SIZE and PROBSIZE.

	KEYWORD ‘size'

	PART_SIZE(I), PROBSIZE(I) - an arbitrary numbers of lines of input, terminated by a blank line.

Optional keyword "dive" is input to indicate that the colloid diversity model is being used. It is followed by FLAG_COL_DAUGHTER, optional keyword "file" and the name of the file containing the CDF table or equation data (a description of the format for the file is provided with the second mptr example below), the TPRPFLAG with optional SIMNUM, optional CDF equation parameters (when "file" is not used), and keyword "irreversible" or "reversible" with FLAG_LOG.

+---------------------------------+----+--------------------------------------------+
| KEYWORD 'dive'                  | or | KEYWORD 'dive                              |
+---------------------------------+----+--------------------------------------------+
| FLAG\_COL\_DAUGHTER             |    | FLAG\_COL\_DAUGHTER                        |
+---------------------------------+----+--------------------------------------------+
| KEYWORD 'file'                  |    | TPRPFLAG                                   |
+---------------------------------+----+--------------------------------------------+
| CDFFILENAME                     |    | or                                         |
+---------------------------------+----+--------------------------------------------+
| TPRPFLAG                        |    | TPRPFLAG, SIMNUM                           |
+---------------------------------+----+--------------------------------------------+
| or                              |    | K\_REV, R\_MIN,R\_MAX, SLOPE\_KF, CINT\_KF |
+---------------------------------+----+--------------------------------------------+
| TPRPFLAG, SIMNUM                |    | KEYWORD ‘irreversible’                     |
+---------------------------------+----+--------------------------------------------+
| KEYWORD ‘irreversible’          |    | or                                         |
+---------------------------------+----+--------------------------------------------+
| or                              |    | KEYWORD ‘reversible’, FLAG_LOG             |
+---------------------------------+----+--------------------------------------------+
| KEYWORD ‘reversible’, FLAG\_LOG |    |                                            |
+---------------------------------+----+--------------------------------------------+

Note that optional KEYWORDs "size" and "dive" are only used when colloid transport is enabled.

* Group 8 - LAYERS

* Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, DIFFMFL

or for simulations using "dfree":

* Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, H2O_DIFF, TORT_DIFF

or for simulations with colloid (``TRANSFLAG < 0``):

* Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, DIFFMFL, KCOLL, RCOLL, FCOLL

or for simulations with colloid using "dfree":

* Group 9 - LAYER_I, TRANSFLAG, KD, RD_FRAC, H2O_DIFF, TORT_DIFF, KCOLL, RCOLL, FCOLL

* Group 10 - NS

* Group 11 - JA, JB, JC, TMPCNSK

Note that because the number of source terms is controlled by the value entered for NS, Group 11 input is not terminated with a blank line.

* Group 12 - PINMASS, T1SK, T2SK 

For transient source terms, Group 12 is repeated for each time interval and terminated with a blank line. Groups 11 and 12 are repeated for each source term (from 1 to NS).

For decay-ingrowth calculations, when the particle injection period is too small (for example, 1.E-4 days) compared to the half-life of the radionuclides and the half-life is large (for example 1.E+9 days), numerical errors in the decay-ingrowth calculation may arise due to truncation error. To get better accuracy, the user should try to increase the length of the injection period.

For particle tracking simulations using the transfer function method (see `See Transfer function curve data input file <Macro49660.html>`_ for input file format), it is sometimes desirable to identify the parameter ranges over which the two- and three-parameter type curves are accessed, so that an assessment can be made regarding the density of transfer function curves in a given part of the parameter space. If the flag output_flag in the transfer function file is set to "out", the code writes the real*8 array param_density to the *.out file in the following format:

For regular parameter spacings, the output is:


.. code::

   i = 1, nump1
      j = 1, nump2
          k = nump3
   
              write(iout.*) param_density(i,j,k)
   
          end do
      end do
   end do


For two-parameter models, only the i and j loops are used. The value of param_density is the number of times any particle passes through any node at those values of the parameters. This allows the user to identify regions in which a greater density of transfer functions may be required. For the option 'free' in which there is no structure to the parameter grid used for the transfer function curves, nump1 is the total number of curves, and nump2 and nump3 are equal to 1. 

+-------------------+--------------+-------------------------------------------------------------------------------------------+
| Input Variable    | Format       | Description                                                                               |
+===================+==============+===========================================================================================+
| NSPECI            | integer      | Number of species in the simulation.                                                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| MAXLAYERS         | integer      | Maximum number of property layers in the model.                                           |
|                   |              | The actual number of layers used in the model must be ≤ MAXLAYERS.                        |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| MAX_PARTICLES     | integer      | Maximum number of particles used for individual species.                                  |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| RIPFEHM           | integer      | | Flag to indicate if simulation is coupled with GoldSim.                                 |
|                   |              | |   RIPFEHM = 0, FEHM standalone simulation                                               |
|                   |              | |   RIPFEHM = 1, GoldSim-FEHM coupling simulation                                         |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| MAX1D             | integer      | Maximum 1-D array size for holding particle tracking information for                      |
|                   |              | all simulated species. The value of MAX1D depends on number of species,                   |
|                   |              | number of time steps, number of radionuclide release bins, number of                      |
|                   |              | species involved in ingrowth, and the length of the decay-ingrowth chain.                 |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| POUT              | integer      | | Flag to specify the concentration output format:                                        |
|                   |              | | 1 -  Concentrations computed as number of particles per unit total volume               |
|                   |              |   (rock and fluid)                                                                        |
|                   |              | | 2 -  Concentrations computed as number of particles per unit fluid volume               |
|                   |              |   (the fluid is liquid for TRAK_TYPE = 1 and gas for TRAK_TYPE = 2).                      |
|                   |              | | 3 -  Concentrations computed as number of particles at a given node point.              |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| PRNT_RST          | integer      | | Flag to specify whether particle information is written to the ".fin",                  |
|                   |              |   ".ptrk_fin", or ".ptrk" files:                                                          |
|                   |              | | If PRNT_RST = 0, Particle information is not written to the output files.               |
|                   |              | | If PRNT_RST = 1, 11, 21, 31, 41 All particle information necessary for a                |
|                   |              |   restart is written to the ".fin" file.                                                  |
|                   |              | | If PRNT_RST = -1, -11, -21, -31, -41 Only particle positions and ages are               |
|                   |              |   written to the ".fin" file.                                                             |
|                   |              | | If ABS (PRNT_RST) = 2, 12, 22, 32, 42 Mass flux values are written to the               |
|                   |              |   ".fin" file followed by particle information.                                           |
|                   |              | | If 10 ≤ ABS(PRNT_RST) < 30 Particle exit locations and count are written                |
|                   |              |   to the ".ptrk_fin" file.                                                                |
|                   |              | | If ABS(PRNT_RST) ≥ 20 Cumulative particle counts versus time are written                |
|                   |              |   to the ".ptrk" file, for variables specified by PRNT_VARNUM (the default                |
|                   |              |   is to output all variables).                                                            |
|                   |              | | If ABS(PRNT_RST) ≥ 40, Cumulative mass output from a FEHM/GoldSim coupled               |
|                   |              |   simulation will be written to file ``FEHM_GSM_Mass_balance.txt``. Note that to          |
|                   |              |   track cumulative mass an additional array of size ``maxparticles*nspeci`` must          |
|                   |              |   be allocated so caution should be used when specifying this option to ensure            |
|                   |              |   sufficient system memory is available.                                                  |
|                   |              | |                                                                                         |
|                   |              | | When particle tracking data or mass fluxes are written to the ``.fin`` file,            |
|                   |              |   the arrays are written after all of the heat and mass simulation information.           |
|                   |              |   The mass fluxes can be read into the code in a subsequent ptrk or mptr simulation       |
|                   |              |   and the code can simulate transport on this steady state flow field (see macro          |
|                   |              |   ``rflo``).The particle information written is sufficient to perform a restart of the    |
|                   |              |   particle tracking simulation and to post-process the data to compile statistics         |
|                   |              |   on the particle tracking run. However, for a large number of particles, this            |
|                   |              |   file can become quite large, so particle tracking information should only be            |
|                   |              |   written when necessary. Thus, 0 should be used for ``PRNT_RST`` unless restarting       |
|                   |              |   or post-processing to obtain particle statistics is required. Selecting the             |
|                   |              |   "-" options allows a subset of the full set of information needed for a                 |
|                   |              |   restart (particle positions and ages) to be written. Restart runs that use              |
|                   |              |   this file as input will only be approximate, since the particle is assumed              |
|                   |              |   to have just entered its current cell. For restart runs, ``PRNT_RST = 1`` is            |
|                   |              |   preferred, while ``PRNT_RST = -1`` is appropriate for output of particle                |
|                   |              |   statistics for post- processing.                                                        |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| PRNT_VARNUM       | integer      | | A list of integers specifying which particle counts should be output. For each          |
|                   |              |   value entered ``PRNT_VAR(PRNT_VARNUM)`` is set to true. If no values are entered        |
|                   |              |   the default is to print all variables.                                                  |
|                   |              | | 1 – Number of particles that have entered the system                                    |
|                   |              | | 2 – Number of particles currently in the system                                         |
|                   |              | | 3 – Number of particles that have left the system                                       |
|                   |              | | 4 – Number of particles that have decayed                                               |
|                   |              | | 5 – Number of particles that have been filtered                                         |
|                   |              | | 6 – Number of particles that left this time interval                                    |
|                   |              |                                                                                           |
|                   |              | | Note: The data found in the ".ptrk" file was previously reported in the                 |
|                   |              |   general output file. From version 2.25 of the code and forward that data                |
|                   |              |   will be reported in the optional, ".ptrk" file unless a coupled GoldSim-FEHM            |
|                   |              |   simulation is being run. In addition, the user has the option of selecting              |
|                   |              |   which statistics parameters are reported. The default is to report all                  |
|                   |              |   statistics parameters.                                                                  |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character    | Optional keyword "tcurve" indicating transfer function curve data should be               |
|                   |              | input to model matrix diffusion. If the keyword is found then NUMPARAMS and               |
|                   |              | FILENAME are entered, otherwise they are omitted.                                         |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| NUMPARAMS         | integer      | Number of parameters that define the transfer function curves being used.                 |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| FFMAX             | real         | The maximum fracture flow fraction used in the transfer function curve                    |
|                   |              | data. Default value: 0.99.                                                                |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TFILENAME         | character    | Name of input file containing the transfer function curve data.                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword ‘zptr' designating zones for breakthrough curves will be                 |
|                   |              | defined. If no keyword is input, IPZONE and IDZONE are also omitted.                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| IPZONE            | integer      | Number of zones for which breakthrough curves are to be output                            |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| IDZONE            | integer      | A list of zones for which particle breakthrough data are required. The code               |
|                   |              | outputs the number of particles that leave the system at each zone IDZONE                 |
|                   |              | at the current time step. This information is written to the ".out" file                  |
|                   |              | at each heat and mass transfer time step.                                                 |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| RSEED             | integer      | 6-digit integer random number seed.                                                       |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| RSEED_RELEASE     | integer      | | 6-digit integer random number seed for particle release location calculation.           |
|                   |              |   If a value is not entered for RSEED_RELEASE it will be set equal to RSEED.              |
|                   |              | |                                                                                         |
|                   |              | | Note that for GoldSim-FEHM coupled simulations the random seeds are controlled          |
|                   |              |   by GoldSIM and the values input in the mptr macro are not used.                         |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword ‘wtri" indicatiing a water table rise calculation should                 |
|                   |              | be performed.                                                                             |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| WATER_TABLE       | real         | | Water table elevation to be used for water table rise calculation.                      |
|                   |              | | Note that for GoldSim-FEHM coupled simulations the water table                          |
|                   |              |   rise calculations are controlled by GoldSIM and the values input                        |
|                   |              |   in the mptr macro are not used and may be omitted.                                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| DAYCS             | real         | Time which the particle tracking solution is enabled (days).                              |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| DAYCF             | real         | Time which the particle tracking solution is disabled (days).                             |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| DAYHF             | real         | Time which the flow solution is disabled (days).                                          |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| DAYHS             | real         | Time which the flow solution is enabled (days).                                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword ‘file' designating alternate transport                                   |
|                   |              | parameter file input for multiple simulation realizations.                                |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| PFILENAME         | character*80 | Name of file from which to read transport parameters.                                     |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword ‘afm' designating the Active Fracture                                    |
|                   |              | Model input for determining fracture spacing should be used.                              |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*5  | Optional keyword ‘dfree' designates that the free water                                   |
|                   |              | diffusion coefficient and tortuosity will be input instead                                |
|                   |              | of the molecular diffusion coefficient.                                                   |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TCLX              | real         | Dispersivity in the x-direction (m). The input value is                                   |
|                   |              | ignored when dispersion is turned off.                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TCLY              | real         | Dispersivity in the y-direction (m). The input value is                                   |
|                   |              | ignored when dispersion is turned off.                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TCLZ              | real         | Dispersivity in the z-direction (m). The input value is                                   |
|                   |              | ignored when dispersion is turned off.                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| APERTUR           | real         | Mean fracture aperture (m). The input value is ignored                                    |
|                   |              | when matrix diffusion is turned off.                                                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| MATRIX_POR        | real         | Porosity of the rock matrix. Used to simulate diffusion                                   |
|                   |              | and sorption in the rock matrix when matrix diffusion                                     |
|                   |              | is invoked, otherwise the input value of MATRIX_POR is ignored.                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| SRESIDUAL         | real         | Residual saturation in the Active Fracture Model used for                                 |
|                   |              | determining the spacing between active fractures.                                         |
|                   |              | This parameter is only needed when the keyword ‘afm' is                                   |
|                   |              | included, in which case the input must be entered.                                        |
|                   |              | However, the model is only used in dual permeability                                      |
|                   |              | simulations at locations where the finite spacing matrix                                  |
|                   |              | diffusion model is invoked.                                                               |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| GAMMA_AFM         | real         | Exponent in the Active Fracture Model used for determining                                |
|                   |              | the spacing between active fractures. See comments for                                    |
|                   |              | SRESIDUAL above.                                                                          |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| ITRC              | integer      | Model number for parameters defined in group 5.                                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| ITH_SPECI         | integer      | Number index of the ith species.                                                          |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TRAK_TYPE         | integer      | | Flag to denote the fluid phase of the particles:                                        |
|                   |              | | 1 - liquid phase particles                                                              |
|                   |              | | 2 - vapor phase particles                                                               |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| HALF_LIFE         | real         | Half-life for irreversible first order decay reaction(s)                                  |
|                   |              | (days). Set HALF_LIFE = 0 for no decay.                                                   |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| IDAUGHTER         | integer      | | Index of the daughter species (i.e., the index number of                                |
|                   |              |   the species to which the current species decays)                                        |
|                   |              | | If IDAUGHTER = 0, there is no decay and no ingrowth,                                    |
|                   |              | | If IDAUGHTER = -1, there is decay but no ingrowth.                                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| CONFACTOR         | real         | | Initial conversion factor for GoldSim-FEHM coupling and                                 |
|                   |              |   FEHM standalone simulations (# of particles/mole).                                      |
|                   |              | |                                                                                         |
|                   |              | | For FEHM stand alone simulations:                                                       |
|                   |              | | If CONFACTOR = 0, no conversion is necessary. The input                                 |
|                   |              |   value of PINMASS is the number of particles.                                            |
|                   |              | |                                                                                         |
|                   |              | | For GoldSim-FEHM coupling:                                                              |
|                   |              | | If CONFACTOR = 0, at each time step, the code selects a                                 |
|                   |              |   conversion factor based on the available memory and the                                 |
|                   |              |   remaining simulation time (end time - current time). The                                |
|                   |              |   code then uses the selected conversion factor to calculate                              |
|                   |              |   the number of particles to be injected at the current time step.                        |
|                   |              | |                                                                                         |
|                   |              | | For both stand alone and GoldSim-FEHM coupling cases:                                   |
|                   |              | | If CONFACTOR > 0, the code assumes the input mass is in moles                           |
|                   |              |   and uses the product of the CONFACTOR and the input mass to calculate                   |
|                   |              |   the input number of particles at each time step.                                        |
|                   |              | | When CONFACTOR >0, FEHM may use an updated conversion factor from                       |
|                   |              |   previous time step(s) as the input for the current time step instead                    |
|                   |              |   of using the original input CONFACTOR for improved results.                             |
|                   |              | |                                                                                         |
|                   |              | | If CONFACTOR < 0, the code uses the product of the absolute value                       |
|                   |              |   of CONFACTOR and the input mass (in moles) to calculate the input                       |
|                   |              |   number of particles at each time step. A CONFACTOR updated from a                       |
|                   |              |   previous time step will not be used.                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| NEWCONFACTOR      | real         | | Replace the initial value of CONFACTOR with that specified by                           |
|                   |              |   NEWCONFACTOR.                                                                           |
|                   |              | | If NEWCONFACTOR = 0, use automatic conversion factors.                                  |
|                   |              | | If NEWCONFACTOR > 0, then use the product of the CONFACTOR                              |
|                   |              |   and the input mass (in moles) to calculate the input number of                          |
|                   |              |   particles at each time step starting from CONFTIME. In this case,                       |
|                   |              |   FEHM may use an updated conversion factor from previous time step(s)                    |
|                   |              |   as a modification to CONFACTOR.                                                         |
|                   |              | | If NEWCONFACTOR < 0, then FEHM uses the product of the absolute value                   |
|                   |              |   of NEWCONFACTOR and the input mass (in moles) to calculate the input                    |
|                   |              |   number of particles at each time step (``CONFACTOR = -NEWCONFACTOR``).                  |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| CONFTIME          | real         | The time at which to change the CONFACTOR value to that specified                         |
|                   |              | by NEWCONFACTOR.                                                                          |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| GMOL              | real         | The molecular weight of the ith species. The code uses GMOL and                           |
|                   |              | CONFACTOR to convert the mass from number of particles to grams                           |
|                   |              | in the final output for GoldSim-FEHM coupling.                                            |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| P_FRACTION        | real         | The decay-ingrowth particle release factor (percentage of the                             |
|                   |              | maximum number of particles released for the current species).                            |
|                   |              | For decay-ingrowth simulations, P_FRACTION is used to reduce                              |
|                   |              | the number of particles released by parent or daughter species,                           |
|                   |              | thus, avoiding memory overflow in the daughter species due to                             |
|                   |              | parent decay-ingrowth where multiple parents decay to the same                            |
|                   |              | daughter species. The normal range of P_FRACTION is from 0 to 1.                          |
|                   |              | The default value is 0.25. A user should select an appropriate                            |
|                   |              | value based on the mass input of parent and daughter species,                             |
|                   |              | half-lives, importance of each species to the transport results,                          |
|                   |              | and simulation time period.                                                               |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| ASTEP             | integer      | Maximum length of array used to hold particle tracking information                        |
|                   |              | for the ith species. Its value depends on number of time steps,                           |
|                   |              | number of release bins, and number of parent species.                                     |
|                   |              | The sum of ASTEP for all species should be equal to or smaller                            |
|                   |              | than MAX1D.                                                                               |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| CFRACTION         | real         | The fraction of the user determined maximum number of particles                           |
|                   |              | (MAX_PARTICLES) to be assigned by mass, (1 – cfraction) will                              |
|                   |              | then be the fraction of particles assigned by time step.                                  |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword ‘size' designating that the colloid size                                 |
|                   |              | distribution model option is being used (combined with                                    |
|                   |              | the interface filtration option in the itfc macro). If the                                |
|                   |              | keyword is not input, PART_SIZE and PROBSIZE are also omitted.                            |
|                   |              | Colloid size is only sampled once for each realization.                                   |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| PART_SIZE         | real         | Colloid particle size for this entry of the particle size                                 |
|                   |              | distribution table (paired with a value of PROBSIZE).                                     |
|                   |              | An arbitrary number of entries can be input, terminated with                              |
|                   |              | a blank line. The code assigns each particle a size based on                              |
|                   |              | this distribution of particle sizes, and decides if particles                             |
|                   |              | are irreversibly filtered based on the pore size distribution                             |
|                   |              | assigned in the itfc macro.                                                               |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| PROBSIZE          | real         | Colloid cumulative probability for the distribution of sizes                              |
|                   |              | (paired with a value of PART_SIZE). See description of                                    |
|                   |              | PART_SIZE above for details. The final entry of the table                                 |
|                   |              | must have PROBSIZE = 1, since the distribution is assumed                                 |
|                   |              | to be normalized to unity.                                                                |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword "dive" signifying that the specie being                                  |
|                   |              | specified is either a colloid species using the colloid                                   |
|                   |              | diversity model or a non-colloid daughter species of a                                    |
|                   |              | colloid species.                                                                          |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| FLAG_COL_DAUGHTER | integer      | When FLAG_COL_DAUGHTER = 1 signals that the species being                                 |
|                   |              | specified is a non-colloid species that can result as a                                   |
|                   |              | daughter product of a colloid parent species. If the species                              |
|                   |              | is not a daughter product or the daughter product is a                                    |
|                   |              | colloid, FLAG_COL_DAUGHTER = 0.                                                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character*4  | Optional keyword ‘file' designating the cumulative probability                            |
|                   |              | distribution function (CDF) retardation parameters for the                                |
|                   |              | colloid diversity model should be read from an external file.                             |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| CDF_FILENAME      | character*80 | | Name of the file containing the cumulative probability                                  |
|                   |              |   distribution function (CDF) (entered if optional keyword                                |
|                   |              |   ‘file' follows keyword ‘dive'). See below for file formats.                             |
|                   |              | | If TPRPFLAG = 11 or 12, Table option                                                    |
|                   |              | | If TPRPFLAG = 13 or 14, Equation option                                                 |
|                   |              | |                                                                                         |
|                   |              | | The following equations are used for :math:`R_{min} \le R \le R_{max}`,                 |
|                   |              | | :math:`R = 1 + K_f / K_{rev}`,                                                          |
|                   |              |   :math:`\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)`                                     |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TPRP_FLAG         | integer      | | Values of TPRPFLAG between 11 and 14 signify that the colloid                           |
|                   |              |   diversity model with equal weight sampling will be used:                                |
|                   |              | | TPRPFLAG = 11: CDF vs retardation factor specified in a table                           |
|                   |              | | TPRPFLAG = 12: similar to 11, but the SQRT(CDF) is used instead                         |
|                   |              |   of CDF for sampling                                                                     |
|                   |              | | TPRPFLAG = 13: CDF vs :math:`K_f` (Attachment rate constant)                            |
|                   |              |   specified as a straight line equation in the log-log space                              |
|                   |              | | TPRPFLAG = 14: similar to 13, but the SQRT(CDF) is used                                 |
|                   |              |   instead of CDF for sampling                                                             |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| SIMNUM            | integer      | Simulation number, used for selecting the table/equation from                             |
|                   |              | the colloid diversity file. For GoldSim-FEHM coupled simulations                          |
|                   |              | or FEHM runs using the ‘msim' option this parameter is passed                             |
|                   |              | to the code. For non-coupled simulations it is an optional                                |
|                   |              | input. (Default value = 1)                                                                |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| K_REV             | real         | Detachment rate constant for reversible filtration of                                     |
|                   |              | irreversible colloids.                                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| R_MIN             | real         | Minimum value of the retardation factor for reversible                                    |
|                   |              | filtration of irreversible colloids.                                                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| R_MAX             | real         | Maximum value of the retardation factor for reversible                                    |
|                   |              | filtration of irreversible colloids                                                       |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| SLOPE_KF          | real         | | Value of the slope (:math:`m`) in the log-log space                                     |
|                   |              |   for the equation:                                                                       |
|                   |              | | :math:`\log_{10}(CDF) = b + m \cdot \log_{10}(K_f)`                                     |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| CINT_KF           | real         | Value of the intercept (:math:`b`) in the log-log space                                   |
|                   |              | for the above equation                                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KEYWORD           | character    | Keyword specifying whether the colloid species is                                         |
|                   |              | ‘irreversible' or ‘reversible'.                                                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| FLAG_LOG          | integer      | | For reversible colloids an average retardation factor is used:                          |
|                   |              | |                                                                                         |
|                   |              | | If FLAG_LOG = 0: a linear average of the distribution is used                           |
|                   |              | | If FLAG_LOG = 1: a log-linear average of the distribution                               |
|                   |              |   is used                                                                                 |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| LAYERS            | integer      | Number of layers in which the transport properties of the                                 |
|                   |              | ith species are to be modified. If no property is altered,                                |
|                   |              | then set layers=0.                                                                        |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| LAYER_I           | integer      | The index number of the ith layer defined in group 5                                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TRANSFLAG         | integer      | | Flag to specify which transport mechanisms apply [abs(TRANSFLAG)]:                      |
|                   |              | | 1 - advection only (no dispersion or matrix diffusion)                                  |
|                   |              | | 2 - advection and dispersion (no matrix diffusion)                                      |
|                   |              | | 3 - advection and matrix diffusion, infinite fracture spacing                           |
|                   |              |   solution (no dispersion)                                                                |
|                   |              | | 4 - advection, dispersion, and matrix diffusion, infinite fracture                      |
|                   |              |   spacing solution                                                                        |
|                   |              | | 5 - advection and matrix diffusion, finite fracture spacing                             |
|                   |              |   solution (no dispersion)                                                                |
|                   |              | | 6 - advection, dispersion, and matrix diffusion, finite fracture                        |
|                   |              |   spacing solution                                                                        |
|                   |              | | 8 - use the the transfer function approach with 3 dimensionless                         |
|                   |              |   parameters and type curves for handling fracture-matrix interactions.                   |
|                   |              | |                                                                                         |
|                   |              | | For TRANSFLAG < 0, transport simulations include colloids.                              |
|                   |              | | For equivalent continuum solutions, the fracture spacing in the                         |
|                   |              |   finite spacing model is determined using                                                |
|                   |              | |                                                                                         |
|                   |              | | :math:`SPACING = APERTURE / POROSITY`                                                   |
|                   |              | |                                                                                         |
|                   |              | | For dual permeability models, the fracture spacing input parameter                      |
|                   |              |   APUV1 in the ``dpdp`` macro is used as the half-spacing between fractures.              |
|                   |              |   If the Active Fracture Model (see keyword ‘afm') is used, APUV1 is the                  |
|                   |              |   geometric fracture half-spacing, and the additional terms SRESIDUAL                     |
|                   |              |   and GAMMA_AFM are used to determine the spacing between active                          |
|                   |              |   fractures (see below).                                                                  |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KD                | real         | Sorption coefficient (linear, reversible, equilibrium sorption).                          |
|                   |              | Units are kg-fluid / kg-rock (these units are equivalent to the                           |
|                   |              | conventional units of cc/g when the carrier fluid is water at                             |
|                   |              | standard conditions). This value applies to the medium as a whole                         |
|                   |              | when matrix diffusion is turned off, whereas for simulations invoking                     |
|                   |              | matrix diffusion, the value applies to the rock matrix.                                   |
|                   |              | For the latter case, sorption in the flowing system (fractures)                           |
|                   |              | is modeled using the RD_FRAC variable.                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| RD_FRAC           | real         | Retardation factor within the primary porosity (fractures) for a                          |
|                   |              | matrix diffusion particle tracking simulation (use 1 for no                               |
|                   |              | sorption on fracture faces). The input value is ignored unless                            |
|                   |              | matrix diffusion is invoked.                                                              |
|                   |              |                                                                                           |
|                   |              |                                                                                           |
|                   |              |                                                                                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| DIFFMFL           | real         | Molecular diffusion coefficient in the rock matrix (m2/s).                                |
|                   |              | The input value is ignored unless matrix diffusion is invoked.                            |
|                   |              |                                                                                           |
|                   |              |                                                                                           |
|                   |              |                                                                                           |
|                   |              |                                                                                           |
|                   |              |                                                                                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| H2O_DIFF          | real         | | Free water diffusion coefficient. The molecular diffusion                               |
|                   |              |   coefficient is calculated as                                                            |
|                   |              |   :math:`H2O\_DIFF \times TORT\_DIFF`                                                     |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TORT_DIFF         | real         | Tortuosity                                                                                |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| KCOLL             | real         | | Colloid distribution parameter, the ratio of contaminant mass                           |
|                   |              |   residing on colloids to the mass present in aqueous form.                               |
|                   |              |   It is used to compute an effective aperture via the following:                          |
|                   |              | | :math:`APWID = APERERTURE \cdot (1 + KCOLL)`                                            |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| RCOLL             | real         | | Colloid retardation factor. Used, in conjunction with kcoll,                            |
|                   |              |   to adjust colloid retardation in fractures using the following                          |
|                   |              |   formula:                                                                                |
|                   |              | | :math:`FRACRD = \frac{RD\_FRAC + KCOLL \cdot RCOLL}{1+KCOLL}`                           |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| FCOLL             | real         | Colloid filtration parameter. Used to compute the probability a colloid                   |
|                   |              | will be irreversibly filtered along the path between two nodes using                      |
|                   |              | the following:                                                                            |
|                   |              | :math:`PROBFILT = 1 - \exp(DISTANCE/FCOLL)` where                                         |
|                   |              | :math:`DISTANCE` iis the leength of the path between nodes.                               |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| NS                |              | Number of spatial source terms for the ith species                                        |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| TMPCNSK           | real         | | Particle injection parameter assigned for nodes defined by JA,                          |
|                   |              |   JB, and JC. Two options are available:                                                  |
|                   |              | |                                                                                         |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
|                   |              | | TMPCNSK > 0. - particles are injected at each node in                                   |
|                   |              |   proportion to the source mass flow rate at the node. This                               |
|                   |              |   boundary condition is equivalent to injecting a solute of                               |
|                   |              |   a given concentration into the system. Note: the source                                 |
|                   |              |   flow rates used to assign the number and timing of particle                             |
|                   |              |   injections are those at the beginning of the particle                                   |
|                   |              |   tracking simulation (time DAYCS). Transient changes in this                             |
|                   |              |   source flow rate during the particle tracking simulation do                             |
|                   |              |   not change the number of particles input to the system.                                 |
|                   |              | | TMPCNSK < 0. - particles are introduced at the node(s),                                 |
|                   |              |   regardless of whether there is a fluid source at the node.                              |
|                   |              | |                                                                                         |
|                   |              | | Default is 0. for all unassigned nodes, meaning that no                                 |
|                   |              |   particles are injected at that node.                                                    |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| PINMASS           | real         | Input mass. If CONFACTOR = 0, PINMASS is the number of particles                          |
|                   |              | to be injected at locations defined by TMPCNSK.                                           |
|                   |              | If CONFACTOR > 0, PINMASS is the input mass expressed in moles.                           |
|                   |              | The code uses CONFACTOR to convert PINMASS into number of particles.                      |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| T1SK              | real         | Time (days) when particle injection begins. Default is 0.                                 |
+-------------------+--------------+-------------------------------------------------------------------------------------------+
| T2SK              | real         | Time (days) when particle injection ends. Default is 0.                                   |
+-------------------+--------------+-------------------------------------------------------------------------------------------+


.. note::
   **Notes on Restarting:** As with all restart runs for FEHM, a ".ini" file is
   specified to be read to set the initial conditions upon restarting. However,
   there are two possibilities for restart calculations with particle tracking
   (mptr or ptrk): 1) the heat and mass transfer solution is being restarted, but
   the particle tracking simulation is initiated during the restart run (it was not
   carried out in the simulation that generated the ".ini" file); or 2) the heat
   and mass transfer solution and the particle tracking simulation are both being
   restarted. If the code does not find the "ptrk" key word at the top of the ".ini"
   file, then the original run did not employ particle tracking, and Case 1 is assumed.
   A common example is a preliminary calculation that establishes a fluid flow steady
   state, followed by a restart simulation of transport.

If "ptrk" was written into the ".ini" file in the original run, the particle data in the ".ini" file are read and used to initialize the particle tracking simulation (Case 2). In this instance, the number of particles (NPART) must be set the same for the restart run as in the original run or the results will be unpredictable.When restarting a particle tracking simulation, certain input data are overwritten by information in the ".ini" file. These parameters include RSEED, RSEED_RELEASE, PCNSK, T1SK, and T2SK. Other input parameters can be set to different values in the restart run than they were in the original run, but of course care must be taken to avoid physically unrealistic assumptions, such as an abrupt change in transport properties (Group 4 input) part way through a simulation.

A final note on restart calculations is in order. A common technique in FEHM restart calculations is to reset the starting time at the top of the ".ini" file to 0 or in the time macro so that the starting time of the restart simulation is arbitrarily 0, rather than the ending time of the original simulation. This is useful for the example of the steady state flow calculation, followed by a restart solute transport calculation. Although this technique is acceptable for particle tracking runs that are initiated only upon restart (Case 1), it is invalid when a particle tracking run is being resumed (Case 2). The reason is that all particle times read from the ".ini" file are based on the starting time of the original simulation during which the particle tracking simulation was initiated.

The following is an example of mptr. A multiple-species decay-chain
(:math:`\rightarrow 2 \rightarrow 3`) is simulated, with decay half lives of the
species equaling 10,000, 3,000, 10,000, and 4,000 years, respectively. In this
simulation a maximum of 3 property layers are specified although only 1 layer is used,
the maximum number of particles is specified to be 1100100, and FEHM is run in stand-alone
mode. Concentrations will be computed as number of particles per unit fluid volume and
no output will be written to the “.fin” file. Use of the ‘zptr’ keyword indicates that
a single zone will be defined for breakthrough curve output which will be written to
the “.out” file. The random number seed is defined to be 244562. The particle tracking
solution is enabled at 0.1 days, and disabled at 3.65e8 days, while the flow solution
is disabled at 38 days and re-enabled at 3.65e8 days. Dispersivity in the X-, Y-, and
Z- directions are defined to be 0.005 m, the mean fracture aperture is 0.0001 m, and
the matrix porosity is 0.3. Particles for species 1 are injected at a constant rate
from 0 to 5,000 years, and species 2, 3, and 4 are formed through the decay reactions,
with no input at the inlet. Advection and dispersion (without matrix diffusion) is being
modeled. The retardation factors for the four species are 1, 1, 1.9, and 1, respectively
(i.e. only species 3 sorbs).

+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| mptr   |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 4      | 3          | 1100100    | 0         |        |    |    |    |     | Group 1  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 2      | 0          |            |           |        |    |    |    |     | Group 2  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| zptr   |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 244562 |            |            |           |        |    |    |    |     | Group 3  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 0.1    | 3.65e8     | 38         | 3.65e8    |        |    |    |    |     | Group 4  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 0.005  | 0.005      | 0.005      | 1.e-4     | 0.3    |    |    |    |     | Group 5  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
|        |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 0          | 0          | 1         |        |    |    |    |     | Group 6  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
|        |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 1          | 3.652485E6 | 2         | 1      | -1 | 1. | 1. | 0.5 | Group 7  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 8  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 2          | 0.         | 1.        | 1.e-14 |    |    |    |     | Group 9  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 10 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 202        | 201        | -1        |        |    |    |    |     | Group 11 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 10000. | 0.         |            | 365.25E2  |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 10000. | 365.25E2   |            | 730.5E2   |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 10000. | 730.5E2    |            | 1.09575E5 |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 10000. | 1.09575E5  |            | 1.461E5   |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| .      |            |            |           |        |    |    |    |     | .        |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| .      |            |            |           |        |    |    |    |     | .        |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| .      |            |            |           |        |    |    |    |     | .        |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 10000. | 1.7532E6   | 1.789725E6 |           |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 10000. | 1.789725E6 | 1.82625E6  |           |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
|        |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 2      | 1          | 1.095745E6 | 3         | 0      | -1 | 1. | 1. | 0.5 | Group 7  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 8  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 2          | 0.         | 1.        | 1.e-14 |    |    |    |     | Group 9  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 10 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 202        | 201        | -1        |        |    |    |    |     | Group 11 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 0      | 0.         | 1.825E6    |           |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
|        |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 3      | 1          | 3.652485E6 | 4         | 0      | -1 | 1. | 1. | 0.5 | Group 7  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 8  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 2          | 0.108      | 1.        | 1.e-14 |    |    |    |     | Group 9  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 10 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 202        | 201        | -1        |        |    |    |    |     | Group 11 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 0      | 0.         | 1.825E6    |           |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
|        |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 4      | 1          | 1.460972E6 | -1        | 0      | -1 | 1. | 1. | 0.5 | Group 7  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 8  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 2          | 0.         | 1.        | 1.e-14 |    |    |    |     | Group 9  |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      |            |            |           |        |    |    |    |     | Group 10 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 1      | 202        | 201        | -1        |        |    |    |    |     | Group 11 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
| 0      | 0.         | 1.825E6    |           |        |    |    |    |     | Group 12 |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+
|        |            |            |           |        |    |    |    |     |          |
+--------+------------+------------+-----------+--------+----+----+----+-----+----------+

In the second example, transfer function data is used along with the active fracture and colloid diversity models. The the cumulative probability distribution function (CDF) retardation parameters for the colloid diversity model are entered in an external file using the table format. The format for the new input files associated with the colloid diversity model are:

For ptrk/sptr/mptr simulations with TPRP_FLAG = 11 or 12:

* Header line indicating the species number (always 1 for ptrk/sptr simulations)

* Multiple tables, each with the following format:

  * One line specifying the realization number 

  * Multiple lines with two columns of data (real) representing rcdiv and probdiv. Note that probdiv should start at 0 and end with 1. 

  * A blank line is used to specify the end of the distribution table.

* A blank line is used to specify the end of the file

For mptr, a distribution will need to be entered for each colloid species. Therefore,
the header line and tables are repeated for each colloid species.

For ptrk/mptr simulations with TPRP_FLAG = 13 or 14:

* Header line containing comments 

* Multiple lines, each line containing realization number, :math:`b, m, k_r, R_{min}, R_{max}`

* A blank line is used to specify the end of the file

For sptr simulations with TPRP_FLAG = 13 or 14:

* Header line containing comments 

* Multiple lines, each line containing realization number, :math:`b, m, k_r, R_{min}, R_{max}, \alpha_L`

* A blank line is used to specify the end of the file

Particle statistics data for the cumulative number of particles that have left the sytem and the number of particles that left during the current timestep are output.

+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| mptr                                         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 3                                            | 100    | 500000 | 0        |          |       |              |                         |          | Group 1  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 0                                            | 30     | 3      | 6        |          |       |              |                         |          | Group 2  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| tcurve                                       |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 3                                            |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| ../colloid_cell/input/uz_tfcurves_nn_3960.in |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 244562                                       |        |        |          |          |       |              |                         |          | Group 3  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 0.                                           | 1.e20  | 0.     | 1.e20    |          |       |              |                         |          | Group 4  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| afm                                          |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1.e-03                                       | 1.e-03 | 1.e-03 | 1.e-03   | 0.2      | 0.01  | 0.6          | // layer 1 tcwm1 zone 1 |          | Group 5  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1.e-03                                       | 1.e-03 | 1.e-03 | 0.00e+00 | 0.2      | 0.00  | 0.0          | // layer 1 tcwm1 zone 2 |          | Group 5  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
|                                              |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            | 10     | 1      | 1        |          |       |              |                         |          | Group 6  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 11                                           | 20     | 1      | 2        |          |       |              |                         |          | Group 6  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
|                                              |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            | 1      | 0      | -1       | 1        | 0     | 1.00E+15243  | 1.0                     | speci1   | Group 7  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| diveresity                                   |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 0                                            |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| file                                         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| ../colloid_cell/input/rcoll_data.dat         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 11 1                                         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| irreversible                                 |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 2                                            |        |        |          |          |       |              |                         |          | Group 8  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            | -2     | 0.0    | 1.0e+00  | 1.00e-30 | 1e+20 | 1            | 1.00                    | #1 tcwM1 | Group 9  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 2                                            | -2     | 0.0    | 1.0e+00  | 1.00e-30 | 1e+20 | 1            | 1.00                    | #2 tcwM2 | Group 9  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
|                                              |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            |        |        |          |          |       |              |                         |          | Group 10 |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            | 1      | 1      | -1       |          |       |              |                         |          | Group 11 |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 100000                                       | 0.     | 0.01   |          |          |       |              |                         |          | Group 12 |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
|                                              |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 3                                            | 1      | 0      | -1       | 1        | 0     | 1.00EE+15243 | 1.0                     | speci3   | Group 7  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| diversity                                    |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 0                                            |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| file                                         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| ../colloid_cell/input/rcoll_data.dat         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 11 3                                         |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| irreversible                                 |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 2                                            |        |        |          |          |       |              |                         |          | Group 8  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            | -2     | 0.0    | 1.0e+00  | 1.00e-30 | 1e+20 | 1            | 1.00                    | #1 tcwM1 | Group 9  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 2                                            | -2     | 0.0    | 1.0e+00  | 1.00e-30 | 1e+20 | 1            | 1.00                    | #2 tcwM2 | Group 9  |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
|                                              |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            |        |        |          |          |       |              |                         |          | Group 10 |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 1                                            | 1      | 1      | -1       |          |       |              |                         |          | Group 11 |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
| 0                                            | 0.     | 0.01   |          |          |       |              |                         |          | Group 12 |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+
|                                              |        |        |          |          |       |              |                         |          |          |
+----------------------------------------------+--------+--------+----------+----------+-------+--------------+-------------------------+----------+----------+

With the file ``rcoll_data.dat`` as:

+---------+----------+--------+--------+
| Test file for 1D importance sampling |
+---------+----------+--------+--------+
| 1       | 0.0      | 0.0    | 0.0    |
+---------+----------+--------+--------+
| 1.0     | 0.0      |        |        |
+---------+----------+--------+--------+
| 2.0     | 0.125    |        |        |
+---------+----------+--------+--------+
| 3.0     | 0.25     |        |        |
+---------+----------+--------+--------+
| 4.0     | 0.375    |        |        |
+---------+----------+--------+--------+
| 5.0     | 0.5      |        |        |
+---------+----------+--------+--------+
| 6.0     | 0.625    |        |        |
+---------+----------+--------+--------+
| 7.0     | 0.75     |        |        |
+---------+----------+--------+--------+
| 8.0     | 0.875    |        |        |
+---------+----------+--------+--------+
| 9.0     | 1.0      |        |        |
+---------+----------+--------+--------+
|         |          |        |        |
+---------+----------+--------+--------+
| .       |          |        |        |
+---------+----------+--------+--------+
| .       |          |        |        |
+---------+----------+--------+--------+
| .       |          |        |        |
+---------+----------+--------+--------+
|         |          |        |        |
+---------+----------+--------+--------+
| 3       | 0.0      | 0.0    | 0.0    |
+---------+----------+--------+--------+
| 1.0     | 0.0      |        |        |
+---------+----------+--------+--------+
| 2.0     | 0.125    |        |        |
+---------+----------+--------+--------+
| .       |          |        |        |
+---------+----------+--------+--------+
| .       |          |        |        |
+---------+----------+--------+--------+
| .       |          |        |        |
+---------+----------+--------+--------+
| 9.0     | 1.0      |        |        |
+---------+----------+--------+--------+
|         |          |        |        |
+---------+----------+--------+--------+



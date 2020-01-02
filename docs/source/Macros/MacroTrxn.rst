========
``trxn``
========

Overview
--------

The ``trxn`` macro is designed to replace the ``trac`` (and optionally ``rxn``) macros with a more convenient, user-friendly input format.  Rather than reading input as sets of nameless numbers, the ``trxn`` macro reads several "blocks" of data in which parameters are defined.  The blocks may be specified in any order (with the exception of the ``ctrl`` and ``lookup`` blocks), and all blocks are optional, except for the ``header`` and ``comp`` blocks.  Any parameters that are not specified will be given default values (usually zero).

The ``trxn`` macro relies heavily on zones, and uses zones for applying all variables that can vary by node.  For this reason, a ``zone`` macro must be supplied in the input file before `trxn` is read.  A ``time`` macro must also be given before ``trxn``.

If the pound character ("#") appears on a line, everything after it on that line will be treated as a comment.  Lines in which the first character is a pound sign are treated as blank lines.  Blocks must not contain any blank or commented lines (although they can contain comments that do not start at the beginning of the line), and blocks must be separated by at least one blank or commented line.  In the ``comp``, ``water``, ``rock``, ``gas``, ``disp``, ``sorp``, and ``assign`` blocks, entire columns can be commented out.  If an asterisk ("*") (separated from surrounding tokens by whitespace) appears in the column header line of one of these blocks, the contents of every column to the right of the asterisk will be ignored.  An entire block (until the next blank line) can be skipped by placing the keyword ``null`` before the block name.

trxn blocks
-----------

The blocks are as follows:

``ctrl``
--------

This block contains control parameters, all on the same line.  If it is supplied, it must be the first block in trxn.  Its options are as follows:

+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Option         | Description                                                                                                                                                         |
+================+=====================================================================================================================================================================+
| ``rxnon``      | Enables reactions, which are by default disabled.  Reactions will not occur in a simulation if `rxnon` is not given, even if reaction-related blocks are specified. |
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``co2_couple`` | Enables CO2 coupling for `rxn`.                                                                                                                                     |
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``debug``      | Enables the output of debugging information.                                                                                                                        |
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Below is an example of the `ctrl` block:

.. code::

   ctrl co2_couple rxnon


In this example, CO2 coupling and reactions are enabled.

``include``
-----------

The `include` block is designed to facilitate the construction of "libraries" of rock/water/gas types, component properties, dispersion models, etc. that can be shared between simulations.  It allows blocks to be included from external files as though these files' contents were placed at the location of the `include` statement.  External include libraries have the exact same syntax as the standard `trxn` macro.  The first non-blank (and non-commented) line of the include file should be `trxn library`.  Following this, as many blocks as desired may be supplied. 
The `end trxn` line should be given at the end of the library file;
otherwise, reading will terminate early. 
Multiple `include` statements referencing several different files may be used in
the same `trxn` macro.  Libraries may include other libraries.

The syntax is as follows:

.. code::

   include /scratch/nts/ms/trxn/test-problems/standard.trxl


This will include the library located at ``/scratch/nts/ms/trxn/test-problems/standard.trxl`` at the current location in the input file.  As an example, a simple library might look something like this:

.. code::

   # Standard library for use with all test problems
   
   trxn library
   
   header
   0	1	1e-6	1
   0	100000	100000	100000
   5	2	1	1	0
   iskip=0, rsdmax=1e-9
   
   disp	lx	ly	lz	vx	vy	vz
   std_x	0.1	1e-9	1e-9	0.1	1e-9	1e-9
   std_y	1e-9	0.1	1e-9	1e-9	0.1	1e-9
   std_z	1e-9	1e-9	0.1	1e-9	1e-9	0.1
   
   sorp	ltype	a1l	a2l	bl	vtype	a1v	a2v	bv
   .std
   *	con	0	0	1	con	0	0	1
   
   diff	l=1e-9, v=1e-9
   
   end trxn


All problems that include this library will be given the standard values for `header` and liquid and vapor diffusion coefficients of 1×10^-9^, and will have one sorption model ("std") and three dispersion models ("std_x", "std_y", and "std_z") available.

Caution should be used when including libraries to prevent multiple definitions of any blocks, as this may cause unexpected behavior.  `trxn` will print a warning message to the error file if a block is supplied more than once, and these warnings should be considered.  Also, care should be taken to avoid placing any problem-specific parameters in the libraries.

``header``
----------

This block contains basic constant values, and is the first three lines of the `trac` macro copied verbatim.

+------+----------------------------------------+
| Line | Variables                              |
+======+========================================+
| 1    | ``AN0 AWC EPC UPWGTA``                 |
+------+----------------------------------------+
| 2    | ``DAYCS DAYCF DAYHF DAYHS``            |
+------+----------------------------------------+
| 3    | ``IACCMX DAYCM DAYCMM DAYCMX NPRTTRC`` |
+------+----------------------------------------+
| 4    | Optional parameters                    |
+------+----------------------------------------+

Please refer to the `trac` section of the FEHM User's Manual for the meanings of the parameters on the first three lines.  On the fourth line, optional parameters may be defined, using the form `variable=value` (no spaces), with commas between these pairs.  The variables that may be set here are `ISKIP`, `RSDMAX`, and `STRAC_MAX`.  Please refer to the `rxn` section of the FEHM User's Manual for the meanings of these variables.  (`STRAC_MAX` is the maximum allowable saturation for vapor and Henry's Law species.  When using `trac`, it is read if provided in the input file; however, it is not documented in the `trac` portion of the FEHM User's Manual.)  `ISKIP` and `RSDMAX` are used only if reactions are enabled (see `rxn` below).  If these optional variables are omitted, they are given the values shown in the example below.

Below is an example of the `header` block.

.. code::

   header
   0.0	1.0	1.0e-6	1.0			# ANO, AWC, EPC, UPWGTA
   0.0	1e20	1e20	1e20			# DAYCS, DAYCF, DAYHF, DAYHS
   10	2.0	1.0	150000.0	1	# IACCMX, DAYCM, DAYCMM, DAYCMX, NPRTTRC
   iskip=0, rsdmax=1e-9, strac_max=0.99


Please note that the keyword ``userc`` is not supported in the ``header`` block.  For ``userc`` support, please refer to the `userc` block below.  Also, unlike in ``trac``, ``NPRTTRC`` may not be omitted from the header.

``userc``
---------

This block invokes the solute transport user subroutine as specified in the ``trac`` section of the user's manual.  On the same line as the macro name should be given the path to a file containing userc parameters.  See the `trac` section of the user's manual for more information on the `userc` input format.

Below is an example of ``userc``:

.. code::

   userc input/userc.dat


In this example, the ``userc`` subroutine is called with an input file located at ``input/userc.dat``.

``comp``
--------

The ``comp`` block is used to define each component present in the simulation.  (A component is any group of compounds, ions, etc., all of the same phase, that the user wishes to be treated as a single entity by the tracer solver.)  It contains one line for each component, and each line consists of a phase designation for that component and the name of the component.  The phase designation is one of "aqueous", "solid", "gas", and "henry"; these may be shortened to the first character to save time.

If kinetic reactions (``rxn`` macro) are being simulated, two additional columns may be included.  The `master` column indicates the "master species" for each aqueous or Henry's Law component.  These master species are arbitrarily chosen forms of the components, by convention the form that is expected to dominate in reactions.  The third column, `guess`, allows the user to specify a guess for the initial uncomplexed concentration of each aqueous and Henry's Law component.  This is not necessary unless the chemical speciation solver has difficulty converging with the default value of 1×10^-9^, in which case specifying a more representative value may help.  If this is not necessary, leave the column out entirely, or place asterisks in the rows of components that do not need help converging.

Below is an example of the ``comp`` block:

.. code::

   comp		master	guess
   aqueous H	H+	*
   aq C_a		HCO3-	8.6e-6
   a Na_a		Na+	*
   a Ca_a		Ca++	*
   a C_a2		CO3--	*
   a Cl_a		Cl-	2.0e-9
   a U238		U02	*
   a Th234		ThO2	*
   solid AlO3	*	*
   s NaCl		*	*
   s NaHCO3	*	*
   s CaCl2		*	*
   gas O2		*	*
   g N2		*	*
   henry C_h2	C6H6	1.2e-8
   h Cl_h		Cl2	*
   h C_h		CO2	*


In this example, there are 17 components.  "H", "C_a", "Na_a", "Ca_a", "C_a2",
"Cl_a", "U238", and "Th234" are aqueous; "AlO3", "NaCl", "NaHCO3", and "CaCl2"
are solid; "O2" and "N2" are gaseous, and "C_h", "Cl_h", and "C_h2" may be liquid
or gas according to Henry's Law.

``water``, ``rock``, ``gas``
----------------------------

These blocks are identical in form, and they are used to assign concentrations in
 the simulation for components of different states.  These blocks allow the user to
 specify different "water types", "rock types", and "gas types", which may consist
 of different combinations of components specified in `comp` in different concentrations.
 See the `moles` block below for an alternative input format.


On the same line as the block name, the names of components specified in `comp`
 are placed, separated by tabs or spaces.  These are column headers.  Below this line,
 one line is given to each "type" desired.  Each line consists of the name of the type,
 followed by numbers representing the concentrations of each of the components given
 in the column headers in that type.  If the columns are separated by tabs,
 this layout forms a neat table with types down the left side, components across
 the top, and the concentration of each component in each type at the intersections
 of the rows and columns.


Only aqueous and Henry's Law components may be included in the `water` block,
only solid components in the `rock` block, and only gaseous and Henry's Law
components in the `gas` block.  In the `water` block, a special column header,
`pH`, may be included.  This is the same as heading the column with "H"
(and may only be done if "H" is specified in `comp` and is aqueous), but
allows the user to enter H+ concentration in terms of pH rather than molarity.
If a concentration in the grid is negative, it is assumed that the value entered
is the base-ten logarithm of the actual concentration, and the value is adjusted
accordingly.


An asterisk ("*") in a space where a number is expected is the same as a 0. If a
component is omitted entirely from the table, it is assumed that that
component is not present in the simulation.
The unit for all numbers in these tables is molal, with the exception of the
`pH` column, if present.


A negative value in the grid for a water or gas inflow type indicates that
the concentration of that solute will be held constant at inflow nodes of that water type.
Negative values in initial condition types should be avoided.


An example of each of the four block types, consistent with the sample `comp`
block above, is given below:

.. code::

   water	pH	C_a	Na_a	C_a2	C_h2
   wt1	7.43	1e-6	1.34e-5	1e-10	0.3
   wt2	5.4	1e-3	0.002	2.4e-6	2.3
   wt3	2.3	0	10	0	0
   
   rock	AlO3	NaCl	NaHCO3
   tuff	3.4e-2	1.24	1.6e-6
   granite	3e-5	1.5e-2	6e-5
   clay	0.3	3.9e-3	0.03
   
   gas	O2	N2	Cl_h	C_h
   vt1	0.23	0.10	1.24	2
   vt2	1.02	0.012	0.2	0


In this example, there are three water types ("wt1", "wt2", and "wt3"), three rock
types ("tuff", "granite", and "clay"), and two vapor types ("vt1" and "vt2").
Water type "wt1" has a pH of 7.3, an HCO3- concentration of 1×10^-6^ //m//, an Na+
concentration of 1.34×10^-5^ //m//, a CO3-- concentration of 1×10^-10^ //m//,
and a C6H6 concentration of 0.3 //m//.

``print``
---------

The `print` block allows the user to specify for which components and complexes information is to be printed to tracer output(`.trc`) files.  This block occupies only one line.  After the keyword `print` is given a list of aqueous component and complex names that are to be printed, delimited by spaces, tabs, or commas.  The keywords `all` and `none` may be given instead of the list, specifying that information is to be printed for all or none of the aqueous components and complexes, respectively.  The default action, if no `print` block is specified, is the same as specifying `print all`.  This information is only used if reactions are enabled.

Below is an example of the `print` block:

.. code::

   print UO2 HCO3- Na_a Ca_a


In this example, information will be printed only about UO2, HCO3-, Na_a, and Ca_a.

``moles``
---------

This block allows the user to specify initial solute conditions as the total number of moles contained in each zone.  FEHM will distribute the specified number of moles evenly throughout the volume of the zone.  The `moles` block conflicts with the `water`, `rock`, and `gas` columns in the `assign` block, and they cannot both appear in the same `trxn` macro.

The `moles` block allows the use of custom zone definitions.  Unlike standard zones, these zones are permitted to overlap, and the concentrations at nodes in overlapping zones become cumulative.  If these custom zones are to be used, a separate file must be created, containing a valid zone macro that defines the desired zones.  Due to limited functionality in the `trac` macro, all zones in this macro must be specified using the `nnum` method and zones must be numbered rather than named.  (See the manual section on the `zone` macro for more information.)  Alternatively, the zones defined in previous `zone` macros may be used.

On the same line as the `moles` keyword in this block, the path to the alternate zone file should be provided.  If the standard zone definitions are to be used, leave the path blank or explicitly turn the alternate zone processing off by providing an asterisk in the place of the path.  The next lines form a table with component names from `comp` across the top and zone numbers down the left-hand side.  At the intersections of rows and columns is placed the number of moles of that component initially in that zone.  Any omitted zones or components are assigned a zero concentration.

Below is an example of the `moles` block:

.. code::

   moles input/moles.zone
   	C_a	Na_a	C_a2	C_h2
   2	0	1.6	0.02	1
   3	0	1.7	0.6	0.89
   5	0	2	0.045	2.2
   4	1.2	0	0.102	1.61


In this example, zone definitions are loaded from the file `input/moles.zone`, which remain in effect for the rest of the `moles` macro.  Zone 2 contains 1.6 moles of Na_a, 0.02 moles of C_a2, and 1 mole of C_h2 distributed evenly throughout it.

``hparam``
----------

The ``hparam`` block is used to assign Henry's Law parameters to different Henry's Law components.  Each Henry's Law component appearing in `comp` must be included in this block.  Below the block title `hparam`, each Henry's Law component is given one line.  Each of these lines consists of the name of the Henry's Law component, the model that will be used to simulate that component, and then several `parameter=value` pairs (separated by commas), specifying parameters for each model.  The models and their parameters are as follows:

+----------------+-------------------------------------------------------------------------------------------------------------------------+
| Option         | Description                                                                                                             |
+================+=========================================================================================================================+
| `1` or `hoff`  | van't Hoff model.  Parameters:  :math:`ah - A_H; dhh - \delta H_H`                                                      |
+----------------+-------------------------------------------------------------------------------------------------------------------------+
| `2` or `multi` | Multi-parameter fit.  Parameters: :math:`ahn - A_Hn` for each n from 1 to 5                                             |
+----------------+-------------------------------------------------------------------------------------------------------------------------+
| `3` or `wvp`   | Use water vapor pressure.  Parameters: :math:`ah-A_H`; hh - Henry's constant modifier :math:`H=P_{wv} \cdot \Delta H_H` |
+----------------+-------------------------------------------------------------------------------------------------------------------------+

Below is an example of an `hparam` block that uses all three of the above models:

.. code::

   hparam
   C6H6	hoff	ah=0.34, dhh=3.2
   Cl2	multi	ah1=4.12, ah2=2.4, ah3=4, ah4=0.18, ah5=1.1
   CO2	wvp	hh=2.8, ah=0.04


In this example, Henry's Law component C6H6 is simulated using the van't Hoff model.  Its A,,H,, value is 0.34, and its ΔH,,H,, value is 3.2.

``diff``
--------

This block sets the molecular diffusion coefficients and diffusion models for liquid
and vapor components.  Due to technical restrictions imposed by FEHM's solver,
there can only be one liquid and one vapor diffusion coefficient for the entire model.
The ``diff`` block consists of one line.  After the keyword are up to four name/value pairs, where the name is `l` to set the liquid diffusion coefficient, `v` to set the vapor diffusion coefficient, `lm` to set the liquid diffusion model, or `vm` to set the liquid diffusion model.  The possible values for the diffusion models are as follows:

+---------------+---------------------------------------------------------------------------------------------------------+
| Option        | Description                                                                                             |
+===============+=========================================================================================================+
| `0` or `con`  | Constant molecular diffusion coefficient                                                                |
+---------------+---------------------------------------------------------------------------------------------------------+
| `1` or `mq`   | Millington-Quirk diffusion model                                                                        |
+---------------+---------------------------------------------------------------------------------------------------------+
| `2` or `cw`   | Conca-Wright diffusion model for liquid, or alternate Millington-Quirk model for vapor                  |
+---------------+---------------------------------------------------------------------------------------------------------+
| `3` or `adif` | Diffusion model calculated as a function of pressure and temperature using tortuosity from `adif` macro |
+---------------+---------------------------------------------------------------------------------------------------------+


Below is an example of the ``diff`` block:

.. code::

   diff	l=1e-9, v=1.6e-8	lm=0, vm=con


In this example, the liquid diffusion coefficient is 1×10^-9^, and the vapor diffusion coefficient is 1.6×10^-8^.  Both liquid and vapor use a constant molecular diffusion coefficient.

``disp``
--------

This block is used to set the dispersivity constants of regions of the simulation.  Dispersivity parameters are applied to dispersivity models, which are applied to zones in the `assign` block.  Model names run down the left side of the block, and parameters run across the top.  Dispersivity can be supplied either for the X, Y, and Z directions, or for the longitudinal and transverse directions.  Parameter names consist of two characters:  the first `l` or `v` for liquid or vapor, and the second `x`, `y`, `z`, `l`, or `t` for X, Y, Z, longitudinal, or transverse, respectively.  Parameter names from the two modes of setting dispersion cannot be mixed.

Below is an example of the ``disp`` block:

.. code::

   disp	lx	ly	lz	vx	vy	vz
   model1	0.2	3	1.5	0.28	3.6	1.92
   model2	0.18	0.9	2.361	1.22	0.56	0.58


In this example, X/Y/Z dispersivity is set.  There are two models, named `model1` and `model2`.  The liquid dispersivity for `model1` is 0.2 in the X direction, 3 in the Y direction, and 1.5 in the Z direction.  The vapor dispersivity for `model1` is 0.28 in the X direction, 3.6 in the Y direction, and 1.92 in the Z direction.

``sorp``
--------

This block sets adsorption parameters for selected components.  On the same line as the block title should appear any or all of the following column headers:

+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Header  | Description                                                                                                                                |
+=========+============================================================================================================================================+
| `ltype` | The type of adsorption model to be used to simulate liquid adsorption.  See below for the possible adsorption models from which to choose. |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `a1l`   | :math:`\alpha_1` parameter in liquid adsorption model.                                                                                     |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `a2l`   | :math:`\alpha_2` parameter in liquid adsorption model.                                                                                     |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `bl`    | :math:`\beta` parameter in liquid adsorption model.                                                                                        |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `vtype` | The type of adsorption model to be used to simulate vapor adsorption.  See below for the possible adsorption models from which to choose.  |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `a1v`   | :math:`\alpha_1` parameter in vapor adsorption model.                                                                                      |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `a2v`   | :math:`\alpha_2` parameter in vapor adsorption model.                                                                                      |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+
| `bv`    | :math:`\beta` parameter in vapor adsorption model.                                                                                         |
+---------+--------------------------------------------------------------------------------------------------------------------------------------------+

If a column header is not specified, it is assumed to be zero for every component.

The following are the available adsorption models to use for ``ltype`` and ``vtype``:

+----------------+---------------------------------------+
| Option         | Description                           |
+================+=======================================+
| `0` or `con`   | Conservative solute                   |
+----------------+---------------------------------------+
| `1` or `lin`   | Linear sorption isotherm              |
+----------------+---------------------------------------+
| `2` or `freu`  | Freundlich sorption isotherm          |
+----------------+---------------------------------------+
| `3` or `mfreu` | Modified Freundlich sorption isotherm |
+----------------+---------------------------------------+
| `4` or `lang`  | Langmuir sorption isotherm            |
+----------------+---------------------------------------+

The :math:`\alpha_1, \alpha_2,` and :math:`\beta` parameters are used differently
according to the adsorption model chosen:

+---------------------+-----------------------------------------------------------+---------------------------------+------------------+-----------------------+
| Model               | Expression                                                | :math:`\alpha_1`                | :math:`\alpha_2` | :math:`\beta`         |
+=====================+===========================================================+=================================+==================+=======================+
| Linear              | :math:`C_r = K_d \cdot C_l`                               | :math:`K_d`                     | 0                | 1                     |
+---------------------+-----------------------------------------------------------+---------------------------------+------------------+-----------------------+
| Freundlich          | :math:`C_r = \Lambda \cdot C_1^\beta`                     | :math:`\Lambda`                 | 0                | :math:`0 < \beta < 1` |
+---------------------+-----------------------------------------------------------+---------------------------------+------------------+-----------------------+
| Modified Freundlich | :math:`C_r / (C_{r,max} - C_r) = \Lambda \cdot C_l^\beta` | :math:`\Lambda \cdot C_{r,max}` | :math:`\Lambda`  | :math:`0 < \beta < 1` |
+---------------------+-----------------------------------------------------------+---------------------------------+------------------+-----------------------+
| Langmuir            | :math:`C_r = (r_b \cdot C_l) / (1 + r \cdot C_l)`         | :math:`r_b`                     | r                | 1                     |
+---------------------+-----------------------------------------------------------+---------------------------------+------------------+-----------------------+

For a more in-depth description of the models used for sorption, please refer to the FEHM Models and Methods Summary.

Models are designated by a line starting with a period ("."), immediately followed by the name of the model and nothing else.  On the next several lines, components specified in `comp` may be given, with their corresponding parameters for the current model in the correct column.  The same components, in the same order, must be given in each model.  If a component is liquid- or vapor-only, then asterisks should be placed in the columns that do not apply to that component.  An asterisk can also be placed in the column containing component names to indicate that all components that have not yet been explicitly given sorption parameters should use the model on the line with the asterisk.  Thus, to assign the same sorption parameters to all components, only one line per model would be supplied, containing an asterisk in the first column.

Below is an example `sorp` block:

.. code::

   sorp	ltype	a1l	a2l	bl	vtype	a1v	a2v	bv
   .model1
   CO3--	freu	3.02	0.061	0.89	*	*	*	*
   C6H6	lin	0.89	3.16	0.2	freu	0.35	4.519	0.688
   CO2	mfreu	1.20	3.31	0.4	con	0.01	2.01	0.61
   .model2
   CO3--	con	0.001	2.02	0.88	*	*	*	*
   C6H6	con	1.20	0.58	1.12	lin	2.2	2.043	2.7
   CO2	mfreu	3.006	1.0	9.8	mfreu	0.229	3.434	2.33


In this example, C6H6 is modeled using the linear sorption isotherm model with a liquid α,,1,, parameter of 0.89 and a vapor α,,1,, parameter of 0.35 in adsorption model `model1`, and modeled using the conservative solute model with a liquid α,,1,, parameter of 1.20 and a vapor α,,1,, parameter of 2.2 in model `model2`.

If ``sorp`` is omitted, it is assumed to contain zeros for all values except the β parameters, which are assumed to be 1.

``cden``
--------

The `cden` block allows the user to input the molecular weights of aqueous and aqueous Henry's Law components and have the code adjust the density of the water according to the concentrations of these components.  It should not be used if `trxn` is preceded by a `cden` macro, as some values may be modified.  `cden` accepts several lines, each consisting of the name of a master species and its molecular weight, separated by spaces, tabs, or commas.

Below is an example of the ``cden`` block:

.. code::

   cden
   HCO3-	61
   CO3--	60
   C6H6	78


In this example, HCO3- is defined to have a molecular weight of 61, CO3-- to have a molecular mass of 60, and C6H6 to have a molecular mass of 78.

If database lookup is enabled (see the `lookup` block below), one of the lines in `cden` may consist of only an asterisk.  If such a line is provided, all components that `lookup` dynamically imports will be inserted into the `cden` block with the appropriate molecular weights.  These imported molecular weights can be overridden by explicitly listing the component and its desired molecular weight on a separate line.

``group``
---------

The information in this block is only used if reactions are enabled (see ``ctrl`` above).  It is used to group aqueous components that take part in rapid kinetic reactions.  On the line below the keyword `group`, place one line for each group.  Each line should contain the names of all aqueous components present in that group, separated by spaces, tabs, or commas.

Below is an example of the `group` block:

.. code::

   group
   H
   C_a Cl_a
   Na_a Ca_a Ca_a2
   U238 Th234


In this example, H is in its own group, C_a and Cl_a are in a group together, Na_a, Ca_a, and Ca_a2 are grouped together, and U238 and Th234 are grouped together.

``cplx``
--------

The `cplx` block allows the user to specify instantaneous reactions that form aqueous complexes from the aqueous master species and non-aqueous components specified in `comp`.  One equation is specified on each line, using a slightly modified version of the standard reaction format detailed in the `rxn` block below.  The left side of the reaction should contain only the name of the aqueous complex, without a number denoting its stoichiometry (the stoichiometry of the aqueous complex must be 1).  The right side should contain stoichiometry/compound pairs as specified by the standard format.  If a compound needs to be removed to make the aqueous complex, negate its stoichiometry.  After the equation, two comma-separated values in the format `variable=value` (not padded by spaces) are required:  `ckeq` for the equilibrium constant of that complex, and `heq` for the enthalpy of the formation of that complex.  On the same line as the keyword `cplx`, the keyword `log10` may be provided to denote that the values for all constants in the block will be given as the base-ten logarithm of the actual values.  Asterisks supplied in place of the equilibrium constant or enthalpy value always signify a zero, even if the `log10` keyword is specified.

If the keyword `equi` is supplied on the same line as the block name, the `equi` block below will be consulted to calculate equilibrium constants as functions of temperature.  In this case, the `log10` keyword and any values given to the right of the equations will be ignored.

If the database lookup option is enabled (see the `lookup` block below), complexes from the `% CPLX` section of the database file may be imported using a line that omits the right-hand side of the equation, and possibly the equilibrium information as well.  If such a line is encountered, the code will search for the named complex in the lookup database, and import the complex's constituents, stoichiometry, and equilibrium-related constants.  The `comp` block will be automatically updated to include all components required for the imported complexes.  Note that the `equi` option cannot be used if database lookup is enabled, as this conflicts with the five-parameter fit used by PHREEQC databases for calculating the equilibrium constant as a function of temperature.  These equilibrium values can be overridden by placing explicit values for ckeq and heq on the same line.

Below is an example of `cplx`:

.. code::

   cplx log10
   CaHCO3+ = 1 Ca++ + 1 HCO3-			ckeq=-13.456, heq=0
   OH- = -1 H+					ckeq=-14, heq=0
   Ca2UO2(CO3)3 = 2 Ca++ + 2 UO2 + 3 HCO3- + -3 H+	ckeq=-20.26, heq=*
   MgHPO4


In this example, there are four complexes.  The complex Ca2UO2(CO3)3 is made by combining 2 Ca++, 2 UO2, and 3 HCO3- and removing 3 H+.  The information for the complex MgHPO4 is automatically imported from the database given in the `lookup` block.  All six equilibrium and enthalpy values are given as the base-ten logarithms of their actual values.  The equilibrium constant for CaHCO3+ is 1.0×10^-13.456^.  The enthalpies of CaHCO3+ and OH- are 1.0×10^0^ = 1, but the enthalpy of Ca2UO2(CO3)3 is zero.

``equi``
--------

This block allows equilibrium constants for aqueous complexes from the `cplx` block above to vary with temperature.  If the keyword `equi` is provided in the `cplx` block above, `trxn` will not read equilibrium constants in `cplx`; instead, the `equi` block will be used to determine the constants as a function of temperature.

Every aqueous complex appearing in `cplx` must be included in `equi`.  For each aqueous complex, a line should be provided that contains only the name of the complex, prefixed by a single period.  On the lines that follow, the equilibrium constants should be specified by placing on each line a temperature in o C, a tab or space, and the equilibrium constant at that temperature.  The number of temperatures need not be the same for each complex.

Below is an example of the `equi` block:

.. code::

   equi
   .CaHCO3+
   0	1e-10
   10	2e-10
   40	4e-9
   .OH-
   0	1e-14
   .CaUO2(CO3)3
   20	1.3e-8
   40	2.64e-6
   60	8.7e-4
   80	2.1e-1


In this example, CaHCO3+ has an equilbirium constant of 1×10^-10^ at 0o C, 2×10^-10^ at 10o C, and 4×10^-9^ at 40o C.  OH- has an equilibrium constant of 1×10^-14^ at all temperatures.

``dist``
--------

This block specifies distribution models that can be used for reaction types 1 and 2 to describe the distribution coefficient as a function of temperature.  For each distribution model, a line beginning with a period and then the name of the model (without a separating space) is given.  This is followed by a set of temperature/distribution coefficient pairs, one per line, for as many lines as desired.  FEHM will perform a piecewise linear interpolation between the values given to determine the value of the distribution coefficient for intervening temperatures.

Below is an example of the `dist` block:

.. code::

   dist
   .model1
   0	1
   10	2
   20	4
   30	8
   40	16
   .model2
   0	20
   50	40
   100	70


There are two models in this example.  The first one, `model1`, has five data points at 10o C intervals from 0 to 40o C.  The distribution coefficient at 0o C is 1, the coefficient at 10o C is 2, and the coefficient at 20o C is 4.  Intermediate temperatures are linearly interpolated, so the distribution coefficient for `model1` at 25o C is 6.

``sol``
-------

The ``sol`` block specifies solubility models that can be used for reaction types 7 and 8.  For each solubility model, a line beginning with a period, immediately followed by the model name, is given.  This is followed by temperature/solubility coefficient pairs, one per line, for as many lines as desired.  FEHM will perform a piecewise linear interpolation between the values given to determine the solubility coefficient for intervening temperatures.

Below is an example of the `sol` block:

.. code::

   sol
   .model1
   0	0
   10	1e-16
   20	1e-4
   30	1e-2
   50	1


This example contains one model, `model1`.  In this model, the solubility coefficient changes from 0 to 1 over a range from 0o C to 50o C.  The solubility coefficient at 10o C is 1×10^-16^, and the coefficient at 20o C is 1×10^-4^.  Because FEHM linearly interpolates between successive values of the solubility coefficient, the coefficient at 40o C is 0.505.

``lookup``
----------

This block enables the dynamic lookup process for mineral dissolution and aqueous complexation.  `lookup` uses a lookup database (generated by `trxndb` from a USGS PHREEQC geochemical database) to determine which reactions occur among a specific set of minerals and aqueous complexes.  The `lookup` block adds data to the `comp`, `cden`, `group`, `cplx`, and `rxn` blocks based on the information provided it.

In order to use `lookup`, a database in `trxn`'s standard format must be supplied.  This format is as follows:

.. code::

   Keyword `% MASTER`\\
   Blank line(s)\\
   Component 1 parameters:  component 1 name, component 1 master species name, component 1 molar mass\\
   Component 2 parameters...\\
   Blank line(s)\\
   Keyword `% CPLX`\\
   Blank line(s)\\
   Complex 1 name\\
   Complex 1 products:  product 1 stoichiometry, product 1 name; product 2 stoichiometry, product 2 name...\\
   Equilibrium constant for this complex\\
   Enthalpy of reaction for this complex\\
   Up to five values defining temperature dependence of equilibrium constant:  A,,1,,, A,,2,,, A,,3,,, A,,4,,, A,,5,,, where log,,10,, K = A,,1,, + A,,2,, · T + A,,3,, / T + A,,4,, · log,,10,, T + A,,5,, / T^2^.\\
   Blank line(s)\\
   Complex 2 parameters...\\
   Blank line(s)\\
   Keyword `% MIN`\\
   Blank line(s)\\
   Mineral 1 name, mineral 1 formula\\
   Mineral 1 products\\
   Equilibrium constant for this mineral\\
   Enthalpy of reaction for this mineral\\
   Up to five values defining temperature dependence of equilibrium constant\\
   Blank line(s)\\
   Mineral 2 parameters...\\
   Blank line(s)\\
   Keyword `% END`\\


Comments can be included anywhere in the input file by using a pound sign ("#").

Below is a brief example of a database in this form.

.. code::

   % MASTER
   
   E	e-	0
   H	H+	1.007942
   Mn(+2)	Mn++	54.938045
   Mn(+3)	Mn+++	54.938045
   F	F-	18.9984032
   Al	AlOH++	26.9815386
   Si	H4SiO4	28.08554
   Mg	Mg++	24.305
   Ca	Ca++	40.0782
   S	HS-	32.0652
   Fe	Fe++	55.845
   
   % CPLX
   
   Al+++
   1 AlOH++ 1 H+ 
   -5.00
   11.49
    -38.253 0.0 -656.27 14.327
   
   Mn++
   1 e- 1 Mn+++ 
   -25.510
   25.800
   
   
   MnF+
   1 Mn++ 1 F- 
   0.840
   
   % IMM
   
   Sepiolite Mg2Si3O7.5OH:3H2O
   3 H4SiO4 -4 H+ 2 Mg++ 
   15.760
   -10.700
   
   Fluorite CaF2
   1 Ca++ 2 F- 
   -10.600
   4.690
   66.348 0.0 -4298.2 -25.271 
   
   
   Pyrite FeS2
   2 HS- -2 e- 1 Fe++ -2 H+ 
   -18.479
   11.300
   
   % END


In this example, there are three complexes and three minerals.  The first complex is Al+++, which is formed by the master species AlOH++ and H+.  Its equilibrium constant is -5.00, its enthalpy of formation is 11.49, and four of five possible temperature-dependence parameters are supplied.  For the complex MnF+, no enthalpy of formation or temperature-dependence parameters are supplied.

A conversion script, ``trxndb``, is available to automate conversion of certain USGS PHREEQC input files to the appropriate format.  The converter script understands a limited subset of the complete syntax used in PHREEQC input files; if it gives improper results or errors, ensure that the input file is in a consistent format and that the keywords `log_k`, `delta_h`, and `-analytic` are used rather than their shortened alternatives.  An example file that can be converted flawlessly by `trxndb` can be found at `/scratch/nts/ms/trxn/geochem/phreeqc.dat`.  The converter script is written in Perl and located at `/scratch/nts/ms/trxn/geochem/trxndb`.  The script should be called with the name of the PHREEQC input file as the only argument.  It will create a file in the same directory with the same name and extension `trxd` containing the `trxn`-compatible database.  For example:

.. code::

   $ ls
   phreeqc.dat
   $ /scratch/nts/ms/trxn/geochem/trxndb phreeqc.dat
   45 master species, 187 solution equations read
   58 mineral equations read
   45 master species, 180 solution equations written
   57 mineral equations written
   $ ls
   phreeqc.dat  phreeqc.trxd
   $


The ``lookup`` block must be the first block in the ``trxn`` macro (excepting the ``ctrl`` block if it is used) and consists of one line, which contains the block name `lookup` and the full path to the database file.  The `lookup` block should be followed by a blank or commented line.

Below is an example of the ``lookup`` block.

.. code::

   lookup /scratch/nts/ms/trxn/geochem/phreeqc.trxd


In this example, the database is located at ``/scratch/nts/ms/trxn/geochem/phreeqc.trxd``.

Once ``lookup`` has been enabled, the ``cden``, ``cplx``, and ``rxn`` blocks may be modified by the user to utilize the information from the lookup database.  Please see the sections for these blocks for the appropriate syntax to take advantage of this information.  If the `debug` option is provided in the `ctrl` block, the final version of each block will be printed to the output (`.out`) file, which may aid in debugging if the results of the run are not as expected.

``rxn``
-------

The ``rxn`` blocks are used to model kinetic reactions between simulated compounds.  Seven types of reactions my be used, each with its own input parameters and input format.  `rxn` is intended to be specified multiple times, once for each reaction that is taking place.  For each reaction, the first line of the block contains the keyword `rxn` and the number representing the type of the reaction (see below).  The next several lines are used for the parameters unique to that reaction type, which are detailed below.  The end of the `rxn` block is signaled by a blank line.

Most of these reactions take one line of input in the standard reaction format, which is used to specify reactants, products, and stoichiometries simultaneously.  In this format, each reactant or product is specified by a number denoting the stoichiometry of that compound, a space, and then the name of the compound as given in `comp`.  Compounds are separated from each other by a plus sign padded on either side with spaces (" + ").  The products and reactants are separated from each other by a token containing an equals sign ("=").  Optionally, directional indicators may be added to the equals sign to indicate the direction of the reaction (e.g., "=>").  The reactants must be placed on the left side of the equals sign, and the products on the right side.  For example:

.. code::

   6 HCl + 2 Al_s => 2 AlCl3 + 3 H2


The reaction types and their parameters are as follows.  A depthier description of the mechanics of each type of reaction can be found at the end of the `rxn` section of the FEHM User's Manual.

+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Number  | Type                               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
+=========+====================================+===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================+
| 1       | Linear kinetic reaction            | This reaction accepts one aqueous reactant and one solid product with 1:1 stoichiometry.  The first parameter line for this reaction is a reaction in standard format, without the stoichiometric coefficients.  The second parameter line in a list of comma-separated name/value pairs, where the acceptable names are `rate` to specify the rate of the reaction and `distcoef` to specify the distribution coefficient.  The value for `distcoef` my be a real number of the coefficient is constant, or the name of a model specified in `dist` (without the beginning period).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2       | Langmuir kinetic reaction          | This reaction's parameter input is identical to the input for reaction type 1 except for an extra available parameter for the name/value pairs.  This parameter is `maxconc`, used to set the maximum concentration that can sorb.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 3       | General kinetic reaction           | This reaction accepts its first line of of input as a generic reaction in the standard format described above.  This is followed by a line of name/value pairs, for which the name can be `forward` to set the forward rate constant for the reaction or `reverse` to set the reverse rate constant.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 4       | Dual Monod biodegradation reaction | The first three to six parameter lines of this reaction consist of a name, a colon, and one or more component/complex names.  The name `substrate` accepts a single immobile component name that is to be the substrate that is degraded.  The name `electronacceptor` accepts the stoichiometry and name (separated by a space) of the aqueous complex that is the electron acceptor for the reaction.  The name `biomass` accepts the stoichiometry and name of the solid component that is the biomass produced by the reaction.  The name `reactants` is optional and accepts the stoichiometries and names of any extra reactants that are participating in the reaction.  Likewise, `products` is optional and accepts the names of any extra products of the reaction.  Only a total of three additional reactants and products can be specified.  If only certain forms of the substrate are biodegradable, those can be listed with the `biodegradable` name.\\The last parameter line contains a list of name/value pairs as follows:  `ks` - the half-maximum-rate concentration of the substrate; `ka` - the half-maximum-rate concentration of the electron acceptor; `decay` - the microbial decay rate (1/hr); `phthreshold` - the pH threshold above which biodegradation will not occur; `qm` - the maximum rate of substrate utilization (mol/kg biomass/hr); `yield` - the microbial yield coefficient (kg biomass/mol substrate); `xminit` - the minimum concentration of biomass (mol/kg rock).                                                                                                                                                          |
+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 5       | Radioactive decay reaction         | This reaction accepts one aqueous component reactant and one aqueous component product with 1:1 stoichiometry.  This first parameter line is a reaction in standard format without stoichiometric coefficients.  The reactant is the component that is decaying, and the product is the decay product.  If the decay product is not being modeled in the simulation, an asterisk may be given in place of the product name.  The second parameter line contains a single name/value pair.  The name is `halflife`, and the value is the half-life of the reaction in years.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 7 and 8 | Precipitation/dissolution reaction |  For reaction type 7, the rates are based on the total aqueous concentration of the components; whereas for reaction type 8, the rates are based on the free-ion concentrations alone.  This reaction accepts one solid component and any number of aqueous master species.  The first line of input for this reaction should be a reaction in standard format, containing only a solid component on one side, and at least one master species on the other.  The second line should contain the following name/value pairs:  `solubility` - the solubility product (either a real number or the name of a solubility model specified in `sol`); `rate` - the reaction rate (mol/m^2^/s); `sarea` - the surface area of the mineral (m^2^/m^3^ rock); `molecularweight` - the molecular weight of the mineral; `density` - the density of the mineral.  `molecularweight` and `density` should only be provided for reaction type 8.\\For reaction types 7 and 8 only, if the database lookup option is enabled (see the `lookup` block above), an alternate form of the reaction can be input.  On the first line of the reaction, provide only the name of a mineral defined in the `% IMM` section of the lookup database.  The information on the reactants and products, as well as the solubility constant, are imported from the database, and any components that are necessary for the reaction are dynamically added to the `comp` block.  However, the `rate` and `sarea` parameters (and `molecularweight` and `density` parameters for reaction type 8) are not contained in the database, and must still be set by the user on the second line of the reaction. |
+---------+------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Please note that in all name/value pairs, the name and value must be separated by only an equals sign that is not padded by spaces.

Here is an example of each reaction type.  Two versions of reaction type 5 are given to demonstrate the optional unsimulated daughter species, and two versions of reaction type 7 are provided to demonstrate the database lookup option.

.. code::

   rxn 1
   Ca_a <=> Ca_s
   rate=10, distcoef=2
   
   rxn 2
   Ca_a <=> Ca_s
   rate=10, maxconc=2, distcoef=model2
   
   rxn 3
   3 H + 2 Ca + 5 U238 <=> Cl + 2 Cl2 + 3 C_a
   forward=2, reverse=1.5
   
   rxn 4
   substrate:  U238
   electronacceptor:  3 H
   biomass:  Ca
   reactants:  2 Cl
   products:  1 Na, 5 Th234
   biodegradable:  UO2
   ks=1.2, ka=1.35, decay=0.69, phthreshold=8.2, qm=0.20, yield=1.2, xminit=0.067
   
   rxn 5
   U238 => Th234
   halflife=20
   
   rxn 5
   U234 => *
   halflife=10
   
   rxn 7
   NaCl <=> Na + Cl
   solubility=0.0231, rate=1.02, sarea=2.2
   
   rxn 7
   Quartz
   rate=0.2, sarea=1
   
   rxn 8
   CaCl2 <=> Ca + 2 Cl
   solubility=model1, rate=0.2, sarea=5, molecularweight=60, density=5.25


``assign``
----------

This block allows the user to assign the parameters stored in models in the above blocks to the zones defined in the `zone` macro.  The `assign` block also allows assignment of some other parameters that are specific to zones.

Zone numbers run down the side of the `assign` block, and parameters run across the top.  The following parameters may be supplied; all are optional:

+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| Option  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | Default Value                                                                                                     |
+=========+===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================+===================================================================================================================+
| `water` | Sets the initial water type filling the zone.  This parameter must be the name of a water type defined in the `water` block.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | Pure water                                                                                                        |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `boun`  | Sets the concentrations of species coming from inflow nodes in this zone.  This parameter may consist of a water type defined in `water` and/or a gas type defined in `gas`.  If both a water type and a gas type are flowing in, they should be separated by a period (not padded by spaces).  If there is no inflow in this zone, give an asterisk for this parameter.  An asterisk can also be used to specify inflow of pure water.                                                                                                                                                                                                                                                                                                                                                                                                                                   | Pure water                                                                                                        |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
|`time`   | Sets the time range during which the inflow nodes in this zone are injecting, in days.  A lone zero ("0") gives no inflow.  An asterisk provides no inflow if the entry in the `boun` column is "*" or if the `boun` column is missing, but provides inflow over the entire simulation otherwise.  A single number other than zero will give inflow for one day starting at the specified day.  A range separated by a greater-than sign (">") not padded by spaces will run injection from the first number specified to the second number.  If the number before or after the greater-than is ommitted, it will default to the beginning or end of the simulation, respectively.  Thus, ">30" will run injection from time 0 through time 30, "30>" will run injection from time 30 to the end of the simulation, and ">" will run injection for the entire simulation. | Inflow for the entire simulation if valid water/gas types are specified in the `boun` column, no inflow otherwise |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `rock`  | Sets the composition of the rock in this zone.  This parameter must be the name of a rock type defined in `rock`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | Rock contains no relevant species                                                                                 |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `gas`   | Sets the initial composition of the gas in this zone.  This parameter must be the name of a gas type defined in `gas`.  If there is no gas in this zone, give an asterisk for this parameter.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | No gas                                                                                                            |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `disp`  | Sets the dispersivity constants for this zone.  This parameter must be the name of a dispersivity model defined in `disp`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | No dispersivity                                                                                                   |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `sorp`  | Sets adsorption parameters for this zone.  This parameter must be the name of an adsorption model defined in `sorp`.  While the models in `sorp` are defined by starting the line with a period, do not include the period in this parameter.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | No adsorption for any components                                                                                  |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `tpor`  | Sets the optional tracer porosity for this zone.  This parameter must be a real number from 0 to 1.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | 0.32                                                                                                              |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| `opt`   | Sets miscellaneous options per-zone.  Options should be separated by periods with no spaces.  The available options are:  `const` causes the concentrations of solutes in the inflow water to be held constant at nodes in the current zone; `accum` enables solute accumulation in the current zone.  Note that `const` and `accum` are mutually exclusive.  An asterisk can be used to specify no options.                                                                                                                                                                                                                                                                                                                                                                                                                                                              | No options                                                                                                        |
+---------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+

An asterisk may be given for any of the parameters above for a given zone, in which case the default values given above are used.  If an entire column is omitted, that parameter will be given the default values shown above for every zone.  If a zone is omitted, it will receive default values for every column.  Nodes that are not in any zone at the time when `trxn` is called will also receive default values.

Below is an example of the ``assign`` block.

.. code::

   assign	water	rock	gas	boun	time	disp	sorp	tpor
   1	wt3	granite	*	vt1.wt1	20>30	model1	model1	0.28
   2	wt2	clay	vt2	*	*	model2	model2	0.69
   3	wt2	granite	vt2	*	*	model1	model2	0.32


In this example, there are three zones.  The first zone is initially filled with water of type `wt3` and no gases, and water of type `wt1` and gas of type `vt1` are flowing into it starting at time 20 days and ending at time 30 days.  Its rock is of type `granite`, it uses the dispersivity parameters defined in `model1` and the adsorption parameters defined in `model1`, and has a tracer porosity of 0.28.  Note also that zones 2 and 3 have no inflow.

After specifying all applicable blocks, use the `end trxn` keyword to end the reading of the `trxn` macro.

Example
-------

Below is a complete, commented working example `trxn` macro, with its accompnaying zone macro for reference.   This example was taken from the `multi_solute` test case in FEHM's standard verification suite. 


.. code::

   zone
   default
   all
   bound
   nnum
   2	1 102
   
   trxn
   
   # In the zone macro above, all nodes are placed into the zone named "default",
   # except for nodes 1 and 202, which are placed in the zone named "bound".
   # Zone "bound" will be used for inflow.
   
   ctrl rxnon # Enable reactions.
   
   # Include header information from trac.
   header
   1.d-80	1.0	1.e-6	0.5
   1.	2000	1.0	2000
   5	5.0	1.e-6	2.8333e-3	0
   iskip=0, rsdmax=1e-10
   
   # There are six components in this simulation:  aqueous cobalt, iron, and EDTA,
   # and solid Co-EDTA, Fe-EDTA, and cobalt.
   comp		master
   a Cobalt	Cobalt_a
   a Iron		Iron_a
   a EDTA		EDTA_a
   s Co-EDTA_s	*
   s Fe-EDTA_s	*
   s Cobalt_s	*
   
   # There is only one type of water, called "inflow", which contains 3.16e-5 M
   # aqueous cobalt, 1e-13 M aqueous iron, and 3.16e-5 M aqueous EDTA.
   water	Cobalt	Iron	EDTA
   inflow	3.16e-5	1e-13	3.16e-5
   
   # There is one sorption model.  It models all components with a linear sorption
   # isotherm, using alpha-1 and alpha-2 parameters of zero, and a beta parameter
   # of 1.
   sorp	ltype	a1l	a2l	bl
   .smod
   Cobalt	lin	0	0	1	
   Iron	lin	0	0	1
   EDTA	lin	0	0	1
   
   # We assign the liquid diffusion coefficient for the simulation to be 1e-9.
   diff	l=1e-9
   
   # There is one model for dispersivity, "dmod".  It sets the dispersivity to
   # 0.05 in the X direction, and 1e-34 in the Y and Z directions.
   disp	lx	ly	lz
   dmod	0.05	1e-34	1e-34
   
   # There are two groups for the coupled solver.  One contains cobalt and EDTA,
   # and the other contains iron.
   group
   Cobalt EDTA
   Iron
   
   # Here, two aqueous complexes are defined:  Co-EDTA and Fe-EDTA.  Co-EDTA is
   # composed of one aqueous cobalt and one aqueous EDTA; Fe-EDTA is composed of
   # one aqueous iron and one aqueous EDTA.  Both complexes have enthalpy zero;
   # Co-EDTA has equilibrium constant 1e18 and Fe-EDTA has equilibrium constant
   # 6.31e27.
   cplx
   Co-EDTA_a = 1 Cobalt_a + 1 EDTA_a	ckeq=1e18, heq=0
   Fe-EDTA_a = 1 Iron_a + 1 EDTA_a		ckeq=6.31e27, heq=0
   
   # There are four reactions taking place in this simulation.
   
   # Reaction 1 is a linear kinetic reaction describing the dissolution of
   # Co-EDTA.  The distribution coefficient is 0.533, and the rate of reaction is
   # 1.
   rxn 1
   Co-EDTA_a <=> Co-EDTA_s
   distcoef=0.533, rate=1
   
   # Reaction 2, also linear kinetic, describes the dissolution of cobalt.  The
   # distribution coefficient is 5.07, and the rate of reaction is again 1.
   rxn 1
   Cobalt_a <=> Cobalt_s
   distcoef=5.07, rate=1
   
   # Reaction 3 is also linear kinetic and describes the dissolution of Fe-EDTA.
   # The distribution coefficient here is 0.427, and the rate of reaction is 1.
   rxn 1
   Fe-EDTA_a <=> Fe-EDTA_s
   distcoef=0.427, rate=1
   
   # Reaction 4 is a general kinetic reaction that describes the complexation of
   # solid cobalt and Fe-EDTA to form Co-EDTA.  The forward rate constant is
   # 1.26e-2, and the reaction never occurs in reverse (the reverse rate constant
   # is zero).
   rxn 3
   Co-EDTA_s = Fe-EDTA_s + Cobalt_s
   for=1.26e-2, rev=0
   
   # Finally, attributes are assigned to zones.  The first zone, "default", con-
   # tains most of the nodes; the zone "bound" contains only the inflow nodes.
   # The initial water filling both zones is pure water, as signified by the
   # asterisks in the water column, and the rock does not contain any species that
   # participate in the reactions, as signified by the lack of a rock column.
   # (The water column could also have been left out entirely, but is included
   # here for clarity.)  No inflow occurs in the default zone, as shown by the
   # asterisk in the boun column, but water of the  "inflow" type defined in the
   # water block is flowing in through the nodes in the bound zone from time 1
   # through time 4.167 of the simulation.  The dispersivity model "dmod" and the
   # sorption model "smod" are applied to both zones.
   assign	water	boun	time	disp	sorp
   default	*	*	0	dmod	smod
   bound	*	inflow	1>4.167	dmod	smod
   
   end trxn


Below are the ``trac`` and ``rxn`` macros that were replaced by the above ``trxn``, for reference and comparison.


.. code::

   trac
   1.d-80 1.0 1.e-6 0.5
   1. 2000 1.0 2000
   5 5.0 1.e-6 2.8333e-3
   6
   1
   1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34
   
   1 202 1 1
   
   
   1 202 101 3.1623e-5 1.0 4.16667
   
   1
   1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34
   
   1 202 1 1
   
   
   1 202 101 1.e-13 1.0 4.16667
   
   1
   1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34
   
   1 202 1 1
   
   
   1 202 101 3.1623e-5 1.0 4.16667
   
   0
   
   
   0
   
   
   0
   
   
   end trac
   rxn
   ** NCPLX, NUMRXN
   2,4
   ** Coupling of the aqueous components (dRi/dUj)
   2
   1 0 1
   0 1 0
   ** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT rows)
   1    Cobalt[aq]     0     0     1.e-9
   2      Iron[aq]     0     0     1.e-9
   3      EDTA[aq]     0     0     1.e-9
   ** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
   101   Co-EDTA[aq]   0
   102   Fe-EDTA[aq]   0
   ** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immoblie spec, NIMM rows)
   1   Co-EDTA[s]     0
   2   Fe-EDTA[s]     0
   3    Cobalt[s]     0
   ** IDVAP(IV), VAPNAM(IM), VAPPRT(IV) (ID # and name of vapor spec, NVAP rows)
   ** Skip nodes
   0
   ** RSDMAX
   1.0e-10
   **** Chemical reaction information for equilibrium reactions ******
   ** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
            0
   ** CKEQ(IX) ,HEQ(IX) (Stability constants and Enthaplys, NCPLX rows)
   1.0e+18    0
   6.31e+27   0
   ** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
   1.0       0.0       1.0
   0.0       1.0       1.0
   ** LINEAR KINETIC REACTION (type 1) **
           1
   ** Where does the reaction take place? **
    1 0 0
   
   ** Aqueous Component/Complex #, Solid Component #
         101      1
   ** Distribution coeffienct (kg water/ kg rock) **  
       0.533
   ** Mass transfer coefficient (1/hr) **
         1.0
   ** LINEAR KINETIC REACTION (type 1) **
           1
   ** Where does the reaction take place? **
    1 0 0
   
   ** Aqueous Component/Complex #, Solid Component #
           1      3
   ** Distribution coeffienct (kg rock/ kg water) **  
        5.07
   ** Mass transfer coefficient (1/hr) **
         1.0
   ** LINEAR KINETIC REACTION (type 1) **
           1
   ** Where does the reaction take place? **
    1 0 0
   
   ** Aqueous Component/Complex #, Solid Component #
         102      2
   ** Distribution coeffienct (kg rock/ kg water) **  
       0.427
   ** Mass transfer coefficient (1/hr) **
         1.0
   ** GENERAL EXCHANGE REACTION (type 3) **
           3
   ** Where does the reaction take place? **
    1 0 0
   
   ** # of solid, liquid and vapor species **
           3   0   0
   ** forward and reverse rate constants (1/hr) **
     1.26e-2 0
   ** Solid Species in reaction **
     1      2      3
   ** Stoichiometry **
     1.0   -1.0  -1.0
   end rxn


Additional Notes
----------------

Removed Features
----------------

Please note that the following features present in `trac` and `rxn` have been removed in `trxn`:

* ``trac``
 * Inflow concentrations for solids
 * Varying molecular diffusion coefficient by either species or location
 * Setting different dispersivities for different species
* ``rxn``
 * Using `IFXCONC` to specify that concentrations are free-ion only (all concentrations in `water`, `rock`, and `gas` must be total aqueous concentrations)
 * Enabling reactions at specific nodes only
 * Reaction type 6
* Both
 * JA JB JC format for entering specific nodes for tracer injection, etc.  (These must be defined by zones.)

If these features are desired, old-style input using `trac` and `rxn` must be used instead.

Further resources and verification
----------------------------------

Several test problems from the standard FEHM test suite that have `trac` and/or `rxn` macros have been converted to the `trxn` format.  These can be found in `/scratch/nts/ms/trxn/test-problems`.  The input files contain the original `trac` and/or `rxn` macros (turned off), along with an equivalent `trxn` macro.  The output for the `trxn` macro has been verified against the output for the original macros in each of these cases; the results of the comparison can be found in the `plot/plot.png` and `plot/plot-orig.png` files in the test directory.  The location of the input file in the test directory is included in parentheses after the test problem name.  The following test problems have been verified to work successfully with `trxn`:

* `baro_trans` (`input/baro_trans.in`)
* `cden_test (`input/static-multi1.dat`)`
* `dissolution` (`input/dissolution.in`)
* `fracture_transport` (`input/tangtestN.in`)
* `henrys_law` (`input/henrytest.in`)
* `multi_solute` (`input/multi_solute.in`)
* `sorption` (`input/sorption.in`)
* `transport3D` (`input/3d_trac.dat`)

In addition, several new test problems have been developed in order to better test the full range of `trxn`'s functionality.  These are also found in `/scratch/nts/ms/trxn/test-problems`, and documentation of the problem setup and expected results can be found in `/scratch/nts/ms/trxn/test-problems/meta/doc`.  The following test problems fall into this category:
* `multirock`
* `inflow`
* `decay`

Bugs, error handling, and limitations
-------------------------------------

* `trxn` will try to print out a nice error message if something goes wrong.  However, this is not guaranteed.
* A `zone` macro with at least one valid zone is required before `trxn`.  This `zone` macro may contain no more than 100 zones.
* If a block is specified multiple times, the values from the last block should be used; however, do not rely on this feature.
* The maximum permitted length of input lines is 200 characters.
* The maximum number of characters in a name (of a component, zone, model, etc.) is 40.
* The maximum number of specifications (models, components, etc.) that any given block may contain is 100.
* The zone macro preceding `trxn` may not contain more than 100 zones.
* Model and species names should be strictly alphanumeric plus the five characters "(", ")", "+", "_", and "-". Use of other characters may cause incorrect behavior.
* If only the `trac` blocks of the macro are given, `trxn` should be compatible with an old `rxn` macro, as long as the `trxn` macro is read before the `rxn` macro.  However, this has not been tested and may not work reliably.  Furthermore, `trxn` is not compatible with the old `trac` macro.

Debug tools
-----------

* The keyword `debug` in the `ctrl` block will enable some informational output that may be useful for debugging problems.  The `stop` keyword in the `ctrl` block will halt FEHM immediately after reading and processing `trxn`.
* `null` blocks are ignored by default, but are printed if debugging output is enabled.

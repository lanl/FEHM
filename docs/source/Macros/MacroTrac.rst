========
``trac``
========

* Group 1 - KEYWORD ‘userc', ANO, AWC, EPC, UPWGTA

  * Optional keyword "file" is used to specify the name of the data file that will contain input for the userc subroutine.

  * KEYWORD ‘file' 

  * USERC_FILENAME

or

* Group 1 - ANO, AWC, EPC, UPWGTA

* Group 2 - DAYCS, DAYCF, DAYHF, DAYHS

* Group 3 - IACCMX, DAYCM, DAYCMM, DAYCMX, NPRTTRC

* Group 4 - KEYWORD ‘tpor'

* Group 5 - JA, JB, JC, PS_TRAC

Tracer porosity is entered only if the Group 4 keyword (‘tpor'), which specifies tracer porosity input, is present, otherwise Groups 4 and 5 are omitted.

* Group 6 - NSPECI

* Group 7 - KEYWORD ‘ldsp'

The Group 7 keyword (‘ldsp') specifies longitudinal / transverse dispersion should be used. If X, Y, Z dispersion is desired Group 7 is omitted, and dispersivities are input in X, Y, Z order in Group 9 or Group 12. When longitudinal / transverse dispersion is invoked the Z-components of dispersivity are omitted from the Group 9 or Group 12 input, and X and Y represent longitudinal and transverse dispersion respectively. Note that an "L" or "V" added to the Group 9 or Group 12 variable names (MFLAG, SEHDIFF, TCX, TCY, TCZ, IADSF, A1ADSF, A2ADSF, BETADF, DIFFM) indicates the value is for the liquid or vapor phase, respectively.

* Group 8 - KEYWORD ‘dspl' or ‘dspv' or ‘dspb'

The Group 8 keyword specifies that the same diffusion coefficient and dispersivities are to be used for all species of the same type (liquid and/or vapor). This will make the calculations more efficient and thus should be used if applicable. If Group 8 is omitted, Groups 9 and 10 are also omitted, and input resumes with Group 11.

If only liquid species are present (keyword ‘dspl') or only vapor species are present (keyword ‘dspv') with no longitudinal / transverse dispersion, Group 9 is defined as follows: 

* Group 9 - MFLAG, SEHDIFF, TCX, TCY, TCZ

Otherwise if both liquid and vapor are present (keyword ‘dspb'), parameters for both must be entered. 

* Group 9 - MFLAGL, SEHDIFFL, TCLX, TCLY, TCLZ, MFLAGV, SEHDIFFV, TCVX, TCVY, TCVZ

* Groups 9 is used to define transport models for which diffusion and dispersion parameters are identical. Group 9 is read in until a blank line is encountered. The model number is incremented by 1 each time a line is read.

* Group 10 -  JA, JB, JC, ITRCDSP

* Group 11 - ICNS [SPNAM]

There are two options for group twelve. If the same diffusion coefficient and dispersivities are to be used for all species of the same type (liquid and/ or vapor - keyword ‘dspl', ‘dspv', or ‘dspb') only sorption parameters are input:

* Group 12 - IADSF, A1ADSF, A2ADSF, BETADF

or for a Henry's Law Species (both liquid and vapor)

* Group 12 - IADSFL, A1ADSFL, A2ADSFL, BETADFL, IADSFV, A1ADSFV, A2ADSFV, BETADFV

In the absence of a Group 8 keyword ("dspl', ‘dspv', or ‘dspb') the following input (for liquid or vapor) which includes the sorption and dispersion parameters is used: 

* Group 12 - IADSF, A1ADSF, A2ADSF, BETADF, MFLAG, DIFFM, TCX, TCY, TCZ

For a Henry's Law Species (both liquid and vapor) if DIFFML ≥ 0

* Group 12 - IADSFL, A1ADSFL, A2ADSFl, BETADFL, MFLAGL, DIFFML, TCLX, TCLY, TCLZ, IADSFV, A1ADSFV, A2ADSFV, BETADFV, MFLAGV, DIFFMV, TCVX, TCVY, TCVZ

* Group 13 - JA, JB, JC, ITRCD

* Group 14 - HENRY_MODEL, HAWWA(1), HAWWA(2), HAWWA(3), HAWWA(4), HAWWA(5) (only input for a Henry's Law species, otherwise omitted)

* Group 15 - JA, JB, JC, ANQO

* Group 16 - JA, JB, JC, CNSK, T1SK, T2SK

Groups 11, 12, 13, 14, 15, and 16 are entered as a unit for each solute. However, for a solid species, only groups 11, 15, and 16 are entered (groups 12, 13, and 14 are not applicable for a solid species). Groups 12 and 13 are used to define transport models for which sorption, diffusion and dispersion parameters are identical. For a liquid or vapor species, only one set of Group 12 parameters should be entered per region. However, for a Henry's Law species, two sets of parameters per region must be entered. For this case, the liquid sorption parameters should be entered on the first line and the vapor sorption parameters on a second line or as a continuation of the first line. Group 12 is read in until a blank line is encountered. The model number is incremented by 1 each time a line is read. Group 13 then assigns a transport model number to every node.

Injection nodes must be specified in control statement **flow**.

+----------------+--------------+-------------------------------------------------------------------------------+
| Input Variable | Format       | Description                                                                   |
+================+==============+===============================================================================+
| KEYWORD        | character*5  | Keyword for invoking a solute transport user subroutine.                      |
|                |              | If the word ``userc`` is placed in this position, then the                    |
|                |              | code invokes a solute transport user subroutine at each time                  |
|                |              | step. Omit this key word if there is no solute user subroutine                |
|                |              | for the simulation.                                                           |
+----------------+--------------+-------------------------------------------------------------------------------+
| KEYWORD        | character*4  | Optional keyword ‘file' designating the name for the user                     |
|                |              | subroutine input transport parameter file will be input.                      |
|                |              | If this keyword and the following line are omitted, the                       |
|                |              | default name will be ``userc_data.dat``.                                      |
+----------------+--------------+-------------------------------------------------------------------------------+
| USERC_FILENAME | character*80 | Name of file from which to read transport parameters for                      |
|                |              | optional user subroutine.                                                     |
+----------------+--------------+-------------------------------------------------------------------------------+
| ANO            | real         | Initial solute concentration, set at all nodes for all                        |
|                |              | species unless overwritten by a restart file input or                         |
|                |              | values in group 14 below (moles/kg fluid).                                    |
+----------------+--------------+-------------------------------------------------------------------------------+
| AWC            | real         | | Implicitness factor for solute solution.                                    |
|                |              | | AWC > 1.0 gives 2nd order solution                                          |
|                |              | | AWC ≤ 1.0 gives 1st order solution                                          |
+----------------+--------------+-------------------------------------------------------------------------------+
| EPC            | real         | Equation tolerance for solute solution. When the square                       |
|                |              | root of the sum of the squared residuals is lower than                        |
|                |              | EPC, the solution is assumed to be converged.                                 |
+----------------+--------------+-------------------------------------------------------------------------------+
| UPWGTA         | real         | | Upstream weighting term for the solute solution.                            |
|                |              | | UPWGTA < 0.5 UPWGTA is set to 0.5                                           |
|                |              | | UPWGTA > 1.0 UPWGTA is set to 1.0                                           |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYCS          | real         | Time which the solute solution is enabled (days).                             |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYCF          | real         | Time which the solute solution is disabled (days).                            |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYHF          | real         | Time which the flow solution is disabled (days).                              |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYHS          | real         | Time which the flow solution is enabled (days).                               |
+----------------+--------------+-------------------------------------------------------------------------------+
| IACCMX         | integer      | Maximum number of iterations allowed in solute                                |
|                |              | solution if time step multiplier is enabled                                   |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYCM          | real         | Time step multiplier for solute solution                                      |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYCMM         | real         | Initial time step for solute solution (days)                                  |
+----------------+--------------+-------------------------------------------------------------------------------+
| DAYCMX         | real         | Maximum time step for solute solution (days)                                  |
+----------------+--------------+-------------------------------------------------------------------------------+
| NPRTTRC        | integer      | Print-out interval for solute information. Data for                           |
|                |              | every NPRTTRC solute time step will be written to the                         |
|                |              | ".trc" file. If this parameter is omitted (for                                |
|                |              | compatibility with old input files) the default                               |
|                |              | value is 1. Note that the first and last solute                               |
|                |              | time step within a heat and mass transfer step                                |
|                |              | automatically get printed.                                                    |
+----------------+--------------+-------------------------------------------------------------------------------+
| KEYWORD        | character*4  | Keyword ‘tpor' specifying optional tracer porosity                            |
|                |              | should be input. If group 4 is omitted, porosities                            |
|                |              | assigned in macro rock are used.                                              |
+----------------+--------------+-------------------------------------------------------------------------------+
| PS_TRAC        | real         | Tracer porosity                                                               |
+----------------+--------------+-------------------------------------------------------------------------------+
| NSPECI         | integer      | Number of solutes simulated.                                                  |
+----------------+--------------+-------------------------------------------------------------------------------+
| KEYWORD        | character*4  | Keyword ‘ldsp' specifying longitudinal / transverse dispersion.               |
|                |              | If x, y, z dispersion is desired group 7 is omitted, and                      |
|                |              | dispersivities are input in x, y, and then z order                            |
|                |              | (group 9 or group 12). Otherwise, if longitudinal                             |
|                |              | / transverse dispersion is desired the keyword                                |
|                |              | ‘ldsp' is entered and dispersivities are instead                              |
|                |              | input in longitudinal and then transverse order                               |
|                |              | with values for the third dimension omitted.                                  |
+----------------+--------------+-------------------------------------------------------------------------------+
| KEYWORD        | character*4  | | Keyword specifying the same diffusion coefficient and                       |
|                |              |   dispersivities are to be used for all species of the same                   |
|                |              |   type (liquid and/or vapor).                                                 |
|                |              | | ``dspl`` indicates that only liquid species exist.                          |
|                |              | | ``dspv`` indicates that only vapor species exist.                           |
|                |              | | ``dspb`` indicates that both liquid and vapor species exist.                |
+----------------+--------------+-------------------------------------------------------------------------------+
| ICNS           | integer      | | Phase designation for the ith solute                                        |
|                |              | | * -2 - Henry's Law species (input and output concentration                  |
|                |              |   values are gas concentrations).                                             |
|                |              | | * -1 - Vapor species.                                                       |
|                |              | | * 0 - Solid species                                                         |
|                |              | | * 1 - Liquid species                                                        |
|                |              | | * 2 - Henry's Law species (input and output concentration                   |
|                |              |   values are liquid concentrations)                                           |
+----------------+--------------+-------------------------------------------------------------------------------+
| SPNAM          | character*20 | For each species, the name of the species (e.g. Sulfate).                     |
|                |              |  This is an optional identifier that may be input when                        |
|                |              |  macro ``rxn`` is not being used.                                             |
+----------------+--------------+-------------------------------------------------------------------------------+
| MFLAG          | integer      | | Flag denoting type of diffusion model to be used                            |
|                |              | | 0 - the molecular diffusion coefficient is a constant.                      |
|                |              | | 1 - Millington Quirk diffusion model for liquid or vapor.                   |
|                |              | | 2 - Conca and Wright diffusion model for liquid, alternate                  |
|                |              |   Millington Quirk diffusion model for vapor.                                 |
|                |              | | 3 - vapor diffusion coefficient is calculated as a function                 |
|                |              |   of pressure and temperature using tortuosity from adif                      |
|                |              |   macro, of the "Models and Methods Summary" of the FEHM Application          |
|                |              |   (Zyvoloski et al. 1999).                                                    |
|                |              | |                                                                             |
|                |              | | FEHM calculates liquid contaminant flux as                                  |
|                |              |   J = (Water Content)x(D*)x(GradC) and vapor contaminant flux as              |
|                |              |   J = (Air Content)x(D*)x(GradC) where D\* is the diffusion                   |
|                |              |   coefficient input in this macro. Water content is defined as                |
|                |              |   porosity x saturation and air content is defined as                         |
|                |              |   porosity x (1 - saturation). For more explanation on the Millington         |
|                |              |   Quirk and Conca/Wright models see Stauffer, PH, JA Vrugt, HJ Turin,         |
|                |              |   CW Gable, and WE Soll (2009) **Untangling diffusion from advection          |
|                |              |   in unsaturated porous media: Experimental data, modeling, and               |
|                |              |   parameter uncertainty assessment**. Vadose Zone Journal,                    |
|                |              |   8:510-522, doi:10.2136/vzj2008.0055.                                        |
|                |              | |                                                                             |
|                |              | | SEHDIFF real Molecular diffusion coefficient (m2/s)                         |
|                |              | |                                                                             |
|                |              | | When MFLAG = 0, the input diffusion coefficient is used                     |
|                |              |   directly in the contaminant flux equations presented above.                 |
|                |              |   However, MFLAG = 1 or 2, the free air or free water diffusion               |
|                |              |   coefficient is input and the correct porous diffusion is                    |
|                |              |   calculated within FEHM. For MFLAG = 3, the code assumes a                   |
|                |              |   free air diffusion coefficient of 2.33e-5 m2/s for water                    |
|                |              |   vapor in air as described in the Models and Methods                         |
|                |              |   Summary Eq. 21.                                                             |
+----------------+--------------+-------------------------------------------------------------------------------+
| SEHDIFF        | real         | Molecular diffusion coefficient (m^2^/s)                                      |
+----------------+--------------+-------------------------------------------------------------------------------+
| TCX            | real         | Dispersivity in x-direction (m)                                               |
+----------------+--------------+-------------------------------------------------------------------------------+
| TCY            | real         | Dispersivity in y-direction (m)                                               |
+----------------+--------------+-------------------------------------------------------------------------------+
| TCZ            | real         | Dispersivity in z-direction (m)                                               |
+----------------+--------------+-------------------------------------------------------------------------------+
| ITRCDSP        | integer      | Region number for dispersion parameters given in group 9                      |
|                |              | (keyword dspl, dspv, or dspv). Default is 1.                                  |
+----------------+--------------+-------------------------------------------------------------------------------+
| IADSF          | integer      | | Adsorption model type for the ith species, ith region                       |
|                |              | | 0 - conservative solute                                                     |
|                |              | | 1 - linear sorption isotherm                                                |
|                |              | | 2 - Freundlich sorption isotherm                                            |
|                |              | | 3 - Modified Freundlich sorption isotherm                                   |
|                |              | | 4 - Langmuir sorption isotherm                                              |
+----------------+--------------+-------------------------------------------------------------------------------+
| A1ADSF         | real         | α1 parameter in adsorption model                                              |
+----------------+--------------+-------------------------------------------------------------------------------+
| A2ADSF         | real         | α2 parameter in adsorption model                                              |
+----------------+--------------+-------------------------------------------------------------------------------+
| BETADF         | real         | β parameter in adsorption model                                               |
+----------------+--------------+-------------------------------------------------------------------------------+
| DIFFM          | real         | Molecular diffusion coefficient (m^2^/s) See discussion                       |
|                |              | for SEHDIFF.                                                                  |
+----------------+--------------+-------------------------------------------------------------------------------+
| ITRCD          | integer      | Region number for group 12 sorption parameters or                             |
|                |              | for sorption and dispersion parameters (no keyword).                          |
|                |              | Default is 1.                                                                 |
+----------------+--------------+-------------------------------------------------------------------------------+
| HENRY_MODEL    | integer      | | Flag denoting which model is to be used for                                 |
|                |              |   defining the temperature dependence of the Henry's law constant             |
|                |              | | 1 - van't Hoff model                                                        |
|                |              | | 2 - Multi-parameter fit to experimental data (used for carbonate system)    |
|                |              | | 3 - Henry's model uses water vapor pressure (:math:`H = P_{wv}`)            |
+----------------+--------------+-------------------------------------------------------------------------------+
| HAWWA(1)       | real         | | Term in Henry's Law temperature dependence model:                           |
|                |              | | For model 1 or 3 - parameter value is :math:`A_H`                           |
|                |              | | For model 2 - parameter value is :math:`A_{H,1}`                            |
|                |              | | For model 3 - not used                                                      |
+----------------+--------------+-------------------------------------------------------------------------------+
| HAWWA(2)       | real         | | Term in Henry's Law temperature dependence model:                           |
|                |              | | For model 1 - parameter value is :math:`\Delta H_H`                         |
|                |              | | For model 2 - parameter value is :math:`A_{H,2}`                            |
|                |              | | For model 3 - Henry's constant modifier,                                    |
|                |              |   :math:`H = P_{wv} \cdot \Delta H_H`                                         |
+----------------+--------------+-------------------------------------------------------------------------------+
| HAWWA(3)       | real         | | Term in Henry's Law temperature dependence model:                           |
|                |              | | For model 1 - not used                                                      |
|                |              | | For model 2 - parameter value is :math:`A_{H,3}`                            |
+----------------+--------------+-------------------------------------------------------------------------------+
| HAWWA(4)       | real         | | Term in Henry's Law temperature dependence model:                           |
|                |              | | For model 1 - not used                                                      |
|                |              | | For model 2 - parameter value is :math:`A_{H,4}`                            |
+----------------+--------------+-------------------------------------------------------------------------------+
| HAWWA(5)       | real         | | Term in Henry's Law temperature dependence model:                           |
|                |              | | For model 1 - not used                                                      |
|                |              | | For model 2 - parameter value is :math:`A_{H,5}`                            |
+----------------+--------------+-------------------------------------------------------------------------------+
| ANQO           | real         | Initial concentration of tracer, which will supersede                         |
|                |              | the value given in group 1. Note that if initial                              |
|                |              | values are read from a restart file, these values                             |
|                |              | will be overwritten. Units are moles per kg vapor or                          |
|                |              | liquid for a liquid, vapor, or Henry's law species,                           |
|                |              | and moles per kg of solid for a solid species.                                |
|                |              | Default is 0.                                                                 |
+----------------+--------------+-------------------------------------------------------------------------------+
| CNSK           | real         | Injection concentration at inlet node (moles per kg                           |
|                |              | liquid or vapor). If fluid is exiting at a node, then                         |
|                |              | the in-place concentration is used. If CNSK < 0, then                         |
|                |              | the concentration at that particular node will be held                        |
|                |              | at a concentration of abs(cnsk) (default is 0 for all                         |
|                |              | unassigned nodes).                                                            |
+----------------+--------------+-------------------------------------------------------------------------------+
| T1SK           | real         | Time (days) when tracer injection begins. Default is 0.                       |
+----------------+--------------+-------------------------------------------------------------------------------+
| T2SK           | real         | Time (days) when tracer injection ends. Default is 0.                         |
|                |              | If T2SK < 0, the absolute value of T2SK is used for this                      |
|                |              | parameter, and the code interprets the negative value as a                    |
|                |              | flag to treat the node as a zero-solute-flux node for cases                   |
|                |              | in which a fluid sink is defined for that node. For this                      |
|                |              | case, the solute will stay within the model at the node                       |
|                |              | despite the removal of fluid at that location.                                |
|                |              | If a fluid source is present at the node, CNSK is                             |
|                |              | the concentration entering with that fluid, as in the                         |
|                |              | normal implementation of a solute source. Note that the                       |
|                |              | code cannot handle the case of T2SK < 0 and CNSK < 0                          |
|                |              | (fixed concentration), as these are incompatible inputs.                      |
|                |              | Therefore, the code prints an error message and stops                         |
|                |              | for this condition.                                                           |
+----------------+--------------+-------------------------------------------------------------------------------+


---
title : Program Considerations
layout : page_getting-started
permalink: /program-specification/program-considerations
hero_height: is-hidden
---

# Program Considerations

## Program Options

The uses and capabilities of FEHM are summarized below with reference to the macro input structure discussed in [Input Data](input-data)

## FEHM Capabilities and Associated Macro Command Statements


* Mass, energy balances in porous media
  * Variable rock properties (```rock```)  
  * Variable permeability (```perm```, ```fper```)  
  * Variable thermal conductivity (```cond```, ```vcon```)  
  * Variable fracture properties, dual porosity, dual porosity/dual permeability (```dual```, ```dpdp,```, ```gdpm```)    
* Multiple components available
  * Air-water isothermal mixture available (```airwater```, ```bous```, ```head```), fully coupled to heat and mass transfer (<strong>ngas</strong>, vapl, adif)  
  * Up to 10 solutes with chemical reactions between each (```trac```, ```rxn```)  
  * Multiple species particle tracking (```ptrk```, ```mptr```, ```sptr```)  
  * Different relative permeability and capillary pressure models (```rlp```, ```exrl```)    
* Equation of state flexibility inherent in code (```eos```)  
* Pseudo-stress and storativity models available
  * Linear porosity deformation (```ppor```)  
  * Gangi stress model (```ppor```)    
* Numerics
  * Finite element with multiple element capabilities (```elem```)  
  * Short form input methods available (```coor```, ```elem```, ```fdm```)  
  * Flexible properties assignment (```zone```, ```zonn```)  
  * Flexible solution methods
    * Upwinding, implicit solution available (```ctrl```)  
    * Iteration control adaptive strategy (```iter```)    
  * Finite volume geometry (```finv```, ```isot```)  
* Flexible time step and stability control (```time```)  
* Time-dependent fixed value and flux boundary conditions (```flow```, ```boun```, ```hflx```)  

## Initialization

The coefficient arrays for the polynomial representations of the density (crl, crv),
enthalpy (cel, cev), and viscosity (cvl, cvv) functions are initialized to the values
enumerated in Appendix of the “Models and Methods Summary” of the
FEHM Application (Zyvoloski et al. 1999), while values for the saturation pressure and
temperature function coefficients are found in Appendix of that document.
All other global array and scalar variables, with the exception of the variables listed in the
table below, are initialized to zero if integer or real, character variables are initialized
to a single blank character, and logical variables are initialized as false.

### Initial (Default) Values

| Variable | Value | Variable | Value | Variable | Value | 
|:---------|:------|:---------|:------|:---------|:------|
| aiaa  | 1.0  | contim  | 1.0e+30  | daymax  | 30.0  |
| daymin  | 1.0e-05  | g1  | 1.0e-06  | g2  | 1.0e-06  |
| g3  | 1.0e-03  | iad_up  | 1000  | iamx  | 500  |
| icons  | 1000  | irlp  | 1  | nbits  | 256  |
| ncntr  | 10000000  | nicg  | 1  | rnmax  | 1.0e+11  |
| str  | 1.0  | strd  | 1.0  | tmch  | 1.0e-09  |
| upwgt  | 1.0  | upwgta  | 1.0  | weight_factor  | 1.0e-3  |

## Restart Procedures

FEHM writes a restart file for each run. The restart output file name may be given in the input control file or as terminal input, or if unspecified will default to fehmn.fin (see [Control File or Terminal I/O Input](input-data)). The file is used on a subsequent run by providing the name of the generated file (via control file or terminal) for the restart input file name. It is recommended that the restart input file name be modified to avoid confusion with the restart output file. For example, by changing the suffix to .ini, the default restart output file, fehmn.fin would be renamed fehmn.ini, and that file name placed in the control file or given as terminal input. Values from the restart file will overwrite any variable initialization prescribed in the input file. The initial time of simulation will also be taken from the restart file unless specified in the macro [time](../Macros/MacroTime.md) input.


## Error Processing

Due to the nonlinearity of the underlying partial differential equations, it is possible to produce an underflow or overflow condition through an unphysical choice of input parameters. More likely the code will fail to converge or will produce results which are out of bounds for the thermodynamic functions. The code will attempt to decrease the time step until convergence occurs. If the time step drops below a prescribed minimum the code will stop, writing a restart file. The user is encouraged to look at the input check file which contains information regarding maximum and minimum values of key variables in the code. All error and warning messages will be output to an output error file or the main output file.

The table below provides additional information on errors that will cause FEHM to terminate.

### Error Conditions Which Result in Program Termination

| Error Condition |   Error Message   | 
|:---------------|:-----------------------|
| I/O file error  | |
| Unable to create / open I/O file | * * * *  Error opening file fileid * * * * ⋅⋅⋅* * * * —————————* * * * * * *  JOB STOPPED * * * * * * —————————* * * *  |
| Coefficient storage file not found | program terminated because coefficient storage file not found |
| Coefficient storage file can not be read | error in parsing beginning of stor file-or-stor file has unrecognized format:quit-or-stor file has neq less than data file:quit  |
| Coefficient storage file already exists  | changing name of new * .stor (old file exists) new file name is fehmn_temp.stor-and-&gt;&gt;&gt; name fehmn_temp.stor is used : stopping  |
| Optional input file not found  | ERROR nonexistant file filenameSTOPPED trying to use optional input file  |
| Unable to open optional input file  | ERROR opening filenameSTOPPED trying to use optional input file  |
| Unable to determine file prefix for AVS output files  | FILE ERROR: nmfil2 file: filename unable to determine contour file prefix  |
| Unable to determine file prefix for pest output files  | FILE ERROR: nmfil15 file: filename unable to determine pest file name-or-FILE ERROR: nmfil16 file: filename unable to determine pest1 file name  |
| Unable to determine file prefix for streamline particle tracking output files  | FILE ERROR: nmfil17 file: filename unable to determine sptr1 file name-or-FILE ERROR: nmfil18 file: filename unable to determine sptr2 file name-or-FILE ERROR: nmfil19 file: filename unable to determine sptr3 file name  |
| Unable to determine file prefix for submodel output file  | FILE ERROR: nmfil24 file: filename unable to determine submodel file name  |
| Input deck errors  | &#160;  |
| Coordinate or element data not found  | * * * *  COOR Required Input * * * * -or-* * * *  ELEM Required Input * * * * * * —————————* * * * * *  JOB STOPPED * * * * * * —————————* * * *   |
| Inconsistent zone coordinates  | inconsistent zone coordinates izone = izone please check icnl in macro CTRL  |
| Invalid AVS keyword read for macro cont  | ERROR:READ_AVS_IOunexpected character string (terminate program execution)Valid options are shown:⋅⋅⋅The invalid string was: string  |
| Invalid keyword or input order read for macro boun  | time change was not first keyword,stop-or-illegal keyword in macro boun, stopping  |
| Invalid keyword read for macro subm  | &gt;&gt;&gt;&gt; error in keyword for macro subm &lt;&lt;&lt;&lt;  |
| Invalid macro read  | * * * *  error in input deck : char * * * *   |
| Invalid parameter values (macros using loop construct)  | Fatal error - for array number arraynummacro - macroGroup number - groupnumSomething other than a real or integer has been specified-or-Line number - lineBad input, check this line-or-Fatal error, too manyreal inputs to initdata2-or-Fatal error, too manyinteger inputs to initdata2  |
| Invalid streamline particle tracking parameter  | ist must be less than or equal to 2  |
| Invalid tracer input  | * *  Using Old InputEnter Temperature Dependency Model Number: 1 - Van Hoff 2 - awwa model, see manual for details * *   |
| Invalid transport conditions  | Fatal error You specified a Henrys Law species with initial concentrations input for the vapor phase (icns = -2), yet the Henrys Constant is computed as 0 for species number speciesnum and node number nodenum. If you want to simulate a vapor-borne species with no interphase transport, then you must specify a gaseous species.  |
| Invalid flag specified for diffusion coefficient calculation  | ERROR – Illegal Flag to concadiffCode Aborted in concadiff  |
| Optional input file name can not be read  | ERROR reading optional input file nameSTOPPED trying to use optional input file  |
| Optional input file contains data for wrong macro  | ERROR –&gt; Macro name in file for macro macroname is wrong_macronameSTOPPED trying to use optional input file  |
| Option not supported  | This option (welbor) not supported.Stop in input-or-  |
| &#160;  | specific storage not available fornon isothermal conditions : stopping-or-gangi model not yet available forair-water-heat conditions : stopping-or-Gencon not yet set for rd1dof.Stop in gencon  |
| Parameter not set  | &gt;&gt;&gt;&gt; gravity not set for head problem: stopping &lt;&lt;&lt;&lt;  |
| Relative permeabilities specified for non-dual or -double porosity model.  | * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * f-m terms but no dpdp : stopping* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   |
| Invalid parameters set  | &#160;  |
| Dual porosity  | * * * *  check fracture volumes,stopping* * * * * * * * *  check equivalent continuum VGs * * * *   |
| Finite difference model (FDM)  | &gt;&gt;&gt;&gt; dimension (icnl) not set to 3 for FDM: stopping &lt;&lt;&lt;&lt;  |
| Maximum number of nodes allowed is less than number of equations  | * * * *  n0(n0) .lt. neq(neq) * * * *  check parameter statements <a href="#id1"><span class="problematic" id="id2">* * </span></a><a href="#id3"><span class="problematic" id="id4">* </span></a>  |
| Node number not in problem domain (macros dvel, flxo, node, nod2, nod3, zone, zonn)  | * * * *  Invalid input: macro macro * * * * ’* * * *  Invalid node specified, value is greater than n0 ( n0 ): stopping <a href="#id5"><span class="problematic" id="id6">* * </span></a><a href="#id7"><span class="problematic" id="id8">* * </span></a>  |
| Noncondensible gas  | cannot input ngas temp in single phase-or-ngas pressure lt 0 at temp and total press givenmax allowable temperature temp-or-ngas pressure gt total pressure i= i-or-ngas pressure lt 0.  |
| Particle tracking  | ERROR: Pcnsk in ptrk must be either always positive or always negative.Code aborted in set_ptrk.f  |
| Relative permeabilities  | cannot have anisotropic perms for rlp model 4 or rlp model 7 with equivalent continuum stopping  |
| Tracer  | ERROR: Can not have both particle tracking (ptrk) and tracer input (trac).Code Aborted in concen.f-or-Gencon not yet set for rd1dof.Stop in gencon-or-ERROR - solute accumulation option cannot be used with cnsk&lt;0-or-* *  On entry to SRNAME parameter number I2 had an illegal value  |
| Insufficient storage  | &#160;  |
| Boundary conditions  | exceeded storage for number of models  |
| Dual porosity  | * * * *  n &gt; n0, stopping * * * *   |
| Generalized dual porosity  | In gdpm macro, ngdpmnodes must be reduced to reduce storage requirementsA value of ngdpm_actual is requiredThe current value set is ngdpmnodes-or-Fatal error in gdpm macroA value of ngdpm_actual is required’The current value set is ngdpmnodesIncrease ngdpmnodes to ngdpm_actual and restart  |
| Geometric coefficients  | program terminated because of insufficient storage  |
| Tracer  | * * * *  memory too small for multiple tracers * * * *   |
| Invalid colloid particle size distribution  | Fatal error, the colloid particle size distribution must end at 1  |
| Invalid particle diffusion  | Fatal errorFor a dpdp simulation, Do not apply the matrix diffusion particle tracking to the matrix nodes, only the fracture nodes  |
| Invalid particle state  | Initial state of particles is invalidParticle number i1  |
| Error computing geometric coefficients  | iteration in zone did not converge, izone = zone_number please check icnl in macro CTRL  |
| Too many negative volumes or finite element coefficients  | too many negative volumes : stopping-or-too many negative coefficients : stopping  |
| Unable to compute local coordinates  | iteration in zone did not converge, izone = zone please check icnl in macro CTRL  |
| Unable to normalize matrix  | cannot normalize  |
| Singular matrix in LU decomposition  | singular matrix in ludcmp  |
| Singular matrix in speciation calculations  | Speciation Jacobian matrix is singular!-or-Scaled Speciation Jacobian matrix is singular!-or-Speciation scaling matrix is singular!  |
| Solution failed to converge  | timestep less than daymin timestep_number current_timestep_size current_simulation_time-or-Tracer Time Step Smaller Than Minimum StepStop in resettrc-or-  |
| &#160;  | Newton-Raphson iteration limit exceeded in speciation subroutine!-or-Newton-Raphson iteration limit exceeded in scaled speciation subroutine!Failure at node i  |

  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
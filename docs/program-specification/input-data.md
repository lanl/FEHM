---
title : Input Data
layout : page_getting-started
permalink: /program-specification/input-data
hero_height: is-hidden
---


# Input Data

## General Considerations

<p><strong>Techniques</strong></p>
<p><strong>Control File or Terminal I/O Startup</strong></p>
<p>The input/output (I/O) file information is provided to the code from an input control file or the terminal. The default control file name is <code class="docutils literal notranslate"><span class="pre">fehmn.files</span></code>. If a control file with the default name is present in the directory from which the code is being executed, no terminal input is required. If the default control file is not present, it is possible to supply the name of the control file on the command line, otherwise input prompts are written to the screen preceded by a short description of the I/O files used by FEHM. It should be noted that a control file name entered on the command line will take precedence over the default control file. The descriptions of the I/O files are elaborated on in <a href="/getting-started/data-files"> Data Files </a>. The initial prompt asks for the name of a control file. It is also If a control file name is entered for that prompt no further terminal input is required. If a control file is not used, the user is then prompted for I/O file names, the tty output flag, and user subroutine number. When the input file name is entered from the terminal the user has the option of letting the code generate the names for the remainder of the auxiliary files using the input file name prefix. The form of the input file name is <code class="docutils literal notranslate"><span class="pre">filen</span></code> or <code class="docutils literal notranslate"><span class="pre">filen.*</span></code> where <code class="docutils literal notranslate"><span class="pre">filen</span></code> is the prefix used by the code to name the auxiliary files and <code class="docutils literal notranslate"><span class="pre">.*</span></code> represents an arbitrary file extension.</p>
<p><strong>Multiple Realization Simulations</strong></p>
<p>The code has an option for performing multiple simulation realizations (calculations) where input (e.g., porosity, permeability, saturation, transport properties or particle distributions) is modified for each realization but the calculations are based on the same geometric model. Multiple realizations are initiated by including a file called fehmn.msim in the directory from which the code is being run. If invoked, a set number of simulations are performed sequentially, with pre- and post-processing steps carried out before and after each simulation. This capability allows multiple simulations to be performed in a streamlined fashion, with processing to change input files before each run and post-processing to obtain relevant results after each run.</p>
<p><strong>Macro Control Structure</strong></p>
<p>The finite element heat and mass transfer code (FEHM) contains a macro control structure for data input that offers added flexibility to the input process. The macro command structure makes use of a set of control statements recognized by the input module of the program. When a macro control statement is encountered in an input file, a certain set of data with a prescribed format is expected and read from the input file. In this way, the input is divided into separate, unordered blocks of data. The input file is therefore a collection of macro control statements, each followed by its associated data block. Blocks of data can be entered in any order, and any blocks unnecessary to a particular problem need not be entered. The macro control statements must appear in the first four columns of a line. The other entries are free format, which adds flexibility, but requires that values be entered for all input variables (no assumed null values).</p>
<p>As an aid to the user, the capabilities of FEHM summarized in <a class="reference external" href="ProgramConsiderations.md">Capabilities of FEHM with Macro Command References</a> refer to applicable macro commands. <a class="reference external" href="../Macros.md">Macro Control Statements for FEHM</a> lists the macro control statements with a brief description of the data associated with each. A more detailed description of each macro control statement and its associated input are found in their respective pages. Macro control statements may be called more than once, if, for example, the user wishes to reset some property values after defining alternate zones. Some statements are required, as indicated in <a class="reference external" href="../Macros.md">Macro Control Statements for FEHM</a>, the others are optional.</p>

## Macro Control Statements for FEHM (Alphabetical List)

| Macro |  Description    | 
|:---------------|:-----------------------|
| [adif](../Macros/MacroAdif.md)| Air-water vapor diffusion|
| [airwater or air](../Macros/MacroAirwater.md)| Isothermal air-water input|
| [anpe](../Macros/MacroAnpe.md)| Anisotropic permeability|
| [boun](../Macros/MacroBoun.md)| Boundary conditions (required for flow problem if macro flow is not used)|
| [bous](../Macros/MacroBous.md)| Boussinesq-type approximation|
| [carb](../Macros/MacroCarb.md)| CO2 input|
| [cden](../Macros/MacroCden.md)| Concentration-dependent density|
| [cflx](../Macros/MacroCflx.md)| Molar flow rate through a zone|
| [cgdp](../Macros/MacroCgdp.md)| Rate-limited gdpm node|
| [chea](../Macros/MacroChea.md)| Output in terms of head, not pressures (non-head problem)|
| [cond](../Macros/MacroCond.md)| Thermal conductivity data (required for non-isothermal problem)|
| [conn](../Macros/MacroConn.md)| Print number of connections for each node and stop|
| [cont](../Macros/MacroCont.md)| Contour plot data|
| [conv](../Macros/MacroConv.md)| Head input conversion for thermal problems|
| [coor](../Macros/MacroCoor.md)| Node coordinate data (required if macro fdm is not used)|
| [cont](../Macros/MacroCont.md)| Program control parameters (required)|
| [dpdp](../Macros/MacroDpdp.md)| Double porosity/double permeability model input|
| [dual](../Macros/MacroDual.md)| Dual porosity model input|
| [dvel](../Macros/MacroDvel.md)| Velocity printout (formerly macro velo)|
| [elem](../Macros/MacroElem.md)| Element node data (required if macro fdm is not used)|
| [eos](../Macros/MacroEos.md)| Simple equation of state data|
| [evap](../Macros/MacroEvap.md)| Evaporation model|
| [exrl](../Macros/MacroExrl.md)| Explicit evaluation of relative permeability|
| [fdm](../Macros/MacroFdm.md)| Finite difference grid generation (required if macro coor and elem are not used)|
| [finv](../Macros/MacroFinv.md)| Finite volume flow coefficients|
| [flgh](../Macros/MacroFlgh.md)| Generalized head boundary conditions (confined aquifer)|
| [flow](../Macros/MacroFlow.md)| Flow data (required for flow problem if macro boun is not used)|
| [flo2](../Macros/MacroFlo2.md)| Alternate format for flow data (input using 3-D planes)|
| [flo3](../Macros/MacroFlo3.md)| Alternate format for flow data (defined for see page faces)|
| [floa](../Macros/MacroFloa.md)| Alternate format for flow data (additive to previous flow definition)|
| [flwt](../Macros/MacroFlwt.md)| Movable source or sink (wtsi only)|
| [flxn](../Macros/MacroFlxn.md)| Write all non-zero source/sink internodal mass flows by node to an output file.|
| [flxo](../Macros/MacroFlxo.md)| Internodal mass flow printout|
| [flxz](../Macros/MacroFlxz.md)| Zone based mass flow output|
| [fper](../Macros/MacroFper.md)| Permeability scaling factor|
| [frlp](../Macros/MacroFrlp.md)| Relative permeability factors for residual air effect|
| [ftsc](../Macros/MacroFtsc.md)| Flux correction for saturations over 1|
| [gdkm](../Macros/MacroGdkm.md)| Generalized dual permeability model|
| [gdpm](../Macros/MacroGdpm.md)| Generalized dual porosity model|
| [grad](../Macros/MacroGrad.md)| Gradient model input|
| [hcon](../Macros/MacroHcon.md)| Set solution to heat conduction only|
| [head](../Macros/MacroHead.md)| Hydraulic head input|
| [hflx](../Macros/MacroHflx.md)| Heat flow input|
| [hist](../Macros/MacroHist.md)| User selected history output|
| [hyco](../Macros/MacroHyco.md)| Hydraulic conductivity input (required if macro perm is not used)|
| [ice](../Macros/MacroIce.md) or [meth](../Macros/MacroIce.md)| Ice phase calculations, methane hydrate input|
| [imex](../Macros/MacroImex.md)| implicit-explicit solution|
| [impf](../Macros/MacroImpf.md)| Time step control based on maximum allowed variable change|
| [init](../Macros/MacroInit.md)| Initial value data (required if macro pres or restart file is not used)|
| [intg](../Macros/MacroIntg.md)| Set integration type for finite element coefs|
| [isot](../Macros/MacroIsot.md)| Isotropic definition of control volume/finite element coefficients|
| [iter](../Macros/MacroIter.md)| Iteration parameters|
| [itfc](../Macros/MacroItfc.md)| Flow and transport between zone interfaces|
| [ittm](../Macros/MacroIttm.md)| Sticking time for phase changes|
| [itup](../Macros/MacroItup.md)| Iterations used with upwinding|
| [iupk](../Macros/MacroIupk.md)| Upwind transmissibility including intrinsic permeability|
| [ivfc](../Macros/MacroIvfc.md)| Enable exponential fracture and volume model|
| [mdnode](../Macros/MacroMdnode.md)| Enables extra connections to be made to nodes|
| [meth](../Macros/MacroIce.md) or [ice](../Macros/MacroIce.md)| Methane hydrate input|
| [mptr](../Macros/MacroMptr.md)| Multiple species particle tracking simulation input|
| [nfinv](../Macros/MacroNfinv.md)| Finite element instead of finite volume calculations|
| [ngas](../Macros/MacroNgas.md)| Noncondensible gas (air) data|
| [nobr](../Macros/MacroNobr.md)| Don’t break connection between nodes with boundary conditions|
| [node](../Macros/MacroNode.md)| Node numbers for output and time histories|
| [nod2](../Macros/MacroNod2.md)| Node numbers for output and time histories, and alternate nodes for terminal output|
| [nod3](../Macros/MacroNod3.md)| Node numbers for output and time histories, alternate nodes for terminal output, and alternate nodes for variable porosity model information|
| [nrst](../Macros/MacroNrst.md)| Stop NR iterations on variable changes|
| [para](../Macros/MacroPara.md)| Parallel FEHM (isothermal only)|
| [perm](../Macros/MacroPerm.md)| Permeability input (required if macro hyco is not used)|
| [pest](../Macros/MacroPest.md)| Parameter estimation routine output|
| [phys](../Macros/MacroPhys.md)| Non-darcy well flow|
| [ppor](../Macros/MacroPpor.md)| Pressure and temperature dependent porosity and permeability|
| [pres](../Macros/MacroPres.md)| Initial pressure, temperature, and saturation data, boundary conditions specification (required if macro init or restart file is not used)|
| [ptrk](../Macros/MacroPtrk.md)| Particle tracking simulation input|
| [renu](../Macros/MacroRenu.md)| Renumbers nodes|
| [rest](../Macros/MacroRest.md)| Manage restart options |
| [rflo](../Macros/MacroRflo.md)| Read in flux values|
| [rich](../Macros/MacroRich.md)| Enable Richards’ equation|
| [rive](../Macros/MacroRive.md) or [well](../Macros/MacroRive.md)| River or implicit well package|
| [rlp](../Macros/MacroRlp.md)| Relative permeability input (required for 2-phase problem if macro rlpm is not used, otherwise optional)|
| [rlpm](../Macros/MacroRlpm.md)| Alternate style relative permeability input (required for 2-phase problem if macro rlp is not used, otherwise optional)|
| [rock](../Macros/MacroRock.md)| Rock density, specific heat, and porosity input (required)|
| [rxn](../Macros/MacroRxn.md)| Chemical reaction rate model input|
| [sol](../Macros/MacroSol.md)| Solver specifications|
| [sptr](../Macros/MacroSptr.md)| Streamline particle tracking simulation input|
| [stea](../Macros/MacroStea.md)| Steady state program termination|
| [stop](../Macros/MacroStop.md)| Signals the end of input (required)|
| [strs](../Macros/MacroStrs.md)| Stress solution enabled|
| [subm](../Macros/MacroSubm.md)| Submodel boundary condition output|
| [svar](../Macros/MacroSvar.md)| Enable pressure-enthalpy variables|
| [text](../Macros/MacroText.md)| Text input to be written to output file|
| [thic](../Macros/MacroThic.md)| Variable thickness input for two-dimensional problems|
| [time](../Macros/MacroTime.md)| Time step and time of simulation data (required)|
| [trac](../Macros/MacroTrac.md)| Solute simulation input|
| [trxn](../Macros/MacroTrxn.md)| A user-friendly replacement for [[wiki:MacroTra|
| [user](../Macros/MacroUser.md)| User subroutine call|
| [vapl](../Macros/MacroVapl.md)| Vapor pressure lowering|
| [vbou](../Macros/MacroVbou.md)| read non-pyhsical gridblock volumes|
| [vcon](../Macros/MacroVcon.md)| Variable thermal conductivity input|
| [weli](../Macros/MacroWeli.md)| Peaceman type well impedance|
| [well](../Macros/MacroRive.md) or [rive](../Macros/MacroRive.md)| Implicit well package or River|
| [wgtu](../Macros/MacroWgtu.md)| Areas, weights (user-defined) for boundary conditions|
| [wflo](../Macros/MacroWflo.md)| Alternate submodel boundary output|
| [wtsi](../Macros/MacroWtsi.md)| Water table, simplified|
| [zneg](../Macros/MacroZneg.md)| Problem solving tool removes negative coefs|
| [zone](../Macros/MacroZone.md)| Geometric definition of grid for input parameter assignment|
| [zonn](../Macros/MacroZonn.md)| Same as zone, except zonn definitions are not overwritten|




| Old Macro |  Original Description, current behavior is undefined.    | 
|:---------------|:-----------------------|
| [alti](../Macros/MacroAlti.md) | Alternate element and coordinate input, grid format out of date.|
| dof|no longer used|
| exuz|no longer used|
| naple or [szna](../Macros/MacroSzna.md)|Isothermal water input, enable simplified napl solution|
| [rflx](../Macros/MacroRflx.md)|Radiation source term|
| solv|no longer used|
| wlbr|no longer used, replaced by [[wiki:MacroRiv|
| [zeol](../Macros/MacroZeol.md)|Zeolite water balance input|


<p>Comments may be entered in the input file by beginning a line with a <code class="docutils literal notranslate"><span class="pre">#</span></code> symbol (the <code class="docutils literal notranslate"><span class="pre">#</span></code> symbol must be found in the first column of the line). Comments may precede or follow macro blocks but may not be found within a block.</p>

<p>Optional input files may be used by substituting a keyword and file name in the main input file (described in detail in <a class="reference internal" href="#optional-input-files">Optional Input Files</a>). The normal macro input is then entered in the auxiliary file.</p>

<p>A macro may be disabled (turned off or omitted from a run) by adding keyword “off” on the macro line and terminating the macro with an end statement of the form <code class="docutils literal notranslate"><span class="pre">end//macro//</span></code> or <code class="docutils literal notranslate"><span class="pre">end</span> <span class="pre">//macro//</span></code> (see <a class="reference internal" href="#option-to-disable-a-macro">Option to disable a macro</a>).</p>

## Input Parameters
<p>Many input parameters such as porosity or permeability vary throughout the grid and need to have different values assigned at different nodes. This is accomplished in two ways. The first uses a nodal loop-type definition (which is the default):</p>

<blockquote>
<div>JA, JB, JC, PROP1, PROP2, …</div></blockquote>
<p>where</p>
<blockquote>
<div><p>JA - first node to be assigned with the properties PROP1, PROP2, …</p>
<p>JB - last node to be assigned with the properties PROP1, PROP2, …</p>
<p>JC -  loop increment for assigning properties PROP1, PROP2, ….</p>
<p>PROP1, PROP2, etc. - property values to be assigned to the indicated nodes.</p>
</div></blockquote>

<p>In the input blocks using this structure, one or more properties are manually entered in the above structure. When a blank line is entered, that input block is terminated and the code proceeds to the next group or control statement. (Note that blank input lines are shaded in the examples shown in <a class="reference internal" href="#individual-input-records-or-parameters">Individual Input Records or Parameters</a>.) The nodal definition above is useful in simple geometries where the node numbers are easily found. Boundary nodes often come at regular node intervals and the increment counter JC can be adjusted so the boundary conditions are easily entered. To set the same property values at every node, the user may set JA and JC to 1 and JB to the total number of nodes, or alternatively set <code class="docutils literal notranslate"><span class="pre">JA</span> <span class="pre">=</span> <span class="pre">1</span></code>, and <code class="docutils literal notranslate"><span class="pre">JB</span> <span class="pre">=</span> <span class="pre">JC</span> <span class="pre">=</span> <span class="pre">0</span></code>.</p>

<p>For dual porosity problems, which have three sets of parameter values at any nodal position, nodes 1 to N [where N is the total number of nodes in the grid (see macro <a class="reference external" href="MacroCoor.md">coor</a>)] represent the fracture nodes, nodes N + 1 to 2N are generated for the second set of nodes, the first matrix material, and nodes 2N + 1 to 3N for the third set of nodes, the second matrix material. For double porosity/double permeability problems, which have two sets of parameter values at any nodal position, nodes 1 to N represent the fracture nodes and nodes N + 1 to 2N are generated for the matrix material.</p>

<p>For more complicated geometries, such as 3-D grids, the node numbers are often difficult to determine. Here a geometric description is preferred. To enable the geometric description the <a class="reference external" href="MacroZone.md">zone</a> control statement <a class="reference external" href="MacroWgtu.md">wgtu</a> is used in the input file before the other property macro statements occur. The input macro <strong>zone</strong> requires the specification of the coordinates of 4-node parallelograms for 2-D problems or 8-node polyhedrons in 3-D. In one usage of the control statement <strong>zone</strong> all the nodes are placed in geometric zones and assigned an identifying number. This number is then addressed in the property input macro commands by specifying a JA &lt; 0 in the definition of the loop parameters given above. For example if JA = -1, the properties defined on the input line would be assigned to the nodes defined as belonging to geometric Zone 1 (JB and JC must be input but are ignored in this case). The control statement <strong>zone</strong> may be called multiple times to redefine geometric groupings for subsequent input. The previous zone definitions are not retained between calls. Up to 1000 zones may be defined. For dual porosity problems, which have three sets of parameter values at any nodal position, Zone <code class="docutils literal notranslate"><span class="pre">ZONE_DPADD</span> <span class="pre">+</span> <span class="pre">I</span></code> is the default zone number for the second set of nodes defined by Zone I, and Zone <code class="docutils literal notranslate"><span class="pre">2*ZONE_DPADD</span> <span class="pre">+</span> <span class="pre">I</span></code> is the default zone number for the third set of nodes defined by Zone I. For double porosity/double permeability problems, which have two sets of parameter values at any nodal position, Zone <code class="docutils literal notranslate"><span class="pre">ZONE_DPADD</span> <span class="pre">+</span> <span class="pre">I</span></code> is the default zone number for the second set of nodes defined by Zone I. The value of <code class="docutils literal notranslate"><span class="pre">ZONE_DPADD</span></code> is determined by the number of zones that have been defined for the problem. If less than 100 zones have been used <code class="docutils literal notranslate"><span class="pre">ZONE_DPADD</span></code> is set to 100, otherwise it is set to 1000. Zones of matrix nodes may also be defined independently if desired.</p>

<p>Alternatively, the <a class="reference external" href="MacroZonn.md">zonn</a> control statement may be used for geometric descriptions. Regions are defined the same as for control statement zone except that previous zone definitions are retained between calls unless specifically overwritten.</p>

<p><strong>Interface</strong></p>
<p>To interface with <a class="reference external" href="http://en.wikipedia.org/wiki/GoldSim">GoldSim</a>, FEHM is compiled as a dynamic link library (DLL) subroutine that is called by the GoldSim code. When FEHM is called as a subroutine from GoldSim, the GoldSim software controls the time step of the simulation, and during each call, the transport step is carried out and the results passed back to GoldSim for processing and/or use as radionuclide mass input to another portion of the GoldSim system, such as a saturated zone transport submodel. The interface version of FEHM is set up only to perform particle tracking simulations of radionuclide transport, and is not intended to provide a comprehensive flow and transport simulation capability for GoldSim. Information concerning the GoldSim user interface may be found in the GoldSim documentation (Golder Associates, 2002).</p>

<p><strong>Consecutive Cases</strong></p>
<p>Consecutive cases can be run using the multiple realizations simulation option (see <a class="reference internal" href="#multiple-realization-simulations">Multiple Realization Simulations</a>). The program retains only the geometric information between runs (i.e., the grid and coefficient information). The values of all other variables are reinitialized with each run, either from the input files or a restart file when used.</p>

<p><strong>Defaults</strong></p>

<p>Default values are set during the initialization process if overriding input is not provided by the user.</p>

## Individual Input Records or Parameters

<p>Other than the information provided through the control file or terminal I/O and the multiple realization simulations file, the main user input is provided using macro control statements in the input file, geometry data file, zone data file, and optional input files. Data provided in the input files is entered in free format with the exception of the macro control statements and keywords which must appear in the first four (or more) columns of a line. Data values may be separated with spaces, commas, or tabs. The primary input file differs from the others in that it begins with a title line (80 characters maximum) followed by input in the form of the macro commands. Each file containing multiple macro commands should be terminated with the <a class="reference external" href="MacroStop.md">stop</a> control statement. In the examples provided in the following subsections, blank input lines are depicted with shading.</p>

<p><strong>Control File or Terminal I/O Input</strong></p>

<p>The file name parameters enumerated below <code class="docutils literal notranslate"><span class="pre">[nmfil(2-13)]</span></code>, are entered in order one per line in the control file (excluding the control file name <code class="docutils literal notranslate"><span class="pre">[nmfil(1)]</span></code> and error file name <code class="docutils literal notranslate"><span class="pre">[nmfil(14)]</span></code>) or in response to a prompt during terminal input. If there is a control file with the name fehmn.files in your local space (current working directory), FEHM will execute using that control file and there will be no prompts. If another name is used for the control file, it can be entered on the command line or at the first prompt.</p>

<p>A blank line can be entered in the control file for any auxiliary files not required, for the “none” option for tty output, and for the “0” option for the user subroutine number.</p>

<p>In version 2.30 an alternate format for the control file has been introduced that uses keywords to identify which input and output files will be used. Please note that the file name input styles may not be mixed.</p>

<blockquote>
<div>Group 1: <code class="docutils literal notranslate"><span class="pre">NMFIL(i)</span></code> (a file name for each i = 2 to 13)</div></blockquote>

<p>or</p>

<blockquote>
<div><p>Group 1: <code class="docutils literal notranslate"><span class="pre">KEYWORD:</span> <span class="pre">NMFIL</span></code></p>
<p>Group 2: <code class="docutils literal notranslate"><span class="pre">TTY_FLAG</span></code></p>
<p>Group 3: <code class="docutils literal notranslate"><span class="pre">USUB_NUM</span></code></p>
</div></blockquote>

<p>Unlike previous versions of the code, if a file name is not entered for the output file, check file, or restart file, the file will not be generated. An error output file will still be generated for all runs (default name <code class="docutils literal notranslate"><span class="pre">fehmn.err</span></code>). However, with the keyword input style the user has the option of naming the error file. File names that do not include a directory or subdirectory name, will be located in the current working directory. With keyword input a root filename may be entered for output files that use file name generation (hist macro output, cont macro avs, surfer or tecplot output, etc.). The data files are described in more detail in <a class="reference external" href="DataFiles.md">Data Files</a>.</p>



| Input Variable |  Format    | Default | Description|
|:---------------|:-----------|:--------|:-----------------------|
| keyword| character*5 | No keywords, old style file format is used | Keyword specifying input or output file type. The keyword is entered followed immediately by a “:” with a “space” preceding the filename. Keywords, which must be entered starting in the first column, are:|
| | | | input    -   Main input filegrid -   Geometry data file; or grida or gridf - Ascii formatted geometry file; or gridu or gridb - Unformatted geometry filezone    -   Initial zone fileoutp   -   Output filersti -   Restart input filersto  -   Restart output filehist -   Simulation history outputtrac   -   Solute history outputcont   -   Contour outputdual, dpdp - Dual porosity, double porosity outputstor    -   Coefficient storage filecheck   –  Input check output filenopf -   Symbolic factorization filecolu -   Column data file for free surface problemserror -   Error output fileroot   -   Root name for output file name generationco2i   -   CO,,2,, parameter data file (default co2_interp_table.txt)look  -   Equation of state data lookup table file (default lookup.in)|
| | | | Keyword file name input is terminated with a blank line. The keywords and file names may be entered in any order unlike the old style input.|
| nmfil| character*100| | Input or output file name|
| | | fehmn.files | nmfil( 1) - Control file name (this file is not included in the old style control file) (optional)|
| | | fehmn.dat | nmfil( 2) - Main input file name (required)|
| | | not used |nmfil( 3) - Geometry data input file name (optional)|
| | | not used | nmfil( 4) - Zone data input file name (optional)|
| | | not used | nmfil( 5) - Main output file name (optional)|
| | | not used | nmfil( 6) - Restart input file name (optional)|
| | | not used | nmfil( 7) - Restart output file name (optional)|
| | | not used | nmfil( 8) - Simulation history output file name (optional)|
| | | not used | nmfil( 9) - Solute history output file name (optional)|
| | | not used | nmfil(10) - Contour plot output file name (optional)|
| | | not used | nmfil(11) - Dual porosity or double porosity / double permeability contour plot output file name (optional)|
| | | not used | nmfil(12) - Coefficient storage file name (optional)|
| | | not used | nmfil(13) - Input check output file name (optional)|
| | | fehmn.err | nmfil(14) - Error output file name (this file is not included in the old style control file). The default name is used if not input.|
|tty_flag | character*4 | none | Terminal output flag: all, some, none|
| usub_num | integer | 0 | User subroutine call number|

<p>The following are examples of the input control file. The first example (left) uses keyword style input, while the second and third examples (right) use the original style control file input form. In the first example, four files are explicitly named, the input file, geometry file, tracer history output file and output error file. A root file name is also provided for file name generation. The “all” keyword indicates that all information should be written to the terminal and the ending “0” indicates that the user subroutine will not be called. In the second example in the center, all input will be found in the current working directory and output files will also be written to that directory. The blank lines indicate that there is no restart initialization file or restart output file, a dual porosity contour plot file is not required, and the coefficient storage file is not used. The “some” keyword indicates that selected information is output to the terminal. The ending “0” indicates that the user subroutine will not be called. In the third example on the right, input will be found in the “groupdir” directory, while output will be written to the current working directory. The “none” keyword indicates that no information should be written to the terminal and the ending “0” indicates that the user subroutine will not be called.</p>
<p>Files <code class="docutils literal notranslate"><span class="pre">fehmn.files</span></code>:</p>

| | | | | |
|:---------------|:-----------|:--------|:----------|:------------|
|input: /groupdir/c14-3| | tape5.dat | | /groupdir/c14-3|
| trac: c14-3.trc | | tape5.dat | | /groupdir/grid-402|
| grid: /groupdir/grid-402 | | tape5.dat | | /groupdir/c14-3|
| root: c14-3 | | tape5.out | | c14-3.out|
| error: c14-3.err | |  | | /groupdir/c14-3.ini|
|  | |  | | c14-3.fin|
| all | | tape5.his | | c14-3.his|
| 0 | | tape5.trc | | c14-3.trc|
|  | | tape 5.con | | c14-3.con|
|   | |   | | c14-3.dp|
|   | |   | | c14-3.stor|
|   | | tape5.chk | | c14-3.chk|
|   | | some | | none|
|   | | 0 | | 0|

## Multiple Realization Simulations

<p>The multiple realization simulations input file (fehmn.msim) contains the number of simulations to be performed and, on UNIX systems, instructions for pre- and post-processing input and output data during a multiple realization simulation. The file uses the following input format:</p>

<p>Line 1  nsim</p>

<p>Lines 2-N single_line</p>

| Input Variable | Format | Description |
|:---------------|:-----------|:--------|
| nsim | integer | Number of simulation realizations to be performed|
| single_line | character*80 | An arbitrary number of lines of UNIX shell script or Windows bat file instructions: lines 2-n:lines which are written to a file called fehmn.pre (UNIX) or fehmn.pre.bat (Windows), which is invoked before each realization using the following command: sh fehmn.pre $1 $2 (UNIX systems) or fehmn.pre.bat $1 $2 (Windows)|
|   |   | line n+1: the keyword ‘post’, placed in the first four columns of the input file, denotes that the previous line is the last line in the fehmn.pre script, and that the data for the post-processing script fehmn.post follows|
|   |   | lines n+2 to N: lines which are written to a file called fehmn.post (UNIX) or fehmn.post.bat (Windows), which is invoked after each realization using the following command: sh fehmn.post $1 $2 (UNIX) or fehmn.post.bat $1 $2 (Windows)|
|   |   | Thus, the files fehmn.pre and fehmn.post are created by the code and are meant to provide the capability to perform complex pre- and post-processing steps during a multiple realization simulation. Script arguments $1 and $2 represent the current simulation number and nsim, the total number of simulations, respectively.|


<p>In the following (UNIX style) example, 100 simulations are performed with pre and post-processing steps carried out. File “fehmn.msim” contains the following:</p>

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span>echo This is run number $1 of $2
rm fehmn.filescurnum=`expr $1`
cp control.run fehmn.filesrm ptrk.inputcp ptrk.np$curnum ptrk.input
echo starting up the run
post
curnum=`expr $1`/home/robinson/fehm/chun_li/ptrk/process1_fujrm
np$curnum.outputmv results.output np$curnum.output
echo finishing the run again
</pre></div></div>

<p>The first line after the number of simulations demonstrates how the current and total number of simulations can be accessed in the fehmn.pre shell script. This line will write the following output for the first realization:</p>

<p>This is run number 1 of 100</p>

<p>The pre-processing steps in this example are to remove the fehmn.files file from the working directory, copy a control file to fehmn.files, copy a particle tracking macro input file to a commonly named file called ptrk.input, and write a message to the screen. The fehmn.files files should be used or else the code will require screen input of the control file name for every realization. One hundred particle tracking input files would have been generated previously, and would have the names ptrk.np1, ptrk.np2, …, ptrk.np100. Presumably, these files would all have different transport parameters, resulting in 100 different transport realizations. The post-processing steps involve executing a post-processor program for the results (process1_fuj). This post-processor code generates an output file called results.output, which the script changes to np1.output, np2.output, …, np100.output, for further processing after the simulation.</p>

<p>One other point to note is that the variable “curnum” in this example is defined twice, once in the pre-processor and again in the post-processor. This is necessary because fehmn.pre and fehmn.post are distinct shell scripts that do not communicate with one another.</p>

## Transfer function curve data input file

<p>In the FEHM particle tracking models, diffusion of solute from primary porosity into the porous matrix is captured using an upscaling procedure designed to capture the small scale behavior within the field scale model. The method is to impart a delay to each particle relative to that of a nondiffusing solute. Each particle delay is computed probabilistically through the use of transfer functions. A transfer function represents the solution to the transport of an idealized system with matrix diffusion at the sub-grid-block scale. After setting up the idealized model geometry of the matrix diffusion system, a model curve for the cumulative distribution of travel times through the small-scale model is computed, either from an analytical or numerical solution. Then, this probability distribution is used to determine, for each particle passing through a given large-scale grid block, the travel time of a given particle. Sampling from the distribution computed from the small scale model ensures that when a large number of particles pass through a cell, the desired distribution of travel times through the model is reproduced. In FEHM, there are equivalent continuum and dual permeability formulations for the model, each of which call for a different set of sub-grid-block transfer function curves. These curves are numerical input to the FEHM, with a data structure described below. Optional input in macros mptr, ptrk, and sptr is used to tell the code when transfer function curves are required and whether 2 or 3 (numparams) parameter curves are to be used.</p>

| Input Variable | Format | Description |
|:---------------|:-----------|:--------|
| DUMMY_STRING | character*4 | If keyword “free” is input at the beginning of the file, the code assumes an irregular distribution of transfer function parameters and performs either a nearest neighbor search to choose the transfer function curve, or a more computationally intensive multi-curve interpolation.|
| NCURVES | integer | Total number of transfer function curves input in the file when keyword “free” is used.|
| LOG_FLAG | integer | Array of flags to denote whether to take the natural log of the parameter before determining the distance from the input parameter to the values for each transfer function when keyword “free” is used. If 0, use the input value, if 1, use ln(value).|
| NORMAL_PARAM | real | Array of values for normalizing the parameters before performing the nearest neighbor or multi-curve interpolation. Each value is divided by NORMAL_PARAM(i) if LOG_FLAG is 0, and ln(LOG_FLAG) if LOG_FLAG is 1.|
| CURVE_STRUCTURE | integer | Flag to denote the type of interpolation to be performed when keyword “free” is used. If 1, simple nearest neighbor search, if &gt; 1, a multi-curve interpolation is performed with curve_strucuture points nearest to the input values of the parameters. It is recommended that values no greater than 4 be used.|
| WEIGHT_FACTOR | real | Optional wieght factor used when CURVE_STRUCTURE &gt; 1. The default value is 1e-3. When determining the interpolated value of time using a multi-curve interpolation there are occasions where the algorithm yields large values of the weights used to compute the particle residence time. In a few such cases numerical errors can make the scheme fail so that the interpolated values for time erroneously get very large. This occurs when the sum of the weights divided by any individual weight is small, that is, large weights of opposite sign cancelling one another out. To prevent this error in the scheme from affecting the results, the code reverts to a nearest neighbor approach to obtain the time. The criterion for this option is that the sum of the weights divided by any individual weight is less than weight_factor. Increasing this value to 1.e-2 or higher can eliminate such occurrences. This parameter is very problem dependent, so this parameter is included for flexibility. It is recommended that the default of 1.e-3 or a higher value of 1.e-2 or so be used.|
| NUMP1 | integer | Number of parameter 1 values used to define transfer function curves.|
| PARAM1 | real | nump1 parameter 1 values defining transfer function curves.|
| NUMP2 | integer | Number of parameter 2 values used to define transfer function curves.|
| PARAM2 | real | nump2 parameter 2 values defining transfer function curves.|
| NUMP3 | integer | Number of parameter 3 values used to define transfer function curves.|
| PARAM3 | real | nump3 parameter 3 values defining transfer function curves.|
| D4 | integer | Fracture-matrix flow interaction flag (d4 = 1, 4). For the three-parameter option, the dual permeability model requires four transfer function curves for each set of parameters. Interactions can occur from fracture-fracture (d4=1), fracture-matrix (d4=2), matrix-fracture (d4=3), and matrix-matrix (d4=4).|
| NUMP_MAX | integer | Maximum number of delay time and concentration values for transfer function curves.|
| NUMP | integer | Number of delay time and concentration values in each transfer function curve (nump1, nump2, nump3, d4).|
| DTIME | real | Transfer function curve delay times (nump1, nump2, nump3, d4, nump).|
| CONC | real | Transfer function curve concentrations (nump1, nump2, nump3, d4, nump).|
| OUTPUT_FLAG | character*3 | If optional keyword “out” is entered at the end of the file the code outputs information on the parameter space encountered during the simulation in the <a href="#id9"><span class="problematic" id="id10">*</span></a>.out file. See [[wiki:MacroMptr|

<p>The transfer function curve data file uses the following format if a regular grid of parameters is input. Please note that parameter values for this format should be input from smallest to largest value:</p>

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">nump1</span>

<span class="n">param1</span> <span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="n">nump1</span>

<span class="n">nump2</span>

<span class="n">param2</span> <span class="p">(</span><span class="n">j</span><span class="p">),</span> <span class="n">j</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="n">nump2</span>
</pre></div>
</div>
<p>If 2 parameter curves are being input</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">nump_max</span>

<span class="n">For</span> <span class="n">each</span> <span class="n">i</span><span class="p">,</span> <span class="n">j</span> <span class="p">(</span><span class="n">nump3</span> <span class="o">=</span><span class="mi">1</span><span class="p">,</span> <span class="n">d4</span> <span class="o">=</span> <span class="mi">1</span><span class="p">)</span>

<span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">),</span> <span class="n">param1</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">param2</span><span class="p">(</span><span class="n">j</span><span class="p">)</span>

<span class="n">followed</span> <span class="n">by</span> <span class="k">for</span> <span class="n">each</span> <span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">)</span>

<span class="n">time</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nump</span><span class="p">),</span> <span class="n">conc</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nump</span><span class="p">)</span>
</pre></div>
</div>
<p>Or if 3 parameter curves are being input</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">nump3</span>

<span class="n">param3</span><span class="p">(</span><span class="n">k</span><span class="p">),</span> <span class="n">k</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="n">nump3</span>

<span class="n">nump_max</span>

<span class="n">For</span> <span class="n">each</span> <span class="n">d4</span><span class="p">,</span> <span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="n">k</span>

<span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="n">k</span><span class="p">,</span> <span class="n">d4</span><span class="p">),</span> <span class="n">param1</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">param2</span><span class="p">(</span><span class="n">j</span><span class="p">),</span> <span class="n">param3</span><span class="p">(</span><span class="n">k</span><span class="p">),</span> <span class="n">d4</span>

<span class="n">followed</span> <span class="n">by</span> <span class="k">for</span> <span class="n">each</span> <span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="n">k</span><span class="p">,</span> <span class="n">d4</span><span class="p">)</span>

<span class="n">time</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="n">k</span><span class="p">,</span> <span class="n">d4</span><span class="p">,</span> <span class="n">nump</span><span class="p">),</span> <span class="n">conc</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">,</span> <span class="n">k</span><span class="p">,</span> <span class="n">d4</span><span class="p">,</span> <span class="n">nump</span><span class="p">)</span>

<span class="n">out_flag</span> <span class="p">(</span><span class="n">optional</span><span class="p">)</span> <span class="o">-</span> <span class="n">keyword</span> <span class="s2">&quot;out&quot;</span>
</pre></div>
</div>
<p>The transfer function curve data file uses the following format for the case in which the transfer functions are input without a regular grid of parameters:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">dummy_flag</span> <span class="o">-</span> <span class="n">keyword</span> <span class="s2">&quot;free&quot;</span>

<span class="n">log_flag</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="n">numparams</span>

<span class="n">normal_param</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="n">numparams</span>

<span class="n">curve_structure</span><span class="p">,</span> <span class="n">weight_factor</span> <span class="p">(</span><span class="n">optional</span><span class="p">)</span>

<span class="n">ncurves</span>

<span class="n">nump_max</span>
</pre></div>
</div>
<p>For “free” form input of transfer function curves <code class="docutils literal notranslate"><span class="pre">(nump1</span> <span class="pre">=</span> <span class="pre">ncurves,</span> <span class="pre">nump2</span> <span class="pre">=</span> <span class="pre">1,</span> <span class="pre">and</span> <span class="pre">nump3</span> <span class="pre">=</span> <span class="pre">1)</span></code></p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">If</span> <span class="mi">2</span> <span class="n">parameter</span> <span class="n">curves</span> <span class="n">are</span> <span class="n">being</span> <span class="nb">input</span>

<span class="n">For</span> <span class="n">each</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span><span class="n">to</span> <span class="n">ncurve</span> <span class="p">(</span><span class="n">d4</span> <span class="o">=</span> <span class="mi">1</span><span class="p">)</span>

<span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">),</span> <span class="n">param1</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">param2</span><span class="p">(</span><span class="n">i</span><span class="p">)</span>

<span class="n">followed</span> <span class="n">by</span> <span class="k">for</span> <span class="n">each</span> <span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">)</span>

<span class="n">time</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nump</span><span class="p">),</span> <span class="n">conc</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nump</span><span class="p">)</span>
</pre></div>
</div>
<p>Or if 3 parameter curves are being input</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">For</span> <span class="n">each</span> <span class="n">d4</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="mi">4</span><span class="p">,</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span> <span class="n">to</span> <span class="n">ncurve</span>

<span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">d4</span><span class="p">),</span> <span class="n">param1</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">param2</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">param3</span><span class="p">(</span><span class="n">i</span><span class="p">),</span> <span class="n">d4</span>

<span class="n">followed</span> <span class="n">by</span> <span class="k">for</span> <span class="n">each</span> <span class="n">nump</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">d4</span><span class="p">)</span>

<span class="n">time</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">d4</span><span class="p">,</span> <span class="n">nump</span><span class="p">),</span> <span class="n">conc</span><span class="p">(</span><span class="n">i</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">d4</span><span class="p">,</span> <span class="n">nump</span><span class="p">)</span>

<span class="n">out_flag</span> <span class="p">(</span><span class="n">optional</span><span class="p">)</span> <span class="o">-</span> <span class="n">keyword</span> <span class="s2">&quot;out&quot;</span>
</pre></div>
</div>
<p>Please note that all fracture-fracture curves are input followed by fracture-matrix curves, followed by matrix-fracture curves, followed by matrix-matrix curves.</p>

## Optional Input Files

<p>The data for any of the FEHM macros (with the exception of coor and elem, where use of a separate geometry input file is handled through control file input) may be entered in an alternate input file. To use this option the keyword ‘file’ must appear on the input line immediately following the control statement (macro name). The line immediately following this keyword will contain the name of the alternate input file. The contents of the alternate input file consist of the regular full macro description: the macro name followed by the data. Note that data from only one macro may be entered per auxilliary file. The entries in the optional input file may be preceded or followed by comments using the “#” designator (see discussion on <a href="#id12"><span class="problematic" id="id13">`Comments`_</span></a> may be entered in the input file by beginning a line with a ‘#’ symbol (the ‘#’ symbol must be found in the first column of the line). Comments may precede or follow macro blocks but may not be found within a block.). As with regular macro input, comments may not be embedded within the data block.</p>

<blockquote>
<div><p>Group 1 - <code class="docutils literal notranslate"><span class="pre">LOCKEYWORD</span></code></p>
<p>Group 2 - <code class="docutils literal notranslate"><span class="pre">LOCFILENAME</span></code></p>
</div></blockquote>

| Input Variable | Format | Description |
|:---------------|:-----------|:--------|
| LOCKEYWORD | character*4 | Keyword ‘file’ to designate an auxiliary input file is used.|
| LOCFILENAME | character*100 | Name of the optional data input file.|

<p>The following illustrate the use of an optional input file and its contents. In this example, optional file “rockfile” is located in the current working directory. Input for macro “rock” is described in Control statement <a class="reference external" href="MacroRock.md">rock</a> (required).</p>

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">rock</span>
<span class="n">file</span>
<span class="n">rockfile</span>
</pre></div>
</div>

<p>File “rockfile”:</p>
<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="c1"># Auxiliary file used for rock macro input</span>
<span class="n">rock</span>
 <span class="mi">1</span>   <span class="mi">140</span>     <span class="mi">1</span>       <span class="mf">2563.</span>   <span class="mf">1010.</span>   <span class="mf">0.3500</span>

<span class="c1"># End of rock macro input</span>
</pre></div>
</div>


## Option to disable a macro

<p>The data from any input macro may be omitted from a simulation by including the “off” keyword on the macro line and terminating the macro with an end macro statement. This also allows the inclusion of an end macro statement for any input macro. The end macro will follow the last specified line of input. If a macro is normally terminated with an end keyword, the macro id is appended to that keyword (an additional end macro line is not added). This facilitates experimentation with input options without the need to comment out or remove unwanted macros from the input file.</p>

<p>In the following example the perm macro is turned off and the hyco macro is used in its place.</p>

<div class="code highlight-default notranslate"><div class="highlight"><pre><span></span><span class="n">perm</span> <span class="n">off</span>
     <span class="mi">1</span>       <span class="mi">0</span>       <span class="mi">0</span>       <span class="mf">1.e-12</span>  <span class="mf">1.e-12</span>  <span class="mf">1.e-12</span>
<span class="n">end</span> <span class="n">perm</span>
<span class="n">hyco</span>
     <span class="mi">1</span>       <span class="mi">0</span>       <span class="mi">0</span>       <span class="mf">3.8e-3</span>  <span class="mf">3.8e-3</span>  <span class="mf">3.8e-3</span>
<span class="n">end</span> <span class="n">hyco</span>
</pre></div>
</div>

  <div role="contentinfo">
    <p>
        
        &copy; Copyright 2018, Los Alamos National Laboratory

    </p>
  </div>
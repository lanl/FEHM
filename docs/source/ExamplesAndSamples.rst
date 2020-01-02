****************************
Examples and Sample Problems
****************************

The following describes execution of the FEHM code. `Constructing an Input File`_ discusses the construction of an input file. `Code Execution`_ illustrates the entire procedure for executing the FEHM code using terminal input. Example 1 describes the setup and results from a simple 2-D heat conduction simulation. The remaining sections provide more complex example problems and deal only with problem setup and expected results.

Constructing an Input File
==========================

FEHM is a very general simulation code. Thus it is preferable to discuss the construction of an input file from a problem oriented point of view. In what follows the needs of the physical problem (initial conditions, boundary conditions, etc.) will be addressed in terms of the macro statements.

**Initial conditions**. These are needed for every problem, even if it is a steady state simulation. If the simulation is comprised of fully saturated water flow or heat conduction only, then the appropriate control statement would be `init <MacroInit.html>`_. The use of ``init`` also allows the specification of initial temperature and pressure (gravity) gradients. If two phase flow is prescribed (thermal or isothermal) then entering the initial conditions through the control statement `pres <MacroPres.html>`_ is more convenient. Initial values for noncondensible gas are handled in the `ngas <MacroNgas.html>`_ control statement. It should be remembered that if a restart file is present, those values will have precedence over values input in control statement ``init`` but not over values input in control statement ``pres`` Solute initial conditions are prescribed through the control statement `trac <MacroTrac.html>`_.

**Boundary conditions**. Fluid and heat flow boundary conditions can be prescribed through control statements ``pres`` `bound <MacroBoun.html>`_, `flow <MacroFlow.html>`_, and `hflx <MacroHflx.html>`_. Boundary conditions are entered with ``pres`` by specifying a negative phase state designation (the code will actually use the absolute value of the phase state designation). In this case the code will keep the variable values constant at whatever value was prescribed in **pres**. Flowing pressures are input with the boun or ``flow`` control statement. Solute boundary conditions are prescribed through the control statement ``trac``

**Material and Energy Balance Equations**. The choice of the coupled system equations is made in control statements `sol <MacroSol.html>`_, ``ngas``, and `air <MacroAir.html>`_.

**Rock or Media Properties**. These are found in the `rock <MacroRock.html>`_ and ``perm`` control statements.

**Fluid Properties**. These are found in control statement `eos <MacroEos.html>`_, which is optional. If **eos** is not invoked, then the properties of water and air included in the code are used. Relative permeabilities, depending on both the fluid and media type, are found in control statement `rlp <MacroRlp.html>`_.

**Mesh Geometry and Nodal Coordinates**. This geometry information is found in control statements `coor <MacroCoor.html>`_ and `elem <MacroElem.html>`_. This information is usually created with a mesh generation program.

**Simulation Time**. The time stepping information including printout intervals and time step sizing is found in control statement `time <MacroTime.html>`_.

**Numerics**. Convergence criteria, upwinding parameters, fill-in for the preconditioned conjugate gradient solver and geometry type (2-D, 3-D, radial) are entered with control statement `ctrl <MacroCtrl.html>`_.

**Advanced Iteration Control**. Reduced degree of freedom methods are invoked with the `iter <MacroIter.html>`_ control statement. One important quantity entered with this statement is the maximum time for the job to run on the computer.

**Sources and Sinks**. These are input with the control statements ``boun`` or ``flow`` Care must be taken as the parameters have different meanings for different physical models.

The following table shows required and optional macros listed by the type of problem being simulated. See the Alphabetic Macro List for all available macros and their definitions.

.. _MacroTable:

Macro Control Statements by Problem Type
----------------------------------------

Heat Conduction
---------------

+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| Required Macros                                                                                         | Optional Macros                                                    |
+=========================================================================================================+====================================================================+
| `title <Macros/MacroTitle.html>`_                                                                       | `cont <Macros/MacroCont.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `boun <Macros/MacroBoun.html>`_ or `flow <Macros/MacroFlow.html>`_ or `hflx <Macros/MacroHflx.html>`_   | `finv <Macros/MacroFinv.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `cond <Macros/MacroCond.html>`_                                                                         | `flo2 <Macros/MacroFlo2.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `coor <Macros/MacroCoor.html>`_                                                                         | `flxo <Macros/MacroFlxo.html>`_ or `flxz <Macros/MacroFlxz.html>`_ |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `ctrl <Macros/MacroCtrl.html>`_                                                                         | `iter <Macros/MacroIter.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `elem <Macros/MacroElem.html>`_                                                                         | `node <Macros/MacroNode.html>`_ or `nod2 <Macros/MacroNod2.html>`_ |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `init <Macros/MacroInit.html>`_ or `pres <Macros/MacroPres.html>`_                                      | `renu <Macros/MacroRenu.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `rock <Macros/MacroRock.html>`_                                                                         | `rflx <Macros/MacroRflx.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `sol <Macros/MacroSol.html>`_                                                                           | text or comments (#)                                               |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `time <Macros/MacroTime.html>`_                                                                         | `user <Macros/MacroUser.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| `stop <Macros/MacroStop.html>`_                                                                         | `vcon <Macros/MacroVcon.html>`_                                    |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
|                                                                                                         | `zone <Macros/MacroZone.html>`_ or `zonn <Macros/MacroZonn.html>`_ |
+---------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+


Water / Water Vapor / Heat Equivalent Continuum, Dual Porosity,Dual Permeability
--------------------------------------------------------------------------------

+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| Required Macros                                                                                                | Optional Macros                                                       |
+================================================================================================================+=======================================================================+
| `title <Macros/MacroTitle.html>`_                                                                              | `cden <Macros/MacroCden.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `boun <Macros/MacroBoun.html>`_ or `flow <Macros/MacroFlow.html>`_ or `hflx <Macros/MacroHflx.html>`_          | `cont <Macros/MacroCont.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `cond <Macros/MacroCond.html>`_                                                                                | `eos <Macros/MacroEos.html>`_                                         |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `coor <Macros/MacroCoor.html>`_                                                                                | `exrl <Macros/MacroExrl.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `ctrl <Macros/MacroCtrl.html>`_                                                                                | `finv <Macros/MacroFinv.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `elem <Macros/MacroElem.html>`_                                                                                | `flo2 <Macros/MacroFlo2.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `init <Macros/MacroInit.html>`_ or `pres <Macros/MacroPres.html>`_                                             | `flxo <Macros/MacroFlxo.html>`_ or `flxz <Macros/MacroFlxz.html>`_    |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `perm <Macros/MacroPerm.html>`_                                                                                | `fper <Macros/MacroFper.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `rlp <Macros/MacroRlp.html>`_                                                                                  | `gdpm <Macros/MacroGdpm.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `rock <Macros/MacroRock.html>`_                                                                                | `hflx <Macros/MacroHflx.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `sol <Macros/MacroSol.html>`_                                                                                  | `iter <Macros/MacroIter.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `time <Macros/MacroTime.html>`_                                                                                | `node <Macros/MacroNode.html>`_ or `nod2 <Macros/MacroNod2.html>`_    |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `stop <Macros/MacroStop.html>`_                                                                                | `ppor <Macros/MacroPpor.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                                | `renu <Macros/MacroRenu.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| dual (* only)                                                                                                  | `rflx <Macros/MacroRflx.html>`_                                       |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| dpdp (** only)                                                                                                 | `rxn <Macros/MacroRxn.html>`_                                         |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                                | text or comments (#)                                                  |
|                                                                                                                +-----------------------------------------------------------------------+
|                                                                                                                | `trac <Macros/MacroTrac.html>`_                                       |
|                                                                                                                +-----------------------------------------------------------------------+
|                                                                                                                | `user <Macros/MacroUser.html>`_ or `userc <Macros/MacroUserc.html>`_  |
|                                                                                                                +-----------------------------------------------------------------------+
|                                                                                                                | `vcon <Macros/MacroVcon.html>`_                                       |
|                                                                                                                +-----------------------------------------------------------------------+
|                                                                                                                | `velo <Macros/MacroVelo.html>`_                                       |
|                                                                                                                +-----------------------------------------------------------------------+
|                                                                                                                | `zone <Macros/MacroZone.html>`_ or `zonn <Macros/MacroZonn.html>`_    |
+----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+


Air / Water / No Heat Equivalent Continuum, Dual Porosity, Dual Permeability 
----------------------------------------------------------------------------

+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| Required Macros                                                                                       | Optional Macros                                                       |
+=======================================================================================================+=======================================================================+
| `title <Macros/MacroTitle.html>`_                                                                     | `adif <Macros/MacroAdif.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `boun <Macros/MacroBoun.html>`_ or `flow <Macros/MacroFlow.html>`_ or `hflx <Macros/MacroHflx.html>`_ | `cden <Macros/MacroCden.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `cond <Macros/MacroCond.html>`_                                                                       | `cont <Macros/MacroCont.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `coor <Macros/MacroCoor.html>`_                                                                       | `eos <Macros/MacroEos.html>`_                                         |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `ctrl <Macros/MacroCtrl.html>`_                                                                       | `finv <Macros/MacroFinv.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `elem <Macros/MacroElem.html>`_                                                                       | `flo2 <Macros/MacroFlo2.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `init <Macros/MacroInit.html>`_ or `pres <Macros/MacroPres.html>`_                                    | `flxo <Macros/MacroFlxo.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `ngas <Macros/MacroNgas.html>`_                                                                       | `fper <Macros/MacroFper.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `perm <Macros/MacroPerm.html>`_                                                                       | `gdpm <Macros/MacroGdpm.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `rlp <Macros/MacroRlp.html>`_                                                                         | `iter <Macros/MacroIter.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `rock <Macros/MacroRock.html>`_                                                                       | `node <Macros/MacroNode.html>`_ or `nod2 <Macros/MacroNod2.html>`_    |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `sol <Macros/MacroSol.html>`_                                                                         | `ppor <Macros/MacroPpor.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `time <Macros/MacroTime.html>`_                                                                       | `renu <Macros/MacroRenu.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| `stop <Macros/MacroStop.html>`_                                                                       | `rflx <Macros/MacroRflx.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `rxn <Macros/MacroRxn.html>`_                                         |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| dual (*only)                                                                                          | `szna <Macros/MacroSzna.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
| dpdp (**only)                                                                                         | text or comments (#)                                                  |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `trac <Macros/MacroTrac.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `user <Macros/MacroUser.html>`_ or `userc <Macros/MacroUserc.html>`_  |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `vapl <Macros/MacroVapl.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `vcon <Macros/MacroVcon.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `velo <Macros/MacroVelo.html>`_                                       |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                                                       | `zone <Macros/MacroZone.html>`_ or `zonn <Macros/MacroZone.html>`_    |
+-------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------+



Air / Water / No Heat Equivalent Continuum, Dual Porosity, Dual Permeability
----------------------------------------------------------------------------

+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| Required Macros                                                      | Optional Macros                                                       |
+======================================================================+=======================================================================+
| `title <Macros/MacroTitle.html>`_                                    | `bous <Macros/MacroBous.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `airwater <Macros/MacroAirwater.html>`_                              | `cont <Macros/MacroCont.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `boun <Macros/MacroBoun.html>`_ or `flow <Macros/MacroFlow.html>`_   | `eos <Macros/MacroEos.html>`_                                         |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `coor <Macros/MacroCoor.html>`_                                      | `exri <Macros/MacroExri.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `ctrl <Macros/MacroCtrl.html>`_                                      | `finv <Macros/MacroFinv.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `elem <Macros/MacroElem.html>`_                                      | `flo2 <Macros/MacroFlo2.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `init <Macros/MacroInit.html>`_ or `pres <Macros/MacroPres.html>`_   | `flxo <Macros/MacroFlxo.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `node <Macros/MacroNode.html>`_ or `nod2 <Macros/MacroNod2.html>`_   | `fper <Macros/MacroFper.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `perm <Macros/MacroPerm.html>`_                                      | `gdpm <Macros/MacroGdpm.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `rock <Macros/MacroRock.html>`_                                      | `head <Macros/MacroHead.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `sol <Macros/MacroSol.html>`_                                        | `iter <Macros/MacroIter.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `time <Macros/MacroTime.html>`_                                      | `ppor <Macros/MacroPpor.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| `stop <Macros/MacroStop.html>`_                                      | `pres <Macros/MacroPres.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | `renu <Macros/MacroRenu.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| dual (*only)                                                         | `rlp <Macros/MacroRlp.html>`_                                         |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
| dpdp (only)                                                          | `rxn <Macros/MacroRxn.html>`_                                         |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | text or comments                                                      |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | `trac <Macros/MacroTrac.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | `user <Macros/MacroUser.html>`_ or `userc <Macros/MacroUserc.html>`_  |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | `vapl <Macros/MacroVapl.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | `velo <Macros/MacroVelo.html>`_                                       |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+
|                                                                      | `zone <Macros/MacroZone.html>`_ or `zonn <Macros/MacroZonn.html>`_    |
+----------------------------------------------------------------------+-----------------------------------------------------------------------+

Code Execution
--------------
To run FEHM, the program executable file name is entered at the system prompt:

.. code::

   <PROMPT> xfehm_v2.30

Note that executable names may vary depending on the system being used. 

The I/O file information is provided to the code from an input control file or the terminal. The default control file name is ``fehmn.files``. If a control file with the default name is present in the directory from which the code is being executed, no terminal input is required. If the default control file is not present, input prompts are written to the screen. A short description of the I/O files used by FEHM precedes the initial prompt. The following assumes the default control file was not found in the execution directory (for this example ``/home/fehm/heat2d``).

After the command ``xfehm_v2.30`` is given, the code queries the user regarding the input files, as follows:

.. code::

   Enter name for iocntl -- default file name: not using
   [(name/na or not using), RETURN = DEFAULT]

This query asks for a control file name. If a control file name is entered no further terminal input is required. The below example shows the control file that would produce the same results as the terminal responses discussed below and illustrated in TerminalQuery_. Files that are not needed for output can be represented with a blank line. If names are not provided for the write file and/or the data check file, the code will use the following defaults: fehmn.fin and fehmn.chk. Following the file names is the flag that controls terminal output. The last line of the file is the user subroutine number. Omitting these values results in no terminal output and no user subroutine call. For now, we assume a carriage return ``<cr>`` is entered and a control file is not being used. The following query will appear:

**Input control file for heat conduction example**

.. code::

   /home/fehm/heat2d/input/heat2d.in
   /home/fehm/heat2d/input/heat2d.in
   /home/fehm/heat2d/input/heat2d.in
   /home/fehm/heat2d/output/heat2d.out
   /home/fehm/heat2d/output/heat2d.fin
   /home/fehm/heat2d/output/heat2d.his
   /home/fehm/heat2d/output/heat2d.chk
   some


   Enter name for inpt -- default file name: fehmn.dat
   [(name/na or not using), RETURN = DEFAULT]

This query asks for an input file name. If a ``<cr>`` is given, the default ``fehmn.dat`` is used for the input file. We shall assume that the input file name entered is

.. code::

   input/heat2d.in

Note that a subdirectory containing the file is also given. If the file did not exist, the code would repeat the prompt for an input file. Next the code would query to determine if the prefix of the input file name (the portion of the name preceding the final "." or first space) should be used for code generated file names.

.. code::

   Do you want all file names of the form input/heat2d.* ? [(y/n), RETURN = y]
   *** Note: If "y" incoor and inzone will equal inpt ***

A ``<cr>`` will produce files with identical prefixes, including the subdirectory. If the response is negative, the code will query for the names of all required files. Assume we enter ``n``.

.. code::

   Enter name for incoor -- default file name: input/heat2d.in
   [(name/na or not using), RETURN = DEFAULT]

(See TerminalQuery_ for the remaining file name queries.) 

Next a query for terminal output appears.

.. code::

   tty output -- show all reference nodes, selected reference nodes, or none:
   [(all/some/none), RETURN = none]

An "all" reply prints out the primary node information to the terminal at every time step. A "some" reply prints a selected subset of the node information. A reply of "none" suppresses all tty output with the exception of error messages printed if code execution is terminated abnormally or when maximum number of iterations are exceeded. Assume we enter "some".

The next query concerns the subroutine USER. This subroutine is used for special purposes and is not available to the general user.

.. code::

   user subroutine number (provided to subroutine USER before every time step):
   [RETURN = none]

Assume a ``<cr>`` is entered.

The code will then print a summary of the I/O files to be used. 

The final query regards the acceptance of the file set just created. A "yes" reply denotes that the user has accepted the file set and the code proceeds with calculations. A "no" reply starts the query sequence again so I/O file names may be reentered or modified. A "stop" reply stops the current computer job.

.. code::

   If data is OK enter yes to continue, no to restart terminal input,
   or stop to end program: [(yes/no/stop), RETURN = yes]


Screen output for this example execution using terminal input and a previous version of the code is shown below. The only difference in the output is that the code version identifier and date are updated for the current version. User responses are shown in *italics*.

.. _TerminalQuery:

**Terminal query for FEHM example run**

.. code::

   <PROMPT> xfehm_v2.30sun
   version FEHM V2.30sun 03-09-15 QA:QA 09/15/2003 11:08:14
   **** Default names for I/O files ****
   control file : fehmn.files
   input file : filen.*
   geometry data file : filen.*
   zone data file : filen.*
   output file : filen.out
   read file (if it exists) : filen.ini
   write file (if it exists) : filen.fin
   history plot file : filen.his
   tracer history plot file : filen.trc
   contour plot file : filen.con
   dual or dpdp contour plot file : filen.dp
   stiffness matrix data read/write file : filen.stor
   input check file : filen.chk
   **** where ****
   "filen.*" may be 100 characters maximum. If a name is not entered
   when prompted for, a default file name is used. "fehmn.dat" is the
   default used for the input file name.
   **** note ****
   A save file and input check file are always written. If you do not
   provide a name for these files, the following defaults will be used:
   fehmn.fin, fehmn.chk
   Enter name for iocntl -- default file name: not using
   [(name/na or not using), RETURN = DEFAULT]
   <cr>
   Enter name for inpt -- default file name: fehmn.dat
   [(name/na or not using), RETURN = DEFAULT]
   input/heat2d.in
   Do you want all file names of the form input/heat2d.* ? [(y/n), RETURN = y]
   *** Note: If "y" incoor and inzone will equal inpt ***
   n
   Enter name for incoor -- default file name: input/heat2d.in
   [(name/na or not using), RETURN = DEFAULT]
   <cr>
   Enter name for inzone -- default file name: input/heat2d.in
   [(name/na or not using), RETURN = DEFAULT]
   <cr>
   Enter name for iout -- default file name: input/heat2d.out
   [(name/na or not using), RETURN = DEFAULT]
   output/heat2d.out
   Enter name for iread -- default file name: input/heat2d.ini
   [(name/na or not using), RETURN = DEFAULT]
   na
   Enter name for isave -- default file name: input/heat2d.fin
   [(name/na or not using), RETURN = DEFAULT]
   output/heat2d.fin
   Enter name for ishis -- default file name: input/heat2d.his
   [(name/na or not using), RETURN = DEFAULT]
   output/heat2d.his
   Enter name for istrc -- default file name: input/heat2d.trc
   [(name/na or not using), RETURN = DEFAULT]
   na
   Enter name for iscon -- default file name: input/heat2d.con
   [(name/na or not using), RETURN = DEFAULT]
   na
   Enter name for iscon1 -- default file name: input/heat2d.dp
   [(name/na or not using), RETURN = DEFAULT]
   na
   Enter name for isstor -- default file name: input/heat2d.stor
   [(name/na or not using), RETURN = DEFAULT]
   na
   Enter name for ischk -- default file name: input/heat2d.chk
   [(name/na or not using), RETURN = DEFAULT]
   output/heat2d.chk
   tty output -- show all reference nodes, selected reference nodes, or none:
   [(all/some/none), RETURN = none]
   some
   user subroutine number (provided to subroutine USER before every time step):
   [RETURN = none]
   <cr>
   First reference output node will be written to tty
   File purpose - Variable - Unit number - File name
   control - iocntl - 0 - not using
   input - inpt - 11 - input/heat2d.in
   geometry - incoor - 11 - input/heat2d.in
   zone - inzone - 11 - input/heat2d.in
   output - iout - 14 - output/heat2d.out
   initial state - iread - 0 - not using
   final state - isave - 16 - output/heat2d.fin
   time history - ishis - 17 - output/heat2d.his
   time his.(tr) - istrc - 18 - not using
   contour plot - iscon - 19 - not using
   con plot (dp) - iscon1 - 20 - not using
   fe coef stor - isstor - 21 - not using
   input check - ischk - 22 - output/heat2d.chk
   Value provided to subroutine user: not using
   If data is OK enter yes to continue, no to restart terminal input,
   or stop to end program: [(yes/no/stop), RETURN = yes]
   <cr>

 
Heat Conduction in a Square
===========================

This simple 2-D problem is used to illustrate input file construction and basic output. Heat conduction in a 1 meter square with an initial temperature, :math:`T_0 = 200 ^{\circ}C`, is modeled after a surface temperature, :math:`T_s = 100 ^{\circ}C`, is imposed at time, :math:`t = 0` (`See Schematic diagram of 2-D heat conduction problem. <FEHM-UM.9.3.html>`_). The input parameters used for the heat conduction problem are defined in `Input Parameters for the 2-D Heat Conduction Problem <FEHM-UM.9.3.htm#71293>`_.The finite element mesh for this problem is shown in `Finite element mesh used for 2-D heat conduction problem. <FEHM-UM.9.3.htm#49464>`_. Only a quarter of the square needs to be modeled because of problem symmetry. 

.. _TerminalQueryTable:

.. figure: FEHM-UM.9.3-1.gif
   :caption: Schematic diagram of 2-D heat conduction problem

+-----------------------------------------+----------------------------------------------+-----------------------------------+
|                                       Input Parameters for the 2-D Heat Conduction Problem                                 |
+=========================================+==============================================+===================================+
|  Rock thermal conductivity              | :math:`\kappa r`                             | :math:`2.7 \frac{W}{m \cdot K}`   |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
|  Rock density                           | :math:`\rho r`                               | :math:`2700 \frac{kg}{m^3}`       |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
|  Rock specific heat                     | :math:`Cr`                                   | :math:`1000 \frac{J}{kg \cdot K}` |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
|  Width                                  | :math:`a`                                    | :math:`0.5 m`                     |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
|  Length                                 | :math:`b`                                    | :math:`0.5 m`                     |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
|  Initial temperature                    | :math:`T_0`                                  | :math:`200 ^{\circ}C`             |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
| Surface temperaturefor all x, y = 0.5 m | :math:`T_s`                                  | :math:`100 ^{\circ}C`             |
+-----------------------------------------+----------------------------------------------+-----------------------------------+
|  Rock thermal diffusivity               | :math:`\kappa = \frac{\kappa_r}{\rho_r C_r}` |                                   |
+-----------------------------------------+----------------------------------------------+-----------------------------------+


.. figure: FEHM-UM.9.3-5.gif
   :caption: Finite element mesh used for 2-D heat conduction problem

The input file (see `FEHM input file for heat conduction example (heat2d.in). <FEHM-UM.9.3.htm#26275>`_) uses optional macro control statement node (output nodes) and the required macro control statements sol (solution specification - heat transfer only), init (initial value data), rock (rock properties), cond (thermal conductivities), perm (permeabilities), time (simulation timing data), ctrl (program control parameters), coor (node coordinates), elem (element node data), and stop. For this problem macro control statement flow is also used to set the temperature boundary conditions. A portion of the output file is reproduced in `FEHM output from the 2-D heat conduction example. <FEHM-UM.9.3.htm#50551>`_.

FEHM input file for heat conduction example (``heat2d.in``)
===========================================================

.. code::

   ***** 2-D Heat Conduction Model (2X2 rectangles) *****node
   2
   7	5
   sol
   
   -1	-1
   init
   
   10.	0.	200.	0.	0.	200.	0.	0.
   rock
   
   1	9	1	2700.	1000.	0. 
   
   cond
   
   1	9	1	2.7e-00	2.7e-00	2.7e-00
   
   perm
   
   1	9	1	1.e-30	1.e-30	1.e-30
   
   flow
   
   1	3	1	 10.00	-100.00	1.e03
   3	9	3	 10.00	-100.00	1.e03
   
   time
   
   0.005	4.00	1000	10	1994	02
   
   ctrl
   
   40	1.e-04	08
   1	9	1	1
   
   1.0	0.0	1.0
   10	1.0	0.00005	0.005
   1	0
   coor Feb 23, 1994 11:39:40 
   
   	9
   
   	1	0.	0.50	0.
   
   	2	0.25	0.50	 0.
   
   	3	0.50	0.50	 0.
   
   	4	0.	0.25	 0.
   
   	5	0.25	0.25	 0.
   
   	6	0.50	0.25	 0.
   
   	7	0.	0.	0.
   
   	8	0.25	 0.	0.
   
   	9	0.50	 0.	0.
   
   
   
   elem
   
   	4	4
   
   	1	4	5	2	1
   
   	2	5	6	3	2
   
   	3	7	8	5	4
   
   	4	8	9	6	5
   
   
   
   stop


FEHM output from the 2-D heat conduction example
================================================

.. code::

   FEHM V2.10 00-06-28 08/07/2000 13:25:08
   ***** 2-D Heat Conduction Model *****
   File purpose - Variable - Unit number - File name
   control - iocntl - 0- not using
   input - inpt - 11- heat2d.in
   geometry- incoor- 11- heat2d.in
   zone - inzone- 11- heat2d.in
   output - iout - 14- heat2d.out
   initial state- iread- 0- not using
   final state- isave- 16- fehmn.fin
   time history- ishis- 17- heat2d.his
   time his.(tr)- istrc- 0- not using
   contour plot- iscon- 0- not using
   con plot (dp)- iscon1- 0- not using
   fe coef stor- isstor- 0- not using
   input check- ischk- 22- fehmn.chk
   Value provided to subroutine user: not using
   **** input title : coor**** incoor = 11 ****
   **** input title : elem**** incoor = 11 ****
   **** input title : stop**** incoor = 11 ****
   **** input title : node**** inpt = 11 ****
   **** input title : sol**** inpt = 11 ****
   **** input title : init**** inpt = 11 ****
   **** input title : rock**** inpt = 11 ****
   **** input title : cond**** inpt = 11 ****
   **** input title : perm**** inpt = 11 ****
   **** input title : flow**** inpt = 11 ****
   **** input title : time**** inpt = 11 ****
   **** input title : ctrl**** inpt = 11 ****
   **** input title : stop**** inpt = 11 ****
   BC to BC connection(s) found(now set=0.0)
   BC to BC connection(s) found(now set=0.0)
   pressures and temperatures set by gradients
   >>>reading nop from file nop.temp.....
   >>>reading nop was succesful.....
   storage needed for ncon43 available 43
   storage needed for nop43 available 46
   storage needed for a matrix33 available 33
   storage needed for b matrix33 available 46
   storage needed for gmres81 available 81
   storage available for b matrix resized to 33<<<<<<
   time for reading input, forming coefficients 0.204E-01
   **** analysis of input data on file fehmn.chk ****
   *********************************************************************
   Time Step 1
   Timing Information
   Years Days Step Size (Days)
   0.136893E-04 0.500000E-02 0.500000E-02
   Cpu Sec for Time Step = 0.8081E-03 Current Total = 0.2650E-02
   Equation Performance
   Number of N-R Iterations: 1
   Avg # of Linear Equation Solver Iterations: 3.0
   Number of Active Nodes: 9.
   Total Number of Newton-Raphson Iterations: 1 , Solver: 3
   Largest Residuals
   EQ1 R= 0.1660E-07 node= 5 x=0.2500 y=0.2500 z= 1.000
   Node Equation 1 Residual Equation 2 Residual
   7 0.111444E-07 0.185894E-01
   5 0.165983E-07 0.135450E+01
   Nodal Information (Water)
   source/sink source/sink
   Node p(MPa) e(MJ) l sat temp(c) (kg/s) (MJ/s)
   7 10.000 0.00 0.000 199.981 0. 0.
   5 10.000 0.00 0.000 198.645 0. 0.
   Global Mass & Energy Balances
   Total mass in system at this time:0.000000E+00 kg
   Total mass of steam in system at this time:0.000000E+00 kg
   Total enthalpy in system at this time:0.105123E+03 MJ
   Water discharge this time step:0.000000E+00 kg (0.000000E+00 kg/s)
   Water input this time step:0.000000E+00 kg (0.000000E+00 kg/s)
   Total water discharge:0.000000E+00 kg (0.000000E+00 kg/s)
   Total water input:0.000000E+00 kg (0.000000E+00 kg/s)
   Enthalpy discharge this time step:0.297800E+02 MJ (0.689352E-01 MJ/s)
   Enthalpy input this time step:0.000000E+00 MJ (0.000000E+00 MJ/s)
   Total enthalpy discharge:0.297800E+02 MJ (0.689352E-01 MJ/s)
   Total enthalpy input:0.297800E+02 MJ (0.689352E-01 MJ/s)
   Net kg water discharge (total out-total in):0.000000E+00
   Net MJ discharge (total out-total in):0.000000E+00
   Conservation Errors: 0.000000E+00 (mass), -0.100326E+01 (energy)
   *********************************************************************
   Time Step 11
   .
   .
   .
   *********************************************************************
   Time Step 801
   Timing Information
   Years Days Step Size (Days)
   0.109515E-01 0.400005E+01 0.500000E-04
   Cpu Sec for Time Step = 0. Current Total = 4.533
   Equation Performance
   Number of N-R Iterations: 1
   Avg # of Linear Equation Solver Iterations: 2.0
   Number of Active Nodes: 9.
   Total Number of Newton-Raphson Iterations: 801 , Solver: 2402
   Largest Residuals
   EQ1 R= 0.9774E-13 node= 7 x= 0.000 y= 0.000 z= 1.000
   Node Equation 1 Residual Equation 2 Residual
   7 0.977369E-13 0.186062E-04
   5 0.621566E-13 0.930309E-05
   Nodal Information (Water)
   source/sink source/sink
   Node p(MPa) e(MJ) l sat temp(c) (kg/s) (MJ/s)
   7 10.000 0.00 0.000 100.230 0. 0.
   5 10.000 0.00 0.000 100.115 0. 0.
   Global Mass & Energy Balances
   Total mass in system at this time:0.000000E+00 kg
   Total mass of steam in system at this time:0.000000E+00 kg
   Total enthalpy in system at this time:0.675565E+02 MJ
   Water discharge this time step:0.000000E+00 kg (0.000000E+00 kg/s)
   Water input this time step:0.000000E+00 kg (0.000000E+00 kg/s)
   Total water discharge:0.000000E+00 kg (0.000000E+00 kg/s)
   Total water input:0.000000E+00 kg (0.000000E+00 kg/s)
   Enthalpy discharge this time step:0.455636E-05 MJ (0.105471E-05 MJ/s)
   Enthalpy input this time step:0.000000E+00 MJ (0.000000E+00 MJ/s)
   Total enthalpy discharge:0.673463E+02 MJ (0.155894E+02 MJ/s)
   Total enthalpy input:0.673463E+02 MJ (0.155894E+02 MJ/s)
   Net kg water discharge (total out-total in):0.000000E+00
   Net MJ discharge (total out-total in):0.000000E+00
   Conservation Errors: 0.000000E+00 (mass), -0.100144E+01 (energy)
   simulation ended: days 4.00 timesteps 801
   total N-R iterations = 801
   total solver iterations = 2402
   total code time(timesteps) = 0.526277
   **** -------------------------------------------------------------- ****
   **** This program for ****
   **** Finite Element Heat and Mass Transfer in porous media****
   **** -------------------------------------------------------------- ****
   **** Version : FEHM V2.10 00-06-28 ****
   **** End Date : 08/07/2000 ****
   **** Time : 13:25:08 ****
   **** --------------------------------------------------------------****


The analytical solution for 2-D heat conduction (Carslaw and Jaeger, 1959) is given by

.. image: FEHM-UM.9.3-6.gif
.. image: FEHM-UM.9.3-7.gif

The below image shows a plot of the simulation results compared to the analytical solution for the selected output nodes at :math:`x = y = 0 m` and :math:`x = y = 0.25 m`. 

.. figure: FEHM-UM.9.3-9.gif
   :caption: Comparison of analytical and model solution for 2-D heat conduction.
 
 
DOE Code Comparison Project, Problem 5, Case A
==============================================

This problem involves multiphase flow in a 2-D horizontal reservoir. The problem is characterized by a moving two-phase region, i.e., the fluid produced at the production well is replaced by cold water recharge over one of the outer boundaries. The problem parameters are given below and the geometry and boundary conditions are shown in the below schematic. Of particular note are the variable initial temperature field, provided to the code through a read file (see `iread <DataFiles.html#Readfileiread>`_), and the prescribed pressure and temperature on the right boundary. A partial listing of the input file is provided in `FEHM input file for DOE problem <FEHM-UM.9.4.htm#30248>`_. In addition to the required macros, macro flow is used to specify the pressure and temperature boundary condition and the production flow rate. Macro rlp is used to set the residual liquid and gas saturations. 

Input Parameters for the DOE Code Comparison Project Problem
------------------------------------------------------------

+-------------------------------+------------------------------------------+-----------------------------------+
| Parameter                     | Symbol                                   | Value                             |
+===============================+==========================================+===================================+
| Reservoir permeability        | :math:`k`                                | :math:`2.5 \cdot 10^{-14} m^2`    |
+-------------------------------+------------------------------------------+-----------------------------------+
| Reservoir porosity            | :math:`\phi`                             | :math:`0.35`                      |
+-------------------------------+------------------------------------------+-----------------------------------+
| Rock thermal conductivity     | :math:`\kappa r`                         | :math:`1 \frac{W}{m \cdot K}`     |
+-------------------------------+------------------------------------------+-----------------------------------+
| Rock density                  | :math:`\rho r`                           | :math:`2563 \frac{kg}{m^3}`       |
+-------------------------------+------------------------------------------+-----------------------------------+
| Rock specific heat            | :math:`C r`                              | :math:`1010 \frac{J}{kg \cdot K}` |
+-------------------------------+------------------------------------------+-----------------------------------+
| Reservoir length              | :math:`x`                                | :math:`300 m`                     |
+-------------------------------+------------------------------------------+-----------------------------------+
| Reservoir thickness           | :math:`y`                                | :math:`200 m`                     |
+-------------------------------+------------------------------------------+-----------------------------------+
| Liquid residual saturation    | slr                                      | :math:`0.3`                       |
+-------------------------------+------------------------------------------+-----------------------------------+
| Gas residual saturation       | sgr                                      | :math:`0.1`                       |
+-------------------------------+------------------------------------------+-----------------------------------+
| Reservoir discharge           | qm                                       | :math:`0.05 \frac{kg}{m \cdot s}` |
+-------------------------------+------------------------------------------+-----------------------------------+
| Initial Pressure              | :math:`P_0`                              | :math:`3.6 MPa`                   |
+-------------------------------+------------------------------------------+-----------------------------------+
| Production well coordinates:  | :math:`x = 62.5 m`, :math:`y = 62.5 m`   |                                   |
+-------------------------------+------------------------------------------+-----------------------------------+
| Observation well coordinates: | :math:`x = 162.5 m`, :math:`y = 137.5 m` |                                   |
+-------------------------------+------------------------------------------+-----------------------------------+


.. figure: FEHM-UM.9.4-4.gif
   :caption: Initial temperature distribution (T in oC, r in m) where r = sqrt(x+y)

.. figure: FEHM-UM.9.4-6.gif
   :caption: Schematic diagram of the geometry and boundary conditions for the DOE code comparison project problem.

FEHM input file for DOE problem
-------------------------------

.. code::

   *** DOE Code Comparison Project, Problem 5, Case A ***
   node
   2
   50 88
   sol
   1 1
   init
   3.6 0. 240. 0. 0. 240. 0. 0.
   rlp
   2 0.3 0.1 0.0 0.0
   1 140 1 1
   rock
   1 140 1 2563. 1010. 0.35
   cond
   1 140 1 1.00e-00 1.00e-00 1.00e-00
   perm
   1 140 1 2.5e-14 2.5e-14 0.e-00
   flow
   88 88 1 0.050 -25.00 0.
   14 140 14 3.600 -160.00 1.
   time
   30.0 3650. 10000 1000 1994 03
   ctrl
   40 1.e-07 08
   1 140 1 1
   1.0 0.0 1.0
   40 1.2 0.1 60.
   1 0
   coor
   140
   ...
   elem
   4 117
   ... stop


There is no analytical solution for this problem, but six researchers produced results for the DOE code comparison project (Molloy, 1980). The reader is referred to this reference for a more detailed discussion of this problem and the code comparison. Results from this problem are compared to those for the other codes, obtained from Molloy (1980), as a check on FEHM. The results for the outlet temperature, shown in `Comparison of FEHM production well temperatures with results from other codes <FEHM-UM.9.4.htm#30437>`_, are in excellent agreement with the other codes. The results for the outlet pressure and pressure at an observation well 125 m distant, `Comparison of FEHM production and observation well pressure drops with results from other codes <FEHM-UM.9.4.htm#89269>`_, are also in good agreement with the other codes. Contour plots of pressure and temperature at the end of the simulation were also generated for this problem and are shown in `Contour plot of pressure at ten years for the DOE problem <FEHM-UM.9.4.htm#99688>`_ and `Contour plot of temperature at ten years for the DOE problem <FEHM-UM.9.4.htm#58720>`_. 

.. figure: FEHM-UM.9.4-7.gif
   :caption: 1. Comparison of FEHM production well temperatures with results from other codes.

.. figure: FEHM-UM.9.4-8.gif
   :caption: 2. Comparison of FEHM production and observation well pressure drops with results from other codes.

.. figure: FEHM-UM.9.4-9.gif
   :caption: 3. Contour plot of pressure at ten years for the DOE problem.

.. figure: FEHM-UM.9.4-10.gif
   :caption: 4. Contour plot of temperature at ten years for the DOE problem.

.. figure: FEHM-UM.9.4-11.gif
   :caption: 5. Comparison of FEHM production and observation well pressure drops with results from other codes.

 
Reactive Transport Example
--------------------------
This one-dimensional example demonstrates the use of the reactive transport module of FEHM. The application of this simulation is the transport of cobalt (Co) in groundwater. Radioactive cobalt is present in the subsurface at several DOE sites. Although its presence as a divalent cation implies that it should sorb strongly to most soils, its migration rate has been shown to be greater than expected due to complexation with EDTA, a decontaminating agent also found in the subsurface of these sites. Much experimental work has gone into studying the transport of Co as CoEDTA, a much less strongly sorbed species. The chemical reactions and equilibrium or rate constants used to perform this simulation are:

.. image: FEHM-UM.9.5-1.gif
.. image: FEHM-UM.9.5-3.gif
.. image: FEHM-UM.9.5-5.gif
.. image: FEHM-UM.9.5-8.gif
.. image: FEHM-UM.9.5-11.gif
.. image: FEHM-UM.9.5-14.gif

.. figure: FEHM-UM.9.5-16.gif
   :caption: Schematic drawing of the geometry and boundary conditions for the cobalt transport problem

+-----------------------------------------------------------------+---------------------+---------------------------------+
|                                             Input Parameters for the Reactive Transport Test Problem                    |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Parameter                                                       | Symbol              | Value                           |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Reactor Length                                                  | :math:`L`           | 10 m                            |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Node spacing                                                    | :math:`\Delta l`    | 0.1 m                           |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Fluid Density                                                   | :math:`\rho_f`      | :math:`1000\:kg/m^3`            |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Bulk Rock Density                                               | :math:`\rho_b`      | :math:`1500\:kg/m^3`            |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Porosity                                                        | :math:`\phi`        | :math:`0.4`                     |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Pore Water Velocity                                             | :math:`u`           | :math:`1\:m/hr`                 |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Dispersivity                                                    | :math:`\alpha`      | :math:`0.05 m`                  |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Time step (tracer)                                              | :math:`\Delta t`    | :math:`0.09 - 360 s`            |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Total elapsed time                                              | :math:`t`           | 7.25 days                       |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Pressure                                                        | :math:`P_0`         | :math:`1.0\:MPa`                |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Co Inlet Concentration                                          | :math:`C_{in\:Co}`  | :math:`3.1623 \cdot 10^{-5}\:M` |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Fe Inlet Concentration                                          | :math:`C_{in\:Fe}`  | :math:`0\:M`                    |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| EDTA Inlet Concentration                                        | :math:`C_{in\:EDTA}`| :math:`3.1623 \cdot 10^{-5}\:M` |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| Boundary conditions:	At l = 0, u = 1 m/hr; At l = 1, P = 1 MPa |                     |                                 |
+-----------------------------------------------------------------+---------------------+---------------------------------+
| ‡Flow rate: :math:`q = 0.5556 kg/s`                             |                     |                                 |
+-----------------------------------------------------------------+---------------------+---------------------------------+

FEHM input file for reactive transport problem
----------------------------------------------

.. code::

   COMPARE FEHMN and PDREACT: Linear Sorption w/ Surface Exchange
   cond
   1 202 1 2.7 2.7 2.7
   ctrl
   50 1e-6 8
   1 202 1 2
   1 0 0.5
   25 2. 1.e-6 1.e-1
   1 0
   flow
   1 202 101 -0.05556 -25 0
   101 202 101 1. -25 -1
   init
   1. 25 25 0 1000 25 0 0
   node
   1
   202
   perm
   1 202 1 5.0e-13 5.0e-30 5.0e-30
   rock
   1 202 1 1500 1000 0.4
   sol
   1 -1
   time
   1.e-6 7.25 1000 10 92 11
   # solute 1: Total Cobalt Concentration
   # solute 2: Total Iron Concentration
   # solute 3: Total EDTA Concentration
   # solute 4: CoEDTA adsorbed concentration
   # solute 5: Co adsorbed concentration
   # solute 6: FeEDTA adsorbed concentration
   trac
   0.0 1.0 1.e-6 0.5
   1. 2000 1.0 2000
   5 5.0 1.e-6 4.1667e-3
   61
   1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34
   1 202 1 1
   1 202 1 0.
   1 202 101 3.1623e-5 1.0 4.16667
   1
   1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34
   1 202 1 1
   1 202 1 0.
   1 202 101 1.e-13 1.0 4.16667
   1
   1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34
   1 202 1 1
   1 202 1 0.
   1 202 101 3.1623e-5 1.0 4.16667
   0
   1 202 1 0.
   0
   1 202 1 0.
   0
   1 202 1 0.0
   rxn
   ** NCPLX, NUMRXN
   2,4
   ** Coupling of the aqueous components (dRi/dUj)
   2
   1 0 1
   0 1 0
   ** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT
   rows)
   1 Cobalt[aq] 0 0 1.e-9
   2 Iron[aq] 0 0 1.e-9
   3 EDTA[aq] 0 0 1.e-9
   ** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
   101 Co-EDTA[aq] 0
   102 Fe-EDTA[aq] 0
   ** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immoblie spec, NIMM rows)
   1 Co-EDTA[s] 0
   2 Fe-EDTA[s] 0
   3 Cobalt[s] 0
   ** IDVAP(IV), VAPNAM(IM), VAPPRT(IV) (ID # and name of vapor spec, NVAP rows)
   ** Skip nodes
   0
   ** RSDMAX
   1.0e-10
   **** Chemical reaction information for equilibrium reactions ******
   ** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
   0
   ** CKEQ(IX) ,HEQ(IX) (Stability constants and Enthaplys, NCPLX rows)
   1.0e+18 0
   6.31e+27 0
   ** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
   1.0 0.0 1.0
   0.0 1.0 1.0
   ** LINEAR KINETIC REACTION (type 1) **
   1
   ** Where does the reaction take place? **
   1 0 0
   ** Aqueous Component/Complex #, Solid Component #
   101 1
   ** Distribution coeffienct (kg water/ kg rock) **
   0.533
   ** Mass transfer coefficient (1/hr) **
   1.0
   ** LINEAR KINETIC REACTION (type 1) **
   1
   ** Where does the reaction take place? **
   1 0 0
   ** Aqueous Component/Complex #, Solid Component #
   1 3
   ** Distribution coeffienct (kg rock/ kg water) **
   5.07
   ** Mass transfer coefficient (1/hr) **
   1.0
   ** LINEAR KINETIC REACTION (type 1) **
   1
   ** Where does the reaction take place? **
   1 0 0
   ** Aqueous Component/Complex #, Solid Component #
   102 2
   ** Distribution coeffienct (kg rock/ kg water) **
   0.427
   ** Mass transfer coefficient (1/hr) **
   1.0
   ** GENERAL EXCHANGE REACTION (type 3) **
   3
   ** Where does the reaction take place? **
   1 0 0
   ** # of solid, liquid and vapor species **
   3 0 0
   ** forward and reverse rate constants (1/hr) **
   1.26e-2 0
   ** Solid Species in reaction **
   1 2 3
   ** Stoichiometry **
   1.0 -1.0 -1.0
   coor n/a
   202
   1 0.000001.000000.00000
   2 0.100001.000000.00000
   3 0.200001.000000.00000
   ...
   2009.800000.000000.00000
   2019.900000.000000.00000
   20210.000000.000000.00000
   elem
   4 100
   1 10210321
   2 10310432
   3 10410543
   ...
   98 1992009998
   99 20020110099
   100201202101100
   stop


FEHM results for this problem are compared to those of PDREACT (Valocchi et al., 1994), a two-dimensional, isothermal, saturated-zone flow and transport code in `Comparison of FEHM and PDREACT for the breakthrough curves of aqueous species <FEHM-UM.9.5.htm#51118>`_ and `Comparison of FEHM and PDREACT for the exit concentration versus time for solid species <FEHM-UM.9.5.htm#30492>`_. 

.. figure: FEHM-UM.9.5-25.gif
   :caption: 1. Comparison of FEHM and PDREACT for the breakthrough curves of aqueous species.

.. figure: FEHM-UM.9.5-26.gif
   :caption: 2. Comparison of FEHM and PDREACT for the exit concentration versus time for solid species.
 

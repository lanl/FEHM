################
System Interface
################

System-Dependent Features
-------------------------
In addition to standard intrinsic math routines only two system routines are required by the FEHM code. The code uses a system call to get the date (subroutine dated) and a system routine to get the CPU clock time (subroutine tyming).

Compiler Requirements
---------------------
FEHM Version 2.30 is written for Fortran 90. FEHM has been successfully compiled and run on SUN, PC, Linux and macOS computers. 

Hardware Requirements
---------------------
No special hardware features or environments are required by the software. The code will run on SUN workstations running Solaris 7 or later (UNIX), PC workstations running Windows 2000 or later, Linux 2.4.18 or later, and macOS Snow Leopard or later. Memory requirements depend on the problem being modeled (based on the number of nodes). It is suggested that the system being used have a minimum of 128 MB of memory.

Control Sequences or Command Files
----------------------------------
None.

Software Environment
--------------------
N/A

#########################
Installation Instructions
#########################

1. Installing a compiled executable
-----------------------------------
Copy the executable to a location on the current search path. Refer to the Installation Test Plan for the FEHM Application Version 2.30 (10086-ITP-2.30-00).

2. Creating the FEHM binary from source (UNIX)
----------------------------------------------
On the system where FEHM is to be installed, make an installation directory, with subdirectories src and objects:

.. code::

  mkdir fehm
  mkdir src objects

Copy all fehm source files (i.e., extract them from a tar file -- ``fehm_src.tar``) into the src directory:

.. code::

  cd fehm/src
  tar xvf fehm_src.tar

A Makefile is included and should be placed in your objects directory. To compile and link FEHM, change into the objects directory and compile the code:

.. code::

  cd fehm/objects
  make -OR- make -f Makefile

The makefile creates an executable called:

.. code::

  xfehm_v2.30

It should be noted that FEHM uses the GZSOLVE Application (ECD-97) reuse components, solve_new, solve_rdof, and slvesu. The GZSOLVE subroutines are compiled directly into this version of FEHM.

3. Installation Verification and Validation
-------------------------------------------

A series of test scripts have been developed to automate the validation procedure for FEHM. They are described in more detail in the FEHM-VTP.APP, of the Validation Test Plan for the FEHM Application Version 2.30 (10086-VTP-2.30-00). See the FEHM VTP for a discussion of the tests performed and their results. 

4. FEHM for YMP
---------------

For official use of the FEHM code on the YMP project, an executable version should be obtained from the project configuration management group. For binary installation instructions, refer to the Installation Test Plan for the FEHM Application Version 2.30 (10086-ITP-2.30-00).

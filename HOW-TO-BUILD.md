## DMFEHM (Distributed Mesh - FEHM) HOW TO BUILD AND USE

This software and all dependencies are open source software available under the BSD-3 or BSD license. 
This dev-version of FEHM is an attempt to run h2o co2 problems ONLY using Domain Decomposition. 
The method of development was to step through the h2o co2 related functions, and change them to a parallel structure.
To change them to the new structure, while still holding the integrity of the master-FEHM functionallity, every editted
file only contains a call to a new .F90 subroutine which consists of these changes.
These new subroutines are only called if the code is compiled with the -DWITHDD=1 compilation flag.
 
### To obtain the latest version of DMFEHM

This version of the code depends on PETSC. I use gcc-9.7.3 to compile PETSC and FEHM, yet any version of gcc > 6.3.0 should work.
The makefile will automatically download and install PETSC
before building FEHM, as well as setting the correct $PETSC_DIR variable, linking (-L), and including (-I)

Use the pre-procesor directive '-DWTIHDD=0' to compile without Domain Decomposition.

Build this version of FEHM with gcc-7.3.0 or greater.

git clone https://github.com/lanl/FEHM.git --branch DMFEHM 

cd FEHM/src

make -f Makefile.fehm_ubuntu

### To build DMFEHM without petsc (Must set $PETSC_DIR and $PETSC_ARCH correctly, must point to a petsc built with same compiler as FEHM):

git clone https://github.com/lanl/FEHM.git --branch DMFEHM

make -f Makefile.fehm_ubuntu xfehm

### To run DMFEHM:

$PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n NPROCS path/to/build/xfehm

NPROCS = number of processes to use for the solver

ex.
$PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n 16 /scratch/ymp/smckinney/fehm_petsc/FEHM/build/xfehm

### Development Outline

The subroutines listed below show the subroutines that required editting, as well as their call order.
Each subroutine will call a new version of itself labeled as dd_*.f90 if compiled for domain decomposition.

fehmn.f -> dd_init.F90

icetrco2.f -> gensco2h2o.f

gensco2h2o.f ->  geneq_h2o_co2.f -> 

gensco2h2o.f -> ther_co2h2o.f ->

gensco2h2o.f -> co2h2o_combine.f ->

gensco2h2o.f -> solve_new.f ->

gensco2h2o.f -> co2h2o_combine -> gensco2h2o.f (end of loop)

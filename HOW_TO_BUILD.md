## FEHM-PETSC_32 HOW TO BUILD AND USE

This software and all dependencies are open source software available under the BSD-3 or BSD license. 
The 32bit version of FEHM-PETSC can handle ~7mil nodes and runs using petsc's parallel solver.
 
### To obtain the latest version of FEHM-PETSC_32

Build this version of FEHM with gcc-7.3.0 or greater.

git clone https://github.com/lanl/FEHM.git --branch PETSC  
cd FEHM/src
make -f Makefile.fehm_ubuntu

### To build FEHM-PETSC without petsc (Must set $PETSC_DIR and $PETSC_ARCH correctly, must point to a petsc built with same compiler as FEHM):

git clone https://github.com/lanl/FEHM.git --branch PETSC
make -f Makefile.fehm_ubuntu xfehm

### To run FEHM-PETSC_32:

$PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n NPROCS path/to/build/xfehm

NPROCS = number of processes to use for the solver

ex.
$PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n 16 /scratch/ymp/smckinney/fehm_petsc/FEHM/build/xfehm

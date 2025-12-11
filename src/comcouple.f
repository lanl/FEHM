      module comcouple
!     comcouple.h
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1
!D1  PURPOSE
!D1  
!D1 Include file containing new parameters related to the solute module.
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY 
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comcouple.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:24   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.1   Tue Jan 16 15:55:12 1996   zvd
!D2 Added prolog.
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  N/A
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4  
!D4  None
!D4  
!***********************************************************************
!D4 Global Variables
!D4
!D4                            COMMON
!D4 Identifier      Type     Block  Description
!D4
!D4 ngroups          int     comrxni.h  The number of groups of 
!D4                                      coupled species
!D4 group            int     comrxni.h  Array indicating the coupling
!D4                                     between species
!D4 pos              int     comrxni.h  Array containing the species
!D4                                     numbers of the coupled species
!D4                                     in each group
!D4 n_couple_species int     comrxni.h  Array containg the number of
!D4                                     coupled species contained in
!D4                                     each grouping
!D4 mdof_sol         int     comrxni.h  The maximum number of degrees
!D4                                     of freedom necessary to solve
!D4                                     the tracer solution
!D4 nb_sol           int    davidi.h    The array for nb in the 
!D4                                     solute transport solution
!D4 nmat_sol         int    davidi.h    Array denoting the positions
!D4                                     of each submatrix in the 
!D4                                     Jacobian for the solute 
!D4                                     transport solution
!D4 nrhs_sol         int    davidi.h    Array denoting the positions
!D4                                     in the residual array for
!D4                                     each degree of freedom
!***********************************************************************

      integer ngroups
      integer, allocatable :: group(:,:)
      integer, allocatable :: pos(:,:)
      integer, allocatable :: n_couple_species(:)
      integer mdof_sol
      integer, allocatable :: nb_sol(:)
      integer, allocatable :: nmat_sol(:)
      integer, allocatable :: nrhs_sol(:)
c zvd 02-Sep-08 Make henry model variables allocatable
      integer, allocatable :: henry_model(:)
c      real*8, allocatable ::  awwa(:,:)
      real*8, allocatable ::  hawwa(:,:)
      real*8, allocatable ::  fzero(:)
      
      end module comcouple

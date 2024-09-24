      subroutine den_vis_spatial(iflg)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To enable a spatially variable density and viscosity that does not change in time
CD1 gaz 121921 this routine reads in density, viscosity, and compressibility of water
CD1 used in isothermal applications, need to make more automatic      
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 012014            G. Zyvoloski   N/A     Initial implementation
CD2
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 iieosd       int     I       Flag denoting if simple
CD3                                  thermodynamics are to be used
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 inpt                   I    File contains input data
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 None
CD4 
CD4
C**********************************************************************
      use comai
      use combi
      use comci
      use comdi
      use comii
      use comdi
      use comdti
      use comki
      implicit none

      integer i,icode
      integer iflg

       if ( iflg .eq. 0 )  then
        macro = 'den '
      
        allocate(den_spatial(n0),vis_spatial(n0),comp_spatial(n0))
        allocate(deng_spatial(n0),visg_spatial(n0))
        den_spatial = 0.0d0
        comp_spatial = 0.0d0
        deng_spatial = 0.0d0
        vis_spatial = 0.0d0
        visg_spatial = 0.0d0
	      
c**** read den, vis for liquids ****
        narrays = 3
        itype(1) = 8
        itype(2) = 8
        itype(3) = 8
        default(1) = 0.
        default(2) = 0.
        default(3) = 0.
        igroup = 1
c gaz 110418 cmacroread(12) to macroread(24)      
        call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(24), macro, igroup, ireturn,
     *     r8_1=den_spatial(1:n0),r8_2=vis_spatial(1:n0),
     *     r8_3 =comp_spatial(1:n0))

c**** read den, vis for gases ****
        narrays = 2
        itype(1) = 8
        itype(2) = 8
        default(1) = 0.
        default(2) = 0.
        igroup = 1
      
        call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(24), macro, igroup, ireturn,
     *     r8_1=deng_spatial(1:n0),r8_2=visg_spatial(1:n0))
c
c  check if values are physical
c
        iden_vis = 0
        ideng_vis = 0
        do i = 1, n0
         if(den_spatial(i).gt.0.0d00.or.vis_spatial(i).gt.0.0d00) 
     *      iden_vis = 1
         if(deng_spatial(i).gt.0.0d00.or.visg_spatial(i).gt.0.0d00) 
     *       ideng_vis = 1
        enddo
      endif
      return
      end

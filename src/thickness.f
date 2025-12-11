      subroutine thickness(iflg)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine modifies volumes and finite element coefficients
CD1  to account for variable thickness.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD 
CD2 Date         Programmer     Number  Comments 
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/thickness.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Fri Jan 12 18:01:24 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.3   Thu Jan 11 11:57:18 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   04/25/95 08:47:34   llt
CD2 retrieved lost log history
CD2 
CD2    Rev 1.1   04/25/95 08:15:12   llt
CD2 added log history information to prolog
CD2
CD2    Rev 1.0   01/28/95 13:56:06   llt
CD2 water balance equation was modified
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

      use combi
      use comdi
      use comgi
      use comrxni
      use comfi
      use comdti
      use comai
      implicit none

      integer iflg,icode,i,j,i2,kb,isx,ipiv
      real*8 thick_avg

c================================================================
c iflg       - parameter for subroutine operations(input,etc)
c neq        - number of nodes
c ndim_sx    - first dimension of coefficient matrix
c i          - do loop index
c j          - do loop index
c kb         - neighbor node 
c isx        - index for coefficient matrix
c ipiv       - do loop index 
c num_neg    - number of negative volumes or coefficients
c iptty      - unit number for output file
c max_num    - maximun number of negative volume or coefficients
c            - before the program stops
c sx1        - volume array
c sx         - flow coefficient array
c sumsx      - intermediate sum of coefficients
c thick      - thickness of reservoir at every node point
c thick_avg  - average thickness between node and its neighbor
c================================================================

      if(iflg.eq.0) then
         allocate(thic(neq))
         
c     read in thickness information for nodes
         call rdthick
      else if(iflg.ne.0.and.ithic.ne.0) then
c     first modify volumes
         do i=1,neq
            sx1(i)=sx1(i)*thic(i)
         enddo
c     next modify permeabilities, assume isotropic conditions)
         do i=1,neq
            cord(i,3)=cord(i,3)*thic(i)
         enddo
         deallocate(thic)
      endif

      return
      end                

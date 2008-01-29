      subroutine indpdp
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
CD1 To modify fracture volume at nodes for dpdp calculations.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/indpdp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Tue Jan 30 11:39:54 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   06/23/95 14:14:34   robinson
CD2 Corrected bug in setting of the fracture volume fraction
CD2 
CD2    Rev 1.1   03/18/94 15:47:36   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:54   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use      Description
CD3 
CD3 None
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier   Type        Description
CD4 
CD4 
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 Identifier   Type        Description
CD4 
CD4 neq, irlp, irlpt, volf1, rp19f, sx1, cord, pnx, pny, pnz, thx,
CD4 thy, thz, nmat, nrhs, nelm, ico2, irdof
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type       Description
CD4 
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop index
CD5 sx1d         real*8      Volume of current node
CD5 neqp1        int         Number of equations plus 1
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 N/A
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.4.9 Double-porosity/double-permeability formulation
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN indpdp
CPS 
CPS FOR each node
CPS   IF the van Genutchen relative permeability model is used
CPS     Set volume fraction for fracture
CPS   ENDIF
CPS   Compute volumes for fracture and matrix nodes
CPS   Set coordinate value in array
CPS ENDFOR
CPS 
CPS FOR each node
CPS   Compute permeability and thermal conductivity values for...
CPS   ... fracture and matrix nodes
CPS ENDFOR
CPS 
CPS Compute values in pointer array
CPS 
CPS IF there is no noncondensible gas present
CPS   FOR each pointer value
CPS     Set the value in the pointer array
CPS   ENDFOR
CPS ELSE there is a noncondensible gas
CPS   IF the solution is without reduced degrees of freedom
CPS     FOR each pointer value
CPS       Set the value in the pointer array
CPS     ENDFOR
CPS   ELSE a reduced degree of freedom model is used
CPS     Set the value in the pointer array
CPS   ENDIF
CPS ENDIF
CPS 
CPS Set values in pointer array for right hand side (solution vector)
CPS 
CPS END indpdp
CPS 
C**********************************************************************

      use davidi
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer i,neqp1
      real*8 sx1d
      logical, save :: first_time = .true.
      if(first_time) then
      do i=1,neq
c     check if van Genutchen rl perms were selected
         if(irlpt(irlp(i)).eq.4.or.irlpt(irlp(i)).eq.7) then
c     calculate volume factor from fracture porosity
            volf1(i)=abs(rp19f(irlp(i)))
         endif
         sx1d=sx1(i)
         sx1(i) =  sx1d*volf1(i)
         sx1(i+neq) =  sx1d*(1.-volf1(i))
         cord(i+neq,3) =  cord(i,3)
      enddo
      end if
      first_time = .false.
c
c**** weight transport properties
c
      do i=1,neq
         pnx(i+neq)=pnx(i+neq)*(1.0-volf1(i))
         pny(i+neq)=pny(i+neq)*(1.0-volf1(i))
         pnz(i+neq)=pnz(i+neq)*(1.0-volf1(i))
         pnx(i)=pnx(i)*volf1(i)
         pny(i)=pny(i)*volf1(i)
         pnz(i)=pnz(i)*volf1(i)
         if(ico2.ge.0) then
           thx(i+neq)=thx(i+neq)*(1.0-volf1(i))
           thy(i+neq)=thy(i+neq)*(1.0-volf1(i))
           thz(i+neq)=thz(i+neq)*(1.0-volf1(i))
           thx(i)=thx(i)*volf1(i)
           thy(i)=thy(i)*volf1(i)
           thz(i)=thz(i)*volf1(i)
         endif
      enddo
c     
c     set up pointer array for the A matrix
c     
      neqp1=neq+1
      nmat(1)=0
      nmat(2)=nelm(neqp1)-neqp1
      if(ico2.le.0) then
c     no non condensible present
         do i=3,16
            nmat(i)=nmat(i-1)+nmat(2)
         enddo
      else
c     we have a non condensible
c     pointer length depends on irdof for this application
         if(irdof.eq.0) then
            do i=3,36
               nmat(i)=nmat(i-1)+nmat(2)
            enddo
         else
c     notice the connection terms are neq long
            nmat(3)=nmat(2)+nmat(2)
            nmat(4)=nmat(3)+nmat(2)
            nmat(5)=nmat(4)+neq
            nmat(6)=nmat(5)+neq
            nmat(7)=nmat(6)+neq
            nmat(8)=nmat(7)+nmat(2)
            nmat(9)=nmat(8)+nmat(2)
            nmat(10)=nmat(9)+nmat(2)
            nmat(11)=nmat(10)+neq
            nmat(12)=nmat(11)+neq
            nmat(13)=nmat(12)+neq
            nmat(14)=nmat(13)+nmat(2)
            nmat(15)=nmat(14)+nmat(2)
            nmat(16)=nmat(15)+nmat(2)
            nmat(17)=nmat(16)+neq
            nmat(18)=nmat(17)+neq
            nmat(19)=nmat(18)+neq
            nmat(20)=nmat(19)+neq
            nmat(21)=nmat(20)+neq
            nmat(22)=nmat(21)+neq
            nmat(23)=nmat(22)+nmat(2)
            nmat(24)=nmat(23)+nmat(2)
            nmat(25)=nmat(24)+nmat(2)
            nmat(26)=nmat(25)+neq
            nmat(27)=nmat(26)+neq
            nmat(28)=nmat(27)+neq
            nmat(29)=nmat(28)+nmat(2)
            nmat(30)=nmat(29)+nmat(2)
            nmat(31)=nmat(30)+nmat(2)
            nmat(32)=nmat(31)+neq
            nmat(33)=nmat(32)+neq
            nmat(34)=nmat(33)+neq
            nmat(35)=nmat(34)+nmat(2)
            nmat(36)=nmat(35)+nmat(2)
         endif
      endif

c     set up pointer array for RHS
      nrhs(1)=0
      nrhs(2)=neq
      nrhs(3)=neq+neq
      nrhs(4)=neq+neq+neq
      nrhs(5)=neq+neq+neq+neq
      nrhs(6)=neq+neq+neq+neq+neq

      return
      end

      subroutine dualtx
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
CD1 To backsubstitute to get tracer solution for dual porosity nodes.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-4-93      G. Zyvoloski   N/A     Initial implementation, but
CD2                                        previous non-YMP versions
CD2                                        of FEHM exist, and the
CD2                                        current version may differ
CD2                                        from these
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dualtx.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 15:45:10 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:50:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:36   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
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
CD4 neq, bp, a21mpf, wb11, rb2mf, tb11, r3mf, a32mpf, an, npn
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
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4 
CD4 
CD4 Global Subprograms
CD4 
CD4 None
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
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 id           int         Do loop index parameter over each node
CD5 idp          int         Number of first matrix node
CD5 idpp         int         Number of second matrix node
CD5 x1           real*8      Correction of current primary unknown
CD5 x2           real*8      Correction of first matrix node
CD5 x3           real*8      Correction of second matrix node
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
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
CD9 2.3.4 Solute-transport equations
CD9 2.4.7 Dual-porosity formulation
CD9 
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN dualtx
CPS 
CPS FOR each node
CPS   Calculate indexes for matrix nodes
CPS   Calculate corrections to concentrations
CPS   Calculate concentrations
CPS ENDFOR each node
CPS 
CPS END dualtx
CPS 
C**********************************************************************

      use combi
      use comdi
      use comei
      use comfi
      use comgi
      use comhi
      use comdti
      use comai
      implicit none
      
      integer id,idp,idpp
      real*8 x1,x2,x3
      
      do id=1,neq
         idp=id+neq
         idpp=id+neq+neq
         
c     calculate node idp solution
c     x2=ab22i*(rb2-a21*x1)
c     find rt=rb2-a21*x1
c     
         x1=bp(id)
         x2=wb11(id)*(rb2mf(id)-a21mpf(id)*x1)
         
c     calculate node idpp solution
c     x3=a33i*(r3-a32*x2)
c     find rt=r3-a32*x2
c     
         x3=tb11(id)*(r3mf(id)-a32mpf(id)*x2)
c     
c     update tracer varible set
c     n-r corrections
         
         an(idp+npn)=an(idp+npn)-x2
         an(idpp+npn)=an(idpp+npn)-x3
      enddo
      
      return
      end

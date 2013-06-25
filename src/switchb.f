      subroutine switchb (bp, nrhs, nrhsb, islord, idof, iswitch, neq)
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To reorder "b" matrix.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/switchb.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jan 17 15:33:10 1996   zvd
CD2 Added to prolog
CD2 
CD2    Rev 1.4   Thu Jan 11 11:49:28 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   Tue Jan 09 14:16:18 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.2   11/15/95 15:26:34   gaz
CD2 fixed bug
CD2 
CD2    Rev 1.1   03/18/94 15:59:56   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:30   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.5.2 Solve nonlinear equation set at each time step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4  This is a general utility routine used in the code whenever 
CD4  necessary.
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5   
CD5   bp              REAL*8   O    Array of Newton-Raphson residuals, after
CD5                                   solution an array of Newton-Raphson
CD5                                   corrections
CD5   idof            INT      I    Number of degrees of freedom per node for
CD5                                   the current problem
CD5   islord          INT      I    Parameter used in the reduced degree of
CD5                                   freedom model
CD5   iswitch         INT
CD5   neq             INT      I    Number of nodes, not including dual
CD5                                   porosity nodes
CD5   nrhs            INT           Array used in the reduced degree of
CD5                                      freedom method
CD5   nrhsb           INT           Array used in the reduced degree of
CD5                                      freedom method
CD5   
CD5 Interface Tables
CD5
CD5   None
CD5
CD5 Files
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6
CD6   None
CD6
CD6 Global Subprograms
CD6
CD6   None
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7   None
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   bpdum2          REAL*8   Temporary storage for bp value when being
CD7                              switched
CD7   bpdum3          REAL*8   Temporary storage for bp value when being
CD7                              switched
CD7   i               INT      Loop index
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN switchb
CPS 
CPS   IF reordering should be done
CPS  
CPS      IF the problem has three degrees of freedom
CPS  
CPS         IF switch 2
CPS            FOR each node
CPS                set residual values
CPS            END FOR
CPS         END IF
CPS  
CPS      ELSE IF the problem has two degrees of freedom
CPS  
CPS         IF switch is 1 and reordering parameter is -1
CPS            set first and second rdof array value from b rdof array
CPS         END IF
CPS         IF switch is 2
CPS            set first rdof array value to zero
CPS            set second rdof array value to neq
CPS         END IF
CPS  
CPS      ELSE IF the problem has four degrees of freedom
CPS  
CPS         IF switch 1
CPS            set first to fourth rdof array value from b rdof array
CPS         END IF
CPS         IF switch 2
CPS            set first rdof array value to zero
CPS            set second rdof array value to neq
CPS            set third rdof array value to 2*neq
CPS            set fourth rdof array value to 3*neq
CPS         END IF
CPS  
CPS      ELSE IF the problem has six degrees of freedom
CPS  
CPS         IF switch 1
CPS            set first to fourth rdof array value from b rdof array
CPS         END IF
CPS         IF switch 2
CPS            set first rdof array value to zero
CPS            set second rdof array value to neq
CPS            set third rdof array value to 2*neq
CPS            set fourth rdof array value to 3*neq
CPS            set fifth rdof array value to 4*neq
CPS            set sixth rdof array value to 5*neq
CPS         END IF
CPS  
CPS      END IF
CPS  
CPS   END IF
CPS   
CPS END switchb
CPS
C***********************************************************************

      implicit none

      integer i, idof, islord, iswitch, neq, nrhs(*), nrhsb(*)
      real*8 bp(*), bpdum1, bpdum2, bpdum3

*     Change if test to eliminate return at start of routine
      if (islord .ne. 0) then
         if (idof .eq. 3) then
            if (iswitch .eq. 2) then
               do i=1,neq
                  bpdum1=bp(i+nrhs(1))
                  bpdum2=bp(i+nrhs(2))
                  bpdum3=bp(i+nrhs(3))
                  bp(i+nrhsb(1))=bpdum1
                  bp(i+nrhsb(2))=bpdum2
                  bp(i+nrhsb(3))=bpdum3
               enddo
            else
c no change:equations ordered naturally, variables arranged diff               
            end if
         else if (idof .eq. 2) then
            if (iswitch .eq. 1 .and. islord .eq. -1) then
               nrhs(1)=nrhsb(1)
               nrhs(2)=nrhsb(2)
            end if
            if (iswitch .eq. 2) then
               nrhs(1)=0         
               nrhs(2)=neq         
            end if
         else if (idof .eq. 4) then
            if (iswitch .eq. 1) then
               nrhs(1)=nrhsb(1)
               nrhs(2)=nrhsb(2)
               nrhs(3)=nrhsb(3)
               nrhs(4)=nrhsb(4)
            end if
            if (iswitch .eq. 2) then
               nrhs(1)=0         
               nrhs(2)=neq         
               nrhs(3)=2*neq
               nrhs(4)=3*neq      
            end if
         else if (idof .eq. 6) then
            if (iswitch .eq. 1) then
               nrhs(1)=nrhsb(1)
               nrhs(2)=nrhsb(2)
               nrhs(3)=nrhsb(3)
               nrhs(4)=nrhsb(4)
               nrhs(5)=nrhsb(5)
               nrhs(6)=nrhsb(6)
            endif
            if (iswitch .eq. 2) then
               nrhs(1)=0         
               nrhs(2)=neq         
               nrhs(3)=2*neq
               nrhs(4)=3*neq      
               nrhs(5)=4*neq
               nrhs(6)=5*neq      
            endif
         end if
      end if

      return
      end

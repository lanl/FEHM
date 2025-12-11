      subroutine switch (nmat, nmatb, islord, idof, iswitch, nsizea)
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
CD1 To reorder the "a" matrix.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/switch.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jan 17 15:33:00 1996   zvd
CD2 Added to prolog
CD2 
CD2    Rev 1.4   Thu Jan 11 11:44:56 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   Tue Jan 09 14:15:44 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.2   11/15/95 15:25:58   gaz
CD2 corrected minor bug
CD2 
CD2    Rev 1.1   03/18/94 15:59:54   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:28   pvcs
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
CD5   idof            INT      I    Number of degrees of freedom per node for
CD5                                   the current problem
CD5   islord          INT      I    
CD5   iswitch         INT      I    
CD5   nmat            INT      O    Array used in the reduced degree of
CD5                                      freedom method      
CD5   nmatb           INT      I    Array used in the reduced degree of
CD5                                      freedom method
CD5   nsizea          INT      I    
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
CPS BEGIN switch
CPS 
CPS   IF reordering should be done 
CPS   
CPS      IF the problem has three degrees of freedom
CPS      
CPS         IF reordering parameter is 1
CPS         
CPS            IF switch 1
CPS               FOR each position
CPS                   set rdof array value
CPS               END FOR
CPS            ELSE IF switch 2
CPS               FOR each  position
CPS                   calculate rdof array value
CPS               END FOR
CPS            END IF
CPS            
CPS         END IF
CPS         
CPS      ELSE IF the problem has two degrees of freedom
CPS      
CPS         IF switch 1
CPS            FOR each  position
CPS                set rdof array value
CPS            END FOR
CPS         ELSE IF switch 2
CPS            FOR each  position
CPS                calculate rdof array value
CPS            END FOR
CPS         END IF
CPS      
CPS      ELSE IF the problem has four degrees of freedom
CPS      
CPS         IF switch 1
CPS            FOR each  position
CPS                set rdof array value
CPS            END FOR
CPS         ELSE IF switch 2
CPS            FOR each  position
CPS                calculate rdof array value
CPS            END FOR
CPS         END IF
CPS      
CPS      ELSE IF the problem has six degrees of freedom
CPS      
CPS         IF switch 1
CPS            FOR each  position
CPS                set rdof array value
CPS            END FOR
CPS         ELSE IF switch 2
CPS            FOR each  position
CPS                calculate rdof array value
CPS            END FOR
CPS         END IF
CPS      
CPS      END IF
CPS
CPS   END IF
CPS 
CPS END switch
CPS
C***********************************************************************

      implicit none

      integer i,idof,islord,iswitch,nmatb(*),nmat(*),nsizea

* Change if test to eliminate return at start of routine
      if (islord .ne. 0) then
         if (idof .eq. 3) then
            if (islord .ne. 0) then
               if (iswitch .eq. 1) then
                  do i = 1, 9
                     nmat(i) = nmatb(i)
                  end do
               else if (iswitch .eq. 2) then
                  do i = 1, 9
                     nmat(i) = (i - 1) * nsizea
                  end do
               end if
            end if
         else if (idof .eq. 2) then
            if (iswitch .eq. 1) then
               do i = 1, 4
                  nmat(i) = nmatb(i)
               end do
            else if (iswitch .eq. 2) then
               do i = 1, 4
                  nmat(i) = (i - 1) * nsizea
               end do
            end if
         else if (idof .eq. 4) then
            if (iswitch .eq. 1) then
               do i = 1, 16
                  nmat(i) = nmatb(i)
               end do
            else if (iswitch .eq. 2) then
               do i = 1, 16
                  nmat(i) = (i - 1) * nsizea
               end do
            end if
         else if (idof .eq. 6) then
            if (iswitch .eq. 1) then
               do i = 1, 36
                  nmat(i) = nmatb(i)
               enddo
            else if (iswitch .eq. 2) then
               do i = 1, 36
                  nmat(i) = (i - 1) * nsizea
               enddo
            endif
         end if 
      end if 
      
      return
      end

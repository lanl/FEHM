      subroutine inflo2
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
CD1 Read flow data input by planes for 3D models.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 21-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/inflo2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:04 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:17:54 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:34:44 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Fri Feb 16 08:59:44 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.2   Tue Jan 30 11:49:02 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:03:28   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:58   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   esk             REAL*8   fdd    Inlet enthalpy associated with a source
CD4   inpt            INT      faai   Unit number for input file
CD4   ka              INT      fbb    Contains boundary type information for
CD4                                     each node
CD4   pflow           REAL*8   fdd    Flowing pressure at each source node
CD4   sk              REAL*8   fdd    Source strength of each node
CD4   wdd1            CHAR     faac   Alternate character input string
CD4   wellim          REAL*8   fdd    Well impedance at each source nodeCD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   null1            LOGICAL  Check for null lines or 0's in lines
CD4 
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   aiped           REAL*8   Temporary variable for impedance parameter
CD5   ij              INT      Loop index
CD5   ja              INT      Loop index
CD5   jb              INT      Loop index
CD5   jc              INT      Loop index
CD5   jd              INT      Loop index
CD5   jk              INT      Loop index
CD5   kl              INT      Loop index
CD5   eskd            REAL*8   Temporary variable for enthalpy
CD5   skd             REAL*8   Temporary variable for source strength
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.3.7 Sources and sinks
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN inflo2
CPS 
CPS   REPEAT
CPS     
CPS     read input line
CPS     
CPS     EXIT IF null line is found
CPS     
CPS     reread input line using flo2 format
CPS     
CPS     IF first input index is zero exit loop
CPS     
CPS     FOR each node in the input range
CPS         calculate index increment
CPS         FOR each node in the calculated range
CPS             set the enthalpy value
CPS             IF the impedance variable is zero
CPS                set the flow and boundary type to 1
CPS             ELSE
CPS                IF the second temporary variable is negative
CPS                   set the boundary type to -2
CPS                ELSE
CPS                   set the boundary type to -1
CPS                END IF
CPS                set the impedance and flowing pressure
CPS             END IF
CPS         END FOR
CPS     END FOR
CPS   
CPS   UNTIL end of flo2 data is found
CPS
CPS END inflo2
CPS
C***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      implicit none
      
      logical null1
      integer ij, ja, jb, jc, jd, jk, kl
      real*8  aiped, eskd, skd
      
 100  continue
      
      read (inpt, '(a80)') wdd1
      if (null1(wdd1)) go to 105 
      backspace inpt
      read(inpt,*) ja, jb, jc, jd, skd, eskd, aiped
      if (ja .eq. 0) go to 105
      jd = max0(jd, 1)
      do jk = ja, jb
         kl = jk - ja
         do ij = ja + kl, jc + kl, jd
            esk(ij) = eskd
            if(esk(ij).gt.0.and.igrav.ne.0) then
             esk(ij) = eskd-grav*cord(ij,igrav)
            endif
            if (abs(aiped) .lt. zero_t) then
               sk(ij) = skd
            else
               if (aiped .lt. 0.0) then
                  ka(ij) = -2
               else
                  ka(ij) = -1
               end if
               wellim(ij) = abs(aiped) * 1.0d+06
               pflow (ij) = skd
            end if
         end do
      end do
      goto  100
 105  continue
      
      end

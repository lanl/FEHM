      subroutine  contrc
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
CD1 To write solute information to the contour file.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/contrc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Mon Jan 29 14:00:30 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   03/18/94 16:15:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/14/94 11:49:34   zvd
CD2 Modified so if contour files are not present a write won't be attempted.
CD2 
CD2    Rev 1.0   01/20/94 10:22:32   pvcs
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
CD3 Name                  Description
CD3 
CD3 File number iscon     Contour plot file
CD3 File number iscon1    Contour plot file for dual or dpdp runs
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 None
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 iscon, altc, nspeci, npn, npt, icns, neq, anv, an, idualp, idpdp,
CD4 iscon1
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
CD5 Identifier   Type        Description
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 iq           int         Do loop index
CD5 iu           int         Do loop index
CD5 idpflag      int         Flag denoting whether dp nodes are to be
CD5                              written
CD5 ifirst       int         First position in array of dp values
CD5 ilast        int         Last position in array of dp values
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following functions are carried out in this routine:
CD6 
CD6   An IF check determines if contour information is to be written,
CD6   and if it is:
CD6   
CD6     The code checks to determine if formatted or unformatted
CD6     writing is to be done, and does the writing accordingly.
CD6     
CD6     To do the writing, the code loops through each solute, and
CD6     writes the entire array of concentrations.
CD6     
CD6   The code then returns to the calling routine.
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
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN contrc
CPS 
CPS   IF concentrations are to be written (output file exists)
CPS   
CPS     IF a formatted write is to be performed
CPS     
CPS        FOR each solute
CPS            IF this is a Henry's Law vapor species
CPS               Write the array of concentrations
CPS            ELSE
CPS               Write the array of concentrations
CPS            ENDIF
CPS        ENDFOR
CPS     
CPS     ELSE an unformatted write is to be performed
CPS     
CPS        FOR each solute
CPS            IF this is a Henry's Law vapor species
CPS               Write the array of concentrations
CPS            ELSE
CPS               Write the array of concentrations
CPS            ENDIF
CPS        ENDFOR
CPS     
CPS     ENDIF
CPS   
CPS   ENDIF
CPS   
CPS   IF dpdp/dual concentrations are to be written (dpdp/dual output . . .
CPS   . . . file exists)
CPS   
CPS      Initialize flag denoting whether dp nodes are to be written
CPS       
CPS      IF this is a dpdp solution
CPS         Set limits for array values to be written for dp nodes
CPS         Set flag to indicate dp nodes are to be written
CPS      ELSEIF it is a dual porosity solution
CPS         Set limits for array values to be written for dp nodes
CPS         Set flag to indicate dp nodes are to be written
CPS      ENDIF
CPS       
CPS      IF dp nodes are to be written
CPS         IF a formatted write is to be performed
CPS     
CPS            FOR each solute
CPS                IF this is a Henry's Law vapor species
CPS                   Write the array of concentrations
CPS                ELSE
CPS                   Write the array of concentrations
CPS                ENDIF
CPS            ENDFOR
CPS     
CPS         ELSE
CPS            
CPS            FOR each solute
CPS                IF this is a Henry's Law vapor species
CPS                   Write the array of concentrations
CPS                ELSE
CPS                   Write the array of concentrations
CPS                ENDIF
CPS            ENDFOR
CPS     
CPS         ENDIF
CPS         
CPS      ENDIF
CPS     
CPS   ENDIF
CPS     
CPS END contrc
CPS 
C***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      use comxi
      implicit none

      integer iq,iu,ifirst,ilast,idpflag

      if ( iscon .gt. 0 )  then
         if(altc .ne. 'fehm') then
            do iq = 1, nspeci
               npn = npt(iq)
               if(icns(iq) .eq. -2) then
                  write(iscon, 6000) (anv(iu+npn), iu = 1, neq)
               else
                  write(iscon, 6000) (an(iu+npn), iu = 1, neq)
               endif
            end do
         else
            do iq = 1, nspeci
               npn = npt(iq)
               if(icns(iq) .eq. -2) then
                  write(iscon) (anv(iu+npn), iu = 1, neq)
               else
                  write(iscon) (an(iu+npn), iu = 1, neq)
               endif
            end do
         endif
      endif

      if (iscon1 .gt. 0) then
         idpflag = 0
         if (idpdp .ne. 0 ) then
            ifirst = neq + 1
            ilast = 2 * neq
            idpflag = 1
         else if(idualp .ne. 0) then
            ifirst = 2 * neq + 1
            ilast = 3 * neq
            idpflag = 1
         end if

         if( idpflag .eq. 1 ) then
            if(altc .ne. 'fehm') then
               do iq = 1, nspeci
                  npn = npt(iq)
                  if (icns(iq) .eq. -2) then
                     write(iscon1 ,6000)
     2                    (anv(iu+npn), iu = ifirst, ilast)
                  else
                     write(iscon1 ,6000)
     2                    (an(iu+npn), iu = ifirst, ilast)
                  end if
               end do
            else
               do iq = 1, nspeci
                  npn = npt(iq)
                  if (icns(iq) .eq .-2) then
                     write(iscon1)
     2                    (anv(iu+npn), iu = ifirst, ilast)
                  else
                     write(iscon1)
     2                    (an(iu+npn), iu = ifirst, ilast)
                  end if
               end do
            end if
         end if
      end if

 6000 format(1x,5f12.8)
      end

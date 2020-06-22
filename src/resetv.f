      subroutine  resetv( ndum )
!***********************************************************************
!  Copyright, 1994, 2004,  The  Regents of the University of California.
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
CD1 To reset variables to old time step values.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/resetv.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:20   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Fri Mar 01 14:21:34 1996   gaz
CD2 changed ( ico2 .gt. 0 ) to  ( ico2 .ne. 0 )  
CD2 
CD2    Rev 1.3   Wed Jan 17 13:00:02 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.2   Thu Jan 11 10:20:14 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:43:20   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:14   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.5.1 Implement time-step mechanism
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
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
CD5 Identifier    Type   Use    Definition
CD5
CD5 ndum          INT     I     Index parameter for dual porosity option
CD5
CD5 Interface Tables
CD5
CD5 None
CD5
CD5 Files
CD5 
CD5 None
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
CD6 phi, pho, s, so, t, to, ieos, ico2, neq, pci, pcio, ice, sii, sio,
CD6 ices
CD6
CD6 Global Subprograms
CD6
CD6 None
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
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   i               INT      Node number
CD7   id              INT      Do loop index
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN resetv
CPS
CPS FOR each node
CPS   Compute unknown number
CPS   Reset pressure, saturation, and temperature to previous value
CPS ENDFOR
CPS
CPS FOR each node
CPS   Compute unknown number
CPS   IF the old fluid state was vapor
CPS     Reset fluid state variable to vapor
CPS   ELSEIF the old fluid state was liquid
CPS     Reset fluid state variable to liquid
CPS   ELSE the fluid was 2 phase
CPS     Reset the fluid state to 2 phase
CPS   ENDIF
CPS
CPS ENDFOR
CPS
CPS IF this is a noncondensible gas simulation
CPS   FOR each node
CPS     Compute unknown number
CPS     Reset capillary pressure to the old value
CPS   ENDFOR
CPS ENDIF
CPS 
CPS IF this is an ice solution
CPS
CPS   FOR each node
CPS     Compute unknown number
CPS     Reset ice state parameter to old value
CPS   ENDFOR
CPS
CPS   FOR each node
CPS
CPS     Compute unknown number
CPS     IF the old ice state was no ice
CPS       Reset ice state variable to no ice
CPS     ELSEIF the old ice state was fully frozen
CPS       Reset ice state variable to fully frozen
CPS     ELSE the system was ice and water
CPS       Reset the ice state to ice and water
CPS     ENDIF
CPS     
CPS   ENDFOR
CPS ENDIF
CPS 
CPS END resetv
CPS
C***********************************************************************

      use comci
      use comdi
      use comfi
      use comdti
      use comai
      use comco2
      use davidi
      implicit none
c gaz 121718   0923219 h2o_crit   moved to comai
c      real*8 pcrit_h2o, tcrit_h2o
c      parameter(pcrit_h2o=22.00d0, tcrit_h2o=373.95)
      integer ndum,id,i
c first check for saturated only

      if(irdof.eq.13) then
         do id=1,neq
            i      =  id+ndum
            phi(i) =  pho(i)
            if (ifree .ne. 0) then
c               s  (i) =  so (i)
                rlxyf(i) =  so (i)
            end if
            deni(i) = 0.0
         enddo
         return
      else
         do id = 1, neq
            i      =  id+ndum
            phi(i) =  pho(i)
            s  (i) =  so (i)
            t  (i) =  to (i)
            deni(i) = 0.0
            denei(i) = 0.0
         end do
      endif
c gaz 121718 modifications to support SC water
      
c**** reset noncondensible gas if appropriate ****
      if ( ico2 .gt. 0 )  then
         do id = 1, neq
            i      =  id+ndum
            pci(i) =  pcio(i)
         end do
      end if
c**** reset fluid state ****
      do id = 1, neq
         i      =  id+ndum
         if(phi(i).ge.pcrit_h2o.and.t(i).ge.tcrit_h2o) then
            ieos(i) = 4
c gaz 102319 insure that pcp = 0.0 and dpcef = 0.0 for SC
            pcp(i) = 0.0
            dpcef(i) = 0.0
         else if ( so(i) .le. 0.0 )  then
            ieos(i) =  3
         else if ( so(i) .ge. 1.0 )  then
            ieos(i) =  1
         else
            ieos(i) =  2
         end if
c gaz 071919 add correction for rich equation (no ieos(i) = 3)
       if(jswitch.ne.0.and.ieos(i).eq.3) ieos(i) = 2
      end do



c**** reset ice solution ****
      if ( ice .ne. 0 )  then
         call icectr(11,0)
      end if
c**** RJP 04/10/07 reset co2 solution ****
      if ( icarb .eq. 1 )  then
         call icectrco2(11,0)
      end if
c**** reset stress solution ****
      if ( istrs .ne. 0 )  then
         call stressctr(-10,0)
      end if
c**** reset time_ieos
      time_ieos = 0.0d0      
      return
      end

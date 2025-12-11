      subroutine  gendat  ( mb,mc,id )
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
CD1 To generate coordinates and element information in simple 
CD1 geometric problems.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gendat.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Thu Feb 15 10:31:12 1996   zvd
CD2 Updated purpose.
CD2 
CD2    Rev 1.4   Wed Jan 17 12:24:52 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.3   Wed Jan 10 11:46:12 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.2   Wed Jan 10 09:54:46 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:43:08   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:06   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
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
CD5 mb            INT     I     First coordinate number
CD5 mc            INT     I     Second coordinate number
CD5 id            INT     I     Flag to decide whether to interpolate
CD5                                on nodes or elements
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
CD6 cord, ns, nelm
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
CD7   mdif            INT      Difference in the two indexes
CD7   mdifm           INT      mdif - 1
CD7   dx              real*8   Spacing in x-direction
CD7   dy              real*8   Spacing in y-direction
CD7   dz              real*8   Spacing in z-direction
CD7   i               INT      Do loop index
CD7   io              INT      Do loop index
CD7   ir              INT      Do loop index
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN gendat
CPS
CPS Compute integer values
CPS
CPS IF we are interpolating on coordinates
CPS
CPS   Compute delta x, y, and z between points
CPS
CPS   FOR each added point
CPS     Compute coordinate value
CPS   ENDFOR
CPS
CPS ELSE we are interpolating on elements
CPS
CPS   FOR each node on this element
CPS      Compute index used later
CPS   ENDFOR
CPS
CPS   FOR each added point
CPS     FOR each node on this element
CPS       Compute element number
CPS     ENDFOR
CPS   ENDFOR
CPS
CPS ENDIF
CPS
CPS END gendat
CPS
C***********************************************************************

      use combi
      use comci
      use comdti
      use comai
      implicit none

      integer nb(20),mb,mc,id,io,i,ir,mdif,mdifm
      real*8 dx,dy,dz

      mdif   =  mb-mc
      mdifm  =  mdif-1
      if (id.ne.2)  then
c     
c**** interpolate on coordinates ****
c
         dx     =  ( cord(mb,1)-cord(mc,1) )/mdif
         dy     =  ( cord(mb,2)-cord(mc,2) )/mdif
         dz     =  ( cord(mb,3)-cord(mc,3) )/mdif
         do io = 1, mdifm
            cord(mc+io,1) =  cord(mc,1)+io*dx
            cord(mc+io,3) =  cord(mc,3)+io*dz
            cord(mc+io,2) =  cord(mc,2)+io*dy
         end do
      else
c
c**** interpolate on elements ****
c
         do i = 1, ns
            nb(i)  =  ( nelm((mb-1)*ns+i)-nelm((mc-1)*ns+i) )/mdif
         end do
         do io = 1, mdifm
            do ir = 1, ns
               nelm((mc+io-1)*ns+ir) =  nelm((mc-1)*ns+ir)+io*nb(ir)
            end do
         end do
      endif

      return
      end

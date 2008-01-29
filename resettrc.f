      subroutine resettrc(daytr)
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
CD1 To cut the time step of the solute transport calculation and reset
CD1 the values of all parameters before reinitiating the calculation.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-09-94      B. Robinson   N/A     Initial Implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/resettrc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:20   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:02 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Fri Jun 28 09:04:36 1996   zvd
CD2 Added timestep error termination write to ierr
CD2 
CD2    Rev 1.3   Wed May 08 13:35:38 1996   hend
CD2 Updated trac time step manipulations
CD2 
CD2    Rev 1.2   Thu Feb 01 15:52:12 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   04/25/95 10:09:48   llt
CD2 retrieved lost log history information
CD2
CD2    Rev 1.0   03/23/94 14:49:44   robinson
CD2 Initial Implementation
CD2
C**********************************************************************
CD3 
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier Type    Use   Description
CD3
CD3 daytr      real*8  I/O   Tracer time step
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
CD4 rc           real*8  fdd1    Source/sink and chemical reaction terms
CD4                              of the tracer mass balance
CD4 dench        real*8  fdd1    Solute mass storage at previous time
CD4                              for each node
CD4 anl          real*8  fdd1    Current concentration at each node
CD4                              for the current tracer
CD4 anlo         real*8  fdd1    Current concentration at each node
CD4                              for the current tracer
CD4 dencj        real*8  fdd1    Solute mass storage at current time
CD4                              for each node
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
CD5 ispecies     int         Do loop index over all solutes
CD5 inode        int         Do loop index over all nodes
CD5 iconc        int         Index of current concentration value
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 None
CD6
C**********************************************************************
CD7
CD7 SPECIAL COMMENTS
CD7
CD7  Requirements from SDN: 10086-RD-2.20-00
CD7    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD7    FEHM Application Version 2.20
CD7
C**********************************************************************
CD8
CD8 REQUIREMENTS TRACEABILITY
CD8
CD8 2.3.4 Solute-transport equations
CD8 2.5.1 Implement time-step mechanism
CD8
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN resettrc
CPS 
CPS Reset the simulation times
CPS Cut time step
CPS 
CPS FOR each species
CPS   FOR each node
CPS     Reset concentration to value at previous time
CPS     Set storage array value to 0
CPS   ENDFOR
CPS ENDFOR
CPS 
CPS END resettrc
CPS 
C**********************************************************************

      use comci
      use comdi
      use comdti
      use comai
      implicit none

      real*8 daytr
      integer ispecies,inode,iconc

      days = days - daytr
      daysi = days
      daytr = daytr/2.
      if (daytr.lt.daymin) then
         if (iout .ne. 0) then
            write(iout, 100) 
            write(iout, 200) 
            write(iout, 300) 
            write(iout, 100) 
         end if
         write(ierr, 100)
         write(ierr, 200) 
         write(ierr, 300) 
         write(ierr, 100)
         if (iptty.gt.0) then
            write(iptty, 100)
            write(iptty, 200) 
            write(iptty, 300)
            write(iptty, 100)
         endif
         stop
      endif
      do ispecies = 1, nspeci
         do inode = 1, neq
            iconc = (ispecies-1)*neq +inode
            an(iconc) = anlo(iconc)
            denci(iconc) = 0.
         end do
      end do
 100  format ('******************************************')
 200  format ('Tracer Time Step Smaller Than Minimum Step')
 300  format ('Stop in resettrc')

      return
      end

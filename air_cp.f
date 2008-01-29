      subroutine  air_cp(temp_celsius, heatcap, dheatcapt)
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
CD1 To compute the heat capacity of air (J/kg-K) and the
CD1 derivative with respect to temperature.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 8-14-94      B. Robinson    22      Initial Implementation
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/air_cp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Jan 29 13:10:42 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.1   Mon Jan 29 09:52:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.0   08/22/94 11:21:04   llt
CD2 Original version
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Name          Type        Description
CD3 
CD3 temp_celsius  real*8      Temperature in deg. C
CD3 heatcap       real*8      Heat capacity of air
CD3 dheatcapt     real*8      Derivative of heat capacity
CD3                              with temperature
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
CD4 None
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4 None
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
CD5 a_air        real*8      Constants in heat capacity
CD5                             calculation
CD5 a_aird       real*8      Constants in heat capacity
CD5                             derivative calculation
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop index
CD5 temp_celsius real*8      Temperature (deg. C)
CD5 tplus        real*8      Running product in heat capacity
CD5                             calculation
CD5 tplusd       real*8      Running product in heat capacity
CD5                             derivative calculation
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
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
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA Sychev, V. V. et al., Thermodynamic Properties of Air, Hemisphere
CDA Publishing Corp., 1988.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN air_cp
CPS 
CPS   Initialize heat capacity and derivative to 0
CPS   Initialize running products
CPS
CPS   FOR each term in series
CPS     Compute running product terms
CPS     Compute next term in series
CPS   ENDFOR
CPS
CPS END air_cp
CPS 
C**********************************************************************
      implicit none
      integer i
      real*8  a_air(4),a_aird(4),temp_celsius,heatcap,dheatcapt,
     &     tplus, tplusd
      data a_air / 1003.7, .025656, .00045457, -2.7107e-7 /
      data a_aird / 0., .025656, .00090914, -8.1321e-7 /
      
      heatcap = a_air(1)
      dheatcapt = 0.0
      tplus = 1.
      
      do i = 2,4
         tplusd = tplus
         tplus = tplus*temp_celsius
         heatcap = heatcap + a_air(i)*tplus
         dheatcapt = dheatcapt + a_aird(i)*tplusd
      enddo

      return
      end




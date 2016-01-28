      subroutine vaporl(tl,pcap,dpcaps,ivapl,delps,ddelt,ddels)
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
CD1 To calculate the vapor pressure lowering contribution to
CD1 saturation pressure.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 ?            G. Zyvoloski   N/A     Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/vaporl.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Fri Feb 02 14:10:08 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:10:16   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:08   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 tl           real*8  I       Temperature
CD3 pcap         real*8  I       Capillary pressure
CD3 dpcaps       real*8  I       Derivative of cap. pressure with
CD3                                  saturation
CD3 ivapl        real*8  I       Control parameter
CD3 delps        real*8  O       Vapor pressure lowering
CD3 ddelt        real*8  O       Derivative of vapor pressure with
CD3                                  temperature
CD3 ddels        real*8  O       Derivative of vapor pressure with
CD3                                  saturation
CD3                                  
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 File with number iout  O    File used to write warning and error
CD3                             messages (screen output if specified
CD3                             as 6)
CD3 File with number iptty O    File used to write warning and error
CD3                             messages (screen output if specified
CD3                             as 6)
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 None
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
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
CD5 dlgs         real*8      Density of water times gas constant
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5 
CD5 None
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
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN vaporl
CPS 
CPS IF vapor pressure lowering is included
CPS   Compute vapor pressure lowering and derivatives
CPS ENDIF
CPS 
CPS 
CPS END vaporl
CPS
C**********************************************************************

      implicit none

      real*8 delps,tl,pcap,dpcaps,dlgs,ddels,ddelt
      integer ivapl
      data dlgs/0.461/
c     delps-vapor pressure lowering
c     tl-temperature
c     pcap-capillary pressure
c     dpcaps-der wrt saturation of the capillary pressure
c     ivapl-control parameter
c     ddels-derivative of vapor pressure wrt saturation
c     ddelt-derivative of vapor pressure wrt temperature
c     dlgs- density of liquid water times the gas constant
      
      if(ivapl.eq.1) then
c     calculate vapor pressure lowering
         delps=exp(-pcap/(dlgs*(tl+273.0)))
         ddels=delps*(-dpcaps/(dlgs*(tl+273.0)))
         ddelt=delps*(pcap/(dlgs*(tl+273.0)**2))
      endif
      
      return
      end
      subroutine vaporl_salt(tl0,salt_con,pv_sc,dsct,dscc)
c
c correction for vapor pressure lowering from salt
c this function uses parameters fitted using PEST
c and data from Sparrow.
c
c tl-temperature(C)
c salt_con- salt concentration (moles)
c water vapor pressure correction (Mpa)
c
       implicit none
       real*8 ms,xf,xm
       real*8 pvwn,pvwn1,pvwn2
       real*8 pvwd,pvwd1,pvwd2
       real*8 pvw0,salt_corr
       real*8 dla0,dlta1,dlta2,dlpa1,dlpa2,bcoef1
       real*8 dlb0,dlpb1,dlpb2,dltb1,dltb2
       real*8 salt_con,tl,tl0,tl_norm,dpvwn2t,dpvwd2t
       real*8 pv_sc,dsct,dscc
       parameter (tl_norm = 360.0d00)

       dla0  =     2.119971900000000 
       dlpa1 =     3.000000000000000 
       dlpa2 =     0.738870490000000 
       dlta1 =     0.000000000000000 
       dlta2 =     1.000000000000000 
       dlb0  =     1.236252800000000 
       dlpb1 =     0.000000000000000 
       dlpb2 =     1.000000000000000 
       dltb1 =    -0.972846800000000 
       dltb2 =     2.377151700000000 
       bcoef1 =   -1.739785000000000 

          tl = tl0/tl_norm
          ms = salt_con*58.55 / 1000.
          xf = ms / ( ms + 1.0 ) 
          xm = xf 
          pvwn1=dla0
          pvwn2=dlta1*tl**dlta2
          pvwn=pvwn1+pvwn2
          pvwd1=dlb0+dlpb1*xm**dlpb2
          pvwd2=dltb1*tl**dltb2
          pvwd=pvwd1+pvwd2
          pvw0=pvwn/pvwd
          
          pv_sc = -dlpa1*xm**dlpa2*(pvw0+bcoef1)
          dpvwn2t = dlta1*dlta2*tl**(dlta2-1.)
          dpvwd2t = dltb1*dltb2*tl**(dltb2-1.)
          dsct =-((dpvwn2t/pvwd - (pvwn/pvwd2**2)*dpvwd2t)*(1./tl_norm))
c tracer is explicitly updated (derivative wrt c = 0)
          dscc = 0.0

       return
       end

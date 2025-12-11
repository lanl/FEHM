      subroutine vgrlp( sl, star, alpha, beta, hmin,
     2     hp, dhp, rl, drls, rv, drvs, iflg )
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
CD1 To compute the liquid relative permeability and derivative for the
CD1 van Genuchten model.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 3-3-94       Bruce Robinson        Initial implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/vgrlp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:23:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:12   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Fri Sep 26 15:14:10 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.5   Fri Feb 02 14:14:54 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   04/25/95 08:47:38   llt
CD2 retrieved lost log history
CD2 
CD2    Rev 1.3   04/25/95 08:15:18   llt
CD2 added log history information to prolog
CD2
CD2    Rev 1.2   01/28/95 13:56:32   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.1   08/23/94 08:23:00   llt
CD2 gaz changes
CD2
CD2    Rev 1.0   03/23/94 14:55:34   robinson 
CD2 Initial implementation
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 sl           real*8   I      Liquid saturation
CD3 star         real*8   I      Normalized liquid saturation
CD3 alpha        real*8   I      Van Ganuchten parameter
CD3 beta         real*8   I      Van Ganuchten parameter
CD3 hmin         real*8   I      minimum capillary pressure
CD3 hp           real*8   O      Capillary pressure
CD3 dhp          real*8   O      Derivative of capillary pressure with
CD3                                  respect to saturation
CD3 rl           real*8   O      Liquid relative permeability
CD3 drls         real*8   O      Derivative of liquid relative
CD3                                  permeability with respect to
CD3                                  saturation
CD3 rv           real*8   O      Vapor relative permeability
CD3 drvs         real*8   O      Derivative of vapor relative
CD3                                  permeability with respect to
CD3                                  liquid saturation
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 NONE
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4 
CD4 NONE
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 NONE
CD4 
CD4 Global Subprograms
CD4 
CD4 NONE
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 NONE
CD5 
CD5 Local Types
CD5
CD5 NONE
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 alamda       real*8      Exponent in correlation
CD5 denom        real*8      Denominator in normalized saturation
CD5                              expression
CD5 hp           real*8      Capillary pressure before unit conversion
CD5 dhp          real*8      Derivative of hp with respect to saturation
CD5 dahp         real*8      Derivative of ahp with respect to
CD5                              capillary pressure
CD5 term1        real*8      Term used in relative permeability
CD5                              calculation
CD5 dterm1       real*8      Derivative of term1
CD5 ahp          real*8      Term used in relative permeability
CD5                              calculation
CD5 al2          real*8      Exponent used in relative permeability
CD5                              calculation
CD5 ratio        real*8      Term used in relative permeability
CD5                              calculation
CD5 term2        real*8      Term used in relative permeability
CD5                              calculation
CD5 terma        real*8      Term used in relative permeability
CD5                              calculation
CD5 dratio       real*8      Derivative of ratio
CD5 dterma       real*8      Derivative of terma
CD5 dterm2       real*8      Derivative of term2
CD5 termhp1      real*8      Term used in relative permeability
CD5 termhp2      real*8      Term used in relative permeability
CD5 ahp1         real*8      Term used in relative permeability
CD5 termahp1     real*8      Term used in relative permeability
CD5 termahp2     real*8      Term used in relative permeability
CD5                              calculation
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 N/A
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
CD8 2.4.4 Relative-permeability and capillary-pressure functions
CD8
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN vgcalc
CPS 
CPS Compute parameters used in calculation
CPS 
CPS IF capillary pressure less than minimum value
CPS
CPS   Set capillary pressure to minimun value
CPS
CPS ENDIF
CPS
CPS IF saturation is within valid range (0 to 1)
CPS 
CPS   Compute liquid relative permeability and derivative
CPS 
CPS ELSEIF saturation is less than 0
CPS 
CPS   Set relative permeability to 0
CPS   Set derivative to 0
CPS 
CPS ELSE saturation is greater than 1
CPS 
CPS   Set relative permeability to 1
CPS   Set derivative to 0
CPS   
CPS ENDIF
CPS 
CPS Compute vapor relative permeability and derivative with respect...
CPS ... to liquid saturation
CPS 
CPS END vgcalc
CPS
C**********************************************************************

      implicit none

      real*8 sl  
      real*8 star
      real*8 alpha
      real*8 beta
      real*8 hp
      real*8 dhp
      real*8 rl
      real*8 drls
      real*8 rv
      real*8 drvs
      real*8 hmin


      real*8 alamda
      real*8 dahp
      real*8 term1
      real*8 ahp
      real*8 al2
      real*8 dterm1
      real*8 ratio
      real*8 term2
      real*8 terma
      real*8 dratio
      real*8 dterma
      real*8 dterm2
      real*8 termhp1
      real*8 termhp2
      real*8 ahp1
      real*8 termahp1
      real*8 termahp2
c
      real*8 term3,term4,term5,term6
      real*8 term7,term8
c
      integer iflg
c
c check for minumum capillary pressure
c
c gaz debug 020824
      if(iflg.le.1) then
c      if(iflg.eq.2) then
         if (hp .lt. hmin ) hp=hmin
         if(sl .gt.0.0.and.sl .lt.1.00) then
c     calculate the relative permeability
            alamda = 1.0-1.0/beta
            al2 = alamda/2.0
            termhp1 = alpha*hp
            termhp2 = termhp1**(beta-1.)
            ahp = termhp1*termhp2
            dahp = beta*termhp2*alpha
            ahp1 = 1. + ahp
            termahp1 = ahp1**al2
            termahp2 = ahp1 * termahp1
            term1 = 1. / termahp1
            dterm1 = -al2 / termahp2
            ratio =ahp/ahp1
            terma=(1.0-ratio**alamda)
            term2 = terma*terma
            dratio = 1.0/ahp1-ahp/ahp1**2
            dterma=-alamda/ratio**(1.0-alamda)*dratio
            dterm2 = 2.0*terma*dterma
            rl = term1*term2
            drls=(dterm1*term2+term1*dterm2)*dhp*dahp
         else if(sl .le.0.0) then
c     lower residual cutoff
            rl = 0.0
            drls= 0.0
c     upper residual cutoff
         else
            rl = 1.0
            drls= 0.0
         endif
         rv = 1.-rl
         drvs=-drls
      else if(iflg.eq.1) then
         if (hp .lt. hmin ) hp=hmin
         if(sl .gt.0.0.and.sl .lt.1.00) then
c     calculate the relative permeability
            alamda = 1.0-1.0/beta
            al2 = alamda/2.0
            termhp1 = alpha*hp
            termhp2 = termhp1**(beta-1.)
            ahp = termhp1*termhp2
            dahp = beta*termhp2*alpha
            ahp1 = 1. + ahp
            termahp1 = ahp1**al2
            termahp2 = ahp1 * termahp1
            term1 = 1. / termahp1
            dterm1 = -al2 / termahp2
            ratio =ahp/ahp1
            terma=(1.0-ratio**alamda)
            term2 = terma*terma
            dratio = 1.0/ahp1-ahp/ahp1**2
            dterma=-alamda/ratio**(1.0-alamda)*dratio
            dterm2 = 2.0*terma*dterma
            rl = term1*term2
            drls=(dterm1*term2+term1*dterm2)*dhp*dahp
         else if(sl .le.0.0) then
c     lower residual cutoff
            rl = 0.0
            drls= 0.0
            rv =1.0
            drvs = 0.0
c     upper residual cutoff                                                            X
         else
            rl = 1.0
            drls= 0.0
            rv = 1.0
            drvs= 0.0
         endif
c     part added by rosangela (sept 2002)
         if(sl .gt.0.0.and.sl .lt.1.00) then
c     calculate the relative permeability
          term3=term1*term1
          term4=(1-term3)**(1/2)
          term5=term3**(1./alamda)
          term6=(1-term5)**(2*alamda)
          term7=term3**(1./beta-1.)
          term8=(1/(1-term3))+4*term7*(1/(1-term5))       
          rv=term4*term6
          drvs=(-0.5)*rv*term8
         endif
      endif
      return
      end

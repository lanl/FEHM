      subroutine vgrlps(iflg, sl,  sl1, sl2, alamda,tol_l,
     2     tol_u, rl, drls, rv, drvs )
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
CD1 To compute the liquid and gas relative permeabilities and derivatives for the
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
CD2 $Log:   /pvcs.config/fehm90/src/vgrlps.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
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
CD7 Requirements from SDN: 10086-RD-2.20-00
CD7   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD7   FEHM Application Version 2.20
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
! calculates van genuchten rlperm as a function of saturation (sl)
! alpha is not used
! iflg = 1: rv=1-rl
! iflg=2 :  more complicated calculation of rv
! requires iflg, sl, star,beta, sl1,sl2,tol_l,tol_u
! returns rl, drls, rv, drvs
C**********************************************************************

      implicit none

      integer iflg
c gaz 021424 modified code to have variable cutoff tolerances
      real*8 sl  
      real*8 star
      real*8 star1
      real*8 starm1
      real*8 alpha
      real*8 beta
      real*8 rl
      real*8 drls
      real*8 rv
      real*8 drvs
      real*8 tol_l
      real*8 tol_u
      real*8 sl1
      real*8 sl2
      real*8 dstar
      real*8 dstarm1

      real*8 alamda
      real*8 alamdai
      real*8 alamda2
      real*8 term1
      real*8 term2
      real*8 term3
      real*8 term4
      real*8 dterm1
      real*8 dterm2
      real*8 dterm3
      real*8 dterm4
      real*8 rlp_l
      real*8 rlp_u
      real*8 rvp_l
      real*8 rvp_u
c
c rlperms as a function of saturation 
c
      star=(sl-sl1)/(sl2-sl1)

      if(iflg.eq.0) then
c calculate interpolated values
      else if(iflg.ne.0) then
      star=(sl-sl1)/(sl2-sl1)
       dstar = 1.0/(sl2-sl1)
       if(star.gt.tol_l.and.star.lt.1.0-tol_u) then
 
c     calculate the relative permeability
c          alamda = 1.0-1.0/beta
          alamdai = 1.0/alamda
          term1 = sqrt(star)
          term2 = star**alamdai
          term3 = (1.0-term2)**alamda
          term4 = (1.0-term3)**2
          dterm1 = 0.5/term1
          dterm2 = alamdai*star**(alamdai-1.0)
          dterm3 = -alamda/(1.0-term2)**(1.0-alamda)*dterm2
          dterm4 = -2.0*(1.0-term3)*dterm3
          rl=term1*term4
          drls = (dterm1*term4+term1*dterm4)*dstar
       else if(star.gt.0.0.and.star.le.tol_l) then
c     lower residual cutoff
          star1= tol_l
c          alamda = 1.0-1.0/beta
          alamdai = 1.0/alamda
          term1 = sqrt(star1)
          term2 = star1**alamdai
          term3 = (1.0-term2)**alamda
          term4 = (1.0-term3)**2
          rlp_l=term1*term4
          rl = rlp_l/tol_l*star
          drls= rlp_l/tol_l*dstar
c     upper residual cutoff
       else if(star.ge.1.0-tol_u.and.star.lt.1.0) then
          star1= 1.0-tol_u
c          alamda = 1.0-1.0/beta
          alamdai = 1.0/alamda
          term1 = sqrt(star1)
          term2 = star1**alamdai
          term3 = (1.0-term2)**alamda
          term4 = (1.0-term3)**2
          rlp_u=term1*term4
          rl = (rlp_u-1.0)/(tol_u)*(1.0-star)+1.0
          drls= -(rlp_u-1.0)/(tol_u)*dstar
       else if(star.lt.0.0) then
c     lower residual cutoff
          rl = 0.0
          drls= 0.0
c     upper residual cutoff
c     gaz 091210
       else if(star.ge.1.0) then
          rl = 1.0
          drls= 0.0
       endif
c gaz 070720 fixed # in column 1  with c      
c now block beginning on line 275 is ended
      endif
c now block beginning on line 270 is ended

      if(iflg.eq.1) then
c rv is 1. -rl
       rv = 1.-rl
       drvs=-drls
      else if(iflg.eq.2) then
c rv has own vg function
       dstar = 1.0/(sl2-sl1)
       starm1 = 1.0-star
       dstarm1 = -dstar
       if(starm1.gt.tol_l.and.starm1.lt.1.0-tol_u) then
c     calculate the relative permeability
c          alamda = 1.0-1.0/beta
          alamdai = 1.0/alamda
          alamda2 = 2.0*alamda
          term1 = sqrt(starm1)
          term2 = star**alamdai
          term3 = (1.0-term2)**alamda2
          dterm1 = -0.5/term1
          dterm2 = alamdai*star**(alamdai-1.0)
          if(alamda2-1.0.lt.0.0) then
           dterm3 = -alamda2/(1.0-term2)**(1.0-alamda2)*dterm2
          else
           dterm3 = -alamda2*(1.0-term2)**(alamda2-1.0)*dterm2
          endif
          rv=term1*term3
          drvs = (dterm1*term3+term1*dterm3)*dstar
       else if(starm1.gt.0.0.and.starm1.le.tol_l) then
c     lower residual cutoff
          star1= tol_l
c          alamda = 1.0-1.0/beta
          alamdai = 1.0/alamda
          alamda2 = 2.0*alamda
          term1 = sqrt(star1)
          term2 = (1.0-star1)**alamdai
          term3 = (1.0-term2)**alamda2
          rvp_l=term1*term3
          rv = rvp_l/tol_l*starm1
          drvs= rvp_l/tol_l*dstarm1
c     upper residual cutoff
       else if(starm1.ge.1.0-tol_u.and.starm1.lt.1.0) then
          star1= 1.0-tol_u
c          alamda = 1.0-1.0/beta
          alamdai = 1.0/alamda
          alamda2 = 2.0*alamda
          term1 = sqrt(star1)
          term2 = (1.0-star1)**alamdai
          term3 = (1.0-term2)**alamda2
          rvp_u=term1*term3
          rv = rvp_u/tol_u*(1.0-starm1)+1.0       
          drvs= -rvp_u/tol_u*dstarm1
       else if(starm1.le.0.0) then
c     lower residual cutoff
          rv = 0.0
          drvs= 0.0
c     upper residual cutoff
       else if(starm1.ge.1.0) then
          rv = 1.0
          drvs= 0.0
       endif
      endif

      end subroutine vgrlps

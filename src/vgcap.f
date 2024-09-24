      subroutine vgcap( sl, slr, smr, alpha, beta, ac1, ac2, ac3, ac4,
     2     smcut, sucut, bc3, bc4, hp, dhp ,ireg)
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
CD1 To compute the capillary pressure and derivatives for the 
CD1 van Genuchten model.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 $Log:   /pvcs.config/fehm90/src/vgcap.f_a  $ 
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:22:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:08   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Fri Sep 26 15:14:08 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.8   Fri Feb 16 13:25:48 1996   zvd
CD2 Spelling correction.
CD2 
CD2    Rev 1.7   Fri Feb 02 14:13:20 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   06/02/95 10:39:52   llt
CD2 gaz changes
CD2 
CD2    Rev 1.5   05/08/95 09:55:24   gaz
CD2 fixed criteria for upper saturation fit (gaz)
CD2 
CD2    Rev 1.4   04/25/95 10:19:44   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.3   03/23/95 19:26:14   gaz
CD2 gaz modified calling sequence for fits to wet region
CD2 
CD2    Rev 1.2   08/23/94 08:22:56   llt
CD2 gaz changes
CD2
CD2    Rev 1.1   06/03/94 15:41:40   gaz
CD2 Made capillary pressure continuous from residual saturation
CD2 to zero saturation
CD2
CD2    Rev 1.0   03/23/94 14:55:32   robinson
CD2 Initial implementation
CD2
CD2 3-3-94       Bruce Robinson        Initial implementation
CD2                                     
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 sl           real*8   I/O    Liquid saturation on input,
CD3                                  normalized liquid saturation on
CD3                                  output
CD3 slr          real*8   I      Residual liquid saturation
CD3 smr          real*8   I      Maximum liquid saturation
CD3 alpha        real*8   I      Van Ganuchten parameter
CD3 beta         real*8   I      Van Ganuchten parameter
CD3 ac1          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac2          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac3          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac4          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 smcut        real*8   I      Lower cutoff saturation below which a
CD3                                  spline fit is used for capillary
CD3                                  pressure versus saturation
CD3 sucut        real*8   I      Upper cutoff saturation above which a
CD3                                  linear fit is used for capillary
CD3                                  pressure versus saturation
CD3 bc3          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 bc4          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 hp           real*8   O      Capillary pressure
CD3 dhp          real*8   O      Derivative of capillary pressure with
CD3                                  respect to saturation
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
CD5 ds           real*8      Reciprocal of demoninator in normalized
CD5                              saturation expression
CD5 denom        real*8      Denominator in normalized saturation
CD5                              expression
CD5 star         real*8      Normalized saturation
CD5 alpi         real*8      Exponent in correlation
CD5 hp           real*8      Capillary pressure before unit conversion
CD5 dhp          real*8      Derivative of hp with respect to saturation
CD5 termstar1    real*8      Term used in capillary pressure calculation
CD5 termstar2    real*8      Term used in capillary pressure calculation
CD5 termb1       real*8      Term used in capillary pressure calculation
CD5 termb2       real*8      Term used in capillary pressure calculation
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
CPS BEGIN vgcap
CPS 
CPS Compute parameters used in calculation
CPS 
CPS IF saturation is within valid range (0 to 1)
CPS 
CPS   IF saturation is within range where normal correlation is valid
CPS     Compute capillary pressure and derivative with respect to...
CPS     ... saturation
CPS   ELSEIF saturation calculation requires low-end cubic spline...
CPS   ... calculation
CPS     Compute capillary pressure and derivative with respect to...
CPS     ... saturation
CPS   ELSE saturation calculation requires high-end cubic spline...
CPS   ... calculation
CPS     Compute capillary pressure and derivative with respect to...
CPS     ... saturation
CPS   ENDIF
CPS   
CPS   Perform units conversion on capillary pressure and derivative
CPS 
CPS ELSEIF saturation is less than 0
CPS 
CPS   Set capillary pressure to maximum value
CPS   Set derivatives to 0
CPS 
CPS ELSE saturation is greater than 1
CPS 
CPS   Set capillary pressure to minimum value
CPS   Set derivatives to 0
CPS   
CPS ENDIF
CPS 
CPS END vgcap
CPS
C**********************************************************************

      implicit none

      integer ireg,ncut
      real*8 sl
      real*8 slr
      real*8 smr
      real*8 alpha
      real*8 beta
      real*8 ac1
      real*8 ac2
      real*8 ac3
      real*8 ac4
      real*8 bc3
      real*8 bc4
      real*8 smcut
      real*8 sucut
      real*8 hp
      real*8 dhp


      real*8 alamda
      real*8 ds
      real*8 denom
      real*8 star
      real*8 alpi
      real*8 termstar1
      real*8 termstar2
      real*8 termb1
      real*8 termb2
      real*8 hp_cut
      parameter(ncut = 1)
      alamda = 1.0-1.0/beta
      denom = smr-slr
      star=(sl-slr)/denom
      ds = 1.0/denom
      ireg = 0
      if(sl.gt.0.0.and.sl.lt.1.0) then
c     check if within spline cutoffs
         if(star.gt.smcut.and.star.lt.sucut) then
c         if(star.gt.smcut.and.star.lt.sucut.and.ireg.ne.5
c     &      .and.ireg.ne.1) then
            alpi = 1.0/alamda
            termstar1 = star**alpi
            termstar2 = star*termstar1
            termb1 = 1./termstar1 - 1.
            termb2 = termb1**(-alamda)
            hp = 1./alpha * termb1*termb2
            dhp = (1.-alamda)/alpha * termb2 * (-alpi/termstar2) * ds
         else  if(star.le.smcut) then
c        else  if(star.le.smcut.or.ireg.eq.1) then
c     use linear interpolation for lower cutoff
            hp = ac1*sl**3 + ac2*sl**2 +ac3*sl + ac4
            dhp=(3.0*ac1*sl**2+2.0*ac2*sl+ac3)
         else  if(star.ge.sucut) then
c        else  if(star.ge.sucut.or.ireg.eq.5) then
c     use linear interpolation for upper cutoff
c            star,sucut = 0.99  hp_cut = 19.9976320885161
c            star,sucut = 0.99925  hp_cut = 4.75925174261405
c            hp = 19.9976320885161
c          hp_cut = 4.75925174261405
c          bc3 = hp_cut/(sucut-1.d0)
c          hp = bc3*(star-1.d0)
c          dhp = bc3*ds
c     use linear interpolation for upper cutoff
            hp = bc3*sl + bc4
            dhp= bc3
            hp = bc3/sl**ncut + bc4
            dhp = -bc3*ncut/sl**(ncut+1)
         endif

      else if(sl.le.0.0) then
c     lower saturation cutoff
         hp = ac4
         dhp= ac3
c     upper saturation cutoff
      else
         hp = bc3*sl + bc4
         dhp= bc3
         hp = 0.0 
         dhp= 0.0
c gaz  02-18-08
         hp = bc3/sl**ncut + bc4
         dhp = -bc3*ncut/sl**(ncut+1)
c ugta version 14-Jun-10, used following
c        dhp = 0.0
         ireg = 100
      endif
      return
      end

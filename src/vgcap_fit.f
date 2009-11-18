      subroutine vgcap_fit(iflg,slr,slm,slcut,smcut,fac,alpha
     &                     ,alamda,c1,c2,c3,c4,hmin)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  To provide linear or cubic fit to capillary pressure data.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD 
CD2 Date         Programmer     Number  Comments 
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/vgcap_fit.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:22:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:10   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Fri Nov 21 16:51:40 1997   gaz
CD2 removed some commented out code
CD2 
CD2    Rev 1.5   Fri Sep 26 15:14:08 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.4   Fri Feb 16 13:28:26 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.3   Thu Jan 11 12:29:54 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   04/25/95 08:47:36   llt
CD2 retrieved lost log history
CD2 
CD2    Rev 1.1   04/25/95 08:15:16   llt
CD2 added log history information to prolog
CD2 
CD2    Rev 1.0   03/23/95 17:53:16 pvcs
CD2 initial release
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.4.4	Relative-permeability and capillary-pressure functions
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
CP5 iflg    int    flag to designate fitting procedure
CP5 slr     real*8 residual liquid saturation
CP5 slm     real*8 maximum liquid saturation
CP5 slcut   real*8 cutoff liquid saturation
CP5 smcut   real*8 normalized cutoff liquid saturation
CP5 fac     real*8 multiple of capillary pressure 
CP5 alpha   real*8 parameter in van Genecthen capillary model
CP5 alamda  real*8 parameter in van Genecthen capillary model
CP5 c1      real*8 cubic coefficient in polynomial fit
CP5 c2      real*8 quadratic coefficient in polynomial fit
CP5 c3      real*8 linear coefficient in polynomial fit
CP5 c4      real*8 constant coefficient in polynomial fit
C***********************************************************************

      implicit none

      integer iflg,ncut 
      real*8 slr,slm,slcut,alpha,alamda,c1,c2,c3,c4 
      real*8 hmax,hmin,hcut,smcut,alpi,ds,dhp,slope
      real*8 fac,fac_min,fac_use,hcut2,smcut2,ancut
      parameter(fac_min=2.0,ncut = 1)
      if(iflg.eq.1) then
c zero slope, zero curvature at sl=0
c h = a*x**3 + d
c this is most common fehm way
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
            ds = 1.0/(slm-slr)
            dhp = 1.0/alpha*(1.0-alamda)/
     *              (1.0/smcut**alpi-1.0)**alamda
     *              *(-alpi/smcut**(alpi+1.0))*ds
            fac=1.0/3.0
       hmax=hcut-fac*dhp*slcut
       c1=dhp/(3.0*slcut**2)
       c2=0.0
       c3=0.0
       c4=hmax
      else if(iflg.eq.2) then
c linear fit from slcut to sl=0.0
c with hmax=fac*hcut
            fac_use=max(fac,fac_min)
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
       hmax=fac_use*hcut
       slope=-(hmax-hcut)/slcut
       c3=slope
       c4=hmax
      else if(iflg.eq.3) then
c linear fit from slcut to sl=0.0
c with hmax: linear extension of slope at slcut
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
            ds = 1.0/(slm-slr)          
            dhp = 1.0/alpha*(1.0-alamda)/
     *              (1.0/smcut**alpi-1.0)**alamda
     *              *(-alpi/smcut**(alpi+1.0))*ds
       hmax=hcut-dhp*slcut
       slope=-(hmax-hcut)/slcut
       c3=slope
       c4=hmax
      else if(iflg.eq.4) then
c linear fit from smcut to star= 1.00
c with hmax: linear extension of slope at slcut
            alpi = 1.0/alamda
            hcut = 1.0/alpha*(1.0/smcut**alpi-1.0)**(1.0-alamda)
       slope=(hcut-hmin)/(1.0-slcut)
       c3=-slope
       c4=slope+hmin
       ancut = (1./slcut**ncut-1.0)
       c3 = hcut/ancut
       c4 = -c3 + hmin
      continue
      endif
      return
      end

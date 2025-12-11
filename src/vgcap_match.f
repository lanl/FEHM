
      subroutine vgcap_match(iflg,smcutm,alpham,alamdam,facm
     &                     ,smcutf,alphaf,alamdaf,facf)
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
CD1 To match the maximum capillary pressures in the fracture and matrix.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 $Log:   /pvcs.config/fehm90/src/vgcap_match.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:23:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
CD2
CD2 Split subroutine out of vgcap_fit 08-Feb-02
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
CP5 smcut   real*8 normalized cutoff liquid saturation
CP5 alpha   real*8 parameter in van Genecthen capillary model
CP5 alamda  real*8 parameter in van Genecthen capillary model
C***********************************************************************

      implicit none

      integer iflg
      real*8 smcutm,alpham,alamdam,facm
      real*8 smcutf,alphaf,alamdaf,facf
      real*8 hcut,alpi,hcutm,hcutf,facfmax
       if(iflg.eq.1) then
c change cutoff saturation so pressures match
            alpi = 1.0d00/alamdam
            hcut = 1.0d00/alpham*(1.0d00/smcutm**alpi-1.0d00)
     &                           **(1.0d00-alamdam)
            smcutf = 1.0d00/((alphaf*hcut)**(1.0d00/(1.0d00-alamdaf)))
     &               **alamdaf 
       else if(iflg.eq.2) then
c change fac so pressures match at s=0
            facfmax = facf
            alpi = 1.0d00/alamdam
            hcutm = 1.0d00/alpham*(1.0d00/smcutm**alpi-1.0d00)
     &                           **(1.0d00-alamdam)
            alpi = 1.0d00/alamdaf
            hcutf = 1.0d00/alphaf*(1.0d00/smcutf**alpi-1.0d00)
     &                           **(1.0d00-alamdaf)
            facf = min(hcutm*facm/hcutf,facfmax)
       endif
      return
      end

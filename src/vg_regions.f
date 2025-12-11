      subroutine vg_regions(iflg,ireg,mi,su_cut)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To determine subregions for vgcap curves.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2
!D2 $Log:   /pvcs.config/fehm90/src/vg_regions.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:22:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2 Split subroutine out of rlperm 08-Feb-02
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.4 Relative-permeability and capillary-pressure functions
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      use comdti
      use comai
      use comdi
      implicit none

      integer iflg,mi,ireg,it,icap0
      real*8 slcut,slcut_l,slcut_u,su_cut
      real*8 srcut,srcut_l,srcut_u
      real*8 fac01,fac02,fac03,fac04
      real*8 stepl_regions
c
      parameter (fac01=0.999d00,fac02=0.975d00)
      parameter (fac03=1.000,fac04=1.000)
      parameter (stepl_regions=1.00d00)
c
                  it=irlp(mi)
                  icap0=icap(mi)
               if(idualp.ne.0.or.idpdp.ne.0) then
c upper fitting region
                  if(iflg.eq.1) then
                   slcut=su_cut*(rp2f(it)-rp1f(it))+rp1f(it)
                   slcut_u=fac01*slcut
                   slcut_l=fac02*slcut
                  if(s(mi).ge.slcut_u.and.icap(mi).ne.5) then
                     icap(mi)=5    
                  elseif(s(mi).le.slcut_l.and.icap(mi).eq.5) then
                     icap(mi)=0    
                  endif
                  else if(iflg.eq.2) then
                   slcut=su_cut*(rp12f(it)-rp11f(it))+rp11f(it)
                   slcut_u=fac01*slcut
                   slcut_l=fac02*slcut
                  if(s(mi).ge.slcut_u.and.icap(mi).ne.5) then
                     icap(mi)=5    
                  elseif(s(mi).le.slcut_l.and.icap(mi).eq.5) then
                     icap(mi)=0    
                  endif
                  endif
c lower fitting region
                  if(iflg.eq.1) then
                   srcut=rp6f(it)*(rp2f(it)-rp1f(it))+rp1f(it)
                   srcut_l=fac03*srcut
                   srcut_u=fac04*srcut
                  if(s(mi).le.srcut_l.and.icap(mi).ne.1) then
                     icap(mi)=1    
                  elseif(s(mi).ge.srcut_u.and.icap(mi).eq.1) then
                     icap(mi)=0    
                  endif
                  else if(iflg.eq.2) then
                   srcut=rp16f(it)*(rp12f(it)-rp11f(it))+rp11f(it)
                   srcut_l=fac03*srcut
                   srcut_u=fac04*srcut
                  if(s(mi).le.srcut_l.and.icap(mi).ne.1) then
                     icap(mi)=1    
                  elseif(s(mi).ge.srcut_u.and.icap(mi).eq.1) then
                     icap(mi)=0    
                  endif
                  endif
               else
c equivalent continuum
c upper fitting region
                  if(iflg.eq.1) then
                   slcut=su_cut*(rp2f(it)-rp1f(it))+rp1f(it)
                   slcut_u=fac01*slcut
                   slcut_l=fac02*slcut
                  if(s(mi).ge.slcut_u.and.icap(mi).ne.5) then
                     if(icap(mi).eq.0) icap(mi)=5    
                     if(icap(mi).eq.6) icap(mi)=7    
                  elseif(s(mi).le.slcut_l.and.icap(mi).ge.5) then
                     if(icap(mi).eq.7) icap(mi)=6    
                     if(icap(mi).eq.5) icap(mi)=0    
                  endif
                  else if(iflg.eq.2) then
                   slcut=su_cut*(rp12f(it)-rp11f(it))+rp11f(it)
                   slcut_u=fac01*slcut
                   slcut_l=fac02*slcut
                  if(s(mi).ge.slcut_u.and.icap(mi).ne.6) then
                     if(icap(mi).eq.0) icap(mi)=6    
                     if(icap(mi).eq.5) icap(mi)=7    
                  elseif(s(mi).le.slcut_l.and.icap(mi).ge.5) then
                     if(icap(mi).eq.7) icap(mi)=5    
                     if(icap(mi).eq.6) icap(mi)=0    
                  endif
                  endif
c lower fitting region
                  if(iflg.eq.1) then
                   srcut=rp6f(it)*(rp2f(it)-rp1f(it))+rp1f(it)
                   srcut_l=1.00*srcut
                   srcut_u=1.00*srcut
                  if(s(mi).le.srcut_l.and.icap(mi).ne.1) then
                     if(icap(mi).eq.2) icap(mi)=3    
                     if(icap(mi).eq.0) icap(mi)=1    
                  elseif(s(mi).ge.srcut_u.and.icap(mi).eq.1
     &                   .or.icap(mi).eq.3) then
                     if(icap(mi).eq.3) icap(mi)=2    
                     if(icap(mi).eq.1) icap(mi)=0    
                  endif
                  else if(iflg.eq.2) then
                   srcut=rp16f(it)*(rp12f(it)-rp11f(it))+rp11f(it)
                   srcut_l=1.00*srcut
                   srcut_u=1.00*srcut
                  if(s(mi).le.srcut_l.and.icap(mi).ne.2) then
                     if(icap(mi).eq.1) icap(mi)=3    
                     if(icap(mi).eq.0) icap(mi)=2    
                  elseif(s(mi).ge.srcut_u.and.icap(mi).eq.2
     &                   .or.icap(mi).eq.3) then
                     if(icap(mi).eq.3) icap(mi)=1    
                     if(icap(mi).eq.2) icap(mi)=0    
                  endif
                  endif
               endif
                  if(icap(mi).ge.1.and.icap(mi).le.3) then
                     ireg=1
                  else if(icap(mi).ge.5) then
                     ireg=5
                  else
                     ireg=0        
                  endif
                  if(icap(mi).ne.icap0) then
                    strd = min(strd,stepl_regions)
                  endif
      return
      end

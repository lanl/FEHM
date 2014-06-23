      subroutine rlperm_wtsi(iflg)
!***********************************************************************
!  Copyright, 2006,  The  Regents  of the  University of California.
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
!D1 To compute the relative permeabilities for water table (wtsi) macro
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 2-14-06, Programmer: G. Zyvoloski
!D2                                     
!D2 $Log:   /pvcs.config/fehm90/src/rlperm_wtsi.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 ?
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
      use combi
      use comci
      use comdi
      use comgi
      use comii 
      use comfi
      use comwt
      implicit none
      integer iflg, inode, mm, ihead_ck, j
      integer iadtest, iad_ck_wtsi
      real*8  hmin, hmax, hfac_l, hfac_h, strd_wtsi
      real*8  rlzf_dum, dis_ck, head_ck1, head_ck2
      real*8  fac_wtsi
      parameter (iad_ck_wtsi= 1)
      parameter (strd_wtsi= 1.d00)
      parameter (hfac_h=1.0d0, hfac_l=1.00001d0, fac_wtsi= 10000.)
      save iadtest
      
      if(iflg.eq.0) then
         
      else if(iflg.eq.1) then
c         rho1grav = crl(1,1)*(9.81d-6)
c     
c     calculate relative perms and other free surface quantities
c     
c     first check and change phase state
c     
         if(iad.lt.iad_up_wtsi) then
            if(iad.eq.0) dry_zone = 0
            if(iad.eq.0) strd = 1.0d0
            do inode=1,n0
               if(izone_free_nodes(inode).eq.3) then
                  hmin=head12(inode,1)
                  hmax=head12(inode,2) 
                  if(phi(inode).lt.hmax*hfac_l
     &                 .and.phi(inode).gt.hmin*hfac_h) then
                     izone_free_nodes(inode)= 2
                     strd = strd_wtsi
                  else if(phi(inode).gt.hmax*hfac_h) then 
                     izone_free_nodes(inode)= 1
                     strd = strd_wtsi
                  endif
               else if(izone_free_nodes(inode).eq.2) then
                  hmin=head12(inode,1)
                  hmax=head12(inode,2) 
                  if(phi(inode).lt.hmin*hfac_l) then
                     izone_free_nodes(inode)= 3
                     strd = strd_wtsi
                  else if(phi(inode).gt.hmax*hfac_h) then 
                     izone_free_nodes(inode)= 1
                     strd = strd_wtsi
                  endif
               else 
                  hmin=head12(inode,1)
                  hmax=head12(inode,2)
                  if(phi(inode).lt.hmax*hfac_l) then
                     izone_free_nodes(inode)= 2
                     strd = strd_wtsi
                  endif
               endif   
            enddo
         endif


         ifree1 = 0
         do inode=1,n0 
            
            if(izone_free_nodes(inode).eq.1) then

               rlxyf(inode) = 1.d0
               drlxyf(inode) = 0.d0               
               
            else if(izone_free_nodes(inode).eq.2) then
               hmin=head12(inode,1)
               hmax=head12(inode,2)
               rlxyf(inode) = (phi(inode)-hmin)/(hmax-hmin) + rlptol
               drlxyf(inode) = 1.d0/(hmax-hmin)
               ifree1 = ifree1 +1
            else if(izone_free_nodes(inode).eq.3) then
               hmin=head12(inode,1)
               hmax=head12(inode,2)
               rlxyf(inode) = 0.d0 + rlptol
               drlxyf(inode) = 1.d0/(hmax-hmin) 
               ifree1 = ifree1 +1
            else
               rlxyf(inode) = 1.d0 + rlptol
               drlxyf(inode) = 0.0d0

            endif
c     
c     new code gaz 12-15-05 for  horizontal-vertical anisotropy
c     

            rlzf_dum = zfac_ani*rlxyf(inode)             
            if(rlzf_dum.le.rlptol) then
               rlzf(inode) = 0.d0
               drlzf(inode)  = zfac_ani*drlxyf(inode)
            else if(rlzf_dum.ge.1.d00) then
               rlzf(inode) = 1.d00
               drlzf(inode)  = 0.0d0
            else
               rlzf(inode) = rlzf_dum 
               drlzf(inode)= zfac_ani*drlxyf(inode)
            endif             


         enddo                 
         

         if(iad.ne.0) then
            head_ck = 0.0
            ihead_ck = 0

            do inode = 1,neq
               if(izone_free_nodes(inode).eq.2.and.ka(inode).ge.0) then
                  hmin=head12(inode,1)
                  hmax=head12(inode,2)
                  dis_ck=(hmax-hmin)/rho1grav
                  head_ck1 = (s(inode)-rlxyf(inode))*dis_ck
                  head_ck2 = abs(head_ck1)
                  if(head_ck2.gt.abs(head_ck)) then
                     ihead_ck = inode
                     head_ck = head_ck1
                  endif 
               else if(izone_free_nodes(inode).eq.1.and.ka(inode).ge.0)
     &                 then
                  head_ck1 = -bp(inode)/rho1grav
                  head_ck2 = abs(head_ck1)
                  if(head_ck2.gt.abs(head_ck)) then
                     ihead_ck = inode
                     head_ck = head_ck1			     
                  endif
               endif
            enddo
!            if (iout .ne. 0) then               
!               write(iout,3000) l, iad, head_ck, ihead_ck,
!     &              (cord(ihead_ck,j),j=1,3),izone_free_nodes(ihead_ck)
!            end if
c            if(iptty.ne.0) then
c               write(iptty,3000) l, iad, head_ck, ihead_ck,
c     &              (cord(ihead_ck,j),j=1,3),izone_free_nodes(ihead_ck)
c            endif
         endif
         
         if(iad.eq.1) then
	    iadtest = 0
	    head_ck_first = abs(head_ck)
         else if(iad.gt.1) then
	    
            if(abs(head_ck).gt.fac_wtsi*head_ck_first) then
               iadtest = iadtest +1
               if(iadtest.eq.iad_ck_wtsi) then
                  iad = abs(maxit)
                  if (iout .ne. 0) write(iout,*)
     &                 'restart time step due to head divergence'
                  if(iptty.ne.0) write(iptty,*)
     &                 'restart time step due to head divergence'
	       endif
            endif
	    
         endif
         s = rlxyf 
         
      else
         
      endif
 3000 format('timestep, iteration, max head diff, node, x,y,z,',
     &     ' grid state',/,1x,i5,1x,i3,1x,1p,g14.5,1x,i9,3(1x,g14.5),
     &     1x,i3) 

      return
      end

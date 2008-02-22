      subroutine uz_wtctr(iz)
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
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
CD1
CD1 PURPOSE
CD1
CD1 Calculate uz correction near water table to make pressures
CD1 continuous
CD1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: Date 2005, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/uz_wtctr.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
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

      use comai
      use combi
      use comdi
      use comdti
      use comii

      implicit none

      integer i,iz,ndummy,mid,it,mi,itp,itperm,itpperm,im
      integer inode, i1, i2, ii,kb
      real*8 cp1,cp3,hmax,hmin,hmid


      if(iwt_uz.ne.0) then
c     uz pressure correction
         if (iz.eq.99) then
	    allocate (head12(n0,2))
	    allocate (dhead12s(n0))
         else if (iz.eq.100) then
c     allocate arrays calculate distances
            
            do inode=1,n0
               i1=nelm(inode)+1
               i2=nelm(inode+1)
               hmid=cord(inode,igrav)
               hmin=0.
               hmax=0.
               do ii =i1,i2
                  kb=nelm(ii)
                  hmax=max(cord(kb,igrav)-hmid,hmax)
                  hmin=min(cord(kb,igrav)-hmid,hmin)
               enddo
               if(ivf.eq.-1) then
     	            head12(inode,1) = max(hmax,abs(hmin))
               else
                  head12(inode,1) = abs(hmax-hmin)/2.   
               endif

            enddo
         else if (iz.eq.101) then
c     addjust boundary files pressues
            
         else if (iz.eq.102) then
c     add saturation component to total pressure
            do kb = 1,neq
               head12(kb,2) = 0.0
               dhead12s(kb) = 0.0
               if(s(kb).le.0.0) then
                  head12(kb,2) = 0.
                  dhead12s(kb) = 0.
               else if(s(kb).ge.1.0) then
                  head12(kb,2) =  head12(kb,1)
                  dhead12s(kb) = head12(kb,1)
               else if(s(kb).gt.0.0.and.s(kb).lt.1.0) then
                  head12(kb,2) = head12(kb,1)*s(kb)
                  dhead12s(kb) = head12(kb,1)
               endif
               head12(kb,2) = head12(kb,1)*s(kb)
               dhead12s(kb) = head12(kb,1)
c     head12(kb,2) = 0.0
c     dhead12s(kb) = 0.0
	    enddo
            

            
         end if
      endif

      return
      end

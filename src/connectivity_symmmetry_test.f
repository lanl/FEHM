      subroutine connectivity_symmetry_test(iflg)   
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
CD1  check the connectivity of the connectivity array (nelm)
CD1  check istrw and sx arrays
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2  $Log:   /pvcs.config/fehm90/src/mdnodes.f_a  $
CD2  initial implementation george zyvoloski 11-29-09
CD2 *********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3  
CD3  2.2 Finite-Element Coefficient Generation
CD3  2.6 Provide Input/Output Data Files
CD3  3.0 INPUT AND OUTPUT REQUIREMENTS
CD3  
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
CD4  
CD4
C***********************************************************************

      use davidi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer iflg,neqp1
      integer i,i1,i2,i3,i4,jj,kk,kb,kc
      integer, allocatable :: ncon_new(:)
      integer, allocatable :: istrw_new(:)
      integer, allocatable :: idum(:)
      
      real*8, allocatable :: area_con(:,:)
      real*8, allocatable :: sx_new(:,:)
      real*8, allocatable :: dum(:)
c gaz 112223 changed ctest to  ctest_sym (avoid confusion with ctest (xnl)     
      logical ctest_sym

c=======================================================================

      if(iflg.eq.0) then
       ctest_sym = .true.
       neqp1 = neq + 1 
       do  i = 1,neq  
        i1 = nelmdg(i)+1
        i2 = nelm(i+1)
        do jj = i1,i2
         kb = nelm(jj)
         i3 = nelm(kb)+1
         i4 = nelmdg(kb)-1
         do kk = i3,i4
          kc = nelm(kk)
          if(kc.eq.i) then
           go to 100
          endif
         enddo
         ctest_sym = .false.
         write (ierr,*)  
     &   'connectivity failed, node (i,kb) = (',i,',',kb,')'
100     continue         
        enddo
       enddo
       if(ctest_sym) then
        if(iptty.ne.0) write(iptty,*)
        if(iptty.ne.0) 
     &   write(iptty,*) ' >>>> passed connectivity symmetry test <<<<'
        if(iptty.ne.0) write(iptty,*)
        if(iatty.ne.0) write(iatty,*)
        if(iatty.ne.0) 
     &   write(iatty,*) ' >>>> passed connectivity symmetry test <<<<'  
        if(iatty.ne.0) write(iatty,*)
       else
        if(iptty.ne.0) write(iptty,*)
        if(iptty.ne.0) 
     &   write(iptty,*) ' >>>> failed connectivity symmetry test <<<<'
        if(iptty.ne.0) write(iptty,*) ' >>>> stopping <<<<'
        if(iptty.ne.0) write(iptty,*) 
        if(iatty.ne.0) write(iatty,*) 
        if(iatty.ne.0) 
     &   write(iatty,*) ' >>>> failed connectivity symmetry test <<<<' 
        if(iatty.ne.0) 
     &   write(iatty,*)  ' >>>> stopping <<<<' 
        if(iatty.ne.0) write(iatty,*)
        stop 
       endif
      else if(iflg.eq.2) then

      endif
      return
      end  
      

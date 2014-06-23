      subroutine wtrise
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************
!D1
!D1  PURPOSE
!D1
!D1  This subroutine modifies liquid saturation and flow field 
!D1  for water table rise situation.
!D1
!***********************************************************************
!D2 
!D2  REVISION HISTORY
!D2
!D2 Initial implementation: 20-Oct-06, Shaoping Chu
!D2    to incorporate functionality of WTRISE V2.0 into FEHM
!D2
!D2  $Log:   /pvcs.config/fehm90/src/wtrise.f
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.3.5 Cell-based particle-tracking module
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!D4 Software Management Report: WTRISE V2.0, 
!D4   Document ID: 10537-SMR-2.0-00, MOL.20030519.0061
!D4
!***********************************************************************

      use comai
      use combi
      use comdi
      use comflow
      use compart
      implicit none
      
      integer i,itot,j,icheck,icheckf
      real*8, parameter :: sink=1.d+10 

      real*8,allocatable :: sf(:)
      real*8,allocatable :: sm(:)
      real*8,allocatable :: ff(:)
      real*8,allocatable :: fm(:)


      itot=nelm(neq+1)-(neq+1)
      

c     get sf,sm,from s; ff,fm from a_axy
      allocate(sf(neq),sm(neq))
      sf =0.0
      sm= 0.0
      do i=1,neq
         sf(i)=s(i)
         if (idpdp .ne. 0) then
            sm(i)=s(neq+i)
         end if
      end do


      allocate(ff(itot),fm(itot))
      ff= 0.0
      fm= 0.0
      do i=1,itot
         ff(i)=a_axy(i)
         if (idpdp .ne. 0) then
            fm(i)=a_axy(itot+i)
         end if
      end do


c     modify liquid saturation below water table
      icheck = 0
      do i=1,neq
         if (cord(i,3) .le. water_table) then
            sf(i)=1.
            sm(i)=1.
            icheck=icheck+1
         end if  
c     assign modified variable back to original array 
         s(i)=sf(i)
         if (idpdp .ne. 0) then
            s(neq+i)=sm(i)
         end if
      end do


c     go through  each connection and modify flow field
c     for those nodes that are below the water table
c     j is the node number
c     i denotes the index in the ncon(i.e. nelm) array starting 
c     with 1 instead of neq+1
      j=1
      icheckf=0
      do i=1,itot 
         if (i.eq.nelmdg(j)-neq-1) then
            if (cord(j,3).le.water_table) then
               ff(i)=sink
               fm(i)=sink
               icheckf=icheckf+1
            endif
            j=j+1
         end if
c     assign modified variable back to original array 
         a_axy(i)=ff(i)
         if (idpdp .ne. 0) then
            a_axy(itot+i)=fm(i)
         end if
      end do

      deallocate(sf,sm)
      deallocate(ff,fm)

      return
      end   

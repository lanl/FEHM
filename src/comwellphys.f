      module comwellphys 
!     comwellphys
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
!D1 PURPOSE
!D1
!D1 well modules
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 25-OCT-08, Programmer: george zyvoloski
!D2
!D2
!***********************************************************************

      integer  nwellphy
      integer, allocatable :: iwellpt(:)
      integer, allocatable :: iwellp(:)
      
      real*8, allocatable ::  rwell1f(:)
      real*8, allocatable ::  rwell2f(:)
      real*8, allocatable ::  rwell3f(:)
      real*8, allocatable ::  rwell4f(:)
      real*8, allocatable ::  rwell5f(:)
      real*8, allocatable ::  rwell6f(:)
      real*8, allocatable ::  rolsvf(:)
      real*8, allocatable ::  rovsvf(:) 
        
      real*8, allocatable ::  mdriftf(:)      
      real*8, allocatable ::  dmdriftp(:)
      real*8, allocatable ::  dmdrifte(:)
      real*8, allocatable ::  dmdriftw(:)                   

      real*8, allocatable ::  vel_m_well2(:)
      real*8, allocatable ::  rolsv(:)
      real*8, allocatable ::  rovsv(:)
      real*8, allocatable ::  vsupl(:)
      real*8, allocatable ::  vsupv(:)
      
      end module comwellphys  

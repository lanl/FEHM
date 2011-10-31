      module comtable
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
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
!D1 Global file for super-critical water.
!D1
!***********************************************************************

C     comtable.h      
C                                        ! phs 4/23/99
C--------------------------------------------------------------
cphs  Common block to hold the Lookup table and lookup table data
cphs  such as the increment in P and T and the num of points
cphs  in each direction 
c---------------------------------------------------------------

           integer inct, incp, numt, nump, tableFLAG

           integer :: iolookup
           real(8) ::  PP(65000,11)
           character(200) :: lookup_file = 'lookup.in'

           end module comtable

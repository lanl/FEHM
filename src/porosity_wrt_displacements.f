      subroutine porosity_wrt_displacements
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

!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To calculate porosity and its derivatives 
!D1 wrt displacements for coupled problems. 
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.00 
!D2 
!D2 Initial implementation: Date 14-Feb-2011, Programmer: sharad kelkar
!D2
!**********************************************************************
      use comai, only : neq
      use comdi, only : psini, ps
      use comsi, only : vol_strain, pore_factor

      implicit none

      integer i
      real*8 vstrain, vstrain_debug

c note: vol_strain < 0 for ocntraction
      do i = 1,neq
         vstrain =  vol_strain(i)
c.....................................
c         ps(i) = psini(i)*(1.+ vstrain)-pore_factor*vstrain
c note that ps calculated above is the actual porosity, but 
c ps calculated below is a factor so that 
c pore volume = ps()*initial bulk volume (=sx1d)
         ps(i) = psini(i) + pore_factor*vstrain
      enddo
      
      return
      
      end subroutine porosity_wrt_displacements

c.....................................................................

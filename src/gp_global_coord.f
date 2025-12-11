      subroutine gp_global_coord(i, j, x, y, z)
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
! 
! Returns the global (x,y,z) coordinates of the j-th gausspoint in the
! i-the element
! 
! Author : Sai Rapaka
!

      use comfem
      use comai, only: ns
      use combi, only: cord

      implicit none
      integer           :: i, j
      real*8            :: x, y, z

      integer           :: node, k

      x = 0.0d0
      y = 0.0d0
      z = 0.0d0

      do k=1,ns
        node = elnode(i, k)
        x = x + cord(node, 1)*Psi(i, j, k)
        y = y + cord(node, 2)*Psi(i, j, k)
        z = z + cord(node, 3)*Psi(i, j, k)
      enddo
        
      end subroutine gp_global_coord

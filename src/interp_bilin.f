      real*8 function interp_bilin (a1, a2, av, b1, b2, bv, t1, t2, t3,
     &     t4)
!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
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
!**********************************************************************
!
! PURPOSE
!
! To calculate particle time delay using type curves.
!
!**********************************************************************
!
! Initial implementation: 09-Nov-10, Programmer: Scott Painter
!
!**********************************************************************
! Return the interpolated time value
 
!argument list same as interp4

      implicit none 

      real*8 a1, a2, av, b1, b2, bv, t1, t2, t3, t4 
      real*8 u, t
      real*8 omt, omu

      t=(av-a1)/(a2-a1) 
      u=(bv-b1)/(b2-b1) 

! if linear in log-space then 
!      t=(dlog10 (av) - dlog10 (a1)) / (dlog10 (a2) - dlog10 (a1)) 
!      u=(dlog10 (bv) - dlog10 (b1)) / (dlog10 (b2) - dlog10 (b1)) 

      omt=1.0d0 - t 
      omu=1.0d0 - u 

      interp_bilin =omt*omu*t1 + t*omu*t3 + t*u*t4 + omt*u*t2 
 
      return 

      end function interp_bilin


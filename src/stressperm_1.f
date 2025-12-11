      subroutine stressperm_1(jpt)
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

c     perm model 1 - not implemented yet (default)
c     not needed unless fully coupled
      use comai
      use combi
      use comdi
      use comsi
      implicit none
      integer jpt

      if(allocated(rlxs)) then
         rlxs(jpt) = 1.0
         rlys(jpt) = 1.0 
         rlzs(jpt) = 1.0
         drlxs(jpt,1) = 0.0
         drlys(jpt,1) = 0.0
         drlzs(jpt,1) = 0.0
         drlxs(jpt,2) = 0.0
         drlys(jpt,2) = 0.0
         drlzs(jpt,2) = 0.0
         drlxs(jpt,3) = 0.0
         drlys(jpt,3) = 0.0
         drlzs(jpt,3) = 0.0
         drlxs(jpt,4) = 0.0
         drlys(jpt,4) = 0.0
         drlzs(jpt,4) = 0.0	  
      endif	 
      
      return
      end
c.....................................................................

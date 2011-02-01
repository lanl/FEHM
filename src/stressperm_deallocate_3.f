      subroutine deallocate_stressperm_3
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

c     
c     deallocate memory for stress derivatives for fully coupled solution
     
      use comai
      use combi
      use comdi
      use comdti
      use comsi
      implicit none

      if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &     allocated(rlxs))then       
         deallocate(rlxs,rlys,rlzs)
         deallocate(drlxs,drlys,drlzs)
         deallocate(idum_str1)
      endif
      
      return
      
      end
c......................................................................

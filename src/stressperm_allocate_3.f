      subroutine allocate_stressperm_3
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

c     only calculate for model 3, model 5, model 7 and model 8
c     initial setup calcs node neighbor information
      use comai
      use combi
      use comdi
      use comdti
      use comsi
      implicit none
c     
c     
      if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &     .not.allocated(rlxs))then  
         allocate(rlxs(n0))
         allocate(rlys(n0))
         allocate(rlzs(n0))
         allocate(drlxs(n0,4))
         allocate(drlys(n0,4))
         allocate(drlzs(n0,4))	  
         allocate (idum_str1(n0))
         idum_str1 = 0 
      endif
      return
      
      end
c......................................................................

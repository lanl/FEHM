      subroutine plasticctr(iflg)
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
! Driver routine for plastic computations using FEHM algorithm
! 
! Author : Sai Rapaka
!

      use comai
      use comsi

      implicit none

      integer iflg
      integer i,j,k,iterm,jterm
      logical null1
      character*10 macro1
      double precision, dimension(1:3, 1:3) :: tmp
      
      if(iflg.eq.initPlastic) then

         if(icnl.eq.0) then
           allocate(row(1:3,1:3,1:9))
           allocate(col(1:3,1:3,1:9))
           tmp(1,1) = 1; tmp(1,2) = 4; tmp(1,3) = 6;
           tmp(2,1) = 4; tmp(2,2) = 2; tmp(2,3) = 5;
           tmp(3,1) = 6; tmp(3,2) = 5; tmp(3,3) = 3;
           do i=1,3
             do j=1,3
               do k=1,9
                 iterm = floor((k-0.5)/3) + 1
                 jterm = k - (iterm-1)*3
                 row(i,j,k) = tmp(i,jterm)
                 col(i,j,k) = tmp(iterm,j)
!                 write(iout,*) i,j,k,row(i,j,k), col(i,j,k)
               enddo
             enddo
           enddo
         else
             allocate(row(1:2,1:2,1:4))
             allocate(col(1:2,1:2,1:4))
         endif

      elseif (iflg.eq.assemblePlastic) then
!        if(icnl.eq.0) then
!          call geneq_plastic_3D()
!        endif
      endif

      end subroutine plasticctr

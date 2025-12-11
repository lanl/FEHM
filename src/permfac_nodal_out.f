      subroutine permfac_nodal_out
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
c s kelkar sep 29 2011
c compute nodal averages of permfactor for output

      use comai, only: neq
      use combi, only: nelm
      use comfem, only:permfactor, permfac_out

      implicit none
      integer i,j,i_begin, i_end 
      integer permfac_out_count

      permfac_out_count = 0
      permfac_out = 1.0

      do i=1, neq
         i_begin = nelm(i) + 1
         i_end = nelm(i+1)
         permfac_out_count = 0
         do j=i_begin, i_end
            if(j.ne.i) then
               permfac_out(i,:) = permfac_out(i,:) + permfactor(j, :)
               permfac_out_count = permfac_out_count + 1
            endif
         enddo
         permfac_out(i,:) = permfac_out(i,:)/permfac_out_count
      enddo

      return

      endsubroutine permfac_nodal_out

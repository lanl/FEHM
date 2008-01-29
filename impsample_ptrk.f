      subroutine impsample_ptrk(model_flag,ith, num_part,
     1           rcoll_div,ret_weight, np_max)
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
!D1 Assign weights for importance sampling based on position in
!D1 CDF curve.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 27-Oct-06, Programmer: S. Kelkar
!D2
!D2 $Log:   /pvcs.config/fehm90/src/impsample.f_a  $
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.6 Streamline particle-tracking module 
!D3
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************
c Oct 25, 06 s kelkar
c Assing retardation factors using the following importance
c sampling scheme: equally divide the range of the retardation factor 
c into num_part intervals, and then assign the weight to the n-th 
c particle as weight(n)=CDF(n)-CDF(n-1) that will be needed for 
c calculating breakthrough curves. The weights are normalized to 1
c between the particles 1 thru num_part. The assignment of the index 
c is randomized at the end

      use comai, only : ierr, iptty
      use comsk, only : r_min, r_max, divs, flag_method

      implicit none

      integer np, model_flag,num_part, num_part_1, ith
      integer i, j, np_max
      integer :: impseed = 41293150

      real*8 cdf_min,cdf_max,dcdf,dret_ran

      save impseed

      real   ran_sp, randum
      real*8 ret_min,ret_max,ret,dret,cdf0,cdf1,sum_weight
      real*8 ret_log,dret_log,ret0,ret1
      real*8 ret_weight(np_max),rcoll_div(np_max)
      real*8 , allocatable :: temp_weight(:)
      real*8 , allocatable :: temp_rcoll(:)

      integer , allocatable :: index (:)
      allocate(index(num_part+1))
      allocate(temp_weight(num_part))
      allocate(temp_rcoll(num_part))

      if(num_part.le.0) then
         write(ierr,*)' num_part<=0 in impsample_ptrk. STOP'
         stop
      endif

      ret_min= r_min(divs(ith))
      ret_max= r_max(divs(ith))
      if(ret_max-ret_min.eq.0)then
         ret_weight = 1./num_part
         rcoll_div = ret_max
      else
         if(flag_method(ith).eq.0) then
c random linear sampling of the retardation factor
            dret=ret_max-ret_min
            do np=1,num_part
               randum=ran_sp(impseed)
               dret_ran=randum*dret
               ret=ret_min+dret_ran
               temp_weight(np)= 1.
               temp_rcoll(np) = ret
            enddo
         elseif(flag_method(ith).eq.1) then
c s kelkar April 11, 07
c     sampling with equal divisions in the CDF space, end-points
c     are not enforced, R value selected randomly for each 
c     subinterval in the log space
               call calculate_cdf(model_flag,ith,divs(ith),ret_min,
     1              cdf_min)
               call calculate_cdf(model_flag,ith,divs(ith),ret_max,
     1              cdf_max)
               dcdf=(cdf_max-cdf_min)/num_part
               ret0=ret_min
               cdf0=cdf_min
               do np=1,num_part
                  cdf1=cdf0+dcdf
                  call calculate_ret(model_flag,ith,divs(ith),cdf1,ret1)
                  dret_log=dlog10(ret1/ret0)
                  randum=ran_sp(impseed)
                  dret_ran=randum*dret_log
                  temp_weight(np)= 1.
                  temp_rcoll(np) = 10.**(dlog10(ret0)+dret_ran)
                  cdf0=cdf1
                  ret0=ret1
               enddo                 
         elseif(flag_method(ith).eq.2) then
c importance sampling with equal divisions in the log(R) space
c     in order to honor the max and min of R, the first and last
c     points are forced to these values, but this makes the two 
c     end intervals half the length of the remaining intervals.
c     np-1 equal divisions in the log space. use midpoints of
c     these intervals to calculate delta(CDF), except for the end 
c     points
            
            if(num_part.gt.1) then
c     dret_log is half the division size in log space
               dret_log=0.5*(dlog10(ret_max/ret_min))/(num_part-1)
c     calculate the first point
               ret=ret_min
               ret_log=dlog10(ret)
               temp_rcoll(1) = ret 
               call calculate_cdf(model_flag,ith,divs(ith),ret,cdf0)
               ret_log=ret_log+dret_log
               ret=10.**(ret_log)
               call calculate_cdf(model_flag,ith,divs(ith),ret,cdf1)
               temp_weight(1)=cdf1-cdf0
               cdf0=cdf1
               
               do np=2,num_part-1
                  ret_log=ret_log+dret_log
                  temp_rcoll(np) =10.**(ret_log) 
c     use the next mid-point for probability
                  ret_log=ret_log+dret_log
                  ret=10.**(ret_log)
                  call calculate_cdf(model_flag,ith,divs(ith),ret,cdf1)
                  temp_weight(np)=cdf1-cdf0
                  cdf0=cdf1
               enddo
c     calculate the last point
               ret_log=ret_log+dret_log
               ret=10.**(ret_log)
               temp_rcoll(num_part) =ret
               call calculate_cdf(model_flag,ith,divs(ith),ret,cdf1)
               temp_weight(num_part)=cdf1-cdf0
               
            else
c     cant sample less than 2 points, output an error message
               if (iptty .ne. 0) write(iptty,*) '# of samples in ',
     &              'impsample_uz table le 1 stop'
               write(ierr,*) '# of samples in impsample_uz table le ',
     &              '1 stop'
               stop
            endif
         endif
            
c randomize the assignment of the indices 1 thru num_part
         do i=1,num_part
            index(i)=i
         enddo
         index(num_part+1)=num_part
         num_part_1=num_part
         do i=1,num_part
            randum = ran_sp(impseed)
            np=1+int(randum*num_part_1)
            ret_weight(i)= temp_weight(index(np))
            rcoll_div(i)= temp_rcoll(index(np))
            num_part_1=num_part_1-1
            do j = np,num_part_1
               index(j) = index(j+1) 
            enddo 
         enddo
      endif
      
c make sure weights are normalized
      sum_weight=0
      do np=1,num_part
         sum_weight=sum_weight+ret_weight(np)
      enddo
      if(sum_weight.le.0.) then
         write(ierr,*)'error in impsample_uz. '
         write(ierr,*)'sum_weight.le.0.. stop.'
         if (iptty .ne. 0) then
            write(iptty,*)'error in impsample_uz. '
            write(iptty,*)'sum_weight.le.0. stop.'
            stop
         endif
      endif
      sum_weight=sum_weight/num_part
      do np=1,num_part
         ret_weight(np)=ret_weight(np)/sum_weight
      enddo

      deallocate(index)
      deallocate(temp_weight)
      deallocate(temp_rcoll)

      return

      end subroutine impsample_ptrk
c...........................................................

      subroutine impsample(model_flag)
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
c Assuming retardation factors using the following importance
c sampling scheme - equally divide the range of the retardation factor 
c into num_part intervals, and then assign the weight to the n-th 
c particle as weight(n)=CDF(n)-CDF(n-1) that will be needed for 
c calculating breakthrough curves

      use comdti
      use comai
      use comsptr
      use comsk

      implicit none

      integer i, j, np, model_flag, num_part_1
      integer :: impseed = 41293150
      real*8 ret_min,ret_max,ret,dret,cdf0,cdf1,sum_weight
      real*8 ret_log,dret_log,ret0,ret1

      save impseed
      real   ran_sp, randum

      real*8 , allocatable :: temp_weight(:)
      real*8 , allocatable :: temp_rcoll(:)

      integer , allocatable :: index (:)
      allocate(index(num_part+1))
      allocate(temp_rcoll(num_part))
      allocate(temp_weight(num_part))
 
      ret_min= r_min(1)
      ret_max= r_max(1)
      if(ret_max-ret_min.eq.0)then
         ret_weight = 1./num_part
         rcoll_div = ret_max
      else

c in order to honor the max and min of R, the first and last
c points are forced to these values, but this makes the two 
c end intervals half the length of the remaining intervals.
c np-1 equal divisions in the log space. use midpoints of
c these intervals to calculate delta(CDF), except for the end 
c points

         if(num_part.gt.1) then
c dret_log is half the division size in log space
            dret_log=0.5*(dlog10(ret_max/ret_min))/(num_part-1)
c calculate the first point
            ret=ret_min
            ret_log=dlog10(ret)
            temp_rcoll(1) = ret 
            call calculate_cdf(model_flag,1,1,ret,cdf0)
            ret_log=ret_log+dret_log
            ret=10.**(ret_log)
            call calculate_cdf(model_flag,1,1,ret,cdf1)
            temp_weight(1)=cdf1-cdf0
            cdf0=cdf1
            
            do np=2,num_part-1
               ret_log=ret_log+dret_log
               temp_rcoll(np) =10.**(ret_log) 
c use the next mid-point for probability
               ret_log=ret_log+dret_log
               ret=10.**(ret_log)
               call calculate_cdf(model_flag,1,1,ret,cdf1)
               temp_weight(np)=cdf1-cdf0
               cdf0=cdf1
            enddo
c calculate the last point
            ret_log=ret_log+dret_log
            ret=10.**(ret_log)
            temp_rcoll(num_part) =ret
            call calculate_cdf(model_flag,1,1,ret,cdf1)
            temp_weight(num_part)=cdf1-cdf0

         else
c cant sample less than 2 points, output an error message
            write(iptty,*)'# of smaples in impsample table le 1. stop'
            write(ierr,*)'# of smaples in impsample table le 1. stop'
            stop
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
            ret_weight(i,1)= temp_weight(index(np))
            rcoll_div(i,1)= temp_rcoll(index(np))
            num_part_1=num_part_1-1
            do j = np,num_part_1
               index(j) = index(j+1) 
            enddo 
         enddo
      endif
      
c make sure weights are normalized to total number of particles
      sum_weight=0
      do np=1,num_part
         sum_weight=sum_weight+ret_weight(np,1)
      enddo
      if(sum_weight.le.0.) then
         write(ierr,*)'error in R_importance_sample. stop.'
         if (iptty .ne. 0) 
     &        write(iptty,*)'error in R_importance_sample. stop.'
         stop
      endif
      sum_weight = sum_weight/num_part
      do np=1,num_part
         ret_weight(np,1)=ret_weight(np,1)/sum_weight
      enddo

      deallocate(index)
      deallocate(temp_weight)
      deallocate(temp_rcoll)

      return

      end subroutine impsample
c...........................................................

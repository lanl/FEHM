      subroutine find_closest_curve(flag, par1v, par2v, par3v,
     2     i, j, k, points, weight) 
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To determine the curve with the closest parameter values to the input
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 
!D2 Initial implementation: 15-DEC-02, Programmer: Bruce Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/find_closest_curve.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************
! Determine the transfer function curve with parameters closest to
! the input values

      use compfrac      
      use comai, only : ierr
      implicit none

      real*8 par1v, par2v, par3v
      real*8 pl1, pl2, pl3
      integer i, j, k
      integer i1
      real*8, allocatable, dimension(:) :: distance
      real*8, allocatable, dimension(:) :: p1
      real*8, allocatable, dimension(:) :: p2
      real*8, allocatable, dimension(:) :: p3
      integer, allocatable, dimension(:) :: iorder
      integer flag
      integer nparamp1, ifind, icol
      integer points(4)
      real weight(4), a(4,4), v(4,4), w(4), wttot
      integer ioindex(4), return_flag
      return_flag = 0

      allocate(distance(nump1), p1(nump1),p2(nump1),p3(nump1))
      allocate(iorder(nump1))
      distance = 0.
      iorder = 0
      p1 = 0.
      p2 = 0.
      p3 = 0.
      nparamp1=numparams+1

c     Load the parameter arrays and convert the appropriate
c     values to log values where specified
      if(log_flag(1).eq.1) then
         pl1 = dlog(par1v)/alog(normal_param(1))
         do i1 = 1, nump1
            p1(i1) = dlog(param1(i1))/alog(normal_param(1))
         end do
      else
         pl1 = par1v/normal_param(1)
         do i1 = 1, nump1
            p1(i1) = param1(i1)/normal_param(1)
         end do
      end if
      if(log_flag(2).eq.1) then
         pl2 = dlog(par2v)/alog(normal_param(2))
         do i1 = 1, nump1
            p2(i1) = dlog(param2(i1))/alog(normal_param(2))
         end do
      else
         pl2 = par2v/normal_param(2)
         do i1 = 1, nump1
            p2(i1) = param2(i1)/normal_param(2)
         end do
      end if
      if(numparams.gt.2) then
         if(log_flag(3).eq.1) then
            pl3 = dlog(par3v)/alog(normal_param(3))
            do i1 = 1, nump1
               p3(i1) = dlog(param3(i1))/alog(normal_param(3))
            end do
         else
            pl3 = par3v/normal_param(3)
            do i1 = 1, nump1
               p3(i1) = param3(i1)/normal_param(3)
            end do
         end if
      else
         pl3 = par3v
         do i1 = 1, nump1
            p3(i1) = pl3
         end do
      end if

c     Compute distance from the parameter to each curve

      do i1 = 1, nump1
         distance(i1) = sqrt((p1(i1)-pl1)**2+
     2        (p2(i1)-pl2)**2+
     3        (p3(i1)-pl3)**2)
      end do
c     Sort from lowest value to highest
      call sort_values(nump1, nump1, iorder, distance)


      if(flag.eq.0) then


c     Index of curve that has smallest deviation from the
c     input curve. j and k are 1 because of the data structure
c     i.e. curves go from 1 to nump1 in i, and j and k are 1 long
         
         i = iorder(1)
         j = 1
         k = 1

      else

c     Determine if any of the parameters are identical over all
c     points and replace with points somewhat farther away.
c     Always use closest two points (unless they are identical)
c     Start replacing 3rd and 4th points as necessary to get
c     some differences in each parameter (otherwise SVD
c     fails)

         call  filter_points

c     Set parameters to do nearest neighbor if an appropriate
c     set of points could not be found for svdcmp

         if(return_flag.ne.0) then
            if(numparams.eq.2) then
               weight(1) = 1.
               weight(2) = 0.
               weight(3) = 0.
               points(1) = iorder(1)
               points(2) = iorder(2)
               points(3) = iorder(3)
            elseif(numparams.gt.2) then
               weight(1) = 1.
               weight(2) = 0.
               weight(3) = 0.
               weight(4) = 0.
               points(1) = iorder(1)
               points(2) = iorder(2)
               points(3) = iorder(3)
               points(4) = iorder(4)
            end if
            return
         end if

c     Least squares fit to numparams+1 nearest curves
c     using singular value decomposition to find weights
c     to be used in interpolation
         a = 0.
         a(1,1) = pl1-p1(iorder(ioindex(1)))
         a(1,2) = pl1-p1(iorder(ioindex(2)))
         a(1,3) = pl1-p1(iorder(ioindex(3)))
         a(2,1) = pl2-p2(iorder(ioindex(1)))
         a(2,2) = pl2-p2(iorder(ioindex(2)))
         a(2,3) = pl2-p2(iorder(ioindex(3)))
         if(numparams.eq.3) then
            a(1,4) = pl1-p1(iorder(ioindex(4)))
            a(2,4) = pl2-p2(iorder(ioindex(4)))
            a(3,1) = pl3-p3(iorder(ioindex(1)))
            a(3,2) = pl3-p3(iorder(ioindex(2)))
            a(3,3) = pl3-p3(iorder(ioindex(3)))
            a(3,4) = pl3-p3(iorder(ioindex(4)))
         end if

c     SVD algorithm finds the weights


         call svdcmp_new(ierr,a,nparamp1,nparamp1,4,4,w,v)

         ifind = 0
         find_soln: do i = 1, nparamp1
            if(w(i).lt.1.e-6) then
               ifind = i
               exit find_soln
            end if
         end do find_soln

         if(ifind.eq.0) then
            if(numparams.eq.2) then
               weight(1) = 1.
               weight(2) = 0.
               weight(3) = 0.
               points(1) = iorder(1)
               points(2) = iorder(2)
               points(3) = iorder(3)
               write(ierr,*) 'Error in particle tracking interp.'
               write(ierr,*) 'Stopping'
            stop
            elseif(numparams.gt.2) then
               weight(1) = 1.
               weight(2) = 0.
               weight(3) = 0.
               weight(4) = 0.
               points(1) = iorder(1)
               points(2) = iorder(2)
               points(3) = iorder(3)
               points(4) = iorder(4)
            end if
            return
         end if
         
         weight = 0.
         wttot = 0.
         do icol = 1, nparamp1
            weight(icol) = v(icol,ifind)
            wttot = wttot + weight(icol)
         end do
         
c     Check and revert to nearest neighbor if poorly behaved

         do icol = 1, nparamp1
            if(abs(weight(icol)).gt.1e-10) then
c     weight_factor is a term that prevents the svd results to be used
c     when the weights themselves are much larger than the sum
c     of the weights. This is a situation in which numerical errors
c     can become a problem in the interpolation scheme. The code
c     reverts to the nearest neighbor if this happens.

               if(abs(wttot/weight(icol)).lt.weight_factor) then
                  wttot = 1.
                  weight(1) = 1.
                  weight(2) = 0.
                  weight(3) = 0.
                  if(numparams.gt.2) then
                     weight(4) = 0.
                  end if
               end if
            end if
         end do
         weight = weight/wttot
         points(1) = iorder(ioindex(1))
         points(2) = iorder(ioindex(2))
         points(3) = iorder(ioindex(3))
         if(numparams.gt.2) then
            points(4) = iorder(ioindex(4))
         end if
      end if

      deallocate(distance, iorder, p1, p2, p3)

      return
      contains

      subroutine filter_points
      implicit none
      real min_difference
      parameter (min_difference = 1.e-3)
      integer ipflag(3), ipts, jpar
      real pmatrix(3,4), denom, comp_value
      real max_p(3)
      logical done

      ioindex(1) = 1
      ioindex(2) = 2
      ioindex(3) = 3
      ioindex(4) = 4

      done = .false.
      do while(.not.done)
         ipflag = 0
c     Load matrix with inital values for checking the parameter values
         max_p = 0.
         do ipts = 1, nparamp1
            pmatrix(1,ipts) = p1(iorder(ioindex(ipts)))
            pmatrix(2,ipts) = p2(iorder(ioindex(ipts)))
            if(numparams.gt.2) then
               pmatrix(3,ipts) = p3(iorder(ioindex(ipts)))
            end if
         end do

         do jpar = 1, numparams

            inner_loop: do ipts = 2, nparamp1
               denom = max(abs(pmatrix(jpar,ipts)),
     2           abs(pmatrix(jpar,ipts-1)))
               if(denom.ne.0) then
                  comp_value = abs(pmatrix(jpar,ipts)-
     2                 pmatrix(jpar,ipts-1))/denom
                  if(comp_value.gt.min_difference) then
                     ipflag(jpar) = 1
                     exit inner_loop
                  else
                  end if
               else
               end if
            end do inner_loop

         end do

         done = .true.
         do jpar = 1, numparams
            if(ipflag(jpar).eq.0) done = .false.
         end do
         if(done) then
         end if
         if(.not.done) then
            if(numparams.eq.2) then
               if(ioindex(3).lt.nump1) then
                  ioindex(3) = ioindex(3) + 1
               else
                  return_flag = 1
                  done = .true.
               end if
            elseif(numparams.gt.2) then
               if(ioindex(4).lt.nump1) then
                  ioindex(3) = ioindex(3) + 1
                  ioindex(4) = ioindex(4) + 1
               else
                  return_flag = 1
                  done = .true.
               end if
            end if
         end if
      end do

      end subroutine filter_points
      end subroutine find_closest_curve

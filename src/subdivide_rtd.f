      subroutine subdivide_rtd(fractional_area)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Subdivide rtd curve for the mixing model.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 17-SEP-04,  Programmer: B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/subdivide_rtd.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
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


c     5. Subdivide f(t) curve

      use comrtd
      implicit none
      integer i
      integer j,isize2
      real*8 fractional_area
      real*8 tfinal, slope_exp, fcurve_final
      real*8 previous_time
      real*8 previous_rtd
      integer jcount
      real*8 delta_time
      real*8 delta_rtd
      real*8 ftotal, itotal
      real*8 previous_integrand
      integer maxdo
      integer subdiv
      parameter(subdiv = -1)

c     subdiv is now a flag that if negative, we are using
c     an equally spaced grid. For this option we want
c     the code to not subdivide the RTD curve, so
c     we are using a max function to set to 1 when it
c     is negative - BAR 6-15-2004

c     Allocate space for final subdivided rtd arrays

      nsubdiv = max(1,subdiv)*rtdcount
      if (.not. allocated(time_subdiv)) then
         allocate(time_subdiv(nsubdiv))
         allocate(rtd_subdiv(nsubdiv))
         allocate(capf(nsubdiv))
         allocate(cumi(nsubdiv))
      else
         isize2=size(time_subdiv,1)
         if(isize2.ne.nsubdiv) then
            deallocate(time_subdiv,rtd_subdiv,capf,cumi)
            allocate(time_subdiv(nsubdiv))
            allocate(rtd_subdiv(nsubdiv))
            allocate(capf(nsubdiv))
            allocate(cumi(nsubdiv))
         end if
      end if

      time_subdiv = 0.
      rtd_subdiv = 0.
      capf = 0.
      cumi = 0.


c     Loop through each original rtd point, subdivide

      previous_time = 0.
      previous_rtd = 0.
      jcount = 0
      if(fractional_area.lt.1.) then
         maxdo = rtdcount-25
      else
         maxdo = rtdcount
      end if
      do i = 1, maxdo
         delta_time = time_rtd(i)-previous_time
         delta_rtd = rtd(i)-previous_rtd
         do j = 1, max(1,subdiv)
            jcount = jcount + 1
            time_subdiv(jcount) = previous_time+ j*delta_time/
     2           max(1,subdiv)
            rtd_subdiv(jcount) = previous_rtd + j*delta_rtd/
     2           max(1,subdiv)
         end do
         previous_time = time_rtd(i)
         previous_rtd = rtd(i)
      end do

c     Loop through each subdivided point, integrate

      ftotal = 0.
      previous_integrand = 0.
      previous_time = 0.
      if(fractional_area.lt.1.) then
         maxdo = rtdcount-25
      else
         maxdo = rtdcount
      end if
      do i = 1, maxdo
         ftotal = ftotal + 0.5*(time_subdiv(i)-previous_time)
     2        *(previous_integrand + rtd_subdiv(i))
         capf(i) = ftotal
         previous_integrand = rtd_subdiv(i)
         previous_time = time_subdiv(i)
      end do

c     Loop through each subdivided point, normalize

      if(fractional_area.lt.1.) then
         maxdo = rtdcount-25
      else
         maxdo = rtdcount
      end if

      do i = 1, maxdo
         rtd_subdiv(i) = fractional_area*
     2        rtd_subdiv(i)/capf(maxdo)
         capf(i) = fractional_area*capf(i)/capf(maxdo)
      end do

      if(fractional_area.lt.1.) then
         tfinal = time_subdiv(nsubdiv-25)
         slope_exp = rtd_subdiv(nsubdiv-25)/(1.-fractional_area)
         delta_time = 1./(5.*slope_exp)
         fcurve_final = rtd_subdiv(nsubdiv-25)
c     Provide extrapolated curve
         do i = nsubdiv-24, nsubdiv
            time_subdiv(i) = time_subdiv(i-1) + delta_time
            rtd_subdiv(i) = fcurve_final*
     2           exp(-slope_exp*(time_subdiv(i)-tfinal))
            capf(i) = capf(i-1)+0.5*delta_time*
     2           (rtd_subdiv(i)+rtd_subdiv(i-1))
         end do

      end if
c     compute mean residence time

      ftotal = 0.
      previous_integrand = 0.
      previous_time = 0.
      if(fractional_area.lt.1.) then
         maxdo = rtdcount-25
      else
         maxdo = rtdcount
      end if
      do i = 1, maxdo
         ftotal = ftotal + 0.5*(time_subdiv(i)-previous_time)
     2        *(previous_integrand + time_subdiv(i)*rtd_subdiv(i))
         previous_integrand = time_subdiv(i)*rtd_subdiv(i)
         previous_time = time_subdiv(i)
      end do
      mean_residence_time = ftotal

c     compute cumulative age distribution

      itotal = 0.
      previous_integrand = 1./mean_residence_time
      previous_time = 0.
      if(fractional_area.lt.1.) then
         maxdo = rtdcount-25
      else
         maxdo = rtdcount
      end if

      do i = 1, maxdo
         itotal = itotal + 0.5*(time_subdiv(i)-previous_time)
     2        *(previous_integrand +
     3        (1.-capf(i))/mean_residence_time)
         cumi(i) = itotal
         previous_integrand = (1.-capf(i))/mean_residence_time
         previous_time = time_subdiv(i)
      end do

      end subroutine subdivide_rtd

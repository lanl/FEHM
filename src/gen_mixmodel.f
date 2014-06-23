      subroutine  gen_mixmodel
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Generate fluxes and flow source/sink flow for mixing model
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.22
!D2 Initial implementation: 17-SEP-04,  Programmer: B. Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inrestart.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai, only : neq, inpt, ierr, iptty, iocntl
      use combi, only : sx1
      use comrtd, only : mean_residence_time, mean_porosity, 
     2     mean_fl_density, maxmix_flag, rtdcount, time_rtd, 
     3     rtd
      implicit none
      integer i
      logical ex
      real*8 sumvol
      character*80 chdum
      character*100 locfilename
      real*8 denom_const
      parameter(denom_const = 2.50662827)
      real*8 mean, std_dev, delta_time
      real*8 meanlt, std_devlt, delta_ltime, logtime
      real*8 tmax, weight, weight_tot
      integer icurve, ncurves
      real*8 fractional_area,pivottime,fact0,factfinal
c     Read in parameters associated with the mixing reactor


      fractional_area = 1.0
      read (inpt, '(a80)') chdum
      if(chdum(1:3) .eq. 'min') then
         maxmix_flag = 1
      elseif(chdum(1:3) .eq. 'max') then
         maxmix_flag = -1
      end if

      read(inpt,*) mean_porosity, mean_fl_density

c     Read in RTD values or compute using normal or lognormal distributions

      read (inpt, '(a80)') chdum
      if(chdum(1:4) .eq. 'file') then

c     Read from file

      read(inpt, '(a100)') locfilename
      read(inpt,'(a80)') chdum
      if(chdum(1:4).eq.'sfac') then
         read(inpt,*) fractional_area
         if(fractional_area.gt.1..or.fractional_area.lt.0.) then
            fractional_area = 1.
         end if
         read(inpt,*) pivottime,fact0,factfinal
      else
         backspace inpt
         fractional_area = 1.
         pivottime = 1.
         fact0 = 1.
         factfinal = 1.
      end if
* Make sure file exists
         inquire (file = locfilename, exist = ex)
         if (.not. ex) then
            write (ierr, 200) locfilename
            if (iptty .gt. 0) write (iptty, 200) locfilename
            stop
 200        format ('ERROR nonexistant file', a100)
         end if

         open(iocntl, file = locfilename)

         call read_rtd(iocntl,fractional_area,
     2        pivottime,fact0,factfinal)
         close (iocntl)

      elseif(chdum(1:3) .eq. 'nor') then

c     Normal distribution of RTDs

         read(inpt,*) mean, std_dev
         rtdcount = 501
         allocate(time_rtd(rtdcount))
         allocate(rtd(rtdcount))
         time_rtd = 0.
         rtd = 0.

c    Compute rtd curve
         time_rtd(1) = 0.
         rtd(1) = 0.
         if(mean-4.*std_dev .gt. 0.) then
            delta_time = 8.*std_dev/(rtdcount-2)
            time_rtd(2) = mean-4.*std_dev
         else
            delta_time = mean/((rtdcount-1)/2)
            time_rtd(2) = delta_time
         end if
            rtd(2) = exp(-(time_rtd(2)-mean)**2/
     2           (2.*std_dev**2))/(denom_const*std_dev)
         do i = 3, rtdcount
            time_rtd(i)=time_rtd(i-1) + delta_time
            rtd(i) = exp(-(time_rtd(i)-mean)**2/
     2           (2.*std_dev**2))/(denom_const*std_dev)
         end do

      elseif(chdum(1:4) .eq. 'lnor') then

c     Lognormal distribution of RTDs

         read(inpt,*) meanlt, std_devlt
         rtdcount = 501
         if (.not. allocated(time_rtd)) then
            allocate(time_rtd(rtdcount))
            allocate(rtd(rtdcount))
         end if
         
         time_rtd = 0.
         rtd = 0.

c    Compute rtd curve
         time_rtd(1) = 0.
         rtd(1) = 0.
         delta_ltime = 8.*std_devlt/(rtdcount-2)
         logtime = meanlt-4.*std_devlt
         time_rtd(2) = exp(logtime)
            rtd(2) = exp(-(logtime-meanlt)**2/
     2        (2.*std_devlt**2))/
     3        (denom_const*std_devlt*time_rtd(2))
         do i = 3, rtdcount
            logtime = logtime + delta_ltime
            time_rtd(i)=exp(logtime)
            rtd(i) = exp(-(logtime-meanlt)**2/
     2           (2.*std_devlt**2))/
     3           (denom_const*std_devlt*time_rtd(i))
         end do


      elseif(chdum(1:4) .eq. 'mnor') then

c     Multiple Normal distribution of RTDs

         read(inpt,*) ncurves
         rtdcount = 1001
         allocate(time_rtd(rtdcount))
         allocate(rtd(rtdcount))
         time_rtd = 0.
         rtd = 0.
         tmax = 0.
         do icurve = 1, ncurves
            read(inpt,*) weight, mean, std_dev
            if(mean+4.*std_dev.gt.tmax) then
               tmax = mean+4.*std_dev
            end if
         end do
         do icurve = 1, ncurves
            backspace inpt
         end do
c     Compute rtd curve
         time_rtd(1) = 0.
         rtd(1) = 0.
         delta_time = tmax/(rtdcount-1)
         weight_tot = 0.
         do icurve=1,ncurves
            time_rtd(1) = 0.
            read(inpt,*) weight, mean, std_dev
            if(icurve.lt.ncurves) then
               weight_tot = weight_tot + weight
            else
               weight = 1.-weight_tot
            end if   
            do i = 2, rtdcount
               time_rtd(i)=time_rtd(i-1) + delta_time
               rtd(i) = rtd(i)+weight*exp(-(time_rtd(i)-mean)**2/
     2              (2.*std_dev**2))/(denom_const*std_dev)
            end do
         end do

      end if

c     Determine RTD arrays

      call subdivide_rtd(fractional_area)




      return
      end

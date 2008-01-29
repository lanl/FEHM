      module comrtd
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
!D1 Include file containing passed parameters and pointers related to
!D1 the mixing model
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 17-SEP-04,  Programmer: B. Robinson
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/comrtd.f_a  $
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

c     rtdcount - integer - number of points in original rtd curve
c     time_rtd - real*8, size rtdcount - array of times in the original
c        rtd curve
c     rtd - real*8, size rtdcount - array of f(t) points in the original
c        rtd curve
c     nsubdiv - integer - number of points in subdivided rtd curve
c     time_subdiv - real*8, size nsubdiv - array of times in the subdivided
c        rtd curve
c     rtd_subdiv - real*8, size nsubdiv - array of f(t) points in the subdivided
c        rtd curve
c     capf - real*8, size nsubdiv - array of cumulative rtd's for the subdivided
c        curve
c     cumi - real*8, size nsubdiv - array of cumulative internal age 
c        distributions
c     sx1 - real*8, size neq - array of finite volumes for each cell in the grid
c     mean_residence_time - mean residence time of reactor, used in 
c        pre-existing grid method for generating the grid
c     neq - integer - number of points in output 1D grid

      integer rtdcount
      real*8, allocatable :: time_rtd(:)
      real*8, allocatable :: rtd(:)
      integer nsubdiv
      real*8, allocatable :: time_subdiv(:)
      real*8, allocatable :: rtd_subdiv(:)
      real*8, allocatable :: capf(:)
      real*8, allocatable :: cumi(:)
      real*8 mean_residence_time, mean_porosity, mean_fl_density
      integer maxmix_flag

      end module comrtd

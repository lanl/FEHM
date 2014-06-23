       module compfrac
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
!D1 Data module for parallel fracture dispersion interpolation calculations.
!D1 
!***********************************************************************
!D2 
!D2 REVISION HISTORY
!D2 
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 Initial implementation: 28-Oct-02, Programmer: Z. Dash
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/compfrac.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:36   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 08:58:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
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
 
      implicit none
      save
 
      integer :: numparams !Number of parameters defining type curves
      integer :: nump1     !Number of param1 values or number of curves
      integer :: nump2     !Number of param2 values
      integer :: nump3     !Number of param3 values
      integer :: d4        !Number of possible fracture-matrix flow interactions
      integer :: curve_structure   ! if 1, free format structure of curves
      integer :: log_flag(100)  ! If 0, use value itself in interpolation
!   if 1, take log of value first
      real :: normal_param(100)
!Number of time/concentration pairs in each curve
      integer, allocatable, dimension(:,:,:,:) :: nump      
      integer, allocatable, dimension(:,:,:) :: itf_curve
      real, allocatable, dimension(:,:,:) :: wt_curve
      real weight_factor
!Maximum fracture flow fraction (upper limit) used in transfer function curves
      real*8 :: ffmax = 0.99
! Parameters to monitor range of values sampled in curves, first value
! will be max or min value in the curves, second value will be the
! max or min value sampled during calculations        
      real*8 :: sigma_low(2), sigma_high(2), omega_low(2), omega_high(2)
      real*8 :: par3_low(2), par3_high(2)
      real*8, allocatable, dimension(:) :: param1 
      real*8, allocatable, dimension(:) :: param2
      real*8, allocatable, dimension(:) :: param3
      real*8, allocatable, dimension(:,:,:,:,:) :: dtime
      real*8, allocatable, dimension(:,:,:,:,:) :: conc
      integer, allocatable, dimension(:,:,:) :: param_density

      logical :: pfrac_read = .false.
      integer :: ipartout = 0

      end module compfrac
 

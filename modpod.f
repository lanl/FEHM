      module modpod
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
!D1 Include file containing passed parameters related to POD basis functions.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20
!D2
!D2 Initial implementation: 2-Sep-04, Programmer: Z. Dash, B. Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comsk.f_a  $
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 N/A
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 None
!D4 
!***********************************************************************

      integer, allocatable :: ntouched(:)
      real*8, allocatable :: sx1(:)
      integer, allocatable :: irray(:,:)
      integer, allocatable :: nelm(:)
      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      real*8, allocatable :: z(:)
      real*8, allocatable :: concentration(:,:)
      real*8, allocatable :: time(:)
      real*8, allocatable :: c(:,:)
      real*8, allocatable :: eigenvec(:,:)
      real*8, allocatable :: eigenval(:)
      real*8, allocatable :: pod_bf(:,:)
      real*8, allocatable :: pod_bf1x(:,:)
      real*8, allocatable :: pod_bf1y(:,:)
      real*8, allocatable :: pod_bf1z(:,:)
      real*8, allocatable :: pod_bf2x(:,:)
      real*8, allocatable :: pod_bf2y(:,:)
      real*8, allocatable :: pod_bf2z(:,:)
      integer neq

      end module modpod

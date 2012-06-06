       module comfem
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
! Module that defines variables and arrays needed for simple, traditional,
! quadrature based integration for the mechanical equations
!
! Author : Sai Rapaka 
!
       integer                                 :: ifem
       integer                                 :: numgausspoints

       integer, allocatable                    :: elnode(:,:)
       real*8,  allocatable                    :: gpcord(:,:)
       real*8,  allocatable                    :: gpweight(:)

       real*8,  allocatable                    :: detJ(:,:)
       real*8,  allocatable                    :: Psi(:,:,:)
       real*8,  allocatable                    :: iPsi(:,:)
       real*8,  allocatable                    :: dPsidX(:,:,:)
       real*8,  allocatable                    :: dPsidY(:,:,:)
       real*8,  allocatable                    :: dPsidZ(:,:,:)

       real*8,  allocatable                    :: fem_stress(:,:,:)
       real*8,  allocatable                    :: fem_strain(:,:,:)
       real*8,  allocatable                    :: conv_strain(:,:,:)
       real*8,  allocatable                    :: conv_pstrain(:,:)
       real*8,  allocatable                    :: pstrain(:)

! Variables needed for permeability-stress modeling
       integer, allocatable                    :: edges(:, :)
       integer, allocatable                    :: edgeNum1(:, :)
       integer, allocatable                    :: edgeNum2(:, :)
       integer, allocatable                    :: numElems(:)
       integer, allocatable                    :: NodeElems(:)
       real*8,  allocatable                    :: permfactor(:,:)
       real*8,  allocatable                    :: permfactor_nodal(:)
       real*8,  allocatable                    :: permtmp(:,:)
c sep 29 2011 s kelkar. Calculate average nodal permfactor for output
       real*8, allocatable                     :: permfac_out(:,:)

       integer                                 :: flag_element_perm

       end module comfem

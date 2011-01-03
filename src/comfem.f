       module comfem

       ! Module that defines variables and arrays needed for simple, traditional,
       ! quadrature based integration for the mechanical equations

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

       end module comfem

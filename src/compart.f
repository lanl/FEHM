	module compart
!     compart
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
!D1  PURPOSE
!D1
!D1 Include file containing passed parameters and pointers related to
!D1 the particle tracking transport option.
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY
!D2
!D2 $Log:   /pvcs.config/fehm90/src/compart.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:57:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:44   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:25:58   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:57:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.4   Tue Jan 16 15:55:14 1996   zvd
!D2 Added prolog.
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  N/A
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4  
!D4  None
!D4  
!***********************************************************************

	integer, allocatable :: num_particles(:)
c bhl_12/8/05	
	integer, allocatable :: num_part_mem(:)
	real*8, allocatable :: cfraction(:)
c bhl_12/8/05
	integer maxlayers, max_particles, nreac
c bhl_11/3/06
	integer min_part
c bhl_11/3/06
	real, allocatable :: timeleft(:,:)
	real, allocatable :: partconc(:)
	real*8, allocatable :: nodepconc(:)
	integer, allocatable :: box(:,:)
	real*8, allocatable :: pcnsk(:)
	real, allocatable :: start_time(:,:)
	real, allocatable :: pconc(:),confactor(:)	
	real, allocatable :: newconfactor(:)
	real, allocatable :: frac_done(:,:)
	real, allocatable :: theta(:,:)
	real*8, allocatable:: flow_ot(:)
	real*8, allocatable:: flow_ottot(:)
	real*8, allocatable :: Rf(:,:)
	real*8, allocatable :: p(:)
	real*8, allocatable :: p_sf(:)
	real*8, allocatable :: mass(:)
	real, allocatable :: kd(:,:)
	real, allocatable :: rd_frac(:,:)
	real, allocatable :: kcoll(:,:)
	real, allocatable :: rcoll(:,:)
	real, allocatable :: fcoll(:,:)
	real, allocatable :: matrix_por(:)
	real, allocatable :: aperture(:)
c zvd 05/21/07 Added parameter for optional input of secondary spacing
	real, allocatable :: secondary(:)
	real, allocatable :: gamma_afm(:)
	real, allocatable :: sresidual(:)
	real, allocatable::  p_fraction(:)   
	integer  pout, prnt_rst,ipzone,ipzone1,ipzone2
	integer rseed, rseed_release
	integer :: ripfehm = 0
	logical restarting, ptrak
	integer, allocatable:: aidex(:), astep(:)
	integer, allocatable:: M_N_region(:), M_N_failed_nodes(:)
	integer, allocatable:: selected_in_region(:)
	integer, allocatable :: dispflag(:,:)
	integer, allocatable :: trak_type(:)
	integer, allocatable :: diffflag(:,:), sumdecayed(:)
	integer, allocatable :: sud_frind(:,:)
	integer, allocatable :: nprevd(:), nsegs(:), nsport(:)
	integer, allocatable :: lsport(:), ori(:), obj(:)
	integer, allocatable :: ptindex(:), dum_p(:), insnode(:)
	integer, allocatable :: idzone(:), ioconfactor(:)
	real*8, allocatable :: kfact(:), ivdt(:), tmsport(:)
	real*8, allocatable :: gmol(:), conftime(:)
	real*8, allocatable :: pcount(:,:), bconfactor(:), pconcmax(:)
        real*8, allocatable :: idcmax(:,:), idcavg(:,:), idflux(:)
        real*8, allocatable :: probsize(:,:), partsize(:,:)
        real*8, allocatable :: porsize(:,:)
        real*8, allocatable :: vg_poresize(:)
        integer, allocatable :: sizes(:)
        integer ncolsizes, ncolspec
	real*8, allocatable :: fracflow(:,:)
        integer, allocatable :: num_exit(:,:)
c***  spchu add water table rise 
	real*8  water_table 
	logical wtrise_flag
c***  add water table rise 
c zvd 11/09/06
	logical :: prnt_var(6)
c bhl_5/15/08
	real*8, allocatable :: bconf_sav(:,:)
	integer ipmbal
	integer:: ibox_offset = 0
c bhl_5/15/08

	end module compart

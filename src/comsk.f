      module comsk
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
!D1 Include file containing passed parameters and pointers related to
!D1 the streamline particle tracking option for OMR grids.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20
!D2
!D2 Initial implementation: 14-Jan-02, Programmer: S Kelkar
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comsk.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:00:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
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

c ...s kelkar  dec 13 01, new OMR stuff using least squares.....

c... april 12 04 3DOMR
      integer, allocatable :: oldnode(:)
      integer, allocatable :: oldnode_rand(:)
      integer, allocatable :: oldnode2(:)
c..........................................................

      integer iomrmax,komrmax,omr_nodes
      parameter ( komrmax = 30 )
      parameter ( iomrmax = 50 )

      integer node_count,komr_count(-3:3)
      integer iomr_count
      integer iomr_neighbour(iomrmax)
      integer isave_omr(-3:3,komrmax)

      integer, allocatable :: node_omr(:)

      integer iall_in,imod_in,iall_out,imod_out,igood_in,igood_out
      integer iall_in_eff,iall_out_eff
      integer listall_in(200),listmod_in(200),listmod_in_dir(200)
      integer listall_out(200),listmod_out(200),listmod_out_dir(200)
      integer colinear_in(200),colinear_out(200)
c...s kelkar  feb 9,04, 3DOMR......................
      integer, allocatable :: iboulist(:,:)
      integer, allocatable :: ipc_save(:)
      integer count_steps_max
      integer, allocatable :: count_steps(:)
c.......................................................

      real*8 mass_in,mass_out
      real*8 vin(3),vout(3),alam
      real*8 distance(200)
      real*8 cmodin(200,3),cmodout(200,3)
      real*8 cgbin(200,3),cgbout(200,3)
      real*8 csm(200,3)
      real*8 flow_all_in(200),flow_all_out(200)
      real*8 area_mod_in(200),area_mod_out(200)

      real*8 epsilonsk
      parameter (epsilonsk = 1.e-8)

      logical omr_flag, save_omr
c.......................................................

      real*8 deltawt
      real*8 smin
      real*8 sirr
      parameter (sirr = 0.01)
c smin is an optional input now
c      parameter (smin = 0.01)

c...s kelkar  aug 4,05, cliff nodes......................
      integer, allocatable :: istep_cases(:)
      logical cliff_flag
      integer , allocatable :: ncase(:)
      integer , allocatable :: ne(:,:)
      integer , allocatable :: node_flag1(:)
      integer , allocatable :: node_save(:)
      integer , allocatable :: list_case(:,:)
      integer , allocatable :: node_global(:)
      integer , allocatable :: node_case(:)
      integer , allocatable :: jtet(:)
      integer , allocatable :: jtetoff(:)
      integer , allocatable :: edge_node_tri(:,:)
      real*8, allocatable :: cords(:,:),tri_med(:,:),tri_norm(:,:)
      real*8, allocatable :: rnodecol(:),relemcol(:)
c....................................................... 
      logical corner_flag
      real*8 epsilon_corner
      integer, allocatable :: templist(:)
      integer, allocatable :: oldlist(:)
      integer, allocatable :: newlist(:)

c.......................................................

c Hari_Kelkar 01-Nov-06 Multispecies ..............................
c if multi-species option is invoked, read in Kds for each particle
c and rock type- 
      logical, allocatable :: flag_diversity(:)
      logical div_flag
      integer maxprobdivs
      integer iread_rcoll, ncoll_daugh
      integer total_colloids, total_rev, total_irrev, total_daughters
      integer, allocatable :: divs(:), flag_col_daughter(:)
      integer, allocatable :: irrevs(:), divs_d(:)
      integer, allocatable :: reves(:)
      integer, allocatable :: nprobdivs(:)
      real*8, allocatable :: probdiv(:,:)
      real*8, allocatable :: rcdiv(:,:)
      real*8, allocatable :: rcoll_div(:,:)
      logical, allocatable :: flag_col_irrev(:)
      real*8, allocatable :: mean_rcoll_reve(:)
      real*8, allocatable :: ret_weight(:,:)
      real*8, allocatable :: ret_weight_daughter(:,:)
      
      real*8, allocatable :: r_min(:)
      real*8, allocatable :: k_rev(:)
      real*8, allocatable :: r_max(:)
      real*8, allocatable :: slope_kf(:)
      real*8, allocatable :: cint_kf(:)
      
      integer, allocatable :: flag_log(:), flag_method(:)
      
c..................................................................
      

      end module comsk

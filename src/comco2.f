      module comco2
!     comco2
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
!D1 Global include file for array variables and pointers (FEHMN application).
!D1 Hydrate variables
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 06-June_2002  G. Zyvoloski           Initial implementation.
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!***********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   None
!D4
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block  Description
!D4
!D4   ***** COMMON Block fdd pointers and associated variables *****
!D4
!D4 Global Subprograms
!D4
!D4   None
!D4
!***********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!***********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!***********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 N/A
!D9
!***********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS None
!PS
!***********************************************************************

      real*8, allocatable ::  denco2i(:)
      real*8, allocatable ::  denco2h(:)
      real*8, allocatable ::  deneco2i(:)
      real*8, allocatable ::  deneco2h(:)
      real*8, allocatable ::  pflowco2(:) 
      real*8, allocatable ::  phico2(:) 
      real*8, allocatable ::  phoco2(:)
      real*8, allocatable ::  tco2(:) 
      real*8, allocatable ::  toco2(:)

      real*8, allocatable ::  skco2_tot(:)
      real*8, allocatable ::  skco2w(:)
      real*8, allocatable ::  qhco2(:)
      real*8, allocatable ::  qhco2w(:)
      real*8, allocatable ::  eskco2(:)
      real*8, allocatable ::  eflowco2(:)
      real*8, allocatable ::  wellco2(:)
      real*8, allocatable ::  rovfco2(:)
      real*8, allocatable ::  envfco2(:)
      real*8, allocatable ::  rolfco2(:)
      real*8, allocatable ::  enlfco2(:)
      real*8, allocatable ::  qhflxco2(:)

      real*8, allocatable :: dpcpw(:)
      real*8, allocatable :: dpcpg(:)
      real*8, allocatable ::  pcg(:)
      real*8, allocatable :: dpcgw(:)
      real*8, allocatable :: dpcgg(:)
      real*8, allocatable ::  dpcpco22(:)
      real*8, allocatable ::  dtpaco2(:)
      real*8, allocatable ::  dtpaeco2(:)

      integer, allocatable ::  kaco2(:)
      
      real*8, allocatable ::  dmycf(:)
      real*8, allocatable ::  dmyaf(:)
      real*8, allocatable ::  deycf(:)
      real*8, allocatable ::  deyaf(:)
c      real*8, allocatable ::  dmxcf(:)
c      real*8, allocatable ::  dexcf(:)

      real*8, allocatable ::  skco2(:)
      real*8, allocatable ::  dskco2w(:)
      real*8, allocatable ::  dskco2ya(:)
      real*8, allocatable ::  dskco2yc(:)
c      real*8, allocatable ::  dskco2xc(:)

      real*8, allocatable ::  dqw(:)
      real*8, allocatable ::  dqya(:)
      real*8, allocatable ::  dqyc(:)
      real*8, allocatable ::  deqya(:)
      real*8, allocatable ::  deqyc(:)
c      real*8, allocatable ::  deqw(:)
c      real*8, allocatable ::  dqxc(:)
c      real*8, allocatable ::  deqxc(:)
      real*8, allocatable ::  dqhw(:)
      real*8, allocatable ::  dqhya(:)
      real*8, allocatable ::  dqhyc(:)
c      real*8, allocatable ::  dqhxc(:)

      real*8, allocatable :: rl_l(:)
      real*8, allocatable :: drl_lw(:)
      real*8, allocatable :: drl_lg(:)
      real*8, allocatable :: rl_v(:)
      real*8, allocatable :: drl_vw(:)
      real*8, allocatable :: drl_vg(:)
      real*8, allocatable :: rl_w(:)
      real*8, allocatable :: drl_ww(:)
      real*8, allocatable :: drl_wg(:)
      real*8, allocatable ::  a44i(:)

      real*8, allocatable :: csalt(:)

      real*8, allocatable :: wat_prop(:)
      real*8, allocatable :: co2_prop(:)
      real*8, allocatable :: dmol(:)

      real*8, allocatable :: diw(:)
	real*8, allocatable :: diwc(:)
      real*8, allocatable :: diwp(:)
      real*8, allocatable :: diwe(:)
      real*8, allocatable :: diww(:)
      real*8, allocatable :: diwyc(:)
      real*8, allocatable :: diwya(:)
      real*8, allocatable :: divyc(:)
      real*8, allocatable :: divya(:)  
      real*8, allocatable :: dilyc(:)
      real*8, allocatable :: dilya(:)
c      real*8, allocatable :: diwxc(:)
c      real*8, allocatable :: divxc(:) 
c      real*8, allocatable :: dilxc(:)

      real*8, allocatable :: dtpsc(:)
      real*8, allocatable ::  diltrac(:)
      real*8, allocatable :: stowat(:)
c     integer, allocatable :: idco2(:)
      real*8, allocatable :: xc(:)
      real*8, allocatable :: xw(:)
      real*8, allocatable :: xa(:)
      real*8, allocatable :: yc(:)
      real*8, allocatable :: yw(:)
      real*8, allocatable :: ya(:)
      real*8, allocatable :: fw(:)
      real*8, allocatable :: fl(:)
      real*8, allocatable :: fg(:)
      real*8, allocatable :: fow(:)
      real*8, allocatable :: fol(:)
      real*8, allocatable :: fog(:)
      real*8, allocatable :: flowco2s(:)
      real*8, allocatable :: fw_tmp(:)
      real*8, allocatable :: fl_tmp(:)
      real*8, allocatable :: fg_tmp(:)
      integer, allocatable ::  inico2flg(:)

      real*8, allocatable :: xoc(:)
      real*8, allocatable :: xow(:)
c      real*8, allocatable :: xoa(:)
      real*8, allocatable :: yoc(:)
      real*8, allocatable :: yow(:)
c      real*8, allocatable :: yoa(:)

      real*8, allocatable :: c_axy(:)
      real*8, allocatable :: c_vxy(:)

      real*8, allocatable :: diff(:)
      real*8, allocatable :: tortco2(:)
      real*8, allocatable :: con_prop(:)
      integer, allocatable :: ico2dis(:)
      integer, allocatable :: ico2diso(:)

	real*8, allocatable :: strd_arr(:)

      real*8  amco20,aeco20,amco2,aeco2
      real*8  qco2_in,qeco2_in,balco2,baleco2
      real*8  qco2,qeco2,qco2ts,qeco2ts
      real*8  qco2ts_in,qeco2ts_in
c     real*8  aihyd                    

      real*8  amco2hyd0,aeco2hyd0,amco2hyd,aeco2hyd
      real*8  qco2hyd_in,qeco2hyd_in,balco2hyd,baleco2hyd
      real*8  qco2hyd,qeco2hyd,qco2hydts,qeco2hydts
      real*8  qco2hydts_in,qeco2hydts_in

      integer idof_co2,ibrine, iprtype, ico2sol, icarb
      integer ico2diff_flg, ico2prop_flg, iwatdis

      end module comco2

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

                logical strong
      real*8, allocatable ::  denco2i(:)
      real*8, allocatable ::  denco2h(:)
      real*8, allocatable ::  deneco2i(:)
      real*8, allocatable ::  deneco2h(:)
      real*8, allocatable ::  pflowco2(:) 
      real*8, allocatable ::  phico2(:) 
      real*8, allocatable ::  phoco2(:)
      real*8, allocatable ::  tco2(:) 
      real*8, allocatable ::  toco2(:)
      real*8, allocatable ::  phih2o(:)
      real*8, allocatable ::  phoh2o(:)
      real*8, allocatable ::  dpsattf(:)

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
      real*8, allocatable ::  rovfmix(:)
      real*8, allocatable ::  drovfmixp(:)
      real*8, allocatable ::  drovfmixe(:)
      real*8, allocatable ::  envfmix(:)
      real*8, allocatable ::  denvfmixp(:)
      real*8, allocatable ::  denvfmixe(:)
      real*8, allocatable ::  rol_co2(:)
      real*8, allocatable ::  drol_co2p(:)
      real*8, allocatable ::  drol_co2e(:)
      real*8, allocatable ::  enl_co2(:)
      real*8, allocatable ::  denl_co2p(:)
      real*8, allocatable ::  denl_co2e(:)

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
      real*8, allocatable :: yc_tmp(:)
      real*8, allocatable :: xco2_vapf(:)
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

      real*8, allocatable :: rate(:)
      real*8  amco20,aeco20,amco2,aeco2
      real*8  qco2_in,qeco2_in,balco2,baleco2
      real*8  qco2,qeco2,qco2ts,qeco2ts
      real*8  qco2ts_in,qeco2ts_in
c     real*8  aihyd                    

      real*8  amco2hyd0,aeco2hyd0,amco2hyd,aeco2hyd
      real*8  qco2hyd_in,qeco2hyd_in,balco2hyd,baleco2hyd
      real*8  qco2hyd,qeco2hyd,qco2hydts,qeco2hydts
      real*8  qco2hydts_in,qeco2hydts_in
      real*8  phi_vap_lmt
      parameter (phi_vap_lmt = 1.d-5)

      integer idof_co2,ibrine, iprtype, ico2sol, icarb, carbon_tracer
      integer ico2diff_flg, ico2prop_flg, iwatdis
c impes-like flow terms  
      integer imped_ex    
      real*8, allocatable :: permsd1_sv_w(:)
      real*8, allocatable :: permsd1_sv_co2(:)
      
      real*8 :: permsd11 = 0., dprmp1 = 0., dprmt1 = 0., dprmw1 = 0.
      real*8 :: permsd12 = 0., dprmp2 = 0., dprmt2 = 0., dprmw2 = 0.
      real*8 :: permsd13 = 0., dprmp3 = 0., dprmt3 = 0., dprmw3 = 0.

      real*8 frac_cl,frac_cg,frac_c,frac_w, yco2,ywat,yair,xco2,xwat
      real*8 rol_h2o,rol_d,emw,drol_dp,drol_dt,drol_dyc,drol_dya,xair
      real*8 drolyc, drolya, roa, droadp, droadt, ena, denadt,denadp
      real*8 dvisadp, denvt, denvp, visca
      real*8 dvisadt, pw, rlw, drlww, drlwg, drlwp, drlwt
      real*8 rll, drllw, drllg, drllp, drllt
      real*8 enw,enl,dhlt,dhlp,dhvt,dhvp
      real*8 drolw, denlt, denlp, denlw, denlya, denlyc
      real*8 visl, dvislya, dvislyc, dvislw
      real*8 row, drowp, drowt, drowya, drowyc, denwp, denwt
      real*8 visw, dviswp, dviswt
      real*8 dprmp, dprmt, dprmw, dprmyc, dprmya
      real*8 damyc, damya, daeyc, daeya
      real*8 dhprdyc, dhprdya, rlv, drlvw, drlvg, drlvp, drlvt
      real*8 drovw, drovya, drovyc, denvw, denvya, denvyc, dvisvya
      real*8 dvisvyc, dvisvw, demwyc, denwyc, denwya
      real*8 vpartial,dvpardt
      real*8 xs, dxsw, dxsg, denwfw, denwfg, s1, s2, ds1dw, ds1dg
c new variables
      real*8 enx,denxp,denxe,denxyc,denxya, vis_tol, xtol
      real*8 pvap_h2o,pvap_co2,pcl0,drow_vp,denv_vapt
      real*8 dxwat_vap,drov_co2p
      real*8 dxco2p,env_h2o,env_co2,denv_co2p,denv_co2t
      real*8 drov_co2w,drov_co2yc,drov_co2ya
      real*8 dvisv_co2w,dvisv_co2yc,dvisv_c2ya
      real*8 denv_co2w,denv_co2yc,denv_co2ya
      real*8 visv_co2,dvisv_co2p,dvisv_co2ya,dvisv_co2t
      real*8 denv_vap,frac_g,drow_vt,dxwat_vat,dxco2t,dyco2p
      real*8 dpvap_h2op,drov_co2t,drow_vtw,drow_vpw,dyco2t
      real*8 daml,damg,dumc(9)
      real*8 daet1,daet2,dael,daeg
      real*8 denv_h2op,denv_h2ot
      real*8 visv_h2o, dvisv_h2op, dvisv_h2ot 
      real*8 dxco2_vap,dxco2_vat
      real*8 drov_h2ot,drov_h2op
      
      real*8 rov, d_rov_1, d_rov_2 ,d_rov_t
      real*8 d_rol_1, d_rol_2 ,d_rol_t
      real*8 d_row_1, d_row_2 ,d_row_t
      real*8 env, d_env_1, d_env_2, d_env_t
      real*8 d_enw_1, d_enw_2, d_enw_t
      real*8 d_enl_1, d_enl_2, d_enl_t
      real*8 visv, d_visv_1, d_visv_2, d_visv_t
      real*8 erock, d_erock_1, d_erock_2, d_erock_t
      real*8 d_am_2, d_ae_2
      real*8 d_por_2, d_pvap_h2o_t
      real*8  rov_h2o, d_rov_h2o_1, d_rov_h2o_2, d_rov_h2o_t 
      real*8  rov_co2, d_rov_co2_1, d_rov_co2_2, d_rov_co2_t 
      real*8 xwat_vap, d_xwat_vap_1, d_xwat_vap_2, d_xwat_vap_t
      real*8 xco2_vap, d_xco2_vap_1, d_xco2_vap_2, d_xco2_vap_t
      real*8  d_env_h2o_1, d_env_h2o_2, d_env_h2o_t 
      real*8  d_env_co2_1, d_env_co2_2, d_env_co2_t 
      real*8  d_env_vap_1, d_env_vap_2, d_env_vap_t     
      real*8 d_visv_h2o_1, d_visv_h2o_2, d_visv_h2o_t
      real*8 d_visv_co2_1, d_visv_co2_2, d_visv_co2_t
      
      real*8 permsd,dtsatpco2
      real*8 drovt,drovp,dvisvt,dvisvp
      real*8 rol,drolt,drolp,dvislt,dvislp
      parameter (xtol = 1.d-16) 
      
      end module comco2

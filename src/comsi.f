      module comsi
!***********************************************************************
! Copyright 2008 Los Alamos National Security, LLC  All rights reserved
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
!D1
!D1 PURPOSE
!D1
!D1 Global include file for array variables (FEHMN application).
!D1 Coupled Fluid Stress application
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

      integer flag_permmodel
      integer incremental_shear_permmodel

      integer icycs,ipsps,iinp,minks,inst,ipchng,nnx,ihms,nvfc 
      integer nomass,noheat,noydis,localx,localy,ifrac,idof_stress
      integer ibp_stress, ibodyforce, ipermstr, initcalc, istrshis 
      integer ispmd,ipermstr1,ipermstr2,ipermstr3
      integer ipermstr4,ipermstr5,ipermstr6,ipermstr7,ipermstr8
      integer ipermstr11,ipermstr21,ipermstr22,ipermstr222
      integer ipermstr23, ipermstr24, ipermstr25
      integer ipermstr91, ipermstr31, ipermstr100
      integer istresscall,ilitho,iad_strs,istresspor,idisp_rel
      integer cnum_stress, ilithgrad, permfile
      
      integer isNonlinear, Nonlin_model_flag, flag_principal
c gaz 031017    
      integer isbiotNonLin
      integer flag_excess_shear
c gaz 061123
      integer cnum_sv
     
      real*8 friction_out,strength_out, pp_fac_out
      real*8 daystr,fpor,fric,wo,fwght,fupwt,fdnwt,tol_stress
      real*8 tptch,pchmin,pchmax,tchmin,tchmax,coftol,bp_stress
      real*8 bpx,bpy,bpz,tol_stress1,abs_tol_stress,bp_update
      real*8 vol_tot_change, daystress, timestress, timestress0
      real*8 strx_min,stry_min,strz_min,perx_m,pery_m,perz_m
      real*8 e10_facx,e10_facy,e10_facz,elev_high,pres_elev
      real*8 tenslie_str_fac,area_str_zone
      parameter (tensile_str_fac = 10.d0)
      integer, allocatable :: kr(:,:)
      integer, allocatable :: npbn(:)
      integer, allocatable :: nvfcl(:)
      integer, allocatable :: ispmt(:)
      integer, allocatable :: ispm(:)
      integer, allocatable :: nskw_stress(:,:)
      integer, allocatable :: frac_flg(:)

      real*8, allocatable ::  elastic_mod(:)
      real*8, allocatable ::  poisson(:)
      real*8, allocatable ::  e1(:)
      real*8, allocatable ::  e2(:)
      real*8, allocatable ::  e3(:)
c s kelkar 12/6/09 axisymmetric anisotropy
      logical stress_anisotropy_in,stress_anisotropy_use
        real*8, allocatable ::  elastic_mod_t(:)
        real*8, allocatable ::  poisson_t(:)
        real*8, allocatable ::  e4(:)
        real*8, allocatable ::  ezz(:)
        real*8, allocatable ::  shearmod_t(:)
c d dempsey 03/11/14 bodyforce vector
      real*8, allocatable ::  bodyforce_x(:)
      real*8, allocatable ::  bodyforce_y(:)
      real*8, allocatable ::  bodyforce_z(:)

c......................................
c sai - variables needed for plasticity
       integer iPlastic, initPlastic, assemblePlastic
       integer NumPlasticModels

       integer, allocatable :: row(:,:,:), col(:,:,:)
 
       integer, allocatable::  plasticModel(:)
       integer, allocatable::  modelNumber(:)

       integer, allocatable::  isPlastic(:,:)
       real*8, allocatable ::  plastic_strain(:,:)
       real*8, allocatable ::  plasticParam1(:)
       real*8, allocatable ::  plasticParam2(:)
       
       real*8, allocatable ::  delta_u(:)
       real*8, allocatable ::  delta_v(:)
       real*8, allocatable ::  delta_w(:)

       real*8, allocatable ::  yield_stress(:)
       real*8, allocatable ::  hardening(:)
 
       real*8, allocatable ::  strain_xx(:)
       real*8, allocatable ::  strain_yy(:)
       real*8, allocatable ::  strain_zz(:)
       real*8, allocatable ::  strain_xy(:)
       real*8, allocatable ::  strain_yz(:)
       real*8, allocatable ::  strain_zx(:)
 
c......................................
c k yoshioka 3/2/10 orthotropy
      logical stress_orthotropy
      real*8, allocatable :: elastic_mod_2(:)
      real*8, allocatable :: elastic_mod_3(:)
      real*8, allocatable :: poisson_13(:)
      real*8, allocatable :: poisson_23(:)
      real*8, allocatable :: shearmod_12(:)
      real*8, allocatable :: shearmod_13(:)
      real*8, allocatable :: shearmod_23(:)
      real*8, allocatable :: e5(:)
      real*8, allocatable :: e6(:)
      real*8, allocatable :: bulk1(:)
      real*8, allocatable :: bulk2(:)
      real*8, allocatable :: bulk3(:)
c.........................................      
      real*8, allocatable ::  bulk(:)
      real*8, allocatable ::  alp(:)
c gaz 042917 added initial thermal and bulk modulus      
      real*8, allocatable ::  bulk0(:)
      real*8, allocatable ::  alp0(:)
      real*8, allocatable ::  du(:) 
      real*8, allocatable ::  dv(:) 
      real*8, allocatable ::  dw(:)
      real*8, allocatable ::  duo(:) 
      real*8, allocatable ::  dvo(:) 
      real*8, allocatable ::  dwo(:)
      real*8, allocatable ::  du_ini(:) 
      real*8, allocatable ::  dv_ini(:) 
      real*8, allocatable ::  dw_ini(:)
      real*8, allocatable ::  du_tot(:) 
      real*8, allocatable ::  dv_tot(:) 
      real*8, allocatable ::  dw_tot(:)
      real*8, allocatable ::  sxx(:) 
      real*8, allocatable ::  syy(:) 
      real*8, allocatable ::  szz(:)
      real*8, allocatable ::  epsxx(:) 
      real*8, allocatable ::  epsyy(:) 
      real*8, allocatable ::  epszz(:)
      real*8, allocatable ::  sxy(:) 
      real*8, allocatable ::  syz(:) 
      real*8, allocatable ::  szx(:)
      real*8, allocatable ::  epsxy(:) 
      real*8, allocatable ::  epsyz(:) 
      real*8, allocatable ::  epszx(:)
      real*8, allocatable ::  str_x(:) 
      real*8, allocatable ::  str_y(:)
      real*8, allocatable ::  str_z(:)
      real*8, allocatable ::  str_x0(:)
      real*8, allocatable ::  str_y0(:)
      real*8, allocatable ::  str_z0(:)
      logical residual_stress
c Bai model      
      real*8, allocatable ::  estr_x0(:)
      real*8, allocatable ::  estr_y0(:)
      real*8, allocatable ::  estr_z0(:)
      real*8, allocatable ::  estr_xy0(:)
      real*8, allocatable ::  estr_xz0(:)
      real*8, allocatable ::  estr_yz0(:)
c Failure + Bai displacement
      real*8, allocatable ::  es_f_x0(:,:)  
      real*8, allocatable ::  es_f_y0(:,:)
      real*8, allocatable ::  es_f_z0(:,:)
      real*8, allocatable ::  s_f_xy0(:,:)
      real*8, allocatable ::  s_f_xz0(:,:)
      real*8, allocatable ::  s_f_yz0(:,:) 
      real*8, allocatable ::  frc_zen(:,:)
      real*8, allocatable ::  frc_azm(:,:)
cccccccccccccccccccccccccccccccccccccccccccccccheck
      real*8,  allocatable ::  check(:)            
c     ___________________________________________      
      real*8, allocatable ::  str_xy(:) 
      real*8, allocatable ::  str_xz(:)
      real*8, allocatable ::  str_yz(:)
      real*8, allocatable ::  str_xy0(:)
      real*8, allocatable ::  str_xz0(:)
      real*8, allocatable ::  str_yz0(:)
      real*8, allocatable ::  asxx(:) 
      real*8, allocatable ::  asyy(:)
      real*8, allocatable ::  disp(:,:)
      real*8, allocatable ::  forc(:,:)
      real*8, allocatable ::  phistr(:) 
      real*8, allocatable ::  tmpstr(:)
      real*8, allocatable ::  dluf(:) 
      real*8, allocatable ::  dlvf(:)
      real*8, allocatable ::  shp(:) 
      real*8, allocatable ::  sht(:)
      real*8, allocatable ::  uxtp(:) 
      real*8, allocatable ::  uytp(:)
      real*8, allocatable ::  uztp(:)
      real*8, allocatable ::  dmadv(:)
      real*8, allocatable ::  deadv(:)
      real*8, allocatable ::  drluf(:)
      real*8, allocatable ::  drlvf(:)
      real*8, allocatable ::  rlxs(:)
      real*8, allocatable ::  rlys(:)
      real*8, allocatable ::  rlzs(:)
      real*8, allocatable ::  drlxs(:,:)
      real*8, allocatable ::  drlys(:,:)
      real*8, allocatable ::  drlzs(:,:)
      real*8, allocatable ::  drlxv(:,:)
      real*8, allocatable ::  drlzv(:,:)
      real*8, allocatable ::  dflv(:)
      real*8, allocatable ::  dvfdv(:,:)
      real*8, allocatable ::  vol_strain(:)
      real*8, allocatable ::  vol_strain0(:)
      real*8, allocatable ::  dvol_strainu(:)
      real*8, allocatable ::  dvol_strainv(:)
      real*8, allocatable ::  dvol_strainw(:)
      real*8, allocatable ::  vol_temp(:)
      real*8, allocatable ::  pnx0(:)
      real*8, allocatable ::  pny0(:)
      real*8, allocatable ::  pnz0(:)
      real*8, allocatable ::  e10(:)
      real*8, allocatable ::  e20(:)
      real*8, allocatable ::  e30(:)
      
      integer, allocatable ::   itstress(:)


      real*8 ftoll
      real*8, allocatable ::  af(:,:)
      real*8, allocatable ::  cf(:,:)
      real*8, allocatable ::  stifi(:,:) 
      real*8, allocatable ::  sigfi(:,:)    
      real*8, allocatable ::  vm(:,:)
      real*8, allocatable ::  dlvy(:) 
      real*8, allocatable ::  dlux(:) 
      real*8, allocatable ::  dsx(:) 
      real*8, allocatable ::  dsy(:) 

c     common/fracturex/dla(n0),dlb(n0),dmpm(n0),dmtm(n0),
c     $    dmu(n0), dmv(n0),depm(n0),detm(n0),deu(n0),dev(n0),
c     $    akkx(n0),akky(n0), aiio(n0),biio(n0),
c     $    ffx(n0),ffy(n0)

      real*8, allocatable ::  dla(:) 
      real*8, allocatable ::  dlb(:) 
      real*8, allocatable ::  dmpm(:) 
      real*8, allocatable ::  dmtm(:) 
      real*8, allocatable ::  dmu(:) 
      real*8, allocatable ::  dmv(:) 
      real*8, allocatable ::  depm(:) 
      real*8, allocatable ::  detm(:)
      real*8, allocatable ::  deu(:) 
      real*8, allocatable ::  dev(:)  
      real*8, allocatable ::  akkx(:) 
      real*8, allocatable ::  akky(:)
      real*8, allocatable ::  aiio(:) 
      real*8, allocatable ::  biio(:)  
      real*8, allocatable ::  ffx(:) 
      real*8, allocatable ::  ffy(:)  

c     dimension a12mpf(n0),a12mef(n0),a12eef(n0)

      real*8, allocatable ::  a12mpf(:)  
      real*8, allocatable ::  a12mef(:) 
      real*8, allocatable ::  a12eef(:)  

      integer, allocatable ::  icxuani(:) 
      integer, allocatable ::  icxvani(:) 
      integer, allocatable ::  icxwani(:)
      integer, allocatable ::  icyuani(:) 
      integer, allocatable ::  icyvani(:) 
      integer, allocatable ::  icywani(:)
      integer, allocatable ::  iczuani(:) 
      integer, allocatable ::  iczvani(:) 			   
      integer, allocatable ::  iczwani(:)

      integer, allocatable ::  ncon_xu1(:) 
      integer, allocatable ::  ncon_xv1(:)
      integer, allocatable ::  ncon_xw1(:)
      integer, allocatable ::  ncon_yu1(:)
      integer, allocatable ::  ncon_yv1(:)
      integer, allocatable ::  ncon_yw1(:)
      integer, allocatable ::  ncon_zu1(:)
      integer, allocatable ::  ncon_zv1(:)
      integer, allocatable ::  ncon_zw1(:)
      
      integer, allocatable ::  ncon_xu2(:) 
      integer, allocatable ::  ncon_xv2(:)
      integer, allocatable ::  ncon_xw2(:)
      integer, allocatable ::  ncon_yu2(:)
      integer, allocatable ::  ncon_yv2(:)
      integer, allocatable ::  ncon_yw2(:)
      integer, allocatable ::  ncon_zu2(:)
      integer, allocatable ::  ncon_zv2(:)
      integer, allocatable ::  ncon_zw2(:)      
      
      integer, allocatable ::  istrws(:)   
      
      real*8, allocatable ::  sx_xu(:) 
      real*8, allocatable ::  sx_xv(:) 
      real*8, allocatable ::  sx_xw(:)
      real*8, allocatable ::  sx_yu(:) 
      real*8, allocatable ::  sx_yv(:) 
      real*8, allocatable ::  sx_yw(:)
      real*8, allocatable ::  sx_zu(:) 
      real*8, allocatable ::  sx_zv(:) 
      real*8, allocatable ::  sx_zw(:)      
      
      real*8, allocatable ::  sx_temp_xu(:,:,:) 
      real*8, allocatable ::  sx_temp_xv(:,:,:) 
      real*8, allocatable ::  sx_temp_xw(:,:,:)
      real*8, allocatable ::  sx_temp_yu(:,:,:) 
      real*8, allocatable ::  sx_temp_yv(:,:,:) 
      real*8, allocatable ::  sx_temp_yw(:,:,:)
      real*8, allocatable ::  sx_temp_zu(:,:,:) 
      real*8, allocatable ::  sx_temp_zv(:,:,:) 
      real*8, allocatable ::  sx_temp_zw(:,:,:)

      integer, allocatable ::  idum_str(:,:)  
      real*8, allocatable ::   dum_str(:,:)   
      integer, allocatable ::  idum_str1(:)   
      
      integer, allocatable :: ipermx(:,:)
      integer, allocatable :: ipermy(:,:)
      integer, allocatable :: ipermz(:,:)      
      real*8, allocatable ::  bp_flow1(:) 
      real*8, allocatable ::  bp_flow2(:)  
      
      real*8, allocatable ::  spm1f(:) 
      real*8, allocatable ::  spm2f(:)  
      real*8, allocatable ::  spm3f(:)  
      real*8, allocatable ::  spm4f(:)  
      real*8, allocatable ::  spm5f(:)
      real*8, allocatable ::  spm6f(:)  
      real*8, allocatable ::  spm7f(:)  
      real*8, allocatable ::  spm8f(:) 
      real*8, allocatable ::  spm9f(:) 
      real*8, allocatable ::  spm10f(:)
      real*8, allocatable ::  spm11f(:)
      real*8, allocatable ::  spm12f(:)
      real*8, allocatable ::  spm13f(:)
      real*8, allocatable ::  spm14f(:)
      real*8, allocatable ::  spm15f(:)
      real*8, allocatable ::  spm16f(:)
      real*8, allocatable ::  spm1f222(:) 
      real*8, allocatable ::  spm2f222(:)  
      real*8, allocatable ::  spm3f222(:)  
      real*8, allocatable ::  spm4f222(:)  
      real*8, allocatable ::  spm5f222(:)
      real*8, allocatable ::  spm6f222(:)  
      real*8, allocatable ::  spm7f222(:)  
      real*8, allocatable ::  spm8f222(:) 
      real*8, allocatable ::  spm9f222(:) 
      real*8, allocatable ::  spm10f222(:)

      real*8, allocatable ::  flitho(:,:) 
      
      integer, allocatable ::  iarea_str(:,:)  
      real*8, allocatable ::  area_str(:,:)  
      
      real*8 frac_vol

c     For Min model
      real*8 ipmd4_fx,ipmd4_br,ipmd4_bmx,ipmd4_alx,ipmd4_aly
      real*8 ipmd4_fdx,ipmd4_dmx,ipmd4_gmx,ipmd4_kc
      real*8 ipmd4_fy,ipmd4_btx,ipmd4_bty, ipmd4_fdy,ipmd4_gmy

c s kelkar      
      integer flag_strainout_first 
      integer flag_strain_output, ifile_strain,n_incr
      integer nentries, nentries_young
      real*8, allocatable ::  e_ini(:), poisson_ini(:)
      real*8, allocatable ::   dNuedt(:)
      real*8, allocatable ::  dEdt(:)
      real*8, allocatable ::  xtan_min(:)        
      real*8, allocatable ::  disp0(:,:)
      real*8, allocatable ::  k_strs91(:,:)
c      real*8, allocatable ::  e_temp91(:,:)
c gaz 042116
c      
      real*8, allocatable ::  e_temp91(:,:,:)
      real*8, allocatable ::  t_non_ref(:)
      integer, allocatable ::  i_tab_youngs(:)
      integer, allocatable ::  iy_tab(:)
      integer, allocatable ::  istr_non_model(:)
      integer  max_y_tab, n_young_table, nentries_young_max
      integer max_non_str, nentries_biot
      parameter (max_non_str = 100)
      parameter (max_y_tab = 100,nentries_young_max = 10000)
c      
c gaz 042116
c      
      real*8, allocatable ::  biot_temp91(:,:,:)
      real*8, allocatable ::  biot_t_non_ref(:)
      integer, allocatable ::  i_tab_biot(:)
      integer, allocatable ::  iy_tab_biot(:)
      integer, allocatable ::  istr_non_model_biot(:)
      integer  max_y_tab_b, n_biot_table, nentries_biot_max
      integer max_non_str_biot
      parameter (max_non_str_biot = 100)
      parameter (max_y_tab_b = 100,nentries_biot_max = 10000)      

      real*8 knx_stressperm, kny_stressperm,knz_stressperm

      real*8, allocatable ::  excess_shear(:)
      real*8, allocatable ::  shear_angle(:)
      real*8 pore_factor

      real*8, allocatable ::  density(:), internal_energy(:)

      integer, allocatable ::  itemp_perm22(:)
      integer, allocatable ::  itemp_perm23(:)

c s karra 17May 2012
      real*8, allocatable ::  plasticParam3(:)
      integer :: flag_pstrain_perm_coupling
c
c s kelkar for permmodel_22 zero out initial excess hsear stress
      real*8, allocatable ::  str_x0_perm(:)
      real*8, allocatable ::  str_y0_perm(:)
      real*8, allocatable ::  str_z0_perm(:)
      real*8, allocatable ::  str_xy0_perm(:)
      real*8, allocatable ::  str_xz0_perm(:)
      real*8, allocatable ::  str_yz0_perm(:)
      real*8, allocatable ::  str_pf0_perm(:)

c d dempsey, for permmodel_25, save fracture orientations from file
      real*8, allocatable ::  str25_dip(:)
      real*8, allocatable ::  str25_azi(:)
      real*8, allocatable ::  str25_xfrac(:,:)
      real*8, allocatable ::  str25_yfrac(:,:)
      real*8, allocatable ::  str25_zfrac(:,:)
      real*8, allocatable ::  perm_mult1(:)
      real*8, allocatable ::  perm_mult2(:)
      real*8, allocatable ::  perm_mult3(:)
      real*8 str25_density(1000)
      integer str25_N_obs
      
      
c gaz 110517    
      real*8 eigenvec(3,3),alambda(3), eigenvec_deg(3)


      end module comsi

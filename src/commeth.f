      module commeth
!     commeth 
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
!D1 Include file for array variables for hydrate variables.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2
!D2 Initial implementation: 06-Jun-02, Programmer: G. Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/commeth.f_a  $
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


      real*8, allocatable ::  denmethi(:)
      real*8, allocatable ::  denmethh(:)
      real*8, allocatable ::  denemethi(:)
      real*8, allocatable ::  denemethh(:)
      real*8, allocatable ::  pcpmeth(:)
      real*8, allocatable ::  pflowmeth(:) 
      real*8, allocatable ::  phimeth(:) 
      real*8, allocatable ::  phometh(:)
      real*8, allocatable ::  smeth(:) 
      real*8, allocatable ::  someth(:) 
      real*8, allocatable ::  tmeth(:) 
      real*8, allocatable ::  tometh(:) 
      real*8, allocatable ::  skmeth(:)
      real*8, allocatable ::  skmeth_tot(:)
      real*8, allocatable ::  skmethw(:)
      real*8, allocatable ::  qhmeth(:)
      real*8, allocatable ::  qhmethw(:)
      real*8, allocatable ::  eskmeth(:)
      real*8, allocatable ::  eflowmeth(:)
      real*8, allocatable ::  wellmeth(:)
      real*8, allocatable ::  rovfmeth(:)
      real*8, allocatable ::  envfmeth(:)
      real*8, allocatable ::  rolfmeth(:)
      real*8, allocatable ::  enlfmeth(:)
      real*8, allocatable ::  fracw(:)
      real*8, allocatable ::  fracwo(:)
      real*8, allocatable ::  qhflxmeth(:)

      real*8, allocatable ::  dpcp3(:)
      real*8, allocatable ::  dpcp4(:)
      real*8, allocatable ::  dpcpmeth2(:)
      real*8, allocatable ::  dtpameth(:)
      real*8, allocatable ::  dtpaemeth(:)

      integer, allocatable ::  kameth(:)
       
      real*8, allocatable ::  dmmf(:)
      real*8, allocatable ::  demf(:)
      real*8, allocatable ::  dilm(:)
      real*8, allocatable ::  divm(:)
      real*8, allocatable ::  deqm(:)
      real*8, allocatable ::  deqpm(:)

      real*8, allocatable ::  dmwf(:)
      real*8, allocatable ::  dewf(:)
      real*8, allocatable ::  dilw(:)
      real*8, allocatable ::  divw(:)
      real*8, allocatable ::  divg(:)
      real*8, allocatable ::  deqw(:)
      real*8, allocatable ::  deqpw(:)

      real*8, allocatable ::  denhydi(:)
      real*8, allocatable ::  denhydh(:)
      real*8, allocatable ::  denehydi(:)
      real*8, allocatable ::  denehydh(:)
      real*8, allocatable ::  frachyd(:)
      real*8, allocatable ::  frachydo(:)
      real*8, allocatable ::  fracgas(:)
      real*8, allocatable ::  fracgaso(:)
      real*8, allocatable ::  enlfhyd(:)
      real*8, allocatable ::  phihyd(:)
      real*8, allocatable ::  permhyd(:)
      real*8, allocatable ::  permhyd0(:)
      real*8, allocatable ::  dpermhyd4(:)
      real*8, allocatable ::  skhyd(:)
      real*8, allocatable ::  qhhyd(:)
      real*8, allocatable ::  skmhyd(:)
      real*8, allocatable ::  qhmhyd(:)
      real*8, allocatable ::  skwhyd(:)
      real*8, allocatable ::  qhwhyd(:)
      real*8, allocatable ::  dskhyd1(:)
      real*8, allocatable ::  dqhhyd1(:)
      real*8, allocatable ::  dskhyd2(:)
      real*8, allocatable ::  dqhhyd2(:)
      real*8, allocatable ::  dskhyd3(:)
      real*8, allocatable ::  dqhhyd3(:)
      real*8, allocatable ::  dskhyd4(:)
      real*8, allocatable ::  dqhhyd4(:)
      real*8, allocatable ::  dskmethw1(:)
      real*8, allocatable ::  dskmethw3(:)
      real*8, allocatable ::  dskmethw4(:)
      real*8, allocatable ::  dqm(:)
      real*8, allocatable ::  dqmt(:)
      real*8, allocatable ::  dqmh(:)
      real*8, allocatable ::  deqmh(:)
      real*8, allocatable ::  dqmw(:)
      real*8, allocatable ::  dqmhyd4(:)
      real*8, allocatable ::  dqd(:,:)
      real*8, allocatable ::  nhyd(:)
c      integer, allocatable ::  nhyd(:)
      integer, allocatable ::  ieoso(:)
      integer, allocatable ::  iceso(:)
      integer, allocatable ::  ihyd(:)
      integer, allocatable ::  ihydo(:)

      real*8, allocatable ::  dq3(:)
      real*8, allocatable ::  dq4(:)

      real*8, allocatable ::  rl_meth(:)
      real*8, allocatable ::  drl_meth3(:)
      real*8, allocatable ::  drl_meth4(:)
      real*8, allocatable ::  rl_h2o(:)
      real*8, allocatable ::  drl_h2o3(:)
      real*8, allocatable ::  drl_h2o4(:)
      real*8, allocatable ::  ahyd_coef(:)
c gaz 2005/August (related to general dissociation models)     
      real*8, allocatable ::  skhyd_temp(:)
      real*8, allocatable ::  dskhydpf(:)
      real*8, allocatable ::  dskhydtf(:)
      real*8, allocatable ::  dskhydmf(:)
      real*8, allocatable ::  dskhydwf(:) 
c hydrate dissociation (formation) model parameters  
      integer max_hyd_model
      parameter (max_hyd_model = 100)
      integer, allocatable ::  idisst(:)
      integer, allocatable ::  idissp(:)
      real*8, allocatable ::  diss1f(:)
      real*8, allocatable ::  diss2f(:)
      real*8, allocatable ::  diss3f(:)
      real*8, allocatable ::  diss4f(:)
      real*8, allocatable ::  diss5f(:)
      real*8, allocatable ::  diss6f(:)
      real*8, allocatable ::  diss7f(:)
      real*8, allocatable ::  diss8f(:)
      real*8, allocatable ::  diss9f(:)
      real*8, allocatable ::  form1f(:)
      real*8, allocatable ::  form2f(:)
      real*8, allocatable ::  form3f(:)
      real*8, allocatable ::  form4f(:)
      real*8, allocatable ::  form5f(:)
      real*8, allocatable ::  form6f(:)
      real*8, allocatable ::  form7f(:)
      real*8, allocatable ::  form8f(:) 
      real*8, allocatable ::  form9f(:)       

      real*8, allocatable ::  a44i(:)
      real*8, allocatable ::  fracw_tot(:)
      real*8, allocatable ::  fracg_tot(:)

c tenma 2005/March
*     real*8, allocatable :: qhhyd_tot(:)
      real*8, allocatable :: qhh2o_tot(:)
      real*8, allocatable :: qhmeth_tot(:)
*     real*8, allocatable :: skhyd_tot(:)
      real*8, allocatable :: skh2o_tot(:)
c tenma 2005/April
      real*8, allocatable :: frac_mudw(:)
      real*8, allocatable :: poro_sand(:)
      real*8, allocatable :: poro_mud(:)
      real*8, allocatable :: z_sand(:)
      real*8, allocatable :: z_mud(:)
      real*8, allocatable :: permhyd_z(:)
      real*8, allocatable :: ka0_xy(:)
      real*8, allocatable :: perm_sand(:)
      real*8, allocatable :: perm_mud(:)
      real*8, allocatable ::  ahyd_coef_1(:)
c tenma 2005/06 => sakamoto model
      real*8, allocatable ::  frachyd_1(:)
      integer, allocatable ::  nhyd_2(:)
      real*8, allocatable :: nhyd_3(:)
      real*8, allocatable :: frachyd_max(:)

      real*8, allocatable :: fracg2(:)
      real*8, allocatable :: fracw2(:)
c
      real*8  ammeth0,aemeth0,ammeth,aemeth
      real*8  qmeth_in,qemeth_in,balmeth,balemeth
      real*8  qmeth,qemeth,qmethts,qemethts
      real*8  qmethts_in,qemethts_in
      real*8  aihyd                    

      real*8  amhyd0,aehyd0,amhyd,aehyd
      real*8  qhyd_in,qehyd_in,balhyd,balehyd
      real*8  qhyd,qehyd,qhydts,qehydts
      real*8  qhydts_in,qehydts_in

      real*8  amh2o0,aeh2o0,amh2o,aeh2o
      real*8  qh2o_in,qeh2o_in,balh2o,baleh2o
      real*8  qh2o,qeh2o,qh2ots,qeh2ots
      real*8  qh2ots_in,qeh2ots_in

      real*8 hdiss01,hdiss11,hdiss02,hdiss12
      real*8 da,oe,pe,le,me,ne,e_act1,e_act2,afhyd
      real*8 dar,oer,per,ler,mer,ner,e_act1r,e_act2r,afhydr

c    porh is pass-through of variable ps(i)
c    volh is pass-through of variable sx1(i)
      real*8 porh,volh,strd_meth,strd_meth0
      real*8 ratw,ratgas,rathyd,fracv,fracl
c    ptc1-ptc4 are patameters for the hydrate pt line
      integer imeth_pt
      real*8 ptc1,ptc2,ptc3,ptc4
      real*8 source_hyd_temp,ds1_tmp,ds2_tmp,ds3_tmp,ds4_tmp

      integer idof_meth, ihydfl, ishyd, ic_hyd, ic_hyd_tot, numhyd
      integer icmw,icmm,icmh
      integer ihyd_diss_type, ihyd_grow_type

      end module commeth

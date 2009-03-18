      module comsi
!*******************************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global include file for array variables and pointers (FEHMN application).
!D1 Coupled Fluid Stress application
!D1
!*******************************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 06-June_2002  G. Zyvoloski           Initial implementation.
!D2
!*******************************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!*******************************************************************************
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
!*******************************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!*******************************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!*******************************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!*******************************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!*******************************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 N/A
!D9
!*******************************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!*******************************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS None
!PS
!*******************************************************************************
      integer icycs,ipsps,iinp,minks,inst,ipchng,nnx,ihms,nvfc
      integer nomass,noheat,noydis,localx,localy,ifrac,idof_stress
	integer ibp_stress, ibodyforce, ipermstr, initcalc, istrshis 

      real*8 daystr,fpor,fric,wo,fwght,fupwt,fdnwt,tol_stress
      real*8 tptch,pchmin,pchmax,tchmin,tchmax,coftol,bp_stress
	real*8 bpx,bpy,bpz,tol_stress1,abs_tol_stress,bp_update
	real*8 vol_tot_change
	real*8 stry_min,strx_min,e10_fac,str_mult
      integer, allocatable :: kr(:,:)
      integer, allocatable :: npbn(:)
      integer, allocatable :: nvfcl(:)
      integer, allocatable :: nskw_stress(:,:)

      real*8, allocatable ::  elastic_mod(:)
      real*8, allocatable ::  poisson(:)
      real*8, allocatable ::  e1(:)
      real*8, allocatable ::  e2(:)
      real*8, allocatable ::  e3(:)
      real*8, allocatable ::  bulk(:)
      real*8, allocatable ::  alp(:)
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
      real*8, allocatable ::  str_xy(:) 
      real*8, allocatable ::  str_xz(:)
      real*8, allocatable ::  str_yz(:)
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

c      common/fracturex/dla(n0),dlb(n0),dmpm(n0),dmtm(n0),
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

      end module comsi




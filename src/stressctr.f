      subroutine stressctr(iflg,ndummy) 
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To manage the (fluid-stress) calculations
!D1 for multi-component systems. 
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20 
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/stressctr.f_a  $
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY 
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comei
      use comfem
      use comfi
      use comgi
      use comii
      use comji
      use comki
      use comsi
      use comsplitts
      use comxi
      use davidi
      
      implicit none
      
      integer iflg,i,ndummy,md,j,k,isstr_temp, neqp1
      integer i1,i2,jj,kb,kc,iforce,nr1,nr2,ieosd
      integer il,ilev,mlev, node, inptorig, num_tab, i_non_strs 
      character*10 macro1, macro2
      real*8 dis_tol,aiter,aminkt
      real*8, allocatable :: stressboun(:)
      integer, allocatable :: kq_dum(:)
      real*8 wgt_tot,wgt_toti
      real*8 wgt_facx,wgt_totx
      real*8 wgt_facy,wgt_toty
      real*8 wgt_facz,wgt_totz
      real*8 xi,yi,zi,xkb,ykb,zkb,disx,disy,disz
      real*8 dui,dvi,dwi,ddu,ddv,ddw
      real*8 e1i,e2i,e3i,e1kb,e2kb,e3kb
      real*8 e1bar,e2bar,e3bar
      real*8 dudx,dudy,dudz
      real*8 dvdx,dvdy,dvdz
      real*8 dwdx,dwdy,dwdz
      real*8 str_x_avg,str_y_avg,str_z_avg
      real*8 str_xy_avg,str_xz_avg,str_yz_avg
      real*8 str_x_tot,str_y_tot,str_z_tot
      real*8 str_xy_tot,str_xz_tot,str_yz_tot
      integer idisx_max,idisy_max,idisz_max
      integer idisx_min,idisy_min,idisz_min
      real*8 disx_max,disy_max,disz_max
      real*8 disx_min,disy_min,disz_min
      real*8 area_facex, area_facey, area_facez, vol_change
      real*8 cordxi, cordyi, cordzi, delx, dely, delz
      real*8 dispxi,dispyi,dispzi,dispx,dispy,dispz
      real*8 str_ratio_xy, str_ratio_xz, str_ratio_yz
      real*8 biot,erat,efac,epi,dpd,shpi,stress_min_prin,lith_min
      real*8 ctherm,eti,shti,lith_damage_fac,str_eff,pterm
      real*8 xdV, ydV, zdV
      integer el
      
      real*8 denf, por, por_new, sx1d
      real*8 denfe, cp_t, ener_fl
      real*8 dporeeui,dporeeukb,ddenieeui
      real*8 dporeevi,dporeevkb,ddenieevi
      real*8 dporeewi,dporeewkb,ddenieewi
      real*8 dmeeui, dmeeukb, deneeui, deneeukb
      real*8 dmeevi, dmeevkb, deneevi, deneevkb
      real*8 dmeewi, dmeewkb, deneewi, deneewkb
      
      real*8, allocatable ::  modulus(:), nu(:)
      integer itmp

      integer jmia,mdkb
      integer iws, max_zone_number
      real*8 sjsix,sjsiy,sjsiz,alpkb,alphab,alpi
      real*8 bulkkb,bulki,bulkb,tdumt,pdumt,elev,sq3
      real*8 termt, termp, bulk_mod
      real*8 grav_save, sdepth, gdepth
      character*24 hist_str
      character*4 stype
      character*6 ptype
      character*9 sdum1
      integer nodestress, icount, ilithod
      logical null1
      
      save isstr_temp
      
      real*8 area_ix, area_iy, area_iz, area_tol
      real*8 dvol_strainu_i, dvol_strainv_i, dvol_strainw_i
      real*8 strx,strxgrad 
      integer open_file, ifile
      
c     s kelkar 12/6/09 axisymmetric anisotropy
      real*8 young_p , young_t, pois_p, pois_t, pois_sq
      real*8 fac1,fac2,fac3, ezzi,ezzkb,ezzbar,efacxy,efacz 
      real*8 shtixy,shtiz,shpixy,shpiz, e4i,e4kb,e4bar
      real*8 efac_ks, betat, efac_betat
      real*8 bulk_tol

      character*100 perm_stress_file, young_temp_file
      integer i91,j91
c     s kelkar nov 5 2010
      real*8  pi
c      real*8 eigenvec(3,3),alambda(3), eigenvec_deg(3)
      real*8 friction, strength, pp_fac
      integer iispmd
      integer ishear

      character*150 input_msg
      integer imsg(16)
      character*32 cmsg(16)
      integer msg(16)
      integer nwds
      real*8 xmsg(16)

      real*8 temp_junk
      real*8 stress_factor_initial
	  
c     d dempsey Feb 2013
      character*80 ltest2
      logical null_new
      integer line_number
c..................................................................
      
      parameter (pi=3.1415926535)
      parameter (area_tol = 1.d-18)
      parameter(dis_tol=1.d-12)
      parameter(bulk_tol=1.d-12)
      parameter(max_zone_number = 1000)
      if(istrs.eq.0) return  
      
      if(iflg.eq.0) then
c     
c     first allocate variables      
c     
         allocate (kr(n0,3))   
c     allocate (npbn(n0))
c     allocate (nvfcl(n0))
         allocate (elastic_mod(n0))
         allocate (poisson(n0))
         allocate (e1(n0))
         allocate (e2(n0))
         allocate (e3(n0))
c     s kelkar 4/20/2010   gaz  042916 allocated later
c         allocate (e_ini(1:n0))
c         allocate (dEdt(1:n0))
c         allocate (dNuedt(1:n0))
c         allocate (poisson_ini(1:n0))
c     d dempsey 3/11/2014
         allocate (bodyforce_x(n0))
         allocate (bodyforce_y(n0))
         allocate (bodyforce_z(n0))
         
c     s kelkar 12/6/09 axisymmetric anisotropy
c     for now, arrays are allocated while reading input,lines ~684-700
c......................................
         allocate (bulk(n0))
         allocate (alp(n0))
         allocate (bulk0(n0))
         allocate (alp0(n0))
         allocate (du(n0)) 
         allocate (dv(n0)) 
         allocate (dw(n0))

c     Sai : allocate delta_(u/v/w) - incremental displacements needed for plasticity
         allocate(delta_u(n0))
         allocate(delta_v(n0))
         allocate(delta_w(n0))
         delta_u = 0.0d0
         delta_v = 0.0d0
         delta_w = 0.0d0

         allocate (duo(n0)) 
         allocate (dvo(n0)) 
         allocate (dwo(n0))
         allocate (du_ini(n0)) 
         allocate (dv_ini(n0)) 
         allocate (dw_ini(n0))
         allocate (du_tot(n0)) 
         allocate (dv_tot(n0)) 
         allocate (dw_tot(n0))
         
         allocate (str_x(n0)) 
         allocate (str_y(n0))
         allocate (str_z(n0))
         allocate (str_xy(n0)) 
         allocate (str_xz(n0))
         allocate (str_yz(n0))
         
         allocate(strain_xx(n0))
         allocate(strain_yy(n0))
         allocate(strain_zz(n0))
         allocate(strain_xy(n0))
         allocate(strain_yz(n0))
         allocate(strain_zx(n0))

         allocate (disp(n0,3))
         allocate (disp0(n0,3))
         allocate (forc(n0,3))
         
         allocate (vol_strain(n0))
         allocate (vol_temp(n0))
         allocate (vol_strain0(n0))
         allocate (idum_str(n0,3))
         allocate (dum_str(max_zone_number,3))
         allocate (iarea_str(n0,3))
         allocate (area_str(n0,3))
c     arrays for stress derivatives in mass and energy equations
         allocate (its21(1,4))
         allocate (its22(100,4))
         allocate (its31(1,4))
         allocate (its32(100,4))
         allocate (its41(1,4))
         allocate (its42(100,4))
         allocate (itstress(200))
         allocate (ts21(1,4))
         allocate (ts22(100,4))
         allocate (ts31(1,4))
         allocate (ts32(100,4))
         allocate (ts41(1,4))
         allocate (ts42(100,4))
         allocate (ts51(100))
         allocate (ts52(100))
         allocate (ts53(100))
c
   
c     zero variables
c     
         incremental_shear_permmodel = 0
         flag_element_perm = 0
         ifem = 0

         flag_strainout_first = 0
         flag_permmodel = 0
c     flag_plastic_message = 0
         ibodyforce = 0
         flag_principal=0
         flag_excess_shear = 0
         flag_strain_output = 0
c     set to constant permeability
         ipermstr =  1	
         initcalc = 0
         istrshis = 0
         istresspor = 0
         itotals_s = 0
         itert_s = 0
         itotal_s = 0
         idisp_rel = 0
         ilithgrad = 0
         kr= 0   
c     npbn = 0
c     nvfcl = 0
         tol_stress = 0.0d0
         
         elastic_mod = 0.0d0
         poisson = 0.0d0
         e1 = 0.0d0
         e2 = 0.0d0
         e3 = 0.0d0
c     axisymmetric anisotropy
         stress_anisotropy_in = .false.
         stress_anisotropy_use = .false.
		 
         bodyforce_x = 0.0d0
         bodyforce_y = 0.0d0
         bodyforce_z = 0.0d0
c......................................
         bulk = 0.0d0
         alp = 0.0d0
         alp0 = 0.0d0
         bulk0 = 0.0d0         
         du = 0.0d0 
         dv = 0.0d0 
         dw = 0.0d0
         duo = 0.0d0 
         dvo = 0.0d0 
         dwo = 0.0d0
         du_ini = 0.0d0 
         dv_ini = 0.0d0 
         dw_ini = 0.0d0
         du_tot = 0.0d0 
         dv_tot = 0.0d0 
         dw_tot = 0.0d0
         
         str_x = 0.0d0 
         str_y = 0.0d0
         str_z = 0.0d0
         str_xy = 0.0d0 
         str_xz = 0.0d0
         str_yz = 0.0d0
         
         disp = 0.0d0
         forc = 0.0d0
         
c     dluf = 0.0d0 
c     dlvf = 0.0d0
c     shp = 0.0d0 
c     sht = 0.0d0
c     uxtp = 0.0d0 
c     uytp = 0.0d0
c     uztp = 0.0d0
c     dmadv = 0.0d0
c     deadv = 0.0d0
c     drluf = 0.0d0
c     drlvf = 0.0d0
         
c     drlzv = 0.0d0
c     dflv = 0.0d0
c     dvfdv = 0.0d0
         vol_strain = 0.0d0
         vol_strain0 = 0.0d0
         
         idum_str = 0
         dum_str = 0.0d0
         iarea_str = 0
         area_str = 0.0
         
         idof_stress = 0
         abs_tol_stress = 1.d-10 
c     Sai: initially set nonlinear flag to zero
         isNonlinear = 0
         isbiotNonLin = 0
         
c     
c     temporary unit number for stress contour
c     
	 hist_str(1:14) = 'history_stress'
         
c     read in input when stress is present
         read (inpt  ,   *) istrs, ihms
         if(ihms.eq.-4) then
            backspace inpt
            read (inpt  ,   *) istrs, ihms, daystress
         endif         
c     istrs = 0 - skip stress solution
c     istrs = 1 - plain strain and 3-D (hookean) solution
c     istrs = 2 - plain stess (hookean) solution (must be 2-D)
c     ihms - identifies the coupling and  stress solution 
c     ihms gt 0 fully coupled
c     ihms = 15 : derivatives of pore volume only wrt displacements 
c                 in heat-mass equations, but not those of permeability
c     (ihms = -15  : allow sequentially coupled porosity changes, 
c                    ihms is reset to -3)
c     ihms = 16 : derivatives of permeability only in heat and 
c                 mass equations, but not those of pore volume
c     ihms = 17 : derivatives of both permeability and pore volume  
c                 in heat-mass equations
c     ihms = -1 only at the end of the simulation
c     ihms = -2 beginning and end of simulation
c     ihms = -3 end of each time step
c     ihms = -4 end of each time segment(or less) requires addtional 
c     input daystress
c     this parameter will not change the flow time step size
c     with the exception of ihms = -3 tini and pini
c     will not be updated after each iteration
c     
         if(ihms.eq.15.or.ihms.eq.-15.or.ihms.eq.16.or.ihms.eq.17) then
            pore_factor = 1
            backspace inpt
            read (inpt, '(a)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
            if(nwds.gt.2) then
               pore_factor = xmsg(3)
            endif
            endif
         if(ihms.eq.15.or.ihms.eq.17) then
            if(.not.allocated(density)) allocate(density(n0))
            if(.not.allocated(internal_energy))
     &           allocate(internal_energy(n0))
         endif

            if(ihms.eq.-15) ihms=-3

            if(ihms.eq.16.or.ihms.eq.17) flag_element_perm =  1

c     
         if(icnl.ne.0.and.ihms.lt.0) then
c     not fully coupled  2D solution    
            idof_stress = 2
            istrs_coupl = ihms
         else if(icnl.eq.0.and.ihms.lt.0) then
c     not fully coupled          
            idof_stress = 3
            istrs_coupl = ihms
         else if(icnl.ne.0.and.ihms.gt.0) then 
c     fully coupled but mat need to be enlarged in startup      
            idof_stress = 4  
c     istrs_coupl = ihms
         else if(icnl.eq.0.and.ihms.gt.0) then  
c     fully coupled but mat need to be enlarged in startup     
            idof_stress = 5
            istrs_coupl = ihms
         endif
         if(istrs.eq.0) then
            idof_stress=0
            istrs_coupl=0
         endif
         if (idof_stress.ne.2.and.istrs.eq.2) then
            write (ierr, 19)
            if (iout .ne. 0) write (iout, 19)
            stop
         endif
         
         
         
 19      format(/, 'plain stress is not 2-D: stopping', /)
 20      format(3x,'**** stress conditions: 1D (uncoupled) ****') 
 21      format(3x,'**** stress conditions: 2D (uncoupled) ****') 
 22      format(3x,'**** stress conditions: 3D (uncoupled) ****') 
 23      format(3x,'**** stress conditions: 1D (coupled) ****')
 24      format(3x,'**** stress conditions: 2D (coupled) ****') 
 25      format(3x,'**** stress conditions: 3D (coupled) ****')
 26      format(3x,'**** not used ****')
         
         macro = "strs"
c     
c     read in "sub macros" for stress
c     
         
         
 100     continue
         read (inpt, '(a80)') wdd1
c     Changed end check for compatibility with macro "off" option
         if (wdd1(1:9) .eq. 'stressend' .or. wdd1(1:8) .eq. 'end strs'
     &        .or. wdd1(1:7) .eq. 'endstrs') go to 200 
         if (wdd1(1:1) .eq. '#') go to 100 
         read (wdd1, '(a10)') macro1
         if (iout .ne. 0) write(iout, 50) macro1
         if (iptty .gt. 0) write(iptty, 50) macro1 
 50      format(3x, '**** stress sub macro : ', a10,' **** ' ) 
         if(macro1.eq.'text      ') then
c     
c     read text here
c     
c     
c     set zone call within stress module
c     
         else if(macro1(1:4).eq.'zone') then
            inptorig = inpt
            macro = macro1(1:4)
            cnum_stress = cnum_stress + 1
            call start_macro(inptorig, inpt, macro)
            call zone(cnum_stress, inpt)
            call done_macro(inpt)
            inpt = inptorig
         else if(macro1(1:4).eq.'zonn') then
            inptorig = inpt
            macro = macro1(1:4)
            cnum_stress = cnum_stress + 1
            call start_macro(inptorig, inpt, macro)
            call zone(cnum_stress, inpt)
            call done_macro(inpt)
            inpt = inptorig
            
         else if(macro1.eq.'bodyforce ') then
c     
c     enable body force (weight of rock)
c     
            backspace inpt
            read (inpt, '(a80)') wdd1
            call parse_string(wdd1, imsg, msg, xmsg, cmsg, nwds)
            ibodyforce = 3
            if(nwds.gt.1) then
              if(cmsg(2).eq.'force') then
                ibodyforce = 4
              else if(cmsg(2).eq.'acceleration') then
                ibodyforce = 5
              endif
				
            igroup = 1
            narrays = 3
            itype(1) = 8
            itype(2) = 8
            itype(3) = 8
            default(1) = 0.
            default(2) = 0.
            default(3) = 0.
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           r8_1 = bodyforce_x(1:n0),r8_2 = bodyforce_y(1:n0),
     &           r8_3 = bodyforce_z(1:n0))
				
            endif

            
         else if(macro1.eq.'reldisp   ') then
c     
c     use relative displacements (for volumestrain, permmodels, contour)
c     
            idisp_rel = 1
            
         else if(macro1.eq.'stresspor ') then
c     
c     change porosity for explicit update
c     
            istresspor = 1        
            
         else if(macro1.eq.'initcalc  ') then
c     
c     calculate initial stress
c     
            initcalc = 1
c     
         else if(macro1.eq.'permmodel ') then
c     
c     enable permeability model
c     ipermstr = 1 is the default
c     
            flag_permmodel = 1
            if(.not.allocated(ispm)) then
               allocate(ispm(n0))
c     set default to model 1          
               ispm = 1
            endif
            
            i = 0
            j = 0
            ex = .false.
            do
               read(inpt,'(a80)') wdd1
               if(null1(wdd1)) exit
               backspace inpt
               read(inpt,*) ispmd
               backspace inpt
               i = i+1           
               if(ispmd .eq.1) then
                  permfile=open_file('permifile.permi', 'unknown')
                  write(permfile,*) 
     &                 'node  time  check  Kyy  Sxx  Syy  Szz  P'
c     default and no input            
                  read(inpt,*) ispmt(i)   
c     model 1 default            
c     model 2 (linear perm variation with stress)
c     model 3 general diaplacement based formulation (fully coupled)
c     model 5 general diaplacement (modified) based form (explicit)
c     Model 3 and 5 can be explicitly coupled or fully coupled
c     model 4 (cubic perm variation with stress)        
c     model 91 table lookup, read table from a file
c     model 6(linear perm variation  with effective stress failure)
c     
               else if (ispmd .eq. 2. or.ispmd.eq.222) then
c     simple explicit tensile stress model
c     spm3f,spm6f,spm9f are ignored for a 2D problem
c     spm1f:strx_min, min tensile stress (x direction) for damage to occur
c     spm2f:stry_min, min tensile stress (y direction) for damage to occur
c     spm3f:stry_min, min tensile stress (z direction) for damage to occur
c     spm4f:e10_facx, damage factor (maximum x) for elastic modulus
c     spm5f:e10_facy, damage factor (maximum y) for elastic modulus
c     spm6f:e10_facz, damage factor (maximum z) for elastic modulus
c     spm7f:str_multx, maximum change in permeability (x direction) allowed 
c     spm8f:str_multy, maximum change in permeability (y direction) allowed 
c     spm9f:str_multz, maximum change in permeability (z direction) allowed 
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i)
               else if (ispmd .eq. 21) then
                  if (.not.allocated(xtan_min)) allocate (xtan_min(n0))
c     mohr-coulomb failure criteria on a user-specified plane
c     spm1f,spm2f,spm3f=direction cosins of the normal to fault plane
c     spm4f:friction coefficient of shear in the fault plane
c     spm5f:shear strength of the fault plane
c     spm6f:factor in effective stress=sigma-(pp_fac*pore pressure)
c     spm7f:range of excess shear stress over which damage is ramped
c     spm8f:maximum multiplier for young's modulus  in x-prime direction
c     spm9f:maximum multiplier for young's modulus  in y-prime direction
c     spm10f:maximum multiplier for young's modulus  in z-prime direction
c     spm11f: maximum multiplier  permeability x-prime direction
c     spm12f: maximum multiplier  permeability y-prime direction
c     spm13f: maximum multiplier  permeability z-prime direction

                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i),spm10f(i),spm11f(i),spm12f(i),spm13f(i)
               else if (ispmd .eq. 22) then
                  if(.not.allocated(itemp_perm22)) 
     &                 allocate (itemp_perm22(neq))
                  itemp_perm22 =0
                  spm11f(i)=1.0
                  spm12f(i)=0.
                  spm13f(i)=0.
                  open(unit=94,file='failed_nodes_perm22.dat')
c     mohr-coulomb failure criteria on the plane that miximizes
c     abs(shear)-friction*normal stress
c     spm1f:friction coefficient of shear in the fault plane
c     spm2f:shear strength of the fault plane
c     spm3f:factor in effective stress=sigma-(pp_fac*pore pressure)
c     spm4f:range of excess shear stress over which damage is ramped
c     spm5f:maximum multiplier for young's modulus  in x-prime direction
c     spm6f:maximum multiplier for young's modulus  in y-prime direction
c     spm7f:maximum multiplier for young's modulus  in z-prime direction
c     spm8f: maximum multiplier  permeability x-prime direction
c     spm9f: maximum multiplier  permeability y-prime direction
c     spm10f: maximum multiplier  permeability z-prime direction
c     spm11f: optional, maximum multiplier  porosity
c     spm12f: optional, permfactor is made a funciton of inceramental
c     spm13f: optional, excess shewar has to exceed this value before
c            permfactor becomes different from 1.
c     here z-prime is along tehnormal to the plane of failure, and
c     y-prime is along the median principal stress
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i),spm10f(i)
c..... s kelkar 18 sep 2012, allow optional porosity damage 
                  backspace inpt
                  read (inpt, '(a)') input_msg
                  call parse_string(input_msg, imsg, msg, xmsg, cmsg
     &                 , nwds)
                  if(nwds.gt.11) then
                     spm11f(i) = xmsg(12)
                  endif
c..... s kelkar 27 Nov 2012, allow permfactor upon failure to depend
c        on increamental shear stress with a user defined factor
c        must have invoked 'initcal' prior to input for this macro
                  backspace inpt
                  read (inpt, '(a)') input_msg
                  call parse_string(input_msg, imsg, msg, xmsg, cmsg
     &                 , nwds)
                  if(nwds.gt.12) then
                     if(.not. residual_stress.and.initcalc.eq.1) then
                        spm12f(i) = xmsg(13)
                        if(spm12f(i).gt.1.0e-18) then
                  if(.not.allocated(str_x0_perm)) then
                     allocate (str_x0_perm(n0))
                     allocate (str_y0_perm(n0))
                     allocate (str_z0_perm(n0))
                     allocate (str_xy0_perm(n0))
                     allocate (str_xz0_perm(n0))
                     allocate (str_yz0_perm(n0))
                     allocate (str_pf0_perm(n0))
                  endif
                  incremental_shear_permmodel = 1
                        endif
                     else
                  write(iptty,*)'must invoke initcalc prior to'
                  write(iptty,*)'this option with stressperm=22.'
                  write(iptty,*)'At present, can not use this with' 
                  write(iptty,*)' residual_stress option'
                  write(iptty,*)'STOP in stressctr'
                  write(ierr,*)'must invoke initcalc prior to'
                  write(ierr,*)'this option with stressperm=22.'
                  write(ierr,*)'At present, can not use this with'
                  write(ierr,*)'residual_stress option'
                  write(ierr,*)'STOP in stressctr'
                  stop
                     endif
                  endif
c..... s kelkar 03 Dec 2012, 
c        on increamental shear stress with a user defined factor
c        must have invoked 'initcal' prior to input for this macro
                  backspace inpt
                  read (inpt, '(a)') input_msg
                  call parse_string(input_msg, imsg, msg, xmsg, cmsg
     &                 , nwds)
                  if(nwds.gt.13) then
                        spm13f(i) = xmsg(14)
                     endif
c...............................
               else if (ispmd .eq. 24) then
c                  if(.not.allocated(itemp_perm22)) 
c     &                 allocate (itemp_perm22(neq))
c                  itemp_perm22 =0
                  incremental_shear_permmodel = 1
c                  spm11f(i)=1.0
c                  spm12f(i)=0.
c                  spm13f(i)=0.
c                  spm14f(i)=0.
                  if(.not.allocated(str_x0_perm)) then
                     allocate (str_x0_perm(n0))
                     allocate (str_y0_perm(n0))
                     allocate (str_z0_perm(n0))
                     allocate (str_xy0_perm(n0))
                     allocate (str_xz0_perm(n0))
                     allocate (str_yz0_perm(n0))
                     allocate (str_pf0_perm(n0))
                  endif
c                  open(unit=94,file='failed_nodes_perm22.dat')
c     David Dempsey, Jan 2013, mohr-coulomb with fracture distribution
c     spm1f:shear fracture toughness 
c     spm2f:static friction coefficient
c     spm3f:dynamic friction coefficient
c     spm4f: number of fractures per control volume
c     spm5f: shear displacement at which perm enhancement begins 
c     spm6f: shear displacement interval to complete perm enhancement
c     spm7f: total perm enhancement (in log(perm))
c     spm8f: fracture cohesion
c     spm9f: fracture density in control volume 
c     here z-prime is along the normal to the plane of failure, and
c     y-prime is along the median principal stress
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i)
c...............................
c     David Dempsey, Feb 2013, mohr-coulomb with fracture distribution
               else if (ispmd .eq. 25) then
                  incremental_shear_permmodel = 1
                  if(.not.allocated(str_x0_perm)) then
                     allocate (str_x0_perm(n0))
                     allocate (str_y0_perm(n0))
                     allocate (str_z0_perm(n0))
                     allocate (str_xy0_perm(n0))
                     allocate (str_xz0_perm(n0))
                     allocate (str_yz0_perm(n0))
                     allocate (str_pf0_perm(n0))
                     allocate (str25_xfrac(n0,1000))
                     allocate (str25_yfrac(n0,1000))
                     allocate (str25_zfrac(n0,1000))
c initial values of 2. indicate these have yet to be assigned
                     do line_number=1,n0
                       str25_xfrac(line_number,1)=2.
                     enddo
                  endif
c read in fracture orientation data from 'fracture_orientations.dat'
                  allocate (str25_dip(10000))
                  allocate (str25_azi(10000))
                  open(unit=92,file='fracture_orientations.dat')
                  str25_N_obs=0
                  str25_density = 0.
                  do line_number=1,10000
                    read(92,'(a80)',end=7012) ltest2
                    if(.not.null_new(ltest2)) then
                      read(ltest2,*,end=7012,err=7012) 
     &                  str25_dip(line_number),str25_azi(line_number)
                      str25_N_obs = str25_N_obs + 1
                    else
                      goto 7012
                    endif
                  enddo					
 7012             continue
c     spm1f: shear fracture toughness
c     spm2f: static friction coefficient
c     spm3f: dynamic friction coefficient
c     spm4f: number of fractures per control volume
c     spm5f: shear displacement at which perm enhancement begins
c     spm6f: shear displacement interval to complete perm enhancement
c     spm7f: total perm enhancement (in log(perm))
c     spm8f: fracture cohesion

c     here z-prime is along the normal to the plane of failure, and
c     y-prime is along the median principal stress
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i)
                  
c................................................
               else if (ispmd .eq. 23) then
                  if(.not.allocated(itemp_perm23)) 
     &                 allocate (itemp_perm23(neq))
                  itemp_perm23 =0
                  open(unit=94,file='failed_nodes_perm23.dat')
c using the hydraulic aperture along the direction of shear movement
c a provided by Stuart Walsh from LANL, 6 Nov 2012.
c direction of shear slip assumed to be vertical. mohr-coulomb failure 
c criteria on the vertical plane 
c     abs(shear)-friction*normal stress
c     spm1f:friction coefficient of shear in the fault plane
c     spm2f:shear strength of the fault plane
c     spm3f:factor in effective stress=sigma-(pp_fac*pore pressure)
c     spm4f:(MPa) empirical modulus to conver post failure excess shear 
c           stress to shear dispalcement for permeability model
c     spm5f:(meters) maximum shear displacement allowed fo rpermeability
c     spm6f:aperture (meters) corrosponding to max shear displacement
c     spm7f:(metrers) lenght scale cubic aperture law
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i), spm7f(i)

c...................................
               else if (ispmd .eq. 222) then
c     combined model  explicit tensile and shear (2 and 22)
c     spm3f,spm6f,spm9f are ignored for a 2D problem
c     spm1f222:strx_min, min tensile stress (x direction) for damage to occur
c     spm2f222:stry_min, min tensile stress (y direction) for damage to occur
c     spm3f222:stry_min, min tensile stress (z direction) for damage to occur
c     spm4f222:e10_facx, damage factor (maximum x) for elastic modulus
c     spm5f222:e10_facy, damage factor (maximum y) for elastic modulus
c     spm6f222:e10_facz, damage factor (maximum z) for elastic modulus
c     spm7f222:str_multx, maximum change in permeability (x direction) allowed 
c     spm8f222:str_multy, maximum change in permeability (y direction) allowed 
c     spm9f222:str_multz, maximum change in permeability (z direction) allowed 
c     mohr-coulomb failure criteria on the plane that miximizes
c     abs(shear)-friction*normal stress
c     spm1f:friction coefficient of shear in the fault plane
c     spm2f:shear strength of the fault plane
c     spm3f:factor in effective stress=sigma-(pp_fac*pore pressure)
c     spm4f:range of excess shear stress over which damage is ramped
c     spm5f:maximum multiplier for young's modulus  in x-prime direction
c     spm6f:maximum multiplier for young's modulus  in y-prime direction
c     spm7f:maximum multiplier for young's modulus  in z-prime direction
c     spm8f: maximum multiplier  permeability x-prime direction
c     spm9f: maximum multiplier  permeability y-prime direction
c     spm10f: maximum multiplier  permeability z-prime direction
c     here z-prime is along tehnormal to the plane of failure, and
c     y-prime is along the median principal stress
                  if(.not.allocated(itemp_perm22)) 
     &                 allocate (itemp_perm22(neq))
                  itemp_perm22 =0
                  open(unit=94,file='failed_nodes_perm22.dat')
                  read(inpt,*)ispmt(i),spm1f222(i),spm2f222(i),
     &                 spm3f222(i),spm4f222(i),spm5f222(i),spm6f222(i),
     &                 spm7f222(i),spm8f222(i),spm9f222(i)
                  read(inpt,*) spm1f(i),spm2f(i),spm3f(i),spm4f(i),
     &                 spm5f(i),spm6f(i),spm7f(i),spm8f(i),spm9f(i),
     &                 spm10f(i)

               else if (ispmd .eq. 3) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i)
               else if (ispmd .eq. 4) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i),spm10f(i),spm11f(i),spm12f(i),spm13f(i)
               else if (ispmd .eq. 5) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i)
               else if (ispmd .eq. 6) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i),spm10f(i),spm11f(i)
               else if(ispmd .eq. 7) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i)
               else if (ispmd .eq. 8) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i),
     &                 spm4f(i),spm5f(i),spm6f(i),spm7f(i),spm8f(i),
     &                 spm9f(i) 
               else if (ispmd .eq.11) then
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i),spm3f(i)     
               else if (ispmd .eq. 100) then
c permmodel for VonMises and Drucker-Prager palsticity. Opens in tension
c  permeability is a ramp as a function of accumulated plastic strain
c  ramp range is spm4f.
c  at present it is isotropic, spm2f and spm3f are not used. 
                  read(inpt,*) ispmt(i),spm1f(i),spm2f(i),spm3f(i)
     &                 ,spm4f(i)
               else if (ispmd. eq. 31) then        
c     model for changing properties when ice is present
                  read(inpt,*)ispmt(i),spm1f(i),spm2f(i)
                  if(.not.allocated(dzrg_new)) allocate (dzrg_new(n0))
               else if (ispmd. eq. 91) then
                  if(.not.allocated(str_x0)) allocate (str_x0(n0)) 
                  if(.not.allocated(str_y0)) allocate (str_y0(n0))
                  if(.not.allocated(str_z0)) allocate (str_z0(n0))
                  read(inpt,*)ispmt(i)
                  read(inpt,'(a100)') perm_stress_file
                  ifile = open_file(perm_stress_file, 'old')
                  read(ifile,*)nentries 
                  allocate (k_strs91(nentries,4))
                  do i91=1,nentries
                     read(ifile,*)(k_strs91(i91,j91),j91=1,4)
                  enddo
                  close (ifile)
               endif
               
               if(ispmt(i).eq.1) ipermstr1 = ispmt(i)
               if(ispmt(i).eq.2) ipermstr2 = ispmt(i)
               if(ispmt(i).eq.3) ipermstr3 = ispmt(i)
               if(ispmt(i).eq.4) ipermstr4 = ispmt(i)
               if(ispmt(i).eq.5) ipermstr5 = ispmt(i)
               if(ispmt(i).eq.6) ipermstr6 = ispmt(i)
               if(ispmt(i).eq.7) ipermstr7 = ispmt(i)
               if(ispmt(i).eq.8) ipermstr8 = ispmt(i)
               if(ispmt(i).eq.11) ipermstr11 = ispmt(i)
               if(ispmt(i).eq.21) ipermstr21 = ispmt(i)
               if(ispmt(i).eq.22) ipermstr22 = ispmt(i)
               if(ispmt(i).eq.23) ipermstr23 = ispmt(i)
               if(ispmt(i).eq.24) ipermstr24 = ispmt(i)
               if(ispmt(i).eq.25) ipermstr25 = ispmt(i)
               if(ispmt(i).eq.31) ipermstr31 = ispmt(i)            
               if(ispmt(i).eq.222) ipermstr222 = ispmt(i)
               if(ispmt(i).eq.91) ipermstr91 = ispmt(i)
               if(ispmt(i).eq.100) ipermstr100 = ispmt(i)
            end do
c     set default         
            ispm = 1
            narrays = 1
            itype(1) = 4
            default(1) = 1
            igroup = 1
            call initdata2( inpt, ischk, n0, narrays,
     $           itype, default, macroread(8), macro, igroup, ireturn,
     $           i4_1=ispm(1:n0))	        
            
            
c     ****** end of input loop
         else if(macro1.eq.'initial  ') then
c     
c     read in initial stress state
c     
            
            igroup = 1
            narrays = 3
            itype(1) = 8
            itype(2) = 8
            itype(3) = 8
            default(1) = 0.
            default(2) = 0.
            default(3) = 0.
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           r8_1=duo(1:n0),r8_2=dvo(1:n0),r8_3 = dwo(1:n0))
            
            
         else if(macro1.eq.'ipini     ') then 
c     
c     use read-in initial pressure and temperature differences 
c     if they exist (from restart file)
c     
            ipini = 1
c.............................................................
         else if(macro1(1:3).eq.'fem') then
            ifem = 1
            backspace inpt
            read (inpt, '(a)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
            if(nwds.gt.1) then
               flag_element_perm = imsg(2)
            endif

c     Copy information about elements from nelm to elnode
            allocate(elnode(nei, ns))
            do i = 1,nei
               do j = 1,ns
                  elnode(i, j) = nelm((i-1)*ns + j)
               enddo
            enddo

c     gauss-point information
            if((icnl.eq.0) .and. (ns.eq.8)) then
               numgausspoints = 8
               allocate(gpcord(numgausspoints, 3))
               allocate(gpweight(numgausspoints))
               allocate(detJ(nei, numgausspoints))
               allocate(Psi(nei, numgausspoints, ns))
               allocate(dPsidX(nei, numgausspoints, ns))
               allocate(dPsidY(nei, numgausspoints, ns))
               allocate(dPsidZ(nei, numgausspoints, ns))
               allocate(fem_stress(nei, numgausspoints, 6))
               allocate(fem_strain(nei, numgausspoints, 6))
               allocate(conv_strain(nei, numgausspoints, 6))
               allocate(iPsi(ns, ns))

               call compute_ipsi()

               fem_stress = 0.0d0
               fem_strain = 0.0d0
               conv_strain = 0.0d0

               sq3 = 1.0d0/sqrt(3.0d0)
               gpcord(1,1) = -sq3
               gpcord(1,2) = -sq3
               gpcord(1,3) = -sq3
               gpweight(1) = 1.0d0

               gpcord(2,1) =  sq3
               gpcord(2,2) = -sq3
               gpcord(2,3) = -sq3
               gpweight(2) = 1.0d0

               gpcord(3,1) =  sq3
               gpcord(3,2) =  sq3
               gpcord(3,3) = -sq3
               gpweight(3) = 1.0d0

               gpcord(4,1) = -sq3
               gpcord(4,2) =  sq3
               gpcord(4,3) = -sq3
               gpweight(4) = 1.0d0

               gpcord(5,1) = -sq3
               gpcord(5,2) = -sq3
               gpcord(5,3) =  sq3
               gpweight(5) = 1.0d0

               gpcord(6,1) =  sq3
               gpcord(6,2) = -sq3
               gpcord(6,3) =  sq3
               gpweight(6) = 1.0d0

               gpcord(7,1) =  sq3
               gpcord(7,2) =  sq3
               gpcord(7,3) =  sq3
               gpweight(7) = 1.0d0

               gpcord(8,1) = -sq3
               gpcord(8,2) =  sq3
               gpcord(8,3) =  sq3
               gpweight(8) = 1.0d0

               call fem_shapefunctions_3r()

c 2D Fem in X-Y plane
            else if((icnl.eq.1) .and. (ns.eq.4)) then
               numgausspoints = 4
               allocate(gpcord(numgausspoints, 2))
               allocate(gpweight(numgausspoints))
               allocate(detJ(nei, numgausspoints))
               allocate(Psi(nei, numgausspoints, ns))
               allocate(dPsidX(nei, numgausspoints, ns))
               allocate(dPsidY(nei, numgausspoints, ns))
               allocate(fem_stress(nei, numgausspoints, 3))
               allocate(fem_strain(nei, numgausspoints, 3))

               fem_stress = 0.0d0
               fem_strain = 0.0d0

               sq3 = 1.0d0/sqrt(3.0d0)

               gpcord(1,1) = -sq3
               gpcord(1,2) = -sq3
               gpweight(1) = 1.0d0

               gpcord(2,1) =  sq3
               gpcord(2,2) = -sq3
               gpweight(2) = 1.0d0

               gpcord(3,1) =  sq3
               gpcord(3,2) =  sq3
               gpweight(3) = 1.0d0

               gpcord(4,1) = -sq3
               gpcord(4,2) =  sq3
               gpweight(4) = 1.0d0

               call fem_shapefunctions_2r()

            endif
c...........................................................
         else if(macro1.eq.'stressboun') then
            ilithod = 0
            ilithgrad = 0	  
            read (inpt, '(a80)') wdd1
            if (wdd1(1:9) .eq. 'lithograd') then     
               ilithgrad = 1	  
               read(wdd1,*) sdum1, sdepth, gdepth
            else if (wdd1(1:11) .eq. 'distributed') then
               iforce = 1
            else if (wdd1(1:11) .eq. 'lithostatic') then
c     if lithostatic then must use initcalc	   
               ilitho = 1
               ilithod = 1
               initcalc = 1
               if(iout.ne.0) write(iout, 101)   
               if(iptty.ne.0) write(iptty, 101)	    
               if(.not.allocated(flitho)) then
                  allocate(flitho(n0,3))
                  flitho = 0.0d0
               endif
 101           format ('initcalc set because lithostatic chosen')
            else
               backspace inpt
               iforce = 0
            endif 
            
            allocate (stressboun(n0))
            allocate (kq_dum(n0))
            kq_dum = 0
            igroup = 1
            narrays = 2
            itype(1) = 8
            itype(2) = 4
            default(1) = 0.
            default(2) = 0
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           r8_1 = stressboun(1:n0),i4_1 = kq_dum(1:n0)) 
            
            if(ilithgrad.ne.0) then    
c     
c     calculate lithostatic gradient
c     
               if(icnl.eq.0) then
c     assumes all gradients are in the vertical direction    
c     z is vertical in 3D (positive upward)     
                  do i = 1,n0
                     if(kq_dum(i).eq.1) then
c     set stress gradient 
                        kr(i,1) = -1
                        strxgrad = stressboun(i)
                        strx = (sdepth+gdepth-cord(i,3))*strxgrad
                        forc(i,1) = strx
                        iarea_str(i,1) = 1
                     endif		      
                     if(kq_dum(i).eq.2) then
c     set stress gradient 
                        kr(i,2) = -2
                        strxgrad = stressboun(i)
                        strx = (sdepth+gdepth-cord(i,3))*strxgrad
                        forc(i,2) = strx
                        iarea_str(i,2) = 2
                     endif		      
                     if(kq_dum(i).eq.3) then
c     set stress gradient 
                        kr(i,3) = -3
                        strxgrad = stressboun(i)
                        strx = (sdepth+gdepth-cord(i,3))*strxgrad
                        forc(i,3) = strx
                        iarea_str(i,3) = 3
                     endif
                  enddo  
               else
c     assumes all gradients are in the vertical direction    
c     y is vertical in 2D  
                  do i = 1,n0
                     if(kq_dum(i).eq.1) then
c     set stress gradient 
                        kr(i,1) = -1
                        strxgrad = stressboun(i)
                        strx = (sdepth+gdepth-cord(i,2))*strxgrad
                        forc(i,1) = strx
                        iarea_str(i,1) = 1
                     endif		      
                     if(kq_dum(i).eq.2) then
c     set stress gradient 
                        kr(i,2) = -2
                        strxgrad = stressboun(i)
                        strx = (sdepth+gdepth-cord(i,2))*strxgrad
                        forc(i,2) = strx
                        iarea_str(i,2) = 2
                     endif
                  enddo
               endif
            else if(ilithod.eq.0) then  
c     
c     apply displacements and forces (distributed forces 
c     are completed later - see iflg=2)
c     
               do i = 1,n0
                  if(kq_dum(i).eq.1) then
c     displacement
                     kr(i,1) = 1
                     disp(i,1) = stressboun(i)
                     disp0(i,1) = disp(i,1)
                  else if(kq_dum(i).eq.-1) then
c     force
                     if(iforce.eq.0) then
                        kr(i,1) = -1
                     else
                        kr(i,1) = -1
                        idum_str(i,1) = izonef(i)
                     endif
                     forc(i,1) = stressboun(i)
                  endif		      
                  if(kq_dum(i).eq.2) then
c     displacement
                     kr(i,2) = 2
                     disp(i,2) = stressboun(i)
                     disp0(i,2) = disp(i,2)
                  else if(kq_dum(i).eq.-2) then
c     force
                     if(iforce.eq.0) then
                        kr(i,2) = -2
                     else
                        kr(i,2) = -2
                        idum_str(i,2) = izonef(i)
                     endif
                     forc(i,2) = stressboun(i)
                  endif   
                  if(kq_dum(i).eq.3) then
c     displacement
                     kr(i,3) = 3
                     disp(i,3) = stressboun(i)
                     disp0(i,3) = disp(i,3)
                  else if(kq_dum(i).eq.-3) then
c     force
                     if(iforce.eq.0) then
                        kr(i,3) = -3
                     else
                        kr(i,3) = -3
                        idum_str(i,3) = izonef(i)
                     endif
                     forc(i,3) = stressboun(i)
                  endif      		    		     
               enddo
            else 
c     new code for lithostatic and principal stresses           
               do i = 1,n0
                  if(kq_dum(i).eq.1) then
c     multiplier for lithostatic load
                     flitho(i,1) = stressboun(i)
                  endif		      
                  if(kq_dum(i).eq.2) then
c     multiplier for lithostatic load
                     flitho(i,2) = stressboun(i)
                  endif   
                  if(kq_dum(i).eq.3) then
c     multiplier for lithostatic load
                     flitho(i,3) = stressboun(i)
                  endif      		    		     
               enddo          
            endif
            deallocate(stressboun,kq_dum)
         else if(macro1.eq.'elastic   ') then
            
            igroup = 1
            narrays = 2
            itype(1) = 8
            itype(2) = 8
            default(1) = 0.
            default(2) = 0.
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           r8_1 = elastic_mod(1:n0),r8_2 = poisson(1:n0))
c     s kelkar 4/20/2010 E and Nue as functions of temperature
         else if(macro1.eq.'nonlinear '.or.macro1.eq.'elastic no') then
c
c  read in  tables or functions and attributes
c
            isNonlinear = 1
            allocate(i_tab_youngs(max_y_tab))
            allocate(iy_tab(n0))
            iy_tab = 0
            allocate (e_temp91(nentries_young_max,3,max_y_tab))
            e_temp91 = 0.0d0
            allocate(e_ini(max_non_str))
            allocate(poisson_ini(max_non_str))
            allocate(t_non_ref(max_non_str))
            allocate(dEdt(max_non_str))
            allocate(dNuedt(max_non_str))
            allocate(istr_non_model(max_non_str))
            e_ini = 0.0d0
            poisson_ini = 0.0d0
            dEdt = 0.0d0
            t_non_ref = 0.0d0
            dNuedt = 0.0d0
            num_tab = 0
 211        continue
            read(inpt,'(a80)') wdd1
            if(.not. null1(wdd1)) then
c gaz 042916 need to put on one line 
               backspace inpt
               read(inpt,*) i_non_strs
               backspace inpt
               num_tab = num_tab + 1
               istr_non_model(num_tab) = i_non_strs
              if(i_non_strs.eq.1)then
               read(inpt,*) i_non_strs,e_ini(num_tab),
     &          poisson_ini(num_tab),t_non_ref(num_tab),dEdt(num_tab),
     &          dNuedt(num_tab)                   
              else
c  read nonlinear youngs table (function of T)            
                 call young_temp_table(0,0,num_tab)
              endif
               go to 211
            endif
            go to 210
c     read in nodal capillary type
 210         continue
c assign tables to nodes
            narrays = 1
            itype(1) = 4
            default(1) = 0
            macro = "etab"
            igroup = 2
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(9), macro, igroup, ireturn,
     2           i4_1=iy_tab(1:n0) )  
c         endif       
c     s kelkar 12/6/09 axisymmetric anisotropy
         else if(macro1.eq.'anisotropy') then
            stress_anisotropy_in = .true.
            stress_anisotropy_use = .true.
            allocate (elastic_mod_t(n0))
            allocate (poisson_t(n0))
            allocate (ezz(n0))
            allocate (shearmod_t(n0))
            allocate (e4(n0))
            elastic_mod_t=0.0d0
            poisson_t=0.0d0
            ezz=0.0d0
            shearmod_t=0.0d0
            e4=0.0d0
            igroup = 1
            narrays = 3
            itype(1) = 8
            itype(2) = 8
            itype(3) = 8
            default(1) = 0.
            default(2) = 0.
            default(3) = 0.
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           r8_1 = elastic_mod_t(1:n0),r8_2 = poisson_t(1:n0),
     &           r8_3 = shearmod_t(1:n0))
c....................................................................
         else if(macro1.eq.'plastic   ') then

c     optional parameter for increamental application of displacement BC
            n_incr=0
            read (inpt, '(a10)') macro2
            if(macro2(1:4).eq.'incr') then
               read(inpt,*)n_incr
            else
               backspace(inpt)
            endif
c     
            iPlastic = 1
            flag_element_perm = 1
            write(iout,*)'................'
            write(iout,*)'plastic model, forcing the use of elemental'
            write(iout,*)'stress-dependant permeability updates'
            write(iout,*)'..................'

            if(ifem.eq.1) then
               allocate(isPlastic(nei, numgausspoints))
               allocate(plastic_strain(nei, numgausspoints))
               allocate(conv_pstrain(nei, numgausspoints))
               allocate(pstrain(1:n0))

               isPlastic = 0
               plastic_strain = 0.0d0
               pstrain = 0.0d0
            endif

c     Define the different flags for the plasticctr subroutine
            initPlastic = 1
            assemblePlastic = 2
            call plasticctr(initPlastic)
            NumPlasticModels = 1
            read(inpt,*) NumPlasticModels
c forcing to have only 1 plastic model for the whole domain
            NumPlasticModels = 1

c            if (iout .ne. 0) write(iout, 102) NumPlasticModels
c 102        format ('Number of plastic models being used : ', i4)

            if(.not.allocated(plasticModel)) then
               allocate(plasticModel(1:NumPlasticModels))
               allocate(modulus(1:NumPlasticModels))
               allocate(nu(1:NumPlasticModels))
               allocate(plasticParam1(1:NumPlasticModels))
               allocate(plasticParam2(1:NumPlasticModels))
               allocate(plasticParam3(1:NumPlasticModels))
            endif
            
            if(.not.allocated(modelNumber)) then
               allocate(modelNumber(1:n0))
            endif

            i = 0
            j = 0
            ex = .false.
            do i=1,NumPlasticModels
               read(inpt,*) itmp
               backspace inpt
               
               if(itmp.eq.1) then
c     Isotropic, linear elastic solid
                  read(inpt,*) plasticModel(i), modulus(i), nu(i)
               else if (itmp.eq.2) then                          
c     von-Mises plasticity model
c     plasticParam1 = yield stress
c     plasticParam2 = hardening parameter
                  read(inpt,*) plasticModel(i), modulus(i), nu(i),
     &                 plasticParam1(i), plasticParam2(i)
               else if (itmp.eq.3) then
                  ! Drucker-Prager plasticity model (without cap)
                  ! plasticParam1 = eta
                  ! plasticParam2 = xi
                  ! plasticParam3 = c
                  read(inpt,*) plasticModel(i), modulus(i), nu(i),
     &                 plasticParam1(i), plasticParam2(i),
     &                 plasticParam3(i)             
               endif
               
               if (iout .ne. 0) write(iout, 103) i, NumPlasticModels,
     &              plasticModel(i)
 103           format ('Model number ', i4, ' of ', i4,
     &              ' is plastic model ', i4)
            end do
c     set default
            modelNumber = 1
c forcing the entire domain to have the same plastic model

c assign model# to each node as specified in the input
c            narrays = 1
c            itype(1) = 4
c            default(1) = 1
c            igroup = 1
c            call initdata2( inpt, ischk, n0, narrays,
c     &           itype, default, macroread(8), macro, igroup, ireturn,
c     &           i4_1=modelNumber(1:n0))ccc
c
            do i=1,n0
               elastic_mod(i) = modulus(modelNumber(i))
               poisson(i) = nu(modelNumber(i))
            enddo
c....................................................................
         else if(macro1(1:10).eq.'biot      ') then
            igroup = 1
            narrays = 2
            itype(1) = 8
            itype(2) = 8
            default(1) = 0.
            default(2) = 0.
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           r8_1 = alp0(1:n0),r8_2 = bulk0(1:n0))            
c....................................................................
       else if(macro1(1:10).eq.'nonlinbiot'.or.
     &         macro1(1:10).eq.'biot nonli') then
c
c  read in  tables of biot parameters as function of temperature
c
            isbiotNonLin = 1
            allocate(i_tab_biot(max_y_tab_b))
            allocate(iy_tab_biot(n0))
            iy_tab_biot = 0
            allocate (biot_temp91(nentries_biot_max,3,max_y_tab_b))
            biot_temp91 = 0.0d0
            allocate(istr_non_model_biot(max_non_str_biot))
            num_tab = 0
 311        continue
            read(inpt,'(a80)') wdd1
            if(.not. null1(wdd1)) then
c gaz 042916 need to put on one line 
               backspace inpt
               read(inpt,*) i_non_strs
               backspace inpt
               num_tab = num_tab + 1
               istr_non_model_biot(num_tab) = i_non_strs
c  read nonlinear biot table (function of T)            
                 call biot_temp_table(0,0,num_tab)
               go to 311
            endif
            go to 310
c     read in nodal capillary type
 310         continue
c assign tables to nodes
            narrays = 1
            itype(1) = 4
            default(1) = 0
            macro = "etab"
            igroup = 2
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(9), macro, igroup, ireturn,
     2           i4_1=iy_tab_biot(1:n0) )          
            
            

         else if(macro1(1:5).eq.'toler') then
            read(inpt,*) tol_stress
c     
c     print out history plot for displacement and stresses 
c     
         else if(macro1.eq.'stresshis') then
c     
            do jj = 1,1000000
               read(inpt,'(a80)') wdd1
               if(null1(wdd1)) go to 37
            enddo
 37         continue
            allocate(nskw_stress(jj-1,2))
            do i = 1,jj
               backspace inpt
            enddo	
            istrshis = jj - 1  
            do i = 1,jj
               read(inpt,'(a80)') wdd1
               if(null1(wdd1)) go to 39
               read(wdd1,'(a4)') stype
               if(stype.eq.'node') then
                  read(wdd1,*) stype, nodestress, ptype
               else
                  read(wdd1,*) stype, xi, yi, zi, ptype
                  call near3 (xi, yi, zi, nodestress, 0)
               endif
               nskw_stress(i,1) = nodestress
               if(ptype.eq.'disx  ')then
                  nskw_stress(i,2) = 1
               else if(ptype.eq.'disy  ')then
                  nskw_stress(i,2) = 2
               else if(ptype.eq.'disz  ')then
                  nskw_stress(i,2) = 3
               else if(ptype.eq.'strx  ')then
                  nskw_stress(i,2) = 4
               else if(ptype.eq.'stry  ')then
                  nskw_stress(i,2) = 5
               else if(ptype.eq.'strz  ')then
                  nskw_stress(i,2) = 6
               else if(ptype.eq.'strxy ')then
                  nskw_stress(i,2) = 7	
               else if(ptype.eq.'strxz ')then
                  nskw_stress(i,2) = 8
               else if(ptype.eq.'stryz ')then
                  nskw_stress(i,2) = 9	            
               endif
            enddo
 39         continue	
         elseif(macro1.eq.'principal') then
            flag_principal=1
            if (.not.allocated(xtan_min)) allocate (xtan_min(n0))
         elseif(macro1.eq.'excess_she') then
            read(inpt,*) friction_out,strength_out, pp_fac_out
            flag_excess_shear=1
            if (.not.allocated(excess_shear)) then
               allocate (excess_shear(n0))
               excess_shear = 0.
               allocate (perm_mult1(n0))
               perm_mult1 = 1.
               allocate (perm_mult2(n0))
               perm_mult2 = 1.
               allocate (perm_mult3(n0))
               perm_mult3 = 1.
            endif
            if (.not.allocated(shear_angle)) then
               allocate (shear_angle(n0))
               shear_angle = 0.
            endif
         elseif(macro1.eq.'strainout') then
            flag_strain_output = 1
            ifile_strain = open_file('strain.out', 'unknown')

c     end input macros
c     
         else
            write (ierr, 104)
            if (iout .ne. 0) write(iout, 104) 
            if (iptty .ne. 0) write(iptty, 104)
            stop
 104        format ('ERROR IN STRESS INPUT(STOPPING)')
         end if
 40      continue
         go to 100
 200     continue
c..............................................................
c     
c     linear isotropic or anisotropic (at present plain strain and 3D)
c     plain stress has different combinations
c     
        call elastic_constants(1) 
         
         macroread(8) = .TRUE. 
         if(iptty.ne.0) then
            write(iptty, 105)
            write(iptty, 106)
            write(iptty, 107)
            write(iptty, 108)
            write(iptty, 109)
            write(iptty, 110)
            write(iptty, 111)
            write(iptty, 105)
         endif
         if (iout .ne. 0) then
            write(iout, 105)
            write(iout, 106)
            write(iout, 107)
            write(iout, 108)
            write(iout, 109)
            write(iout, 110)
            write(iout, 111)
            write(iout, 105)
         end if
 105     format ('*********************************')
 106     format (' NOTE displacements are zeroed after every time step')
 107     format (' NOTE macro "nobr" always set with stress (disabled)')
 108     format (' NOTE macro "isot" not allowed with stress')
 109     format (' NOTE input changed (thermal expansion now volume ',
     &        '(=linear*3)')
 110     format (' NOTE input changed (Biot term now = a in the term ',
     &        'a*(1-K/Kb)')
 111     format (' NOTE ifinv=1 (finite volume disabled when stress ',
     &        'calls present)')
         inobr = 1
         ivf = 0
         mlz = 0
c
c check that all elastic parameters are defined
c  
         j = 0
         do i = 1, n0
          if(elastic_mod(i). gt. 0.1) then
           go to 811
          elseif (isNonlinear.eq.1) then
           if(iy_tab(i).ne.0) go to 811
          else
           j = j + 1
          endif
811      continue
         enddo   
        if(j.gt.0) then
         if(iout.ne.0) write(iout,812) j
         if(iptty.ne.0) write(iptty,812) j
         if(iout.ne.0) write(iout,105) 
         if(iptty.ne.0) write(iptty,105) 
        endif  
812     format(' Warning: ',i8,' nodes without elastic params defined')
      else if(iflg.eq.2) then
c     
c     sort out applied forces            
c     
         do i = 1,n0
            if(idum_str(i,1).ne.0) then
               dum_str(idum_str(i,1),1) = dum_str(idum_str(i,1),1) + 
     &              sx1(i)
            endif
            if(idum_str(i,2).ne.0) then
               dum_str(idum_str(i,2),2) = dum_str(idum_str(i,2),2) + 
     &              sx1(i)
            endif
            if(idum_str(i,3).ne.0) then
               dum_str(idum_str(i,3),3) = dum_str(idum_str(i,3),3) + 
     &              sx1(i)
            endif
         enddo
         
         do i = 1,n0
            if(idum_str(i,1).ne.0) then
               forc(i,1) = forc(i,1)*sx1(i)/dum_str(idum_str(i,1),1)
            endif
            if(idum_str(i,2).ne.0) then
               forc(i,2) = forc(i,2)*sx1(i)/dum_str(idum_str(i,2),2)
            endif
            if(idum_str(i,3).ne.0) then
               forc(i,3) = forc(i,3)*sx1(i)/dum_str(idum_str(i,3),3)
            endif
         enddo
c     
c     sort out applied stress (turn into forces)             
c     
c     
         call geom_stress(1)
c     before loop, stress at depth - after loop, applied force    
c     
         do i = 1,n0
            if(iarea_str(i,1).ne.0) then
               forc(i,1) = forc(i,1)*area_str(i,1)
            endif
            if(iarea_str(i,2).ne.0) then 
               forc(i,2) = forc(i,2)*area_str(i,2)
            endif
            if(iarea_str(i,3).ne.0) then
               forc(i,3) = forc(i,3)*area_str(i,3)
            endif
         enddo
         
         deallocate(idum_str,dum_str,iarea_str,area_str)
         
      else if(iflg.eq.3) then
c     
c     calculate boundary conditions based on lithostatic stresses
c     need to be called after the stress calculation
c     flitho(i,1) contains the fraction of lithostatic 
c     remove fixed zero displacement
c     
         if(ilitho.ne.0) then
            if(icnl.eq.0) then
               do i = 1,n0
                  if(flitho(i,1).ne.0.0.or.flitho(i,2).ne.0.0
     &                 .or.flitho(i,3).ne.0.0) then
                     grav_save = grav
                     grav = 0.0
                     call geneq_stress_uncoupled_3D(i)
                     grav = grav_save
                     if(flitho(i,1).ne.0.0) then
                        kr(i,1) = -1
                        str_ratio_xz = flitho(i,1)
                        forc(i,1) = bp(i+nrhs(1))*str_ratio_xz
                     endif
                     if(flitho(i,2).ne.0.0) then
                        kr(i,2) = -2
                        str_ratio_yz = flitho(i,2)
                        forc(i,2) = bp(i+nrhs(2))*str_ratio_yz
                     endif
                     if(flitho(i,3).ne.0.0) then
                        kr(i,3) = -3
                        forc(i,3) = bp(i+nrhs(3))*flitho(i,3)
                     endif
                  endif
               enddo 
            else
               do i = 1,n0
                  if(flitho(i,1).ne.0.0.or.flitho(i,2).ne.0.0) then
                     call geneq_stress_uncoupled_2D(i)
                     if(flitho(i,1).ne.0.0) then
                        kr(i,1) = -1
                        str_ratio_xy = flitho(i,1)
                        forc(i,1) = bp(i+nrhs(1))*str_ratio_xy
                     endif
                     if(flitho(i,2).ne.0.0) then
                        kr(i,2) = -2
                        forc(i,2) = bp(i+nrhs(2))*flitho(i,2)
                     endif
                  endif
               enddo        
            endif
            deallocate(flitho)  
         endif   
c     
c     If Bai model(model 7) is chosen,store the initial effective stresses
c     If  model 2 is chosen, store the initial effective stresses
c     If  model 6 is chosen, store the initial effective stresses
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(ipermstr7.ne.0.or.ipermstr5.ne.0.or.ipermstr8.ne.0) then
            if (.not. allocated(check)) allocate(check(n0))
         endif

         if(ipermstr2.ne.0.or.ipermstr6.ne.0.or.ipermstr7.ne.0) then
            if(icnl.ne.0) then
               if(.not.allocated(estr_x0)) then
               allocate(estr_x0(n0))
               allocate(estr_y0(n0))
               endif
               if(.not.allocated(str_xy0)) then
               allocate(str_xy0(n0))
               endif
            elseif(icnl.eq.0) then
               if(.not.allocated(estr_x0)) then
               allocate(estr_x0(n0))
               allocate(estr_y0(n0))
               allocate(estr_z0(n0))           
               endif
               if(.not.allocated(str_xy0)) then
               allocate(str_xy0(n0))  
               allocate(str_xz0(n0))
               allocate(str_yz0(n0))
            endif
            endif
            if(icnl.ne.0) then  
               do i = 1,n0
                  estr_x0(i) = str_x(i)-pho(i)
                  estr_y0(i) = str_y(i)-pho(i)
                  str_xy0(i) = str_xy(i)
               enddo
            else
               do i = 1,n0
                  estr_x0(i) = str_x(i)-pho(i)
                  estr_y0(i) = str_y(i)-pho(i)
                  estr_z0(i) = str_z(i)-pho(i)
                  str_xy0(i) = str_xy(i)
                  str_xz0(i) = str_xz(i)
                  str_yz0(i) = str_yz(i)
               enddo
            endif
         endif
         
         if(ipermstr8.ne.0 .and. .not. (allocated(es_f_x0))) then
            if(icnl.ne.0) then
               allocate(es_f_x0(n0,2))
               allocate(es_f_y0(n0,2))
               allocate(s_f_xy0(n0,2))
            else 
               allocate(es_f_x0(n0,3))
               allocate(es_f_y0(n0,3))
               allocate(es_f_z0(n0,3))           
               allocate(s_f_xy0(n0,3))  
               allocate(s_f_xz0(n0,3))
               allocate(s_f_yz0(n0,3))
               allocate(frc_zen(n0,3))
               allocate(frc_azm(n0,3))
            endif
         endif  
c     
c     if pore pressure damage model,  save the initial lithostsic stress
c     this next part is to remove stresses caused by initial temperature 
c     and pressure fields       
         
         if(ipermstr6.ne.0) then
c     
c     for effective stress we need the lowest pressure 
c     at highest elevation  
c     
            elev_high = -1.e30
            do i = 1,n0
               if(cord(i,igrav).gt.elev_high) then
                  elev_high = cord(i,igrav)
                  pres_elev = pho(i)
               endif
            enddo  
            pres_elev = pres_elev-0.1   
c     
c     analysis of lowest pressure
c     
            elev = pres_elev/rho1grav*1.e-3	
c     
 112        format ('>>>> Warning for effective stress (pressure)',
     &           ' usage <<<<')
 113        format ('Lowest reservoir pressure is: ', g16.9, 'MPa')
 114        format ('Equivalent depth (ie the top of reservoir) is: ',
     &           g16.9, ' m')
 115        format ('>>Should be consistent with lithostatic stress<<')
            if(iout.ne.0) then
               write(iout,*) 
               write(iout, 112) 
               write(iout, 113) pres_elev
               write(iout, 114) elev
               write(iout, 115)
            endif 
            if(iptty.ne.0) then
               write(iptty,*) 
               write(iptty, 112) 
               write(iptty, 113) pres_elev
               write(iptty, 114) elev
               write(iptty, 115)
            end if
            do i = 1,n0
c     stress due to temp and pore pressure changes
               dpd=phi(i)-phini(i)
               tdumt=t(i)-tini(i)
               biot=bulk(i)
               ctherm=alp(i)
               e1i = e1(i)
               e2i = e2(i)
               e3i = e3(i)
               erat=e2i/e1i
c     s kelkar 12/5/09 axisymmetric anisotropy.........................
               if(stress_anisotropy_use) then
                  ezzi= ezz(i)
                  e4i=e4(i)
                  efacxy= (e1i+e2i+e4i)
                  efacz = 2.0d0*e4i + ezzi
c     pore pressure terms. See Keita's notes dated 2/25/2010, Here 
c     efacxy is defined as efacxy=3*Hp=(C11+C12+C13) and  
c     bulk()=beta_p/3Hp or  beta_p=bulk*efacxy. Then
c     Ks=Hp/(1-beta_p)=(efacxy/3)/(1-efacxy*bulkb). then beta_t 
c     calculated from
c     betat=1-Ht/Ks=1-(efacz/3)/Ks where efacz=3*Ht=(2C13+C33)
c     also define efac_betat=betat/efacz for consistancy of notation
                  efac_ks=(efacxy/3.0)/(1.0-efacxy*biot)
                  betat=1.0-(efacz/3.0)/efac_ks
                  efac_betat=1.0/efacz-1.0/(3.0*efac_ks)
                  shtixy=efacxy*ctherm*tdumt
                  shpixy=efacxy*biot*dpd
                  shtiz=efacz*ctherm*tdumt
                  shpiz=efacz*efac_betat*dpd
                  estr_x0(i) = str_x(i)-pho(i)-shtixy-shpixy
                  estr_y0(i) = str_y(i)-pho(i)-shtixy-shpixy
                  estr_z0(i) = str_z(i)-pho(i)-shtiz-shpiz
c............................................................
               else
                  efac=3.d0*e2i+2.d0*e3i
                  epi=efac*biot
                  eti=efac*ctherm
                  shti=(eti*tdumt)
                  shpi=(epi*dpd)
                  if(icnl.eq.0) then
                     estr_x0(i) = str_x(i)-pho(i)-shti-shpi
                     estr_y0(i) = str_y(i)-pho(i)-shti-shpi
                     estr_z0(i) = str_z(i)-pho(i)-shti-shpi
                  else
                     estr_x0(i) = str_x(i)-pho(i)-shti-shpi
                     estr_y0(i) = str_y(i)-pho(i)-shti-shpi
                  endif
               endif

            enddo
         endif
         if(idisp_rel.ne.0) then
            du_ini = du
            dv_ini = dv
            if(icnl.eq.0) dw_ini = dw
         endif
c     iflg=4 calculate nonlinear material properties
      else if(iflg.eq.4) then 
         call stress_mech_props(0,Nonlin_model_flag,0)
         
c     iflg=5 calculate stress fluid interaction properties
      else if(iflg.eq.5) then
         call stress_fluid_mech_props(0,0)           
         
c     update volumetric strain
      else if(iflg.eq.-6) then 
c     vol_strain0 = vol_strain    
c     calculate volumetric strain
      else if(iflg.eq.6) then
c     
         if(icnl.eq.0) then
c     3-d volumetric change and volume strain
            neqp1 = neq+1
            k = nelm(neqp1)-neqp1
            
            if(.not.allocated(dvol_strainu)) then
               allocate(dvol_strainu(k))
               allocate(dvol_strainv(k))
               allocate(dvol_strainw(k))
               dvol_strainu = 0.0
               dvol_strainv = 0.0
               dvol_strainw = 0.0
            endif
            
c     vol_tot_change = 0.0
            do i = 1,neq
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               cordxi = cord(i,1)
               cordyi = cord(i,2)
               cordzi = cord(i,3)
               dispxi = du(i)-du_ini(i)
               dispyi = dv(i)-dv_ini(i)
               dispzi = dw(i)-dw_ini(i)
               area_ix = 0.0d0
               area_iy = 0.0d0
               area_iz = 0.0d0
               dvol_strainu_i = 0.0d0
               dvol_strainv_i = 0.0d0
               dvol_strainw_i = 0.0d0
               vol_change = 0.0
               md = nelmdg(i)
               do jj = i1,i2
c     area calculation 
                  j=istrw(jj-neqp1)
                  kb = nelm(jj)
                  if(kb.ne.i) then
                     delx = (cord(kb,1)-cordxi)
                     dely = (cord(kb,2)-cordyi)
                     delz = (cord(kb,3)-cordzi)
                     area_facex = -sx(j,1)*delx
                     area_facey = -sx(j,2)*dely
                     area_facez = -sx(j,3)*delz
                     area_ix = area_ix + area_facex
                     area_iy = area_iy + area_facey
                     area_iz = area_iz + area_facez
c     
c     dispx = 0.5d0*(dispxi + du_tot(kb)) 
c     dispy = 0.5d0*(dispyi + dv_tot(kb)) 	 
c     dispz = 0.5d0*(dispzi + dw_tot(kb))  
c     
                     dispx = 0.5d0*(dispxi + (du(kb)-du_ini(kb))) 
                     dispy = 0.5d0*(dispyi + (dv(kb)-dv_ini(kb))) 	 
                     dispz = 0.5d0*(dispzi + (dw(kb)-dw_ini(kb)))    
c     
                     vol_change = vol_change + area_facex*dispx + 
     &                    area_facey*dispy + area_facez*dispz 
                     if(kr(kb,1).ne.1) then
                        dvol_strainu(jj-neqp1) = area_facex*0.5d0
                     else
                        dvol_strainu(jj-neqp1) = 0.
                     endif
                     if(kr(kb,2).ne.2) then
                        dvol_strainv(jj-neqp1) = area_facey*0.5d0
                     else
                        dvol_strainv(jj-neqp1) = 0.
                     endif 
                     if(kr(kb,3).ne.3) then
                        dvol_strainw(jj-neqp1) = area_facez*0.5d0
                     else
                        dvol_strainw(jj-neqp1) = 0.
                     endif
                  endif
               enddo
               if(abs(area_ix).ge.area_tol) then
                  vol_change = vol_change - area_ix*dispxi
                  if(kr(i,1).ne.1) then
                     dvol_strainu(md-neqp1) = -area_ix*0.5d0
                  else
                     dvol_strainu(md-neqp1) = 0.0
                  endif
               endif
               if(abs(area_iy).ge.area_tol) then
                  vol_change = vol_change - area_iy*dispyi
                  if(kr(i,2).ne.2) then 
                     dvol_strainv(md-neqp1) = -area_iy*0.5d0
                  else
                     dvol_strainv(md-neqp1) = 0.0
                  endif
               endif 
               if(abs(area_iz).ge.area_tol) then
                  vol_change = vol_change - area_iz*dispzi
                  if(kr(i,3).ne.3) then
                     dvol_strainw(md-neqp1) = -area_iz*0.5d0
                  else
                     dvol_strainw(md-neqp1) = 0.0
                  endif
               endif 	 
               vol_strain(i) = vol_change/sx1(i)               
c.............................................................
c     vol_tot_change = vol_tot_change + vol_change
            enddo
c     vol_strain = 0.0
c     dvol_strainw = 0.0
c     dvol_strainv = 0.0
c     dvol_strainu = 0.0
         else
c     2-d implimentation
            
c     2-d volumetric change (unit thickness) and volume strain
            neqp1 = neq+1
            if(.not.allocated(dvol_strainu)) then
               k = nelm(neqp1)-neqp1
               allocate(dvol_strainu(k))
               allocate(dvol_strainv(k))
            endif
            do i = 1,neq
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               cordxi = cord(i,1)
               cordyi = cord(i,2)
               dispxi = du(i)-du_ini(i)
               dispyi = dv(i)-dv_ini(i)
               area_ix = 0.0d0
               area_iy = 0.0d0
               dvol_strainu_i = 0.0d0
               dvol_strainv_i = 0.0d0
c     vol_change = vol_strain0(i)*sx1(i)
               vol_change = 0.0
               md = nelmdg(i)
               do jj = i1,i2
c     area calculation 
                  j=istrw(jj-neqp1)
                  kb = nelm(jj)
                  if(kb.ne.i) then
                     delx = (cord(kb,1)-cordxi)
                     dely = (cord(kb,2)-cordyi)
                     area_facex = -sx(j,1)*delx
                     area_facey = -sx(j,2)*dely
                     area_ix = area_ix + area_facex
                     area_iy = area_iy + area_facey
c     
c     dispx = 0.5d0*(dispxi + du_tot(kb)) 
c     dispy = 0.5d0*(dispyi + dv_tot(kb)) 	   
c     
                     dispx = 0.5d0*(dispxi + (du(kb)-du_ini(kb))) 
                     dispy = 0.5d0*(dispyi + (dv(kb)-dv_ini(kb))) 	  
c     
                     vol_change = vol_change + area_facex*dispx + 
     &                    area_facey*dispy 
                     if(kr(kb,1).ne.1) then
                        dvol_strainu(jj-neqp1) = area_facex*0.5d0
                        dvol_strainu_i = dvol_strainu_i + 
     $                       area_facex*0.5d0
                     endif
                     if(kr(kb,2).ne.2) then
                        dvol_strainv(jj-neqp1) = area_facey*0.5d0
                        dvol_strainv_i = dvol_strainv_i + 
     $                       area_facey*0.5d0
                     endif 
                     
                  endif
               enddo
               if(abs(area_ix).ge.area_tol) then
                  vol_change = vol_change - area_ix*dispxi
                  if(kr(i,1).ne.1) then
                     dvol_strainu(md-neqp1) = dvol_strainu_i  - area_ix
                  else
                     dvol_strainu(md-neqp1) = 0.0
                  endif
               endif
               if(abs(area_iy).ge.area_tol) then
                  vol_change = vol_change - area_iy*dispyi
                  if(kr(i,2).ne.2) then 
                     dvol_strainv(md-neqp1) = dvol_strainv_i  - area_iy
                  else
                     dvol_strainv(md-neqp1) = 0.0
                  endif
               endif 
               
               vol_strain(i) = vol_change/sx1(i)
c     vol_tot_change = vol_tot_change + vol_change
            enddo
         endif
      else if(iflg.eq.-7) then
c     
c     porosity changes (associate with volume changes)
c     
         if(istresspor.ne.0) then
            do i = 1,neq
               por = psini(i)
               ps(i) = por*(1.0 + vol_strain(i))
            enddo
	 endif
         
c     
c     derivative of porosity wrt displacements
c     
      else if(iflg.eq.7) then
c     
c     porosity changes
c     
         if(icnl.eq.0) then
c     3-d volumetric change

            neqp1 = neq+1

            do i = 1,neq
               i1 = nelm(i)+1
               i2 = nelm(i+1) 
               dispxi = du(i)
               dispyi = dv(i)
               dispzi = dw(i)
               md = nelmdg(i)
               por = psini(i)
c     
c     define porosity multiplier term (mass equation)
c     
               denf = (denh(i) + deni(i)*dtot)/por/dtot
c     
c     define porosity multiplier term (energy equation)
c     
               cp_t=denr(i)*cpr(i)*t(i)
               denfe = (deneh(i) + denei(i)*dtot)
               ener_fl = (denfe - (1.0-por)*cp_t)/por
               denfe = (ener_fl-cp_t)/dtot
c     
               sx1d = sx1(i)
               por_new = por*(1.0 + vol_strain(i))
               ps(i) = por_new
c     identify the diagonal terms   
               jmia = md - neqp1     
               dporeeui = por*dvol_strainu(jmia)/sx1d
               dporeevi = por*dvol_strainv(jmia)/sx1d
               dporeewi = por*dvol_strainw(jmia)/sx1d
               if(kr(i,1).eq.1) dporeeui = 0.0
               if(kr(i,2).eq.2) dporeevi = 0.0
               if(kr(i,3).eq.3) dporeewi = 0.0
c     
c     must get nmat numbers correct (done in stress_combine)
c     
c     conservation of water mass equation (x-y-z directions)
c     
               a(jmia+nmat(5))=a(jmia+nmat(5))+sx1d*denf*dporeeui 
               a(jmia+nmat(6))=a(jmia+nmat(6))+sx1d*denf*dporeevi
               a(jmia+nmat(7))=a(jmia+nmat(7))+sx1d*denf*dporeewi
c     
c     conservation of energy equation (x-y-z directions)
c     
               a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*denfe*dporeeui
               a(jmia+nmat(9))=a(jmia+nmat(9))+sx1d*denfe*dporeevi
               a(jmia+nmat(10))=a(jmia+nmat(10))+sx1d*denfe*dporeewi
c     
               do jj = i1,i2
                  kb = nelm(jj)
                  if(kb.ne.i) then
                     mdkb = jj-neqp1

                     dmeeukb = por*dvol_strainu(mdkb)/sx1d
                     dmeevkb = por*dvol_strainv(mdkb)/sx1d
                     dmeewkb = por*dvol_strainw(mdkb)/sx1d
                     if(kr(kb,1).eq.1) dmeeukb = 0.0
                     if(kr(kb,2).eq.2) dmeevkb = 0.0
                     if(kr(kb,3).eq.3) dmeewkb = 0.0
c     conservation of mass equation
                     a(mdkb+nmat(5))=a(mdkb+nmat(5))+sx1d*denf*dmeeukb 
                     a(mdkb+nmat(6))=a(mdkb+nmat(6))+sx1d*denf*dmeevkb 
                     a(mdkb+nmat(7))=a(mdkb+nmat(7))+sx1d*denf*dmeewkb 
c     conservation of energy equation 
                     a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*denfe*dmeeukb
                     a(jmia+nmat(9))=a(jmia+nmat(9))+sx1d*denfe*dmeevkb
                     a(jmia+nmat(10))=a(jmia+nmat(10))+sx1d*denfe*
     &                    dmeewkb

                  endif
               enddo
            enddo
c     
c     apply boundary conditions (fixed displacements)
c     
c     call stress_boun3(3,0) (accounted for above)
c     
         else
c     2-d volumetric change

            neqp1 = neq+1

            do i = 1,neq
               i1 = nelm(i)+1
               i2 = nelm(i+1) 
               dispxi = du(i)
               dispyi = dv(i)

               md = nelmdg(i)
               por = psini(i)
c     
c     define porosity multiplier term (mass equation)
c     
               denf = (denh(i) + deni(i)*dtot)/por/dtot
c     
c     define porosity multiplier term (energy equation)
c     
               cp_t=denr(i)*cpr(i)*t(i)
               denfe = (deneh(i) + denei(i)*dtot)
               ener_fl = (denfe - (1.0-por)*cp_t)/por
               denfe = (ener_fl-cp_t)/dtot
c     
               sx1d = sx1(i)
               por_new = por*(1.0 + vol_strain(i))
               ps(i) = por_new
c     identify the diagonal terms   
               jmia = md - neqp1     
               dporeeui = por*dvol_strainu(jmia)/sx1d
               dporeevi = por*dvol_strainv(jmia)/sx1d

               if(kr(i,1).eq.1) dporeeui = 0.0
               if(kr(i,2).eq.2) dporeevi = 0.0

c     
c     must get nmat numbers correct (done in stress_combine)
c     
c     conservation of water mass equation (x-y-z directions)
c     
               a(jmia+nmat(5))=a(jmia+nmat(5))+sx1d*denf*dporeeui 
               a(jmia+nmat(6))=a(jmia+nmat(6))+sx1d*denf*dporeevi

c     
c     conservation of energy equation (x-y-z directions)
c     
               a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*denfe*dporeeui
               a(jmia+nmat(9))=a(jmia+nmat(9))+sx1d*denfe*dporeevi

c     
               do jj = i1,i2
                  kb = nelm(jj)
                  if(kb.ne.i) then
                     mdkb = jj-neqp1

                     dmeeukb = por*dvol_strainu(mdkb)/sx1d
                     dmeevkb = por*dvol_strainv(mdkb)/sx1d

                     if(kr(kb,1).eq.1) dmeeukb = 0.0
                     if(kr(kb,2).eq.2) dmeevkb = 0.0

c     conservation of mass equation
                     a(mdkb+nmat(5))=a(mdkb+nmat(5))+sx1d*denf*dmeeukb 
                     a(mdkb+nmat(6))=a(mdkb+nmat(6))+sx1d*denf*dmeevkb 

c     conservation of energy equation 
                     a(jmia+nmat(8))=a(jmia+nmat(8))+sx1d*denfe*dmeeukb
                     a(jmia+nmat(9))=a(jmia+nmat(9))+sx1d*denfe*dmeevkb

                  endif
               enddo
            enddo
         endif         
c     generate derivatives wrt p and T for stresss eqs
      else if(iflg.eq.15) then   
c     return
         neqp1 = neq+1
	 do i = 1,neq
            e1i = e1(i)
            e2i = e2(i)
            e3i = e3(i)
            if(stress_anisotropy_use) then
               ezzi = ezz(i)
               e4i = e4(i)
            endif
            alpi=alp(i)
            bulki=bulk(i)
            i1 = nelm(i)+1
            i2 = nelm(i+1)
            do jj = i1,i2
               kb = nelm(jj)
               jmia = jj-neqp1
               iws = istrws(jj-neqp1)

c     x term for pore pressure and thermal expansion term
               sjsix=sxs(iws,7)
c     y term for pore pressure and thermal expansion term
               sjsiy=sxs(iws,8)
c     z term for pore pressure and thermal expansion term
               sjsiz=sxs(iws,9)        
               e1kb = e1(kb)
               e2kb = e2(kb)
               e3kb = e3(kb)
               e1bar=2.*e1i*e1kb/(e1i+e1kb + dis_tol)
               e2bar=2.*e2i*e2kb/(e2i+e2kb + dis_tol)
               e3bar=2.*e3i*e3kb/(e3i+e3kb + dis_tol)
               alpkb=alp(kb)
               alphab=2.*alpi*alpkb/(alpi+alpkb + dis_tol)
c     biot term
               bulkkb=bulk(kb)
               bulkb=2.*bulkkb*bulki/(bulkkb+bulki + dis_tol)
               efac = 3.d0*e2bar + 2.d0*e3bar

               
               tdumt=t(kb)-tini(kb)
               pdumt=phi(kb)-phini(kb)
c     
c     tdumx=sjsix*(tdumt*alphab+pdumt*bulkb)*efac
c     tdumy=sjsiy*(tdumt*alphab+pdumt*bulkb)*efac
c     tdumz=sjsiz*(tdumt*alphab+pdumt*bulkb)*efac
c     s kelkar 12/5/09 axisymmetric anisotropy.........................
               if(stress_anisotropy_use) then
                  ezzkb = ezz(kb)
                  e4kb = e4(kb)
                  ezzbar= 2.*ezzi*ezzkb/(ezzi+ezzkb + dis_tol)
                  e4bar= 2.*e4i*e4kb/(e4i+e4kb + dis_tol)
                  efacxy= (e1bar+e2bar+e4bar)
                  efacz = 2.0d0*e3bar + ezzbar
c     pore pressure terms. See Keita's notes dated 2/25/2010, Here 
c     efacxy is defined as efacxy=3*Hp=(C11+C12+C13) and  
c     bulk()=beta_p/3Hp or  beta_p=bulk*efacxy. Then
c     Ks=Hp/(1-beta_p)=(efacxy/3)/(1-efacxy*bulkb). then beta_t 
c     calculated from
c     betat=1-Ht/Ks=1-(efacz/3)/Ks where efacz=3*Ht=(2C13+C33)
c     also define efac_betat=betat/efacz for consistancy of notation
                  efac_ks=(efacxy/3.0)/(1.0-efacxy*bulkb)
                  betat=1.0-(efacz/3.0)/efac_ks
                  efac_betat=1.0/efacz-1.0/(3.0*efac_ks)
c     termt = efac*alphab
c     termp = efac*bulkb
                  a(jmia+nmat(10))=a(jmia+nmat(10))+sjsix*efacxy*bulkb 
                  a(jmia+nmat(12))=a(jmia+nmat(12))+sjsiy*efacxy*bulkb 
                  a(jmia+nmat(14))=a(jmia+nmat(14))+sjsiz*efacz*
     &                 efac_betat
                  a(jmia+nmat(11))=a(jmia+nmat(11))+sjsix*efacxy*alphab 
                  a(jmia+nmat(13))=a(jmia+nmat(13))+sjsiy*efacxy*alphab 
                  a(jmia+nmat(15))=a(jmia+nmat(15))+sjsiz*efacz*alphab
c.....................................................................
               else
                  termt = efac*alphab
                  termp = efac*bulkb
                  a(jmia+nmat(10))=a(jmia+nmat(10))+sjsix*termp 
                  a(jmia+nmat(12))=a(jmia+nmat(12))+sjsiy*termp 
                  a(jmia+nmat(14))=a(jmia+nmat(14))+sjsiz*termp  
                  
                  a(jmia+nmat(11))=a(jmia+nmat(11))+sjsix*termt 
                  a(jmia+nmat(13))=a(jmia+nmat(13))+sjsiy*termt 
                  a(jmia+nmat(15))=a(jmia+nmat(15))+sjsiz*termt
               endif     

            enddo
	 enddo
c     
c     apply boundary conditions (fixed displacements)
c     
         call stress_boun3(2,0)
c     
      else if(iflg.eq.8) then
c     generate thermo, equations and call solver
         if(icnl.ne.0)  then
            call gensl_stress_2D
         else 
            
            call gensl_stress_3D
            
         endif
 3333    continue
         
         

c     update the variables
      else if(iflg.eq.9) then
         if(idof_stress.eq.5) then

c     strd is passed through common
            nr1=nrhs(4)
            nr2=nrhs(5)
            do i=1,neq
               i1=i+nr1
               i2=i+nr2
               ieosd=ieos(i)
               if(ps(i).eq.0.0.or.ieosd.eq.0) then
                  t(i)=t(i)-bp(i2)*strd
               elseif(ieosd.eq.1) then
                  phi(i)=phi(i)-bp(i1)*strd
                  t(i)=t(i)-bp(i2)*strd
               elseif(ieosd.eq.2) then
                  phi(i)=phi(i)-bp(i1)*strd
                  s(i)=s(i)-bp(i2)*strd
               elseif(ieosd.eq.3) then
                  phi(i)=phi(i)-bp(i1)*strd
                  t(i)=t(i)-bp(i2)*strd
               endif
            enddo    
            bp_update = 0.0d0
            do i = 1,neq
               du(i) = du(i)-bp(i+nrhs(1))
               dv(i) = dv(i)-bp(i+nrhs(2))
               dw(i) = dw(i)-bp(i+nrhs(3))
               bp_update = max(bp_update, abs(bp(i+nrhs(1))),
     &              abs(bp(i+nrhs(2))),abs(bp(i+nrhs(3))))		 
            enddo
         else if(icnl.ne.0) then
            bp_update = 0.0d0
            do i = 1,neq
               du(i) = du(i)-bp(i+nrhs(1))
               dv(i) = dv(i)-bp(i+nrhs(2))
               bp_update = max(bp_update, abs(bp(i+nrhs(1))),
     &              abs(bp(i+nrhs(2))))			 
            enddo

         else if(icnl.eq.0) then
            bp_update = 0.0d0
            do i = 1,neq
               du(i) = du(i)-bp(i+nrhs(1))
               dv(i) = dv(i)-bp(i+nrhs(2))
               dw(i) = dw(i)-bp(i+nrhs(3))
               bp_update = max(bp_update, abs(bp(i+nrhs(1))),
     &              abs(bp(i+nrhs(2))),abs(bp(i+nrhs(3))))		
            enddo

         endif    
c     total displacements
      else if(iflg.eq.10) then
	 if(icnl.eq.0) then
            du_tot = du_tot+du
            dv_tot = dv_tot+dv
            dw_tot = dw_tot+dw  
	 else
            du_tot = du_tot+du
            dv_tot = dv_tot+dv	 
	 endif  
c     reset displacements
      else if(iflg.eq.-10) then
	 if(icnl.eq.0) then
c     du = 0.0
c     dv = 0.0
c     dw = 0.0
            du = duo
            dv = dvo
            dw = dwo
	 else
            du = duo
            dv = dvo
	 endif  
	 vol_strain = vol_strain0      

c     call stress output subroutines
      else if(iflg.eq.11)then
         iad_strs    =  max0( 1,iad_strs )
         aiter  =  dfloat( itert_s )/dfloat( iad_strs )
         aminkt =  neq

         if(iptty.ne.0) then 
            write(iptty,6116) 
            write(iptty,6117) 
            write(iptty,773)
            write(iptty,75) iad_strs
            write(iptty,76) aiter
            write(iptty,704) itotal_s,itotals_s
         endif
         if(ntty.eq.2) then
            write(iout,6116) 
            write(iout,6117) 
            write(iout,773)
            write(iout,75) iad_strs
            write(iout,76) aiter
            write(iout,704) itotal_s,itotals_s          
         endif
         if(icnl.ne.0) then            
c     
****  printout displacement information   ****
c     
            ibp_stress=max(ibp_stress,1)
            if(tol_stress.ne.0.0d0) then 
               if(mlz.eq.-1) then
                  if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) 
     &                 '>>> stress solution did not reach tolerance <<<'
               endif	 
               if(ntty.eq.2) write(iout,6118) ibp_stress,
     &              bp(ibp_stress+nrhs(1)),bp(ibp_stress+nrhs(2)),
     &              0.
            endif
            if (iout .ne. 0) then
            if(ntty.eq.2) write(iout,*) ' ' 
            if(ntty.eq.2) write(iout,*) 'displacements and forces ' 
            if(ntty.eq.2) write(iout,6119) 
            end if
            if(tol_stress.ne.0.0d0) then 
               if(mlz.eq.-1) then
                  if(iptty.gt.0) write(iptty,*) 
     &                 '>>> stress solution did not reach tolerance <<<'
               endif
               if(iptty.gt.0) write(iptty,6118) ibp_stress,
     &              bp(ibp_stress+nrhs(1)),bp(ibp_stress+nrhs(2)),
     &              0.
            endif
            if(iptty.gt.0) write(iptty,*) ' ' 
            if(iptty.gt.0) write(iptty,*)'displacements and forces '
            if(iptty.gt.0) write(iptty,6119) 
            
c     
            do   i=1,m
               md     =  nskw(i)               
               if(ntty.eq.2 .and. iout .ne. 0) write(iout,6120)  
     &              md, du(md), dv(md), 0., (forc(md,j),j=1,2),0.,
     &              vol_strain(md)
               if(iptty.gt.0) write(iptty,6120) 
     &              md, du(md), dv(md), 0., (forc(md,j),j=1,2), 0.,
     &              vol_strain(md)
c     
            end  do
            
c     
c     ****   printout average stress information   ****
c     
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) ' ' 
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) 
     &           'Stresses (Convention:Compression Positive)' 
            if(ntty.eq.2 .and. iout .ne. 0) then
               if(flag_principal.eq.1) then
                  write(iout,8119) 
               elseif(flag_excess_shear.eq.1) then
                  write(iout,8118) 
               else                 
                  write(iout,7119)
               endif
            endif
            if(iptty.gt.0) write(iptty,*) ' ' 
            if(iptty.gt.0) write(iptty,*) 
     &           'Stresses (Convention:Compression Positive)' 
            if(flag_excess_shear.eq.1) then
               do ishear=1,n0
                  stress_factor_initial = 0.
                  if(flag_permmodel.eq.1) then
                     iispmd = ispm(ishear)    
                     ispmd = ispmt(iispmd) 
                     if(ispmd.eq.22.or.ispmd.eq.23) then
                        friction = spm1f(iispmd)
                        strength = spm2f(iispmd)
                        pp_fac=spm3f(iispmd)
                        if(incremental_shear_permmodel.eq.1) then
                           stress_factor_initial = spm12f(iispmd)
                        endif
                     else
                        friction = friction_out
                        strength = strength_out
                        pp_fac=pp_fac_out
                     endif
                  else
                     friction = friction_out
                     strength = strength_out
                     pp_fac=pp_fac_out
                  endif
                  call principal_incremental_stress_3D(ishear,alambda,
     &                 eigenvec,stress_factor_initial)
c                  call principal_stress_3D(ishear,alambda,eigenvec)
                  call max_excess_shear(ishear,friction,strength,
     &                 alambda,pp_fac,stress_factor_initial,eigenvec)
                  shear_angle(md) = shear_angle(md)*180./pi
               enddo
            endif
            if(iptty.gt.0) then
               if(flag_principal.eq.1) then
                  write(iptty,8119) 
               elseif(flag_excess_shear.eq.1) then
                  write(iptty,8118)
               else                 
                  write(iptty,7119)
               endif
            endif
            do   i=1,m
               md     =  nskw(i)               
               if(flag_principal.eq.1) then
c     s kelkar nov 5,2010 output plane of max excess shear
                  call principal_stress_3D(md,alambda,eigenvec)
                  eigenvec_deg(1)=dacos(eigenvec(3,3))*180./pi   
                  eigenvec_deg(2)=dacos(eigenvec(1,3))*180./pi   
                  eigenvec_deg(3)=dacos(eigenvec(1,1))*180./pi   
                  if(ntty.eq.2 .and. iout .ne. 0) write(iout,7120)  
     &                 md, alambda(3),alambda(2),alambda(1)
     &                 ,eigenvec_deg(1),eigenvec_deg(2),eigenvec_deg(3)
                  if(iptty.gt.0) write(iptty,7120) 
     &                 md, alambda(3),alambda(2),alambda(1)
     &                 ,eigenvec_deg(1),eigenvec_deg(2),eigenvec_deg(3)
c               elseif(flag_excess_shear.eq.1) then
c     s kelkar nov 5,2010 output plane of max excess shear
c                  call max_excess_shear(md,friction,strength,alambda,
c     &                 pp_fac,stress_factor_initial,eigenvec)
c                  shear_angle = shear_angle*180./pi
                  if(ntty.eq.2 .and. iout .ne. 0) write(iout,8120)  
     &                 md, str_x(md), str_y(md), str_z(md), 
     &                 str_xy(md), str_xz(md), str_yz(md),ps(md)
     &                 ,excess_shear(md),shear_angle(md)
                  if(iptty.gt.0) write(iptty,8120) 
     &                 md, str_x(md), str_y(md), str_z(md), 
     &                 str_xy(md), str_xz(md), str_yz(md),ps(md)
     &                 ,excess_shear(md),shear_angle(md)
                  
               else
                  if(ntty.eq.2 .and. iout .ne. 0) write(iout,7120)  
     &                 md, str_x(md), str_y(md), 0.0, 
     &                 str_xy(md), 0.0, 0.0, ps(md)
                  if(iptty.gt.0) write(iptty,7120) 
     &                 md, str_x(md), str_y(md), 0.0, 
     &                 str_xy(md), 0.0, 0.0, ps(md)
               endif
            enddo
            
c     add total volume change
            vol_tot_change = 0.0
            do i = 1, neq
               vol_tot_change = vol_tot_change + vol_strain(i)*sx1(i)
            enddo
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*)  
            if(iptty.gt.0) write(iptty,*)
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,7121)  
     &           vol_tot_change, vtot
            if(iptty.gt.0) write(iptty,7121) vol_tot_change, vtot
         else
c     
c     ****   printout displacement/force information   ****
c     
            ibp_stress=max(ibp_stress,1)
            
            if(tol_stress.ne.0.0d0) then  
               if(mlz.eq.-1) then
                  if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) 
     &                 '>>> stress solution did not reach tolerance <<<'
               endif
               if(ntty.eq.2 .and. iout .ne. 0) write(iout,6118) 
     &              ibp_stress,bp(ibp_stress+nrhs(1)),
     &              bp(ibp_stress+nrhs(2)), bp(ibp_stress+nrhs(3))
            endif 
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) ' ' 
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) 
     &           'Displacements, Forces, Volume Strain ' 
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,6119) 
            
            
            if(tol_stress.ne.0.0d0) then 
               if(mlz.eq.-1) then
                  if(iptty.gt.0) write(iptty,*) 
     &                 '>>> stress solution did not reach tolerance <<<'
               endif
               if(iptty.gt.0) write(iptty,6118) ibp_stress,
     &              bp(ibp_stress+nrhs(1)),bp(ibp_stress+nrhs(2)),
     &              bp(ibp_stress+nrhs(3))
            endif
            if(iptty.gt.0) write(iptty,*) ' ' 
            if(iptty.gt.0) write(iptty,*) 
     &           'Displacements, Forces, Volume Strain ' 
            if(iptty.gt.0) write(iptty,6119) 
            
c     
            do   i=1,m
               md     =  nskw(i)               
               if(ntty.eq.2 .and. iout .ne. 0) write(iout,6120)  
     &              md,du(md),dv(md),dw(md)
     &              ,(forc(md,j),j=1,3),vol_strain(md)
               if(iptty.gt.0) write(iptty,6120) 
     &              md, du(md), dv(md), dw(md)
     &              , (forc(md,j),j=1,3), vol_strain(md)
c     
            end  do
            
c     
c     ****   printout average stress information   ****
c     
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) ' ' 
            if(iptty.gt.0) write(iout,*)
     &           'Stresses(Convention: Compression is Positive (sign))' 
            if(flag_excess_shear.eq.1) then
               if(ntty.eq.2 .and. iout .ne. 0) write(iout,8118)
               do ishear= 1,n0
                  stress_factor_initial = 0.
                  if(flag_permmodel.eq.1) then
                     iispmd = ispm(ishear)    
                     ispmd = ispmt(iispmd) 
                     if(ispmd.eq.22.or.ispmd.eq.23) then
                        friction = spm1f(iispmd)
                        strength = spm2f(iispmd)
                        pp_fac = spm3f(iispmd)
                        if(incremental_shear_permmodel.eq.1) then
                           stress_factor_initial = spm12f(iispmd)
                        endif
                     else
                        friction = friction_out
                        strength = strength_out
                        pp_fac = pp_fac_out
                     endif
                  else
                     friction = friction_out
                     strength = strength_out
                     pp_fac = pp_fac_out
                  endif
                  call principal_incremental_stress_3D(ishear,alambda,
     &                 eigenvec,stress_factor_initial)
c                  call principal_stress_3D(ishear,alambda,eigenvec)
                  call max_excess_shear(ishear,friction,strength,
     &                 alambda,pp_fac,stress_factor_initial,eigenvec)
                  shear_angle(ishear) = shear_angle(ishear)*180./pi
               enddo
            endif                 
            if(ntty.eq.2 .and. iout .ne. 0) then
               if(flag_principal.eq.1) then
                  write(iout,8119) 
               elseif(flag_excess_shear.eq.1) then
                  write(iout,8118)
               else                 
                  write(iout,7119)
               endif
            endif
            if(iptty.gt.0) write(iptty,*) ' ' 
            if(iptty.gt.0) write(iptty,*)
     &           'Stresses(Convention: Compression is Positive (sign))'  
            if(iptty.gt.0) then
               if(flag_principal.eq.1) then
                  write(iptty,8119) 
               elseif(flag_excess_shear.eq.1) then
                  write(iptty,8118) 
               else                 
                  write(iptty,7119)
               endif
            endif
            do   i=1,m
               md     =  nskw(i) 
               
               if(flag_principal.eq.1) then
c     s kelkar nov 5,2010 output plane of max excess shear
                  call principal_stress_3D(md,alambda,eigenvec)
                  eigenvec_deg(1)=dacos(eigenvec(3,3))*180./pi   
                  eigenvec_deg(2)=dacos(eigenvec(1,3))*180./pi   
                  eigenvec_deg(3)=dacos(eigenvec(1,1))*180./pi   
                  if(ntty.eq.2) write(iout,7120)  
     &                 md, alambda(3),alambda(2),alambda(1)
     &                 ,eigenvec_deg(1),eigenvec_deg(2),eigenvec_deg(3)
                  if(iptty.gt.0) write(iptty,7120) 
     &                 md, alambda(3),alambda(2),alambda(1)
     &                 ,eigenvec(3,3),eigenvec(1,3),eigenvec(1,1)
               elseif(flag_excess_shear.eq.1) then
c     s kelkar nov 5,2010 output plane of max excess shear
                  if(flag_permmodel.eq.1) then
                     stress_factor_initial = 0.
                     iispmd = ispm(md)    
                     ispmd = ispmt(iispmd) 
                     if(ispmd.eq.22.or.ispmd.eq.23) then
                        friction = spm1f(iispmd)
                        strength = spm2f(iispmd)
                        pp_fac=spm3f(iispmd)
                        if(incremental_shear_permmodel.eq.1) then
                           stress_factor_initial = spm12f(iispmd)
                        endif
                     else
                        friction = friction_out
                        strength = strength_out
                        pp_fac=pp_fac_out
                     endif
                  else
                     friction = friction_out
                     strength = strength_out
                     pp_fac=pp_fac_out
                  endif
                  call principal_incremental_stress_3D(md,alambda,
     &                 eigenvec,stress_factor_initial)
c                  call principal_stress_3D(md,alambda,eigenvec)
                  call max_excess_shear(md,friction,strength,alambda,
     &                 pp_fac,stress_factor_initial,eigenvec)
                  shear_angle = shear_angle*180./pi
                  if(ntty.eq.2) write(iout,8120)  
     &                 md, str_x(md), str_y(md), str_z(md), 
     &                 str_xy(md), str_xz(md), str_yz(md),ps(md)
     &                 ,excess_shear(md),shear_angle(md)
                  if(iptty.gt.0) write(iptty,8120) 
     &                 md, str_x(md), str_y(md), str_z(md), 
     &                 str_xy(md), str_xz(md), str_yz(md),ps(md)
     &                 ,excess_shear(md),shear_angle(md)
                  
               else
                  if(ntty.eq.2) write(iout,7120)  
     &                 md, str_x(md), str_y(md), str_z(md), 
     &                 str_xy(md), str_xz(md), str_yz(md),ps(md)
                  if(iptty.gt.0) write(iptty,7120) 
     &                 md, str_x(md), str_y(md), str_z(md), 
     &                 str_xy(md), str_xz(md), str_yz(md),ps(md)
c     
               endif
            end  do


                                !!!! SAI : PRINTOUT FOR STRAINS !!!!

            if(ifem.eq.0) call compute_strains()

            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*) ' '
            if(iptty.gt.0) write(iout,*)
     &           'Strains'
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,7119)
            if(iptty.gt.0) write(iptty,*) ' '
            if(iptty.gt.0) write(iptty,*)
     &           'Strains'
            if(iptty.gt.0) write(iptty,7119)
            do   i=1,m
               md     =  nskw(i)
               if(ntty.eq.2 .and. iout .ne. 0) write(iout,7120)
     &              md, strain_xx(md), strain_yy(md), strain_zz(md),
     &              strain_xy(md), strain_zx(md), strain_yz(md),ps(md)
               if(iptty.gt.0) write(iptty,7120)
     &              md, strain_xx(md), strain_yy(md), strain_zz(md),
     &              strain_xy(md), strain_zx(md), strain_yz(md),ps(md)
c     
            end  do

            if(flag_strain_output.eq.1) then
               write(ifile_strain,7191) l,days
 7191          format(1x,'Time Step',i10,'Days',g16.9)
               if(flag_strainout_first.eq.0) then
                  flag_strainout_first = 1 
                  do md=1,neq
                     write(ifile_strain,7192)
     &                    cord(md,1), cord(md,2), cord(md,3), md,
     &                    strain_xx(md), strain_yy(md), strain_zz(md),
     &                    strain_xy(md), strain_zx(md), strain_yz(md)
                  enddo
               else
                  do md=1,neq
                     write(ifile_strain,7193) md,
     &                    strain_xx(md), strain_yy(md), strain_zz(md),
     &                    strain_xy(md), strain_zx(md), strain_yz(md)
                  enddo
               endif
            endif

 6116       format(/,
     &           ' ***************************************************',
     &           '***************************',/)
 6117       format(1x,' OUTPUT FOR STRESS SOLUTION')

 773        format(/,20x,'Equation Performance')
 75         format(1x,'Number of N-R Iterations: ',1i10)
 76         format(1x,'Avg # of Linear Equation Solver Iterations: ',
     &           1f5.1)
 77         format(1x,'Number of Active Nodes: ',1f10.0)
 704        format(1x,'Total Number of Iterations, N-R: ',i10,
     &           ' , Solver: ',i10)
 6118       format(1x,' Max residual Eqs.: ',
     &           'node ',i7,1p,' x  ',g12.5,' y  ',g12.5,
     &           ' z  ',g12.5) 
 6119       format(3x,'Node',3x,'u disp',7x,'v disp',7x,'w disp',
     &           4x,' force x   force y   force z  vol strain')
 6120       format(i7,1p,3(2x,g11.4),4(1x,g9.2))
 7119       format(3x,'Node',4x,'x str ',7x,'y str ',7x,'z str ',
     &           5x,' xy str    xz str    yz str  ',' porosity')
 8118       format(3x,'Node',4x,'x str ',7x,'y str ',7x,'z str ',
     &           5x,' xy str    xz str    yz str  ',' porosity',
     &           '  excess_shear  ', '  shear_angle  ')
 8119       format(3x,'Node',4x,'sigmax ',7x,'sig-2 ',7x,'sig-min ',
     &           5x,' angle-sigmax_z    sigmax-x    sigmin-x  ')
 7120       format(i7,1p,3(1x,g12.5),3(1x,g9.2),1x,0p,f9.5)
 7192       format(3(1x,g12.5),i7,1p,9(1x,g12.5))
 7193       format(i7,1p,9(1x,g12.5))
 8120       format(i7,1p,3(1x,g12.5),3(1x,g9.2),1x,0p,f9.5,2(1x,g9.2))
c     add total volume change
            vol_tot_change = 0.0
            do i = 1, neq
               vol_tot_change = vol_tot_change + vol_strain(i)*sx1(i)
	    enddo
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,*)  
            if(iptty.gt.0) write(iptty,*)
            if(ntty.eq.2 .and. iout .ne. 0) write(iout,7121)  
     &           vol_tot_change, vtot
            if(iptty.gt.0) write(iptty,7121) vol_tot_change, vtot
 7121       format (1x,'*** Total Volume Change = ', 1p,g11.4,
     &           1x,'Total Volume ',g11.4,1x,'***')

         endif

      else if(iflg.eq.12) then
         if(icnl.eq.0) then
            do i = 1,neq
               duo(i) = du(i)
               dvo(i) = dv(i)
               dwo(i) = dw(i)
            enddo  
         else 
            do i = 1,neq
               duo(i) = du(i)
               dvo(i) = dv(i)
            enddo           	            
         endif
c     vol_strain0 = vol_strain
      else if(iflg.eq.13) then
c     
c     calculate stresses (now improved at the cost of more storage 
c     gaz - 111106)
c     
         if(icnl.eq.0) then 
c     3D  
            if(ifem.eq.1) then
               call stress_3D_post_fem()
            else
               do i = 1,neq
                  call stress_3D_post(i)
               enddo
            endif

         else
c     2D
            if(ifem.eq.1) then
               call stress_2D_post_fem()
            else
               do i = 1,neq
                  call stress_2D_post(i)
               enddo
            endif
         end if
      else if(iflg.eq.14) then
c     
c     write special history plot for stress 
c     create one file for each time series
c     
         if(istrshis.eq.0) return
         if(l.eq.0) then
            isstr_temp = 80
            icount = 1000
            do i = 1, istrshis 
               isstr_temp = isstr_temp + 1
               icount = icount +1
               write(hist_str(15:18),'(i4)') icount
               hist_str(19:22) = '.txt'
               open(isstr_temp, file = hist_str, status='unknown',
     &              form = 'formatted')
               
               jj = nskw_stress(i,1)
               xi = cord(jj,1)
               yi = cord(jj,2)
               if(icnl.eq.0) then
                  zi = cord(jj,3)
               else
                  zi = 0.0
               endif
               if(nskw_stress(i,2).eq. 1)then
                  ptype='disx '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype
               else if(nskw_stress(i,2).eq. 2)then
                  ptype='disy '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype  
               else if(nskw_stress(i,2).eq. 3)then
                  ptype='disz '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype 	  
               else if(nskw_stress(i,2).eq. 4)then
                  ptype='strx '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype  	  
               else if(nskw_stress(i,2).eq. 5)then
                  ptype='stry '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype
               else if(nskw_stress(i,2).eq. 6)then
                  ptype='strz '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype  
               else if(nskw_stress(i,2).eq. 7)then
                  ptype='strxy '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype  
               else if(nskw_stress(i,2).eq. 8)then
                  ptype='strxz '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype       
               else if(nskw_stress(i,2).eq. 9)then
                  ptype='stryz '
                  write(isstr_temp,1200) jj,xi,yi,zi,ptype      
               endif          
 1200          format(1x,'node = ',i8,' xyz = ',1p,3g14.5,' type = ',a5)
            enddo
            return
         endif
         isstr_temp = 80
         do i = 1, istrshis 
            isstr_temp = isstr_temp + 1
            jj =nskw_stress(i,1)
            if(nskw_stress(i,2).eq.1)then
               write(isstr_temp,*) days, du(jj)
            else if(nskw_stress(i,2).eq.2)then
               write(isstr_temp,*) days, dv(jj)
            else if(nskw_stress(i,2).eq.3.and.icnl.eq.0)then
               write(isstr_temp,*) days, dw(jj)  
            else if(nskw_stress(i,2).eq.4)then
               write(isstr_temp,*) days, str_x(jj)  
            else if(nskw_stress(i,2).eq.5)then
               write(isstr_temp,*) days, str_y(jj)  
            else if(nskw_stress(i,2).eq.6.and.icnl.eq.0)then
               write(isstr_temp,*) days, str_z(jj)  
            else if(nskw_stress(i,2).eq.7)then
               write(isstr_temp,*) days, str_xy(jj)       
            else if(nskw_stress(i,2).eq.8.and.icnl.eq.0)then
               write(isstr_temp,*) days, str_xz(jj)   
            else if(nskw_stress(i,2).eq.9.and.icnl.eq.0)then
               write(isstr_temp,*) days, str_yz(jj)   	      
            endif               
         enddo
      else if(iflg.eq.16) then
c     manage perm changes with displacement
         call stress_perm(0,0)
      else if(iflg.eq.17) then 
c     save residual information for flow
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
         if(.not.allocated(bp_flow1)) then
            if(idualp.ne.0) then
               ilev=3
               mlev=m/3
            else if(idpdp.ne.0) then
               ilev=2
               mlev=m/2
            else
               ilev=1
               mlev=m
            endif
            k = 0
            do il=1,ilev
               do i=1,mlev 
                  k = k + 1
               enddo
            enddo  
            allocate (bp_flow1(k))
            allocate (bp_flow2(k))         
         endif
         if(idualp.ne.0) then
            ilev=3
            mlev=m/3
         else if(idpdp.ne.0) then
            ilev=2
            mlev=m/2
         else
            ilev=1
            mlev=m
         endif
         k = 0
         do il=1,ilev
            do i=1,mlev
               md=nskw(i+(il-1)*mlev)
               k = k+1
               bp_flow1(k) = bp(md)
               bp_flow2(k) = bp(md+neq)
            enddo
         enddo  
      else if(iflg.eq.18) then 
c     repopulate residual information for flow
         if(idualp.ne.0) then
            ilev=3
            mlev=m/3
         else if(idpdp.ne.0) then
            ilev=2
            mlev=m/2
         else
            ilev=1
            mlev=m
         endif
         k = 0
         do il=1,ilev
            do i=1,mlev
               md=nskw(i+(il-1)*mlev)
               k = k+1
               bp(md)= bp_flow1(k)
               bp(md+neq)= bp_flow2(k)
            enddo
         enddo   
      else if(iflg.eq.19) then  
c     calculate subsidence
         call subsidence(1)              
      else if(iflg.ge.100) then
c     user subroutines. not enabled
         
      endif
      

      return
      end
      subroutine geom_stress(iflg) 
c     
c     routine to calculate area, forces, etc. for stress bcs
c     
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comji
      use comki
      use comxi
      use davidi
      use comsi  
      
      implicit none
      
      integer iflg,izone,idir,i,i1,i2,jj,kb
      real*8 x1,x2,y1,y2,areai
      real*8 disi,dist1,dist2,gridblock_dem
      
      if(iflg.eq.0) then
      else if(iflg.eq.1) then
c     
c     calculate the area in a direction for a zone
c     for use with stress BC
c     
         area_str_zone = 0.0
         do i = 1, n0
c     check x direction       
            idir = iarea_str(i,1)
            if(idir.ne.0) then
               disi = cord(i,idir)
               dist1 = disi
               dist2 = disi
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do jj = i1,i2
                  kb = nelm(jj)
                  dist1 = max(cord(kb,idir),dist1)
                  dist2 = min(cord(kb,idir),dist2)
               enddo
               if(dist1-disi.ge.disi-dist2) then
                  gridblock_dem = 0.5*(dist1-dist2)
               else
                  gridblock_dem = -0.5*(dist1-dist2)
               endif
               areai = sx1(i)/gridblock_dem
               area_str(i,idir) = areai
            endif
c     check y direction       
            idir = iarea_str(i,2)
            if(idir.ne.0) then
               disi = cord(i,idir)
               dist1 = disi
               dist2 = disi
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do jj = i1,i2
                  kb = nelm(jj)
                  dist1 = max(cord(kb,idir),dist1)
                  dist2 = min(cord(kb,idir),dist2)
               enddo
               if(dist1-disi.ge.disi-dist2) then
                  gridblock_dem = 0.5*(dist1-dist2)
               else
                  gridblock_dem = -0.5*(dist1-dist2)
               endif
               areai = sx1(i)/gridblock_dem
               area_str(i,idir) = areai
            endif
c     check z direction       
            idir = iarea_str(i,3)
            if(idir.ne.0) then
               disi = cord(i,idir)
               dist1 = disi
               dist2 = disi
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do jj = i1,i2
                  kb = nelm(jj)
                  dist1 = max(cord(kb,idir),dist1)
                  dist2 = min(cord(kb,idir),dist2)
               enddo
               if(dist1-disi.ge.disi-dist2) then
                  gridblock_dem = 0.5*(dist1-dist2)
               else
                  gridblock_dem = -0.5*(dist1-dist2)
               endif
               areai = sx1(i)/gridblock_dem
               area_str(i,idir) = areai
            endif                
         enddo    
      else if(iflg.eq.2) then
      endif
      return
      end
      


      subroutine wellphysicsctr(iflg,ndummy)
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
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
!D1  PURPOSE
!D1
!D1      Well physics including the drift flux model
!D1      
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 10-20-08, Programmer: G. Zyvoloski
!D2
!D2 $Log:
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.2 Finite-Element Coefficient Generation
!D3
C***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 GAZ initial implementation 1-17-05
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
C***********************************************************************
C  
C  
C Notes:
C
C
C
C***********************************************************************
C
C     INPUT ARGUMENTS -
C
C
C#######################################################################
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comki
      use comrxni
      use comii
      use comflow
      use comxi
      use davidi
      use comwellphys
      use comco2
      use commeth
c     
      implicit none
      integer ndummy, iz, iwell_type, mi, ieosd, it, iwelld
      integer iflg,i,j,i1,i2,ik,izone,inode,ii,ij,kb,iafile,neqp1_old
      integer neqp1,iii,node,neq_new, id, nbc, kk
      integer imr,ipr,imrl,iprl,kb2,j1,j2,i3, ncont
      integer neqpi, iexpfl, icoef,nr_old,kbmin,kbmax
      integer open_file, inwel,max_seg_div,max_con,max_well_con
      integer maxiarea,iareap, jj, n_ncon, irnode, isnode, icesd

      real*8  vol0, vol2, pi, pi2, disa, raddif,disr,rad2,amult,area_tol
      real*8  sx_max,area_dum,area_new,sx_dum
      real*8  seg,seg2,dis,segx,segy,segz,segx2,segy2,segz2
      real*8  disx,disy,disz
      real*8  vg,vl,vvkb,sg,sl,covg,vvol,vsgf,sur_ten
      real*8  aparam,bparam,gamma,beta,fv,vdrift,vdrift1,vdrift2
      real*8  ar1,ar2,vv_max,vl_max,vb_max,akuta,vcar
      real*8  well1, well2 , well3, well4, well5, well6, well_rad  
      real*8 term1a,dterm1a,term1,dterm1,rlgas,drlgass
c      real*8 rov,drovt,drovp,dvisvt,dvisvp
c      real*8 rol,drolt,drolp,dvislt,dvislp  
      real*8 vismix,dvismixp,dvismixt,dvismixs  
c      real*8 sw,visv,visl,divipl,tref,tempc
c gaz 110819 removed tref, pref (now global)      
      real*8 sw,divipl,tempc
      real*8 xvisl0,xvisl,drocp0
      real*8 pl,xvisv,agrav

      real*8, allocatable  :: dum(:)
      real*8, allocatable  :: dum1(:,:)
      real*8, allocatable  :: seg_min(:)

      integer, allocatable :: idum(:)
      integer, allocatable :: idum1(:,:)
      integer, allocatable :: idum2(:)
      integer, allocatable :: ncon_new(:)
      
c     
      parameter (area_tol=1.e-10, pi= 3.1415927, pi2=6.2831850)
      parameter (max_seg_div = 100,max_con = 10,max_well_con= 4)
c      parameter(pcl0 = 0.101325)
c      parameter(roc0 = 1.292864)
c     

      logical null1

      character*80 dum_string

      integer i4,ic1,node_seg,iseg,ntot2,nodew
      if (.not. allocated(seg_min)) allocate(seg_min(max_con))
      

C     
C#######################################################################
C     
c     
c     return if no (non-darcy) well physics  
c   
      pcl0 = 0.101325 
c gaz debug roc0 now in comai
c      roc0 = 1.292864
      if(nwellphy.eq.0) return
      agrav = abs(grav)
      if(iflg.eq.0) then
c     
c     any required initial calcs
c     
c     allocate arrays

         allocate(iwellpt(neq))
         allocate(iwellp(neq))
         
         allocate(rwell1f(neq))
         allocate(rwell2f(neq))
         allocate(rwell3f(neq))
         allocate(rwell4f(neq))
         allocate(rwell5f(neq))
         allocate(rwell6f(neq))

         allocate(rolsvf(neq))
         allocate(rovsvf(neq)) 
         
         allocate(mdriftf(neq))      
         allocate(dmdriftp(neq))
         allocate(dmdrifte(neq))
         allocate(dmdriftw(neq))                   

         allocate(vel_m_well2(neq))
         allocate(rolsv(neq))
         allocate(rovsv(neq))
c     
c     read in data (like rlperm.f)
c     
c     check for read from other file
c     
         iwellp = 0
         i = 0
         j = 0
         ex = .false.
 10      continue
         read(inpt,'(a80)') wdd1
         if(.not. null1(wdd1)) then
            backspace inpt
            read(inpt,*) iwell_type
            backspace inpt
            i = i+1
            
            if (iwell_type .eq. 1.or.iwell_type .eq. 2) then
c     simplified drift flux model of vapor flow
c     for use with CO2 
c     rwell1f is constant term ar1
c     rwell2f is constant term ar2
c     rwell3f is linear term aparam
c     rwell4f is linear term bparam
c     rwell5f is well radius    
               read(inpt,*) iwellpt(i),rwell1f(i),rwell2f(i),rwell3f(i),
     &              rwell4f(i), rwell5f(i), rwell6f(i)
            endif
         else
            go to 20
         endif   
         go to 10   
 20      continue
c     count the number of models
         nwellphy = i
         narrays = 1
         itype(1) = 4
         default(1) = 1
         macro = "wellphysics "
         igroup = 2
         call initdata2( inpt, ischk, n0, narrays,
     2        itype, default, macroread(7), macro, igroup, ireturn,
     3        i4_1=iwellp(1:n0) )

         macroread(7) = .TRUE.
         

      else if(iflg.eq.1) then
c     
c     calculate drift flux flow parameters and derivatives
c     
         do i = 1,neq
            mi = i+ndummy
            
            if (iwellp_flag .eq. 0) then
               it = 0
            else
               it = iwellp(mi)
            end if
            if(it.eq.0) then
               iwelld=0
            else
               iwelld = iwellpt(it)
            endif

            if (iwelld .ne. 0) then
               well1 = rwell1f(it)
               well2 = rwell2f(it)
               well3 = rwell3f(it)
               well4 = rwell4f(it)
               well5 = rwell5f(it) 
               well6 = rwell6f(it)                  
            else
!     Use linear model as default
               well1 = 0.
               well2 = 0.
               well3 = 1.
               well4 = 1.
               well5 = 1.
               well6 = 1.
               iwelld = 1
            end if
            
            
            if(iwelld.eq.1) then

c     CO2 model
c     calculate parameters for drift-flux model
c     vg = covg*vvol+vdrift
c     vg = gas phase velocity
c     vl = liquid phase velocity
c     covg = paramater in drft flux equation (calculated)
c     vvol = average velocity
c     vdrift = gas drift velocity
c     sur_ten = surface tension  (kg/m**2)
c     sl = liquid saturation
c     sg = gas saturation
c     
c     fixed parameters in drift flux model (for calculating covg)
c     
               icesd = ices(mi)       
               sw = fw(mi)
               sl = fl(mi)
               sg = fg(mi)
               if(icesd.eq.3) sl = fg(mi)         
               ar1 = well1
               ar2 = well2
               aparam = well3
               bparam = well4
               fv = well5 
               well_rad = well6            
c     for now overwrite with fixed values      
c     ar1 = 0.2
c     ar2 = 0.4
c     aparam = 1.2
c     bparam = 0.3
c     fv = 1.0
c     fix surface tension      
               sur_ten = 0.7
c     
c     sg is the CO2 gas fraction         
c     for now set vsgf = vvol (see discussion after eq 5, Shi et al,2005)   
c              vvol = vel_m_well2(ii)   
               vvol = 0.0
               vsgf = 1.0
               beta = max(sg, fv*sg*abs(vvol)/vsgf)        
               gamma = (beta - bparam)/(1.0 - bparam)
               covg = aparam/(1.0 + (aparam - 1)*gamma**2)      
c     
c     determine region where gas saturation lies
c     for now all set to the same expression see discusion after eq 15, Shi et al,2005)
               if(sg.le.ar1) then
                  akuta = 1.53/covg
               else if(sg.ge.ar2) then
                  akuta = 1.53/covg
               else
                  akuta = 1.53/covg 
               endif 
c     calculate the characteristic velocity
c     formulation below evaluates drift velocity at last time step
c     = explicit update
               if(iad.eq.0) then
                  rol = co2_prop(mi)
                  rov = co2_prop(9*neq+mi)
                  rolsv(mi) = rol
                  rovsv(mi) = rov
               else
                  rol = rolsv(mi)
                  rov = rovsv(mi)         
               endif
c     set surface tension (water-air) to 0.7 kg/m**2  
c     need a table - now hardwired     
               vcar = (sur_ten*grav*(rol-rov)/rol**2)**0.25        
c     now set drift velocity
               term1 = (1.0-sg*covg)
               vdrift1 = term1*covg*akuta*vcar
               vdrift2 = sg*covg*sqrt(rov/rol)+term1
               vdrift = vdrift1/vdrift2    
c     
               
c     
c     apparent gas phase rel perm
c     assumption (set linear rel perm in rlperm_co2,no cap pressure!)
c     gas velocity is a function of volumetric flux 
c     vgas is now related to vl
c     vgas = covg*(1.0-sg)/1.0-covg*sg)*vl + vdrift/(1.0-covg*sg)
c     leading sg is for the mass fraction
c     term1 should go to 1 as sg goes to 1 
               
               term1a =  1./(1.0-covg*sg)
               dterm1a = -covg/(1.0-covg*sg)**2  
               term1 = covg*(1.0-sg)*term1a  
               dterm1 = covg*(1.0-sg)*dterm1a - covg*term1a
               rlgas = sg*term1  
               drlgass = term1 + sg*dterm1        
c     
               rov = co2_prop(9*neq+mi)
               drovt=co2_prop(10*neq+mi)
               drovp=co2_prop(11*neq+mi)
               visv=co2_prop(15*neq+mi)
               dvisvt=co2_prop(16*neq+mi)
               dvisvp=co2_prop(17*neq+mi)                 
               rol = co2_prop(mi)
               drolt=co2_prop(neq+mi)
               drolp=co2_prop(2*neq+mi)
               visl=co2_prop(6*neq+mi)
               dvislt=co2_prop(7*neq+mi)
               dvislp=co2_prop(8*neq+mi)    
               vismix = sg*visv + (1.0-sg)*visl
               dvismixp = sg*dvisvp + (1.0-sg)*dvislp 
               dvismixt = sg*dvisvt + (1.0-sg)*dvislt 
               dvismixs = visv-visl                
c     this definition below gives mass flux of gas with velocity of liquid          
               div(mi) = rov/vismix*rlgas
c     
               divp(mi) = (drovp/vismix-rov/vismix**2*dvismixp)*rlgas
               dive(mi) = (drovt/vismix-rov/vismix**2*dvismixt)*rlgas
               divw(mi) = (-rov/vismix**2*dvismixs)*rlgas +
     &              rov/vismix*drlgass
c     
c     mass associated with drift velocity (for now vdrift is explicitly updated so no derivatives)
c     
               mdriftf(mi) = rov*vdrift*term1a 
               dmdriftp(mi) = drovp*vdrift*term1a
               dmdrifte(mi) = drovt*vdrift*term1a     
               dmdriftw(mi) = rov*vdrift*dterm1a                      
c     
            else if(iwelld.eq.2) then

c     aiw-water model (isothermal)
c     calculate parameters for drift-flux model
c     vg = covg*vvol+vdrift
c     vg = gas phase velocity
c     vl = liquid phase velocity
c     covg = paramater in drft flux equation (calculated)
c     vvol = average velocity
c     vdrift = gas drift velocity
c     sur_ten = surface tension  (kg/m**2)
c     sl = liquid saturation
c     sg = gas saturation
c     
c     fixed parameters in drift flux model (for calculating covg)
c     
               pl = phi(mi)
               sg = 1.0-s(mi)
               ar1 = well1
               ar2 = well2
               aparam = well3
               bparam = well4
               ar1 = well1
               ar2 = well2
               aparam = well3
               bparam = well4
               fv = well5 
               well_rad = well6            
c     for now overwrite with fixed values      
c     ar1 = 0.2
c     ar2 = 0.4
c     aparam = 1.2
c     bparam = 0.3 
c     fv = 0.3
c     fix surface tension      
               sur_ten = 0.7
c     
c     sg is the CO2 gas fraction         
c     for now set vsgf = vvol (see discussion after eq 5, Shi et al,2005)   
c     vvol = vel_m_well2(ii)   
               vvol = 0.0
               vsgf = 1.0
               beta = max(sg, fv*sg*abs(vvol)/vsgf)        

c     
c     constrain gamma to 0 <= gamma <= 1.0 
c     
               if(beta.lt.bparam.and.beta.gt.bparam) then
                  gamma = (beta - bparam)/(1.0 - bparam) 
               else if(beta.ge.bparam) then
                  gamma = 1.0
               else if(beta.le.bparam) then 
                  gamma = 0.0           
               endif    
               covg = aparam/(1.0 + (aparam - 1)*gamma**2)      
c     
c     determine region where gas saturation lies
c     for now all set to the same expression see discusion after eq 15, Shi et al,2005)
               if(sg.le.ar1) then
                  akuta = 1.53/covg
               else if(sg.ge.ar2) then
                  akuta = 1.53/covg
               else
                  akuta = 1.53/covg 
                  
               endif 
c     calculate the characteristic velocity
c     formulation below evaluates drift velocity at last time step
c     = explicit update
               if(iad.eq.0) then
                  rol = rolf(mi)
                  rov = rovf(mi)
                  rolsv(mi) = rol
                  rovsv(mi) = rov
               else
                  rol = rolsv(mi)
                  rov = rovsv(mi)         
               endif
c     set surface tension (water-air) to 0.7 kg/m**2  
c     need a table - now hardwired     
               vcar = (sur_ten*agrav*(rol-rov)/rol**2)**0.25        
c     now set drift velocity
               term1 = (1.0-sg*covg)
               vdrift1 = term1*covg*akuta*vcar
               vdrift2 = sg*covg*sqrt(rov/rol)+term1
               vdrift = vdrift1/vdrift2    
c     
c     
c     apparent gas phase rel perm
c     assumption (set linear rel perm in rlperm_co2,no cap pressure!)
c     gas velocity is a function of volumetric flux 
c     vgas is now related to vl
c     vgas = covg*(1.0-sg)/1.0-covg*sg)*vl + vdrift/(1.0-covg*sg)
c     leading sg is for the mass fraction
c     term1 should go to 1 as sg goes to 1 
c     remember  dsg = -ds (sg = 1-s(mi)) 
c     
               term1a =  1./(1.0-covg*sg)
               dterm1a = -covg/(1.0-covg*sg)**2  
               term1 = covg*(1.0-sg)*term1a  
               dterm1 = covg*(1.0-sg)*dterm1a - covg*term1a
               rlgas = sg*term1  
               drlgass = term1 + sg*dterm1        
c     
c     misc. constants  (like at the top of thrair.f)
c  
c gaz 110819 pref, tref (global) read in scanin                 
c               tref = crl(6,1)
               tempc=(273.0)/(tref+273.0)
               drocp0=roc0*tempc/pcl0
               xvisl0=crl(2,1)
               xvisl = xvisl0          
c               pref=crl(4,1)
               xvisv=crl(5,1)
               rov = drocp0*pl
               drovp = drocp0
               vismix = sg*xvisv + (1.0-sg)*xvisl
               dvismixp = 0.0
               dvismixt = 0.0
               dvismixs = xvisv-xvisl                
c     this definition below gives mass flux of gas with velocity of liquid          
               div(mi) = rov/vismix*rlgas
c     
               divp(mi) = (drovp/vismix-rov/vismix**2*dvismixp)*rlgas
c     dive has derivative wrt saturation (-1 from sg = 1.-sl)
               dive(mi) = ((-rov/vismix**2*dvismixs)*rlgas +
     &              rov/vismix*drlgass)*(-1.)
c     
c     mass associated with drift velocity (for now vdrift is explicitly updated so no derivatives)
c     
               mdriftf(mi) = rov*vdrift*term1a 
               dmdriftp(mi) = drovp*vdrift*term1a
c     dmdrifte has derivative wrt saturation (-1 from sg = 1.-sl) 
               dmdrifte(mi) = (rov*vdrift*dterm1a)*(-1.)                   
c     
               mdriftf(mi) = 0.0
               dmdriftp(mi) = 0.0
               dmdrifte(mi) = 0.0           
            endif
         enddo
         continue
      else if(iflg.eq.3) then
c     
c     add BC for wellphysics
c     
      else if(iflg.eq.4) then
c     
c     apply N-R correction
c     
      else if(iflg.eq.6) then	
c     
c     organize output 
c     
      endif
      end

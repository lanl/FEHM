subroutine rlp_cap(ndummy)
!*************************************************************************
! Copyright  2015.   Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos
! National  Laboratory  (LANL),  which is operated by  Los Alamos National
! Security, LLC  for the U. S. Department of Energy.  The U. S. Government
! has rights to use, reproduce, and distribute this software.  Neither the
! U. S. Government nor Los Alamos National Security, LLC or persons acting
! on their behalf,  make any warranty,  express or implied, or assumes any
! liability for the accuracy, completeness, or usefulness of the software,
! any information pertaining to the software,  or  represents that its use
! would not infringe privately owned rights.

! The  software  being licensed  may  be Export Controlled.  It may not be
! distributed  or  used by individuals  or entities prohibited from having
! access to the software package, pursuant to United States export control
! laws and regulations. An export control review and determination must be
! completed before LANS will provide access to the identified Software.
!*************************************************************************

  use comai
  use combi
  use comcomp, only: ioil, pcpow, pcpog, rl_g
  use comci
  use comco2
  use comdi
  use comdti
  use comei
  use comfi
  use comhi
  use comki
  use comrlp
  implicit none
  real*8 scutm,hmin,tol_l,tol_u,su_cut
  real*8 rl_assign,ds_assign,drls_assign
  parameter(scutm = 1.d-03)
  parameter(hmin = 1.d-18)
  parameter(tol_l  = 1.d-5)
  parameter(tol_u  = 1.d-5)
  parameter(su_cut = 0.99d00)
  !      parameter(su_cut = 0.70d00)

  integer i, j, k, ndummy, mrlp, it, ir, iphase, lphase, phase(3),ii
  integer ieosd, mi, icouple, ir2, irf, i2, aphase, itbl, iclmn,mm
  integer rphase, cap_select,pair(20,2),ity
  integer :: ireg = 1
  real*8 alpha,beta,alamda,alpi,smcut,slcut,ds,dhp,dhp1
  real*8 rp1,rp2,rp3,rp4,denom,star,hp,hp1,b1,s_star
  real*8 akf,akm,porf,permb,sl,sl1,sl2,rp5,rp6,rp7
  real*8 smcutm,smcutf,alpham,alamdam,facf,sw,sg,snw
  real*8 rl, rv, drls, drvs, rl1, rv1, drls1, drvs1, cp, dpcp
  real*8 drlw, drlg, drlk, drvw, drvg, drvk, drtw, drtg, drtk
  real*8 skip,dum1,dum2,dum3,dummy,rlc,drlc,rt,drt,drl_wk  ,pn
  !      real*8 drllg, rlv,drlvw,drlvg
  real*8 drls2,drvs2,blank(10)
  logical null1,ex
  character model_type*10
  !     initialize everything as if single wetting phase
  if(icarb.ne.0) then
     pcp=0.
     pcg=0.
     dpcgw=0.
     dpcgg=0.
     dpcpw=0
     dpcpg=0 
     rl_w=0.;drl_ww=0.;drl_wg=0.;rl_l=0.;drl_lw=0.;drl_lg=0.
     rl_v=0.;drl_vw=0.;drl_vg=0.
  else if (ioil .ne. 0) then
     pcpow = 0.; pcpog = 0.
     rl_w=1.; rl_l=0.; rl_g = 0.
  else
     if (ndummy .eq. 0) then
        rlf = 0.; rvf = 0.; drlef = 0.; drvef = 0.
        pcp=0.
        dpcef=0.
     end if
  endif
  !      
  ! loop over every node
  do i = 1, neq
     mi = i + ndummy
     
     ieosd = ieos(mi)
     !     set it to the relative permeability/cap pressure group number
     it = irlp(mi)

     ! based on what type of problem (oil, co2, other) determine which phases are present 
     ! assign values to  k   
     k = 100
     ! k=1-8 single phase     
     ! k=25 air-water or k=28 water-vapor 
     ! 10-water/co2 liquid  11-water/co2 gas  12: liquid co2/gas co2 
     ! 13-water/oil   14:water/gas  15:oil/gas
     ! 16:three phase

     if (ioil .eq. 1) then
        call oil_phase(mi,it,k,sw) 
     else if (icarb .eq. 0) then
        call water_air_phase(ieosd,k)            
     else
        call co2_phase(mi,k) 
     endif

     ! define wetting phase saturation (sw) and non-wetting saturation  (snw)   
	  if(k<9) then
	  		sw=s(mi)
      elseif(k.eq.25 .or. k.eq.28) then   ! water/air or water/vapor
           sw = s(mi)
           snw = 1 - sw
           !                  snw = ?????
        elseif(k==22)  then  ! co2 liquid/gas
           sw=fl(mi)
           snw=fg(mi)
        elseif(k==21) then  ! water/co2_gas
           sw=fw(mi)
           snw=fg(mi)
        else   ! co2 or oil/gas with water
           sw= fw(mi)
           snw=fl(mi)
        endif



     !  single phase  (need to finish this)
     !  I think a cap calculation should be done single phase or not so modify conditional
     !if(k<9)	then
!     write(*,*) 'phase is ',k
! 1=air 2=water (not co2 problem), 3=water (co2 problem), 4=co2_gas, 5=vapor, 6=co2_liquid
! 7=oil 8=gas, 9=co2_sc , 11=methane_hydrate
     if(k<9)	then     
     	select case(k)
        ! need air/water case
        ! these two phase cases never seem to use all 4 derivatives
        case(1)  !air
    		if(icarb>0) then
        	rl_w(mi)=0
        	rl_l(mi)=1
        	else
        	rlf(mi)=0.;rvf(mi)=1.
        	endif
        
     	case(2,3)   ! water
     		if(icarb>0) then
        	rl_w(mi)=1
        	rl_l(mi)=0
        	else
        	rlf(mi)=1.;rvf(mi)=0.
        	endif
     	case(4,6)  !co2_liquid or co2 gas
        	rl_w(mi)=0.     
        	rl_l(mi)=1.
        	rl_v(mi)=0.
! need to set arrays properly for other single phase cases
		end select
        !
     else  ! multi-phase 
     
        !*****CAPILLARY PRESSURE CALCULATION ***********************************   	     
        !			    j=(iphase-1)*max_cp	        
        ! skip this if there is no cap pressure model defined for this type of phase couple			
        ! determine the wetting phase saturation (sw) for this phase couple   (works for 2-phase only)         
!        call wetting_phase(mi,k,sw)
! Assign wetting phase saturation in call to determine phase (k) above
        ! determine the type of cap pressure model
        cap_select=cap_type2(it,k)
	    if(cap_select>0) then 
!	cap_select=cap_type(it)
        	iphase = cap_pos(it,k)
        	j = (iphase - 1) * max_cp
        	if (iphase.le.0) then
           		write(ierr,*) 'no cap pressure model specified for ',k,pn(k)
           		stop
        	end if
! need to insert a special call for fracture case , with flag for fracture  
      	
        	call cap_pressure(cap_select, sw, it,  iphase, mi, cp, dpcp)
!		endif        	
         	hp = cp  ! cap_pressure returns cp in MPa
       	 	dhp = dpcp
        	pcp(mi)=cp
        	if(icarb==0) then
           		dpcef(mi) = dpcp
			else				
           		dpcpw(mi)=dpcp
           ! why do we care about the gas derivative 	?				
           		dpcpg(mi)=-1.*dpcp 	
        	endif
        ! do oil later
        	if(k==27) then
           !     oil/gas
           		pcpog(mi) = cp
           		pcpog(mi+neq) = 0.d0
           		pcpog(mi+2*neq) = dpcp
           		pcpog(mi+3*neq) = 0.d0
        	elseif(k==23) then
           !     oil/water
           		pcpow(mi) = cp
           		pcpow(mi+neq) = dpcp
           		pcpow(mi+2*neq) = 0.d0
           		pcpow(mi+3*neq) = 0.d0
        	end if
     	endif
      !*****RELATIVE PERMEABILITY CALCULATION ***********************************  
! k will always be >=20
! first time, assume all phases will use the same rlp model.  phase_comp(k,1) is the wetting phase
!	write(*,*) 'blah ',it,k
!	do mm=1,2
!	write(*,*) mm, phase_comp(k,mm),rlp_pos(it,phase_comp(k,mm))
!	end do
!	stop
        call model(k,sw,snw,phase_comp(k,1),it,rlw,drlww,drlwg,rll,drllw,drllg,    &
             rlv,drlvw,drlvg,cp,dpcp,mi)

     ! now place results in correct arrays for wetting and non-wetting phases
     select case(k)
        ! need air/water case
        ! these two phase cases never seem to use all 4 derivatives
     case(20,21)   !  water/co2 liquid , water/co2 gas , or liquid co2/gas co2
        rl_w(mi)=rlw
        rl_l(mi)=rll
        drl_ww(mi)=drlww
        drl_wg(mi)=drlwg
        drl_lw(mi)=drllw
        drl_lg(mi)=drllg
     case(22)  !co2_liquid/co2_gas (but icarb>0 ??)
        rl_l(mi)=rlw
        rl_v(mi)=rll
        drllg=drlww
        drlvg=drllw
     case(25, 28)    ! water/air, water/vapor
        rlf(mi) = rlw
        rvf(mi) = rll
        drlef(mi) = drlww
        drvef(mi) = drllw
        drlpf(mi) = 0.0
        drvpf(mi) = 0.0

     case(23,24)  ! oil/water gas/water
        write(ierr,*) 'doing nothing for rel perm of oil/gas'
        ! oil and gas here (13/14)
        !
     case(26)  !coupled water/co2_liquid/co2_gas
        rl_w(mi)=rlw
        rl_l(mi)=rll
        rl_v(mi)=rlv			
        drl_ww(mi)=drlww
        drl_lw(mi)=drllw		
        drl_wg(mi)=drlwg
        drl_lg(mi)=drllg
        drl_vw(mi)=drlvw
        drl_vg(mi)=drlvg
 
     end select
!     if(mi.eq.3) write(*,'(a20,18e15.3)') 'after model ',days, rlf(mi)   !,rvf(mi),drlef(mi),drvef(mi),rlw,rll,drlww,drllw 
	 ity=100
     if(rlp_pos(it,phase_comp(k,2))>0) then  
     ity=rlp_type2(it,phase_comp(k,2))
         if(ity<4) then
         call model(k,1.-sw,1.-snw,phase_comp(k,2),it,rlw,drlww,drlwg,rll,drllw,drllg,    &
             rlv,drlvw,drlvg,cp,dpcp,mi)	
         drlww=-1.*drlww
         else
         call model(k,sw,snw,phase_comp(k,2),it,rlw,drlww,drlwg,rll,drllw,drllg,    &
             rlv,drlvw,drlvg,cp,dpcp,mi)	
         endif
     endif
     ! now place results in correct arrays for non-wetting phase
!	if(ity.eq.0) ity=rlp_type2(it,phase_comp(k,1))  ! this is for the case that no rlp model was specified for the nw phase 
     select case(k)
        ! need air/water case
        ! these two phase cases never seem to use all 4 derivatives
     case(20,21)   !  water/co2 liquid , water/co2 gas , or liquid co2/gas co2  co2_test is k=20
     	if(ity<4) then  ! this logic works for co2_test which invokes a second model call
        	rl_l(mi)=rlw
        	drl_lw(mi)=drlww
        	drl_lg(mi)=drlwg
        else
        	rl_l(mi)=rll
        	drl_lw(mi)=drllw
        	drl_lg(mi)=drllg
		endif        	
     case(22)  !water/co2_gas (but icarb>0 ??)
     	if(ity<4) then
        	rl_v(mi)=rlw
        	drlvg=drlww
        else
        	rl_v(mi)=rlv
        	drlvg=drlvw
        endif
        
     case(25, 28)    ! water/air, water/vapor
     	if(ity<4) then   
        	rvf(mi) = rlw
        	drvef(mi) = drlww
        else
        	rvf(mi) = rll
        	drvef(mi) = drllw
		endif
        drlpf(mi) = 0.0
        drvpf(mi) = 0.0
     case(23,24)  ! oil/water gas/water
        write(ierr,*) 'doing nothing for rel perm of oil/gas'
        ! oil and gas here (13/14)
        !
     case(26)  !water/co2_liquid/co2_gas
		write(ierr,*) ' 2 different relperm models cannot be used for a 3 phase system'
		stop
       
     end select
!	 write(*,*) 'before end of rlp_cap ',mi,rvf(mi)
   
     endif     ! close single / multiphase endif block
     ! proceed to next node      

103 end do
  return
end subroutine rlp_cap
!***************************************************************      

subroutine co2_phase(mi,k)
  ! if co2 problem (1,2, or 3 phase)
  ! input:  node # (mi), rel perm group (it)    
  ! returns:  phase array, k
  use comco2, only : fw,fl,rl_w,rl_l,rl_v,fg,icarb	
  use comrlp, only : rlp_pos	
  use comdi, only : ices
  implicit none      
  integer phase(3),k,mi,it  ,iff,if2
  real*8 sw
  k=1  	
  iff=ices(mi)
  sw=fw(mi)
  if (fw(mi) .ge. 1.0) then
     !     single phase (water)
     k=3
     !				
  else if (fl(mi) .ge. 1.0) then
     !     single phase (liquid/sc CO2)

     k=6  ! does not distinguish between 6 and 9
  else if (fg(mi) .ge. 1.0) then
     !     single phase (gaseous CO2)
     k=4
     ! keating 2012
  else if (fw(mi) .gt. 0.0) then
     !     co2 with water 
     if (ices(mi) .eq. 1 .or. ices(mi) .eq. 4) then
        !     water/co2 liquid 
        k = 20
     else if (ices(mi) .eq. 3) then
        !     water/co2 gas
        k = 21
     else
        !     we're 3-phase
        k=26
     endif
  else
     !     no water (liquid co2 / gas co2)
     k = 22
  endif
!  if(mi.le.20) write(53,566) mi,ices(mi),fw(mi),k
566 format(2i10,f10.3,i10)     	             
  return
end subroutine co2_phase
subroutine oil_phase(mi,it,k,sw) 
  use comai, only : ierr
  use comrlp, only : rlp_pos
  use comco2, only : fw,rl_w,rl_l,fl,fg
  use comcomp, only: rl_g
  implicit none      
  integer phase(3),k,mi,it
  real*8 sw
  k=0
  !     oil/gas
  if (fw(mi) .ge. 1.0) then
     !     single phase (water)
     k=1
     rl_w(mi) = 1.
     sw = fw(mi)
  else if (fl(mi) .ge. 1.0) then
     !     single phase oil
     k=7
     rl_l(mi) = 1.
  else if (fg(mi) .ge. 1.0) then
     !     single phase gas
     k=8
     rl_g(mi) = 1.
  else if (fw(mi) .gt. 0.0) then
     !     oil with water
     sw = fw(mi)
     if (fl(mi).gt.0. .and. fg(mi).eq.0.) then
        !     water/oil 
        k = 23
     else if (fl(mi).eq.0. .and. fg(mi).gt.0.) then
        !     water/gas
        k = 24
     else if (fl(mi).eq.0. .and. fg(mi).eq.0.) then
        !     single phase (water)
        rl_w(mi) = 1.
     else
        !     we're 3-phase
        write(ierr,*) 'do not know how to handle 3-phase oil problem'
        stop
     endif
  endif
  return
end subroutine oil_phase

subroutine water_air_phase(ieosd,k)
  ! for a water/air problem
  ! inputs:  node (mi), relperm group (it), current value of ieosd
  ! returns:  phase array and k
  ! if single phase, will return rlf or rvf	
  use comci, only : rlf,rvf
  use comdi, only : s
  use comrlp, only : rlp_pos
  use comai, only: ico2		
  implicit none      	
  integer phase(3),k,mi,it,ieosd
  real*8 sw
  if(ico2.lt.0) then
    k=25  ! water/air  
  else if (ieosd.eq.1) then
     !     single phase (water)
     k=2
  else if (ieosd .eq. 3) then
     !     single phase (air or vapor)
     k=1
  else if (ieosd.eq.2) then
    if(ico2.eq.0) then
    	k=28   ! water/vapor
    else
  		k=25  ! water/air  
    end if
 endif
  return
end subroutine water_air_phase

subroutine wetting_phase(mi,k,sw)
  ! for a given node (mi) and phase couple (k), returns the wetting phase saturation (sw)		
  ! k=1 single phase     
  ! k=11 air-water or water-vapor    
  ! k:  co2:   2-water/co2 liquid  3-water/co2 gas  4: liquid co2/gas co2 5:three phase
  ! k:  oil:   7-water/oil   9:water/gas   
  use comdi, only : s
  use comco2, only : fw,fl
  use comrlp, only : cap_coupling
  implicit none
  integer icouple,mi,k
  real*8 sw
  !     Assign saturation for wetting phase based on phase couple present

  if(k.eq.25 .or. k.eq.28) then
     ! 2-phase only 
     !  with water		
     !     air/water or vapor/water water/co2
     sw = s(mi)
  else if(k.eq.2.or.k.eq.3)  then
     !     co2/water (liquid or gas) 
     sw = fw(mi)        
  else if (k.eq.4) then
     !     co2 gas/co2 water
     sw = fl(mi)
  else if(k.eq.7) then
     !     gas/liquid (oil problem)
     sw = fl(mi)
  endif
  return
end subroutine wetting_phase
!  ****************************************** MODEL ***********************************************
subroutine model(k,sw, ss, phase, it, rl_w, drl_ww, drl_wg, rl_l, drl_lw, &
     drl_lg, rl_v, drl_vw, drl_vg, hp, dhp, mi)
  ! this routine calculates rel perms and their derivatives for wetting and non-wetting phases, assuming sw is the saturation of the wetting phase
  ! ss is saturation of non-wetting phase, group it, the phase number (phase), and rlp model for that phase (rlp_type2(it,phase)    
  ! 
  use comrlp , only: rlp_pos, rlp_type2, rlp_param	, phase_comp,max_rp
  use comco2 , only: fg	
  use comai, only: ierr	 
  implicit none
  real*8 ss,sw,rl_w,drl_ww,drl_wg,rl_l,drl_lw,drl_lg,rl_v,drl_vw,drl_vg
  real*8 rl,drls,drlg,drlk,drlc,drvw,drvk,drlw,rv,drvg,rt,drtw,drtg,rlc,drtk
  real*8 brl_w,bdrl_ww,bdrl_wg,drl_l,bdrl_lw,bdrl_lg,brl_l,hmin,hp,dhp,hp1,dhp1
  integer iflag,wet,ir,it,k,itbl,irf,mi,m,itype,phase,mm
  ! make sure rlp_pos(it,k) is correctly set in inrlp for the table case, and that the first col is always wetting phase
	phase=abs(phase)
  ir=rlp_pos(it,phase)  ! this points to the correct rlp model for this phase
  itype=rlp_type2(it,phase)  ! this is the rlp model type for this phase
  if(ir>0.and.itype>0) then
     select case(itype)
     case(4,5)   !corey , brooks-corey
        !  ss should be set to co2 liquid sat if 3-phase  	
        call rlp2(itype,sw,ss,it,ir,rl_w,drl_ww,drl_wg,rl_l,drl_lw,drl_lg,    &
             rl_v,drl_vw,drl_vg)
     case(1,2,3)  !  constant, linear, exponential, or 3-phase 
        call rlp1(itype,it,ir,sw,rl_w,drl_ww,rl_l,drl_lw)
        drl_wg = 0;drl_lg=0

     case(6,7,8,9,10)  ! Van Genuchten
        call rlp3(itype, it, ir, mi, sw, hp/9.8e-3, dhp/9.8e-3,  rl_w, drl_ww, rl_l, drl_lw)
        drl_wg = 0; drl_lg=0
     case(11)   ! table
        !   ***********  values from table, called with wetting phase saturation, table number *******
        ! is sw set correctly in 2-phase no water situation?
        irf=2	
        itbl=rlp_pos(it,k)
        call rlp_cap_table (sw, itbl,irf, rl_w, drl_ww)
        irf=3
        call rlp_cap_table (sw, itbl,irf, rl_l, drl_lw)
        if(k==26) then
        	write(ierr,*) 'table format rlpm not available for three phase '
        	stop
        endif
        
     case(12)    
        !*****************     oil_water  ******************************************
        itbl = ir
        irf=3
        call rlp_cap_table (sw, itbl, irf, rl_w, drl_ww)
        drlg = 0.d0
        drlk = 0.d0
        !     oil relative perm in presence of connate water only
        sw = rlp_param(it, ir + 3)
        call rlp_cap_table (sw, itbl, 3, rlc, drlc)	
        ! do not have a place to put rlc   (maybe it is input for stone)              			
        !     gas saturation
        sw = fg(mi)
        !     oil_gas
        itbl = ir
        call rlp_cap_table (sw, itbl, 3, rl_v, drl_vw)
        drvw = 0.d0
        drvk = 0.d0
        !	call stone model
        call stone_rlp(rlc, drlc, mi,    &
             rl, drlw, drlg, drlk, rv, drvw, drvg, drvk,   &
             rt, drtw, drtg, drtk)
        rl_l = rt
        ! passing back rl_l(mi+neq) in  drl_wg   
        !               rl_l(mi+2*neq) in drl_lg
        !               rl_l(mi+3*neq) in drl_vg
        drl_wg = drtw
        drl_lg = drtg
        drl_vg = drtk                                   		
     end select
  else
	call error_v(phase,it)
  endif
  return
end subroutine model

subroutine rlp1(itype,it,ir,sw,rw,drww,rnw,drnw_w)
  use comai, only : ierr
  use comrlp, only : rlp_param,max_rp,rlp_type2
  implicit none
  integer itype,it,ir,iwet,k,ik
  real*8 sw,rl,drls,p1,p2,st1,st2,st3,dummy,drww,drnw_w
  real*8 rw,rnw,p3,p4, akf, akm, porf
  ! this routine handles constant, linear, and exponential rel perm models
  ! all three models return rw,rnw,drww, and drnw
  ! sw is wetting phase saturation
  ! rw - rel perm of wetting phase
  ! rnw - rel perm of non-wetting phase
  ! drww - derivatiove of rw with respect to sw
  ! drnw_w - derivative of rnw with respect to sw

  k=(ir-1)*max_rp          
  p1=rlp_param(it, k+1 )
  p2=rlp_param(it, k+2 )
  select case(itype)

  case(1)  ! constant
     ! only uses p1, only returns value and derivative for one phase
     ! places result in both wetting and non-wetting variables             
     rw=p1
     drww=0.
     rnw=p1
     drnw_w=0.
          
  case(2)
     ! linear               
     ! p1, p2 are residual wetting and residual non-wetting saturations
     call linear(sw, p1, p2, rw, drww)                 
     rnw=1.-rw
     drnw_w=-1.*drww
!	  write(230,*) 'in rlp1 ',sw,rw,rnw
  case(3)
     ! exponential   - provides both wetting and non-wetting values and derivatives
     ! p1, p2 are residual wetting and max wetting saturations 
     ! p3 is exponent; p4 is optional multiplier
     p3=rlp_param(it, k+3 )
     p4=rlp_param(it, k+4 )

     call exponential(sw, p1,p2,p3,p4, rw,drww,rnw,drnw_w)  

  case default
     write(ierr,*) 'cannot use rlp1 subroutine to calculate ',itype,' rel perm type'
     stop  
  end select
  return
end subroutine rlp1
subroutine rlp2(itype,sw,ss,it,ir,rw,drlww,drlwg,rl,drllw,drllg,rv,drlvw,drlvg)   
  use comai, only : ierr
  use comrlp, only : rlp_param,max_rp
  implicit none
  integer itype,it,ir,ir2,k
  real*8 sw,rw,drlww,drlwg,rl,drllw,drllg,rv,drlvw,drlvg,ss,drws
  ! we need to figure out if drws needs to be set              
  ! this subroutine handles 2-phase corey and 2- and 3-phase brooks_corey 
  ! this routine returns rel perm of each phase   
  ! if 2-phase drlwg=0 and drllg=0
  drlwg=0; drllg=0  
  k=(ir-1)*max_rp    
  select case(itype)           
  case(4)               			
     call corey(sw, rlp_param(it, k + 1),    &
          1.-rlp_param(it, k + 2), rw, drlww, rl, drllw)
  case(5)
     call brooks_corey_rlp(sw, rlp_param(it, ir + 1),   &
          1.-rlp_param(it, ir + 2),rlp_param(it, ir + 3),   &
          rlp_param(it, ir + 4),rl, drlww, rl, drllw)
  case(7)
     call brooks_corey_rlp(sw, rlp_param(it, ir + 1),   &
          1.-rlp_param(it, ir + 2),rlp_param(it, ir + 3),   &
          rlp_param(it, ir + 4),rl, drlww, rl, drllw)
     ! 3-phase call still needs work  
     call brooks_corey_rlp3(sw, ss,     &
          rlp_param(it, ir + 1), rlp_param(it, ir + 2),      &
          rlp_param(it, ir  + 3), rl,drllw,drllg,rv,drlvw,drlvg)
     ! not sure if rl,drls,drls2 will be different from the first call                           
  end select

end subroutine rlp2

subroutine rlp3(itype,it,ir,mi,sw,hp,dhp,rw,drww,rnw,drnw_w)
  use comai, only : idpdp, idualp, ierr, neq,days,iupk
  use comco2, only : icarb
  use comrlp, only : rlp_param,max_rp,max_rpf,rlp_type2,rlp_fparam,cap_type,cap_param,cap_fparam,max_cp
  use comdi, only: pnx,pny,pnz,pcp,xfperm,yfperm,zfperm
  implicit none
  integer itype,it,ir,iwet,k,kf,ik, mi, cap_select,mm,iflag,kc
  real*8 sw,rl,drls,rv,drvs,p1,p2,st1,st2,st3,dummy,drww,drnw_w, permb,slr,dummy1
  real*8 rw,rnw,p3,p4, akf, akm, porf, hp, dhp, hp1, dhp1, star, smax,lambda,n
  real*8 hmin, tol_l, tol_u,alphag,nf,slrf,smaxf,drww1,rnw1,rw1,darcyf
  parameter(hmin = 1.d-18)
  !     tol_l and tol_u are lower and upper cutoff saturations
  parameter(tol_l  = 1.d-5)
  parameter(tol_u  = 1.d-5)
    parameter(darcyf = 1.d12)
  ! sw is wetting phase saturation
  ! rw - rel perm of wetting phase
  ! rnw - rel perm of non-wetting phase
  ! drww - derivatiove of rw with respect to sw
  ! drnw_w - derivative of rnw with respect to sw

  k=(ir-1)*max_rp
  kf=(ir-1)*max_rpf
  kc=(ir-1)*max_cp
  slr=rlp_param(it, k+1 )  ! slr
  smax=rlp_param(it, k+2 )  ! smax
  n=rlp_param(it, k+3 )
  lambda=(1.-1./n)  ! convert n to lambda
  if(rlp_fparam(it,kf+1).ge.0.) then ! fracture parameters set
			akf=rlp_fparam(it,kf+4)
			akm=rlp_fparam(it,kf+5)
			porf=rlp_fparam(it,kf+6)
  endif  
  ! Van Genuchten
 select case(itype)
 case (6, 7, 10)   ! Van Genuchten using saturations
       if(itype.eq.6) then  !rnw will be 1-rw
       	iflag=1   
       else
       	iflag=2
       endif
    if (rlp_fparam(it,kf+1).lt.0.) then  ! this will be -99 if fracture parameters were not specified
       ! No fracture model 
! gaz 041424   
       	call vgrlps(iflag, sw,  slr, smax, lambda, tol_l, tol_u, &  
            rw, drww, rnw, dummy)
         drnw_w = -1.*dummy            
		if(itype.eq.10) then
    ! if vg_corey use corey for non-wetting   
          ! don't overwrite wetting phase values, call corey to get non-wetting phase values
          call corey(sw, slr, 1.-smax, rl, drls, rnw, drnw_w)
       	end if
    else   ! fracture parameters were specified
       if (idpdp .ne. 0 .or. idualp .ne. 0) then
          ! Fracture model (double or dual porosity simulation)
          if( mi .le. neq ) then  ! fracture
            call vgrlps(iflag, sw,  rlp_fparam(it, kf + 1), &
                  rlp_fparam(it, kf + 2), 1.-1./rlp_fparam(it, kf + 3), &
                   tol_l, tol_u, rw, drww,rnw,dummy)
             permb = akf * porf    
             if(rlp_fparam(it,kf+7)>0) then ! fracture interaction term was specified
             	call rlp_frac2(1, mi, sw, rlp_fparam(it, kf + 1),rlp_fparam(it, kf + 2),rw, drww, 1.-rw, -drww, &
                  rlp_fparam(it, kf  + 7))
          	 end if
          else  ! matrix
            call vgrlps(iflag, sw,  rlp_param(it, k + 1), &
                  rlp_param(it, k + 2), 1.-1./rlp_param(it, k + 3), &
                   tol_l, tol_u, rw, drww,rnw,dummy)
              permb = akm * (1. - porf)
             if(rlp_fparam(it,kf+7)>0) then ! fracture interaction term was specified
             	call rlp_frac2(1, mi, sw, rlp_param(it, k + 1),rlp_param(it, k + 2), rw, drww, 1.-rw, -drww, &
                  rlp_fparam(it, kf  + 7))
          	 end if
          endif
		else  
            call vgrlps(2, sw,  rlp_param(it, k + 1), &
                  rlp_param(it, k + 2), 1.-1./rlp_param(it, k + 3), &
                   tol_l, tol_u, rw, drww,rnw,dummy)
            call vgrlps(2, sw,  rlp_fparam(it, kf + 1), &
                  rlp_fparam(it, kf + 2), 1.-1./rlp_fparam(it, kf + 3), &
                   tol_l, tol_u, rw1, drww1,rnw1,dummy1)
            permb = akf * porf + akm * (1. - porf)
            rw=(akf*rw1*porf+akm*rw*(1.0-porf))/permb
            drww=(akf*drww1*porf+akm*drww*(1.0-porf))/permb
            rnw = 1. - rw
            drnw_w = -drww
        end if
 
     endif   
  case (8,9)
    if(itype.eq.8) then  !rnw will be 1-rw
       	iflag=0   
       else
       	iflag=1
    endif
    ! Van Genuchten rlp calculation using capillary pressures
    ! this assumes cap pressures (hp,dhp) are already available
    ! if alpha g was specified on rlp it was ignored.  3 parameter model (smin, smax, N)
    if (rlp_fparam(it, kf + 1 ).ge.0.) then
       ! fracture model
 455	format(a15,1x,20(e10.3,1x))	
	  	if (idpdp .ne. 0 .or. idualp .ne. 0) then
			if( mi .le. neq ) then    ! fractures
! vgrlp2 assumes hmin, hp,dhp are defined
       			call vgrlp2(sw, cap_fparam(it, kc + 3 ), rlp_fparam(it, kf + 3 ), hmin, hp, dhp, rw, drww, rnw, drnw_w, iflag)
              permb = akf * porf    
		 	else  ! matrix
       			call vgrlp2(sw, cap_param(it, kc + 3 ), rlp_param(it, k + 3 ), hmin, hp, dhp, rw, drww, rnw, drnw_w, iflag)
              permb = akm * (1. - porf)

			endif 
 	    else
       ! Fracture model, continuum
       		call vgrlp2(sw, cap_param(it, kc + 3), rlp_param(it, k + 3), &
            hmin, hp, dhp, rw, drww, rnw, dummy, 0)
       		call vgrlp2(sw, cap_fparam(it, kc + 3), rlp_fparam(it, kf + 3), &
            hmin, hp, dhp, rl, drls, rv, drvs, 0)
       		permb = akf * porf + akm * (1. - porf)
       		rw = (akf*rl*porf + akm*rw*(1.0-porf))/permb
       		drww = (akf*drls*porf + akm*drww*(1.0-porf))/permb
       		rnw = 1. - rw
       		drnw_w = -drww
     	end if
     else ! no fracture model
       		call vgrlp2(sw, cap_param(it, kc + 3), rlp_param(it, k + 3), &
            hmin, hp, dhp, rw, drww, rnw, drnw_w, iflag)
     
     end if
 case default
     write(ierr,*) 'cannot use rlp3 subroutine to calculate ',itype,' rel perm type'
     stop  
  end select
   if (rlp_fparam(it, kf + 1 ).ge.0.) then
!     upstream weighting of intrinsic permeability
                  if(iupk.gt.0) then
                     permb = permb*darcyf
                     rw = rw*permb
                     rnw = rnw*permb
                     drww = drww*permb
                     drvs = drvs*permb
                     permb = 1.0/darcyf
!     set saturated permeabilities(assume isotropic)
                     pnx(mi) = permb*1.e6
                     pnz(mi) = pnx(mi) * zfperm(it)
                     pny(mi) = pnx(mi) * yfperm(it)
                     pnx(mi) = pnx(mi) * xfperm(it)
                  else
!     set saturated permeabilities(assume isotropic)
                     if(porf.ge.0.0) then
                        pnx(mi) = permb*1.e6
                        pny(mi) = pnx(mi) * yfperm(it)
                        pnz(mi) = pnx(mi) * zfperm(it)
                        pnx(mi) = pnx(mi) * xfperm(it)
                     endif
                  endif
    endif
 	return
end subroutine rlp3

subroutine error_v(ip,it)
  use comai, only : ierr
  implicit none
  integer k,it,ip
  character phase(25)*30
  data phase/'air','water(non co2 prob)','water(co2 prob)','co2 gas','h20 vapor','co2 liquid '   &
       ,'oil','gas','co2 sc','10','methane hydrate',8*'0','water/co2 liquid','water/co2 gas',       &   
       'liquid co2/gas co2','water/oil','water/gas','3 phase'/
  write(ierr,*) phase(ip),' rel perm model is not defined in rlpm group# ' ,it
  stop
end subroutine error_v



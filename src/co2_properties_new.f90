subroutine co2_properties_new(iflg,iphase,var1,var2,var3,istate,var4,var5,var6)
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
!***********************************************************************

! Nov 2  2015 gaz modified from Rajesh Pawar's co2_properties
! gaz 070720  new routine modified from Pawars h2o_properties
      use property_interpolate_2
      use comco2, only: co2_prop
      use comai, only : itsat 
      use comrxni, only : cden_flag

  implicit none
  real*8 mol, mco2, mco22
  real*8 var1,var2,var3,var4,var5(9),var6
  real*8 a1,a2,a3,a4,a5,t1,t2,t3,t4,t5,sum,sum2,t11,t21  
  real*8 a6,x,x2,x3,x4,x5,fn,dfn, ps, p1, t, p, nu, ps1, ps2
  real*8 c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15
  real*8 tc,pc,fg,dfgdp,dfgdt,liq_cp,dliq_cpdp,dliq_cpdt
  real*8 lambdaco2_na, dlambdaco2_nadp, dlambdaco2_nadt
  real*8 tauco2_na_cl, dtauco2_na_cldp, dtauco2_na_cldt
  real*8 temperature, pressure,rhs,drhsdp,drhsdt,dmco2dp,dmco2dt
  real*8 mco21, dmco21dp,dmco21dt,dmco21dxc,dmco2dx,dvar4dt,dvar4dp
  real*8 :: cden_correction
  integer iflg, iphase, iphase1, ifail, istate, icode(9)
  character*200 interpfile, amessage
  !     units are P : MPa
  !     T : degree celsius
  !     h(enthalpy) : MJ/Kg
  !     dhdt : MJ/kg/C
  !     dhdp : MJ/Kg/MPa
  !     rho(density) : kg/m3
  !     dddt : kg/m3/C
  !     dddp : kg/m3/MPa
  !     visc(viscosity) : Pa.s
  !     dvdt: Pa.s/C
  !     dvdp : Pa.s/Mpa

  temperature = var2
  pressure = var1
  var5 = 0.0
  
  if (iflg .eq. 1) then
     ! determine air state
     call get_air_state(ifail,temperature,pressure,istate)
     if(ifail.ne.0) go to 9890
  elseif(iflg.eq.2) then
     ! determine sat temperature given pressure
     call get_air_sat_temperature(ifail,pressure,var4,var5(1))
  elseif(iflg.eq.3) then
     ! determine sat pressure given temprature
     call get_air_sat_pressure(ifail,temperature,var4,var5(1))
  elseif(iflg.eq.4) then
     ! determine properties and derivatives for liquid,vapor,supercritical phases
     icode=1   ! an array
     if(iphase.eq.4) iphase=1
     call get_air_properties(ifail,iphase,9,icode,temperature,pressure,var5)
     ! change units to density-kg/m3 (done) , enthalpy-Mj/kg, and vis-cp
     var5(4)=var5(4)*1.d-3
     var5(5)=var5(5)*1.d-3
     var5(6)=var5(6)*1.d-3
     var5(7)=var5(7)
     var5(8)=var5(8)
     var5(9)=var5(9)
     if(ifail.ne.0) go to 9890
  elseif(iflg.eq.5) then
     ! determine properties and derivatives for liquid phase of 2-phase
     icode=1   ! an array
     call get_air_sat_line_props_pressure(ifail,1,9,icode,pressure,var5)
     var5(4)=var5(4)*1.d-3
     var5(5)=var5(5)*1.d-3
     var5(6)=var5(6)*1.d-3
     var5(7)=var5(7)
     var5(8)=var5(8)
     var5(9)=var5(9)
  elseif(iflg.eq.6) then
     icode=1   ! an array
     call get_air_sat_line_props_pressure(ifail,2,9,icode,pressure,var5)
     var5(4)=var5(4)*1.d-3
     var5(5)=var5(5)*1.d-3
     var5(6)=var5(6)*1.d-3
     var5(7)=var5(7)
     var5(8)=var5(8)
     var5(9)=var5(9)
     ! determine properties and derivatives for vapor phase of 2-phase
  elseif(iflg.eq.7) then
     var5(1) = 0.0
     var5(2) = 0.0
     var5(3) = 0.0
  elseif(iflg.eq.8) then
     if(iphase.eq.1) then
        var5(1) = 1.-var1
        var5(2) = 0.0
        var5(3) = -1.
     else
        var5(1) = var1
        var5(2) = 0.0
        var5(3) = 1.
     endif
  else if (iflg.eq.10) then

  endif

9890 continue

end subroutine co2_properties_new



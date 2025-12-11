subroutine co2_properties(iflg,iphase,var1,var2,var3,istate,var4,var5,var6)
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

  use property_interpolate
  use comco2, only: co2_prop
  use comdi, only: ices
  use comai, only: neq
  use comrxni, only : cden_flag

  implicit none
  real*8 mol, mco2, mco22
  real*8 var1,var2,var3,var4,var5(9)
  real*8 a1,a2,a3,a4,a5,t1,t2,t3,t4,t5,sum,sum2,t11,t21  
  real*8 a6,x,x2,x3,x4,x5,fn,dfn, ps, p1, t, p, nu, ps1, ps2
  real*8 c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15
  real*8 tc,pc,fg,dfgdp,dfgdt,liq_cp,dliq_cpdp,dliq_cpdt
  real*8 lambdaco2_na, dlambdaco2_nadp, dlambdaco2_nadt
  real*8 tauco2_na_cl, dtauco2_na_cldp, dtauco2_na_cldt
  real*8 temperature, pressure,rhs,drhsdp,drhsdt,dmco2dp,dmco2dt
  real*8 mco21, dmco21dp,dmco21dt,dmco21dxc,dmco2dx,dvar4dt,dvar4dp
  real*8 :: cden_correction
  integer iflg, iphase, iphase1, ifail, istate, icode(9), var6
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

  if (iflg .eq. 1) then
     ! determine co2 state
     call get_co2_state(ifail,temperature,pressure,istate)
     if(ifail.ne.0) go to 9890
  elseif(iflg.eq.2) then
     ! determine sat temperature given pressure
     call get_co2_sat_temperature(ifail,pressure,var4,var5(1))
  elseif(iflg.eq.3) then
     ! determine sat pressure given temprature
     call get_co2_sat_pressure(ifail,temperature,var4,var5(1))
  elseif(iflg.eq.4) then
     ! determine properties and derivatives for liquid,vapor,supercritical phases
     icode=1   ! an array
     if(iphase.eq.4) iphase=1
     call get_co2_properties(ifail,iphase,9,icode,temperature,pressure,var5)
     var5(4)=var5(4)*1.d-3
     var5(5)=var5(5)*1.d-3
     var5(6)=var5(6)*1.d-3
     var5(7)=var5(7)*1.d-6
     var5(8)=var5(8)*1.d-6
     var5(9)=var5(9)*1.d-6
     if(ifail.ne.0) go to 9890
  elseif(iflg.eq.5) then
     ! determine properties and derivatives for liquid phase of 2-phase
     icode=1   ! an array
     call get_co2_sat_line_props_pressure(ifail,1,9,icode,pressure,var5)
     var5(4)=var5(4)*1.d-3
     var5(5)=var5(5)*1.d-3
     var5(6)=var5(6)*1.d-3
     var5(7)=var5(7)*1.d-6
     var5(8)=var5(8)*1.d-6
     var5(9)=var5(9)*1.d-6
  elseif(iflg.eq.6) then
     icode=1   ! an array
     call get_co2_sat_line_props_pressure(ifail,2,9,icode,pressure,var5)
     var5(4)=var5(4)*1.d-3
     var5(5)=var5(5)*1.d-3
     var5(6)=var5(6)*1.d-3
     var5(7)=var5(7)*1.d-6
     var5(8)=var5(8)*1.d-6
     var5(9)=var5(9)*1.d-6
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
     ! calculate the dissolved CO2 concentration based on Duan (2003) eqn.
     ! need to calculate fugacity first based on the density.
     c1 = 0.d0
     c2 = 0.d0
     c3 = 0.d0
     c4 = 0.d0
     c5 = 0.d0
     c6 = 0.d0
     c7 = 0.d0
     c8 = 0.d0
     c9 = 0.d0
     c10 = 0.d0
     c11 = 0.d0
     c12 = 0.d0
     c13 = 0.d0
     c14 = 0.d0
     c15 = 0.d0

     !     added fugacity calculation from Duan 2006 paper formulation
     !     equation 2.
     if (var1 .le. 1.d-3) var1 = 1.d-3
     t = var2 + 273.15
     p = var1*10.
     !	  t = var2
     !	  p = var1
     if (cden_flag .eq. 2) then
        ! Sum concentrations to get total moles/kg-water
        mol = cden_correction(var6)
     else
        ! Convert ppm of salt to moles/kg-water
        mol = var3/58.443d3
     end if
     tc = 304.1282d0
     pc = 7.3773d0
     if ((t.gt.273.d0).and.(t.lt.573.d0)) then
        if (t.lt.tc) then
           a1 = -7.0602087d0
           a2 = 1.9391218d0
           a3 = -1.6463597d0
           a4 = -3.2995634
           nu = 1.d0-(t/tc)
           ps1 = a1*nu+(a2*(nu**1.5d0))+(a3*(nu**2.d0))+(a4*(nu**4.d0))
           ps2 = a1+(1.5d0*a2*(nu**0.5d0))+(2.0d0*a3*nu)+(4.d0*a4*(nu**3.d0))
           ps = ps1*tc/t
           ps = dexp(ps)*pc
           p1 = ps
        else if ((t.gt.tc).and.(t.lt.405.d0)) then
           p1 = 75.d0+(t-tc)*1.25d0
        else if (t.gt.405.d0) then
           p1 = 200.d0
        endif
        if( p.lt.p1) then
           c1 = 1.d0
           c2 = 4.7586835d-3
           c3 = -3.3569963d-6
           c5 = -1.3179356
           c6 = -3.8389101d-6
           c8 = 2.2815104d-3
        elseif (t.lt.340.d0) then
           if(p.lt.1000) then
              c1 = -7.1734882d-1
              c2 = 1.5985379d-4
              c3 = -4.9286471d-7
              c6 = -2.7855285d-7
              c7 = 1.1877015d-9
              c12 = -96.539512
              c13 = 4.4774938d-1
              c14 = 101.81078d0
              c15 = 5.3783879d-6
           else
              c1 = -6.5129019d-2
              c2 = -2.1429977d-4
              c3 = -1.144493d-6
              c6 = -1.1558081d-7
              c7 = 1.195237d-9
              c12 = -221.34306d0
              c14 = 71.820393d0
              c15 = 6.6089246d-6
           endif
        elseif (t.lt.435.d0) then
           if (p.lt.1000) then
              c1 = 5.0383896d0
              c2 = -4.4257744d-3
              c4 = 1.9572733d0
              c6 = 2.4223436d-6
              c8 = -9.3796135d-4
              c9 = -1.502603d0
              c10 = 3.027224d-3
              c11 = -31.377342d0
              c12 = -12.847063d0
              c15 = -1.5056648d-5
           else
              c1 = -16.063152d0
              c2 = -2.705799d-3
              c4 = 1.4119239d-1
              c6 = 8.1132965d-7
              c8 = -1.1453082d-4
              c9 = 2.3895671d0
              c10 = 5.0527457d-4
              c11 = -17.76346d0
              c12 = 985.92232d0
              c15 = -5.4965256d-7
           endif
        else
           c1 = -1.569349d-1
           c2 = 4.4621407d-4
           c3 = -9.1080591d-7
           c6 = 1.0647399d-7
           c7 = 2.4273357d-10
           c9 = 3.5874255d-1
           c10 = 6.3319710d-5
           c11 = -249.89661d0
           c14 = 888.768d0
           c15 = -6.6348003d-7
        endif
     endif

	if(ices(var6).eq.3) then
		call co2_fugacity(var1,t,co2_prop(9*neq+var6),fg,dfgdp,dfgdt)
	else
		call co2_fugacity(var2,t,co2_prop(var6),fg,dfgdp,dfgdt)
	endif
     liq_cp = 28.9447706d0 + (-0.0354581768d0*t) + (-4770.67077d0/t)  &
          +(1.02782768d-5*t*t) + (33.8126098/(630-t)) + (9.0403714d-3*p)  &
          +(-1.14934031d-3*p*dlog(t)) + (-0.307405726*p/t)  &
          + (-0.0907301486*p/(630.d0-t))   &
          + (9.32713393d-4*p*p/((630.d0-t)**2.d0))

     dliq_cpdp = 9.0403714d-3-1.14934031d-3*dlog(t)-(0.307405726/t) &
          -0.0907301486/(630.d0-t)+ (2.d0*9.32713393d-4*p/((630.d0-t)**2.d0))

     dliq_cpdt = -0.0354581768d0 + (4770.67077d0/(t*t))  &
          +(2.d0*1.02782768d-5*t) + (33.8126098/(630-t)**2)  &
          +(-1.14934031d-3*p/t) + (0.307405726*p/(t*t))  &
          + (-0.0907301486*p/(630.d0-t)**2)   &
          + (2.d0*9.32713393d-4*p*p/((630.d0-t)**3.d0))
! changed the above line not multiplied by 2 before
     lambdaco2_na = -0.411370585d0 + (6.07632013d-4*t) &
          + (97.5347708d0/t) + (-0.0237622469d0*p/t)  &
          + (0.0170656236*p/(630.d0-t)) + (1.41335834d-5*t*dlog(p))

     dlambdaco2_nadp =  (-0.0237622469d0/t)  &
          + (0.0170656236/(630.d0-t)) + (1.41335834d-5*t/p)

     dlambdaco2_nadt =  (6.07632013d-4) &
          + (-97.5347708d0/(t*t)) + (+0.0237622469d0*p/(t*t))  &
          + (0.0170656236*p/(630.d0-t)**2.) + (1.41335834d-5*dlog(p))

     tauco2_na_cl = 3.36389723d-4 + (-1.9829898d-5*t)  &
          + (2.1222083d-3*p/t) + (-5.24873303d-3*p/(630.d0-t))


     dtauco2_na_cldp = (2.1222083d-3/t) + (-5.24873303d-3/(630.d0-t))


     dtauco2_na_cldt = (-1.9829898d-5)  &
          + (-2.1222083d-3*p/(t*t)) + (-5.24873303d-3*p/(630.d0-t)**2.)

     rhs = liq_cp + (2.d0*lambdaco2_na*mol) + (tauco2_na_cl*mol*mol)

     drhsdp = dliq_cpdp + 2.d0*dlambdaco2_nadp*mol &
          + dtauco2_na_cldp*mol*mol

     drhsdt = dliq_cpdt + 2.d0*dlambdaco2_nadt*mol &
          + dtauco2_na_cldt*mol*mol

!     drhsdp = drhsdp*10.d0

!     dvar4dp = var5(2)
!     dvar4dt = var5(3)
     dvar4dp = 0.d0
     dvar4dt = 0.d0

     mco21 = 	var4*fg*p*44.d-3
     mco22 = var4*fg*p*dexp(-rhs)
     dmco21dp = 44.d-3*(var4*(fg+p*dfgdp)+fg*p*dvar4dp)
     dmco21dt = 44.d-3*p*(var4*dfgdt+fg*dvar4dt)
     dmco21dxc = fg*p*44.d-3
     mco2 = mco21/(mco21+dexp(rhs))

     dmco2dp = (dexp(rhs)*(dmco21dp-mco21*drhsdp))/ &
          ((mco21+dexp(rhs))*(mco21+dexp(rhs)))
     dmco2dt = (dexp(rhs)*(dmco21dt-mco21*drhsdt))/ &
          ((mco21+dexp(rhs))*(mco21+dexp(rhs)))
     dmco2dx = (dexp(rhs)*dmco21dxc)/ &
          ((mco21+dexp(rhs))*(mco21+dexp(rhs)))


!	 mco2 = 0.016*var1
!	 dmco2dp = 0.0016d0
!	 dmco2dt = 0.00d0
!	 dmco2dx = 0.d0
     var5(1) = mco2
     var5(2) = dmco2dp*10.d0
     var5(3) = dmco2dt
     var5(4) = dmco2dx
     var5(5) = fg
     var5(6) = rhs

  endif
9890 continue
  !      call get_error_message(amessage)
  !9891  write(6,*)
  !      call writmess(6,amessage)
  !      go to 9900




end subroutine co2_properties

subroutine writmess(iunit,amessage)

  implicit none

  integer iunit,jend,i,nblc,junit,leadblank,itake,j
  character*(*) amessage
  character (len=20) ablank

  ablank=' '
  itake=0
  j=0
  junit=iunit

  amessage=' '//trim(amessage)
  if(amessage.eq.' ')then
     write(junit,*)
     return
  end if
  do i=1,min(20,len(amessage))
     if(amessage(i:i).ne.' ')go to 21
20 end do
21 leadblank=i-1
  nblc=len_trim(amessage)
5 jend=j+78-itake
  if(jend.ge.nblc) go to 100
  do i=jend,j+1,-1
     if(amessage(i:i).eq.' ') then
        if(itake.eq.0) then
           write(junit,'(a)') amessage(j+1:i)
           itake=2+leadblank
        else
           write(junit,'(a)') ablank(1:leadblank+2)//    &
                amessage(j+1:i)
        end if
        j=i
        go to 5
     end if
  end do
  if(itake.eq.0)then
     write(junit,'(a)') amessage(j+1:jend)
     itake=2+leadblank
  else
     write(junit,'(a)') ablank(1:leadblank+2)//     &
          amessage(j+1:jend)
  end if
  j=jend
  go to 5
100 jend=nblc
  if(itake.eq.0)then
     write(junit,'(a)') amessage(j+1:jend)
  else
     write(junit,'(a)') ablank(1:leadblank+2)//    &
          amessage(j+1:jend)
  end if
  return


end subroutine writmess



subroutine writint(atemp,ival)

  !       Subroutine WRITINT writes an integer to a character variable.

  integer*4 ival
  character*6 afmt
  character*(*) atemp

  afmt='(i   )'
  write(afmt(3:5),'(i3)') len(atemp)
  write(atemp,afmt)ival
  atemp=adjustl(atemp)
  return
end subroutine writint



SUBROUTINE REALREAD(IFAIL,CLINE,RTEMP)

  ! -- Subroutine REALREAD reads a real number from a string.

  INTEGER IFAIL
  real RTEMP
  CHARACTER*8 AFMT
  CHARACTER*(*) CLINE

  IFAIL=0
  AFMT='(F   .0)'
  WRITE(AFMT(3:5),'(I3)') len_trim(CLINE)
  READ(CLINE,AFMT,ERR=100) RTEMP
  RETURN

100 IFAIL=1
  RETURN
END SUBROUTINE REALREAD

SUBROUTINE intREAD(IFAIL,CLINE,iTEMP)

  ! -- Subroutine intREAD reads a real number from a string.

  INTEGER IFAIL
  integer iTEMP
  CHARACTER*6 AFMT
  CHARACTER*(*) CLINE

  IFAIL=0
  AFMT='(i   )'
  WRITE(AFMT(3:5),'(I3)') LEN(CLINE)
  READ(CLINE,AFMT,ERR=100) iTEMP
  RETURN

100 IFAIL=1
  RETURN
END SUBROUTINE intREAD

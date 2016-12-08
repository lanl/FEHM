      subroutine ther_co2_h2o(iflg,ndummy)
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

!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To calculate the equation coeffients and derivatives for a
!D1 air-h2o-CO2 simulation.supercritical-liquid-gas system.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30 
!D2 
!D2 Initial implementation: Date 05-Feb-07, Programmer: Rajesh Pawar
!**********************************************************************
!D3 
!D3
!D3 
!D3
!D3
!D3
!D3 
!**********************************************************************

      use comai
      use combi
      use comci
      use comco2
      use comdi
      use comdti
      use comei
      use comfi
      use comgi
      use comii
      use commeth
      use comriv
      use comrlp, only : rlpnew
      use comrxni
      implicit none

      integer ndummy,mid,mi,ieosd,kq,icesd
      real*8 dtin,dporpl,dportl
      real*8 tl,pl,dpsatt
c      real*8 enw,enl
c      real*8 dhlt,dhlp
c      real*8 rol
c      real*8 drolp,drolt
      real*8 xvisl
c      real*8 dvislp
c      real*8 dvislt
c      real*8 env
c      real*8 dhvt,dhvp
c      real*8 rov
      real*8 xvisv
      real*8 qdis,dtd
      real*8 cp,por,vol,pldif,eskd,eqdum
      real*8 damp,daep,daep1,damh,daeh
      real*8 dql,hprod,dhprdp,dhprde,htc,tbound,hflux
      real*8 dhflxp,sbound,cprd,edif,denrd
      real*8 dpldt,daep2
      real*8 dqv,dhflxe,dtps,permsd1
      real*8 dhflxem,dhflxpm,daeh1,daeh2,pcpww
      real*8, allocatable :: sto1(:)

      real*8 psatd
      real*8 dum1,dum2,dum3,dum5,dum6
      real*8 damw,daew,dhprdw                        
      real*8 dhflxew
      real*8 ensc,denscp,densct,rosc,droscp,drosct,xvisc,dviscp,dvisct
      real*8 dtol
      parameter(dtol=1.d-10)
c     
c     gaz took it out of loop (can be specified source as well)
c     will over-ride above co2 production
      real*8 fracwmin
      real*8 value(9),dumb
      
c      real*8 drolyc, drolya, roa, droadp, droadt, ena, denadt,denadp
      real*8   ycp, dycp
c      real*8  rlw, drlww, drlwg, drlwp, drlwt
c      real*8 rll, drllw, drllg, drllp, drllt
c      real*8 drolw, denlt, denlp, denlw, denlya, denlyc
c      real*8 visl, dvislya, dvislyc, dvislw
c      real*8 row, drowp, drowt, drowya, drowyc, denwp, denwt
c      real*8 visw, dviswp, dviswt
c      real*8 dprmp, dprmt, dprmw, dprmyc, dprmya
c      real*8 damyc, damya, daeyc, daeya
c      real*8 dhprdyc, dhprdya, rlv, drlvw, drlvg, drlvp, drlvt
c      real*8 drovw, drovya, drovyc, denvw, denvya, denvyc, visv, dvisvya
c      real*8 dvisvyc, dvisvw, demwyc, denwyc, denwya
c      real*8 vpartial,dvpardt
c      real*8 xs, dxsw, dxsg, denwfw, denwfg, s1, s2, ds1dw, ds1dg
c new variables
c      real*8 enx,denxp,denxe,denxyc,denxya, vis_tol
      real*8 :: cden_correction, mol
c      real*8 :: permsd11 = 0., dprmp1 = 0., dprmt1 = 0., dprmw1 = 0.
c      real*8 :: permsd12 = 0., dprmp2 = 0., dprmt2 = 0., dprmw2 = 0.
c      real*8 :: permsd13 = 0., dprmp3 = 0., dprmt3 = 0., dprmw3 = 0.
c      parameter (fracwmin=0.1,vis_tol = 1.d-12)
      integer iflg,duma
      
c      save dprmya

      allocate(sto1(n0*2))
c     
c     rol  -  density liquid
c     rov  -  density vapour
c     enl  -  enthalpy liquid
c     env  -  enthalpy vapour
c     visl -  viscosity of liquid
c     visv -  viscosity of vapour
c     rl   -  relative permeability of liquid phase
c     rv   -  relative permeability of vapour phase
c     tfun -  temperature
c     sw   -  saturation liquid
c
      fracwmin=0.1     
      vis_tol = 1.d-12
      dtin=1.0/dtot
      mi = l

c     ****************************************************************
      if(iflg.eq.0) then
c     calculations that are common to co2, co2 hydrate, and water
c     calculate pressure dependant porosity and derivatives

         if(iporos.ne.0) call porosi(1)

c     RJP 02/08/07 Modified the array names below. 
c     relative permeability for multi-phase CO2/water case
c     rl_l represents relative permeability for liquid CO2 
c     in 2-phase CO2 case,or for super-critical CO2
c     derivatives wrt water saturation are in w while wrt CO2 (g/sc)
c     are in g
         allocate(rl_l(n0))
         allocate(drl_lw(n0))
         allocate(drl_lg(n0))
         allocate(rl_v(n0))
         allocate(drl_vw(n0))
         allocate(drl_vg(n0))
         allocate(rl_w(n0))
         allocate(drl_ww(n0))
         allocate(drl_wg(n0))
         rl_w = 0.0
         rl_l = 0.0
         rl_v = 0.0
         drl_lw = 0.0
         drl_lg = 0.0
         drl_vw = 0.0
         drl_vg = 0.0
         drl_ww = 0.0
         drl_wg = 0.0
         
         if (rlpnew) then
            call rlp_cap(ndummy)
            do mid=1,neq
               mi=mid+ndummy
               phi(mi)=phico2(mi)-pcp(mi)
            end do
         else
            do mid=1,neq
               mi=mid+ndummy
c     calculate multi-phase relative perms.
               call rlperm_co2(ndummy,0,mi,rl_w(mi),
     &           drl_ww(mi),drl_wg(mi),rl_l(mi),drl_lw(mi),drl_lg(mi),
     &           rl_v(mi),drl_vw(mi),drl_vg(mi))
c     calculate multi-phase cap. pres.
               call rlperm_co2(ndummy,1,mi,pcp(mi),
     &              dpcpw(mi),dpcpg(mi),pcg(mi),dpcgw(mi),dpcgg(mi),
     &              dum1,dum2,dum3)
               phi(mi)=phico2(mi)-pcp(mi)
            enddo
         end if
      else if(iflg.eq.11) then
c     
c     RJP 02/09/07. All the thermodynamic properties are now calculated 
c     here and stored in arrays for later.	 
c     
c     con_prop are user-defined properties read in as consts in icectrco2.f
c
         if(ico2prop_flg.eq.1) then
            rosc=con_prop(1)
            droscp=con_prop(2)
            drosct=con_prop(3)
            ensc=con_prop(4)
            denscp=con_prop(5)
            densct=con_prop(6)
            xvisc=con_prop(7)
            dviscp=con_prop(8)
            dvisct=con_prop(9)

            row=con_prop(10)
            drowp=con_prop(11)
            drowt=con_prop(12)
            enw=con_prop(13)
            denwp=con_prop(14)
            denwt=con_prop(15)
            visw=con_prop(16)
            dviswp=con_prop(17)
            dviswt=con_prop(18)
         endif

         mi=ndummy
c     phase information now contained in array ices
         ieosd=ieos(mi)
         pl=phico2(mi)
         tl=t(mi)
         dtps=dtpsc(mi)
c     
c     terms to account for space taken up with water-rich phase & co2-rich phase
c     
         frac_w = fw(mi)
         frac_cl=fl(mi)
         frac_cg=fg(mi)
         yco2 = yc(mi)
         ywat = yw(mi)
         yair = ya(mi)
         xco2 = xc(mi)
         xwat = xw(mi)
         xair = xa(mi)
         dporpl=dporp(mi)
         dportl=dport(mi)
            
         if(ico2prop_flg.eq.1) then
            wat_prop(mi) = row
            wat_prop(neq+mi) = drowp
            wat_prop(2*neq+mi) = drowt
            wat_prop(3*neq+mi) = 0.d0
            wat_prop(4*neq+mi) = 0.d0
            wat_prop(5*neq+mi) = enw
            wat_prop(6*neq+mi) = denwp
            wat_prop(7*neq+mi) = denwt
            wat_prop(8*neq+mi) = visw
            wat_prop(9*neq+mi) = dviswp
            wat_prop(10*neq+mi) = dviswt
            call co2_properties(1,duma,pl,tl,dum1,icesd,dumb,value,
     &           duma)
            if(icesd.eq.3) then
               co2_prop(9*neq+mi) = rosc
               co2_prop(10*neq+mi) = drosct
               co2_prop(11*neq+mi) = droscp
               co2_prop(12*neq+mi) = ensc
               co2_prop(13*neq+mi) = densct
               co2_prop(14*neq+mi) = denscp
               co2_prop(15*neq+mi) = xvisc
               co2_prop(16*neq+mi) = dvisct
               co2_prop(17*neq+mi) = dviscp
            else
               co2_prop(mi) = rosc
               co2_prop(neq+mi) = drosct
               co2_prop(2*neq+mi) = droscp
               co2_prop(3*neq+mi) = ensc
               co2_prop(4*neq+mi) = densct
               co2_prop(5*neq+mi) = denscp
               co2_prop(6*neq+mi) = xvisc
               co2_prop(7*neq+mi) = dvisct
               co2_prop(8*neq+mi) = dviscp
            endif

         else
c     vapor enthalpy and derivatives

            call h2o_properties(1,3,pl,tl,dum2,dum3,
     &           env,dhvp,dhvt,dum5,dum6)
               
c     vapor density and derivatives
               
            call h2o_properties(2,3,pl,tl,dum2,dum3,
     &           rov,drovp,drovt,dum5,dum6)

c     vapor viscosity and derivatives
               
            call h2o_properties(3,3,pl,tl,dum2,dum3,
     &           xvisv,dvisvp,dvisvt,dum5,dum6)

               
            if (cden_flag .eq. 2) then
               mol = cden_correction(mi)
c     liquid enthalpy and derivatives
               call h2o_properties(1,2,pl,tl,mol,dum3,
     &              enl,dhlp,dhlt,dum5,dum6)
c     liquid density and derivatives
               call h2o_properties(2,2,pl,tl,mol,dum3,
     &              rol,drolp,drolt,dum5,dum6)
c     liquid viscosity and derivatives               
               call h2o_properties(3,2,pl,tl,mol,dum3,
     &              xvisl,dvislp,dvislt,dum5,dum6)
            else
c     liquid enthalpy and derivatives
               call h2o_properties(1,2,pl,tl,csalt(mi),dum3,
     &              enl,dhlp,dhlt,dum5,dum6)
c     liquid density and derivatives
               call h2o_properties(2,2,pl,tl,csalt(mi),dum3,
     &              rol,drolp,drolt,dum5,dum6)
c     liquid viscosity and derivatives
               call h2o_properties(3,2,pl,tl,csalt(mi),dum3,
     &              xvisl,dvislp,dvislt,dum5,dum6)
            end if

               

c     Calculate density of CO2 dissolved brine. It is assumed the enthalpy
c     & viscosity of brine do not change due to CO2 dissolution. It is 
c     also assumed that the density does not change due to air dissolution.
c     zvd 07-Aug-08 added following lines
            drolyc=0.d0
            drolya=0.d0
c     zvd 07-Aug-08
            if(iprtype.ge.4) then
               vpartial=37.51d0-(9.585d-2*tl)+(8.74d-4*(tl**2.d0))-
     &              (5.044d-7*(tl**3.d0))
               vpartial=vpartial*1d-6
               dVpardt = -9.585d-2+(2.d0*8.74d-4*tl)-
     &              (3.d0*5.044d-7*(tl**2.d0)) 
               dVpardt=dVpardt*1d-6
               rol_h2o = rol
c     calculate effective molecular weight
               emw = yco2*44.d-3+ywat*18.d-3
               demwyc=26.d-3
               ycp = yco2/44.d-3/(yco2/44.d-3+ywat/18.d-3)
               dycp = ywat/(44.d-3*18.d-3)/(yco2/44.d-3+ywat/18.d-3)**2
               rol_d = 1.d0-((44.d-3-rol_h2o*vpartial)*ycp/emw)
               drol_dp = drolp*vpartial*ycp/emw
               drol_dt = ycp*(drolt*vpartial+rol_h2o*dVpardt)/emw
               drol_dyc = -(44.d-3-rol_h2o*vpartial)*(emw*dycp-ycp*
     &              demwyc)/(emw*emw)
c     if(ico2dis(mi).eq.1) then
c     drol_dp=drol_dp+drol_dyc*dmol(mi)
c     drol_dt=drol_dt+drol_dyc*dmol(mi+neq)
c     endif
               drol_dya = 0.0
               rol = rol_h2o/rol_d
               drolp = (rol_d*drolp-rol_h2o*drol_dp)/(rol_d*rol_d)
               drolt = (rol_d*drolt-rol_h2o*drol_dt)/(rol_d*rol_d)
               drolyc = -(rol_h2o*drol_dyc)/(rol_d*rol_d)
               drolya = 0.0
            endif           
c     modify derivatives for 2-phase (should now depend on p, not t)

            if(abs(ieosd).eq.2) then
               drolp=drolp+drolt*dtps
               drovp=drovp+drovt*dtps
               dhlp=dhlp+dhlt*dtps
               dhvp=dhvp+dhvt*dtps
               dvisvp=dvisvp+dvisvt*dtps
               dvislp=dvislp+dvislt*dtps
               dporpl=dporpl+dportl*dtps
               drolt=0.0
               drovt=0.0
               dhlt=0.0
               dhvt=0.0
               dvisvt=0.0
               dvislt=0.0
               dportl=0.0
               dtd=1.0
            endif

            wat_prop(mi) = rol
c Hari
            rolf(mi) = rol
            wat_prop(neq+mi) = drolp
            wat_prop(2*neq+mi) = drolt
            wat_prop(3*neq+mi) = drolyc
            wat_prop(4*neq+mi) = drolya
            wat_prop(5*neq+mi) = enl
            wat_prop(6*neq+mi) = dhlp
            wat_prop(7*neq+mi) = dhlt
            wat_prop(8*neq+mi) = xvisl
            wat_prop(9*neq+mi) = dvislp
            wat_prop(10*neq+mi) = dvislt
c     
c     air properties
c     
            roa = 1.292864d0*(273.15d0/(tl+273.15d0))*(pl/0.101325d0)
            droadp = 1.292864d0*(273.15d0/(tl+273.15d0))/0.101325d0
            droadt = -1.292864d0*(273.15d0/((tl+273.15d0)*
     &           (tl+273.15d0)))*(pl/0.101325d0)
            ena = (1003.7d0+0.025656d0*tl+0.00045457d0*tl*
     &           tl-2.7107d-7*tl*tl*tl)*tl*1.d-6
            denadt = (1003.7d0+2.d0*0.025656d0*tl+3.d0*
     &           0.00045457d0*tl*tl-4.d0*2.7107d-7*tl*tl*tl)*1.d-6
            denadp = 0.d0
            visca = 1.82d-5
            dvisadp = 0.d0
            dvisadt = 0.d0
c     
c     RJP 02/09/07 Calculate CO2-rich phase properties. We are assuming
c     the properties of CO2-rich phase are mol fraction weighted properties
c     of CO2+air+water vapor
c     

            icesd=ices(mi)
            tl=tco2(mi)
            if (fg(mi) .eq. 0.) then
               xs = 0.d0
               dxsg = 0.d0
               dxsw = 0.d0
            else
               xs = fg(mi)/(fg(mi)+fl(mi))
               dxsg = 1.d0/(fg(mi)+fl(mi))+fg(mi)/(fg(mi)+fl(mi))**2
               if (fw(mi) .eq. 1.) then
                  dxsw = 0.
               else
                  dxsw = fg(mi)/((1-fw(mi))*(1-fw(mi)))
               end if
            end if
            if(icesd.eq.2) then
c     two phase (gas/liquid) conditions
c     calculate phase-change temperature and dt/dp
               call co2_properties(2,icesd,pl,dumb,dum1,duma,tl,
     &              value,duma)
               dtps=value(1)
               dtpsc(mi)=dtps
c     RJP 01/22/09 added following for modifying derivatives for
c     dissolved CO2 mass fraction for two phase
               dmol(mi) = dmol(mi)+dmol(mi+neq)*dtps
               dmol(mi+neq) = 0.d0
c     tco2(mi) = tl
            endif
               
            if(icesd.ne.2) then
               call co2_properties(4,icesd,pl,tl,dum1,duma,dumb,
     &              value,duma)
               if(icesd.eq.3) then
                  co2_prop(9*neq+mi) = value(1)
                  co2_prop(10*neq+mi) = value(2)
                  co2_prop(11*neq+mi) = value(3)
                  co2_prop(12*neq+mi) = value(4)
                  co2_prop(13*neq+mi) = value(5)
                  co2_prop(14*neq+mi) = value(6)
                  co2_prop(15*neq+mi) = value(7)
                  co2_prop(16*neq+mi) = value(8)
                  co2_prop(17*neq+mi) = value(9)
               else
                  co2_prop(mi) = value(1)
                  co2_prop(neq+mi) = value(2)
                  co2_prop(2*neq+mi) = value(3)
                  co2_prop(3*neq+mi) = value(4)
                  co2_prop(4*neq+mi) = value(5)
                  co2_prop(5*neq+mi) = value(6)
                  co2_prop(6*neq+mi) = value(7)
                  co2_prop(7*neq+mi) = value(8)
                  co2_prop(8*neq+mi) = value(9)
               endif
               dtpaco2(mi)=0.d0
               dtpaeco2(mi)=1.d0
c     RJP 03/02/08 added mixture enthalpy for water. added the enthalpy derivative wrt yc and ya
c     to wat_prop array. As the value of wat_prop(5*neq+mi) is getting modified, the derivative
c     is calculated first.
               wat_prop(11*neq+mi) = value(4)-wat_prop(5*neq+mi)
               wat_prop(12*neq+mi) = -wat_prop(5*neq+mi)
               wat_prop(5*neq+mi) = ywat*wat_prop(5*neq+mi)+yco2*
     &              value(4)
               wat_prop(6*neq+mi) = ywat*wat_prop(6*neq+mi)+yco2*
     &              value(6)
               wat_prop(7*neq+mi) = ywat*wat_prop(7*neq+mi)+yco2*
     &              value(5)
c     wat_prop(11*neq+mi) = 0.d0
c     wat_prop(12*neq+mi) = 0.d0
            else
c     RJP 2/12/07. It is assumed in case of two phase CO2 (l+g) in equilibrium
c     with water, the liquid phase rich is CO2 is pure CO2 and water vapor and air
c     are part of the gas phase.
               call co2_properties(5,icesd,pl,tl,dum1,duma,dumb,
     &              value,duma)
               co2_prop(mi) = value(1)
               co2_prop(neq+mi) = 0.d0
               co2_prop(2*neq+mi) = value(3)
               co2_prop(3*neq+mi) = value(4)
               co2_prop(4*neq+mi) = 0.d0
               co2_prop(5*neq+mi) = value(6)
               co2_prop(6*neq+mi) = value(7)
               co2_prop(7*neq+mi) = 0.d0
               co2_prop(8*neq+mi) = value(9)
               call co2_properties(6,icesd,pl,tl,dum1,duma,dumb,
     &              value,duma)
C     RJP 04/08/07 assumed density of co2-rich phase is independent of mixed water & air.
               co2_prop(9*neq+mi) = value(1)
               co2_prop(10*neq+mi) = 0.d0
               co2_prop(11*neq+mi) = value(3)
               co2_prop(12*neq+mi) = value(4)
               co2_prop(13*neq+mi) = 0.d0
               co2_prop(14*neq+mi) = value(6)
               co2_prop(15*neq+mi) = value(7)
               co2_prop(16*neq+mi)=0.d0
               co2_prop(17*neq+mi)=value(9)	
c     RJP 03/02/08 added mixture enthalpy for water
               wat_prop(11*neq+mi) = value(4)-wat_prop(5*neq+mi)
               wat_prop(12*neq+mi) = -wat_prop(5*neq+mi)
               wat_prop(13*neq+mi) = dxsg*yco2*
     &              (value(4)-co2_prop(3*neq+mi))
               wat_prop(14*neq+mi) = dxsw*yco2*
     &              (value(4)-co2_prop(3*neq+mi))
               wat_prop(5*neq+mi) = ywat*wat_prop(5*neq+mi)+yco2*
     &              (value(4)*xs+(1.d0-xs)*co2_prop(3*neq+mi))
               wat_prop(6*neq+mi) = ywat*wat_prop(6*neq+mi)+yco2*
     &              (value(6)*xs+(1.d0-xs)*co2_prop(5*neq+mi))
               wat_prop(7*neq+mi) = ywat*wat_prop(7*neq+mi)+yco2*0.d0
               wat_prop(neq+mi) = wat_prop(neq+mi)+
     &              wat_prop(2*neq+mi)*dtpsc(mi)
               wat_prop(2*neq+mi) = 0.d0
               wat_prop(6*neq+mi) = wat_prop(6*neq+mi)+
     &              wat_prop(7*neq+mi)*dtpsc(mi)
               wat_prop(7*neq+mi) = 0.d0
               wat_prop(9*neq+mi) = wat_prop(9*neq+mi)+
     &              wat_prop(10*neq+mi)*dtpsc(mi)
               wat_prop(10*neq+mi) = 0.d0
               dtpaco2(mi)=dtpsc(mi)
               dtpaeco2(mi)=0.d0

c     wat_prop(11*neq+mi) = 0.d0
c     wat_prop(12*neq+mi) = 0.d0
            endif
         endif

c     ****************************************************************
c     
      else if(iflg.eq.-1) then
c     release memory
         deallocate(rl_l)
         deallocate(drl_lw)
         deallocate(drl_lg)
         deallocate(rl_v)
         deallocate(drl_vw)
         deallocate(drl_vg)
         deallocate(rl_w)
         deallocate(drl_ww)
         deallocate(drl_wg)

c     ****************************************************************
      else if(iflg.eq.1) then
c     
c     calculations for water
c     
         do mid=1,neq
            mi=mid+ndummy
            ieosd=ieos(mi)
            pw=phi(mi)
            tl=t(mi)
            pl=phico2(mi)
            pcpww=dpcpw(mi)
            icesd = ices(mi)
            frac_w = fw(mi)
            frac_c=fl(mi)
            if((icesd.eq.3).or.(icesd.eq.2)) frac_c=fg(mi)
            yco2 = yc(mi)
            ywat = yw(mi)
            yair = ya(mi)
            xco2 = xc(mi)
            xwat = xw(mi)
            xair = xa(mi)		
               diw(mi)=0.0
               diwp(mi)=0.0
               diwe(mi)=0.0
               diww(mi)=0.0
               diwyc(mi)=0.0
               diwya(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divw(mi)=0.0
               divyc(mi)=0.0
               divya(mi)=0.0
            cp=denr(mi)*cpr(mi)
            por=ps(mi)
            dporpl=dporp(mi)
            dportl=dport(mi)
            vol=volume(mi)
            dtps=dtpsc(mi)
c     rel-perms  
            rlw = rl_w(mi)
            drlww = drl_ww(mi)
            drlwg = drl_wg(mi)
            drlwp = 0.0
            drlwt = 0.0
            rll = rl_l(mi)
            drllw = drl_lw(mi)
            drllg = drl_lg(mi)
            drllp = 0.0
            drllt = 0.0

            dmpf(mi)=0.d0
            dmef(mi)=0.d0
            dmwf(mi)=0.d0
            dmycf(mi)=0.d0
            dmyaf(mi)=0.d0
            depf(mi)=0.d0
            deef(mi)=0.d0
            dewf(mi)=0.d0
            deycf(mi)=0.d0
            deyaf(mi)=0.d0

            if(icesd.eq.2) then
               rll = rl_v(mi)
               drllw = drl_vw(mi)
               drllg = drl_vg(mi)
               drllp = 0.0
               drllt = 0.0
            endif
c     
c     The following represents properties for co2-rich phase. It is assumed that
c     if 2 phase CO2 is present the dissolved water is in CO2 gas phase.
c     
            if((icesd.eq.1).or.(icesd.eq.4))then 
               rol = co2_prop(mi)
               drolt=co2_prop(neq+mi)
               drolp=co2_prop(2*neq+mi)
               drolw=0.d0
               drolya=0.d0
               drolyc=0.d0
               enl=co2_prop(3*neq+mi)
               denlt=co2_prop(4*neq+mi)
               denlp=co2_prop(5*neq+mi)
               denlw=0.d0
               denlya=0.d0
               denlyc=0.d0
               visl=co2_prop(6*neq+mi) + vis_tol
               dvislt=co2_prop(7*neq+mi)
               dvislp=co2_prop(8*neq+mi)
               dvislya=0.d0
               dvislyc=0.d0
               dvislw=0.d0
            else
               rol = co2_prop(9*neq+mi)
               drolt=co2_prop(10*neq+mi)
               drolp=co2_prop(11*neq+mi)
               drolw=0.d0
               drolya=0.d0
               drolyc=0.d0
               enl=co2_prop(12*neq+mi)
               denlt=co2_prop(13*neq+mi)
               denlp=co2_prop(14*neq+mi)
               denlw=0.d0
               denlya=0.d0
               denlyc=0.d0
               visl=co2_prop(15*neq+mi)+ vis_tol
               dvislt=co2_prop(16*neq+mi)
               dvislp=co2_prop(17*neq+mi)
               dvislya=0.d0
               dvislyc=0.d0
               dvislw=0.d0
            endif				
            row = wat_prop(mi)
            drowp=wat_prop(neq+mi)
            drowt=wat_prop(2*neq+mi)
            drowyc=wat_prop(3*neq+mi)
            drowya=wat_prop(4*neq+mi)	
            enw=wat_prop(5*neq+mi)
            denwp=wat_prop(6*neq+mi)
            denwt=wat_prop(7*neq+mi)
            visw=wat_prop(8*neq+mi) + vis_tol
            dviswp=wat_prop(9*neq+mi)
            dviswt=wat_prop(10*neq+mi)
            denwyc=wat_prop(11*neq+mi)
            denwya=wat_prop(12*neq+mi)
            denwfg=wat_prop(13*neq+mi)
            denwfw=wat_prop(14*neq+mi)

c     source terms and its derivatives
c     sk : water mass production rate, originally declared in comdi, old variable
c     dq : derivative of sk wrt P, originally declared in comci, old variable
c     dqt: derivative of sk wrt T, originally declared in comci, old variable
c     dqw: derivative of sk wrt fw, originally declared in comco2, new variable
c     dqxc: derivative of sk wrt xco2, originally declared in comco2, new variable
c     dqyc: derivative of sk wrt yco2, originally declared in comco2, new variable
c     dqya: derivative of sk wrt ya, originally declared in comco2, new variable   
c     qh : water energy production rate, originally declared in comdi, old variable
c     dqh : derivative of qh wrt P, originally declared in comci, old variable
c     deqh: derivative of qh wrt T, originally declared in comci, old variable
c     dqhw: derivative of qh wrt fw, originally declared in comco2, new variable
c     dqhxc: derivative of qh wrt xco2, originally declared in comco2, new variable
c     dqhyc: derivative of qh wrt yco2, originally declared in comco2, new variable
c     dqhya: derivative of qh wrt ya, originally declared in comco2, new variable   
c     the algorithm used for calculating source/sink terms depends on whether
c     wellbore model is used or not.
            
c     if(iriver.ne.1) then
            dq(mi)=0.0
            dqt(mi)=0.0
            dqw(mi) = 0.0
            dqyc(mi)=0.0
            dqya(mi)=0.0
            qh(mi)=0.0
            dqh(mi)=0.0
            deqh(mi)=0.0
            dqhw(mi)=0.0
            dqhyc(mi)=0.0
            dqhya(mi)=0.0
            
            qdis=sk(mi)
            kq=ka(mi)
            if(kq.lt.0) then
               qh(mi)=0.
               sk(mi)=0.0
               permsd=wellim(mi)
c     form pressure dependent flow term
c     RJP 02/13/07 using the formulation used by Coats, Thomas & Pierson, 
c     SPE 29111, 1995.
c     RJP 4/7/07 It is assumed that for CO2-water problem pflow(mi) will always
c     be greater than zero.
               pldif=pl-pflow(mi)
c     if(kaco2(mi).lt.0) then
c         if((kaco2(mi).eq.-1).or.(kaco2(mi).eq.-4).or.
c    &		(kaco2(mi).eq.-2)) then
c               permsd1=permsd
c                dprmp = 0.d0
c               dprmt = 0.d0
c               dprmw = 0.d0
c               dprmyc = 0.d0
c               dprmya = 0.d0
c           else
               permsd11=permsd*rlw*row*ywat/visw
               dprmp1=permsd*rlw*ywat*(drowp/visw-row*dviswp/visw*visw)
               dprmt1=permsd*rlw*ywat*(drowt/visw-row*dviswt/
     &              visw*visw)
               if(icesd.eq.2) dprmt1=permsd*drlwg*ywat*row/visw
               dprmw1=permsd*drlww*ywat*row/visw

               permsd12=permsd*rll*rol*xwat/visl
               dprmp2=permsd*rll*xwat*(drolp/visl-rol*dvislp/visl*visl)
               dprmt2=permsd*rll*xwat*(drolt/visl-rol*dvislt/
     &              visl*visl)
               if(icesd.eq.2) dprmt2=permsd*drllg*xwat*rol/visl
               dprmw2=permsd*drllw*xwat*rol/visl

               dprmyc=permsd*(rlw*(-row+ywat*drowyc)/visw)
               dprmya=permsd*(rlw*(-row+ywat*drowya)/visw)

               if(pldif.le.0.0d00) then
                  if(kq.eq.-2) then
                     permsd1=0.
                     dprmp = 0.d0
                     dprmt = 0.d0
                     dprmw = 0.d0
                     dprmyc = 0.d0
                     dprmya = 0.d0
                  endif
               endif	
               if(idof_co2.eq.1) then
                  permsd1=wellim(mi)
                  dprmp = 0.d0
                  dprmt = 0.d0
                  dprmw = 0.d0
                  dprmyc = 0.d0
                  dprmya = 0.d0
               endif
			 permsd1=permsd11+permsd12
			 dprmp=dprmp1+dprmp2
			 dprmt=dprmt1+dprmt2
			 dprmw=dprmw1+dprmw2

               qdis=permsd1*pldif
               sk(mi)=qdis
               dq(mi)=permsd1+dprmp*pldif
               dqt(mi)=dprmt*pldif
               dqw(mi)=dprmw*pldif
               dqyc(mi)=dprmyc*pldif
               dqya(mi)=dprmya*pldif
               if(iprtype.ge.4) then
                  if(ico2dis(mi).eq.0) then
                     dqw(mi)=dqyc(mi)
                  else
                     if(icesd.eq.2) then
                        dq(mi)=dq(mi)+dqyc(mi)*dmol(mi)
     &                       +dqyc(mi)*dmol(mi+neq)
                     else
                        dq(mi)=dq(mi)+dqyc(mi)*dmol(mi)
                        dqt(mi)=dqt(mi)+dqyc(mi)*dmol(mi+neq)
                     endif
                  endif
               endif

			if(pldif.gt.0.d0) then
				qh(mi)=(permsd11*enw+permsd12*enl)*pldif
				dqh(mi)=(dprmp1*enw+permsd11*denwp+dprmp2*enl+
     &					permsd12*denlp)*pldif+permsd11*enw+
     &					permsd12*enl
				deqh(mi)=(dprmt1*enw+permsd11*denwt+dprmt2*enl+
     &					permsd12*denlt)*pldif
				dqhw(mi)=(dprmw1*enw+dprmw2*enl)*pldif
				dqhyc(mi)=(dprmyc*enw+permsd11*denwyc)*pldif
				dqhya(mi)=(dprmya*enw+permsd11*denwya)*pldif
				if(iprtype.ge.4) then
					if(ico2dis(mi).eq.0) then
						dqhw(mi)=dqhyc(mi)
					else
						if(icesd.eq.2) then
							dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
     &                         +dqhyc(mi)*dmol(mi+neq)
	                    else
							dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
							deqh(mi)=deqh(mi)+dqhyc(mi)*dmol(mi+neq)
						endif
					endif
                  endif
			endif
            endif
c     endif

c     compressed liquid
            eqdum=rol*frac_c*xwat+row*frac_w*ywat
            den=por*eqdum
            damp=dtin*(dporpl*eqdum+
     &           por*(drolp*frac_c*xwat+drowp*frac_w*ywat))
            damh=dtin*(dportl*eqdum+
     &           por*(drolt*frac_c*xwat+drowt*frac_w*ywat))
            damw=por*dtin*(-rol*xwat+row*ywat)
            if(icesd.eq.2) damh=por*dtin*rol*xwat
            if(icesd.eq.2) damw=por*dtin*row*ywat
            damyc=por*dtin*frac_w*(drowyc*ywat-row)
            damya=por*dtin*frac_w*(drowya*ywat-row)
            eqdum=frac_c*(rol*enl-pl)*xwat+
     &           frac_w*(row*enw-pl)*ywat
            dene=por*eqdum
            daep=dtin*(dporpl*eqdum+
     &           por*(frac_c*xwat*(drolp*enl+rol*denlp-1)+
     &           frac_w*ywat*(drowp*enw+row*denwp-1)))
            daeh=dtin*(dportl*eqdum+
     &           por*(frac_c*xwat*(drolt*enl+rol*denlt)+
     &           frac_w*ywat*(drowt*enw+row*denwt)))
            if(icesd.eq.2) daeh=por*dtin*(xwat*(rol*enl-pl)+
     &           frac_w*row*denwfg*ywat)
            if(iprtype.eq.1) then
c     RJP For water-only problem, add heat capacity of rock
               dene = dene+(1.-por)*cp*tl
               daep2=-dporpl*cp*dtin
               daep=daep+daep2
               daeh=daeh+(1.d0-por)*cp*dtin
            endif
            daew=por*dtin*(-xwat*(rol*enl-pl)+
     &           ywat*(row*enw-pl))
            if(icesd.eq.2) daew=por*dtin*ywat*((row*enw-pl)+
     &           frac_w*row*denwfw)
c     daeyc=por*dtin*frac_w*(ywat*(drowyc*enw+row*denwyc)-(row*enw-pl))
            daeyc=por*dtin*frac_w*(ywat*(drowyc*enw+row*denwyc)
     &           -(row*enw-pl))
            daeya=por*dtin*frac_w*(ywat*(drowya*enw+row*denwya)
     &           -(row*enw-pl))
            if(ps(mi).ne.0.0) then          
c     store derivatives of accumulation terms
               dmpf(mi)=damp
               dmef(mi)=damh
               dmwf(mi)=damw
               dmycf(mi)=damyc
               dmyaf(mi)=damya
               depf(mi)=daep
               deef(mi)=daeh
               dewf(mi)=daew
               deycf(mi)=daeyc
               deyaf(mi)=daeya
               if(iprtype.ge.4) then
                  if(ico2dis(mi).eq.0) then
                     dmwf(mi)=damyc
                     dewf(mi)=daeyc
                  else
                     if(icesd.eq.2) then
                        dmpf(mi)=dmpf(mi)+damyc*dmol(mi)+
     &                       damyc*dmol(mi+neq)
                        depf(mi)=depf(mi)+daeyc*dmol(mi)+
     &                       daeyc*dmol(mi+neq)
                     else
                        dmpf(mi)=dmpf(mi)+damyc*dmol(mi)
                        dmef(mi)=dmef(mi)+damyc*dmol(mi+neq)
                        depf(mi)=depf(mi)+daeyc*dmol(mi)
                        deef(mi)=deef(mi)+daeyc*dmol(mi+neq)
                     endif
                  endif
               endif
               t(mi)=tl
c     save accumulation terms for possible volume changes
               sto1(mi)=den
               sto1(mi+neq)=dene
c     zvd 07-Aug-08 added following line
               rov = co2_prop(9*neq+mi)
c     zvd 07-Aug-08
               dstm(mi)=por*rov*sv*vol
               diw(mi)=0.0
               diwp(mi)=0.0
               diwe(mi)=0.0
               diww(mi)=0.0
               diwyc(mi)=0.0
               diwya(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divw(mi)=0.0
               divyc(mi)=0.0
               divya(mi)=0.0
c     
c     modify flow terms for new upwind scheme
c     .First calculate transfer through water-rich phase
c     RJP 03/26/07 introduced 'w' to refer to water-rich
c     phase. 'l' will now refer to CO2-rich (liquid) phase 
c     for 2 phase CO2 (3-phase problem) and v refers to CO2-rich
c     (gas/liquid/SC) phase or air for single phase CO2 (2-phase
c     problem) or CO2-rich gas phase for 2 phase CO2 (3-phase 
c     problem). 
c     
               dql=rlw*row*ywat/visw
               diw(mi)=dql
               diwp(mi)=rlw*(drowp*ywat/visw-
     &              row*ywat*dviswp/visw**2)
               diwe(mi)=rlw*(drowt*ywat/visw-
     &              row*ywat*dviswt/visw**2)
               if(icesd.eq.2) diwe(mi)=drlwg*row*ywat/visw
               diww(mi)=drlww*row*ywat/visw
               diwya(mi)=rlw*(-row+drowya*ywat)/visw
               diwyc(mi)=rlw*(-row+drowyc*ywat)/visw
c     
c     Next calculate transfer through CO2-rich phase
c     
               dqv=rll*rol*xwat/visl
               div(mi)=dqv
               divp(mi)=rll*(drolp*xwat/visl-
     &              rov*xwat*dvislp/visl**2)
               dive(mi)=rll*(drolt*xwat/visl-
     &              rov*xwat*dvislt/visl**2)
               if(icesd.eq.2) dive(mi)=drllg*rol*xwat/visl
               divw(mi)=drllw*rol*xwat/visl
               divya(mi)=0.d0
               divyc(mi)=0.d0

               if(iprtype.ge.4) then
                  if(ico2dis(mi).eq.0) then
                     diww(mi)=diwyc(mi)
                     divw(mi)=divyc(mi)
                  else
                     if(icesd.eq.2) then
                        diwp(mi)=diwp(mi)+diwyc(mi)*dmol(mi)
     &                       +diwyc(mi)*dmol(mi+neq)
                        divp(mi)=divp(mi)+divyc(mi)*dmol(mi)
     &                       +divyc(mi)*dmol(mi+neq)
                     else
                        diwp(mi)=diwp(mi)+diwyc(mi)*dmol(mi)
                        diwe(mi)=diwe(mi)+diwyc(mi)*dmol(mi+neq)
                        divp(mi)=divp(mi)+divyc(mi)*dmol(mi)
                        dive(mi)=dive(mi)+divyc(mi)*dmol(mi+neq)
                     endif
                  endif
               endif

               if(qdis.le.0.) then
                  eskd=eflow(mi)
                  qh(mi)=qdis*eflow(mi)
                  if((kaco2(mi).ne.2).and.(kaco2(mi).ne.3)) then
                     dqh(mi)=eskd*dq(mi)
                     deqh(mi)=eskd*dqt(mi)
                     dqhw(mi)=eskd*dqw(mi)
                     dqhyc(mi)=eskd*dqyc(mi)
                     dqhya(mi)=eskd*dqya(mi)
                     if(iprtype.ge.4) then
                        if(ico2dis(mi).eq.0) then
                           dqhw(mi)=dqhyc(mi)
                        else
                           if(icesd.eq.2) then
                              dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
     &                          +dqhyc(mi)*dmol(mi+neq)
                           else
                              dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
                              deqh(mi)=deqh(mi)+dqhyc(mi)*dmol(mi+neq)
                           endif
                        endif
                     endif
                  endif
               else
                  qh(mi) = qdis*enw
                  dqh(mi) = qdis*denwp + enw*dq(mi)
                  deqh(mi) = qdis*denwt + enw*dqt(mi)
                  dqhw(mi) = 0.d0
                  dqhyc(mi) = 0.d0
                  dqhya(mi) = 0.d0
               endif
            endif
c     
c     add intercomponent heat flux
c     this term transfers heat from water to co2 (and back)
c     also note this term only contains derivative wrt water   
c     derivative with respect to co2 can also be derived 
c     
            if(qhflxco2(mi).ne.0.0) then
               htc = qhflxco2(mi)
               hflux=htc*(t(mi)-tco2(mi))
               if(abs(ieosd).ne.2.and.abs(icesd).ne.2) then
                  dhflxp=0.0
                  dhflxe=htc
               else
                  dhflxp=htc*dtpa(mi)
                  dhflxe=0.0
               endif
c     this next if block contains derivative wrt co2 component
               if(abs(ices(mi)).ne.2) then
                  dhflxew=-htc
               else
                  dhflxpm=-htc*dtpaco2(mi)
                  dhflxem=0.0
               endif
               qh(mi) = qh(mi) + hflux
               dqh(mi) = dqh(mi) + dhflxp
               deqh(mi) = deqh(mi) + dhflxe
               deqpm(mi) = deqpm(mi) + dhflxpm
               deqm(mi) = deqm(mi) + dhflxem
            endif
c     
c     add heat source term
c     
            if(qflux(mi).ne.0.0) then
               if(qflxm(mi).gt.0.0) then
                  htc=qflxm(mi)
                  tbound=qflux(mi)
                  hflux=htc*(tl-tbound)
                  if(abs(ieosd).ne.2.and.abs(icesd).ne.2) then
                     dhflxp=0.0
                     dhflxe=htc
                  else
                     dhflxp=htc*dtps
                     dhflxe=0.0
                  endif
                  qh(mi)=qh(mi)+hflux
                  dqh(mi)=dqh(mi)+dhflxp
                  deqh(mi)=deqh(mi)+dhflxe
               else if(qflxm(mi).lt.0.0) then
                  htc=abs(qflxm(mi))
                  sbound=qflux(mi)
                  hflux=htc*(tl-sbound)
                  if(abs(ieosd).ne.2.and.abs(icesd).ne.2) then
                     dhflxp=0.0
                     dhflxe=0.0
                  else
                     dhflxp=0.0
                     dhflxe=htc
                  endif
                  qh(mi)=qh(mi)+hflux
                  dqh(mi)=dqh(mi)+dhflxp
                  deqh(mi)=deqh(mi)+dhflxe
               else
                  qh(mi)=qh(mi)+qflux(mi)
               endif
            endif

            if(ps(mi).eq.0.0) then
c     if(idof_co2.eq.1) then
c     heat conduction only
               cprd=cpr(mi)
c     if(kq.lt.0) then
c     eskd=eflow(mi)
c     edif=cprd*(tl-eskd)
c     permsd=wellim(mi)
c     qh(mi)=permsd*edif
c     deqh(mi)=permsd*cprd
c     dqhw(mi)=0.d0
c     dqhyc(mi)=0.d0
c     dqhya(mi)=0.d0
c     endif
               if(kq.ge.0.and.qflux(mi).eq.0.0) deqh(mi)=0.
               dtpae(mi)=1.
               denrd=denr(mi)*cprd
               sto1(mi)=0.0
               sto1(mi+neq)=denrd*tl
               deef(mi)=denrd*dtin
c     endif
            endif
         enddo
         
         do     mid=1,neq
            mi=mid+ndummy
            deni(mi)=(sto1(mi)-denh(mi))*dtin
            denei(mi)=(sto1(mi+neq)-deneh(mi))*dtin
         enddo    
c     Hari store diw in diwc for tracer solution   
         diwc = diw 
         
c     ****************************************************************
      else if(iflg.eq.2) then
c     
c     calculations for co2
c     
         do mid=1,neq
            mi=mid+ndummy
c     phase information now contained in array ices
            icesd=ices(mi)
            pl=phico2(mi)
            tl=tco2(mi)
               diw(mi)=0.0
               diwp(mi)=0.0
               diwe(mi)=0.0
               diww(mi)=0.0
               diwyc(mi)=0.0
               diwya(mi)=0.0
               dil(mi)=0.0
               dilp(mi)=0.0
               dile(mi)=0.0
               dilw(mi)=0.0
               dilyc(mi)=0.0
               dilya(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divw(mi)=0.0
               divyc(mi)=0.0
               divya(mi)=0.0           
            dtps=dtpsc(mi)
            frac_w= fw(mi)
            frac_cl=fl(mi)
            frac_cg=fg(mi)
            yco2 = yc(mi)
            ywat = yw(mi)
            yair = ya(mi)
            xco2 = xc(mi)
            xwat = xw(mi)
            xair = xa(mi)
            cp=denr(mi)*cpr(mi)
            por=ps(mi)
            dporpl=dporp(mi)
            dportl=dport(mi)
            vol=volume(mi)
            dtps=dtpsc(mi)
c     rel-perms  
            rlw = rl_w(mi)
            drlww = drl_ww(mi)
            drlwg = drl_wg(mi)
            drlwp = 0.0
            drlwt = 0.0
            rll = rl_l(mi)
            drllw = drl_lw(mi)
            drllg = drl_lg(mi)
            drllp = 0.0
            drllt = 0.0
            rlv = rl_v(mi)
            drlvw = drl_vw(mi)
            drlvg = drl_vg(mi)
            drlvp = 0.0
            drlvt = 0.0

            dmpf(mi)=0.d0
            dmef(mi)=0.d0
            dmwf(mi)=0.d0
            dmycf(mi)=0.d0
            dmyaf(mi)=0.d0
            depf(mi)=0.d0
            deef(mi)=0.d0
            dewf(mi)=0.d0
            deycf(mi)=0.d0
            deyaf(mi)=0.d0

c     
c     The following represents properties for co2-rich phase. It is assumed that
c     if 2 phase CO2 is present the dissolved water is in CO2 gas phase.
c     
            rol = co2_prop(mi)
            drolt=co2_prop(neq+mi)
            drolp=co2_prop(2*neq+mi)
            drolw=0.d0
            drolya=0.d0
            drolyc=0.d0
            enl=co2_prop(3*neq+mi)
            denlt=co2_prop(4*neq+mi)
            denlp=co2_prop(5*neq+mi)
            denlw=0.d0
            denlya=0.d0
            denlyc=0.d0
            visl=co2_prop(6*neq+mi)
            dvislt=co2_prop(7*neq+mi)
            dvislp=co2_prop(8*neq+mi)
            dvislya=0.d0
            dvislyc=0.d0
            dvislw=0.d0

            rov = co2_prop(9*neq+mi)
            drovt=co2_prop(10*neq+mi)
            drovp=co2_prop(11*neq+mi)
            drovw=0.d0
            drovya=0.d0
            drovyc=0.d0
            env=co2_prop(12*neq+mi)
            denvt=co2_prop(13*neq+mi)
            denvp=co2_prop(14*neq+mi)
            denvw=0.d0
            denvya=0.d0
            denvyc=0.d0
            visv=co2_prop(15*neq+mi) + vis_tol
            dvisvt=co2_prop(16*neq+mi)
            dvisvp=co2_prop(17*neq+mi)
            dvisvya=0.d0
            dvisvyc=0.d0
            dvisvw=0.d0
            
            row = wat_prop(mi)
            drowp=wat_prop(neq+mi)
            drowt=wat_prop(2*neq+mi)
            drowyc=wat_prop(3*neq+mi)
            drowya=wat_prop(4*neq+mi)	
            enw=wat_prop(5*neq+mi)
            denwp=wat_prop(6*neq+mi)
            denwt=wat_prop(7*neq+mi)
            visw=wat_prop(8*neq+mi) + vis_tol
            dviswp=wat_prop(9*neq+mi)
            dviswt=wat_prop(10*neq+mi)
            denwyc=wat_prop(11*neq+mi)
            denwya=wat_prop(12*neq+mi)
            denwfg=wat_prop(13*neq+mi)
            denwfw=wat_prop(14*neq+mi)

c     if(iriver.ne.1) then
            dq(mi)=0.0
            dqt(mi)=0.0
            dqw(mi)=0.0
            dqyc(mi)=0.0
            dqya(mi)=0.0
            qhco2(mi)=0.0
            dqh(mi)=0.0
            deqh(mi)=0.0
            dqhw(mi)=0.0
            dqhyc(mi)=0.0
            dqhya(mi)=0.0
            
            kq=kaco2(mi)

            if(kq.lt.0) then
c     form pressure dependent flow term
c     gaz 2-19-04 added gas relative perms
               qhco2(mi)=0.
               skco2(mi)=0.0
               permsd=wellco2(mi)
c
c gaz debug 110713
c pflowco2 should equal pflow
c               pldif=pl-pflowco2(mi)
c
               pldif=pl-pflow(mi)
               if((kq.eq.-4)) then
                  permsd1=permsd
                  dprmp=0.d0
                  dprmt=0.d0
                  dprmw=0.d0
                  dprmyc=0.d0
                  dprmya = 0.d0
               else
                  if(icesd.ne.2) then
                     if(icesd.eq.3) then
                        permsd11=permsd*rlw*row*yco2/visw
                        permsd12=permsd*rll*rov*xco2/visv

                        dprmp1=permsd*rlw*yco2*(drowp/visw-row*dviswp/
     &                       visw*visw)
                        dprmp2=permsd*rll*xco2*(drovp/visv-rov*dvisvp/
     &                       visv*visv)

                        dprmt1=permsd*rlw*yco2*(drowt/visw-row*dviswt/
     &                       visw*visw)
                        dprmt2=permsd*rll*xco2*(drovt/visv-rov*dvisvt/
     &                       visv*visv)

                        dprmw1=permsd*drlww*row*yco2/visw
                        dprmw2=permsd*drllw*rov*xco2/visv

                        dprmyc=permsd*rlw*(row+yco2*drowyc)/visw
                        dprmya=0.0

					  permsd1=permsd11+permsd12
					  dprmp=dprmp1+dprmp2
					  dprmt=dprmt1+dprmt2
					  dprmw=dprmw1+dprmw2

                     else
                        permsd11=permsd*rlw*row*yco2/visw
                        permsd12=permsd*rll*rol*xco2/visl

                        dprmp1=permsd*rlw*yco2*(drowp/visw-row*dviswp/
     &                       visw*visw)
                        dprmp2=permsd*rll*xco2*(drolp/visl-rol*dvislp/
     &                       visl*visl)
                        dprmt1=permsd*rlw*yco2*(drowt/visw-row*dviswt/
     &                       visw*visw)
                        dprmt2=permsd*rlw*yco2*(drowt/visw-row*dviswt/
     &                       visw*visw)

                        dprmw1=permsd*drlww*row*yco2/visw
                        dprmw2=permsd*drllw*rol*xco2/visl

                        dprmyc=permsd*rlw*(row+yco2*drowyc)/visw
                        dprmya=0.0

					  permsd1=permsd11+permsd12
					  dprmp=dprmp1+dprmp2
					  dprmt=dprmt1+dprmt2
					  dprmw=dprmw1+dprmw2

                     endif
                  else
                     permsd11=permsd*rlw*row*yco2/visw
                     permsd12=permsd*rll*rol/visl
                     permsd13=permsd*rlv*rov*xco2/visv

                     dprmp1=permsd*rlw*yco2*(drowp/visw-row*dviswp/
     &                    visw*visw)
                     dprmp2=permsd*rll*(drolp/visl-rol*dvislp/
     &                    visl*visl)
                     dprmp3=permsd*rlv*xco2*(drovp/visv-rov*dvisvp/
     &                    visv*visv)

                     dprmt1=permsd*drlwg*yco2*row/visw
                     dprmt2=permsd*drllg*rol/visl
                     dprmt3=permsd*drlvg*xco2*rov/visv

                     dprmw1=permsd*drlww*row*yco2/visw
                     dprmw2=permsd*drllw*rol/visl
                     dprmw3=permsd*drlvw*rov*xco2/visv

                     dprmyc=permsd*rlw*(row+yco2*drowyc)/visw
                     dprmya=0.0

				   permsd1=permsd11+permsd12+permsd13
            dprmp=dprmp1+dprmp2+dprmp3
				   dprmt=dprmt1+dprmt2+dprmt3
				   dprmw=dprmw1+dprmw2+dprmw3

                  endif
               endif

               if((pldif.le.0.0d00).and.(kq.eq.-2)) then
                  permsd1=0.d0
                  dprmp=0.d0
                  dprmt=0.d0
                  dprmw=0.d0
                  dprmya=0.d0
                  dprmyc=0.d0
               endif
c gaz 09-02-2010               
               if((pldif.le.0.0d00).and.(kq.eq.-1)) then
                  permsd1=wellco2(mi)
                  dprmp=0.d0
                  dprmt=0.d0
                  dprmw=0.d0
                  dprmya=0.d0
                  dprmyc=0.d0
                  qdis=permsd1*pldif
                  skco2(mi)=qdis
                  dq(mi)=permsd1+dprmp*pldif
               endif               
               if(idof_co2.eq.2) then
                  permsd1=wellco2(mi)
                  dprmp = 0.d0
                  dprmt = 0.d0
                  dprmw = 0.d0
                  dprmyc = 0.d0
                  dprmya = 0.d0
               endif

               qdis=permsd1*pldif
               skco2(mi)=qdis
               dq(mi)=permsd1+dprmp*pldif
               dqt(mi)=dprmt*pldif
               dqw(mi)=dprmw*pldif
               dqyc(mi)=dprmyc*pldif
               dqya(mi)=dprmya*pldif
               if((kq.eq.-4).or.(kq.eq.-5)) then
                  if(icesd.eq.3) then					
                     skco2(mi) = skco2(mi)+permsd1*(fg(mi)-
     &                    flowco2s(mi))
                     dqw(mi) = dqw(mi)-permsd1
                     dqyc(mi) = dqyc(mi)+permsd1*pldif
                  elseif(icesd.eq.2) then
                     skco2(mi) = skco2(mi)+permsd1*(fg(mi)-
     &                    flowco2s(mi))
                     dqt(mi) = dqt(mi)+permsd1
                     dqyc(mi) = dqyc(mi)+permsd1*pldif
                  else
                     skco2(mi) = skco2(mi)+permsd1*(fl(mi)-
     &                    flowco2s(mi))
                     dqw(mi) = dqw(mi)-permsd1
                     dqyc(mi) = dqyc(mi)+permsd1*pldif
                  endif
               endif
               if(iprtype.ge.4) then
                  if(ico2dis(mi).eq.0) then
                     dqw(mi)=dqyc(mi)
                  else
                     if(icesd.eq.2) then
                        dq(mi)=dq(mi)+dqyc(mi)*dmol(mi)
     &                       +dqyc(mi)*dmol(mi+neq)
                     else
                        dq(mi)=dq(mi)+dqyc(mi)*dmol(mi)
                        dqt(mi)=dqt(mi)+dqyc(mi)*dmol(mi+neq)
                     endif
                  endif
               endif

			if(pldif.gt.0.d0) then
				if(icesd.eq.3) then
				qhco2(mi)=(permsd11*enw+permsd12*env)*pldif
				dqh(mi)=(dprmp1*enw+permsd11*denwp+dprmp2*env+
     &					permsd12*denvp)*pldif+permsd11*enw+
     &					permsd12*env
				deqh(mi)=(dprmt1*enw+permsd11*denwt+dprmt2*env+
     &					permsd12*denvt)*pldif
				dqhw(mi)=(dprmw1*enw+dprmw2*env)*pldif
				dqhyc(mi)=(dprmyc*enw+permsd11*denwyc)*pldif
				dqhya(mi)=(dprmya*enw+permsd11*denwya)*pldif
				elseif(icesd.eq.2) then
				qhco2(mi)=(permsd11*enw+permsd12*enl+permsd13*env)
     &			*pldif
				dqh(mi)=(dprmp1*enw+permsd11*denwp+dprmp2*enl+
     &					permsd12*denlp+dprmp3*env+permsd13*denvp)
     &					*pldif+permsd11*enw+permsd12*enl+permsd13*env
				deqh(mi)=(dprmt1*enw+permsd11*denwt+dprmt2*enl+
     &					permsd12*denlt+dprmt3*env+permsd13*denvt)
     &					*pldif+permsd11*enw+permsd12*enl+permsd13*env
				dqhw(mi)=(dprmw1*enw+dprmw2*enl+dprmw3*env)*pldif
				dqhyc(mi)=(dprmyc*enw+permsd11*denwyc)*pldif
				dqhya(mi)=(dprmya*enw+permsd11*denwya)*pldif
				else
				qhco2(mi)=(permsd11*enw+permsd12*enl)*pldif
				dqh(mi)=(dprmp1*enw+permsd11*denwp+dprmp2*enl+
     &					permsd12*denlp)*pldif+permsd11*enw+
     &					permsd12*enl
				deqh(mi)=(dprmt1*enw+permsd11*denwt+dprmt2*enl+
     &					permsd12*denlt)*pldif
				dqhw(mi)=(dprmw1*enw+dprmw2*enl)*pldif
				dqhyc(mi)=(dprmyc*enw+permsd11*denwyc)*pldif
				dqhya(mi)=(dprmya*enw+permsd11*denwya)*pldif
				endif

				if(iprtype.ge.4) then
					if(ico2dis(mi).eq.0) then
						dqhw(mi)=dqhyc(mi)
					else
						if(icesd.eq.2) then
							dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
     &                         +dqhyc(mi)*dmol(mi+neq)
	                    else
							dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
							deqh(mi)=deqh(mi)+dqhyc(mi)*dmol(mi+neq)
						endif
					endif
                  endif
			endif

            endif
c     identify flux term
c     endif
            qdis=skco2(mi)
c     accumulation terms
            if(icesd.eq.2) then
c     two phase condition
               eqdum=row*frac_w*yco2+frac_cg*rov*xco2+frac_cl*rol
               den=por*eqdum
               damp=dtin*(dporpl*eqdum+por*(frac_w*yco2*drowp+
     &              frac_cg*xco2*drovp+frac_cl*drolp))
               damh=dtin*por*(xco2*rov-rol)
               damw=dtin*por*(yco2*row-rol)
               damyc=dtin*por*frac_w*(row+yco2*drowyc)
               damya=0.0
               eqdum=frac_w*(row*enw-pl)*yco2+frac_cg*(rov*env-pl)*
     &              xco2+frac_cl*(rol*enl-pl)
               dene=por*eqdum+(1.-por)*cp*tl
               daep1=por*dtin*(frac_w*yco2*(drowp*enw+row*denwp-1.d0)
     &              +frac_cg*xco2*(drovp*env+rov*denvp-1.d0)+frac_cl*
     &              (drolp*enl+rol*denlp-1.d0))
               daep2=(-dporpl*cp*tl+(1-por)*cp*dtps)*dtin
               daep=daep1+daep2
               daeh=dtin*por*((xco2*(rov*env-pl)-(rol*enl-pl))+
     &              yco2*frac_w*row*denwfg)
               daew=dtin*por*(yco2*(frac_w*row*denwfw+(row*enw-pl))
     &              -(rol*enl-pl))
               daeyc=dtin*por*frac_w*((row*enw-pl)+yco2*(drowyc*enw+
     &              row*denwyc))
               daeya=0.0
            elseif(icesd.eq.3) then
               eqdum=row*frac_w*yco2+frac_cg*rov*xco2
               den=por*eqdum
               damp=dtin*(dporpl*eqdum+por*(frac_w*yco2*drowp+
     &              frac_cg*xco2*drovp))
               damh=dtin*(dportl*eqdum+por*(frac_w*yco2*drowt+
     &              frac_cg*xco2*drovt))
               damw=dtin*por*(row*yco2-rov*xco2)
               damyc=dtin*por*frac_w*(row+yco2*drowyc)
               damya=0.0
               eqdum=frac_w*(row*enw-pl)*yco2+
     &              frac_cg*(rov*env-pl)*xco2
               dene=por*eqdum+(1.-por)*cp*tl
               daep1=por*dtin*(frac_w*yco2*(drowp*enw+row*denwp-1.d0)
     &              +frac_cg*xco2*(drovp*env+rov*denvp-1.d0))
               daep2=-dporpl*cp*tl*dtin
               daep=daep1+daep2
               daeh1=dtin*por*(frac_w*(drowt*enw+row*denwt)*yco2+
     &              frac_cg*xco2*(drovt*env+rov*denvt))
               daeh2=(1.d0-por)*cp*dtin
               daeh=daeh1+daeh2
               daew=dtin*por*((row*enw-pl)*yco2-xco2*(rov*env-pl))
               daeyc=dtin*por*frac_w*((row*enw-pl)+yco2*(drowyc*enw+
     &              row*denwyc))
               daeya=0.0
            else
               eqdum=row*frac_w*yco2+frac_cl*rol*xco2
               den=por*eqdum
               damp=dtin*(dporpl*eqdum+por*(frac_w*yco2*drowp+
     &              frac_cl*xco2*drolp))
               damh=dtin*(dportl*eqdum+por*(frac_w*yco2*drowt+
     &              frac_cl*xco2*drolt))
               damw=dtin*por*(row*yco2-rol*xco2)
               damyc=dtin*por*frac_w*(row+yco2*drowyc)
               damya=0.0
               eqdum=frac_w*(row*enw-pl)*yco2+
     &              frac_cl*(rol*enl-pl)*xco2
               dene=por*eqdum+(1.-por)*cp*tl
               daep1=por*dtin*(frac_w*yco2*(drowp*enw+row*denwp-1.d0)
     &              +frac_cl*xco2*(drolp*enl+rol*denlp-1.d0))
               daep2=-dporpl*cp*tl*dtin
               daep=daep1+daep2
               daeh1=dtin*por*(frac_w*(drowt*enw+row*denwt)*yco2+
     &              frac_cl*xco2*(drolt*enl+rol*denlt))
               daeh2=(1.d0-por)*cp*dtin
               daeh=daeh1+daeh2
               daew=dtin*por*((row*enw-pl)*yco2-xco2*(rol*enl-pl))
               daeyc=dtin*por*frac_w*((row*enw-pl)+yco2*(drowyc*enw+
     &              row*denwyc))
               daeya=0.0
            endif
            if(ps(mi).ne.0.0) then
c     derivatives of accumulation terms 
c     store derivatives of accumulation terms
               dmef(mi)=damh
               dmpf(mi)=damp
               dmwf(mi)=damw
               dmycf(mi)=damyc
               dmyaf(mi)=damya
               depf(mi)=daep
               deef(mi)=daeh
               dewf(mi)=daew
               deycf(mi)=daeyc
               deyaf(mi)=daeya
               if(iprtype.ge.4) then
                  if(ico2dis(mi).eq.0) then
                     dmwf(mi)=damyc
                     dewf(mi)=daeyc
                  else
                     if(icesd.eq.2) then
                        dmpf(mi)=dmpf(mi)+damyc*dmol(mi)
     &                       +damyc*dmol(mi+neq)
                        depf(mi)=depf(mi)+daeyc*dmol(mi)
     &                       +daeyc*dmol(mi+neq)
                     else
                        dmpf(mi)=dmpf(mi)+damyc*dmol(mi)
                        dmef(mi)=dmef(mi)+damyc*dmol(mi+neq)
                        depf(mi)=depf(mi)+daeyc*dmol(mi)
                        deef(mi)=deef(mi)+daeyc*dmol(mi+neq)
                     endif
                  endif
               endif
c     save accumulation terms for possible volume changes
               sto1(mi)=den
               sto1(mi+neq)=dene
c     dtpaco2(mi)=dtps
c     dtpaeco2(mi)=dtd
               dstm(mi)=por*rov*sv*vol
c     RJP 03/27/07. diw notes transport through water-rich phase
c     dil notes transport through co2-rich liquid phase (only present
c     for 3-phase problem, 2-phase CO2 & water)
c     div notes transport through co2-rich vapor phase
               diw(mi)=0.0
               diwp(mi)=0.0
               diwe(mi)=0.0
               diww(mi)=0.0
               diwyc(mi)=0.0
               diwya(mi)=0.0
               dil(mi)=0.0
               dilp(mi)=0.0
               dile(mi)=0.0
               dilw(mi)=0.0
               dilyc(mi)=0.0
               dilya(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divw(mi)=0.0
               divyc(mi)=0.0
               divya(mi)=0.0
               dql=rlw*yco2*row/visw
c     transport in water-rich phase.
               diw(mi)=dql
               diwp(mi)=rlw*yco2*(drowp/visw-
     &              row*dviswp/visw**2)
               diwe(mi)=rlw*yco2*(drowt/visw-
     &              row*dviswt/visw**2)
               if(icesd.eq.2) diwe(mi)=drlwg*yco2*row/visw
               diww(mi)=drlww*yco2*row/visw
               diwyc(mi)=rlw*(row+yco2*drowyc)/visw
               diwya(mi)=0.0
c     transport in CO2-rich phase.
               if(icesd.eq.2) then
                  dqv=rlv*xco2*rov/visv
                  div(mi)=dqv
                  divp(mi)=rlv*xco2*(drovp/visv-
     &                 rov*dvisvp/visv**2)
                  dive(mi)=drlvg*xco2*rov/visv
                  divw(mi)=drlvw*xco2*rov/visv
                  divyc(mi)=0.0
                  divya(mi)=0.0

                  dqv=rll*rol/visl
                  dil(mi)=dqv
                  dilp(mi)=rll*(drolp/visl-
     &                 rol*dvislp/visl**2)
                  dile(mi)=drllg*rol/visl
                  dilw(mi)=drllw*rol/visl
                  dilyc(mi)=0.0
                  dilya(mi)=0.0
               elseif(icesd.eq.3) then
                  dqv=rll*xco2*rov/visv
                  div(mi)=dqv
                  divp(mi)=rll*xco2*(drovp/visv-
     &                 rov*dvisvp/visv**2)
                  dive(mi)=rll*xco2*(drovt/visv-
     &                 rov*dvisvt/visv**2)
                  divw(mi)=drllw*rov*xco2/visv
                  divyc(mi)=0.0
                  divya(mi)=0.0
               else
                  dqv=rll*xco2*rol/visl
                  dil(mi)=dqv
                  dilp(mi)=rll*xco2*(drolp/visl-
     &                 rol*dvislp/visl**2)
                  dile(mi)=rll*xco2*(drolt/visl-
     &                 rol*dvislt/visl**2)
                  dilw(mi)=drllw*rol*xco2/visl
                  dilyc(mi)=0.0
                  dilya(mi)=0.0
               endif

               if(iprtype.ge.4) then
                  if(ico2dis(mi).eq.0) then
                     diww(mi)=diwyc(mi)
                     dilw(mi)=dilyc(mi)
                     divw(mi)=divyc(mi)
                  else
                     if(icesd.eq.2) then
                        diwp(mi)=diwp(mi)+diwyc(mi)*dmol(mi)
     &                       +diwyc(mi)*dmol(mi+neq)
                        dilp(mi)=dilp(mi)+dilyc(mi)*dmol(mi)
     &                       +dilyc(mi)*dmol(mi+neq)
                        divp(mi)=divp(mi)+divyc(mi)*dmol(mi)
     &                       +divyc(mi)*dmol(mi+neq)
                     else
                        diwp(mi)=diwp(mi)+diwyc(mi)*dmol(mi)
                        diwe(mi)=diwe(mi)+diwyc(mi)*dmol(mi+neq)
                        dilp(mi)=dilp(mi)+dilyc(mi)*dmol(mi)
                        dile(mi)=dile(mi)+dilyc(mi)*dmol(mi+neq)
                        divp(mi)=divp(mi)+divyc(mi)*dmol(mi)
                        dive(mi)=dive(mi)+divyc(mi)*dmol(mi+neq)
                     endif
                  endif
               endif

               if(qdis.le.0.) then
                  eskd=eflowco2(mi)
                  qhco2(mi)=eskd*qdis
				if((kaco2(mi).ne.2).and.(kaco2(mi).ne.3)) then
                  dqh(mi)=eskd*dq(mi)
                  deqh(mi)=eskd*dqt(mi)
                  dqhw(mi)=eskd*dqw(mi)
                  dqhyc(mi)=eskd*dqyc(mi)
                  dqhya(mi)=eskd*dqya(mi)
                  if(iprtype.ge.4) then
                     if(ico2dis(mi).eq.0) then
                        dqhw(mi)=dqhyc(mi)
                     else
                        if(icesd.eq.2) then
                           dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
     &                          +dqhyc(mi)*dmol(mi+neq)
                        else
                           dqh(mi)=dqh(mi)+dqhyc(mi)*dmol(mi)
                           deqh(mi)=deqh(mi)+dqhyc(mi)*dmol(mi+neq)
                        endif
                     endif
                  endif							
               endif
            endif
            endif

            if(ps(mi).eq.0.0) then
c     heat conduction only
               if(idof_co2.eq.2) then
                  cprd=cpr(mi)
c     if(kq.lt.0) then
c     eskd=eflowco2(mi)
c     edif=cprd*(tl-eskd)
c     permsd=wellco2(mi)
c     qh(mi)=permsd*edif
c     deqh(mi)=permsd*cprd
c     dqhw(mi)=0.d0
c     dqhyc(mi)=0.d0
c     dqhya(mi)=0.d0
c     endif
                  if(kq.ge.0.and.qflux(mi).eq.0.0) deqh(mi)=0.
                  denrd=denr(mi)*cprd
                  sto1(mi)=0.0
                  sto1(mi+neq)=denrd*tl
                  deef(mi)=denrd*dtin
               else
                  sto1(mi)=0.d0
                  sto1(mi+neq)=0.d0
               endif
               dtpaco2(mi)=0.d0
               dtpaeco2(mi)=1.d0
            endif
         enddo

         do mid=1,neq
            mi=mid+ndummy
            denco2i(mi)=(sto1(mi)  -denco2h(mi))*dtin
            deneco2i(mi)=(sto1(mi+neq)-deneco2h(mi))*dtin
         enddo    

c     ****************************************************************
      else if (iflg.eq.3) then
c     
c     calculations for air 
c     
         do mid=1,neq
            mi=mid+ndummy
            ieosd=ieos(mi)
            pw=phi(mi)
            tl=t(mi)
            pl=phico2(mi)
            frac_w = fw(mi)
            frac_cl=fl(mi)
            frac_cg=fg(mi)
            yco2 = yc(mi)
            ywat = yw(mi)
            yair = ya(mi)
            xco2 = xc(mi)
            xwat = xw(mi)
            xair = xa(mi)		
            cp=denr(mi)*cpr(mi)
            por=ps(mi)
            dporpl=dporp(mi)
            dportl=dport(mi)
            vol=volume(mi)
            dtps=dtpsc(mi)
            icesd = ices(mi)
c     rel-perms  
            rlw = rl_w(mi)
            drlww = drl_ww(mi)
            drlwg = drl_wg(mi)
            drlwp = 0.0
            drlwt = 0.0
            if(icesd.eq.2) then
               rll = rl_l(mi)
               drllw = drl_lw(mi)
               drllg = drl_lg(mi)
               drllp = 0.0
               drllt = 0.0
            else
               rll = rl_v(mi)
               drllw = drl_vw(mi)
               drllg = drl_vg(mi)
               drllp = 0.0
               drllt = 0.0
            endif
c     
c     The following represents properties for co2-rich phase. It is assumed that
c     if 2 phase CO2 is present the dissolved water is in CO2 gas phase.
c     
            if(icesd.eq.2) then 
               rol = co2_prop(mi)
               drolt=co2_prop(neq+mi)
               drolp=co2_prop(2*neq+mi)
               drolw=0.d0
               drolya=0.d0
               drolyc=0.d0
               enl=co2_prop(3*neq+mi)
               denlt=co2_prop(4*neq+mi)
               denlp=co2_prop(5*neq+mi)
               denlw=0.d0
               denlya=0.d0
               denlyc=0.d0
               visl=co2_prop(6*neq+mi) + vis_tol
               dvislt=co2_prop(7*neq+mi)
               dvislp=co2_prop(8*neq+mi)
               dvislya=0.d0
               dvislyc=0.d0
               dvislw=0.d0
            else
               rol = co2_prop(9*neq+mi)
               drolt=co2_prop(10*neq+mi)
               drolp=co2_prop(11*neq+mi)
               drolw=0.d0
               drolya=0.d0
               drolyc=0.d0
               enl=co2_prop(12*neq+mi)
               denlt=co2_prop(13*neq+mi)
               denlp=co2_prop(14*neq+mi)
               denlw=0.d0
               denlya=0.d0
               denlyc=0.d0
               visl=co2_prop(15*neq+mi) + vis_tol
               dvislt=co2_prop(16*neq+mi)
               dvislp=co2_prop(17*neq+mi)
               dvislya=0.d0
               dvislyc=0.d0
               dvislw=0.d0
            endif				
            row = wat_prop(mi)
            drowp=wat_prop(neq+mi)
            drowt=wat_prop(2*neq+mi)
            drowyc=wat_prop(3*neq+mi)
            drowya=wat_prop(4*neq+mi)	
            enw=wat_prop(5*neq+mi)
            denwp=wat_prop(6*neq+mi)
            denwt=wat_prop(7*neq+mi)
            visw=wat_prop(8*neq+mi) + vis_tol
            dviswp=wat_prop(9*neq+mi)
            dviswt=wat_prop(10*neq+mi)

c     modify derivatives for 2-phase (should now depend on p, not t)
            if((abs(ieosd).eq.2).or.(abs(icesd).eq.2)) then
               drolp=drolp+drolt*dtps
               drovp=drovp+drovt*dtps
               dhlp=dhlp+dhlt*dtps
               dhvp=dhvp+dhvt*dtps
               dvisvp=dvisvp+dvisvt*dtps
               dvislp=dvislp+dvislt*dtps
               dporpl=dporpl+dportl*dtps
               drolt=0.0
               drovt=0.0
               dhlt=0.0
               dhvt=0.0
               dvisvt=0.0
               dvislt=0.0
               dportl=0.0
               dtd=1.0
            endif

c     source terms and its derivatives
c     sk : water mass production rate, originally declared in comdi, old variable
c     dq : derivative of sk wrt P, originally declared in comci, old variable
c     dqt: derivative of sk wrt T, originally declared in comci, old variable
c     dqw: derivative of sk wrt fw, originally declared in comco2, new variable
c     dqxc: derivative of sk wrt xco2, originally declared in comco2, new variable
c     dqyc: derivative of sk wrt yco2, originally declared in comco2, new variable
c     dqya: derivative of sk wrt ya, originally declared in comco2, new variable   
c     qh : water energy production rate, originally declared in comdi, old variable
c     dqh : derivative of qh wrt P, originally declared in comci, old variable
c     deqh: derivative of qh wrt T, originally declared in comci, old variable
c     dqhw: derivative of qh wrt fw, originally declared in comco2, new variable
c     dqhxc: derivative of qh wrt xco2, originally declared in comco2, new variable
c     dqhyc: derivative of qh wrt yco2, originally declared in comco2, new variable
c     dqhya: derivative of qh wrt ya, originally declared in comco2, new variable   
c     the algorithm used for calculating source/sink terms depends on whether
c     wellbore model is used or not.
            
            if(iriver.ne.1) then
               dq(mi)=0.0
               dqt(mi)=0.0
               dqw(mi) = 0.0
               dqyc(mi)=0.0
               dqya(mi)=0.0
               qh(mi)=0.0
               dqh(mi)=0.0
               deqh(mi)=0.0
               dqhw(mi)=0.0
               dqhyc(mi)=0.0
               dqhya(mi)=0.0
               
               qdis=sk(mi)
               kq=ka(mi)
               if(kq.lt.0) then
                  qh(mi)=0.
                  sk(mi)=0.0
                  permsd=wellim(mi)
c     form pressure dependent flow term
c     RJP 02/13/07 using the formulation used by Coats, Thomas & Pierson, 
c     SPE 29111, 1995.
                  if(pflow(mi).gt.0.0) then
                     pldif=pw-pflow(mi)
                     dpldt=0.0
                  else
                     pldif=pw-(psatd-pflow(mi))
                     dpldt=-dpsatt
                  endif
                  if(pldif.le.0.0d00) then
                     if(kq.eq.-2) then
                        permsd=0.
                     else
                        permsd=permsd*(rlw*row*yair/visw+
     &                       rll*rol*xair/visl)
                        dprmp=permsd*(rlw*yair*(drowp/visw-
     &                       row*dviswp/visw*visw)+rll*xair*
     &                       (drolp/visl-rol*dvislp/visl*visl))
                        dprmt=permsd*(rlw*yair*(drowt/visw-row*dviswt/
     &                       visw*visw)+rll*xair*(drolt/visl-rol*dvislt/
     &                       visl*visl))
                        dprmw=permsd*(drlww*yair*row/visw+
     &                       drllw*xair*rol/visl)
                        dprmyc=permsd*rlw*yair*drowyc/visw
                        dprmya=permsd*(rlw*(row+yair*drowya)/visw)
                     endif
                  endif
                  qdis=permsd*pldif
                  sk(mi)=qdis
                  dq(mi)=permsd+dprmp*pldif
                  dqt(mi)=dprmt*pldif
                  dqw(mi)=dprmw*pldif
                  dqyc(mi)=dprmyc*pldif
                  dqya(mi)=dprmya*pldif
               endif
            endif

c     compressed liquid
            eqdum=rol*frac_cl*xair+row*frac_w*yair
            den=por*eqdum
            damp=dtin*(dporpl*eqdum+
     &           por*(drolp*frac_cl*xair+drowp*frac_w*yair))
            damh=dtin*(dportl*eqdum+
     &           por*(drolt*frac_cl*xair+drowt*frac_w*yair))
            damw=por*dtin*(-rol*xair+row*yair)
            damyc=por*dtin*frac_w*drowyc*yair
            damya=por*dtin*frac_w*(drowya*yair+row)
            eqdum=frac_cl*(rol*enl-pl)*xair+
     &           frac_w*(row*enw-pl)*yair
            dene=por*eqdum
            daep=dtin*(dporpl*eqdum+
     &           por*(frac_cl*xair*(drolp*enl+rol*denlp-1)+
     &           frac_w*yair*(drowp*enw+row*denwp-1)))
            daeh=dtin*(dportl*eqdum+
     &           por*(frac_cl*xair*(drolt*enl+rol*denlt)+
     &           frac_w*yair*(drowt*enw+row*denwt)))
            daew=por*dtin*(-xair*(rol*enl-pl)+
     &           yair*(row*enw-pl))
            daeyc=por*dtin*frac_w*yair*drowyc*enw
            daeya=por*dtin*frac_w*(yair*drowya*enw+(row*enw-pl))
            if(ps(mi).ne.0.0) then          
c     store derivatives of accumulation terms
               dmpf(mi)=damp
               dmef(mi)=damh
               dmwf(mi)=damw
               dmycf(mi)=damyc
               dmyaf(mi)=damya
               depf(mi)=daep
               deef(mi)=daeh
               dewf(mi)=daew
               deycf(mi)=daeyc
               deyaf(mi)=daeya
               dtpa(mi)=dtps
               dtpae(mi)=dtd
               t(mi)=tl
c     save accumulation terms for possible volume changes
               sto1(mi)=den
               sto1(mi+neq)=dene
               dstm(mi)=por*rov*sv*vol
               diw(mi)=0.0
               diwp(mi)=0.0
               diwe(mi)=0.0
               diww(mi)=0.0
               diwyc(mi)=0.0
               diwya(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divw(mi)=0.0
               divyc(mi)=0.0
               divya(mi)=0.0
c     
c     modify flow terms for new upwind scheme
c     .First calculate transfer through water-rich phase
c     Note here l refers to water-rich phase and v refers to CO2-rich
c     (gas/liquid/SC) phase or air. 
c     
               dql=rlw*row*yair/visw
               diw(mi)=dql
               diwp(mi)=rlw*(drowp*yair/visw-
     &              row*yair*dviswp/visw**2)
               diwe(mi)=rlw*(drowt*yair/visw-
     &              row*yair*dviswt/visw**2)
               diww(mi)=drlww*row*yair/visw
               diwya(mi)=rlw*(row+drowya*yair)/visw
               diwyc(mi)=rlw*drowyc*yair/visw
c     
c     Next calculate transfer through CO2-rich phase
c     
               dqv=rll*rol*xair/visl
               div(mi)=dqv
               divp(mi)=rll*(drolp*xair/visl-
     &              rov*xair*dvislp/visl**2)
               dive(mi)=rll*(drolt*xair/visl-
     &              rov*xair*dvislt/visl**2)
               divw(mi)=drllw*rol*xair/visl
               divya(mi)=0.d0
               divyc(mi)=0.d0
               if(qdis.gt.0.0) then
                  
c     organize source terms and derivatives

                  hprod=xair*enl+yair*enw
                  dhprdp=denlp*xair+yair*denwp
                  dhprde=denlt*xair+yair*denwp
                  dhprdw=0.d0
                  dhprdyc=0.d0
                  dhprdya=-enw
                  qh(mi)=hprod*qdis
                  dqh(mi)=dhprdp*qdis+hprod*dq(mi)
                  deqh(mi)=dhprde*qdis+hprod*dqt(mi)
                  dqhw(mi)=dhprdw*qdis+hprod*dqw(mi)
                  dqhyc(mi)=dhprdyc*qdis+hprod*dqyc(mi)
                  dqhya(mi)=dhprdya*qdis+hprod*dqya(mi)
               endif
               if(qdis.le.0.) then
                  eskd=eflow(mi)
                  qh(mi)=qdis*eflow(mi)
                  dqh(mi)=eskd*dq(mi)
                  deqh(mi)=eskd*dqt(mi)
                  dqhw(mi)=eskd*dqw(mi)
                  dqhyc(mi)=eskd*dqyc(mi)
                  dqhya(mi)=eskd*dqya(mi)
               endif
            endif
c     
c     add intercomponent heat flux
c     this term transfers heat from water to co2 (and back)
c     also note this term only contains derivative wrt water   
c     derivative with respect to co2 can also be derived 
c     
            if(qhflxco2(mi).ne.0.0) then
               htc = qhflxco2(mi)
               hflux=htc*(t(mi)-tco2(mi))
               if(abs(ieosd).ne.2.and.abs(icesd).ne.2) then
                  dhflxp=0.0
                  dhflxe=htc
               else
                  dhflxp=htc*dtpa(mi)
                  dhflxe=0.0
               endif
c     this next if block contains derivative wrt co2 component
               if(abs(ices(mi)).ne.2) then
                  dhflxew=-htc
               else
                  dhflxpm=-htc*dtpaco2(mi)
                  dhflxem=0.0
               endif
               qh(mi) = qh(mi) + hflux
               dqh(mi) = dqh(mi) + dhflxp
               deqh(mi) = deqh(mi) + dhflxe
               deqpm(mi) = deqpm(mi) + dhflxpm
               deqm(mi) = deqm(mi) + dhflxem
            endif
c     
c     add heat source term
c     
            if(qflux(mi).ne.0.0) then
               if(qflxm(mi).gt.0.0) then
                  htc=qflxm(mi)
                  tbound=qflux(mi)
                  hflux=htc*(tl-tbound)
                  if(abs(ieosd).ne.2.and.abs(icesd).ne.2) then
                     dhflxp=0.0
                     dhflxe=htc
                  else
                     dhflxp=htc*dtps
                     dhflxe=0.0
                  endif
                  qh(mi)=qh(mi)+hflux
                  dqh(mi)=dqh(mi)+dhflxp
                  deqh(mi)=deqh(mi)+dhflxe
               else if(qflxm(mi).lt.0.0) then
                  htc=abs(qflxm(mi))
                  sbound=qflux(mi)
                  hflux=htc*(tl-sbound)
                  if(abs(ieosd).ne.2.and.abs(icesd).ne.2) then
                     dhflxp=0.0
                     dhflxe=0.0
                  else
                     dhflxp=0.0
                     dhflxe=htc
                  endif
                  qh(mi)=qh(mi)+hflux
                  dqh(mi)=dqh(mi)+dhflxp
                  deqh(mi)=deqh(mi)+dhflxe
               else
                  qh(mi)=qh(mi)+qflux(mi)
               endif
            endif
            if(ps(mi).eq.0.0) then
               
c     heat conduction only
               cprd=cpr(mi)
               if(kq.lt.0) then
                  eskd=eflow(mi)
                  edif=cprd*(tl-eskd)
                  permsd=wellim(mi)
                  qh(mi)=permsd*edif
                  deqh(mi)=permsd*cprd
                  dqhw(mi)=0.d0
                  dqhyc(mi)=0.d0
                  dqhya(mi)=0.d0
               endif
               if(kq.ge.0.and.qflux(mi).eq.0.0) deqh(mi)=0.
               dtpae(mi)=1.
               denrd=denr(mi)*cprd
               sto1(mi)=0.0
               sto1(mi+neq)=denrd*tl
               deef(mi)=denrd*dtin
            endif
         enddo
         
         do     mid=1,neq
            mi=mid+ndummy
            deni(mi)=(sto1(mi)-denh(mi))*dtin
            denei(mi)=(sto1(mi+neq)-deneh(mi))*dtin
         enddo    
         
      endif
      
      deallocate(sto1)
     
      return
      end

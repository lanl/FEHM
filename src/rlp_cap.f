      subroutine rlp_cap(ndummy)
!***********************************************************************
! Copyright 2009 Los Alamos National Security, LLC  All rights reserved
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
!D1 Calculate relative permeabilities and capillary pressures
!D1 and derivatives.
!D1
!***********************************************************************

      use comai
      use combi
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

      real*8 scutm,hmin,darcyf,tol_l,tol_u,su_cut
      parameter(scutm = 1.d-03)
      parameter(hmin = 1.d-18)
      parameter(darcyf = 1.d12)
      parameter(tol_l  = 1.d-5)
      parameter(tol_u  = 1.d-5)
      parameter(su_cut = 0.99d00)
c      parameter(su_cut = 0.70d00)

      integer i, j, k, ndummy, mrlp, it, ir, iphase, lphase, phase(3)
      integer ieosd, mi, icouple, ir2, irf, i2, aphase, itbl
      integer :: ireg = 1
      real*8 alpha,beta,alamda,alpi,smcut,slcut,ds,dhp,dhp1
      real*8 rp1,rp2,rp3,rp4,denom,star,hp,hp1,b1,s_star
      real*8 akf,akm,porf,permb,sl,sl1,sl2,rp5,rp6,rp7
      real*8 smcutm,smcutf,alpham,alamdam,facf,cp1,cp3,sw,sg
      real*8 rl, rv, drls, drvs, rl1, rv1, drls1, drvs1, cp, dpcp
      real*8 skip,dum1,dum2,dum3,dummy,rlc,drlc,rt,drt,drl_wk     
      logical null1,ex
      character model_type*10

c     zero everything
      if(icarb.ne.0) then
         pcp=0.
         pcg=0.
         dpcgw=0.
         dpcgg=0.
         dpcpw=0
         dpcpg=0 
         rl_w=0.;drl_ww=0.;drl_wg=0.;rl_l=0.;drl_lw=0.;drl_lg=0.
         rl_v=0.;drl_vw=0.;drl_vg=0.		   
      else
         if (ndummy .eq. 0) then
            rlf = 0.; rvf = 0.; drlef = 0.; drvef = 0.
            pcp=0.
            dpcef=0.
         end if
      endif

      do i = 1, neq
         mi = i + ndummy
         ieosd = ieos(mi)
c     set it to the relative permeability model number
         it = irlp(mi)
         if (it .eq. 0) then
c No relative permeability model is specified
            mrlp = 0
         else
            if (rlp_phase(it,3) .eq. 0) then
               mrlp = 2
            else
               mrlp = 3
            end if
         end if

         
         phase(2) = 0 
         phase(3) = 0
c     determine which phases are active
         if (icarb .eq. 0) then
            if (ieosd.eq.1) then
c     single phase (water)
               phase(1) = rlp_pos(it, 1)
            else if (ieosd .eq. 3) then
c     single phase (air)
               phase(1) = rlp_pos(it, 2)
               if (phase(1) .eq. 0) phase(1) = rlp_pos(it, 5)
             else
c     **  air/water or water/vapor	 
               k = 1
               sw = s(mi)
               phase(1) = rlp_pos(it, 1)
               phase(2) = rlp_pos(it, 2)
               if (phase(2) .eq. 0) phase(2) = rlp_pos(it, 5)
            end if
         else
c     co2
            if (fw(mi) .eq. 1.0) then
c     single phase (water)
               phase(1) = rlp_pos(it, 1)
            else if (fl(mi) .eq. 1.0) then
c     single phase (liquid/sc CO2)
               phase(1) = rlp_pos(it, 3)
            else if (fg(mi) .eq. 1.0) then
c     single phase (gaseous CO2)
               phase(1) = rlp_pos(it, 4)
            else if (fw(mi) .gt. 0.0) then
c     co2 with water 
               sw = fw(mi)
               phase(1) = rlp_pos(it, 1)
               if (ices(mi) .eq. 1 .or. ices(mi) .eq. 4) then
c     water/co2 liquid 
                  k = 2
                  phase(2) = rlp_pos(it, 3)
               else if (ices(mi) .eq. 3) then
c     water/co2 gas
                  k = 3
                  phase(2) = rlp_pos(it, 4)
               else
c     we're 3-phase
                  k = 0
                  phase(2) = rlp_pos(it, 3)
                  phase(3) = rlp_pos(it, 4)
               endif
            else
c     no water (liquid co2 / gas co2)
               sw = fl(mi)
               k = 4
               phase(1) = rlp_pos(it, 3)
               phase(2) = rlp_pos(it, 4)
            endif
         end if

         if (phase(2) .eq. 0) then
c     single phase conditions
            select case (rlp_phase(it, phase(1)))
            case (1)
c     single phase (water)
               rlf(mi) = 1.
            case (2)
c     single phase (air)
               rvf(mi) = 1.
            case (-1)
c     single phase (water - co2 problem)
               rl_w(mi) = 1.
            case (3)
c     single phase (liquid/sc CO2)
               rl_l(mi) = 1.
            case (4)
c     single phase (gaseous CO2)
               rl_v(mi) = 1.
            end select

         else if (mrlp .ne. 0) then
c     ** first handle two-phase cases 
c     phase(1) wetting; phase(2) non-wetting
            if (phase(3) .eq. 0) then
               aphase = 2
            else
               aphase = 3
            end if
c ***************
c Need to make sure we are using the parameters assigned to the active phase if we have a 3 phase problem and only 2-phases are active
c**************
c     calculate capillary pressure
            do icouple = 1, aphase - 1
               if (k .ne. 0) then
                  iphase = cap_pos(it, k)
                  if (iphase .eq. 0) exit
               else
               end if
               j = (iphase - 1) * max_cp
               select case (cap_type(it, iphase))
               case (1)
c     Linear
                  call linear(sw, cap_param(it, j + 1),
     &                    cap_param(it, j + 2), cp, dpcp)
                  select case (cap_coupling(it, iphase))
                  case (1, 5)
c     air/water or water/vapor 
                     pcp(mi) = cp
                     dpcef(mi) = dpcp
                  case (2)
c     co2_liquid/water
                     pcp(mi) = cp
                     dpcpw(mi) = dpcp
                     dpcpg(mi)=-1.*dpcpw(mi)
                  case (3)
c     co2_gas/water
                     pcg(mi) = cp
                     dpcgg = dpcp
                     dpcgw(mi)=-1.*dpcgg(mi)
                  end select
               case (2)
c     Linear Forsythe
c     
c     linear forsythe(1988) model  cp=f(S)
c     
                  pcp(mi)=cap_param(it, j + 1)*(cap_param(it, j + 2)-sw)
                  if (icarb .eq. 0) then
                     if (pcp(mi).lt.0.0) then
                        pcp(mi) = 0.0
                        dpcef(mi) = 0.
                     else
                        dpcef(mi) = -cap_param(it, j + 1)
                     end if
                  else
                     if(pcp(mi).lt.0.0) then
                        pcp(mi)=0.0
                        dpcpw(mi)=0.0
                     else
                        dpcpw(mi)=-cap_param(it, j + 1)
                     end if
                  endif
               case (3)
c     Exponential
                  call exponential(sw, cap_param(it, j + 1),
     &                 cap_param(it, j + 2), cap_param(it, j + 3),
     &                 cp, dpcp)
                  select case (cap_coupling(it, iphase))
                  case (1, 5)
c     air/water or water/vapor
                     pcp(mi) = cp
                     dpcef(mi) = dpcp
                  case (2)
c     co2_liquid/water
                     pcp(mi) = cp
                     dpcpw(mi) = dpcp
                     dpcpg(mi)=-1.*dpcpw(mi)
                  case (3)
c     co2_gas/water
                     pcg(mi) = cp
                     dpcgg = dpcp
                     dpcgw(mi)=-1.*dpcgg(mi)
                  end select
               case (4)
c     Brooks-Corey
                  if (icouple .eq. 1) then
                     call brooks_corey_cap(sw, cap_param(it, j + 1),
     &                    cap_param(it, j + 2), cap_param(it, j + 3),
     &                    cap_param(it, j + 4), cap_param(it, j + 5),
     &                    cap_param(it, j + 6), pcp(mi), dpcp)
                     if (icarb .eq. 0) then
                        dpcef(mi) = dpcp
                     else
                        dpcpw(mi) = dpcp
                        dpcpg(mi)=-1.*dpcpw(mi)
                     end if
                  else if (icouple .eq. 2) then
                     call brooks_corey_cap3(sg, cap_param(it, j + 1),
     &                    cap_param(it, j + 2), cap_param(it, j + 4),
     &                    cap_param(it, j + 3), pcg(mi), dpcgw(mi), 
     &                    dpcgg(mi))
                  end if
               case (5, 6, 7)
c     Van Genuchten
                  if(rlp_fparam(it, 7) .eq. 0. .or. 
     &                 icarb .ne. 0) then
                
                     call vgcap(sw, cap_param(it, j + 1), 
     &                    cap_param(it, j + 2), cap_param(it, j + 3), 
     &                    cap_param(it, j + 4), vg1(it, iphase), 
     &                    vg2(it, iphase), vg3(it, iphase), 
     &                    vg4(it, iphase), cap_param(it, j + 6),
     &                    su_cut, cp1s(it, iphase), cp2s(it, iphase),
     &                    pcp(mi), dpcp, ireg)
                     hp = pcp(mi)
                     dhp = dpcp
                     if (icarb .eq. 0) then
                        dpcef(mi) = dpcp
                     else
                        dpcpw(mi) = dpcp
                     end if
                  else
c     check for fractures, only for air/water case
                     if( idpdp .ne. 0 .or. idualp .ne. 0 ) then
                        if(mi.le.neq) then
                           call vgcap(sw, cap_fparam(it, j + 1),
     &                          cap_fparam(it, j + 2), 
     &                          cap_fparam(it, j + 3),
     &                          cap_fparam(it, j + 4), vg5(it, iphase), 
     &                          vg6(it, iphase), vg7(it, iphase), 
     &                          vg8(it, iphase), cap_fparam(it, j + 6),
     &                          su_cut, cp3s(it, iphase), 
     &                          cp4s(it, iphase), 
     &                          pcp(mi), dpcef(mi), ireg)
                        else
                           call vgcap(sw,  cap_param(it, j + 1),
     &                          cap_param(it, j + 2), 
     &                          cap_param(it, j + 3),
     &                          cap_param(it, j + 4), vg1(it, iphase), 
     &                          vg2(it, iphase), vg3(it, iphase), 
     &                          vg4(it, j + 1), cap_param(it, j + 6),
     &                          su_cut, cp1s(it, iphase), 
     &                          cp2s(it, iphase), 
     &                          pcp(mi), dpcef(mi) ,ireg)
                        end if
                        hp = pcp(mi)
                        dhp = dpcef(mi)
                     else
c     we're just using matrix cap pressure here
                        call vgcap(sw,  cap_param(it, j + 1),
     &                       cap_param(it, j + 2), cap_param(it, j + 3),
     &                       cap_param(it, j + 4), vg1(it, iphase), 
     &                       vg2(it, iphase), vg3(it, iphase), 
     &                       vg4(it, iphase), cap_param(it, j + 6), 
     &                       su_cut, cp1s(it, iphase), 
     &                       cp2s(it, iphase), hp, dhp, ireg)
                        pcp(mi) = hp
                        dpcef(mi) = dhp
                        if (cap_type(it, iphase) .eq. 7) then
                           call vgcap(sw, cap_fparam(it, j + 1),
     &                          cap_fparam(it, j + 2), 
     &                          cap_fparam(it, j + 3),
     &                          cap_fparam(it, j + 4), vg5(it, iphase), 
     &                          vg6(it, iphase), vg7(it, iphase), 
     &                          vg8(it, iphase), cap_fparam(it, j + 6),
     &                          su_cut, cp3s(it, iphase), 
     &                          cp4s(it, iphase), 
     &                          hp1, dhp1, ireg)
                        else
                           hp1 = hp
                           dhp1 = dhp
                        end if
                     end if
                  end if
c     conversion to MPa
                  pcp(mi) = 9.8e-3 * pcp(mi)
                  if (icarb .eq. 0) then
                     dpcef(mi) = 9.8e-3 * dpcef(mi)
                  else
                     dpcpw(mi) = 9.8e-3 * dpcpw(mi)
                     dpcpg(mi)= -1. * dpcpw(mi)
                  end if
               case (8)
c     values from table, called with wetting phase saturation, table number
                  itbl = rlp_param(it, j + 1)
                  call rlp_cap_table (sw, itbl, 4, cp, dpcp)

               end select
            end do
c     now do rlp *******************************************
            lphase = 0
            do i2 = 1, aphase
               iphase = phase (i2)
               ir = (iphase - 1) * max_rp
               irf = (iphase - 1) * max_rpf
               select case (rlp_phase(it, iphase))
               case (1)
c     Water
                  sw = s(mi)
               case (-1)
c     Water (co2 problem)
                  sw = fw(mi)
               case(2, 5)
c     Air
                  sw = 1.0 - s(mi)
               case (3)
c     Liquid or super-critical CO2
                  sw = fl(mi)
               case (4)
c     Gaseous co2
                  sw = fg(mi)
               end select

               select case (rlp_type(it, iphase))
               case (1)
c     Constant
                  select case (rlp_phase(it, iphase))
                  case (1)
                     rlf(mi) = rlp_param(it, ir + 1)
                  case(-1)
                     rl_w(mi) = rlp_param(it, ir + 1)
                  case (2)
                     rvf(mi) = rlp_param(it, ir + 1)
                  case (3)
                     rl_l(mi) = rlp_param(it, ir + 1)
                  case (4)
                     rl_v(mi) = rlp_param(it, ir + 1)
                  case (5)
                     rvf(mi) = rlp_param(it, ir + 1)
                  case (6)                    
                  end select
               case (2, 3)
c     Linear or exponential
                  if (rlp_type(it, iphase) .eq. 2) then
                     call linear(sw, rlp_param(it, ir + 1), 
     &                    rlp_param(it, ir + 2), rl, drls)
                  else if (rlp_type(it, iphase) .eq. 3) then
                     call exponential(sw, rlp_param(it, ir + 1),
     &                    rlp_param(it, ir + 2), rlp_param(it, ir + 3),
     &                    rl, drls)
                  end if
                  select case (rlp_phase(it, iphase))
                  case (1)
                     rlf(mi) = rl
                     drlef(mi) = drls
                  case(-1)
                     rl_w(mi) = rl
                     drl_ww(mi) = drls
                 case (2)
                     rvf(mi) = rl
                     drvef(mi) = -drls
                  case (3)
                     rl_l(mi) = rl
                     drl_lw(mi) = drls
                  case (4)
                     rl_v(mi) = rl
                     drl_vw(mi) = drls
                  case (5)
                     rvf(mi) = rl
                     drvef(mi) = -drls
                  case (6)                    
                  end select
               case (4)
c     Corey
                  if (i2 .eq. 1 .or. (rlp_type(it, 1) .ne. 
     &                 rlp_type(it, 2))) then
                     call corey(sw, rlp_param(it, ir + 1), 
     &                    rlp_param(it, ir + 2),rl, drls, rv, drvs)
                     if (rlp_type(it, 1) .eq. rlp_type(it, 2)) then
                        if (k .eq. 1) then
                           rlf(mi) = rl
                           drlef(mi) = drls
                           rvf(mi) = rv
                           drvef(mi) = drvs
                        else if (k .eq. 2) then
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                           rl_l(mi) = rv
                           drl_lw(mi) = drvs
                        else if (k .eq. 3) then
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs
                        else if (k .eq. 4) then
                           rl_l(mi) = rl
                           drl_lw(mi) = drls
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs
                        endif
                        exit
                     else 
                        select case (rlp_phase(it, iphase))
                        case (1)
c     water
                           rlf(mi) = rl
                           drlef(mi) = drls
                        case (-1)
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                        case (2)
c     air/vapor
                           rvf(mi) = rv
                           drvef(mi) = drvs
                        case (3)
c     liquid co2
                           if (k .eq. 4) then
                              rl_l(mi) = rl
                              drl_lw(mi) = drls
                           else
                              rl_l(mi) = rv
                              drl_lw(mi) = drvs
                           end if
                        case (4)
c     gaseous co2
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs 
                        end select 
                     end if
                  end if
               case (5)
c     Brooks Corey
                  if (aphase .eq. 3) then
c     Special case, all three phases active
                     do j = 1, 3
                        if (phase(1) .eq. rlp_phase(it, j)) then
                           iphase = j
                        else if  (phase(2) .eq. rlp_phase(it, j)) then
                           lphase = j
                        end if
                     end do
                     ir = (iphase - 1) * max_rp
                     ir2 = (lphase - 1) * max_rp
                     sl = fl(mi)
                     call brooks_corey_rlp(sw, rlp_param(it, ir + 1),
     &                    rlp_param(it, ir + 3), rl_w(mi), drl_ww(mi), 
     &                    rl_l(mi), drl_lw(mi))
                     call brooks_corey_rlp3(sw, sl, 
     &                    rlp_param(it, ir + 1), rlp_param(it, ir + 3),
     &                    rlp_param(it, ir2 + 1), rl_l(mi), drl_lw(mi),
     &                    drl_lg(mi), rl_v(mi), drl_vw(mi), drl_vg(mi))
                     exit
                  else
                     if (i2 .eq. 1) then
                        lphase = iphase
                     else if (rlp_type(it, iphase) .eq. 
     &                       rlp_type(it, lphase)) then
                        exit
                     end if
                     call brooks_corey_rlp(sw, rlp_param(it, ir + 1),
     &                    rlp_param(it, ir + 3), rl, drls, dum1, rv, 
     &                    drvs, dum2)
                     if (rlp_type(it, iphase) .eq. 
     &                    rlp_type(it, lphase)) then
                        if (k .eq. 1) then
                           rlf(mi) = rl
                           drlef(mi) = drls
                           rvf(mi) = rv
                           drvef(mi) = drvs
                        else if (k .eq. 2) then
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                           rl_l(mi) = rv
                           drl_lw(mi) = drvs
                        else if (k .eq. 3) then
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs
                        else if (k .eq. 4) then
                           rl_l(mi) = rl
                           drl_lw(mi) = drls
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs
                        endif
                     else 
                        select case (rlp_phase(it, iphase))
                        case (1)
c     water
                           rlf(mi) = rl
                           drlef(mi) = drls
                        case (-1)
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                        case (2, 5)
c     air/vapor
                           rvf(mi) = rv
                           drvef(mi) = drvs
                        case (3)
c     liquid co2
                           if (k .eq. 4) then
                              rl_l(mi) = rl
                              drl_lw(mi) = drls
                           else
                              rl_l(mi) = rv
                              drl_lw(mi) = drvs
                           end if
                        case (4)
c     gaseous co2
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs 
                        end select 
                     end if                                             
                  endif
               case (6 , 7)
c     Van Genuchten
                  akf = rlp_fparam(it, irf + 5)
                  akm = rlp_fparam(it, irf + 6)
                  porf = rlp_fparam(it, irf + 7)

                  if (i2 .eq. 1 .or. (rlp_type(it, 1) .ne. 
     &                 rlp_type(it, 2))) then 
                     if (rlp_type(it, iphase) .eq. 6) then
c Van Genuchten using saturations
                        if (porf .eq. 0. .or. icarb .ne. 0) then
c     No fracture model
                           star = (sw - rlp_param(it, ir + 1)) /
     &                          (rlp_param(it, ir + 2) - 
     &                          rlp_param(it, ir + 1))
                           call vgrlps(2, sw, star, 
     &                          rlp_param(it, ir + 3),
     &                          rlp_param(it, ir + 4), 
     &                          rlp_param(it, ir + 1),
     &                          rlp_param(it, ir + 2), tol_l, tol_u,  
     &                          rl, drls, rv, drvs)
                        else
                           if (idpdp .ne. 0 .or. idualp .ne. 0) then
c     Fracture model (double or dual porosity simulation)
                              if( mi .le. neq ) then
                                 star = (sw - rlp_fparam(it, ir + 1)) /
     &                                (rlp_fparam(it, ir + 2) -
     &                                rlp_fparam(it, ir + 1))
                                 call vgrlps(2, sw, star, 
     &                                rlp_fparam(it, irf + 3),
     &                                rlp_fparam(it, irf + 4), 
     &                                rlp_fparam(it, irf + 1),
     &                                rlp_fparam(it, irf + 2), tol_l, 
     &                                tol_u, rl, drls, rv, drvs)
                                 permb = akf * porf    
                              else
                                 star = (sw - rlp_param(it, ir + 1)) /
     &                                (rlp_param(it, ir + 2) - 
     &                                rlp_param(it, ir + 1))
                                 call vgrlps(2, sw, star, 
     &                                rlp_param(it, ir + 3),
     &                                rlp_param(it, ir + 4), 
     &                                rlp_param(it, ir + 1),
     &                                rlp_param(it, ir + 2), tol_l, 
     &                                tol_u, rl, drls, rv, drvs)
                                 permb = akm * (1. - porf)
                              end if
                              if (rlp_fparam(it, irf + 8) .ne. 0.) then
                                 call rlp_frac2(1, mi, sw, star,
     &                                rl, drls, 1.-rl, -drls,
     &                                rlp_fparam(it, irf + 8))
                              end if
                           else
c     Fracture model, continuum
                              star = (sw - rlp_param(it, ir + 1)) /
     &                             (rlp_param(it, ir + 2) - 
     &                             rlp_param(it, ir + 1))
                              call vgrlps(2, sw, star, 
     &                             rlp_param(it, ir + 3),
     &                             rlp_param(it, ir + 4), 
     &                             rlp_param(it, ir + 1),
     &                             rlp_param(it, ir + 2), tol_l,
     &                             tol_u, rl, drls, rv, drvs)
                              star = (sw - rlp_fparam(it, ir + 1)) /
     &                             (rlp_fparam(it, ir + 2) -
     &                             rlp_fparam(it, ir + 1))
                              call vgrlps(2, sw, star, 
     &                             rlp_fparam(it, irf + 3),
     &                             rlp_fparam(it, irf + 4), 
     &                             rlp_fparam(it, irf + 1),
     &                             rlp_fparam(it, irf + 2), tol_l, 
     &                             tol_u, rl1, drls1, rv1, drvs1)
                              permb = akf * porf + akm * (1. - porf)
                              rl = (akf*rl1*porf + akm*rl*(1.0 - porf))
     &                             /permb
                              drls = (akf*drls1*porf + akm*drls*
     &                             (1.0-porf))/permb
                              rv = 1. - rl
                              drvs = -drls
                           end if
                        end if
                     else if (rlp_type(it, iphase) .eq. 7) then
c     Van Genuchten using capillary pressures
c     this assumes cap pressures (hp,dhp) are already available
                        if (porf .eq. 0. .or. icarb .ne. 0) then
c     No fracture model
                           call vgrlp2(sw, rlp_param(it, ir + 3),
     &                       rlp_param(it, ir + 4), hmin, hp, dhp,
     &                       rl, drls, rv, drvs, 0)
                        else if (idpdp .ne. 0 .or. idualp .ne. 0) then
c     Fracture model (double or dual porosity simulation)
                           if( mi .le. neq ) then
                              call vgrlp2(sw, rlp_fparam(it, irf + 3),
     &                             rlp_fparam(it, irf + 4), hmin, hp, 
     &                             dhp, rl, drls, rv, drvs, 0)
                              permb = akf * porf                        
                           else
                              call vgrlp2(sw, rlp_param(it, ir + 3),
     &                             rlp_param(it, ir + 4), hmin, hp, dhp,
     &                             rl, drls, rv, drvs, 0)
                              permb = akm * (1. - porf)
                           end if
                        else
c     Fracture model, continuum
                           call vgrlp2(sw, rlp_param(it, ir + 3),
     &                          rlp_param(it, ir + 4), hmin, hp, dhp,
     &                          rl, drls, rv, drvs, 0)
                           call vgrlp2(sw, rlp_fparam(it, irf + 3),
     &                          rlp_fparam(it, irf + 4), hmin, hp1, 
     &                          dhp1, rl1, drls1, rv1, drvs1, 0)
                           permb = akf * porf + akm * (1. - porf)
                           rl = (akf*rl1*porf + akm*rl*(1.0-porf))/permb
                           drls = (akf*drls1*porf + akm*drls*(1.0-porf))
     &                          /permb
                           rv = 1. - rl
                           drvs = -drls
                        end if
                     end if
                     if (porf .ne. 0) then
                        if(iupk.gt.0) then
c     upstream weighting of intrinsic permeability
                           permb = permb*darcyf
                           rl = rl*permb
                           rv = rv*permb
                           drls = drls*permb
                           drvs = drvs*permb
                           permb = 1.0/darcyf
c     set saturated permeabilities(assume isotropic)
                           pnx(mi) = permb*1.e+6
                           pnz(mi) = pnx(mi)
                           pny(mi) = pnx(mi) 
                        else
c     set saturated permeabilities(assume isotropic)
                           if(porf.ge.0.0) then
                              pnx(mi) = permb*1.e+6
                              pny(mi) = pnx(mi)
                              pnz(mi) = pnx(mi)
                           endif
                        endif
                     end if
                     if (rlp_type(it, 1) .eq. rlp_type(it, 2)) then
                        if (k .eq. 1) then
                           rlf(mi) = rl
                           drlef(mi) = drls
                           rvf(mi) = rv
                           drvef(mi) = drvs
                        else if (k .eq. 2) then
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                           rl_l(mi) = rv
                           drl_lw(mi) = drvs
                        else if (k .eq. 3) then
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs
                        else if (k .eq. 4) then
                           rl_l(mi) = rl
                           drl_lw(mi) = drls
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs
                        endif
                        exit
                     else 
                        select case (rlp_phase(it, iphase))
                        case (1)
c     water
                           rlf(mi) = rl
                           drlef(mi) = drls
                        case (-1)
                           rl_w(mi) = rl
                           drl_ww(mi) = drls
                        case (2, 5)
c     air/vapor
                           rvf(mi) = rv
                           drvef(mi) = drvs
                        case (3)
c     liquid co2
                           if (k .eq. 4) then
                              rl_l(mi) = rl
                              drl_lw(mi) = drls
                           else
                              rl_l(mi) = rv
                              drl_lw(mi) = drvs
                           end if
                        case (4)
c     gaseous co2
                           rl_v(mi) = rv
                           drl_vw(mi) = drvs 
                        end select 
                     endif
                  end if
               case (8)
c     values from table, called with wetting phase saturation, table number
                  itbl = rlp_param(it, j + 1)
                  call rlp_cap_table (sw, itbl, 2, rt, drt)
                  select case (rlp_phase(it, iphase))
                  case (1)
                     rl_w(mi) = rt
                     drl_ww(mi) = drt
                  case (7)
c                     rl_g(mi) = rt
                     drl_wg(mi) = drt
                  end select

               case (9)
                  itbl = rlp_param(it, j + 1)
                  call rlp_cap_table (sw, itbl, 3, rl, drls)
                  call rlp_cap_table (rlp_param(it, j+3), itbl, 3, rlc, 
     &                 drlc)
                  itbl = rlp_param(it, j + 2)
                  call rlp_cap_table (sw, itbl, 3, rv, drvs)
                  call stone_rlp(rlc, drlc, rl_w(mi), rl_v(mi), rl, 
     &                 drls, rv, drvs, rt, drt)
                  rl_l(mi) = rt
c                  drl_wk(mi) = drt
               end select
            end do
         else
c     No model was defined
         end if
         if (icarb .eq. 0) then
            drlpf(mi) = 0.
            drvpf(mi) = 0.
         end if
      end do

      end subroutine rlp_cap


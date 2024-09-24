      subroutine  wrtout (tassem,tas,totalflin,totalein,curinflow,
     &     cureinflow,is_ch,is_ch_t)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.s
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine handles output.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/wrtout.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:30 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.14   Mon Mar 31 08:45:50 1997   gaz
CD2 changes for head output
CD2 
CD2    Rev 1.12   Wed May 29 16:07:00 1996   hend
CD2 Added output kg/s and MJ/s
CD2 
CD2    Rev 1.11   Mon May 20 14:38:28 1996   hend
CD2 Added to Screen Output iptty
CD2 
CD2    Rev 1.10   Thu May 16 12:44:04 1996   hend
CD2 Added enthalphy info to screen output
CD2 
CD2    Rev 1.9   Tue May 14 14:32:58 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.8   Wed May 08 15:06:10 1996   hend
CD2 Fixed nodal Energy Out column data -- qh not rqhd
CD2 
CD2    Rev 1.7   Wed May 08 14:17:24 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.6   Thu Apr 04 12:23:22 1996   hend
CD2 Tabbed code for readability
CD2 
CD2    Rev 1.5   Thu Jan 11 12:37:34 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   11/27/95 17:00:22   gaz
CD2 format change
CD2 
CD2    Rev 1.1   03/18/94 15:43:34   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:36   pvcs
CD2 original version in process of being certified
CD2 
c 17-mar-94 gaz
c got rid of b calls
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use davidi
      use comdti
      use comai
      use comxi
      use comii
      use comwt
      use comflow
      use comco2
      implicit none

      real*8 totalflin,totalein,curinflow,cureinflow
      real*8 phod, dumconv, dumconv1, rolconv
      real*8 aiter, aminkt, years, dayold, sl, eqd, rhomd
      real*8 hmd, rqd, tas, tassem, tdum, pdum, pres_out, qtmp
      integer ilev, mlev, il, i, md
      integer dummyint, m2lev
      integer izone, inode, inneq, iconn, indexa_axy
      integer addnode, idummy, i1, i2, is_ch, is_ch_t
      integer isalt_save
      real*8 sumfout, sumsink, sumsource, sumboun 
      logical matrix_node
      character*20 message_ts

      if (istrs.ge.0)  then
         
         if (ntty.eq.2) write(iout,7) 
 7       format('****************************************************'
     &        ,'*****************')
         if (ntty.eq.2)  write(iout,6000)  l
         if (iptty.gt.0) write(iptty,777)  l
 6000    format(1x,'Time Step',i10)
 777     format(/,1x,'Time Step',i10)

         iad    =  max0( 1,iad )
         aiter  =  dfloat( itert )/dfloat( iad )
         aminkt =  dfloat( minkt )/dfloat( iad )
         years  =  days  /  365.25d00
         dayold =  dtotdm/86400.00d00
         
c     zero out enthalpy for air water system
         
         if (ihf.ne.ihs .or. .not. compute_flow)  then
            if (ntty.eq.2) then
               write(iout,781) 
               write(iout,782)
               write(iout,783) years,days,dayold
            endif
            if (iptty.gt.0) then
               write(iptty,781)
               write(iptty,782)
               write(iptty,783) years,days,dayold
            endif
 781        format(/,20x,'Timing Information')
 782        format(11x,'Years',14x,'Days',9x,'Step Size (Days)')
 783        format(4x, g16.9, 4x, g16.9, 3x, e16.9,/,1x,
     &           'Heat and Mass Solution Disabled')
            goto  40
         endif
         if (ntty.eq.2) then
            write(iout,772) 
            write(iout,778)
            write(iout,78) years,days,dayold
            write(iout,703) tassem,tas
            write(iout,773)
            write(iout,75) iad
            write(iout,76) aiter
            write(iout,77) aminkt
            write(iout,704) itotal,itotals
            write(iout,705) is_ch, is_ch_t
            write(iout,706) nphase_liq, dnphase_liq, nphase_2, 
     &            dnphase_2, nphase_gas, dnphase_gas, nphase_sc,
     &            dnphase_sc 
         endif
         if (iptty.gt.0) then
            write(iptty,772)
            write(iptty,778)
            write(iptty,78) years,days,dayold
            write(iptty,703) tassem,tas
            write(iptty,773)
            write(iptty,75) iad
            write(iptty,76) aiter
            write(iptty,77) aminkt
            write(iptty,704) itotal,itotals
            write(iptty,705) is_ch, is_ch_t
            write(iptty,706) nphase_liq, dnphase_liq, nphase_2, 
     &            dnphase_2, nphase_gas, dnphase_gas, nphase_sc,
     &            dnphase_sc 
         endif
         if(fimp.le.1.0d00) then
            message_ts = '                    '
         else
            message_ts = '   Time Step Limited'
         endif
         if(impf.eq.1) then
            if (ntty.eq.2) write(iout,5996) fimp*delpt,message_ts
            if(iptty.gt.0) write(iptty,5996) fimp*delpt,message_ts
         else if(impf.eq.2) then
            if (ntty.eq.2) write(iout,5997) fimp*deltt,message_ts
            if(iptty.gt.0) write(iptty,5997) fimp*deltt,message_ts
         else if(impf.eq.3) then
            if (ntty.eq.2) write(iout,5998) fimp*delst,message_ts
            if(iptty.gt.0) write(iptty,5998) fimp*delst,message_ts
         else if(impf.eq.4) then
            if (ntty.eq.2) write(iout,5999) fimp*delat,message_ts
            if(iptty.gt.0) write(iptty,5999) fimp*delat,message_ts
         endif
 5996    format(1x,'Largest pressure change = ',g16.9,a20)
 5997    format(1x,'Largest temperature change = ',f10.4,a20)
 5998    format(1x,'Largest saturation change = ',f10.4,a20)
 5999    format(1x,'Largest air pressure change = ',g16.9,a20)
 772     format(/,20x,'Timing Information')
 778     format(11x,'Years',14x,'Days',9x,'Step Size (Days)')
 78      format(4x,g16.9,4x,g16.9,3x,g16.9)
 703     format(1x,'Cpu Sec for Time Step = ',g11.4,
     &        ' Current Total = ',g13.4)
 773     format(/,20x,'Equation Performance')
 75      format(1x,'Number of N-R Iterations: ',1i10)
 76      format(1x,'Avg # of Linear Equation Solver Iterations: ',
     &        1f5.1)
 77      format(1x,'Number of Active Nodes: ',1f10.0)
 704     format(1x,'Total Number of Iterations, N-R: ',i10,
     &        ' , Solver: ',i10)
 705     format(1x,'Phase Changes This Time Step: ',i8,' Total ',i11)
706     format(1x,'Nodes Liq Phase: ',i8,' change ',i8,/,
     &   ' Nodes Two Phase: ',i8,' change ',i8,/,' Nodes Gas Phase: ',
     &     i9,' change ', i8,/,' Nodes SC Phase: ',i8,' change ',i8)
         if(ifree.ne.0) then
          if(ntty.eq.2) then
            write(iout,*) 'Number of partially filled cells ', ifree1
          endif
          if(iatty.gt.0) then
            write(iatty,*) 'Number of partially filled cells ', ifree1
          endif
         endif
         if(ntty.eq.2) then
            write(iout,*) 'Number of restarted time steps ', nrestart_ts
         endif
         if(iatty.gt.0) then
            write(iatty,*)'Number of restarted time steps ', nrestart_ts
         endif
c     insert diagnostics code here
         if (ihf.eq.ihs .and. compute_flow) then
            call diagnostics(1)
            call diagnostics(2)
         endif
         if(iccen.ne.0) then
          resid_out = .false.
         endif
         if ((ntty.eq.2).and.(m.gt.0).and.resid_out) then
            write(iout,774)
            if (iatty.gt.0)  write(iatty,774)
 774        format(3x,'Node',3x,'Equation 1 Residual',3x,
     &           'Equation 2 Residual')
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
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
            do il=1,ilev
               if(il.ne.1) then
                  write(iout,776) il
                  if (iatty.gt.0) write(iatty,776) il
 776              format(3x,'Matrix Level = ',i1)
               endif
               do i=1,mlev
                  md=nskw(i+(il-1)*mlev)
                  write(iout,775) md,bp(md),bp(md+neq)
                  if (iatty.gt.0) write(iatty,775) md,bp(md),bp(md+neq)
 775              format(i7,4x,e14.6,9x,e14.6)
               enddo
            enddo
         endif
         
         if ( ntty .eq. 2 )  then
            if ( m .gt. 0 .and. icarb .ne. 1)  then
c zvd 02-27-09 Output nodal water information in icectrco2 for CO2 problems
               write(iout,779)
               if (iatty.gt.0) write(iatty,779)
 779           format(/,20x,'Nodal Information (Water)')
c Isothermal
               if(ico2.lt.0) then
                  if(ihead.eq.0.and.ice.eq.0.and.ichead.eq.0) then
                     write(iout,6230)
                     if (iatty.gt.0)  write(iatty,6230)
                  else
                     write(iout,6130)
                     if (iatty.gt.0)  write(iatty,6130)
                  end if
               else
                  if(ichead.eq.0) then
                     write(iout,6030)
                     if (iatty.gt.0)  write(iatty,6030)
                  else
                     write(iout,6330)
                     if (iatty.gt.0)  write(iatty,6330)
                  end if
               endif
c gaz 110121 added State (ieos()_ info
c Heat and mass w/o head output
 6030          format(50x, 'source/sink', 2x, 'E source/sink', /, 3x,
     &              'Node', 2x, ' P (MPa)', 3x, ' E (MJ)', 4x, 'S liq',
     &     5x, 'Temp (C)', 4x, '(kg/s)', 7x, '(MJ/s)',2x,'State',
     &     12x,'Zone')
c Heat and mass w/ head output
 6330          format(61x, 'source/sink', 2x, 'E source/sink', /,
     &              3x, 'Node', 2x, 'Head (m)', 3x, 'P (MPa)', 4x,
     &              'E (MJ)', 4x, 'S liq', 5x, 'Temp (C)', 5x,
     &              '(kg/s)', 7x, '(MJ/s)',2x,'State',
     &                12x,'Zone')
c Isothermal output w/ head
 6130          format(62x,'source/sink', /, 3x, 'Node', 2x, 'Head (m)',
     &              4x, 'P (MPa)', 4x, 'E (MJ)', 4x, '  S liq',3x,
     &              'Temp (C)', 4x, '(kg/s)','     State','  Zone' )
c Isothermal output w/o head
 6230          format(52x, 'source/sink', /, 3x,'Node',2x,' P (MPa) ',
     &              3x, ' E (MJ)', 4x, 'S liq', 5x, 'Temp (C)', 5x,
     &              '(kg/s)','      State','  Zone')
               
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
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
               if (ichead .ne. 0) then
                  ihead=1
                  dumconv = crl(1,1)
                  dumconv1 = crl(4,1)
                  pdum = pres0+rol0*head0*(-grav)
                  tdum = temp0        
                  call water_density(tdum,pdum,rolconv)
                  crl(1,1)=rolconv
                  crl(4,1)=pres0
                  rho1grav = rolconv*9.81d-6                  
               end if
               do il=1,ilev
                  if(il.ne.1) then
                     write(iout,780) il
                     if (iatty.gt.0) write(iatty,780) il
 780                 format(2x,'Matrix Level = ',i1)
                  endif
                  do i=1,mlev
                     md=nskw(i+(il-1)*mlev)
c gaz 110123
                     if(idoff.eq.-1)then 
                      sl =1.0
c                     else
c                      sl = s(md)
                     endif
                     if (irdof .eq. 13 .and. ifree .ne. 0) then
                        if(izone_free_nodes(md).gt.1) then        
                           sl=min(s(md)-rlptol,1.00d0)     
                        else
c	                   sl = s(md) 
                         if(abs(sl).lt.1.d-98) sl = 1.d-98
                        endif
                     else if(irdof .ne. 13) then
                        sl = s(md)
                        if(abs(sl).lt.1.d-98) sl = 1.d-98
                     else   
                        sl = 1.0d0
                     endif
                     sv=1.0-sl
                     if(ico2.lt.0.and.ice.eq.0) then
                        eqd=0.0
                     else
                       if(ieos(md).ne.0.and.ps(md).ne.0.) then
                        eqd=0.0
                        rhomd = sl * rolf(md) + sv * rovf(md)
                        hmd=sl*rolf(md)*enlf(md)+sv*rovf(md)*envf(md)
                        if (abs(rhomd) .gt. zero_t)  eqd = hmd/rhomd
                       else
                         if(ivrock.ne.0) then
                          eqd = urock(md)/denr(md)
                         else
                          eqd = cpr(md)*t(md)  
                         endif
                       endif
                     endif
                     rqd    =  sk(md)
c     rqhd   =  0.0
c     if ( abs( rqd ) .gt. zero_t )  rqhd =  qh(md)/rqd
c     CHANGE ABOVE TO JUST PRINT OUT qh ARRAY
                     if(ico2.lt.0) then
	                pres_out = pho(md)-crl(1,1)*head0*(-grav)
                     else
	                pres_out = pho(md)                                
                     endif
                     if(ihead.eq.0.and.ichead.eq.0) then
                        phod=pho(md)
                        if (ico2.lt.0) then
c gaz 062723 changed izonef(:) to zones_char(:)
                           write(iout, 6029)  md ,
     *                          phod , eqd , sl , t(md) , rqd, ieos(md),
     &                          zones_char(md)                      
                           if ( iatty .gt. 0 )  write(iatty ,6029) md ,
     *                          phod , eqd , sl , t(md) , rqd, ieos(md),
     &                          zones_char(md)                      
                        else
                           write(iout, 6031)  md ,
     *             phod , eqd , sl , t(md) , rqd , qh(md), ieos(md),
     &             zones_char(md) 
                           if ( iatty .gt. 0 )  write(iatty ,6031) md ,
     *             phod , eqd , sl , t(md) , rqd , qh(md), ieos(md),
     &             zones_char(md) 
                        end if 
                     else
                        if(ichead.eq.0) then
                           call headctr(4,md,pho(md),phod)
                           if(sl.lt.rlptol+sattol) phod=head_id
                        else
                           call headctr(4,md,pho(md),phod)  
                        endif
c     phod is head with offset removed
                        if (ico2 .lt. 0) then
                          write(iout, 6031)  md , phod , pres_out ,
     &                    eqd , sl , t(md) , rqd , ieos(md),
     &                    zones_char(md)(1:30)
                        if ( iatty .gt. 0 )  write(iatty ,6031)  md ,
     *                       phod, pres_out , eqd , sl , t(md) , rqd,
     &                        ieos(md), zones_char(md)(1:30)
                        else
                           write(iout, 6032)  md , phod , pres_out ,
     &                          eqd , sl , t(md) , rqd , qh(md), 
     &                          ieos(md), zones_char(md)(1:30) 
                           if ( iatty .gt. 0 )  write(iatty ,6032) md ,
     *                          phod, pres_out , eqd , sl , t(md) , 
     &                          rqd , qh(md), ieos(md), 
     &                          zones_char(md)(1:30)
                        end if
                     endif
c gaz 062723 changed 6031 format i7 to a30 
 6029                format(i7,2x,g11.4,1x,g9.3,1x,f9.5,1x,
     *                    g11.4,1x,g10.3,1x,i5,5x,a30)
 6031                format(i7,2x,g11.4,1x,g9.3,1x,g9.3,1x,f9.5,1x,
     *                    g11.4,1x,g10.3,1x,i5,5x,a30)
 6032                format(i7,1x,g11.4,1x,g9.3,1x,g9.3,1x,g9.3,1x,
     *                    f8.3,3x,g11.4,1x,g11.4,1x,i5,5x,a30)
                  enddo
               enddo
               if (ichead .ne. 0) then
                  crl(1,1)= dumconv
                  crl(4,1)= dumconv1
                  ihead=0
                  if(ico2.lt.0) then
                     rho1grav = crl(1,1)*(9.81d-6)
                  else
                     rho1grav = rol0*9.81d-6
                  endif
               end if
c     
c**** call varible porosity output ****
c     
               if (iporos.ne.0) then
                 call porosi (2)
                 call saltctr(8,0,0.0d00,0.0d00)
               elseif(isalt. ne.0) then
                call saltctr(7,0,0.0d00,0.0d00)
                call saltctr(8,0,0.0d00,0.0d00)
               elseif(cden) then
                isalt_save = isalt
                 isalt = 1
                 call saltctr(8,0,0.0d00,0.0d00)
                isalt = isalt_save
               endif
      
c     
c**** output for co2 and ngas****
c     
               if (ico2.gt.0) call  co2ctr  ( 5 )
               if (ico2.lt.0.and.ice.eq.0) then
                  call  airctr  ( 5,0 )
               else if (ico2.lt.0.and.ice.ne.0) then
                  call  icectr  ( 5,0 )
                  call  icectr  (-5,0 )
               endif
            else if(icarb.eq.1) then
             
               call icectrco2(5,0)
               call icectrco2(-5,0)      

            endif
c     
c     compute and printout fluxes
c     
            call flxo(2)
c
c     compute and printout fluxes for CO2
c
            if(icarb.eq.1)  call flxo(3)

c     Compute flux passing thorugh a zone

            if(nflxz.ne.0) then

               call flxz (1, 0.d0)
               call eflxz(1, 0.d0)

            end if
                    
c     
c**** printout global mass and energy flows ****
c     
            if(ico2.ge.0.or.ice.ne.0) then
               if (ntty.eq.2) write(iout,784)
               if (iatty.gt.0) write(iatty,784)
            else if(ico2.lt.0) then
               if (ntty.eq.2) write(iout,796)
               if (iatty.gt.0) write(iatty,796)
            endif
 784        format(/,20x,'Global Mass & Energy Balances')
 796        format(/,20x,'Global Water & Gas Balances')
            if(ico2.ge.0.or.ice.ne.0) then
               if (ntty.eq.2) write(iout,785) amass,asteam,aener
               if (iatty.gt.0) write(iatty,785) amass,asteam,aener
            else if(ico2.lt.0) then
               if (ntty.eq.2) write(iout,797) amass,asteam,aener
               if (iatty.gt.0) write(iatty,797) amass,asteam,aener
            endif
 785        format(1x,'Total mass in system at this time:          ',
     &           e14.6,' kg',/,1x,'Total mass of steam in system at ', 
     &           'this time: ',e14.6,' kg',/,1x,
     &           'Total enthalpy in system at this time:      ',
     &           e14.6,' MJ')
 797        format(1x,'Total water in system at this time:         ',
     &           e14.6,' kg',/,1x,'Total mass of steam in system at ', 
     &           'this time: ',e14.6,' kg',/,1x,
     &           'Total Air(gas) in system at this time:      ',
     &           e14.6,' kg')
            if(ntty.eq.2) 
     &           write(iout,792) qtoti,qtoti/(day*8.64e4),curinflow,
     &           curinflow/(day*8.64e4),toutfl,toutfl/(days*8.64e4),
     &           totalflin,totalflin/(days*8.64e4)
            if (iatty.gt.0) write(iatty,792) qtoti,
     &           qtoti/(day*8.64e4),curinflow,
     &           curinflow/(day*8.64e4),toutfl,toutfl/(days*8.64e4),
     &           totalflin,totalflin/(days*8.64e4)
 792        format(/,1x,'Water discharge this time step: ',e14.6,' kg ',
     &           '(',e12.6,' kg/s)',/,1x,'Water input this time ',
     &           'step:     ',e14.6,' kg (',e12.6,' kg/s)',/,1x,
     &           'Total water discharge:          ',e14.6,' kg ',
     &           '(',e12.6,' kg/s)',/,1x,'Total water ',
     &           'input:              ',e14.6,' kg (',e12.6,' kg/s)')
            if(ico2.ge.0.or.ice.ne.0) then
               if (ntty.eq.2)
     &              write(iout,791) qtotei,qtotei/(day*8.64e4),
     &              cureinflow,cureinflow/(day*8.64e4),teoutf,
     &              teoutf/(days*8.64e4),totalein,totalein/(days*8.64e4)
               if (iatty.gt.0.) write(iatty,791) qtotei,
     &              qtotei/(day*8.64e4),cureinflow,
     &              cureinflow/(day*8.64e4),teoutf,teoutf/(days*8.64e4),
     &              totalein,totalein/(days*8.64e4)
            else if(ico2.lt.0) then
               if (ntty.eq.2)
     &              write(iout,793) qtotei,qtotei/(day*8.64e4),
     &              cureinflow,cureinflow/(day*8.64e4),teoutf,
     &              teoutf/(days*8.64e4),totalein,totalein/(days*8.64e4)
               if (iatty.gt.0.) write(iatty,793) qtotei,
     &              qtotei/(day*8.64e4),cureinflow,
     &              cureinflow/(day*8.64e4),teoutf,teoutf/(days*8.64e4),
     &              totalein,totalein/(days*8.64e4)
            endif
 791        format(/,1x,'Enthalpy discharge this time step: ',e13.6,
     &           ' MJ (',e12.6,' MJ/s)',/,1x,'Enthalpy input this ',
     &           'time step:     ',e13.6,' MJ (',e12.6,' MJ/s)',/,1x,
     &           'Total enthalpy discharge:          ',e13.6,
     &           ' MJ (',e12.6,' MJ/s)',/,1x,'Total enthalpy ',
     &           'input:              ',e13.6,' MJ (',e12.6,' MJ/s)')
 793        format(/,1x,'Air(gas) discharge this time step: ',e13.6,
     &           ' kg (',e12.6,' kg/s)',/,1x,'Air(gas) input this ',
     &           'time step:     ',e13.6,' kg (',e12.6,' kg/s)',/,1x,
     &           'Total air(gas) discharge:          ',e13.6,
     &           ' kg (',e12.6,' kg/s)',/,1x,'Total air(gas) ',
     &           'input:              ',e13.6,' kg (',e12.6,' kg/s)')
            if(ico2.ge.0.or.ice.ne.0) then
               if (ntty.eq.2) write(iout,787) qt,qte
               if (iatty.gt.0) write(iatty,787) qt,qte
            else if(ico2.lt.0) then
               if (ntty.eq.2) write(iout,794) qt,qte
               if (iatty.gt.0) write(iatty,794) qt,qte
            endif
 787        format(/,1x,'Net kg water discharge (total out-total in): '
     &           ,e14.6,/,1x,
     &           'Net MJ energy discharge (total out-total in): ',e14.6)
 794        format(/,1x,'Net kg water discharge (total out-total in): '
     &           ,e14.6,/,1x,
     &           'Net kg air discharge   (total out-total in): ',e14.6)
            if(ico2.ge.0.or.ice.ne.0) then
               if (ntty.eq.2) write(iout,788) difm,dife
               if (iatty.gt.0) write(iatty,788) difm,dife
            else if(ico2.lt.0) then
               if (ntty.eq.2) write(iout,795) difm,dife
               if (iatty.gt.0) write(iatty,795) difm,dife
            endif
 788        format(1x,'Conservation Errors: ',e14.6,' (mass), ',e14.6,
     &           ' (energy)')
 795        format(1x,'Conservation Errors: ',e14.6,' (water), ',e14.6,
     &           ' (air)')
            if(fdum.eq.-999.) then
c     write out flux discrepency (kg/s)
c  gaz 110519 changed g1 to tmch                 
	       if (ntty.eq.2) write(iout,979) abs(tmch),fdum1
               if (iatty.gt.0)	write(iatty,979) abs(tmch),fdum1	
 979           format(' >>> End on flowrate error (kg/s) ',1p,g12.4,
     &              ' Max error ',g12.4,' <<<')
c     now write in terms of (m**3/s), use density of water
               if (ntty.eq.2) write(iout,978) abs(tmch)/997.,fdum1/997.
              if (iatty.gt.0) write(iatty,978) abs(tmch)/997.,fdum1/997.	
 978          format(' >>> End on volume flux error (m**3/s) ',1p,g12.4,
     &              ' Max error ',g12.4,' <<<')
            endif    
         endif
         
         if (ntty.eq.1 .and. iptty.gt.0)  then
            if ( m2 .gt. 0 )  then
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
               if(idualp.ne.0) then
                  ilev=3
                  m2lev=m2/3
               else if(idpdp.ne.0) then
                  ilev=2
                  m2lev=m2/2
               else
                  ilev=1
                  m2lev=m2
               endif
               
               if(ihead.eq.0) then
                  if(iptty.gt.0) write(iptty ,6030)
               else
                  if(iptty.gt.0) write(iptty ,6130)
               endif
               do il=1,ilev
                  if(il.ne.1) then
                     if(iptty.gt.0) write(iptty,*) 'matrix level = ', il
                  endif
                  do i=1,m2lev
                     md=nskw2(i+(il-1)*m2lev)
                     if (irdof .ne. 13 .or. ifree .ne. 0) then
                        sl=s(md)
                     else
                        sl = 1.0d0
                     endif
                     sv=1.0-sl
                     if(ico2.lt.0.and.ice.eq.0) then
                        eqd=0.0
                     else
                        eqd=0.0
                        rhomd = sl * rolf(md) + sv * rovf(md)
                        hmd=sl*rolf(md)*enlf(md)+sv*rovf(md)*envf(md)
                        if (abs(rhomd) .gt. zero_t)  eqd = hmd/rhomd
                     endif
                     rqd=sk(md)
c     rqhd=0.0
c     if ( abs( rqd ) .gt. zero_t )  rqhd =  qh(md)/rqd
c     CHANGE ABOVE TO JUST PRINT OUT qh ARRAY

                     phod=pho(md)
                     if(ihead.ne.0) then
                        call headctr(4,md,pho(md),phod)
                        if(sl.lt.1.0) phod=0.0
                     endif
                     if(iptty.gt.0) write(iptty,6031) md,phod   ,eqd,
     *                    sl,t(md),rqd,qh(md)
                  enddo
               enddo
c     
c**** output for co2, air, or hydrate ****
c     
               if (ico2.gt.0) call  co2ctr  ( 5 )
               if (ico2.lt.0.and.ice.eq.0) then
                  call  airctr  ( 5,0 )
               else if (ico2.lt.0.and.ice.ne.0) then
                  call  icectr  ( 5,0 )
                  call  icectr  (-5,0 )
               endif
c     
c**** call varible porosity output ****
c     
c     dont call variable porosity for terminal monitor
c     
            endif
c     
c**** printout global mass and energy flows ****
c     
            if (iptty.gt.0) then
               if(ico2.ge.0.or.ice.ne.0) then
                  write(iptty,784)
                  write(iptty,785) amass,asteam,aener
               else if(ico2.lt.0) then
                  write(iptty,796)
                  write(iptty,797) amass,asteam,aener
               endif
               write(iptty,792) qtoti,qtoti/(day*8.64e4),curinflow,
     &              curinflow/(day*8.64e4),toutfl,toutfl/(days*8.64e4),
     &              totalflin,totalflin/(days*8.64e4)
               if(ico2.ge.0.or.ice.ne.0) then
                  write(iout,791) qtotei,qtotei/(day*8.64e4),cureinflow,
     &                 cureinflow/(day*8.64e4),teoutf,
     &                 teoutf/(days*8.64e4),
     &                 totalein,totalein/(days*8.64e4)
               else if(ico2.lt.0) then
                  write(iout,793) qtotei,qtotei/(day*8.64e4),cureinflow,
     &                 cureinflow/(day*8.64e4),teoutf,
     &                 teoutf/(days*8.64e4),
     &                 totalein,totalein/(days*8.64e4)
               endif
               if(ico2.ge.0.or.ice.ne.0) then
                  write(iptty,787) qt,qte
                  write(iptty,788) difm,dife
               else if(ico2.lt.0) then
                  write(iptty,794) qt,qte
                  write(iptty,795) difm,dife
               endif
            endif
            if (ntty.eq.2.or.ntty.eq.1) then
               write(iout,792) qtoti,qtoti/(day*8.64e4),curinflow,
     &              curinflow/(day*8.64e4),toutfl,toutfl/(days*8.64e4),
     &              totalflin,totalflin/(days*8.64e4)
               if(ico2.ge.0.or.ice.ne.0) then
                  write(iout,791) qtotei,qtotei/(day*8.64e4),cureinflow,
     &                 cureinflow/(day*8.64e4),teoutf,
     &                 teoutf/(days*8.64e4),
     &                 totalein,totalein/(days*8.64e4)
               else if(ico2.lt.0) then
                  write(iout,793) qtotei,qtotei/(day*8.64e4),cureinflow,
     &                 cureinflow/(day*8.64e4),teoutf,
     &                 teoutf/(days*8.64e4),
     &                 totalein,totalein/(days*8.64e4)
               endif
            endif
            
c     powavg =  qtote/(days*86400.0)
c     if(iptty.gt.0) write(iptty ,6042)  qtoti , qtotei , pow
c     if(iptty.gt.0) write(iptty ,6043)  qtot , qtote , powavg
c     if(iptty.gt.0) write(iptty ,6044)  ichng
         endif
      endif
      
 40   continue
      
      call   river_ctr (6)      
      call  concen  ( 2,0 )
      call  stress  ( 2 )
      call  paractr (3)

      call flush (iout)
      return
      end
      

      subroutine icectrco2(iflg,ndummy)
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
!D1 To manage the CO2(Supercritical/Liquid/Gas)-Water(Brine)-Air calculations
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 02/05/07 R. J. Pawar wrote the new subroutine 
!D2 capable of simulating air/water/CO2 flow taking into account dissolution 
!D2 of water in CO2 and CO2 in water. 
!D2 Refer to C:\rajesh\FEHM_CO2\FEHM_V230\combined_co2_omrwellb\co2_specific_files\equations_co2.doc 
!D2 for equations.
!**********************************************************************
!D3 
!D3 
!D3 
!D3 
!D3 
!D3 
!**********************************************************************


      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comki
      use comxi
      use davidi
      use comco2
      use comriv
      use commeth
      use comsplitts 

      implicit none

      integer iflg,ndummy,i,nr1,nr2,nr3,nr4,nr5,mi      
      integer i1,i2,i3,i4,i5,i6,ilev,mlev,il,md,k,j, kb, iphase
      integer ii,ij,icedc,iced,ieosd,icesd,idco2,duma,neqp1,ncont
      real*8, allocatable :: aiped(:)
      real*8, allocatable :: sktmp(:)
      real*8, allocatable :: esktmp(:)
      integer, allocatable :: iflg_flowmac(:)
      real*8 tliquid,dummyreal
      real*8 tw,pc,tc,ensrc,delsrc,phase_frac,skmd,qhmd
      real*8 skmd1,skmd10,skmd2,skmd3,amco2w,amco2w1,amco21
      real*8 amco2w0,skmd21,skmd31,amco2f,amco2f0,amco2f1
      real*8 dtps,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real*8 strd_co2,strd1
      real*8 eostol
      real*8 eosmg
      real*8 eosml
      real*8 deltc, delpc, tem(9)
      parameter(deltc=31.1d0,delpc=3.78d0)
      real*8 amaxflx, dums1,tolw, zero_e
      real*8 teinfl, flemax, amaxener
      real*8 fh, power, dpowerh, px,py,pz
      real*8 pcrit, tcrit,dumb
      real*8 psatd, tl1,tl2,dtps1,dtps2,dpsatt
      real*8 sl,sg,sw, pl, tl, tsolid
      real*8 denc,dencp,denct,enc,encp,enct,visc,viscp,visct
      real*8 denw,enwp,enwt,viswp,viswt
      real*8 ycmax, xwp, xcp, xcpb, mwc, mww, xc_prime, tmp1
      integer ico2d, ico2dc
      integer imped_ex_0
      parameter(pcrit=7.377d0, tcrit=30.98)
      character*8 macro1
      parameter(eostol=0.0001d00)
      parameter(eosmg=1.0001d00)
      parameter(eosml=0.9999d0)
      parameter(strd1=1.0d0)
      parameter(tolw=1.d-90, zero_e=1.e-10)
      parameter(mwc=44.d0,mww=18.d0)
      parameter(imped_ex_0 = 3)
      integer myfile, open_file
      save myfile
      logical it_is_open
      save amco2f, amco2f0, amco2f1,  amco2w, amco2w0, amco2w1, amco21
      save skmd1, skmd2, skmd21, skmd3, skmd31
      if(nr_stop.eq.0) then
       strd_co2=0.7d0
      else
       strd_co2 = strd_iter
      endif

      if(icarb.ne.0) then  

         if(iflg.eq.0) then

c     RJP 04/15/07 Changed the dof definitions
            idof=5
            if(iprtype.eq.1) then
               
               idof_co2 = 1
               if (iout .ne. 0) write(iout, 20) 
               if (iptty .gt. 0) write(iptty, 20)  
            else if(iprtype.eq.2) then

               idof_co2 = 2
               if (iout .ne. 0) write(iout, 21) 
               if (iptty .gt. 0) write(iptty, 21)  
            else if(iprtype.eq.3) then
               
               idof_co2 = 3
               if (iout .ne. 0) write(iout, 22) 
               if (iptty .gt. 0) write(iptty, 22)  
            else if(iprtype.eq.-3) then

               idof_co2 = 3
               if (iout .ne. 0) write(iout, 222) 
               if (iptty .gt. 0) write(iptty, 222)  

            else if(iprtype.eq.4) then

               idof_co2 = 3
               if (iout .ne. 0) write(iout, 23) 
               if (iptty .gt. 0) write(iptty, 23)  
            else if(iprtype.eq.5) then

               idof_co2 = 4
               if (iout .ne. 0) write(iout, 24) 
               if (iptty .gt. 0) write(iptty, 24)  
            else 
               if (iout .ne. 0) write(iout, *) 
     &              'Stopping : iprtype not valid' 
               if (iptty .gt. 0) write(iptty, *)
     &              'Stopping : iprtype not valid'  
            endif

 20         format(3x,'**** water only problem: ',
     &           'P and T (1-phase Water) ****') 
 21         format(3x,'**** CO2 only problem: ',
     &           'P and T (1-phase CO2) ****') 
 22         format(3x,'**** water-CO2 problem: ',
     &           'P, T and fw (2-phase CO2, with-out dissolution) ****')
 222        format(3x,'**** water-CO2 isothermal problem: ',
     &           'P and fw (water-CO2, with-out dissolution) ****') 
 23         format(3x,'**** water-CO2 problem: ',
     &           'P, T and fw (2-phase CO2, with dissolution) ****') 
 24         format(3x,'**** water-CO2-air problem: ',
     &           'P, T, fw, xco2 (2-phase CO2, dissolution)  ****') 
c     
c     allocate memory for co2 arrays
c     
            allocate(denco2i(n0))
            allocate(denco2h(n0))
            allocate(deneco2i(n0))
            allocate(deneco2h(n0))
            allocate(pflowco2(n0))
            allocate(flowco2s(n0))
            allocate(eflowco2(n0))
            allocate(phico2(n0))
            allocate(phoco2(n0))
            allocate(tco2(n0))
            allocate(toco2(n0))
            allocate(fw(n0))
            allocate(fow(n0))
            allocate(fg(n0))
            allocate(fog(n0))
            allocate(fl(n0))
            allocate(fol(n0))
            allocate(yc(n0))
            allocate(yw(n0))
            allocate(xc(n0))
            allocate(xw(n0))
            allocate(yoc(n0))
            allocate(yow(n0))
            allocate(xoc(n0))
            allocate(xow(n0))
            allocate(ya(n0))
            allocate(xa(n0))
            allocate(skco2(n0))
            allocate(eskco2(n0))
            allocate(qhco2(n0))
            allocate(wellco2(n0))
            allocate(dmpf(n0))
            allocate(dmef(n0))
            allocate(dmwf(n0))
            allocate(dmycf(n0))
            allocate(dmyaf(n0))
            allocate(depf(n0))
            allocate(deef(n0))
            allocate(dewf(n0))
            allocate(deycf(n0))
            allocate(deyaf(n0))
            allocate(diw(n0))
            allocate(diwc(n0))
            allocate(diwp(n0))
            allocate(diwe(n0))
            allocate(diww(n0))
            allocate(diwyc(n0))
            allocate(diwya(n0))
            allocate(dilp(n0))
            allocate(dile(n0))
            allocate(dilw(n0))
            allocate(dilyc(n0))
            allocate(dilya(n0))
            allocate(divp(n0))
            allocate(dive(n0))
            allocate(divw(n0))
            allocate(divyc(n0))
            allocate(divya(n0))

            allocate(dqw(n0))
            allocate(dqyc(n0))
            allocate(dqya(n0))
            allocate(dqh(n0))
            allocate(deqh(n0))
            allocate(dqhw(n0))
            allocate(dqhyc(n0))
            allocate(dqhya(n0))
            allocate(deqw(n0))

            allocate(qhflxco2(n0))
            allocate(dtpaco2(n0))
            allocate(dtpaeco2(n0))
            allocate(kaco2(n0))
            allocate(iceso(n0))
            allocate(ieoso(n0))
            allocate(dpcpw(n0))
            allocate(dpcpg(n0))
            allocate(pcg(n0))
            allocate(dpcgg(n0))
            allocate(dpcgw(n0))
            allocate(diltrac(n0))
            allocate(csalt(n0))
            allocate(stowat(n0*2))
            allocate(dtpsc(n0))
            
            allocate(dskco2w(n0))
            allocate(dskco2yc(n0))
            allocate(deqyc(n0))
            allocate(deqya(n0))

            allocate(wat_prop(n0*15))
            allocate(co2_prop(n0*18))
            allocate(dmol(n0*14))
            allocate(diff(n0))
            allocate(tortco2(n0))
            allocate(ico2dis(n0))
            allocate(ico2diso(n0))

            allocate(fw_tmp(n0))
            allocate(fl_tmp(n0))
            allocate(fg_tmp(n0))
            allocate(inico2flg(n0))

            allocate(strd_arr(n0))
c     
c     zero out arrays for co2        
c     
            stowat = 0.d0
            denco2i = 0.0d00
            denco2h = 0.0d00
            deneco2i = 0.0d00
            deneco2h = 0.0d00
            pcp = 0.0d00
            pflowco2 = 0.0d00
            phico2 = 0.0d00
            phoco2 = 0.0d00
            fg = 0.0d00
            fog = 0.0d00
            fl = 0.0d00
            fol = 0.0d00
            tco2 = 0.0d00
            toco2 = 0.0d00
            skco2 = 0.0d00
            qhco2 = 0.0d00
            eskco2 = 0.0d00
            eflowco2 = 0.0d00
c            rovfco2 = 0.0d00
c            envfco2 = 0.0d00
            wellco2 = 0.0d00
c            rolfco2 = 0.0d00
c            enlfco2 = 0.0d00
c            dpcpco22 = 0.0d00
c            dmmf = 0.0d00
c            demf = 0.0d00
c            dilm = 0.0d00
c            divm = 0.0d00
c            deqm = 0.0d00
c            deqpm = 0.0d00
            dmwf = 0.0d00
            dewf = 0.0d00
            dilw = 0.0d00
            divw = 0.0d00
            deqw = 0.0d00
c            deqpw = 0.0d00
            fw = 0.0d00
            fow = 0.0d00
            qhflxco2 = 0.0d00
            dtpaco2 = 0.0d00
            dtpaeco2 = 0.0d00
            skmd1 = 0.0d00
c            dq3 = 0.d00
c            dq4 = 0.d00
            csalt = 0.d00
            co2_prop = 0.d00
            wat_prop = 0.d00
            diltrac = 0.d00
c           dpcp3 = 0.00
c           dpcp4 = 0.00
            dtpsc = 0.00
            kaco2 = 0
            iceso = 0
            ieoso = 0
            idco2 = 0
            diff = 0.d0
            tortco2 = 0.d0
            ico2diff_flg = 0
            ico2prop_flg = 0
            ico2dis = 0
            imped_ex = 10000000
            yoc = 0.0
            xoc = 0.0
            yow = 0.0
            xow = 0.0
		  dmol = 0.0
            inico2flg = 0
            macro = "carb "
c     
c     read in "sub macros" for CO2
c     

            
 100        continue
            read (inpt, '(a80)') wdd1
c Changed end check for compatibility with macro "off" option
            if (wdd1(1:6) .eq. 'co2end' .or. wdd1(1:7) .eq. 'endcarb'
     &           .or. wdd1(1:8) .eq. 'end carb') go to 200 
            if (wdd1(1:1) .eq. '#') go to 100 
            read (wdd1, '(a8)') macro1
            if (iout .ne. 0) write(iout, 50) macro1
            if (iptty .gt. 0) write(iptty, 50) macro1 
 50         format(3x, '**** CO2 sub macro : ', a8,' **** ' ) 
            if(macro1.eq.'co2pres') then
c     
c     read in initial pressure, temperature & CO2 phase state
c     
               allocate(sktmp(n0))

               igroup = 1
               narrays = 3
               itype(1) = 8
               itype(2) = 8
               itype(3) = 4
               default(1) = 0.
               default(2) = 0.
               default(3) = 0
               call initdata2( inpt, ischk, n0, narrays, itype, 
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1=phico2(1:n0),r8_2=sktmp(1:n0),i4_1=ices(1:n0))
               
               do i = 1, n0
                  if (sktmp(i) .ne. default(2)) then
                     if (abs (ices(i)) .eq. 1) then
                        tco2(i) =  sktmp(i)
                     else if (abs (ices(i)) .eq. 2) then
                        pl=phico2(i)
                        iced = ices(i)
                        call co2_properties(2,iced,pl,dumb,dum1,duma,
     &                       tliquid,dumc,duma)
                        tco2(i)=tliquid
                     else if (abs (ices(i)) .eq. 3) then
                        tco2(i) =  sktmp(i)
                     else if(abs(ices(i)).eq.4) then
                        tco2(i) =  sktmp(i)
                     end if
                  end if
               end do
               
               deallocate(sktmp)
            elseif (macro1.eq.'co2diff') then
               ico2diff_flg = 1
               
               igroup = 4
               narrays = 2
               itype(1) = 8
               itype(2) = 8
               default(1) = 0.d00
               default(2) = 0.d00
c     
c     read in CO2 diffusivity in water
c     

               call initdata2( inpt, ischk, n0, narrays, itype, default,
     &              macroread(8), macro, igroup, ireturn,
     &              r8_1 = diff(1:n0), r8_2 = tortco2(1:n0))

            elseif (macro1.eq.'userprop') then
               
               ico2prop_flg = 1
               allocate(con_prop(18))
c     read user defined properties for CO2 & brine
c     first read CO2 properties with derivatives order is
c     density, density der with p, density der with t
c     enthalpy, enthalpy der with p, enthalpy der with t
c     viscosity, viscosity der with p, viscosity der with t
               read(inpt, *) denc,dencp,denct,enc,encp,enct,visc,viscp,
     &              visct	     
               con_prop(1) = denc
               con_prop(2) = dencp
               con_prop(3) = denct
               con_prop(4) = enc
               con_prop(5) = encp
               con_prop(6) = enct
               con_prop(7) = visc
               con_prop(8) = viscp
               con_prop(9) = visct

c     Next read brine properties with derivatives
               read(inpt, *) denw,denwp,denwt,enw,enwp,enwt,visw,viswp,
     &              viswt
               con_prop(10) = denw
               con_prop(11) = denwp
               con_prop(12) = denwt
               con_prop(13) = enw
               con_prop(14) = enwp
               con_prop(15) = enwt
               con_prop(16) = visw
               con_prop(17) = viswp
               con_prop(18) = viswt
            else if(macro1.eq.'co2frac') then
               
c     
c     read in initial CO2, air, water/brine saturation. fw represents water-rich
c     liquid saturation (volume fraction), fg represents CO2/air-rich gas
c     saturation (volume fraction), fl represent CO2-rich super-critical/liquid  
c     phase saturation (volume fraction). In addition, brine salt concentration 
c     is also specified in ppm.
c     RJP 01/11/08 modified following to specify initial CO2 mass fraction in
c     CO2 dominant phase, xc.
c     RJP 03/09/08 removed qhflxco2 present in earlier versions.
               igroup = 1
               narrays = 5
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               itype(4) = 8
               itype(5) = 4
               default(1) = 0.
               default(2) = 0.
               default(3) = 0.
               default(4) = 0.
               default(5) = 0
               call initdata2( inpt, ischk, n0, narrays,itype, default,
     &              macroread(8), macro, igroup, ireturn,
     &              r8_1 = fw(1:n0),r8_2 = fl(1:n0),
     &              r8_3 = yc(1:n0),r8_4 = csalt(1:n0),
     &              i4_1 = inico2flg(1:n0))
               
               do i = 1, n0
                  fg(i) = 1.d0-fw(i)-fl(i)
                  if(inico2flg(i).eq.1) then
                     fw_tmp(i) = fw(i)
                     fl_tmp(i) = fl(i)
                     fg_tmp(i) = fg(i)
                  endif
               enddo


            else if (macro1.eq.'brine') then

               ibrine = 1
               if (iout .ne. 0) write(iout, 25) 
               if (iptty .gt. 0) write(iptty, 25)  
 25            format(3x,'**** brine option invoked  ****') 
            else if(macro1.eq.'co2flow') then
               igroup = 4
               narrays = 4
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               itype(4) = 4
               default(1) = 0.0d00
               default(2) = 0.0d00
               default(3) = 0.0d00
               default(4) = 0

c     
c     read in initial CO2 flow and boundary data
c     
               allocate(aiped(n0),esktmp(n0),sktmp(n0),iflg_flowmac(n0))
               iflg_flowmac = 0
               call initdata2( inpt, ischk, n0, narrays, itype, default,
     &              macroread(8), macro, igroup, ireturn,  r8_1 = 
     &              sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0),
     &              i4_1 = iflg_flowmac)
               do i=1,n0
                  if(iflg_flowmac(i).eq.1) then
c     constant pressure boundary condition with inflow allowed. aiped is user specified
                     kaco2(i) = -1
                     pflowco2(i) = sktmp(i)
                     eskco2(i) = esktmp(i)
                     wellco2(i) = abs(aiped(i)) * 1.e+6
                  elseif(iflg_flowmac(i).eq.2) then
c     constant pressure boundary condition with only outflow allowed. aiped is user specified
                     kaco2(i) = -2
                     pflowco2(i) = sktmp(i)
                     eskco2(i) = esktmp(i)
                     wellco2(i) = abs(aiped(i)) * 1.e+6
                  elseif(iflg_flowmac(i).eq.3) then
c     constant pressure boundary condition. aiped is calculated in the code based on block geometric params
                     kaco2(i) = -3
                     pflowco2(i) = sktmp(i)		
                     eskco2(i) = esktmp(i)
                  elseif(iflg_flowmac(i).eq.4) then
c     constant pressure and constant saturation boundary condition. aiped is user specified
                     kaco2(i) = -4
                     pflowco2(i) = sktmp(i)
                     flowco2s(i) = esktmp(i)
                     wellco2(i) = abs(aiped(i)) * 1.e+6
                  elseif(iflg_flowmac(i).eq.5) then
c     constant pressure and constant saturation boundary condition. aiped is calculated in the code based
c     on block geometric params
                     kaco2(i) = -5
                     pflowco2(i) = sktmp(i)
                     flowco2s(i) = esktmp(i)
                  elseif(iflg_flowmac(i).eq.6) then
c     constant co2 mass flow rate boundary condition
                     kaco2(i) = 1
                     skco2(i) = sktmp(i)
                     eskco2(i) = esktmp(i)
                  elseif(iflg_flowmac(i).eq.7) then
c     constant dissolved co2 mass flow rate boundary condition
                     kaco2(i) = 2
                     skco2(i) = sktmp(i)
                     eskco2(i) = esktmp(i)
                     flowco2s(i) = aiped(i)
                  elseif(iflg_flowmac(i).eq.8) then
c     constant dissolved co2 mass flow rate boundary condition
                     kaco2(i) = 3
                     flowco2s(i) = sktmp(i)
                     eskco2(i) = esktmp(i)
                     sk(i) = aiped(i)
                  elseif(iflg_flowmac(i).eq.9) then
c     partial explicit update of nonliniar part of co2 constant pressure
                     kaco2(i) = -1
                     imped_ex = imped_ex_0
                     kaco2(i) = -1
                     pflowco2(i) = sktmp(i)
                     eskco2(i) = esktmp(i)
                     wellco2(i) = abs(aiped(i)) * 1.e+6                  
                  endif
               enddo
               if(imped_ex.eq.imped_ex_0) then
                allocate(permsd1_sv_w(n0))
                allocate(permsd1_sv_co2(n0))
               endif
               deallocate(sktmp,esktmp,aiped,iflg_flowmac)
               
            else
               if (iout .ne. 0) write(iout,*) 
     &              'ERROR IN CO2 INPUT(STOPPING)'
               write(ierr,*) 'ERROR IN CO2 INPUT(STOPPING)'
            end if
 40         continue
            go to 100
 200        continue
            
            do i=1,n0
               phoco2(i) = phico2(i)
               toco2(i) = tco2(i)
               fow(i)=fw(i)
               fog(i)=fg(i)
               fol(i)=fl(i)
               iceso(i)  = ices(i)
            enddo
            
            macroread(8) = .TRUE.
            
         elseif(iflg.eq.-2) then
c     
c     set tempertures of components equal for certain equilibrium conditions
c     
            if (.not. co2_read) then
               do i=1,n0
                  to(i) = tco2(i)
                  t(i) = tco2(i)
               enddo         
c     set phase of water to 1 (liquid)
               do i=1,n0
                  ieos(i) = 1
                  ieoso(i)  = ieos(i)
               enddo
            end if
c
c   set water BC (constant pressure) in CO2 BC set
c
c          go to 1599
               do i=1,n0
                  if(kaco2(i).lt.0.and.ka(i).eq.0) then
                   ka(i) = -1 
                   wellim(i) = 1.e-13
                   pflow(i) = pflowco2(i)
                  endif
               enddo
c1599     continue
c
         else if(iflg.eq.-1) then
c     
c     determine phase state for water solid-liquid-gas system
c     should count phase changes
c     h20 phase changes not allowed for co2
          if(itsat.le.10.and.icarb.eq.0) then           
            do ii=1,neq
               ij=ii+ndummy
               iced=ieos(ij)
               icedc=iced
               pl=phi(ij)
               if(iced.eq.3) then
c     gas only conditions
c     tl=t(ij)
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  if(tl.le.tliquid) then
c     gas-liquid can form
                     icedc=2
                     s(ij)=eostol
                     t(ij)=tliquid
                  endif
               elseif(iced.eq.2) then
c     gas-liquid conditions
                  sl=s(ij)
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  t(ij)=tliquid
                  if(sl.ge.1.0) then
c     liquid only can form                  
                     icedc=1
                     s(ij)=1.0
                     t(ij)=tliquid*eosmg
                  else if(sl.le.0.0) then
c     gas only can form    
                     icedc=3
                     s(ij)=1.0
                     t(ij)=tliquid*eosml
                  endif
               elseif(iced.eq.1) then
c     liquid only conditions
                  tl=t(ij)
                  call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  if(tl.ge.tliquid) then
c     gas-liquid can form                
                     icedc=2
                     s(ij)=eosml
                     t(ij)=tliquid
                  else if(tl.le.tsolid) then
c     liquid-solid can form                
                     icedc=-2
                     s(ij)=eostol
                     t(ij)=tsolid
                  endif
               elseif(iced.eq.-2) then
c     liquid-solid conditions
                  sl=s(ij)
                  call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  if(sl.ge.1.0) then
c     liquid only can form                
                     icedc=1
                     s(ij)=1.0   
                     t(ij)=tsolid*eosmg
                  else if(sl.le.0.0) then
c     solid only can form                
                     icedc=-3
                     s(ij)=0.0    
                     t(ij)=tsolid*eosml
                  endif
               elseif(iced.eq.-3) then
c     solid only conditions
                  tl=t(ij)
                  call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  if(tl.ge.tsolid) then
c     liquid-solid can form                
                     icedc=-2
                     s(ij)=eosml 
                     t(ij)=tsolid
                  endif
               endif
               if(icedc.ne.iced) strd = strd_co2
               ieos(ij)=icedc
            enddo
          else if(itsat.gt.10) then
            ieos = 1
          endif
c     
         else if(iflg.eq.1) then
c     
c     determine phase state for co2 supercritical-liquid-gas system
c     should count phase changes
c     
            
            do ii=1,neq
               ij=ii+ndummy
               iced=ices(ij)
               icedc=iced
               ico2d=ico2dis(ij)
               ico2dc=ico2d
               pl=phico2(ij)
               tl=tco2(ij)
               if(iced.eq.4) then
                  call co2_properties(1,iced,pl,tl,dum1,icedc,dumb,dumc,
     &			duma)
                  if(icedc.eq.3) then
                     tco2(ij)=tl
                     fg(ij) = fl(ij)
                     fl(ij) = 0.d0
                  elseif(icedc.eq.1) then
                     tco2(ij)=tl
                  endif
               elseif(iced.eq.3) then
                  call co2_properties(1,iced,pl,tl,dum1,icedc,dumb,dumc,
     &			duma)
                  if(icedc.eq.4) then
                     fl(ij) = fg(ij)
                     fg(ij) = 0.d0
                  elseif(tl.lt.tcrit) then
                     if(icedc.eq.1) then
                        if(pl.lt.pcrit) then
                           call co2_properties(2,iced,pl,dumb,dum1,duma,
     &                          tliquid,dumc,duma)
                           if(fw(ij).lt.1.d0) then
                              icedc=2
                              fg(ij) = fg(ij)*0.9d0
                              fl(ij)=1.d0-fg(ij)-fw(ij)
                              tco2(ij)=tliquid
                           else
                              fg(ij) = 0.d0
                              fl(ij)=1.d0-fg(ij)-fw(ij)
                           endif
                        endif
                     endif
                  endif
               elseif(iced.eq.2) then
c     gas-liquid conditions
                  sg=fg(ij)
                  sl=fl(ij)
                  call co2_properties(2,iced,pl,dumb,dum1,duma,tliquid,
     &                 dumc,duma)
                  tco2(ij) = tliquid
                  if(sg.le.0.0) then
c     liquid only can form                  
                     icedc=1
                     fg(ij) = 0.d0
                     fl(ij)=1.d0-fg(ij)-fw(ij)
c                     tco2(ij)=tliquid*eosml
                  else if(sl.le.0.0) then
c     gas only can form
c  gaz debug 040115 
                     icedc=3
                     fl(ij) = 0.d0
                     fg(ij) = 1. - fw(ij)
c                     tco2(ij)=tliquid*eosmg
                  endif
                  t(ij)=tco2(ij)
               elseif(iced.eq.1) then
c     liquid only conditions
                  call co2_properties(1,iced,pl,tl,dum1,icedc,dumb,dumc,
     &			duma)
                  if(pl.lt.pcrit) then
                     if(icedc.eq.3) then
                        if(tl.lt.tcrit) then
                           call co2_properties(2,iced,pl,dumb,dum1,duma,
     &                          tliquid,dumc,duma)
                           if(fw(ij).lt.1.d0) then
                              icedc=2
                              fl(ij)=fl(ij)*0.99d0
                              fg(ij)=1.d0-fl(ij)-fw(ij)
c gaz debug 040115
c                              tco2(ij)=tliquid
                           else
                              fl(ij) = 0.d0
                              fg(ij)=1.d0-fl(ij)-fw(ij)
                           endif
                        endif
                     endif
                  endif
               endif
c     check if dissolved co2 is reached max.
c     if water dissolving in co2 is activated, calculate partial pressure
c     of co2
c RJP added below 03/14/2011. This should be activated only for 
c carb option '4' where dissolution comes into play.
               ices(ij)=icedc

			call ther_co2_h2o(11,ij)
			if(iprtype.ge.4) then			
               if(iwatdis.eq.1) then
c     calculate phase-change pressure and dp/dt
                  call h2o_properties(5,2,tl,dum1,dum2,dum3,
     &                 psatd,dum4,dpsatt,dum5,dum6)
                  xwp = psatd/pl
                  xcp = 1.d0-xwp		  
               else
                  xwp = 0.d0
                  xcp = 1.d0-xwp		  
               endif
               xcpb=(xcp/mwc)+(xwp/mww)
               xc_prime = xcp/(mwc*xcpb)
               call co2_properties(10,1,phico2(ij),tco2(ij),csalt(ij),
     &              1,xc_prime,tem,ij)
               ycmax = tem(1)
               if(ico2d.eq.0) then
                  if(yc(ij).ge.ycmax) then
                     ico2dis(ij)=1
				   ico2dc=ico2dis(ij)
                     yc(ij) = ycmax
                     yw(ij) = 1.d0 - yc(ij)
                     if(icedc.eq.2) then
                        fg(ij) = 0.000001d0
                        fw(ij) = fw(ij) - 0.000001d0
                        fl(ij) = 1.d0 - fw(ij) - fg(ij)
                     elseif(icedc.eq.3) then
                        fg(ij) = 0.000001d0
                        fw(ij) = fw(ij) - 0.000001d0
                     else
                        fl(ij) = 0.000001d0
                        fw(ij) = fw(ij) - 0.000001d0
                     endif
                     dmol(ij) = tem(2)
                     dmol(ij+neq) = tem(3)
                  endif
               else
                  if(yc(ij).lt.ycmax) then
                     if(fw(ij).lt.1.d0) then
                        yc(ij) = ycmax
                        yw(ij) = 1.d0 - yc(ij)
                        dmol(ij) = tem(2)
                        dmol(ij+neq) = tem(3)
                     else
                        ico2dis(ij) = 0
				   ico2dc=ico2dis(ij)
c     fw(ij)= 1.d0
c     fl(ij)= 0.d0
                        yc(ij) = ycmax*0.999999d0
                        yw(ij) = 1.d0 - yc(ij)
                     endif
                  else
                     yc(ij) = ycmax
                     yw(ij) = 1.d0 - yc(ij)
                     dmol(ij) = tem(2)
                     dmol(ij+neq) = tem(3)
                  endif
               endif
			 endif
               if((icedc.ne.iced).or.(ico2dc.ne.ico2d)) then
                  strd = strd_co2
               else
                  strd = strd1
               endif
               strd_arr(ij) = strd
            enddo
c     check initial phase state and saturation, this is only called in 
c     startup
         else if(iflg.eq.14) then
c     
c     determine phase state for co2 supercritical-liquid-gas system
c     should count phase changes
c     
            
            do ii=1,neq
               ij=ii+ndummy
               pl=phico2(ij)
               tl=tco2(ij)
               iced = ices(ij)
               if(iced.ne.2) then
                  call co2_properties(1,iced,pl,tl,dum1,icedc,dumb,dumc,
     &			duma)
                  ices(ij)=icedc
               else
                  call co2_properties(2,iced,pl,dumb,dum1,duma,
     &                 tliquid,dumc,duma)
                  tco2(ij)=tliquid
               endif
c     
c     check if dissolved co2 is reached max.
c     if water dissolving in co2 is activated, calculate partial pressure
c     of co2
            enddo
         elseif(iflg.eq.12) then
c     RJP 11/09/06 Introduced new state check routine for initialization
            do ii=1, neq
               pl = phico2(ii)
               tl = tco2(ii)
               call co2_properties(1,duma,pl,tl,dum1,ices(ii),dumb,dumc,
     &		 duma)
            enddo
c     
         elseif(iflg.eq.2) then
c     
c     update variables with N-R corrections           
c     
c     
            if(nr_stop.ne.0.and.iad.ge.1)then
             call nr_stop_ctr(1)
            endif
            if(iprtype.eq.1) then
c     water-only problem
               nr1=nrhs(1)
               nr2=nrhs(2)
               do i = 1,neq
                  i1=i+nr1
                  i2=i+nr2
                  if(ps(i).eq.0.0) then
                     tco2(i)=tco2(i)-bp(i2)*strd
                  else
                     phico2(i) = phico2(i)-bp(i1)*strd
                     tco2(i) = tco2(i)-bp(i2)*strd
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  endif
                  t(i) = tco2(i)
               enddo
            elseif(iprtype.eq.2) then
c     co2-only problem, includes phase-change
               nr1=nrhs(1)
               nr2=nrhs(2)
               do i = 1,neq
                  i1=i+nr1
                  i2=i+nr2
                  icesd=ices(i)
                  if(ps(i).eq.0.0) then
                     tco2(i)=tco2(i)-bp(i2)*strd
                  elseif(icesd.ne.2)then
                     phico2(i) = phico2(i)-bp(i1)*strd
                     tco2(i) = tco2(i)-bp(i2)*strd
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  else
                     phico2(i)=phico2(i)-bp(i1)*strd
                     fg(i) = fg(i) - bp(i2)*strd
                     fl(i) = 1.d0-fg(i)-fw(i)
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  endif
                  t(i) = tco2(i)
               enddo
            elseif(iprtype.eq.3) then
c     co2-water problem, includes phase-change, no CO2 dissolution
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
               do i = 1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  icesd=ices(i)
                  if(ps(i).eq.0.0) then
                     tco2(i)=tco2(i)-bp(i2)*strd
                  elseif(icesd.ne.2)then
                     phico2(i) = phico2(i)-bp(i1)*strd
                     tco2(i) = tco2(i)-bp(i2)*strd
                     fw(i) = fw(i) - bp(i3)*strd
                     if (icesd.eq.3) then
                        fg(i) = 1.d0-fw(i)
                        fl(i) = 0.0
                     else
                        fl(i) = 1.d0-fw(i)
                        fg(i) = 0.0
                     endif
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  else
                     phico2(i)=phico2(i)-bp(i1)*strd
                     fg(i) = fg(i) - bp(i2)*strd
                     fw(i) = fw(i) - bp(i3)*strd
                     fl(i) = 1.d0-fg(i)-fw(i)
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  endif
                  t(i) = tco2(i)
               enddo
            elseif(iprtype.eq.-3) then
c     Isothermal co2-water problem no CO2 dissolution
               nr1=nrhs(1)
               nr2=nrhs(2)
               do i = 1,neq
                  i1=i+nr1
                  i2=i+nr2
                  if(ps(i).eq.0.0) then
                     tco2(i)=tco2(i)
                  else
                     phico2(i)=phico2(i)-bp(i1)*strd
                     fw(i) = fw(i) - bp(i2)*strd
                     if (icesd.eq.3) then
                        fg(i) = 1.d0-fw(i)
                     else
                        fl(i) = 1.d0-fw(i)
                     endif
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  endif
                  tco2(i) = tco2(i)
                  t(i) = tco2(i)
               enddo
            elseif(iprtype.eq.4) then
c     co2-water problem, includes phase-change,  CO2 dissolution
c     RJP 02/24/08 added a new array 'ico2ds' to track co2 dissolution
c     state.
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
               do i = 1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  icesd=ices(i)
                  ico2d=ico2dis(i)
                  strd = strd_arr(i)
                  if(ps(i).eq.0.0) then
                     tco2(i)=tco2(i)-bp(i2)*strd
                  elseif(icesd.ne.2)then
                     phico2(i) = phico2(i)-bp(i1)*strd
                     tco2(i) = tco2(i)-bp(i2)*strd
                     if(ico2d.eq.1) then
                        fw(i) = fw(i) - bp(i3)*strd
                        if (icesd.eq.3) then
                           fg(i) = 1.d0-fw(i)
                        else
                           fl(i) = 1.d0-fw(i)
                        endif
                     else
                        yc(i) = yc(i) - bp(i3)*strd
                     endif
c                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
c     &                    phico2(i) = 0.1
                  else
                     phico2(i)=phico2(i)-bp(i1)*strd
                     fg(i) = fg(i) - bp(i2)*strd
                     if(ico2d.eq.1) then
                        fw(i) = fw(i) - bp(i3)*strd
c     if (icesd.eq.3) then
c     fg(i) = 1.d0-fw(i)
c     else
c     fl(i) = 1.d0-fw(i)
c     endif
                     else
                        yc(i) = yc(i) - bp(i3)*strd
                     endif
                     fl(i) = 1.d0-fg(i)-fw(i)
c                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
c     &                    phico2(i) = 0.1
                  endif
                  t(i) = tco2(i)
                  if(iwatdis.eq.1) then
c     calculate phase-change pressure and dp/dt
                     call h2o_properties(5,2,tco2(i),dum1,dum2,dum3,
     &                    psatd,dum4,dpsatt,dum5,dum6)
                     xwp = psatd/phico2(i)
                     xcp = 1.d0-xwp		  
                  else
                     xwp = 0.d0
                     xcp = 1.d0-xwp		  
                  endif
                  xcpb=(xcp/mwc)+(xwp/mww)
                  xc_prime = xcp/(mwc*xcpb)
                  xc(i) = xc_prime
                  xw(i) = 1.d0 - xc(i)
                  yw(i) = 1.d0 - yc(i)
CHari 
                  s(i) = fw(i)
               enddo
            elseif(iprtype.eq.5) then
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
               nr4=nrhs(4)
               do i = 1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  i4=i+nr4
                  icesd=ices(i)
                  if(ps(i).eq.0.0) then
                     tco2(i)=tco2(i)-bp(i2)*strd
                  elseif(icesd.ne.2)then
                     phico2(i) = phico2(i)-bp(i1)*strd
                     tco2(i) = tco2(i)-bp(i2)*strd
                     fw(i) = fw(i) - bp(i3)*strd
                     if(icesd.eq.3) then
                        fg(i) = 1.d0-fw(i)
                     else
                        fl(i) = 1.d0-fw(i)
                     endif
                     xc(i) = xc(i) - bp(i4)*strd
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  else
                     phico2(i)=phico2(i)-bp(i1)*strd
                     fg(i) = fg(i) - bp(i2)*strd
                     fw(i) = fw(i) - bp(i3)*strd
                     fl(i) = 1.d0-fg(i)-fw(i)
                     xc(i) = xc(i) - bp(i4)*strd
                     if (phico2(i) .ge. 0. .and. phico2(i) .lt. 0.1)
     &                    phico2(i) = 0.1
                  endif
                  t(i) = tco2(i)
               enddo
            endif
         elseif(iflg.eq.-34) then
c     call EOS routine to allocate space
            call ther_co2_h2o(0,ndummy)
         elseif(iflg.eq.-35) then
c     call EOS routine to deallocate space
            call ther_co2_h2o(-1,ndummy)
         elseif(iflg.eq.-3) then
c     call EOS routines co2
            call ther_co2_h2o(2,ndummy)
         elseif(iflg.eq.3) then
c     call EOS routines water
            call ther_co2_h2o(1,ndummy)
         elseif(iflg.eq.-4) then
c     first ckeck for temperature-specified source terms
c     need to comment out call in main routine for water if meth is activated
c           if(iad.eq.0) then
c     only for first iteration
               do ii=1,neq
                  ij=ii+ndummy
                  if(eskco2(ij).lt.0.0) then
                     if(ps(ij).eq.0.d0) then
                        eflowco2(ij) = cpr(ij)*(-eskco2(ij))
                     else
                        if(abs(kaco2(ij)).gt.0) then
                           if((kaco2(ij).eq.2).or.(kaco2(ij).eq.3)) then
                              call co2_properties(4,ices(ij),
     &                             flowco2s(ij),-eskco2(ij),dum1,duma,
     &							 dumb,dumc,duma)
	                        eflowco2(ij)=dumc(4)
                              call h2o_properties(9,1,phi(ij),
     &                             -eskco2(ij),csalt(ij),
     &                             dum3,ensrc,dum4,dum5,dum6,dum7)
                              eflow(ij)=ensrc
                              if(iwatdis.eq.1) then
c     calculate phase-change pressure and dp/dt
                                 call h2o_properties(5,2,-eskco2(ij),
     &                           dum1,dum2,dum3,psatd,dum4,dpsatt,dum5,
     &							dum6)
                                 xwp = psatd/phico2(ij)
                              else
                                 xwp = 0.d0
                              endif
                              xcp = 1.d0-xwp		  
                              xcpb=(xcp/mwc)+(xwp/mww)
                              xc_prime = xcp/(mwc*xcpb)
							if(kaco2(ij).eq.2) then
                              call co2_properties(10,1,phico2(ij),
     &                             -eskco2(ij),csalt(ij),1,xc_prime,tem,
     &							 ij)
                              if((flowco2s(ij).lt.tem(1)).and.
     &                             (flowco2s(ij).ne.0d0))then
                                 sk(ij) = skco2(ij)*(1-flowco2s(ij))
                                 skco2(ij) = skco2(ij)*flowco2s(ij)
                              else
                                 sk(ij) = skco2(ij)*(1-tem(1))
                                 skco2(ij) = skco2(ij)*tem(1)
                              endif
                           else
                              call co2_properties(10,1,flowco2s(ij),
     &                             -eskco2(ij),csalt(ij),1,xc_prime,tem,
     &							 ij)
                              skco2(ij) = sk(ij)*tem(1)
                              sk(ij) = sk(ij)*(1-tem(1))
							endif
                           else
                              call co2_properties(4,ices(ij),phico2(ij),
     &                             -eskco2(ij),dum1,duma,dumb,dumc,duma)
                           endif
                        else
                           call co2_properties(4,ices(ij),phico2(ij),
     &                          -eskco2(ij),dum1,duma,dumb,dumc,duma)
                        endif
                        eflowco2(ij)=dumc(4)
                     endif
                  elseif((kaco2(ij).eq.-4).or.(kaco2(ij).eq.-5)) then
                     call co2_properties(4,ices(ij),phico2(ij),
     &                    tco2(ij),dum1,duma,dumb,dumc,duma)
                     eflowco2(ij)=dumc(4)
                  else						
                     eflowco2(ij)=eskco2(ij)
                  endif
                  if((kaco2(ij).ne.2).and.(kaco2(ij).ne.3)) then
                     if(esk(ij).lt.0.0) then
                        if(ps(ij).eq.0.d0) then
                           eflow(ij) = cpr(ij)*(-eskco2(ij))
                        else
                           call h2o_properties(9,1,phi(ij),-esk(ij),
     &                          csalt(ij),dum3,ensrc,dum4,dum5,dum6,
     &                          dum7)
                           eflow(ij)=ensrc
                        endif
                     else
                        eflow(ij)=esk(ij)
                     endif
                  endif
               enddo
c           endif                
         elseif(iflg.eq.4) then
c     call equation generation and load jacobian array
            
            call gensco2h2o      
            
         elseif(iflg.eq.5) then
            if(ntty.eq.2) then
c     
c     output for co2 information
c     
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
               if(m.gt.0) then

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
c     Water                  
                  if (iout .ne. 0) write(iout,800)
                  if(iatty.ne.0) write(iatty,800)
                  do il=1,ilev
                     if(il.ne.1) then
                        if (iout .ne. 0) write(iout,600) il
                        if (iatty.gt.0) write(iatty,600) il
                     endif
                     do i=1,mlev
                        md=  nskw(i+(il-1)*mlev)
                        if (ico2dis(md) .eq. 1) then
                           if (ices(md) .eq. 2) then
                              iphase = 3
                           else
                              iphase = 2
                           end if
                        else
                           iphase = 1
                        end if
                           
                        if (iout .ne. 0) write(iout, 801)  md , pho(md),
     &                       tco2(md),
     &                       fw(md), yc(md), iphase, sk(md) , qh(md) 
                        if ( iatty .gt. 0 )  write(iatty, 801) md ,
     &                       pho(md), tco2(md), fw(md), yc(md), iphase, 
     &                       sk(md), qh(md)  
                     enddo
                  enddo
c     CO2
                  if (iout .ne. 0) then
                     write(iout,802)
                     write(iout,803)
                  end if
                  if(iatty.ne.0) then
                     write(iatty,802)
                     write(iatty,803)
                  end if
                  do il=1,ilev
                     if(il.ne.1) then
                        if (iout .ne. 0) write(iout,600) il
                        if (iatty.gt.0) write(iatty,600) il
                     endif
                     do i=1,mlev
                        md=  nskw(i+(il-1)*mlev)
                        skmd = skco2(md)
                        qhmd = qhco2(md)
                       if (iout .ne. 0) 
     &                       write(iout,804) md, phico2(md),  fl(md),
     &                       fg(md), ices(md), skmd, qhmd       
                        if(iatty.ne.0)
     &                       write(iatty,804) md, phico2(md), fl(md), 
     &                       fg(md), ices(md), skmd, qhmd       
                     enddo
                  enddo
               endif
 600           format(2x,'Matrix Level = ',i1)
 800           format(/,20x,'Nodal Information (Water)', /,
     &              28x, 'Vol Frac', 2x, 'Mass Frac', 2x, 'Phase', 
     &              2x, 'W source/sink',2x,'E source/sink',/,
     &              3x,'Node',1x,'P(MPa)',4x,'Tco2(C)',3x,'Water',
     &              5x,'Diss CO2',3x,'State',5x,'(kg/s)',9x,'(MJ/s)')
 801           format(i7, 1x, g10.4, 1x, g9.3, 1x, g9.3, 3x, g9.3, 2x,
     &              i1, 5x, g11.3, 4x, g11.3)
 802           format(/, 20x, 'Nodal Information (CO2)', /, 21x,
     &              'Volume  Fraction', 6x, 'CO2', 4x, 
     &              'CO2 source/sink', 1x, 'E source/sink' )
 803           format(3x, 'Node', 1x, 'Pco2(MPa)', 1x, 'Liquid CO2',
     &              1x, 'Gaseous CO2', 2x, 'phase', 7x, '(kg/s)', 9x,
     &              '(MJ/s)')
 804           format(i7, 1x, g9.4, 1x, g9.3, 2x, g9.3, 6x, i1, 7x,
     &              g11.3, 4x, g11.3)
c     
c**** printout global mass and energy flows ****
c**** printout global mass and energy balances ****
c     
               if (iout .ne. 0) then
                  write(iout,700)
                  write(iout,709) amco20
                  write(iout,729) amco2f0
                  write(iout,719) amco2w0
                  write(iout,701) amco2
                  write(iout,722) amco2f
                  write(iout,712) amco2w
                  write(iout,711) amco21
                  write(iout,732) amco2f1
                  write(iout,713) amco2w1
                  write(iout,707) (skmd3)
                  write(iout,717) abs(skmd2)
                  write(iout,727) (skmd31)
                  write(iout,737) abs(skmd21)
                  write(iout,708) skmd1
                  write(iout,702)
                  write(iout,703) 
     &                 qco2,abs(qco2ts-qco2ts_in)/dtotdm,
     &                 abs(qco2ts_in)/dtotdm,
     &                 qeco2,abs(qeco2ts-qeco2ts_in)/dtotdm,
     &                 abs(qeco2ts_in)/dtotdm
                  write(iout,704)
               end if
               if (iatty.gt.0) then
                  write(iatty,700)
                  write(iatty,709) amco20
                  write(iatty,729) amco2f0
                  write(iatty,719) amco2w0
                  write(iatty,701) amco2
                  write(iatty,722) amco2f
                  write(iatty,712) amco2w
                  write(iatty,711) amco21
                  write(iatty,732) amco2f1
                  write(iatty,713) amco2w1
                  write(iatty,707) (skmd3)
                  write(iatty,717) abs(skmd2)
                  write(iatty,727) (skmd31)
                  write(iatty,737) abs(skmd21)
                  write(iatty,708) skmd1
                  write(iatty,702)
                  write(iatty,703) 
     &                 qco2,abs(qco2ts-qco2ts_in)/dtotdm,
     &                 abs(qco2ts_in)/dtotdm,
     &                 qeco2,abs(qeco2ts-qeco2ts_in)/dtotdm,
     &                 abs(qeco2ts_in)/dtotdm
                  write(iatty,704)
               end if
               if(baleco2.ne.-999999.) then
                  if (iout.ne.0) write(iout,705) balco2,baleco2
                  if (iatty.gt.0) write(iatty,705) balco2,baleco2
               else 
                  if (iout.ne.0) write(iout,706) balco2   
                  if (iatty.gt.0) write(iatty,706) balco2   
               endif

               inquire(file='co2inout.his',opened=it_is_open)
               if(.not. it_is_open) then
                  myfile=open_file('co2inout.his','unknown')
               endif
               write(myfile,'(g19.8,1x,g19.8,1x,g19.8)') days/365.25, 
     &              abs(skmd21), skmd31
               if(days.ge.tims) close(myfile)

 700           format(/,20x,'Global Mass & Energy for CO2')

 709           format(/,1x,'Initial Total CO2 (kg) in the System: ',
     &              e14.6)

 729           format(/,1x,'Initial Free CO2 (kg) in the System: ',
     &              e14.6)

 719           format(/,1x,'Initial Dissolved CO2 (kg) in the System: ',
     &              e14.6)

 701           format(1x,'Total CO2 (kg) at this time:      ',
     &              e14.6,' kg')
               
 722           format(1x,'Total Free CO2 (kg) at this time:      ',
     &              e14.6,' kg')
               
 712           format(1x,'Dissolved CO2 (kg) in system at this time: ',
     &              e14.6)
               
 711           format(1x,'Total CO2 (kg) in system at last time step:',
     &              e14.6)
               
 732           format(1x,'Free CO2 (kg) in system at last time step:',
     &              e14.6)
               
 713           format(1x,'Dissolved CO2(kg) in system at last time ',
     &              'step: ', e14.6)
               
 707           format(1x,'Total co2 produced (kg) at This Time Step ',
     &              e14.6)
               
 717           format(1x,'Total co2 injected (kg) at This Time Step '
     &              ,e14.6)
               
 727           format(1x,'Total co2 produced (kg) '
     &              ,e14.6)
               
 737           format(1x,'Total co2 injected (kg) '
     &              ,e14.6)

 708           format(/,1x,'Net kg co2 discharge (total out-total in): '
     &              ,e14.6)
               
 702           format(/,20x,'Global Mass & Energy flows for CO2')
 703           format(1x,
     &              'Mass flux:       Total      Time Step(in)',
     &              '       Time Step(out) ',/,
     &              9x,e14.6,5x,e14.6,7x,e14.6,' kg/s',/,1x,
     &              'Energy  flux:    Total      Time Step(in)',
     &              '       Time Step(out) ',/,
     &              9x,e14.6,5x,e14.6,7x,e14.6,' MJ/s')      
 704           format(/,20x,'Global Mass & Energy balances for CO2')
 705           format(1x,'Mass balance:',1x,e14.6,
     &              5x,'Energy balance:',1x,e14.6)
 706           format(1x,'Mass balance:',1x,e14.6,
     &              5x,'Energy balance: N/A')
            endif
         elseif(iflg.eq.6) then
c     initialize variables
            
            do i = 1, neq
			if((kaco2(i).eq.-3).or.(kaco2(i).eq.-5)) then
				tmp1 = 0.d0
				i1 = nelm(i)+1
				i2 = nelm(i+1)
				do j = i1, i2
					kb = nelm(j)
					if((ka(kb).eq.0).or.(kaco2(kb).eq.0))then
						ii = istrw(j-neq-1)
						if(ii.ne.0) then
						tmp1 = max(tmp1,abs(sx(ii,isox)+sx(ii,isoy)+
     &                       sx(ii,isoz)))
						endif
                     endif
                  enddo
				wellim(i) = tmp1*(pnx(i)+pny(i)+pnz(i))/3.d0
                  wellco2(i) = 10.d0*wellim(i)
               endif
            enddo
            if (.not. co2_read) then
               do i=1,n
                  if(ices(i).eq.2) then
                     call co2_properties(2,ices(i),pl,dumb,dum1,duma,
     &                    tliquid,dumc,duma)
                     tco2(i)=tliquid
                  else
                     tco2(i)= toco2(i)
                  endif
                  phoco2(i) = phico2(i)
                  fog(i)= fg(i)
                  if(iprtype.eq.1) then
                     xc(i)=0.d0
                     yw(i)=1.d0
                     xw(i)=1.d0
                     yc(i)=0.d0
                  elseif(iprtype.eq.2) then
                     xc(i)=1.d0
                     yw(i)=0.d0
                     xw(i)=0.d0
                     yc(i)=1.d0
                  elseif(abs(iprtype).eq.3) then
                     xc(i)=1.d0
                     yw(i)=1.d0
                     xw(i)=0.d0
                     yc(i)=0.d0
                  else
                     if(iwatdis.eq.1) then
c     calculate phase-change pressure and dp/dt
                        call h2o_properties(5,2,tco2(i),dum1,dum2,dum3,
     &                       psatd,dum4,dpsatt,dum5,dum6)
                        xwp = psatd/phico2(i)
                        xcp = 1.d0-xwp		  
                     else
                        xwp = 0.d0
                        xcp = 1.d0-xwp		  
                     endif
                     xcpb=(xcp/mwc)+(xwp/mww)
                     xc_prime = xcp/(mwc*xcpb)
                     xc(i) = xc_prime
                  endif				
                  xw(i) = 1.d0 - xc(i)
                  yw(i) = 1.d0 - yc(i)            
               enddo
               do i=1,n
                  if(ieos(i).eq.-2) then
                     call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                    tsolid,dtps,dum4,dum5,dum6)
                     t(i)=tsolid
                  else if(ieos(i).eq.2) then
                     call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                    tliquid,dtps,dum4,dum5,dum6)
                     t(i)=tliquid
                  else
                     t(i)= to(i)
                  endif
                  phi(i) = pho(i)
                  so(i)= s(i)
               enddo
            else
               do i=1,n
                  if(ices(i).eq.2) then
                     call co2_properties(2,ices(i),pl,dumb,dum1,duma,
     &                    tliquid,dumc,duma)
                     tco2(i)=tliquid
                  else
                     tco2(i)= toco2(i)
                  endif
                  phoco2(i) = phico2(i)
                  fog(i)= fg(i)
                  if(iprtype.eq.1) then
                     xc(i)=0.d0
                     xw(i)=1.d0
                  elseif(iprtype.eq.2) then
                     xc(i)=1.d0
                     xw(i)=0.d0
                  elseif(abs(iprtype).eq.3) then
                     xc(i)=1.d0
                     xw(i)=0.d0
                  else
                     if(iwatdis.eq.1) then
c     calculate phase-change pressure and dp/dt
                        call h2o_properties(5,2,tco2(i),dum1,dum2,dum3,
     &                       psatd,dum4,dpsatt,dum5,dum6)
                        xwp = psatd/phico2(i)
                        xcp = 1.d0-xwp		  
                     else
                        xwp = 0.d0
                        xcp = 1.d0-xwp		  
                     endif
                     xcpb=(xcp/mwc)+(xwp/mww)
                     xc_prime = xcp/(mwc*xcpb)
                     xc(i) = xc_prime
                  endif				
                  xw(i) = 1.d0 - xc(i)
                  yw(i) = 1.d0 - yc(i)
               end do               
            end if 
         elseif(iflg.eq.7) then
            
c     calulate initial mass and energy of co2
c     calulate initial mass and energy of co2 hydrate 
            
            amco20 = 0.0
            amco2w0 = 0.0
            aeco20 = 0.0
            qco2 = 0.0
            qeco2 = 0.0
            qco2_in = 0.0
            qeco2_in = 0.0
            skmd1 = 0.0
            do i=1,n
               denco2h(i) = denco2i(i)*dtot
               deneco2h(i) = deneco2i(i)*dtot
               denco2i(i) = 0.0 
               deneco2i(i) = 0.0 
               amco20 = amco20 + denco2h(i)*volume(i)
               amco2w0 = amco2w0 + denh(i)*volume(i)*yc(i)
               aeco20 = aeco20 + deneco2h(i)*volume(i)
            enddo
c     Total free CO2 = Total CO2 - Dissolved CO2
            amco2f0=amco20-amco2w0
c     calulate initial mass and energy of water   
            
            amh2o0 = 0.0
            aeh2o0 = 0.0
            qh2o = 0.0
            qeh2o = 0.0
            qh2o_in = 0.0
            qeh2o_in = 0.0
            do i=1,n
               amh2o0 = amh2o0 + denh(i)*volume(i)
               aeh2o0 = aeh2o0 + deneh(i)*volume(i)
            enddo
c     initialize space for co2 arrays
            
            nmat(37) = nmat(36)+nelm(neq+1)-(neq+1)
            nrhs(7)  = nrhs(6) +neq
            
         elseif(iflg.eq.8) then
            
c     calculate current mass and energy of co2
c     calculate current fluxes in and out of model  
c     note : for the current timestep iflg=8 
c     must be called after iflg=9
c     skmd1 = net flux of CO2 for this time step
c     skmd10 = net flux of CO2 for previous time step
c     skmd2 = Total CO2 injected in this time step
c     skmd3 = Total CO2 produced in this time step
            skmd10 = skmd1
            amco21 = amco2
            amco2w1 = amco2w
            amco2f1 = amco2f
            amco2 = 0.0
            amco2w = 0.0
            skmd2 = 0.0
            skmd3 = 0.0
            aeco2 = 0.0
            qco2ts = 0.0
            qeco2ts = 0.0
            qco2ts_in = 0.0
            qeco2ts_in = 0.0
            do i=1,n
               amco2 = amco2 + denco2h(i)*volume(i)
               amco2w = amco2w + denh(i)*volume(i)*yc(i)
               skmd1 = skmd1 + skco2(i)*dtotdm
               if (skco2(i).le.0.d0) then
                  skmd2 = skmd2 + skco2(i)*dtotdm
               else
                  skmd3 = skmd3 + skco2(i)*dtotdm
               endif
               if(sk(i).gt.0.d0) then
                  skmd3 = skmd3 + sk(i)*dtotdm*yc(i)
               endif
               dums1 = skco2(i)*dtotdm
               qco2 = qco2 + dums1
               aeco2 = aeco2 + deneco2h(i)*volume(i)
               qco2ts = qco2ts + dums1
               if(skco2(i).lt.0.0) then
                  qco2_in = qco2_in + dums1           
                  qco2ts_in = qco2ts_in + dums1                
               endif
               dums1 = qhco2(i)*dtotdm
               qeco2 = qeco2 + dums1                
               qeco2ts = qeco2ts + dums1                
               if(qhco2(i).lt.0.0) then
                  qeco2_in = qeco2_in + dums1                
                  qeco2ts_in = qeco2ts_in + dums1                
               endif
            enddo
c     Total free CO2 = Total CO2 - Dissolved CO2
            amco2f=amco2-amco2w
            skmd21=skmd2+skmd21
            skmd31=skmd3+skmd31
            
c     calculate current mass and energy of water   
c     calculate current fluxes in and out of model  
c     note : for the current timestep iflg=8 
c     must be called after iflg=9
            
            amh2o = 0.0
            aeh2o = 0.0
            qh2ots = 0.0
            qeh2ots = 0.0
            qh2ots_in = 0.0
            qeh2ots_in = 0.0
            do i=1,n
               amh2o = amh2o + denh(i)*volume(i)
               aeh2o = aeh2o + deneh(i)*volume(i)
               dums1 = sk(i)*dtotdm
               qh2o = qh2o + dums1           
               qh2ots = qh2ots + dums1           
               if(sk(i).lt.0.0) then
                  qh2o_in = qh2o_in + dums1           
                  qh2ots_in = qh2ots_in + dums1                
               endif
               dums1 = qh(i)*dtotdm
               qeh2o = qeh2o + dums1                
               qeh2ots = qeh2ots + dums1                
               if(qh(i).lt.0.0) then
                  qeh2o_in = qeh2o_in + dums1                
                  qeh2ots_in = qeh2ots_in + dums1                
               endif
            enddo
            
            
c     calculate mass and energy balance for co2
            
            amaxflx = max(abs(qco2_in),abs(qco2-qco2_in))
            if(amaxflx.gt.zero_t) then
               balco2 = (amco2 - amco20 + qco2) / amaxflx
            else if(amco20.gt.0.0) then
               balco2 = (amco2 - amco20 + qco2) / amco20
            else 
c     set co2 mass balance to n/a
               balco2 = -999999.0                                    
            endif

            amaxflx = max(abs(qeco2_in),abs(qeco2-qeco2_in))
            if(amaxflx.gt.zero_t) then
               baleco2 = (aeco2 - aeco20 + qeco2) / amaxflx
            else if(aeco20.gt.0.0) then
               baleco2 = (aeco2 - aeco20 + qeco2) / aeco20
            else 
c     set co2 energy balance to n/a
               baleco2 = -999999.0                                    
            endif
            
c     calculate mass and energy balance for water
            
            amaxflx = max(abs(qh2o_in),abs(qh2o-qh2o_in))
            if(amaxflx.gt.zero_t) then
               balh2o = (amh2o - amh2o0 + qh2o) / amaxflx
            else if(amh2o0.gt.0.0) then
               balh2o = (amh2o - amh2o0 + qh2o) / amh2o0
            else 
c     set h2o mass balance to n/a
               balh2o = -999999.0                                    
            endif
            difm = balh2o

            amaxflx = max(abs(qeh2o_in),abs(qeh2o-qeh2o_in))
            if(amaxflx.gt.zero_t) then
               baleh2o = (aeh2o - aeh2o0 + qeh2o) / amaxflx
            else if(aeh2o0.gt.0.0) then
               baleh2o = (aeh2o - aeh2o0 + qeh2o) / aeh2o0
            else 
c     set h2o energy balance to n/a
               baleh2o = -999999.0                                    
            endif
            
c     correct mass and energy for water in combining equations
            
            if(idof_co2.eq.3) then
               
               amaxflx = max(abs(qeh2o_in+qeco2_in),
     &              abs((qeh2o+qeco2)-(qeh2o_in+qeco2_in)))
               if(amaxflx.gt.zero_t) then
                  baleh2o = ((aeh2o+aeco2) - 
     &                 (aeh2o0+aeco20) + (qeh2o+qeco2)) / amaxflx
               else if(aeh2o0+aeco20.gt.0.0) then
                  baleh2o = ((aeh2o+aeco2) - (aeh2o0+aeco20)
     &                 + (qeh2o+qeco2)) / (aeh2o0+aeco20)
               else 
c     set h2o energy balance to n/a
                  baleh2o = -999999.0
               endif
               
c     correct h2o energy balance
               
               dife = baleh2o
               
            endif
            if(idof_co2.eq.5) then
               
c     set co2 energy balance to n/a
c     baleh2o represents energy balance of water, co2, and hydrate
               
               
               amaxflx = max(abs(qeh2o_in+qeco2_in+qeco2hyd_in),
     &              abs((qeh2o+qeco2+qeco2hyd)-
     &              (qeh2o_in+qeco2_in+qeco2hyd_in)))
               amaxener = (aeh2o+aeco2+aeco2hyd)-(aeh2o0+aeco20+
     &              aeco2hyd0)
c     if(amaxflx.gt.zero_t.and.amaxener.gt.zero_e) then
               if(amaxener*amaxflx.gt.zero_e) then
                  baleh2o = ((aeh2o+aeco2+aeco2hyd) - 
     &                 (aeh2o0+aeco20+aeco2hyd0)+(qeh2o+qeco2+qeco2hyd))
     &                 /amaxflx
               else if(aeh2o0+aeco20+aeco2hyd0.gt.0.0) then
                  baleh2o = ((aeh2o+aeco2+aeco2hyd) - 
     &                 (aeh2o0+aeco20+aeco2hyd0) + (qeh2o+qeco2+
     &                 qeco2hyd)) / (aeh2o0+aeco20+aeco2hyd0)
               else 
c     set h2o energy balance to n/a
                  baleh2o = -999999.0
               endif
               
c     correct h2o energy balance
               
               dife = baleh2o
            endif
         elseif(iflg.eq.9) then
            
c     advance variable in transient simulation
c     
            do i = 1,n
               phoco2(i) = phico2(i)
               toco2(i) = tco2(i)
               fow(i) = fw(i)
CHari
               so(i) = fw(i)
               fog(i) = fg(i)
               fol(i) = fl(i)
               denco2h(i) = denco2h(i) + denco2i(i)*dtot
               deneco2h(i) =  deneco2h(i) + deneco2i(i)*dtot
               denco2i(i) = 0.0
               deneco2i(i) = 0.0
               xoc(i) = xc(i)
               yoc(i) = yc(i)
               xow(i) = xw(i)
               yow(i) = yw(i) 
               ieoso(i) = ieos(i)
               iceso(i) = ices(i)
               ico2diso(i) = ico2dis(i)
            enddo
            
         elseif(iflg.eq.11) then
            
c     reset variables
            
            do i = 1,n
               phico2(i) = phoco2(i)
               phi(i) = pho(i)
               t(i) = to(i)
               tco2(i) = toco2(i)
               fw(i) = fow(i)
               fg(i) = fog(i)
               fl(i) = fol(i)
               ieos(i)= ieoso(i)
               ices(i)= iceso(i)
               denco2i(i) = 0.0
               deneco2i(i) = 0.0 
               xc(i) = xoc(i)
               yc(i) = yoc(i)
               xw(i) = xow(i)
               yw(i) = yow(i) 
               ieos(i) = ieoso(i)
               ices(i) = iceso(i) 
               ico2dis(i) = ico2diso(i)
            enddo         
            
         elseif(iflg.eq.-20) then
c     write descriptors
c     if (ico2.gt.0) write(isave,'(a4)') 'ngas'
c     if (ico2.lt.0) write(isave,'(a4)') 'air '
c     if (ico2.eq.0) write(isave,'(a4)') 'h20 '
c     if (iccen.ne.0) write(isave,'(a4)') 'trac'
c     if (istrs.eq.0) write(isave,'(a4)') 'nstr'
c     if (istrs.ne.0) write(isave,'(a4)') 'strs'
c     if (idpdp.eq.0) write(isave,'(a4)') 'ndpd'
c     if (idpdp.ne.0) write(isave,'(a4)') 'dpdp'
c     if (idualp.eq.0) write(isave,'(a4)') 'ndua'
c     if (idualp.ne.0) write(isave,'(a4)') 'dual'
            
c     write out restart file
            write(isave, 6002) 'carb'
            write(isave ,6002)  (max(to(mi),tolw),   mi=1,n )
            write(isave ,6002)  (max(so(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(pho(mi),tolw),  mi=1,n )
            
            write(isave ,6002)  (max(toco2(mi),tolw),   mi=1,n )
            write(isave ,6002)  (max(phoco2(mi),tolw),  mi=1,n )
            
            write(isave ,6002)  (max(fow(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(fog(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(fol(mi),tolw),    mi=1,n )

            write(isave ,6003)  (ieoso(mi),    mi=1,n )
            write(isave ,6003)  (iceso(mi),    mi=1,n )
c     if (iccen.ne.0) call  concen  ( 4,0,dummyreal )
 6002       format(4g25.16)
 6003       format(25i4)
         elseif(iflg.eq.20) then
            
c     read(iread,'(a4)') wdd1(1:4)
c     read(iread,'(a4)') wdd1(5:8)
c     read(iread,'(a4)') wdd1(9:12)
c     read(iread,'(a4)') wdd1(13:16)
c     read(iread,'(a4)') wdd1(17:20)
c     read restart file (will not work for double porosity)
            read(iread,'(a4)') wdd1(1:4)
            if(wdd1(1:3).eq.'h20') then
               read(iread,'(a4)') wdd1(5:8)
               read(iread,'(a4)') wdd1(9:12)
               read(iread,'(a4)') wdd1(13:16)
               read(iread,'(a4)') wdd1(17:20)
               if(iriver.eq.1) then
                  n = n-npoint_riv
                  read(iread ,*)  (to(mi), mi=1,n )
                  read(iread ,*)  (so(mi), mi=1,n )
                  read(iread ,*)  (pho(mi), mi=1,n )
                  do i = n+1, n+npoint_riv
                     j = iriver_con_node(i-n)
                     to(i) = to(iriver_con_node(i-n))
                     so(i) = so(iriver_con_node(i-n))
                     pho(i) = pho(iriver_con_node(i-n))
                  enddo
                  n = n+npoint_riv
               else
                  read(iread ,*)  (to(mi), mi=1,n )
                  read(iread ,*)  (so(mi), mi=1,n )
                  read(iread ,*)  (pho(mi), mi=1,n )
               endif
               do i = 1, n
                  if(inico2flg(i).eq.1) then
                     fow(i) = fw_tmp(i)
                     fol(i) = fl_tmp(i)
                     fog(i) = fg_tmp(i)
                  else
                     fow(i) = 1.d0
                     fog(i) = 0.d0
                     fol(i) = 0.d0
                  endif
               enddo

               t = to
               s = so
               phi = pho
               toco2 = to
               phoco2 = pho
               tco2 = toco2
               phico2 = phoco2
c     fow = 1.d0
c     fog = 0.d0
c     fol = 0.d0
               fw = fow
               fg = fog
               fl = fol
               do ii=1, n
                  pl = phico2(ii)
                  tl = tco2(ii)
                  call co2_properties(1,duma,pl,tl,dum1,
     &                 ices(ii),dumb,dumc,duma)
               enddo
               iceso = ices
               ieos = 1
               ieoso = ieos
            else
c     if(iriver.eq.1) then
               if(wdd1(1:3).eq.'car') then
                  if(iriver.eq.1) then
                     n = n-npoint_riv
                     read(iread ,*)  (to(mi), mi=1,n )
                     read(iread ,*)  (so(mi), mi=1,n )
                     read(iread ,*)  (pho(mi), mi=1,n )
                     
                     read(iread ,*)  (toco2(mi),  mi=1,n )
                     read(iread ,*)  (phoco2(mi), mi=1,n )
                     
                     read(iread ,*)  (fow(mi),  mi=1,n )
                     read(iread ,*)  (fog(mi),  mi=1,n )
                     read(iread ,*)  (fol(mi),  mi=1,n )
                     
                     read(iread ,*)  (ieoso(mi),  mi=1,n )
                     read(iread ,*)  (iceso(mi),  mi=1,n )

                     do i = n+1, n+npoint_riv
                        j = iriver_con_node(i-n)
                        to(i) = to(iriver_con_node(i-n))
                        so(i) = so(iriver_con_node(i-n))
                        pho(i) = pho(iriver_con_node(i-n))
                        toco2(i) = toco2(iriver_con_node(i-n))
                        phoco2(i) = phoco2(iriver_con_node(i-n))
                        fow(i) = fow(iriver_con_node(i-n))
                        fog(i) = fog(iriver_con_node(i-n))
                        fol(i) = fol(iriver_con_node(i-n))
                        ieoso(i) = ieoso(iriver_con_node(i-n))
                        iceso(i) = iceso(iriver_con_node(i-n))
                     enddo
                     n = n+npoint_riv
                  else
                     read(iread ,*)  (to(mi), mi=1,n )
                     read(iread ,*)  (so(mi), mi=1,n )
                     read(iread ,*)  (pho(mi), mi=1,n )
                     
                     read(iread ,*)  (toco2(mi),  mi=1,n )
                     read(iread ,*)  (phoco2(mi), mi=1,n )
                     
                     read(iread ,*)  (fow(mi),  mi=1,n )
                     read(iread ,*)  (fog(mi),  mi=1,n )
                     read(iread ,*)  (fol(mi),  mi=1,n )
                     
                     read(iread ,*)  (ieoso(mi),  mi=1,n )
                     read(iread ,*)  (iceso(mi),  mi=1,n )

                  endif
                  t = to
                  s = so
                  phi = pho
                  tco2 = toco2
                  phico2 = phoco2
                  fw = fow
                  ieos = ieoso
                  ices = iceso
                  fg = fog
                  fl = fol
               endif
            endif
         elseif(iflg.eq.22) then
            do i = 1, n
               ieosd = ieos(i)
               icesd = ices(i)
               
               if(ieosd.eq.2) then
                  pl = phi(i)
c     two phase (gas/liquid) conditions
                  
c     calculate phase-change temperature and dt/dp
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
c     calculate phase-change pressure and dp/dt
                  call h2o_properties(5,2,tl,dum1,dum2,dum3,
     &                 psatd,dum4,dpsatt,dum5,dum6)
c     
               else if(ieosd.eq.-2) then
                  pl = phi(i)
c     two phase (liquid/solid) conditions
                  
c     calculate temperature and dt/dp
                  call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                 tl,dtps,dum4,dum5,dum6)
c     calculate pressure  and dp/dt
                  call h2o_properties(5,1,tl,dum1,dum2,dum3,
     &                 psatd,dum4,dpsatt,dum5,dum6)
               else
                  tl=t(i)
               endif
               tl1 = tl
               dtps1 = dtps
               if(icesd.eq.2) then
                  pl = phico2(i)
c     two phase (gas/liquid) conditions
                  
c     calculate phase-change temperature and dt/dp
                  call co2_properties(2,icesd,pl,dumb,dum1,duma,tl,dumc,
     &			duma)
                  dtps=dumc(1)
                  call co2_properties(3,icesd,dumb,tl,dum1,duma,pl,dumc,
     &			duma)
                  psatd=dumc(1)
               else
                  tl=tco2(i)
               endif
               tl2 = tl
               dtps2 = dtps
c     RJP 04/10/05 Enforce thermal equilibrium
               if ((ieosd.eq.2).and.(icesd.ne.2)) then
                  tco2(i) = tl1
                  t(i) = tl1
                  dtpsc(i) = dtps1
               endif
               if ((icesd.eq.2).and.(ieosd.ne.2)) then
                  tco2(i) = tl2
                  t(i) = tl2
                  dtpsc(i) = dtps2
               endif
            enddo

         elseif(iflg.eq.21) then
            
c     write surfer or tecplot contour files
            
            write(iscon,790)
 790        format('         ','X',12x,'Y',12x,'Z',6x,'   pres',12x,'t',
     &           8x,'fracw',8x,'fco2l',7x,'fco2g',7x,'skco2')
            j = ndummy
            do k=1,neq_primary
               if(izone_surf_nodes(k).eq.j) then
                  
                  write(iscon,855) 
     &                 cord(k,1),cord(k,2),cord(k,3),phi(k),
     &                 t(k),fw(k),fl(k),fg(k),skco2(k)
               endif
            enddo
 855        format(1x,10(1x,g12.5))
         endif
         
      endif
      
      return
      end

      subroutine icectr(iflg,ndummy)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To manage the (solid-liquid-gas) calculations for multi-component 
!D1 systems.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/icectr.f_a  $
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
      use commeth
      use comsplitts 

      implicit none

      integer iflg,ndummy,i,nr1,nr2,nr3,nr4,nr5,mi      
      integer i1,i2,i3,i4,i5,ilev,mlev,il,md,k,j
      integer ii,ij,icedc,iced,ieosd,ihydd,ihyddc,np_hyd
c ich_max limits the number of phase changes/timestep (gaz 7-17-06)
c ich_max is variable and in comai
c 2-3 seems pretty good for MH (still testing)
      real*8, allocatable :: aiped(:)
      real*8, allocatable :: sktmp(:)
      real*8, allocatable :: esktmp(:)

      real*8 tsolid,tliquid,aihyd1
      real*8 pl,tl,sl,ensrc,phase_frac,skmd,qhmd,skmd1,skmd10 
      real*8 dtps,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real*8 eostol
      real*8 eosmg
      real*8 eosml
      real*8 amaxflx, dums1,tolw, zero_e
      real*8 teinfl, flemax, amaxener         
      real*8 hyd_frac,w_frac,pdis1,tdis1,frach
      real*8 fh, power, dpowerh, px,py,pz,frach1,frach2,fracg       
c  tenma 04/12/2005 => AIST parallel layer (Sand & Mud)
      real*8 fracw_temp, frachyd_temp, para_mix
      real*8 temp_kz_s, temp_kz_m, temp_z, temp_poroz, temp_fh
	real*8 normal_hyd
c  tenma 04/12/2005 => AIST parallel layer (Sand & Mud)
      character*8 macro1
      parameter(eostol=0.0001d00)
      parameter(eosmg=1.0001d00)
      parameter(eosml=0.999d00)
      parameter(tolw=1.d-40, zero_e=1.e-05)
C      parameter(tolw=1.d-90, zero_e=1.e-04)
      save skmd1

      if(ice.ne.0) then  
c
c these values (ie large) seemto be best for MH
c
         ich_m1 = 100
         ich_m2 = 100

         if(iflg.eq.0) then

c     read in input when meth is present
            read (inpt  ,   *) ice, idof_meth
            if(idof_meth.eq.1) then
               idof=4
               if (iout .ne. 0) write(iout, 20) 
               if (iptty .gt. 0) write(iptty, 20)  
            else if(idof_meth.eq.2) then
               idof=4
               if (iout .ne. 0) write(iout, 21) 
               if (iptty .gt. 0) write(iptty, 21)  
            else if(idof_meth.eq.3) then
               idof=6
               if (iout .ne. 0) write(iout, 22) 
               if (iptty .gt. 0) write(iptty, 22)  
            else if(idof_meth.eq.0) then
               idof=6
               if (iout .ne. 0) write(iout, 23) 
               if (iptty .gt. 0) write(iptty, 23)  
            else if(idof_meth.eq.4) then
               idof=6
               if (iout .ne. 0) write(iout, 24) 
               if (iptty .gt. 0) write(iptty, 24)  
            else if(idof_meth.eq.5) then
               idof=6
               if (iout .ne. 0) write(iout, 25) 
               if (iptty .gt. 0) write(iptty, 25) 
            else if(idof_meth.eq.7) then
               idof=6
               if (iout .ne. 0) write(iout, 25) 
               if (iptty .gt. 0) write(iptty, 26)   
            else 
               if (iout .ne. 0) write(iout, *) 
     &              'Stopping : idof_meth not valid' 
               if (iptty .gt. 0) write(iptty, *)
     &              'Stopping : idof_meth not valid'  
            endif
 20         format(3x,'**** equilibrium condition: ',
     &           'P and T (2-phase M) ****') 
 21         format(3x,'**** equilibrium condition: ',
     &           'P and T (1-phase M) ****') 
 22         format(3x,'**** equilibrium condition: T (1-phase M) ****') 
 23         format(3x,'**** equilibrium condition: ',
     &           'Components independent')
 24         format(3x,'**** equilibrium meth_frac: ',
     &           'Components independent')
 25         format(3x,'**** variables: frac_m,frac_w,phi_m,temperature')
 26         format(3x,'**** equilibrium cals: rate equations = 0.0')
c     
c     allocate memory for methane arrays
c     
            allocate(denmethi(n0))
            allocate(denmethh(n0))
            allocate(denemethi(n0))
            allocate(denemethh(n0))
            allocate(pcpmeth(n0))
            allocate(pflowmeth(n0))
            allocate(phimeth(n0))
            allocate(phometh(n0))
            allocate(smeth(n0))
            allocate(someth(n0))
            allocate(tmeth(n0))
            allocate(tometh(n0))
            allocate(skmeth(n0))

            allocate(rl_h2o(n0))
            allocate(rl_meth(n0))
c
c     save storage gaz 1-13-04
c     allocate(skmethw(n0))
            allocate(eskmeth(n0))
            allocate(eflowmeth(n0))
            allocate(qhmeth(n0))
            allocate(wellmeth(n0))
            allocate(rovfmeth(n0))
            allocate(envfmeth(n0))
            allocate(rolfmeth(n0))
            allocate(enlfmeth(n0))
            allocate(dpcpmeth2(n0))
            allocate(dmmf(n0))
            allocate(demf(n0))
            allocate(dilm(n0))
            allocate(divm(n0))
            allocate(deqm(n0))
            allocate(deqpm(n0))
            allocate(dmwf(n0))
            allocate(dewf(n0))
            allocate(dilw(n0))
            allocate(divw(n0))
            allocate(deqw(n0))
            allocate(dq3(n0))
            allocate(dq4(n0))
            allocate(deqpw(n0))
            allocate(fracw(n0))
            allocate(fracwo(n0))
            allocate(qhflxmeth(n0))
            allocate(dtpameth(n0))
            allocate(dtpaemeth(n0))
            allocate(kameth(n0))
            allocate(denhydi(n0))
            allocate(denhydh(n0))
            allocate(denehydi(n0))
            allocate(denehydh(n0))
            allocate(frachyd(n0))
            allocate(frachydo(n0))
            allocate(fracgas(n0))
            allocate(fracgaso(n0))
            allocate(enlfhyd(n0))
            allocate(phihyd(n0))
            allocate(permhyd(n0))
            allocate(permhyd0(n0))
            allocate(skhyd(n0))
            allocate(qhhyd(n0))
            allocate(skmhyd(n0))
            allocate(qhmhyd(n0))
            allocate(skwhyd(n0))
            allocate(qhwhyd(n0))
            allocate(dskhyd1(n0))
            allocate(dqhhyd1(n0))
            allocate(dskhyd2(n0))
            allocate(dqhhyd2(n0))
            allocate(dskhyd3(n0))
            allocate(dqhhyd3(n0))
            allocate(dskhyd4(n0))
            allocate(dqhhyd4(n0))
            allocate(dpermhyd4(n0))
            allocate(nhyd(n0))
            allocate(iceso(n0))
            allocate(ihydo(n0))
            allocate(ieoso(n0))
            allocate(dqm(n0))
            allocate(dqd(n0,9))
            allocate(dskmethw1(n0))
            allocate(dskmethw3(n0))
            allocate(dskmethw4(n0))
            allocate(dqmw(n0))
            allocate(dqmhyd4(n0))
            allocate(dqmt(n0))
            allocate(dqmh(n0))
            allocate(deqmh(n0))
            allocate(ihyd(n0))
            allocate(dpcp3(n0))
            allocate(dpcp4(n0))
c  tenma 04/12/2005 => AIST parallel layer (Sand & Mud)
            allocate(frac_mudw(n0))
            allocate(poro_sand (n0))
            allocate(poro_mud (n0))
            allocate(z_sand (n0))
            allocate(z_mud (n0))
            allocate(permhyd_z(n0))
            allocate(ka0_xy(n0))
            allocate(perm_sand (n0))
            allocate(perm_mud (n0))
c     allocate(frac_mudwo(n0))
c tenma 06/22/2005 => sakamoto model Sh1
            allocate(frachyd_1(n0))
            allocate(nhyd_2(n0))
            allocate(nhyd_3(n0))
            allocate(frachyd_max(n0))

            allocate(fracg2(n0))
            allocate(fracw2(n0))
c            
c     zero out arrays for methane        
c     
            ihyd = 0
            
            denmethi = 0.0d00
            denmethh = 0.0d00
            denemethi = 0.0d00
            denemethh = 0.0d00
            pcpmeth = 0.0d00
            pflowmeth = 0.0d00
            phimeth = 0.0d00
            phometh = 0.0d00
            smeth = 0.0d00
            someth = 0.0d00
            tmeth = 0.0d00
            tometh = 0.0d00
            skmeth = 0.0d00
            qhmeth = 0.0d00
            eskmeth = 0.0d00
            eflowmeth = 0.0d00
            rovfmeth = 0.0d00
            envfmeth = 0.0d00
            wellmeth = 0.0d00
            rolfmeth = 0.0d00
            enlfmeth = 0.0d00
            dpcpmeth2 = 0.0d00
            dmmf = 0.0d00
            demf = 0.0d00
            dilm = 0.0d00
            divm = 0.0d00
            deqm = 0.0d00
            deqpm = 0.0d00
            dmwf = 0.0d00
            dewf = 0.0d00
            dilw = 0.0d00
            divw = 0.0d00
            deqw = 0.0d00
            deqpw = 0.0d00
            fracw = 0.0d00
            fracwo = 0.0d00
            qhflxmeth = 0.0d00
            dtpameth = 0.0d00
            dtpaemeth = 0.0d00
            skmd1 = 0.0d00
            frachyd = 0.0d00
            frachydo = 0.0d00
            fracgas = 0.0d00
            fracgaso = 0.0d00
            denhydi = 0.0d00
            denhydh = 0.0d00
            denehydi = 0.0d00
            denehydh = 0.0d00
            enlfhyd = 0.0d00
            phihyd = 0.0d00
            permhyd = 0.0d00
            permhyd0 = 0.0d00
            skhyd = 0.0d00
            qhhyd = 0.0d00
            skmhyd = 0.0d00
            qhmhyd = 0.0d00
            skwhyd = 0.0d00
            qhwhyd = 0.0d00
            dskhyd1 = 0.0d00
            dqhhyd1 = 0.0d00
            dskhyd2 = 0.0d00
            dqhhyd2 = 0.0d00
            dskhyd3 = 0.0d00
            dqhhyd3 = 0.0d00
            dskhyd4 = 0.0d00
            dqhhyd4 = 0.0d00
            dpermhyd4 = 0.0d00
            dq3 = 0.d00
            dq4 = 0.d00
            dpcp3 = 0.00
            dpcp4 = 0.00
            kameth = 0
            nhyd = 0
            iceso = 0
            ieoso = 0
            ihydo = 0
c     initialize  da,oe,pe,le,me,ne,e_act,permh
            ihyd_grow_type = 2
            ihyd_diss_type = 1
            da = .375d-6
            oe = 1.d0
            pe = 1.d0
            le = 0.d0
            me = 0.d0
            ne = 0.d0 
            e_act1 = 9400.d0
            e_act2 = 0.d0
            afhyd = 1.24d05
            dar = 1.d05
            oer = -0.572d0
            per = .422d0
            ler = 0.d0
            mer = 1.d0
            ner = 2.441d0
            e_act1r = 2.158E04
            e_act2r = -87.58
            afhydr = 1.d10
            hdiss01 = 215.59d-3
            hdiss11 = -394.945d-6
            hdiss02 = 446.12d-3
            hdiss12 = -132.638d-6
            ptc1 = 0.1091d0
            ptc2 = -29.173d0
            ptc3 = 0.0
            ptc4 = 0.0
            imeth_pt = 1
            ic_hyd_tot = 0
            strd_meth = 0.95d0
            strd_meth0 = 1.d0
c  tenma 04/12/2005 => AIST parallel layer (Sand & Mud)
            frac_mudw = 0.0d00
            poro_sand = 0.0d00
            poro_mud = 0.0d00
            z_sand = 0.0d00
            z_mud = 0.0d00
            perm_sand = 0.0d00
            perm_mud = 0.0d00
c	  frac_mudwo = 0.0d00  
            frachyd_1 = 0.0d00
c        nhyd_2 = 0
            nhyd_2 = -1      

            rathyd = 7.46875d0
            ratgas = 1.d0
            ratw = 6.46875d0
            
            fracv=0.1339 
            fracl=0.8661

            macro = "ice "
c     
c     read in "sub macros" for clatharte
c     

            
 100        continue
            read (inpt, '(a80)') wdd1
            if (wdd1(1:7) .eq. 'methend') go to 200 
            if (wdd1(1:1) .eq. '#') go to 40 
            read (wdd1, '(a8)') macro1
            if (iout .ne. 0) write(iout, 50) macro1
            if (iptty .gt. 0) write(iptty, 50) macro1 
 50         format(3x, '**** clathrate sub macro : ', a8,' **** ') 
            if(macro1.eq.'methpres') then
c     
c     read in initial clathrate pressure
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
c     
c     read in initial clathrate saturation and phase state
c     set hydrate phase state here
c     
               call initdata2( inpt, ischk, n0, narrays, itype, 
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1=phimeth(1:n0),r8_2=sktmp(1:n0),i4_1=ices(1:n0))

               do i = 1, n0
                  if (sktmp(i) .ne. default(2)) then
                     if (abs (ices(i)) .eq. 1) then
                        tmeth(i) =  sktmp(i)
                        smeth(i) =  1.0
                     else if (abs (ices(i)) .eq. 2) then
                        smeth(i) =  sktmp(i)
                        pl=phimeth(i)
                        iced = ices(i)
                        call methane_properties(4,iced,pl,dum1,dum2,
     &                       dum3,tl,dtps,dum4,dum5,dum6)
                        tmeth(i)=tl
                     else if (abs (ices(i)) .eq. 3) then
                        tmeth(i) =  sktmp(i)
                        smeth(i) =  0.0
                     end if
                  end if
                  
               end do

               deallocate(sktmp)

            else if(macro1.eq.'methfrac') then

               igroup = 2
               narrays = 3
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               default(1) = 0.
               default(2) = 0.
               default(3) = 0.
c     
c     read in initial clathrate mass fraction
c     
               call initdata2( inpt, ischk, n0, narrays, itype, 
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = fracw(1:n0),r8_2 = frachyd(1:n0),
     &              r8_3 = qhflxmeth(1:n0))
c
c
            frachyd_1(1) = 0
            do i = 2,n0
              frachyd_1(i) = frachyd_1(1)
            enddo
c     
c     added by RJP 10/30/03
c     
c     
c     assume no initial gas present
c     
               ihydfl=0
c     
c     if there is gas turn flag to 1, if only water no hydrate
c     turn flag to 1 as well
c     
               do i = 1, n0
                  if (fracw(i).eq.1.d0) then
                     ihydfl=1
                  elseif ((1.d0-fracw(i)-frachyd(i)).gt.0.d0) then
                     ihydfl=1 
                  elseif (frachyd(i).eq.0.d0) then
                     ihydfl=1
                  endif
               enddo

c  temma add 2005/06/22
c  Sakamoto-model 
c
            else if(macro1.eq.'methfra2') then
	 
               igroup = 2
               narrays = 4
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               itype(4) = 8
               default(1) = 0.
               default(2) = 0.
               default(3) = 0.
               default(4) = 0.
c     
c     read in initial clathrate mass fraction
c 
               call initdata2( inpt, ischk, n0, narrays, itype,
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = fracw(1:n0),r8_2 = frachyd_1(1:n0),
     &              r8_3 = frachyd(1:n0),r8_4 = qhflxmeth(1:n0))
c
c     added by RJP 10/30/03
c
c
c     assume no initial gas present
c
               ihydfl=0
c
c     if there is gas turn flag to 1, if only water no hydrate
c     turn flag to 1 as well
c
               do i = 1, n0
                  if (fracw(i).eq.1.d0) then
                     ihydfl=1
                  elseif ((1.d0-fracw(i)-frachyd(i)).gt.0.d0) then
                     ihydfl=1 
                  elseif (frachyd(i).eq.0.d0) then
                     ihydfl=1
                  endif
               enddo
c
cc
c   tenma 04/12/2005 => AIST parallel layer (Sand & Mud)

            else if(macro1.eq.'methpara') then

               igroup = 2
               narrays = 4
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               itype(4) = 8
               default(1) = 0.
               default(2) = 0.
               default(3) = 0.
               default(4) = 0.
c     
c     read in initial clathrate mass fraction
c 
               call initdata2( inpt, ischk, n0, narrays, itype,
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = fracw(1:n0),r8_2 = frachyd(1:n0),
     &              r8_3 = qhflxmeth(1:n0),r8_4 = frac_mudw(1:n0))
c
               do i = 1, n0
                  fracw_temp = fracw(i) 
	          frachyd_temp = frachyd(i)
                  para_mix = z_sand(i)*poro_sand(i) + 
     &                 z_mud(i)*poro_mud(i)
                  frachyd(i)= z_sand(i)*poro_sand(i)*frachyd_temp/
     &                 para_mix
                  fracw(i) = ( z_sand(i)*poro_sand(i)*fracw_temp +
     &                 z_mud(i)*poro_mud(i)*frac_mudw(i)) /
     &                 para_mix
               enddo  
c
c     added by RJP 10/30/03
c
c
c     assume no initial gas present
c
               ihydfl=0
c
c     if there is gas turn flag to 1, if only water no hydrate
c     turn flag to 1 as well
c
               do i = 1, n0
                  if (fracw(i).eq.1.d0) then
                     ihydfl=1
                  elseif ((1.d0-fracw(i)-frachyd(i)).gt.0.d0) then
                     ihydfl=1 
                  elseif (frachyd(i).eq.0.d0) then
                     ihydfl=1
                  endif
               enddo

            else if(macro1.eq.'methperm') then
               
               igroup = 3
               narrays = 2
               itype(1) = 8
c               itype(2) = 4
               itype(2) = 8
               default(1) = 1.
               default(2) = 0
c     
c     read in initial clathrate mass fraction
c     
               call initdata2( inpt, ischk, n0, narrays, itype, 
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = permhyd0(1:n0),r8_2 = nhyd(1:n0))

            else if(macro1.eq.'methper2') then

               igroup = 3
               narrays = 3
               itype(1) = 8
c            itype(2) = 4
               itype(2) = 8
               itype(3) = 4
               default(1) = 1.
               default(2) = 0
               default(3) = 0
c            
c            nhyd_2(1) = -1
c            do i = 2,n0
c              nhyd_2(i) = nhyd_2(1)
c            enddo            
c     
c     read in initial clathrate mass fraction
c 
               call initdata2( inpt, ischk, n0, narrays, itype,
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = permhyd0(1:n0),r8_2 = nhyd(1:n0),
     &              i4_1 = nhyd_2(1:n0))

            else if(macro1.eq.'methper3') then	 

               igroup = 3
               narrays = 3
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               default(1) = 1.
               default(2) = 1.
               default(3) = 0
c     
c     read in initial clathrate mass fraction
c 
               call initdata2( inpt, ischk, n0, narrays, itype,
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = permhyd0(1:n0),r8_2 = frachyd_max(1:n0),
     &              r8_3 = nhyd_3(1:n0))

            else if(macro1.eq.'methpowe') then

               read (inpt, '(a80)') wdd1
               do i=1,80
                  if(wdd1(i:i).ne.' ') go to 899
               enddo
               go to 900
 899           continue
               read(wdd1,*) np_hyd
               backspace inpt     
               nhyd(1) = -np_hyd
               allocate(ahyd_coef(np_hyd+1))
               read(inpt,*) np_hyd,permhyd0(1),
     &              (ahyd_coef(i),i=1,np_hyd+1)
 900           continue
               do i = 2,n0
                  permhyd0(i) = permhyd0(1)
               enddo

c tenma 04/19/2005 Absolute permeability 
c  approximation formula K/Ko = 207.7891*Sh^6 - 522.2554*Sh^5 + 526.1379*Sh^4
c                              -272.9801*Sh^3 + 78.3724*Sh^2 - 12.4814*Sh + 0.9835
            else if(macro1.eq.'methappr') then
c modified gaz 080105
               read (inpt, '(a80)') wdd1
               do i=1,80
                  if(wdd1(i:i).ne.' ') go to 901
               enddo
               go to 902
 901           continue
               backspace inpt     
               nhyd(1) = -6
               allocate(ahyd_coef(7))
               read(inpt,*) permhyd0(1), (ahyd_coef(i),i=1,7)
 902           continue
               do i = 2,n0
                  permhyd0(i) = permhyd0(1)
               enddo
c tenma 04/19/2005 

            else if(macro1.eq.'methflow') then
               igroup = 4
               narrays = 3
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               default(1) = 1.d50
               default(2) = 0.0d00
               default(3) = 0.0d00
c     
c     read in initial clathrate flow and boundary data
c     
               allocate(aiped(n0),esktmp(n0),sktmp(n0))

               call initdata2( inpt, ischk, n0, narrays, itype, 
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1=sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0))
               do i=1,n0
                  if(sktmp(i).eq.default(1)) then
                     skmeth(i)=0.0
                  elseif(sktmp(i).eq.0.0) then
                     skmeth(i) = zero_e
                  else
                     eskmeth(i) = esktmp(i)
                     if (abs(aiped(i)) .lt. zero_e) then
                        skmeth(i) = sktmp(i)
                        kameth(i) = 1
                     else
                        if (aiped(i) .lt. 0.) then
                           kameth(i) = -2
                        else
                           kameth(i) = -1
                        end if
                        wellmeth(i) = abs(aiped(i)) * 1.0e+06
                        pflowmeth(i) = sktmp(i)
                     end if
                  endif
               enddo

               deallocate(sktmp,esktmp,aiped)

            else if(macro1.eq.'methdisv') then
c
c alternate input for hydrate rate equations
c
               call hyddiss(0,0)

            else if(macro1.eq.'methdiss') then

c     da - grain area term
c     oe - exponent of grain area term
c     pe - exponent of fugacity difference
c     le - exponent of hydrate saturation term
c     me - exponent of water saturation term
c     ne - exponent of gas saturation term
c     e_act - activation energy  (two terms)
c     afhydh - activity fequency term
               read (inpt,*) ihyd_diss_type
               if(ihyd_diss_type.eq.4) then
                  backspace inpt
                  read (inpt,*) ihyd_diss_type,da
               else if(ihyd_diss_type.eq.1) then
                  backspace inpt
                  read (inpt,*) ihyd_diss_type,da,e_act1,afhyd
               else
                  backspace inpt
                  read (inpt,*) ihyd_diss_type,da,oe,pe,le,me,ne,
     &                 e_act1,e_act2,afhyd
                  if(oe.lt.0.0) then
                     oe=abs(oe)
                     da = 1.0/da
                  endif
               endif

c     dar - grain area term (reformation)
c     oer - exponent of grain area term (reformation)
c     per - exponent of fugacity difference (reformation)
c     ler - exponent of hydrate saturation term (reformation)
c     mer - exponent of water saturation term (reformation)
c     ner - exponent of gas saturation term (reformation)
c     e_actr - activation energy (reformation or growth)
c     afhydr - activity fequency term (reformation or growth)

               read (inpt,*) ihyd_grow_type
               if(ihyd_grow_type.eq.1) then
                  backspace inpt
                  read (inpt,*) ihyd_grow_type,dar,e_act1r,afhydr

               else
                  backspace inpt
                  read (inpt,*) ihyd_grow_type,dar,oer,per,ler,mer,ner,
     &                 e_act1r,e_act2r,afhydr
                  if(oer.lt.0.0) then
                     oer=abs(oer)
                     dar = 1.0/dar
                  endif
               endif 

c     hdiss01 - constant term in diss enthalpy (low temp)
c     hdiss11 - linear term in diss enthalpy (low temp)
c     hdiss02 - constant term in diss enthalpy (high temp)
c     hdiss12 - linear term in diss enthalpy (high temp)

               read (inpt,*) hdiss01,hdiss11,hdiss02,hdiss12

            else if(macro1.eq.'methiter') then

               read (inpt,*) strd_meth0, strd_meth

            else if(macro1.eq.'methpt') then

               read (inpt,*) imeth_pt, ptc1, ptc2

c   tenma 04/12/2005 => AIST parallel layer (Sand & Mud)

            else if(macro1.eq.'methlaye') then	 

               igroup = 5
               narrays = 4
               itype(1) = 8
               itype(2) = 8
               itype(3) = 8
               itype(4) = 8
               default(1) = 0.
               default(2) = 0.
               default(3) = 0.
               default(4) = 0.
c     
c     read in initial clathrate mass fraction
c 
               call initdata2( inpt, ischk, n0, narrays, itype,
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = poro_sand(1:n0),r8_2 = z_sand(1:n0),
     &              r8_3 = poro_mud(1:n0),r8_4 = z_mud(1:n0))

            else if(macro1.eq.'permsamu') then

               igroup = 6
               narrays = 2
               itype(1) = 8
               itype(2) = 8
               default(1) = 0.
               default(2) = 0.
c     
c     read in initial clathrate mass fraction
c 
               call initdata2( inpt, ischk, n0, narrays, itype,
     &              default, macroread(8), macro, igroup, ireturn,
     &              r8_1 = perm_sand(1:n0),r8_2 = perm_mud(1:n0))

c   tenma 04/12/2005 => AIST parallel layer (Sand & Mud)

            else
               if (iout .ne. 0) write(iout,*) 
     &              'ERROR IN HYDRATE INPUT(STOPPING)'
               if (iptty .ne. 0) write(iptty,*)  
     &              'ERROR IN HYDRATE INPUT(STOPPING)'
               write(ierr,*) 'ERROR IN HYDRATE INPUT(STOPPING)'
            end if
 40         continue
            go to 100
 200        continue

            if(idof_meth.ne.7) then
               do i = 1,n0
c     set initial hydrate phase state
                  pl=phimeth(i)
                  tl = tmeth(i)
                  call hydrate_properties(6,1,pl,dum1,dum2,dum3,
     &                 tdis1,dum5,dum6,dum7,dum8)
                  if(tl.lt.tdis1.or.frachyd(i).le.0.0) then
                     if(fracw(i).ne.0.0) then
                        ihyd(i) = -1
                     else
                        ihyd(i) = -2
                     endif
                  elseif(tl.ge.tdis1.and.frachyd(i).gt.0.0) then 
                     ihyd(i) = 0
                  else
                     ihyd(i) = 1
                  endif
               enddo
            else  if(idof_meth.eq.7) then
c     equilibrium formulation (idof_meth = 7)	      
               ratgas = 0.0
               ratw = 0.0
               rathyd = 0.0
               do i = 1,n0
c     set initial hydrate phase state
c     note that we do not input frachyd
                  pl=phimeth(i)
                  tl = tmeth(i)
                  w_frac = fracw(i)
                  frach = frachyd(i)
                  fracg = 1.0d0-w_frac-frach
                  fracgas(i)= fracg
                  w_frac = fracw(i)
                  call hydrate_properties(6,1,pl,dum1,dum2,dum3,
     &                 tdis1,dum5,dum6,dum7,dum8)

                  
                  if(tl.lt.tdis1) then
c     set initial gas or water fraction
c     can have free gas or water but not both
c     calculate frachyd
                     if(w_frac.gt.0.0d0) then
                        fracgas(i) = frach*fracv
                        fracw(i) = 1.- fracgas(i)
                        ihyd(i) = -1
                     else
                        w_frac = frach*fracl
                        fracw(i) = w_frac
                        fracgas(i) = 1.- w_frac
                        ihyd(i) = -2
                     endif
                     
                  elseif(tl.ge.tdis1) then 
                     if(frachyd(i).ge.1.0) then
c     all hydrate
                        ihyd(i) = 1
                     else if(frachyd(i).lt.1.0.and.frachyd(i).gt.0.0) 
     &                       then
                        tmeth(i) = tdis1
                        ihyd(i) = 0
c     set initial gas or water fraction
c     can have free gas and/or water 
                        fracw(i) = w_frac + fracl*frach         
                        fracgas(i) = fracg + fracv*frach
                        
                     else
                        ihyd(i) = 1
                     endif
                  endif
               enddo
            endif
            
            numhyd = 0
            do i=1,n0
               phometh(i) = phimeth(i)
               tometh(i) = tmeth(i)
               someth(i) = smeth(i)
               iceso(i)  = ices(i)
               ihydo(i) = ihyd(i)
               fracwo(i) = fracw(i)
               frachydo(i) = frachyd(i)
               fracgaso(i) = fracgas(i)
               if(frachyd(i).gt.tolw) then
                  numhyd = numhyd + 1 
               endif
            enddo

            macroread(8) = .TRUE.

         elseif(iflg.eq.-2) then
c     
c     set tempertures of components equal for certain equilibrium conditions
c     
            
            if(idof_meth.eq.3.or.idof_meth.ge.5) then
               do i=1,n0
                  pho(i) = phimeth(i)
                  to(i) = tmeth(i)
                  phi(i) = phimeth(i)
                  t(i) = tmeth(i)
               enddo
            endif

c     1-15-04 gaz 
c     set phase of water to 1 (liquid)
            do i=1,n0
               ieos(i) = 1
               ieoso(i)  = ieos(i)
            enddo

         else if(iflg.eq.-1) then
c     
c     determine phase state for water solid-liquid-gas system
c     should count phase changes
c     
c     ieos_ch is the number of phase changes for each node
c     zero at iad = 0, this sould bbe sufficient for W,M and H  
c
            if(iad.eq.0) then
               if(l.eq.1) then
                  ieos_ch=-100
               else
                  ieos_ch = 0
               endif
            endif
            icmw = 0
            do ii=1,neq
               ij=ii+ndummy
               iced=ieos(ij)
               icedc=iced
               pl=phi(ij)
               if(iced.eq.3) then
c     gas only conditions
                  tl=t(ij)
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  if(tl.le.tliquid.and.ieos_ch(ij).lt.ich_max) then
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
                  if(sl.ge.1.0.and.ieos_ch(ij).lt.ich_max) then
c     liquid only can form                  
                     icedc=1
                     s(ij)=1.0
                     t(ij)=tliquid*eosml
                  else if(sl.le.0.0.and.ieos_ch(ij).lt.ich_max) then
c     gas only can form    
                     icedc=3
                     s(ij)=1.0
                     t(ij)=tliquid*eosmg
                  endif
               elseif(iced.eq.1) then
c     liquid only conditions
                  tl=t(ij)
                  call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  if(tl.ge.tliquid.and.ieos_ch(ij).lt.ich_max) then
c     gas-liquid can form                
                     icedc=2
                     s(ij)=eosml
                     t(ij)=tliquid
                  else if(tl.le.tsolid.and.ieos_ch(ij).lt.ich_max) then
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
                  if(sl.ge.1.0.and.ieos_ch(ij).lt.ich_max) then
c     liquid only can form                
                     icedc=1
                     s(ij)=1.0   
                     t(ij)=tsolid*eosmg
                  else if(sl.le.0.0.and.ieos_ch(ij).lt.ich_max) then
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
                  if(tl.ge.tsolid.and.ieos_ch(ij).lt.ich_max) then
c     liquid-solid can form                
                     icedc=-2
                     s(ij)=eosml 
                     t(ij)=tsolid
                  endif
               endif
               if(icedc.ne.iced) icmw = 1
               ieos_ch(ij) = ieos_ch(ij) +1
               ieos(ij)=icedc
            enddo
c     
         else if(iflg.eq.1) then
c     
c     determine phase state for methane system
c     should count phase changes
c     
            icmm = 0
            do ii=1,neq
               ij=ii+ndummy
               iced=ices(ij)
               icedc=iced
               pl=phimeth(ij)
               if(iced.eq.3) then
c     gas only conditions
                  tl=tmeth(ij)
                  call methane_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  if(tl.le.tliquid.and.ieos_ch(ij).lt.ich_max) then
c     gas-liquid can form
                     icedc=2
                     smeth(ij)=eostol
                     tmeth(ij)=tliquid
                  endif
               elseif(iced.eq.2) then
c     gas-liquid conditions
                  sl=smeth(ij)
                  call methane_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  tmeth(ij)=tliquid
                  if(sl.ge.1.0.and.ieos_ch(ij).lt.ich_max) then
c     liquid only can form                  
                     icedc=1
                     smeth(ij)=1.0
                     tmeth(ij)=tliquid*eosml
                  else if(sl.le.0.0.and.ieos_ch(ij).lt.ich_max) then
c     gas only can form    
                     icedc=3
                     smeth(ij)=1.0
                     tmeth(ij)=tliquid*eosmg
                  endif
               elseif(iced.eq.1) then
c     liquid only conditions
                  tl=tmeth(ij)
                  call methane_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  call methane_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  if(tl.ge.tliquid.and.ieos_ch(ij).lt.ich_max) then
c     gas-liquid can form                
                     icedc=2
                     smeth(ij)=eosml
                     tmeth(ij)=tliquid
                  else if(tl.le.tsolid.and.ieos_ch(ij).lt.ich_max) then
c     liquid-solid can form                
                     icedc=-2
                     smeth(ij)=eostol
                     tmeth(ij)=tsolid
                  endif
               elseif(iced.eq.-2) then
c     liquid-solid conditions
                  sl=smeth(ij)
                  call methane_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  if(sl.ge.1.0.and.ieos_ch(ij).lt.ich_max) then
c     liquid only can form                
                     icedc=1
                     smeth(ij)=1.0   
                     tmeth(ij)=tsolid*eosmg
                  else if(sl.le.0.0.and.ieos_ch(ij).lt.ich_max) then
c     solid only can form                
                     icedc=-3
                     smeth(ij)=0.0    
                     tmeth(ij)=tsolid*eosml
                  endif
               elseif(iced.eq.-3) then
c     solid only conditions
                  tl=tmeth(ij)
                  call methane_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  if(tl.ge.tsolid.and.ieos_ch(ij).lt.ich_max) then
c     liquid-solid can form                
                     icedc=-2
                     smeth(ij)=eosml 
                     tmeth(ij)=tsolid
                  endif
               endif
               if(icedc.ne.iced) then
                  icmm = 1
                  ieos_ch(ij) = ieos_ch(ij) +1
               end if
               ices(ij)=icedc
            enddo
c     
c     
         elseif(iflg.eq.-6) then
c     check hydrate state
            if(idof_meth.ne.7) then
               icmh = 0
               do ii=1,neq
                  ij=ii+ndummy
                  ihydd=ihyd(ij)
                  ihyddc=ihydd
                  tl=tmeth(ij)
                  pl=phimeth(ij)
                  w_frac = fracw(ij)
                  frach = frachyd(ij)
                  fracg = fracgas(ij)
                  call hydrate_properties(6,1,pl,dum1,dum2,dum3,
     &                 tdis1,dum5,dum6,dum7,dum8)
                  if(tl.lt.tdis1) then
c     hydrate will form if products are available

                     if(ihydd.eq.-1) then
c     only free water 

c     check for product usage 
                        if(fracg.ge.1.0d0.and.ieos_ch(ij).lt.ich_max) 
     &                       then
c     change to only free gas
                           ihyddc = -2
                        else if(fracg.ge.zero_e .and.
     &                          ieos_ch(ij).lt.ich_max) then
                           ihyddc = 0
                        endif 
                     else if(ihydd.eq.-2) then
                        
                        if(fracw(ij).ge.1.0d0 .and.
     &                       ieos_ch(ij).lt.ich_max) then
c     change to only free gas
                           ihyddc = -1
                        else if(fracw(ij).ge.zero_e .and.
     &                          ieos_ch(ij).lt.ich_max) then
                           ihyddc = 0
                        endif  
                     else if(ihydd.eq.1) then
c     check for product usage 
                        if(frach.gt.zero_e.and.
     &                       ieos_ch(ij).lt.ich_max) then
                           ihyddc= 0
                        endif
                        
                     endif
                  else if(tl.gt.tdis1) then
c     hydrate will dissociate if hydrate is available
                     ihyddc = 1

                     if(frach.gt.zero_e.and.ieos_ch(ij).lt.ich_max) then
                        
                        ihyddc= 0
                     endif

                  else if(tl.eq.tdis1) then
c     hydrate will dissociate if hydrate is available
                     ihyddc = 0
                     if(frach.le.zero_e.and.ieos_ch(ij).lt.ich_max) then
                        ihyddc = 1
                        frachyd(ij) = 0.0
                        fracgas(ij) = 1.0d0-fracw(ij)
                     endif
                  endif
                  if(ihyddc.ne.ihydd) then
                     icmh = 1
                     ieos_ch(ij) = ieos_ch(ij) +1
                  endif
                  ihyd(ij)=ihyddc
               enddo	
            else
c     equilibrium formulation
               icmh = 0
               do ii=1,neq
                  ij=ii+ndummy
                  ihydd=ihyd(ij)
                  ihyddc=ihydd
                  tl=tmeth(ij)
                  pl=phimeth(ij)
                  w_frac = fracw(ij)
                  frach = frachyd(ij)
                  fracg = fracgas(ij)
                  call hydrate_properties(6,1,pl,dum1,dum2,dum3,
     &                 tdis1,dum5,dum6,dum7,dum8)
                  if(tl.lt.tdis1) then
c     hydrate will form if products are available

                     if(ihydd.eq.-1) then
c     only free water 
                        fracg = 1.d0 - w_frac
c     check for product usage 
                        call hydrate_properties(9,ihyddc,pl,dum1,w_frac,
     &                       fracg,frach,dum5,dum6,dum7,dum8)
                        
                     else if(ihydd.eq.-2) then
c     check for product usage 
                        call hydrate_properties(9,ihyddc,pl,dum1,w_frac,
     &                       fracg,frach,dum5,dum6,dum7,dum8)
                     else if(ihydd.eq.1) then
c     check for product usage 
                        if(frach.gt.tolw.and.ieos_ch(ij).lt.ich_max) 
     &                       then
                           ihyddc= 0
                        endif
                        
                     endif
                  else if(tl.gt.tdis1) then
c     hydrate will dissociate if hydrate is available
                     ihyddc = 1

                     if(frach.gt.tolw.and.ieos_ch(ij).lt.ich_max) then
                        tmeth(ij) = tdis1
                        ihyddc= 0
                     endif

                  else if(tl.eq.tdis1) then
c     hydrate will dissociate if hydrate is available
                     ihyddc = 0
                     if(frach.le.tolw.and.ieos_ch(ij).lt.ich_max) then
                        ihyddc = 1
                        frachyd(ij) = 0.0
                        fracgas(ij) = 1.0d0-fracw(ij)
                     else if(fracgas(ij)-fracv*frachyd(ij).le.tolw
     &                       .and.ieos_ch(ij).lt.ich_max) then
                        ihyddc = -1
                     else if(fracw(ij)-fracl*frachyd(ij).le.tolw
     &                       .and.ieos_ch(ij).lt.ich_max) then
                        ihyddc = -2 
                     endif
                     
                  endif
                  if(ihyddc.ne.ihydd) then
                     icmh = 1
                     ieos_ch(ij) = ieos_ch(ij) +1
                  endif
                  ihyd(ij)=ihyddc
               enddo		   	     	      
            endif
         elseif(iflg.eq.2) then
c     
c     update variables with N-R corrections           
c     
c     
c     strd is passed through common 
c     
            if(max(icmw,icmm,icmh).le.0) then
               strd = strd_meth0   
            else
               strd = strd_meth 
            endif
            if(iad.gt.abs(maxit)/2) then
               strd = strd_meth
            endif  
c     
            if(idof_meth.lt.5) then
c     component #1         

               nr1=nrhs(1)
               nr2=nrhs(2)
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  ieosd=ieos(i)
                  if(ps(i).eq.0.0) then
                     t(i)=t(i)-bp(i2)*strd
                  elseif(ieosd.eq.1) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                  elseif(abs(ieosd).eq.2) then
                     phi(i)=phi(i)-bp(i1)*strd
                     s(i)=s(i)-bp(i2)*strd
                  elseif(abs(ieosd).eq.3) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                  endif
               enddo

c     
c     component #2         
c     
c     strd is passed through common
               nr1=nrhs(3)
               nr2=nrhs(4)
               nr3=nrhs(5)
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  ieosd=ices(i)
                  if(ps(i).eq.0.0) then
                     tmeth(i)=tmeth(i)-bp(i2)*strd
                  elseif(ieosd.eq.1) then
                     phimeth(i)=phimeth(i)-bp(i1)*strd
                     tmeth(i)=tmeth(i)-bp(i2)*strd
                  elseif(abs(ieosd).eq.2) then
                     phimeth(i)=phimeth(i)-bp(i1)*strd
                     smeth(i)=smeth(i)-bp(i2)*strd
                  elseif(abs(ieosd).eq.3) then
                     phimeth(i)=phimeth(i)-bp(i1)*strd
                     tmeth(i)=tmeth(i)-bp(i2)*strd
                  endif
                  frachyd(i) = frachyd(i) - bp(i3)*strd
               enddo
            else if(idof_meth.eq.5) then
c     
c     both component done together (special case)
c     idof_meth = 5
c     
c     strd is passed through common
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
               nr4=nrhs(4)
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  i4=i+nr4
                  if(ps(i).eq.0.0) then
                     tmeth(i)=tmeth(i)-bp(i2)*strd
                  else
                     phimeth(i)=phimeth(i)-bp(i1)*strd
                     phi(i) = phimeth(i)
                     tmeth(i)=tmeth(i)-bp(i2)*strd
                     t(i) = tmeth(i) 
                     fracw(i) = fracw(i) - bp(i3)*strd
c gaz 7-25-2006
                     fracw(i) = min(1.0d0,max(0.0d0,fracw(i)))
                     frachyd(i) = frachyd(i) - bp(i4)*strd
c gaz 7-25-2006
                     frachyd(i) = min(1.0d0,max(0.0d0,frachyd(i)))
                  endif
               enddo
            else if(idof_meth.eq.7) then
c     
c     equilibrium model
c     idof_meth = 7
c     
c     strd is passed through common
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
               nr4=nrhs(4)
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  i4=i+nr4
                  if(ps(i).eq.0.0) then
                     tmeth(i)=tmeth(i)-bp(i2)*strd
                  else
                     if(ihyd(i).eq.-1) then
                        phimeth(i)=phimeth(i)-bp(i1)*strd
                        phi(i) = phimeth(i)
                        tmeth(i)=tmeth(i)-bp(i2)*strd
                        t(i) = tmeth(i) 
                        fracw(i) = fracw(i) - bp(i3)*strd
                        fracgas(i) = 1.d0 - fracw(i)
c     calculate hydrate
                        call hydrate_properties(8,-1,pl,tl,fracw(i),
     &                       fracgas(i),frach,dum5,dum6,dum1,dum2) 
                        frachyd(i) = frach
                     else if(ihyd(i).eq.-2) then
                        phimeth(i)=phimeth(i)-bp(i1)*strd
                        phi(i) = phimeth(i)
                        tmeth(i)=tmeth(i)-bp(i2)*strd
                        t(i) = tmeth(i) 
                        fracgas(i) = fracgas(i) - bp(i3)*strd
                        fracw(i) = 1.d0 - fracgas(i)
c     calculate hydrate
                        call hydrate_properties(8,-2,pl,tl,fracw(i),
     &                       fracgas(i),frach,dum5,dum6,dum1,dum2) 
                        frachyd(i) = frach
                     else if(ihyd(i).eq.0) then
                        phimeth(i)=phimeth(i)-bp(i1)*strd
                        phi(i) = phimeth(i)
                        fracw(i) = fracw(i) - bp(i3)*strd
c     fracw(i) = min(1.d0,fracw(i))
                        fracgas(i) = 1.d0 - fracw(i)
c     calculate temperature
                        call hydrate_properties(6,0,phi(i),
     &                       0.,0.,0.,tmeth(i),dums1,0.,0.,0.)
                        t(i) = tmeth(i)
                        frachyd(i) = frachyd(i) - bp(i2)*strd
                     else if(ihyd(i).eq.1) then
                        phimeth(i)=phimeth(i)-bp(i1)*strd
                        phi(i) = phimeth(i)
                        tmeth(i)=tmeth(i)-bp(i2)*strd
                        t(i) = tmeth(i) 
                        fracw(i) = fracw(i) - bp(i3)*strd
                        frachyd(i) = 0.0
                     endif
                  endif
               enddo
            endif

         elseif(iflg.eq.-33) then
c     call EOS routines methane hydrate
            call ther_meth_h2o(3,ndummy)
         elseif(iflg.eq.-34) then
c     call EOS routine to allocate space
            call ther_meth_h2o(0,ndummy)
         elseif(iflg.eq.-35) then
c     call EOS routine to deallocate space
            call ther_meth_h2o(-1,ndummy)
         elseif(iflg.eq.-3) then
c     call EOS routines methane
            call ther_meth_h2o(2,ndummy)
         elseif(iflg.eq.3) then
c     call EOS routines water
            call ther_meth_h2o(1,ndummy)
         elseif(iflg.eq.-4) then
c     first ckeck for temperature-specified source terms
c     need to comment out call in main routine for water if meth is activated
            if(iad.eq.0) then
c     only for first iteration
               do ii=1,neq
                  ij=ii+ndummy
                  if(eskmeth(ij).lt.0.0) then
                     call methane_properties(9,1,phimeth(ij),
     &                    -eskmeth(ij),dum2,dum3,ensrc,dum4,dum5,dum6,
     &                    dum7)
                     eflowmeth(ij)=ensrc
                  else
                     eflowmeth(ij)=eskmeth(ij)
                  endif
                  if(esk(ij).lt.0.0) then
                     call h2o_properties(9,1,phi(ij),-esk(ij),dum2,dum3,
     &                    ensrc,dum4,dum5,dum6,dum7)
                     eflow(ij)=ensrc
                  else
                     eflow(ij)=esk(ij)
                  endif 
               enddo
            endif                
         elseif(iflg.eq.4) then
c     call equation generation and load jacobian array
            if(idof_meth.ne.7) then
               call gensmethh2o    
            else
               call gensmethh2o_equil    
            endif  
         elseif(iflg.eq.5) then
            if(ntty.eq.2) then
c     
c     output for methane information
c     
               if (iout .ne. 0) write(iout,803)
               if(iatty.ne.0) write(iatty,803)
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
                     if (iout .ne. 0) write(iout,600) il
                     if (iatty.gt.0) write(iatty,600) il
 600                 format(2x,'Matrix Level = ',i1)
                  endif
                  do i=1,mlev
                     md=  nskw(i+(il-1)*mlev)
                     if(idof_meth.ne.7) then
                        w_frac = fracw(md) + fracl*frachyd(md)
                        fracg = 1.d0 - w_frac 
                        phase_frac = fracg - fracv*frachyd(md)
                     else
                        if(ihyd(md).eq.0) then
                           w_frac = fracw(md)
                           fracg=1.d0-w_frac
                           phase_frac=fracg-fracv*frachyd(md)
                        else if(ihyd(md).eq.-2) then
                           fracg = fracgas(md)
                           phase_frac=fracg-fracv*frachyd(md)
                           w_frac = 1.d0-fracg
                        else if(ihyd(md).eq.-1) then
                           w_frac = fracw(md)	                  
                           call hydrate_properties(8,ihyd(md),pl,tl,
     &                          w_frac,fracg,frach,dum5,dum6,dum7,dum8) 
                           fracg = 1. - w_frac	                    
                           phase_frac = 0.d0
                        else if(ihyd(md).eq.1) then
                           w_frac = fracw(md)
                           fracg=1.d0-w_frac
                           phase_frac=fracg-fracv*frachyd(md)
                        endif
                     endif

                     skmd = skmeth(md)
                     qhmd = qhmeth(md)
                     if (iout .ne. 0) write(iout,804) 
     &                    md,w_frac,fracg,phase_frac,ices(md),
     &                    skmd,qhmd       
                     if(iatty.ne.0)
     &                    write(iatty,804) 
     &                    md,w_frac,fracg,phase_frac,ices(md), 
     &                    skmd,qhmd       
                  enddo
               enddo
 803           format(/,20x,'Nodal Information (Methane)  ',//,2x,'Node'
     &              ,3x,'water(tot)',2x,'meth(tot)',1x,' meth(free)',
     &              1x,'phase     src/snk    E src/snk')
 804           format(i6,2x,f9.3,1x,f8.3,1x,f10.3,4x,i3,1x,
     &              g11.3,2x,g11.3)
c     
c**** printout global mass and energy flows ****
c**** printout global mass and energy balances ****
c     
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,700)
               if (iatty.gt.0) write(iatty,700)
 700           format(/,20x,'Global Mass & Energy for Methane')
               if (ntty.eq.2.and. iout .ne. 0) write(iout,701) 
     &              ammeth, aemeth
               if (iatty.gt.0) write(iatty,701) ammeth, aemeth      
 701           format(1x,'Total mass in system at this time:          ',
     &              e14.6,' kg',/,1x,
     &              'Total enthalpy in system at this time:      ',
     &              e14.6,' MJ')
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,702)
               if (iatty.gt.0) write(iatty,702)
 702           format(/,20x,'Global Mass & Energy flows for Methane')
               if (ntty.eq.2 .and. iout .ne. 0)write(iout,703) 
     &              qmeth,abs(qmethts-qmethts_in)/dtotdm,
     &              abs(qmethts_in)/dtotdm,
     &              qemeth,abs(qemethts-qemethts_in)/dtotdm,
     &              abs(qemethts_in)/dtotdm
               if (iatty.gt.0) write(iatty,703) 
     &              qmeth,abs(qmethts-qmethts_in)/dtotdm,
     &              abs(qmethts_in)/dtotdm,
     &              qemeth,abs(qemethts-qemethts_in)/dtotdm,
     &              abs(qemethts_in)/dtotdm
 703           format(1x,
     &              'Mass flux:       Total      Time Step(in)',
     &              '       Time Step(out) ',/,
     &              9x,e14.6,5x,e14.6,7x,e14.6,' kg/s',/,1x,
     &              'Energy  flux:    Total      Time Step(in)',
     &              '       Time Step(out) ',/,
     &              9x,e14.6,5x,e14.6,7x,e14.6,' MJ/s')      
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,704)
               if (iatty.gt.0) write(iatty,704)
 704           format(/,20x,'Global Mass & Energy balances for Methane')
               if(balemeth.ne.-999999.) then
                  if (ntty.eq.2 .and. iout .ne. 0)
     &                 write(iout,705) balmeth,balemeth   
                  if (iatty.gt.0) write(iatty,705) balmeth,balemeth
               else 
                  if (ntty.eq.2 .and. iout .ne. 0)
     &                 write(iout,706) balmeth   
                  if (iatty.gt.0) write(iatty,706) balmeth   
               endif
 705           format(1x,'Mass balance:',1x,e14.6,
     &              5x,'Energy balance:',1x,e14.6)
 706           format(1x,'Mass balance:',1x,e14.6,
     &              5x,'Energy balance: N/A')
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,709) ammeth0
               if (iatty.gt.0) write(iatty,709) ammeth0
 709          format(/,1x,'Initial Mass of Methane (kg) in the System: '
     &              ,e14.6)
               
               if (ntty.eq.2 .and. iout .ne. 0)
     &             write(iout,707) (skmd1-skmd10)
               if (iatty.gt.0) write(iatty,707) (skmd1-skmd10)
 707           format(/,1x,'Total methane produced (kg) at This Time ', 
     &              'Step ', e14.6)
               
               if (ntty.eq.2 .and. iout .ne. 0) write(iout,708) skmd1
               if (iatty.gt.0) write(iatty,708) skmd1
 708           format(/,1x,'Net kg methane discharge ',
     &              '(total out-total in): ', e14.6)

               if (ntty.eq.2 .and. iout .ne. 0)
     &              write(iout,905)ic_hyd, ic_hyd_tot
               if(iatty.gt.0) 
     &              write(iatty,905)ic_hyd, ic_hyd_tot
 905           format(/,' Hydrate dissociation changes this Time Step ',
     &              i7, ' Total ', i7)
               if(ntty.eq.2 .and. iout .ne. 0) write(iout,906) numhyd
               if(iatty.gt.0) write(iatty,906) numhyd
 906           format(/,' Number of nodes with hydrate persent ',i7)
            endif
         elseif(iflg.eq.-5) then
            if(ntty.eq.2) then
c     
c     output for hydrate information
c     
               if (iout .ne. 0) write(iout,903)
               if(iatty.ne.0) write(iatty,903)
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
                     if (iout .ne. 0) write(iout,600) il
                     if (iatty.gt.0) write(iatty,600) il
                  endif
                  do i=1,mlev
                     md=  nskw(i+(il-1)*mlev)
                     hyd_frac = frachyd(md)
                     w_frac = fracw(md)
                     skmd = rathyd*skhyd(md)
                     qhmd = rathyd*qhhyd(md)
                     ihyddc = ihyd(md)
                     if(idof_meth.ne.7) then
                        if(ihyd(md).eq.-1) then
                           w_frac  =fracw(md)  
                        else if(ihyd(md).eq.-2) then
                           w_frac = fracw(md) 
                        else if(ihyd(md).eq.0) then
                           w_frac = fracw(md) 
                        else if(ihyd(md).eq.1) then
                           w_frac = fracw(md) 
                        endif
                     else 
                        if(ihyd(md).eq.-1) then 
                           w_frac =fracw(md)-hyd_frac*fracl
                        else if(ihyd(md).eq.-2) then
                           w_frac = 0.0d00
                        else if(ihyd(md).eq.0) then
                           w_frac = fracw(md)-hyd_frac*fracl
                        else if(ihyd(md).eq.1) then
                           w_frac = fracw(md)-hyd_frac*fracl
                        endif
                     endif                      
                     if (iout .ne. 0) write(iout,904) 
     &                    md,w_frac,hyd_frac,skmd,qhmd,ihyddc
                     if(iatty.ne.0)
     &                    write(iatty,904) 
     &                    md,w_frac,hyd_frac,skmd,qhmd,ihyddc
                  enddo
               enddo
 903           format(/,20x,'Nodal Information (hydrate)  ',//,2x,'Node'
     &              ,2x,'water(free)',1x,'hydrate frac ',1x,
     &              'hydrate   src/snk    E src/snk',' hydrate state')
 904           format(i6,4x,f9.3,4x,f9.3,12x,g11.3,2x,g11.3,4x,i6)
c     
            endif
         elseif(iflg.eq.6) then
c     initialize variables
            do i=1,n
               if(ices(i).eq.-2) then
                  call methane_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  tmeth(i)=tsolid
               else if(ices(i).eq.2) then
                  call methane_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  tmeth(i)=tliquid
               else
                  tmeth(i)= tometh(i)
               endif
               phometh(i) = phimeth(i)
               someth(i)= smeth(i)
            enddo

            do i=1,n
               if(ieos(i).eq.-2) then
                  call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &                 tsolid,dtps,dum4,dum5,dum6)
                  t(i)=tsolid
               else if(ieos(i).eq.2) then
                  call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &                 tliquid,dtps,dum4,dum5,dum6)
                  t(i)=tliquid
               else
                  t(i)= to(i)
               endif
               phi(i) = pho(i)
               so(i)= s(i)
            enddo

         elseif(iflg.eq.7) then

c     calulate initial mass and energy of methane 
c     calulate initial mass and energy of methane hydrate 

            ammeth0 = 0.0
            aemeth0 = 0.0
            qmeth = 0.0
            qemeth = 0.0
            qmeth_in = 0.0
            qemeth_in = 0.0
            skmd1 = 0.0
            amhyd0 = 0.0
            aehyd0 = 0.0
            qhyd = 0.0
            qehyd = 0.0
            qhyd_in = 0.0
            qehyd_in = 0.0
            do i=1,n
               denmethh(i) = denmethi(i)*dtot
               denemethh(i) = denemethi(i)*dtot
               denmethi(i) = 0.0 
               denemethi(i) = 0.0 
               ammeth0 = ammeth0 + denmethh(i)*volume(i)
               aemeth0 = aemeth0 + denemethh(i)*volume(i)
               denhydh(i) = denhydi(i)*dtot
               denehydh(i) = denehydi(i)*dtot
               amhyd0 = amhyd0 + denhydh(i)*volume(i)
               aehyd0 = aehyd0 + denehydh(i)*volume(i)
            enddo

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
c     
c     modify global balance terms for equilibrium
c     
            if(idof_meth.eq.7) then
               ammeth0 = ammeth0 + fracv*amhyd0
               aemeth0 = aemeth0 + fracv*aehyd0
               amh2o0 = amh2o0 + fracl*amhyd0
               aeh2o0 = aeh2o0 + fracl*aehyd0
               aehyd0 = 0.0d0
            endif
c     initialize space for methane arrays

            nmat(37) = nmat(36)+nelm(neq+1)-(neq+1)
            nrhs(7)  = nrhs(6) +neq

         elseif(iflg.eq.8) then

c     calculate current mass and energy of methane 
c     calculate current fluxes in and out of model  
c     note : for the current timestep iflg=8 
c     must be called after iflg=9
            skmd10 = skmd1
            ammeth = 0.0
            aemeth = 0.0
            qmethts = 0.0
            qemethts = 0.0
            qmethts_in = 0.0
            qemethts_in = 0.0
            do i=1,n
               ammeth = ammeth + denmethh(i)*volume(i)
               skmd1 = skmd1 + skmeth(i)*dtotdm
               dums1 = (skmeth(i)+skmhyd(i))*dtotdm
               qmeth = qmeth + dums1
               aemeth = aemeth + denemethh(i)*volume(i)
               qmethts = qmethts + dums1
               if(skmeth(i)+skmhyd(i).lt.0.0) then
                  qmeth_in = qmeth_in + dums1           
                  qmethts_in = qmethts_in + dums1                
               endif
               dums1 = (qhmeth(i)+qhmhyd(i))*dtotdm
               qemeth = qemeth + dums1                
               qemethts = qemethts + dums1                
               if(qhmeth(i)+qhmhyd(i).lt.0.0) then
                  qemeth_in = qemeth_in + dums1                
                  qemethts_in = qemethts_in + dums1                
               endif
            enddo

c     calculate current mass and energy of hydrate     
c     calculate current fluxes in and out of model  
c     note : for the current timestep iflg=8 
c     must be called after iflg=9

            amhyd = 0.0
            aehyd = 0.0
            qhydts = 0.0
            qehydts = 0.0
            qhydts_in = 0.0
            qehydts_in = 0.0
            do i=1,n
               amhyd = amhyd + denhydh(i)*volume(i)
               aehyd = aehyd + denehydh(i)*volume(i)
               dums1 = skhyd(i)*dtotdm
               qhyd = qhyd + dums1           
               qhydts = qhydts + dums1           
               if(skhyd(i).lt.0.0) then
                  qhyd_in = qhyd_in + dums1           
                  qhydts_in = qhydts_in + dums1                
               endif
               dums1 = qhhyd(i)*dtotdm
               qehyd = qehyd + dums1                
               qehydts = qehydts + dums1                
               if(qhhyd(i).lt.0.0) then
                  qehyd_in = qehyd_in + dums1                
                  qehydts_in = qehydts_in + dums1                
               endif
            enddo

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
               dums1 = (sk(i)+skwhyd(i))*dtotdm
               qh2o = qh2o + dums1           
               qh2ots = qh2ots + dums1           
               if(sk(i)+skwhyd(i).lt.0.0) then
                  qh2o_in = qh2o_in + dums1           
                  qh2ots_in = qh2ots_in + dums1                
               endif
               dums1 = (qh(i)+qhwhyd(i))*dtotdm
               qeh2o = qeh2o + dums1                
               qeh2ots = qeh2ots + dums1                
               if(qh(i)+qhwhyd(i).lt.0.0) then
                  qeh2o_in = qeh2o_in + dums1                
                  qeh2ots_in = qeh2ots_in + dums1                
               endif
            enddo

c     
c     modify global balance terms for equilibrium
c     
            if(idof_meth.eq.7) then
               ammeth = ammeth + fracv*amhyd
               aemeth = aemeth + fracv*aehyd
               amh2o = amh2o + fracl*amhyd
               aeh2o = aeh2o + fracl*aehyd
               aehyd = 0.0d0
            endif

c     calculate mass and energy balance for methane

            amaxflx = max(abs(qmeth_in),abs(qmeth-qmeth_in))
            if(amaxflx.gt.zero_e) then
               balmeth = (ammeth - ammeth0 + qmeth) / amaxflx
            else if(ammeth0.gt.0.0) then
               balmeth = (ammeth - ammeth0 + qmeth) / ammeth0
            else 
c     set methane mass balance to n/a
               balmeth = -999999.0                                    
            endif

            amaxflx = max(abs(qemeth_in),abs(qemeth-qemeth_in))
            if(amaxflx.gt.zero_e) then
               balemeth = (aemeth - aemeth0 + qemeth) / amaxflx
            else if(abs(aemeth0).gt.zero_e) then
               balemeth = (aemeth - aemeth0 + qemeth) / aemeth0
            else 
c     set methane energy balance to n/a
               balemeth = -999999.0                                    
            endif

c     calculate mass and energy balance for water

            amaxflx = max(abs(qh2o_in),abs(qh2o-qh2o_in))
            if(amaxflx.gt.zero_e) then
               balh2o = (amh2o - amh2o0 + qh2o) / amaxflx
            else if(amh2o0.gt.0.0) then
               balh2o = (amh2o - amh2o0 + qh2o) / amh2o0
            else 
c     set h2o mass balance to n/a
               balh2o = -999999.0                                    
            endif
            difm = balh2o

            amaxflx = max(abs(qeh2o_in),abs(qeh2o-qeh2o_in))
            if(amaxflx.gt.zero_e) then
               baleh2o = (aeh2o - aeh2o0 + qeh2o) / amaxflx
            else if(aeh2o0.gt.0.0) then
               baleh2o = (aeh2o - aeh2o0 + qeh2o) / aeh2o0
            else 
c     set h2o energy balance to n/a
               baleh2o = -999999.0                                    
            endif

c     correct mass and energy for water in combining equations

            if(idof_meth.eq.3) then

c     set methane energy balance to n/a
               balemeth = -999999.0                                    

               amaxflx = max(abs(qeh2o_in+qemeth_in),
     &              abs((qeh2o+qemeth)-(qeh2o_in+qemeth_in)))
               if(amaxflx.gt.zero_e) then
                  baleh2o = ((aeh2o+aemeth) - 
     &                 (aeh2o0+aemeth0) + (qeh2o+qemeth)) / amaxflx
               else if(aeh2o0+aemeth0.gt.0.0) then
                  baleh2o = ((aeh2o+aemeth) - (aeh2o0+aemeth0) 
     &                 + (qeh2o+qemeth))/(aeh2o0+aemeth0)
               else 
c     set h2o energy balance to n/a
                  baleh2o = -999999.0
               endif

c     correct h2o energy balance

               dife = baleh2o

            endif
            if(idof_meth.ge.5) then

c     set methane energy balance to n/a
c     baleh2o represents energy balance of water,methane, and hydrate

               balemeth = -999999.0                                    

               amaxflx = max(abs(qeh2o_in+qemeth_in+qehyd_in),abs((qeh2o
     &              +qemeth+qehyd)-(qeh2o_in+qemeth_in+qehyd_in)))
               amaxener = (aeh2o+aemeth+aehyd) - (aeh2o0+aemeth0+aehyd0)
c     if(amaxflx.gt.zero_e.and.amaxener.gt.zero_e) then
               if(amaxener*amaxflx.gt.zero_e) then
                  baleh2o = ((aeh2o+aemeth+aehyd) - (aeh2o0+aemeth0
     &                 +aehyd0) + (qeh2o+qemeth+qehyd)) / amaxflx
               else if(abs(aeh2o0+aemeth0+aehyd0).gt.zero_e) then
                  baleh2o = ((aeh2o+aemeth+aehyd) - 
     &                 (aeh2o0+aemeth0+aehyd0) + (qeh2o+qemeth+qehyd)) 
     &                 / (aeh2o0+aemeth0+aehyd0)
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
c     check for dissociation changes
c     count cells with hydrate present
c     
            ic_hyd = 0
            numhyd = 0
            do i = 1,n0
               if(ihyd(i).ne.ihydo(i)) then
                  ic_hyd = ic_hyd + 1
               endif
               if(frachyd(i).gt.tolw) then
                  numhyd = numhyd + 1 
               endif
            enddo
            ic_hyd_tot =  ic_hyd_tot + ic_hyd

            do i = 1,n
               phometh(i) = phimeth(i)
               tometh(i) = tmeth(i)
               someth(i) = smeth(i)
               fracwo(i) = fracw(i)
               frachydo(i) = frachyd(i)
               fracgaso(i) = fracgas(i)
               denmethh(i) = denmethh(i) + denmethi(i)*dtot
               denemethh(i) =  denemethh(i) + denemethi(i)*dtot
               denmethi(i) = 0.0
               denemethi(i) = 0.0 
               ieoso(i) = ieos(i)
               iceso(i) = ices(i)
               ihydo(i) = ihyd(i)
            enddo

            do i = 1,n
               denhydh(i) = denhydh(i) + denhydi(i)*dtot
               denehydh(i) =  denehydh(i) + denehydi(i)*dtot
            enddo

         elseif(iflg.eq.10) then

c     permeability reduction factor
c     temma add k=ko(1-SH/SH_max)^N 2005/11/04

            if(nhyd(1).eq.0) then
               do i=1, n
                  if (nhyd(i).eq. 0 .and. permhyd0(i).eq.-3) then
                     normal_hyd = frachyd(i)/frachyd_max(i)
                     if (normal_hyd.le.0.0) then
                        normal_hyd = 0.0
		        permhyd(i) = 1.0     
                        dpermhyd4(i) = -nhyd_3(i)/frachyd_max(i)
                     else if (normal_hyd.ge.1.0) then
                        normal_hyd = 1.0
                        permhyd(i) = 0.0
                        dpermhyd4(i) = 0.0
                     else
		        permhyd(i) = (1.0-normal_hyd)**nhyd_3(i)     
                        dpermhyd4(i) = -nhyd_3(i)*(1/frachyd_max(i))*
     &                       (1.0-normal_hyd)**(nhyd_3(i)-1)
                     endif	         
                  end if

               end do
	    end if
ccc                  
            if(nhyd(1).gt.0) then
               do i = 1,n
c                  if(nhyd(i).gt.0) then
c   tenma 04/14/2005 => AIST parallel layer (Sand & Mud)
c   methperm ( permhyd0 = 1.0 )
                  if(nhyd(i).gt.0 .and. permhyd0(i).ne.-2) then
                     if (nhyd_2(i).eq.-1) then
                        permhyd(i) = permhyd0(i)*(1.0-frachyd(i))**
     &                       nhyd(i)
                        dpermhyd4(i) = -nhyd(i)*permhyd0(i)*
     &                       (1.0-frachyd(i))**(nhyd(i)-1)
c     else if(nhyd_2(i).ne.0) then
c temma should check this condition
                     else
                        if(frachyd(i) .gt. frachyd_1(i)) then
                           permhyd(i) = permhyd0(i)*((1.0-frachyd_1(i))
     &                          **nhyd(i))*((1.0-frachyd(i))/
     &                          (1.0-frachyd_1(i)))**nhyd_2(i)      
                           dpermhyd4(i) = -nhyd_2(i)*permhyd0(i)*
     &                          ((1.0-frachyd_1(i))**(nhyd(i)-1))
     &                          *((1.0-frachyd(i))/(1.0-frachyd_1(i)))
     &                          **(nhyd_2(i)-1)
                        else 
                           permhyd(i) = permhyd0(i)*(1.0-frachyd(i))
     &                          **nhyd(i)     
                           dpermhyd4(i) = -nhyd(i)*permhyd0(i)*
     &                          (1.0-frachyd(i))**(nhyd(i)-1) 
                        endif
                     end if
c    parallel layer permeability reduction factor (Sand & Mud )
c    methperm (permhyd0 = -2.0 )
                  else if(nhyd(i).gt.0 .and. permhyd0(i).eq.-2) then
                     temp_kz_s = perm_sand(i) * z_sand(i)
                     temp_kz_m = perm_mud(i) * z_mud(i)
                     temp_z = z_sand(i) + z_mud(i)
                     temp_poroz = poro_sand(i) * z_sand(i)
                     para_mix = z_sand(i)*poro_sand(i) + 
     &                    z_mud(i)*poro_mud(i)
                     temp_fh=(temp_poroz-para_mix*frachyd(i))/temp_poroz

c
                     ka0_xy(i) = (temp_kz_s + temp_kz_m ) / temp_z
                     permhyd(i) = (temp_kz_s * (temp_fh)**nhyd(i)
     &                    + temp_kz_m ) / ( temp_kz_s + temp_kz_m )     
                     dpermhyd4(i) = -nhyd(i)*(para_mix/temp_poroz)*
     &                    (temp_kz_s/(temp_kz_s + temp_kz_m ))*
     &                    (temp_fh)**(nhyd(i)-1)
                     permhyd_z(i) = 1.0/ permhyd(i) * ((z_sand(i) + 
     &                    z_mud(i)) * perm_sand(i) * perm_mud(i) *
     &                    (temp_fh)**nhyd(i))/( z_sand(i)*perm_mud(i) 
     &                    + z_mud(i)*perm_sand(i)*(temp_fh)**nhyd(i)) 
c   tenma 04/14/2005 => AIST parallel layer (Sand & Mud)
                  else if(nhyd(i).le.0) then 
                     permhyd(i) = 1.
                     dpermhyd4(i) =  0.0
                  endif
               enddo
            else if(nhyd(1).lt.0) then
               np_hyd = abs(nhyd(1)) 
               do i = 1,n
                  fh = frachyd(i)
                  power = ahyd_coef(1)
                  do j= 1 , np_hyd
                     power = power + ahyd_coef(j+1)*fh**j
                  enddo
                  dpowerh = ahyd_coef(2)
                  do j= 2 , np_hyd
                     dpowerh = dpowerh + j*ahyd_coef(j+1)*fh**(j-1)
                  enddo
                  permhyd(i) = permhyd0(i)*(1.0-frachyd(i))**power  
                  dpermhyd4(i) = -power*permhyd0(i)*
     &                 (1.0-frachyd(i))**(power-1)*dpowerh
               enddo
            endif

         elseif(iflg.eq.11) then

c     reset variables

            do i = 1,n
               phimeth(i) = phometh(i)
               phi(i) = pho(i)
               t(i) = to(i)
               tmeth(i) = tometh(i)
               smeth(i) = someth(i)
               s(i) = so(i)
               fracw(i) = fracwo(i)
               frachyd(i) = frachydo(i)
               fracgas(i) = fracgaso(i)
               denmethi(i) = 0.0
               denemethi(i) = 0.0 
               ieos(i)= ieoso(i)
               ices(i)= iceso(i)
               ihyd(i)= ihydo(i)
            enddo
            

         elseif(iflg.eq.-20) then

c     write out restart file

            write(isave ,6002)  (max(to(mi),tolw),   mi=1,n )
            write(isave ,6002)  (max(so(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(pho(mi),tolw),  mi=1,n )
            
            write(isave ,6002)  (max(tometh(mi),tolw),   mi=1,n )
            write(isave ,6002)  (max(someth(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(phometh(mi),tolw),  mi=1,n )

            write(isave ,6002)  (max(fracwo(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(fracgaso(mi),tolw),    mi=1,n )
            write(isave ,6002)  (max(frachydo(mi),tolw),    mi=1,n )

            write(isave ,6003)  (ieoso(mi),    mi=1,n )
            write(isave ,6003)  (iceso(mi),    mi=1,n )
            write(isave ,6003)  (ihydo(mi),    mi=1,n )
            
 6002       format(4g25.16)
 6003       format(25i4)
         elseif(iflg.eq.20) then

c     read restart file (will not work for double porosity)

            read(iread ,*)  (to(mi), mi=1,n )
            read(iread ,*)  (so(mi), mi=1,n )
            read(iread ,*)  (pho(mi), mi=1,n )

            read(iread ,*)  (tometh(mi),  mi=1,n )
            read(iread ,*)  (someth(mi),  mi=1,n )
            read(iread ,*)  (phometh(mi), mi=1,n )

            read(iread ,*)  (fracwo(mi),  mi=1,n )
            read(iread ,*)  (fracgaso(mi),  mi=1,n )
            read(iread ,*)  (frachydo(mi),  mi=1,n )

            read(iread ,*)  (ieoso(mi),  mi=1,n )
            read(iread ,*)  (iceso(mi),  mi=1,n )
            read(iread ,*)  (ihydo(mi),  mi=1,n )
            t = to
            s = so
            phi = pho
            smeth = someth
            tmeth = tometh
            phimeth = phometh
            fracw = fracwo
            frachyd = frachydo
            fracgas = fracgaso
            ieos = ieoso
            ices = iceso
            ihyd = ihydo
            numhyd = 0
            do i=1,n0
               if(frachyd(i).gt.tolw) then
                  numhyd = numhyd + 1 
               endif
            enddo 
            
         elseif(iflg.eq.21) then

c     write surfer or tecplot contour files

            write(iscon,790)
 790        format('         ','X',12x,'Y',12x,'Z',6x,'   pres',12x,'t',
     &           8x,'fracg',8x,'frach',11x,'sk',7x,'skmeth')
            j = ndummy
            if(idof_meth.ne.7) then
               do k=1,neq_primary
                  if(izone_surf_nodes(k).eq.j) then
                     write(iscon,855) cord(k,1),cord(k,2),cord(k,3),
     &                    phi(k),t(k),1.-frachyd(k)-fracw(k),
     &                    frachyd(k),sk(k),skmeth(k)
                  endif
               enddo
            else
               do k=1,neq_primary
                  if(izone_surf_nodes(k).eq.j) then 
                     write(iscon,855) 
     &                    cord(k,1),cord(k,2),cord(k,3),phi(k),
     &                    t(k),fracgas(k)-fracv*frachyd(k),frachyd(k),
     &                    sk(k),skmeth(k)
                  endif
               enddo	        
            endif
 855        format(1x,9(1x,g12.5))
         endif

      endif

      return
      end

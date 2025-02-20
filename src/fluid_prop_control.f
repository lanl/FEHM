      subroutine fluid_props_control(iflg, istart, iend, 
     &           fluid_type, prop, phase)
c      
c subroutine to manage fluid properties with 
c gaz 070321 initial coding
c
c gaz 052522 define local variables
c iflg = 0, allocate memory
c iflg = 1, loop for property evaluation
c iflg = -1, initialize metrics  
c iflg = -2, finish metrics   
c istart = node number to start property calculation
c iend = node number to finish  property calculation 
c fluid_type = h2o, air, co2      
c prop = property to calculate (den,vis,enth,all) for ngas always 'all'
c phase =  phase, ngas always is 'gas'               
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, iwater_table,
     &       iair_table, num_eos_table
          implicit none
             
          integer iflg, istart, iend, ii
          character*9 prop, phase, fluid_type
c gaz 052522 describe  subroutine pass-through variables
c iflg : flag passed to fluid prop routines          
          if(fluid_type(1:3).eq.'h2o') then
            if(iwater_table.eq.0) then
             call h2o_props_polynomials(iflg, istart, iend, 
     &       prop, phase)
            else
             call h2o_props_table(iflg, istart, iend, 
     &       prop, phase)
            endif
          else if(fluid_type(1:3).eq.'air') then
            if(iair_table.eq.0) then
             call air_props(iflg, istart, iend, prop, phase)
            else
             call air_props_table(iflg, istart, iend, prop, phase)
            endif              
          else if(fluid_type(1:3).eq.'co2') then
             call co2wh_props_table(iflg, istart, iend, prop, phase)
          endif
       return 
       end
          
      subroutine h2o_props_polynomials(iflg, istart, iend, prop, phase)
c      
c subroutine to evaluate water properties using rational polynomials
c gaz 060621 initial coding
c      
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, iout, p_tol, t_tol,
     &     pc_tol,roc0, iptty, l, iad, eval_test_h2o
          use combi, only: cord
          use comci, only: cden
          use comdi, only: phi, t, iieos, ieos
          use comdti, only: n0
          use comii, only: cel, cev, crl, crv, cvl, cvv
          use comfi, only: pci

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate, n_allocate_iso
          character*9 prop, phase

          integer ndummy,iieosl,mid,mi,ieosd,iieosd,kq

      	  real*8 psatl,dtsatp,dpsats
      	  real*8 dtin,dporpl,dportl,xrl,xrv,drl,drv,drlp,drvp,ela0,elpa1
      	  real*8 elpa2,elpa3,elta1,elta2,elta3,elpta,elp2ta,elpt2a
      	  real*8 elb0,elpb1,elpb2,elpb3,eltb1,eltb2,eltb3,elptb
      	  real*8 elp2tb,elpt2b,dla0,dlpa1,dlpa2,dlpa3,dlta1,dlta2,dlta3
      	  real*8 dlpta,dlp2ta,dlpt2a,dlb0,dlpb1,dlpb2,dlpb3,dltb1,dltb2
      	  real*8 dltb3,dlptb,dlp2tb,dlpt2b,vla0,vlpa1,vlpa2,vlpa3,vlta1
      	  real*8 vlta2,vlta3,vlpta,vlp2ta,vlpt2a,vlb0,vlpb1,vlpb2,vlpb3
      	  real*8 vltb1,vltb2,vltb3,vlptb,vlp2tb,vlpt2b,eva0,evpa1,evpa2
      	  real*8 evpa3,evta1,evta2,evta3,evpta,evp2ta,evpt2a,evb0,evpb1
      	  real*8 evpb2,evpb3,evtb1,evtb2,evtb3,evptb,evp2tb,evpt2b,dva0
      	  real*8 dvpa1,dvpa2,dvpa3,dvta1,dvta2,dvta3,dvpta,dvp2ta,dvpt2a
      	  real*8 dvb0,dvpb1,dvpb2,dvpb3,dvtb1,dvtb2,dvtb3,dvptb,dvp2tb
      	  real*8 dvpt2b,vva0,vvpa1,vvpa2,vvpa3,vvta1,vvta2,vvta3,vvpta
      	  real*8 vvp2ta,vvpt2a,vvb0,vvpb1,vvpb2,vvpb3,vvtb1,vvtb2,vvtb3
      	  real*8 vvptb,vvp2tb,vvpt2b,tl,pl,x,x2,x3,x4,dpsatt,tl2,tl3
      	  real*8 tlx,tl2x,tlx2,enwn1,enwn2,enwn3,enwn,enwd1,enwd2,enwd3
      	  real*8 enwd,enw,enl,dhwpn1,dhwpn2,dhwpn,dhwpd,dhwp,dhwtn1
      	  real*8 dhwtn2,dhwtn,dhwtd,dhwt,dhlt,dhlp,rnwn1,rnwn2,rnwn3
      	  real*8 rnwd1,rnwd2,rnwd3,rnwn,rnwd,rnw,rol,drlpn1,drlpn2
      	  real*8 drlpn,drolpd,drolp,drlen1,drlen2,drlen,droled,drolt
      	  real*8 viln1,viln2,viln3,viln,vild1,vild2,vild3,vild,vil
      	  real*8 dvlpn1,dvlpn2,dvlpn,dvilpd,dvislp,dvlen1,dvlen2,dvlen
      	  real*8 dviled,dvislt,ensn1,ensn2,ensn3,ensn,ensd1,ensd2,ensd3
      	  real*8 ensd,ens,env,dhvp1,dhvp2,dhvpn,dhvpd,dhvt1,dhvt2
      	  real*8 dhvtn,dhvtd,dhvt,dhvp,rnsn1,rnsn2,rnsn3,rnsd1,rnsd2
      	  real*8 rnsd3,rnsn,rnsd,rns,rov,drspn1,drspn2,drspn,drospd
      	  real*8 drsen1,drsen2,drsen,drostd,visn1,visn2,visn3,visn,visd1
      	  real*8 visd2,visd3,visd,vis,xvisv,dvspn1,dvspn2,dvspn,dvispd
      	  real*8 xvisl,dvisvp,dvsen1,dvsen2,dvsen,dvised,dvisvt  
            real*8 drovp,drovt,ros
      
      	  real*8 cden_correction, cden_cor
            real*8 p_energy, pcl 
            real*8 xv,xv2,xv3,tlxv,tlxv2,tl2xv
            real*8 pv
c gaz 120421 moved eval_test_h2o to comai
c            integer eval_test_h2o
c            parameter(eval_test_h2o = 1) 
            character*3 eval_dum
        if(iflg.eq.0) then
c gaz 120421 used code from tables     
c initialize and allocate memory
c gaz 110324 added n_allocate_iso
         if(.not.allocated(xv_h2o)) then
          if(ico2.lt.0) then
           n_allocate_iso = n0
           n_allocate = 1
          else
           n_allocate_iso = n0
           n_allocate = n0
          endif
          allocate(den_h2o(n_allocate_iso,6))
          allocate(enth_h2o(n_allocate,6))
          allocate(visc_h2o(n_allocate_iso,6))
          allocate(den_h2o_old(n_allocate_iso,6))
          allocate(enth_h2o_old(n_allocate,6))
          allocate(visc_h2o_old(n_allocate_iso,6))
          allocate(psat_h2o(n_allocate,4))
          allocate(humid_h2o(n_allocate,3))
          allocate(xv_h2o(n_allocate))
          allocate(phi_old(n_allocate))
          allocate(pci_old(n_allocate))         
          allocate(t_old(n_allocate))
          allocate(ieos_old(n_allocate))
          allocate(ieval_flag(n_allocate))
          ieval_flag = 0
          phi_old = 0.0d0
          pci_old = 0.0d0
          t_old = 0.0d0
          ieosd_old = 0
         endif
         pl_last  = 0.0d0
         pcl_last = 0.0d0
         tl_last  = 0.0d0  
         ieosd_last = 0
         return
        else if(iflg.eq.1) then
c             
c loop structure             
c 
        iieosl = 0
        do ii = istart,iend
         iieosd = iieos(ii)
         if(iieosd.ne.iieosl) then
c    check for change in coefficient set          
c          
c        liquid phase rational polynomial coefficients
c      
c       liquid enthalpy
c       numerator coefficients
            ela0=cel(1,iieosd)
            elpa1=cel(2,iieosd)
            elpa2=cel(3,iieosd)
            elpa3=cel(4,iieosd)
            elta1=cel(5,iieosd)
            elta2=cel(6,iieosd)
            elta3=cel(7,iieosd)
            elpta=cel(8,iieosd)
            elp2ta=cel(9,iieosd)
            elpt2a=cel(10,iieosd)
c       denomenator coefficients
            elb0=cel(11,iieosd)
            elpb1=cel(12,iieosd)
            elpb2=cel(13,iieosd)
            elpb3=cel(14,iieosd)
            eltb1=cel(15,iieosd)
            eltb2=cel(16,iieosd)
            eltb3=cel(17,iieosd)
            elptb=cel(18,iieosd)
            elp2tb=cel(19,iieosd)
            elpt2b=cel(20,iieosd)
c      
c       liquid density
c       numerator coefficients
            dla0=crl(1,iieosd)
            dlpa1=crl(2,iieosd)
            dlpa2=crl(3,iieosd)
            dlpa3=crl(4,iieosd)
            dlta1=crl(5,iieosd)
            dlta2=crl(6,iieosd)
            dlta3=crl(7,iieosd)
            dlpta=crl(8,iieosd)
            dlp2ta=crl(9,iieosd)
            dlpt2a=crl(10,iieosd)
c       denomenator coefficients
            dlb0=crl(11,iieosd)
            dlpb1=crl(12,iieosd)
            dlpb2=crl(13,iieosd)
            dlpb3=crl(14,iieosd)
            dltb1=crl(15,iieosd)
            dltb2=crl(16,iieosd)
            dltb3=crl(17,iieosd)
            dlptb=crl(18,iieosd)
            dlp2tb=crl(19,iieosd)
            dlpt2b=crl(20,iieosd)
c      
c      
c       liquid viscosity
c       numerator coefficients
            vla0=cvl(1,iieosd)
            vlpa1=cvl(2,iieosd)
            vlpa2=cvl(3,iieosd)
            vlpa3=cvl(4,iieosd)
            vlta1=cvl(5,iieosd)
            vlta2=cvl(6,iieosd)
            vlta3=cvl(7,iieosd)
            vlpta=cvl(8,iieosd)
            vlp2ta=cvl(9,iieosd)
            vlpt2a=cvl(10,iieosd)
c       denomenator coefficients
            vlb0=cvl(11,iieosd)
            vlpb1=cvl(12,iieosd)
            vlpb2=cvl(13,iieosd)
            vlpb3=cvl(14,iieosd)
            vltb1=cvl(15,iieosd)
            vltb2=cvl(16,iieosd)
            vltb3=cvl(17,iieosd)
            vlptb=cvl(18,iieosd)
            vlp2tb=cvl(19,iieosd)
            vlpt2b=cvl(20,iieosd)
c      
c      
c        vapor phase coefficients
c      
c       vapor enthalpy
c       numerator coefficients
            eva0=cev(1,iieosd)
            evpa1=cev(2,iieosd)
            evpa2=cev(3,iieosd)
            evpa3=cev(4,iieosd)
            evta1=cev(5,iieosd)
            evta2=cev(6,iieosd)
            evta3=cev(7,iieosd)
            evpta=cev(8,iieosd)
            evp2ta=cev(9,iieosd)
            evpt2a=cev(10,iieosd)
c       denomenator coefficients
            evb0=cev(11,iieosd)
            evpb1=cev(12,iieosd)
            evpb2=cev(13,iieosd)
            evpb3=cev(14,iieosd)
            evtb1=cev(15,iieosd)
            evtb2=cev(16,iieosd)
            evtb3=cev(17,iieosd)
            evptb=cev(18,iieosd)
            evp2tb=cev(19,iieosd)
            evpt2b=cev(20,iieosd)
c      
c       vapor density
c       numerator coefficients
            dva0=crv(1,iieosd)
            dvpa1=crv(2,iieosd)
            dvpa2=crv(3,iieosd)
            dvpa3=crv(4,iieosd)
            dvta1=crv(5,iieosd)
            dvta2=crv(6,iieosd)
            dvta3=crv(7,iieosd)
            dvpta=crv(8,iieosd)
            dvp2ta=crv(9,iieosd)
            dvpt2a=crv(10,iieosd)
c       denomenator coefficients
            dvb0=crv(11,iieosd)
            dvpb1=crv(12,iieosd)
            dvpb2=crv(13,iieosd)
            dvpb3=crv(14,iieosd)
            dvtb1=crv(15,iieosd)
            dvtb2=crv(16,iieosd)
            dvtb3=crv(17,iieosd)
            dvptb=crv(18,iieosd)
            dvp2tb=crv(19,iieosd)
            dvpt2b=crv(20,iieosd)
c      
c       vapor viscosity
c       numerator coefficients
c      
            vva0=cvv(1,iieosd)
            vvpa1=cvv(2,iieosd)
            vvpa2=cvv(3,iieosd)
            vvpa3=cvv(4,iieosd)
            vvta1=cvv(5,iieosd)
            vvta2=cvv(6,iieosd)
            vvta3=cvv(7,iieosd)
            vvpta=cvv(8,iieosd)
            vvp2ta=cvv(9,iieosd)
            vvpt2a=cvv(10,iieosd)
c       denomenator coefficients
            vvb0=cvv(11,iieosd)
            vvpb1=cvv(12,iieosd)
            vvpb2=cvv(13,iieosd)
            vvpb3=cvv(14,iieosd)
            vvtb1=cvv(15,iieosd)
            vvtb2=cvv(16,iieosd)
            vvtb3=cvv(17,iieosd)
            vvptb=cvv(18,iieosd)
            vvp2tb=cvv(19,iieosd)
            vvpt2b=cvv(20,iieosd)
            iieosl=iieosd     
      endif
c
c
c  check for no-change from last variable list
c 
       if(igrav.ne.0) then
          p_energy = -grav*cord(ii,igrav)
       else
          p_energy = 0.0d0
       endif 
      pl = phi(ii)
      pcl = pci(ii)
      tl = t(ii)
      ieosd = ieos(ii)
      pl_old = phi_old(ii)
      pcl_old = pci_old(ii)
      tl_old = t_old(ii)
      ieosd_old = ieos_old(ii)
      phi_old(ii) =  pl
      pci_old(ii) =  pcl
      t_old(ii) = tl
      ieos_old(ii) = ieosd
      var_last_flag = .false.
      var_old_flag = .false.
      if(abs(pl-pl_last).le.p_tol
     & .and.abs(tl-tl_last).le.t_tol
     & .and.abs(pcl-pcl_last).le.p_tol.and.
     &   ieosd.eq.ieosd_last) then
       var_last_flag = .true.
       else if(abs(pl-pl_old).le.p_tol
     & .and.abs(tl-tl_old).le.t_tol
     & .and.abs(pcl-pcl_old).le.p_tol.and.
     &   ieosd.eq.ieosd_old) then     
       var_old_flag = .true.
       endif
c gaz 120821 moved pl_last etc to property storage
c      pl_last = pl 
c      pcl_last = pcl
c      tl_last = tl
c      ieosd_last = ieosd
c gaz 101521 added pure water and heat sat pressure      
      if(ico2.eq.0) then
c water_vapor_calc ids table-lookup or polynomial evaluation            
       call water_vapor_calc(1, ii, ii, '', '')
       xv = pl
       pv =xv
      else
c water_vapor_calc_ngas ids table-lookup or polynomial evaluation  
       call water_vapor_calc_ngas(1, ii, ii, '', '')
       xv = xv_h2o(ii)  
       pv =xv
      endif       
c      
      itot_calls = itot_calls + 1
c gaz debug 120321  printout test 
c var_last_flag = .true. - last variables match  
c var_old_flag = .true. - old variables match    
c  eval_test_h2o = 0 or 2 don't evaluate properties if variable are close to previous values 
c  eval_test_h2o = 1 evaluate properties always     
      if((var_last_flag).or.(var_old_flag)) then
c don't evaluate properties          
          if(eval_test_h2o.ne.1) then          
              eval_dum ='no '
          else
              eval_dum ='yes'
          endif
      else 
c evaluate properties                
          eval_dum = 'yes'
      endif
      if(l.eq.0.and.ii.eq.1.and.eval_test_h2o.eq.2) then
       write(ierr,1500)
1500   format('l  iad ii',t12,'pl',t25,'pl_last',t40,'pl_old',t55,'tl'
     &  ,t70,'tl_last',t85,'tl_old',t100,'ieosd ieosd_last ieosd_old')
      else if(l.gt.0.and.eval_test_h2o.eq.2) then
       write(ierr,1501) l,iad,ii,pl,pl_last,pl_old,tl,tl_last,tl_old,
     &  ieosd, ieosd_last, ieosd_old, eval_dum
1501   format(1x,i5,i3,i5,1p,t15,g14.7,t29,g14.7,t43,g14.7,t57,g14.7,
     &  t71, g14.7,t85,g14.7,t100,i3,i3,i3,1x,a3)
      endif
      ieval_flag(ii) = 0
c gaz 120321 simplified
      if(eval_dum.eq.'yes') then
       ieval_flag(ii) = 1
       ic_eval = ic_eval + 1

      x=pl
      x2=x*x
      x3=x2*x                
      xv2=xv*xv
      xv3=xv2*xv                
      tl2=tl*tl
      tl3=tl2*tl
      tlx=x*tl
      tl2x=tl2*x
      tlx2=tl*x2                
      tl2=tl*tl
      tl3=tl2*tl
      tlxv=xv*tl
      tl2xv=tl2*xv
      tlxv2=tl*xv2  
      
c liquid enthalpy

      if(prop(1:3).eq.'all'.or. (prop(1:8).eq.'enthalpy'.and.
     & phase(1:6).eq.'liquid'))	 then
       enwn1=ela0+elpa1*x+elpa2*x2+elpa3*x3
       enwn2=elta1*tl+elta2*tl2+elta3*tl3
       enwn3=elpta*tlx+elpt2a*tl2x+elp2ta*tlx2
       enwn=enwn1+enwn2+enwn3
       enwd1=elb0+elpb1*x+elpb2*x2+elpb3*x3
       enwd2=eltb1*tl+eltb2*tl2+eltb3*tl3
       enwd3=elptb*tlx+elpt2b*tl2x+elp2tb*tlx2
       enwd=enwd1+enwd2+enwd3
       enw=enwn/enwd
       enl=enw + p_energy
c
c derivatives of enthalpy
c
       dhwpn1=elpa1+2*elpa2*x+3*elpa3*x2+elpta*tl
       dhwpn1=enwd*(dhwpn1+elpt2a*tl2+elp2ta*2*tlx)
       dhwpn2=elpb1+2*elpb2*x+3*elpb3*x2+elptb*tl
       dhwpn2=enwn*(dhwpn2+elpt2b*tl2+elp2tb*2*tlx)
       dhwpn=dhwpn1-dhwpn2
       dhwpd=enwd**2
       dhwp=dhwpn/dhwpd
       dhwtn1=elta1+2*elta2*tl+3*elta3*tl2+elpta*x
       dhwtn1=enwd*(dhwtn1+elpt2a*2*tlx+elp2ta*x2)
       dhwtn2=eltb1+2*eltb2*tl+3*eltb3*tl2+elptb*x
       dhwtn2=enwn*(dhwtn2+elpt2b*2*tlx+elp2tb*x2)
       dhwtn=dhwtn1-dhwtn2
       dhwtd=enwd**2
       dhwt=dhwtn/dhwtd
       dhlt=dhwt
       dhlp=dhwp
c       dhlpc=0.0
      endif
c
c liquid density
c
      if(prop(1:3).eq.'all'.or. (prop(1:7).eq.'density'.and.
     & phase(1:6).eq.'liquid'))	 then
       rnwn1=dla0+dlpa1*x+dlpa2*x2+dlpa3*x3
       rnwn2=dlta1*tl+dlta2*tl2+dlta3*tl3
       rnwn3=dlpta*tlx+dlpt2a*tl2x+dlp2ta*tlx2
       rnwn=rnwn1+rnwn2+rnwn3
       rnwd1=dlb0+dlpb1*x+dlpb2*x2+dlpb3*x3
       rnwd2=dltb1*tl+dltb2*tl2+dltb3*tl3
       rnwd3=dlptb*tlx+dlpt2b*tl2x+dlp2tb*tlx2
       rnwd=rnwd1+rnwd2+rnwd3
       rnw=rnwn/rnwd
       rol=rnw

c
c  derivatives of density
c
       drlpn1=dlpa1+2*dlpa2*x+3*dlpa3*x2+dlpta*tl
       drlpn1=rnwd*(drlpn1+2*dlp2ta*tlx+dlpt2a*tl2)
       drlpn2=dlpb1+2*dlpb2*x+3*dlpb3*x2+dlptb*tl
       drlpn2=rnwn*(drlpn2+2*dlp2tb*tlx+dlpt2b*tl2)
       drlpn=drlpn1-drlpn2
       drolpd=rnwd**2
       drolp=drlpn/drolpd
       drlen1=dlta1+2*dlta2*tl+3*dlta3*tl2+dlpta*x
       drlen1=rnwd*(drlen1+dlp2ta*x2+2*dlpt2a*tlx)
       drlen2=dltb1+2*dltb2*tl+3*dltb3*tl2+dlptb*x
       drlen2=rnwn*(drlen2+dlp2tb*x2+2*dlpt2b*tlx)
       drlen=drlen1-drlen2
       droled=rnwd**2
       drolt=drlen/droled
c       drolpc=0.0
      endif
c
c liquid viscosity
c
      if(prop(1:3).eq.'all'.or. (prop(1:9).eq.'viscosity'.and.
     & phase(1:6).eq.'liquid'))	then
       viln1=vla0+vlpa1*x+vlpa2*x2+vlpa3*x3
       viln2=vlta1*tl+vlta2*tl2+vlta3*tl3
       viln3=vlpta*tlx+vlpt2a*tl2x+vlp2ta*tlx2
       viln=viln1+viln2+viln3
       vild1=vlb0+vlpb1*x+vlpb2*x2+vlpb3*x3
       vild2=vltb1*tl+vltb2*tl2+vltb3*tl3
       vild3=vlptb*tlx+vlpt2b*tl2x+vlp2tb*tlx2
       vild=vild1+vild2+vild3
       vil=viln/vild
       xvisl=vil
c
c      derivatives of liquid viscosity
c
       dvlpn1=vlpa1+2*vlpa2*x+3*vlpa3*x2+vlpta*tl
       dvlpn1=vild*(dvlpn1+2*vlp2ta*tlx+vlpt2a*tl2)
       dvlpn2=vlpb1+2*vlpb2*x+3*vlpb3*x2+vlptb*tl
       dvlpn2=viln*(dvlpn2+2*vlp2tb*tlx+vlpt2b*tl2)
       dvlpn=dvlpn1-dvlpn2
       dvilpd=vild**2
       dvislp=dvlpn/dvilpd
       dvlen1=vlta1+2*vlta2*tl+3*vlta3*tl2+vlpta*x
       dvlen1=vild*(dvlen1+vlp2ta*x2+2*vlpt2a*tlx)
       dvlen2=vltb1+2*vltb2*tl+3*vltb3*tl2+vlptb*x
       dvlen2=viln*(dvlen2+vlp2tb*x2+2*vlpt2b*tlx)
       dvlen=dvlen1-dvlen2
       dviled=vild**2
       dvislt=dvlen/dviled
c       dvlpc=0.0
      endif
c water vapor properties      
c gaz 060421 need to define xv  
c      
c     water vapor enthalpy

      if(prop(1:3).eq.'all'.or. (prop(1:8).eq.'enthalpy'.and.
     & phase(1:5).eq.'vapor'))	then
       ensn1=eva0+evpa1*xv+evpa2*xv2+evpa3*xv3
       ensn2=evta1*tl+evta2*tl2+evta3*tl3
       ensn3=evpta*tlxv+evpt2a*tl2xv+evp2ta*tlxv2
       ensn=ensn1+ensn2+ensn3
       ensd1=evb0+evpb1*xv+evpb2*xv2+evpb3*xv3
       ensd2=evtb1*tl+evtb2*tl2+evtb3*tl3
       ensd3=evptb*tlxv+evpt2b*tl2xv+evp2tb*tlxv2
       ensd=ensd1+ensd2+ensd3
       ens=ensn/ensd
       env=ens + p_energy
c
c        derivatives of water vapor enthalpy
c
       dhvp1=evpa1+2*evpa2*xv+3*evpa3*xv2+evpta*tl
       dhvp1=ensd*(dhvp1+2*evp2ta*tlxv+evpt2a*tl2)
       dhvp2=evpb1+2*evpb2*xv+3*evpb3*xv2+evptb*tl
       dhvp2=ensn*(dhvp2+2*evp2tb*tlxv+evpt2b*tl2)
       dhvpn=dhvp1-dhvp2
       dhvpd=ensd**2
       dhvp=dhvpn/dhvpd
       dhvt1=evta1+2*evta2*tl+3*evta3*tl2+evpta*xv
       dhvt1=ensd*(dhvt1+evp2ta*xv2+2*evpt2a*tlxv)
       dhvt2=evtb1+2*evtb2*tl+3*evtb3*tl2+evptb*xv
       dhvt2=ensn*(dhvt2+evp2tb*xv2+2*evpt2b*tlxv)
       dhvtn=dhvt1-dhvt2
       dhvtd=ensd**2
       dhvt=dhvtn/dhvtd
c       dhvpc=-dhvp
      endif
c
c water vapor density
c
      if(prop(1:3).eq.'all'.or. (prop(1:7).eq.'density'.and.
     & phase(1:5).eq.'vapor'))	then
       rnsn1=dva0+dvpa1*xv+dvpa2*xv2+dvpa3*xv3
       rnsn2=dvta1*tl+dvta2*tl2+dvta3*tl3
       rnsn3=dvpta*tlxv+dvpt2a*tl2xv+dvp2ta*tlxv2
       rnsn=rnsn1+rnsn2+rnsn3
       rnsd1=dvb0+dvpb1*xv+dvpb2*xv2+dvpb3*xv3
       rnsd2=dvtb1*tl+dvtb2*tl2+dvtb3*tl3
       rnsd3=dvptb*tlxv+dvpt2b*tl2xv+dvp2tb*tlxv2
       rnsd=rnsd1+rnsd2+rnsd3
       rns=rnsn/rnsd
       ros=rns
       rov=ros
c
c        water derivatives of vapor density
c
       drspn1=dvpa1+2*dvpa2*xv+3*dvpa3*xv2+dvpta*tl
       drspn1=rnsd*(drspn1+2*dvp2ta*tlxv+dvpt2a*tl2)
       drspn2=dvpb1+2*dvpb2*xv+3*dvpb3*xv2+dvptb*tl
       drspn2=rnsn*(drspn2+2*dvp2tb*tlxv+dvpt2b*tl2)
       drspn=drspn1-drspn2
       drospd=rnsd**2
       drovp=drspn/drospd
       drsen1=dvta1+2*dvta2*tl+3*dvta3*tl2+dvpta*xv
       drsen1=rnsd*(drsen1+dvp2ta*xv2+2*dvpt2a*tlxv)
       drsen2=dvtb1+2*dvtb2*tl+3*dvtb3*tl2+dvptb*xv
       drsen2=rnsn*(drsen2+dvp2tb*xv2+2*dvpt2b*tlxv)
       drsen=drsen1-drsen2
       drostd=rnsd**2
       drovt=drsen/drostd
c       drovpc=-drovp
      endif
c
c water vapor viscosity
c
      if(prop(1:3).eq.'all'.or. (prop(1:7).eq.'viscosity'.and.
     & phase(1:5).eq.'vapor'))	then
       visn1=vva0+vvpa1*xv+vvpa2*xv2+vvpa3*xv3
       visn2=vvta1*tl+vvta2*tl2+vvta3*tl3
       visn3=vvpta*tlxv+vvpt2a*tl2xv+vvp2ta*tlxv2
       visn=visn1+visn2+visn3
       visd1=vvb0+vvpb1*xv+vvpb2*xv2+vvpb3*xv3
       visd2=vvtb1*tl+vvtb2*tl2+vvtb3*tl3
       visd3=vvptb*tlxv+vvpt2b*tl2xv+vvp2tb*tlxv2
       visd=visd1+visd2+visd3
       vis=visn/visd
       xvisv=vis
c
c      derivatives of water vapor viscosity
c
       dvspn1=vvpa1+2*vvpa2*xv+3*vvpa3*xv2+vvpta*tl
       dvspn1=visd*(dvspn1+2*vvp2ta*tlxv+vvpt2a*tl2)
       dvspn2=vvpb1+2*vvpb2*xv+3*vvpb3*xv2+vvptb*tl
       dvspn2=visn*(dvspn2+2*vvp2tb*tlxv+vvpt2b*tl2)
       dvspn=dvspn1-dvspn2
       dvispd=visd**2
       dvisvp=dvspn/dvispd
       dvsen1=vvta1+2*vvta2*tl+3*vvta3*tl2+vvpta*xv
       dvsen1=visd*(dvsen1+vvp2ta*xv2+2*vvpt2a*tlxv)
       dvsen2=vvtb1+2*vvtb2*tl+3*vvtb3*tl2+vvptb*xv
       dvsen2=visn*(dvsen2+vvp2tb*xv2+2*vvpt2b*tlxv)
       dvsen=dvsen1-dvsen2
       dvised=visd**2
       dvisvt=dvsen/dvised      
c       dvvpc=-dvisvp 
      endif
c      
      den_h2o(ii,1)= rol
      den_h2o(ii,2)= drolp
      den_h2o(ii,3)= drolt
      den_h2o(ii,4)= rov
      den_h2o(ii,5)= drovp
      den_h2o(ii,6)= drovt
c      
      enth_h2o(ii,1)= enl
      enth_h2o(ii,2)= dhlp
      enth_h2o(ii,3)= dhlt
      enth_h2o(ii,4)= env
      enth_h2o(ii,5)= dhvp
      enth_h2o(ii,6)= dhvt
c
      visc_h2o(ii,1)= xvisl
      visc_h2o(ii,2)= dvislp
      visc_h2o(ii,3)= dvislt
      visc_h2o(ii,4)= xvisv
      visc_h2o(ii,5)= dvisvp
      visc_h2o(ii,6)= dvisvt          

      den_h2o_old(ii,1) = rol
      den_h2o_old(ii,2) = drolp
      den_h2o_old(ii,3) = drolt
      den_h2o_old(ii,4) = rov
      den_h2o_old(ii,5) = drovp
      den_h2o_old(ii,6) = drovt
      enth_h2o_old(ii,1)= enl
      enth_h2o_old(ii,2)= dhlp
      enth_h2o_old(ii,3)= dhlt
      enth_h2o_old(ii,4)= env
      enth_h2o_old(ii,5)= dhvp
      enth_h2o_old(ii,6)= dhvt
      visc_h2o_old(ii,1)= xvisl
      visc_h2o_old(ii,2)= dvislp
      visc_h2o_old(ii,3)= dvislt
      visc_h2o_old(ii,4)= xvisv
      visc_h2o_old(ii,5)= dvisvp
      visc_h2o_old(ii,6)= dvisvt 

c gaz make sure variable match properties (in time)  
         pl_last = pl 
         pcl_last = pcl
         tl_last = tl
         ieosd_last = ieosd      
      
      	rol_last	  =     rol	   
      	drolp_last    =     drolp	   
      	drolt_last    =     drolt	   
      	rov_last	  =     rov	   
      	drovp_last    =     drovp	   
      	drovt_last    =     drovt	   
      	enl_last	  =     enl	   
      	dhlp_last     =     dhlp	   
      	dhlt_last     =     dhlt	   
      	env_last	  =     env	   
      	dhvp_last     =     dhvp	   
      	dhvt_last     =     dhvt	   
      	xvisl_last    =     xvisl	   
      	dvislp_last   =     dvislp	   
      	dvislt_last   =     dvislt	   
      	xvisv_last    =     xvisv	   
      	dvisvp_last   =     dvisvp	   
      	dvisvt_last   =     dvisvt 	   

      elseif(var_old_flag) then  
c use old unchanged values  if p,t(node) = p,t(node-1) 
c last properties rol-dvisvt stored on	com_prop_data
c or if p,t(node) = p,t(node,last iter)
      
       den_h2o(ii,1) =den_h2o_old(ii,1) 
       den_h2o(ii,2) =den_h2o_old(ii,2) 
       den_h2o(ii,3) =den_h2o_old(ii,3) 
       den_h2o(ii,4) =den_h2o_old(ii,4) 
       den_h2o(ii,5) =den_h2o_old(ii,5) 
       den_h2o(ii,6) =den_h2o_old(ii,6) 
       enth_h2o(ii,1)=enth_h2o_old(ii,1)
       enth_h2o(ii,2)=enth_h2o_old(ii,2)
       enth_h2o(ii,3)=enth_h2o_old(ii,3)
       enth_h2o(ii,4)=enth_h2o_old(ii,4)
       enth_h2o(ii,5)=enth_h2o_old(ii,5)
       enth_h2o(ii,6)=enth_h2o_old(ii,6)
       visc_h2o(ii,1)=visc_h2o_old(ii,1)
       visc_h2o(ii,2)=visc_h2o_old(ii,2)
       visc_h2o(ii,3)=visc_h2o_old(ii,3)
       visc_h2o(ii,4)=visc_h2o_old(ii,4)
       visc_h2o(ii,5)=visc_h2o_old(ii,5)
       visc_h2o(ii,6)=visc_h2o_old(ii,6) 
       
      elseif(var_last_flag) then
          
       den_h2o(ii,1) =rol_last
       den_h2o(ii,2) =drolp_last
       den_h2o(ii,3) =drolt_last
       den_h2o(ii,4) =rov_last
       den_h2o(ii,5) =drovp_last
       den_h2o(ii,6) =drovt_last
       enth_h2o(ii,1)=enl_last
       enth_h2o(ii,2)=dhlp_last
       enth_h2o(ii,3)=dhlt_last
       enth_h2o(ii,4)=env_last
       enth_h2o(ii,5)=dhvp_last
       enth_h2o(ii,6)=dhvt_last
       visc_h2o(ii,1)=xvisl_last
       visc_h2o(ii,2)=dvislp_last
       visc_h2o(ii,3)=dvislt_last
       visc_h2o(ii,4)=xvisv_last
       visc_h2o(ii,5)=dvisvp_last
       visc_h2o(ii,6)=dvisvt_last       
      endif 
c
      enddo
      else if(iflg.eq.-1) then
c start counting and other metrics  
          ic_eval = 0
          itot_calls = 0
      else if(iflg.eq.-2) then
c end counting and other metrics 
       if(ic_eval.gt.1) then
         ratio_eval_tot =  ic_eval/(itot_calls + 1.e-9)  
        if(iout.ne.0) write(iout,100)ic_eval,itot_calls,ratio_eval_tot
        if(iptty.ne.0) write(iptty,100)
     &                ic_eval,itot_calls,ratio_eval_tot 
c gaz 032922 increase the integer format i9 to i12            
100    format(/,' number of unique property calcs',1x,i12,1x,/,
     &    ' total number of prossible calls',1x,i12,1x,/,
     &    ' fraction of calcs to total calls',1x,f8.3) 
      endif
 
      endif
      return
      end
      subroutine h2o_props_table(iflg, istart, iend, prop, phase)      
c      
c subroutine to evaluate water properties using tables
c gaz 060621 initial coding
c      
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, iout, p_tol, t_tol,
     &     pc_tol,roc0, iptty, iad, l, eval_test_h2o
          use combi, only: cord
          use comci, only: cden
          use comdi, only: phi, t, ieos
          use comdti, only: n0
          use comii, only: cel, cev, crl, crv, cvl, cvv
          use comfi, only: pci

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate, n_allocate_iso
          character*9 prop, phase

          integer mid,mi,ieosd,kq
c gaz 102621
          integer ieosd_sc
      	  real*8 pl,tl,psatl,dtsatp,dpsats,xv
            real*8 vis,xvisv,xvisl
            real*8 xrl,xrv
            real*8 rov,rol,ros
            real*8 enw,enl,env,ens
            real*8 dvisvt,dvisvp,dvislt,dvislp
            real*8 drovt,drovp,drovpc,dvvpc
            real*8 drvp,drv
            real*8 drolt,drolp
            real*8 drlp,drl
            real*8 dhvt,dhvp,dhvpc,dhlt,dhlp      
c      	  real*8 cden_correction, cden_cor
            real*8 p_energy, pcl 

c gaz 070321            
            real*8 dum1,dumb,dumc,value(9), value_a(9)
            real*8 pv, dhlpc,drolpc,dvlpc
            integer istate, ifail, isol_print   
c gaz 120421 moved eval_test_h2o to comai
c            integer eval_test_h2o
c            parameter(eval_test_h2o = 1) 
            character*3 eval_dum
        if(iflg.eq.0) then
c initialize and allocate memory
c gaz 110324
         if(.not.allocated(xv_h2o)) then
          if(ico2.lt.0) then
           n_allocate_iso = n0
           n_allocate = 1
          else
           n_allocate_iso = n0
           n_allocate = n0
          endif
          allocate(den_h2o(n_allocate_iso,6))
          allocate(enth_h2o(n_allocate,6))
          allocate(visc_h2o(n_allocate_iso,6))
          allocate(den_h2o_old(n_allocate_iso,6))
          allocate(enth_h2o_old(n_allocate,6))
          allocate(visc_h2o_old(n_allocate_iso,6))
          allocate(psat_h2o(n_allocate,4))
          allocate(humid_h2o(n_allocate,3))
          allocate(xv_h2o(n_allocate))
          allocate(phi_old(n_allocate))
          allocate(pci_old(n_allocate))         
          allocate(t_old(n_allocate))
          allocate(ieos_old(n_allocate))
          allocate(ieval_flag(n_allocate))
          ieval_flag = 0
          phi_old = 0.0d0
          pci_old = 0.0d0
          t_old = 0.0d0
          ieosd_old = 0
         endif
         pl_last  = 0.0d0
         pcl_last = 0.0d0
         tl_last  = 0.0d0  
         ieosd_last = 0
         return
        else if(iflg.eq.1) then
c             
c loop structure             
c 

        do ii = istart,iend
        
c          
c        liquid phase table data
c      

       if(igrav.ne.0) then
          p_energy = -grav*cord(ii,igrav)
       else
          p_energy = 0.0d0
       endif 
      pl = phi(ii)
      pcl = pci(ii)
      tl = t(ii)
      ieosd = ieos(ii)
      pl_old = phi_old(ii)
      pcl_old = pci_old(ii)
      tl_old = t_old(ii)
      ieosd_old = ieos_old(ii)
      phi_old(ii) =  pl
      pci_old(ii) =  pcl
      t_old(ii) = tl
      ieos_old(ii) = ieosd
      var_last_flag = .false.
      var_old_flag = .false.
      if(abs(pl-pl_last).le.p_tol
     & .and.abs(tl-tl_last).le.t_tol
     & .and.abs(pcl-pcl_last).le.p_tol.and.
     &   ieosd.eq.ieosd_last) then
       var_last_flag = .true.
       else if(abs(pl-pl_old).le.p_tol
     & .and.abs(tl-tl_old).le.t_tol
     & .and.abs(pcl-pcl_old).le.p_tol.and.
     &   ieosd.eq.ieosd_old) then     
       var_old_flag = .true.
      endif
c gaz 120821 moved to property assignment
c      pl_last = pl 
c      pcl_last = pcl
c      tl_last = tl
c      ieosd_last = ieosd
c gaz 101521 added pure water and heat sat pressure      
      if(ico2.eq.0) then
c water_vapor_calc ids table-lookup or polynomial evaluation            
       call water_vapor_calc(1, ii, ii, '', '')
       xv = pl
       pv =xv
      else
c water_vapor_calc_ngas ids table-lookup or polynomial evaluation  
       call water_vapor_calc_ngas(1, ii, ii, '', '')
       xv = xv_h2o(ii)  
       pv =xv
      endif       
c      
      itot_calls = itot_calls + 1
c gaz debug 120321  printout test 
c var_last_flag = .true. - last variables match  
c var_old_flag = .true. - old variables match    
c  eval_test_h2o = 0 or 2 don't evaluate properties if variable are close to previous values 
c  eval_test_h2o = 1 evaluate properties always       
      if((var_last_flag).or.(var_old_flag)) then
c don't evaluate properties          
          if(eval_test_h2o.ne.1) then
              eval_dum ='no '
          else
              eval_dum ='yes'
          endif
      else 
c evaluate properties                
          eval_dum = 'yes'
      endif
      if(l.eq.0.and.ii.eq.1.and.eval_test_h2o.eq.2) then
       write(ierr,1500)
1500   format('l  iad ii',t12,'pl',t25,'pl_last',t40,'pl_old',t55,'tl'
     &  ,t70,'tl_last',t85,'tl_old',t100,'ieosd ieosd_last ieosd_old')
      else if(l.gt.0.and.eval_test_h2o.eq.2) then
       write(ierr,1501) l,iad,ii,pl,pl_last,pl_old,tl,tl_last,tl_old,
     &  ieosd, ieosd_last, ieosd_old, eval_dum
1501   format(1x,i5,i3,i5,1p,t15,g14.7,t29,g14.7,t43,g14.7,t57,g14.7,
     & t71,g14.7,t85,g14.7,t100,i3,i3,i3,1x,a3)
      endif
      ieval_flag(ii) = 0
c gaz 120321 simplified
      if(eval_dum.eq.'yes') then
       ieval_flag(ii) = 1
       ic_eval = ic_eval + 1
c
c gaz 102621 grouping of SC phase is different for pure water (ico2=0) and ngas (ico2>0)
c               
             ieosd = ieos(ii) 
             if(ieosd.eq.4.and.ico2.eq.0) then
              ieosd_sc = 1
             else if(ieosd.eq.4.and.ico2.gt.0) then
              ieosd_sc = 3
             else
              ieosd_sc =ieosd
             endif
              value = 0.0

               if(ieosd_sc.eq.1)then
                  call h2o_properties_new(4,ieosd,pl,tl,dum1,istate,
     &                 dumb,value,dumc)
			      
               elseif (ieosd_sc.eq.2) then                  
c call with ieosd = 1 and use chain rule 
                if(ico2.eq.0) then
                  call h2o_properties_new(5,ieosd,pl,tl,dum1,istate,
     &                 dumb,value,dumc)                     
c don't need chain rule for pure water                  
c                drolt = 0.0
c                dhlt = 0.0
c                dvislt = 0.0  
                 value(2) = 0.0
                 value(5) = 0.0
                 value(8) = 0.0                     
                else                    
                  call h2o_properties_new(4,1,pl,tl,dum1,istate,
     &                 dumb,value,dumc)
                endif
               endif     
                   rol = value(1)
                   drolt = value(2)
                   drolp = value(3)
                   enl = value(4) + p_energy
c gaz 120518 enw has slight (grav) correction                   
                   enw = enl
                   dhlt = value(5)
                   dhlp = value(6)
                   xvisl = value(7)
                   dvislt = value(8)
                   dvislp = value(9)  
                   
                   dhlpc=0.0
                   drolpc=0.0 
                   dvlpc=0.0
c chain rule not needed for table (d/dp already includes)

c-----------------------------------------------------------------------       
     
c water vapor properties      

                 value = 0.0
c gaz 072420 add ieosd.eq.4 for SC                
                if(ieosd_sc.eq.3)then
                  call h2o_properties_new(4,ieosd,pv,tl,dum1,istate,
     &                 dumb,value,dumc)
c gaz 072520 force ieosd to stay = 3, changes from 4 to 1 in h2o_properties_new
                  ieosd = 3
                  xrv = 1
                  xrl = 1
                  drl = 0
                  drv = 0
                  drlp = 0
                  drvp = 0
                elseif (ieosd_sc.eq.2) then
c gaz 120718  derivatives match rat polynomials  
c gaz pure water uses slightly different approach   
                 if(ico2.eq.0) then
                  call h2o_properties_new(6,ieosd,pl,tl,dum1,istate,
     &                 dumb,value,dumc)  
c don't need chain rule for pure water                  
c                drovt = 0.0
c                dhvt = 0.0
c                dvisvt = 0.0  
                 value(2) = 0.0
                 value(5) = 0.0
                 value(8) = 0.0
                 else 
                  call h2o_properties_new(4,3,pv,tl,dum1,istate,
     &                 dumb,value,dumc)
                 endif
                endif    
                   rov = value(1)               
                   ros = rov
                   drovt = value(2)
                   drovp = value(3)
                   env = value(4) + p_energy                 
                   ens = env
                   dhvt = value(5)
                   dhvp = value(6)
                   xvisv = value(7)
                   vis = xvisv
                   dvisvt = value(8)
                   dvisvp = value(9)  
                   dhvpc=-dhvp
                   drovpc=-drovp
                   dvvpc=-dvisvp     
     
c      
      den_h2o(ii,1)= rol
      den_h2o(ii,2)= drolp
      den_h2o(ii,3)= drolt
      den_h2o(ii,4)= rov
      den_h2o(ii,5)= drovp
      den_h2o(ii,6)= drovt
c      
      enth_h2o(ii,1)= enl
      enth_h2o(ii,2)= dhlp
      enth_h2o(ii,3)= dhlt
      enth_h2o(ii,4)= env
      enth_h2o(ii,5)= dhvp
      enth_h2o(ii,6)= dhvt
c
      visc_h2o(ii,1)= xvisl
      visc_h2o(ii,2)= dvislp
      visc_h2o(ii,3)= dvislt
      visc_h2o(ii,4)= xvisv
      visc_h2o(ii,5)= dvisvp
      visc_h2o(ii,6)= dvisvt          
c save values for node ii for comparison with next update
      den_h2o_old(ii,1) = rol
      den_h2o_old(ii,2) = drolp
      den_h2o_old(ii,3) = drolt
      den_h2o_old(ii,4) = rov
      den_h2o_old(ii,5) = drovp
      den_h2o_old(ii,6) = drovt
      enth_h2o_old(ii,1)= enl
      enth_h2o_old(ii,2)= dhlp
      enth_h2o_old(ii,3)= dhlt
      enth_h2o_old(ii,4)= env
      enth_h2o_old(ii,5)= dhvp
      enth_h2o_old(ii,6)= dhvt
      visc_h2o_old(ii,1)= xvisl
      visc_h2o_old(ii,2)= dvislp
      visc_h2o_old(ii,3)= dvislt
      visc_h2o_old(ii,4)= xvisv
      visc_h2o_old(ii,5)= dvisvp
      visc_h2o_old(ii,6)= dvisvt 
c save values for node ii for comparison with ii+1
c if abs(tl-tl_tast) < tol_t and    abs(pl-pl_tast) < tol_p use previous
c densities,enthalpies,viscosities   
c gaz 022020 make sure variable match properties (in time)  
           pl_last = pl 
           pcl_last = pcl
           tl_last = tl
           ieosd_last = ieosd    
           
      	rol_last	  =     rol	   
      	drolp_last    =     drolp	   
      	drolt_last    =     drolt	   
      	rov_last	  =     rov	   
      	drovp_last    =     drovp	   
      	drovt_last    =     drovt	   
      	enl_last	  =     enl	   
      	dhlp_last     =     dhlp	   
      	dhlt_last     =     dhlt	   
      	env_last	  =     env	   
      	dhvp_last     =     dhvp	   
      	dhvt_last     =     dhvt	   
      	xvisl_last    =     xvisl	   
      	dvislp_last   =     dvislp	   
      	dvislt_last   =     dvislt	   
      	xvisv_last    =     xvisv	   
      	dvisvp_last   =     dvisvp	   
      	dvisvt_last   =     dvisvt 	   

      elseif(var_old_flag) then  
c use old unchanged values  if p,t(node) = p,t(node-1) 
c last properties rol-dvisvt stored on com_prop_data
c or if p,t(node) = p,t(node,last iter)
      
       den_h2o(ii,1) =den_h2o_old(ii,1) 
       den_h2o(ii,2) =den_h2o_old(ii,2) 
       den_h2o(ii,3) =den_h2o_old(ii,3) 
       den_h2o(ii,4) =den_h2o_old(ii,4) 
       den_h2o(ii,5) =den_h2o_old(ii,5) 
       den_h2o(ii,6) =den_h2o_old(ii,6) 
       enth_h2o(ii,1)=enth_h2o_old(ii,1)
       enth_h2o(ii,2)=enth_h2o_old(ii,2)
       enth_h2o(ii,3)=enth_h2o_old(ii,3)
       enth_h2o(ii,4)=enth_h2o_old(ii,4)
       enth_h2o(ii,5)=enth_h2o_old(ii,5)
       enth_h2o(ii,6)=enth_h2o_old(ii,6)
       visc_h2o(ii,1)=visc_h2o_old(ii,1)
       visc_h2o(ii,2)=visc_h2o_old(ii,2)
       visc_h2o(ii,3)=visc_h2o_old(ii,3)
       visc_h2o(ii,4)=visc_h2o_old(ii,4)
       visc_h2o(ii,5)=visc_h2o_old(ii,5)
       visc_h2o(ii,6)=visc_h2o_old(ii,6) 
       
      elseif(var_last_flag) then
          
       den_h2o(ii,1) =rol_last
       den_h2o(ii,2) =drolp_last
       den_h2o(ii,3) =drolt_last
       den_h2o(ii,4) =rov_last
       den_h2o(ii,5) =drovp_last
       den_h2o(ii,6) =drovt_last
       enth_h2o(ii,1)=enl_last
       enth_h2o(ii,2)=dhlp_last
       enth_h2o(ii,3)=dhlt_last
       enth_h2o(ii,4)=env_last
       enth_h2o(ii,5)=dhvp_last
       enth_h2o(ii,6)=dhvt_last
       visc_h2o(ii,1)=xvisl_last
       visc_h2o(ii,2)=dvislp_last
       visc_h2o(ii,3)=dvislt_last
       visc_h2o(ii,4)=xvisv_last
       visc_h2o(ii,5)=dvisvp_last
       visc_h2o(ii,6)=dvisvt_last       
      endif 
c      
      enddo
      else if(iflg.eq.-1) then
c start counting and other metrics  
          ic_eval = 0
          itot_calls = 0
      else if(iflg.eq.-2) then
c end counting and other metrics 
       ratio_eval_tot =  ic_eval/(itot_calls + 1.e-9)  
c       write(ierr,100) ic_eval, itot_calls, ratio_eval_tot
      if(iout.ne.0) write(iout,100) ic_eval,itot_calls,ratio_eval_tot
      if(iptty.ne.0) write(iptty,100) ic_eval,itot_calls,ratio_eval_tot
c gaz 032922 increase the integer format i9 to i12          
100    format(/,' number of unique property calcs',1x,i12,1x,/,
     &    ' total number of prossible calls',1x,i12,1x,/,
     &    ' fraction of calcs to total calls',1x,f8.3) 
c      stop 
      endif
      return
      end            
      
      subroutine water_vapor_calc_ngas(iflg, istart, iend, prop, phase)
c      
c subroutine to evaluate water vapor properties and saturation pressure
c gaz 060621 initial coding
c for ngas      
c
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, pcrit_h2o,
     &     tcrit_h2o, l
          use combi, only: cord
          use comci, only: cden, dpcef
          use comdi, only: phi, t, iieos, an, pcp, ieos
          use comdti, only: n0
          use comii, only: cel, cev, crl, crv, cvl, cvv
          use comfi, only: pci
          use comgi, only: dfpcp

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate, n_allocate_iso
          integer ipv_tol
          character*9 prop, phase

          integer ndummy,iieosl,mid,mi,ieosd,iieosd,kq
          real*8 dpsatt, dpct, psatl_100, dpsats
          real*8 pv_tol, pa_tol, pv, tv, xv, pcl, dtdp
          real*8 psatl, dpsatt_100, pl, tl
          parameter(pv_tol=1.d-9, pa_tol = 1.d-16)  
c gaz 070621 debug  
c note that psatl evaluates either table  or polynomial as needed       
        if(iflg.eq.0) then
c allocate space and initialize       

        else if(iflg.eq.1) then   
c property evaluation         
         do ii = istart, iend          
c
c two phase conditions
c 
          pv = 0.0   
          dpsatt=0.0
          dpct = 0.0
          psatl_100 = 0.0
          if(ieos(ii).eq.2) then
c gaz 101521, gaz 103021 (back again)
            pv=psatl(t(ii),pcp(ii),dpcef(ii),dpsatt,dpsats,
     &                   0,an(ii))
c gaz 103021              
c            tv=psatl(phi(ii),pcp(ii),dpcef(ii),dpsatt,dpsats,
c     &                   1,an(ii))            
            pl = phi(ii)
            pcl= pl-pv
            pcl=max(pa_tol,pcl)
            pcl=min(pl,pcl)
            pci(ii)=pcl            
            dpct=-dpsatt
            dtdp = 1./dpsatt
            psat_h2o(ii,1) = pv
            psat_h2o(ii,2) = dtdp
            psat_h2o(ii,3) = dpct
            psat_h2o(ii,4) = dpsats
          else if(ieos(ii).eq.3) then
c calculate water vapor pressure at 100 % humidity
c needed for rel perm calc
c gaz 070221
          pl = phi(ii)
          pcl = pci(ii)           
           if(t(ii).ge.tcrit_h2o) then
              psatl_100 = pcrit_h2o
           else
            psatl_100 = psatl(t(ii),pcp(ii),dpcef(ii),
     &      dpsatt_100,dpsats,0,an(ii))
           endif
            psat_h2o(ii,1) = psatl_100             
            psat_h2o(ii,2) = 0.0
            psat_h2o(ii,3) = 0.0
            psat_h2o(ii,4) = 0.0
          endif
c 
c calculate vapor pressure and xv (for rat polys)
c
           pv = pl-pcl
           if(pv.lt.pv_tol) then
            ipv_tol = 0
            xv= pv_tol
           else
            xv= pv
            ipv_tol = 0
           endif 
           xv_h2o(ii) = xv 
c
c  humidity calculations         
c
         if(ieos(ii).eq.3) then          
          humid_h2o(ii,1) = xv/max(psatl_100,1.d-15)
          humid_h2o(ii,2) = 1./max(psatl_100,1.d-15)
c          dhumidpc = -humid_h2o(ii,2)
          humid_h2o(ii,3) = (-xv/max(psatl_100,1.d-15)**2)*dpsatt_100
         else 
          humid_h2o(ii,1) = 1.  
          humid_h2o(ii,2) = 0.0
          humid_h2o(ii,3) = 0.0  
c         dhumidpc = 0.0
         endif               
       enddo
      endif 
      return 
      end

      subroutine water_vapor_calc(iflg, istart, iend, prop, phase)
c      
c subroutine to evaluate water vapor properties and saturation pressure
c gaz 060621 initial coding
c for pure water and heat      
c
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, pcrit_h2o,
     &     tcrit_h2o, l
          use combi, only: cord
          use comci, only: cden, dpcef
          use comdi, only: phi, t, iieos, an, pcp, ieos
          use comdti, only: n0
          use comii, only: cel, cev, crl, crv, cvl, cvv
          use comfi, only: pci
          use comgi, only: dfpcp

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate, n_allocate_iso
          integer ipv_tol
          character*9 prop, phase

          integer ndummy,iieosl,mid,mi,ieosd,iieosd,kq
          real*8 dpsatt, dpct, psatl_100, dpsats
          real*8 pv_tol, pa_tol, pv, tv, xv, pcl, dtdp
          real*8 psatl, dpsatt_100, pl, tl
          parameter(pv_tol=1.d-9, pa_tol = 1.d-16)  
c gaz 070621 debug  
c note that psatl evaluates either table  or polynomial as needed       
        if(iflg.eq.0) then
c allocate space and initialize       

        else if(iflg.eq.1) then   
c property evaluation         
         do ii = istart, iend          
c
c two phase conditions
c 
          pl = phi(ii)
          pv = pl   
          dpsatt=0.0
          dpct = 0.0
          psatl_100 = 0.0
          if(ieos(ii).eq.2) then
c gaz 101521
c     calculate temperature and dt/dp
c matches coding in sub thermw.f 
           tl=psatl(pl,pcp(ii),dpcef(ii),
     2           dtdp,dpsats,1,an(ii))

            psat_h2o(ii,1) = pl 
            psat_h2o(ii,2) = dtdp
            psat_h2o(ii,3) = 0.0
            psat_h2o(ii,4) = 0.0

          else if(ieos(ii).eq.3) then
c calculate water vapor pressure at 100 % huiidity
c needed for rel perm calc
c gaz 070221
          pl = phi(ii)
          pcl = pci(ii)           
           if(t(ii).ge.tcrit_h2o) then
              psatl_100 = pcrit_h2o
           else
            psatl_100 = psatl(t(ii),pcp(ii),dpcef(ii),
     &      dpsatt_100,dpsats,0,an(ii))
           endif
            psat_h2o(ii,1) = psatl_100             
            psat_h2o(ii,2) = 0.0
            psat_h2o(ii,3) = 0.0
            psat_h2o(ii,4) = 0.0
          endif
c 
c calculate vapor pressure and xv (for rat polys)
c
           pv = pl
           if(pv.lt.pv_tol) then
            ipv_tol = 0
            xv= pv_tol
           else
            xv= pv
            ipv_tol = 0
           endif 
           xv_h2o(ii) = xv 
       enddo
      endif 
      return 
      end

      subroutine air_props(iflg, istart, iend, prop, phase)
c      
c subroutine to evaluate air properties using tables
c gaz 070921 initial coding
c air property calls assume water property calls are made earlier
c call optimization probably workd only for nodal sequence (all water and air props for m=1,2,..)  
c note that the derivatives wrt t and f(t,p) dfdt, dfdp (value_a(1:3) wheeas the storage for
c den_ngas(:,:) is f(t,p) dfdp, dfdt   
c      
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, iout, p_tol, t_tol,
     &     pc_tol,roc0, iptty
          use combi, only: cord
          use comdi, only: phi, t, ieos
          use comdti, only: n0
          use comfi, only: pci
          use comii, only: pmin_air_tabl

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate, n_allocate_iso
          character*9 prop, phase

          integer mid,mi,ieosd,kq

      	  real*8 pl,tl,psatl,dtsatp,dpsats,xv
            real*8 vis,xvisv,xvisl
            real*8 xrl,xrv
            real*8 rov,rol,ros
            real*8 enw,enl,env,ens
            real*8 dvisvt,dvisvp,dvislt,dvislp
            real*8 drovt,drovp,drovpc,dvvpc
            real*8 drvp,drv
            real*8 drolt,drolp
            real*8 drlp,drl
            real*8 dhvt,dhvp,dhvpc,dhlt,dhlp   
            real*8 pv, dhlpc,drolpc,dvlpc
            real*8 p_energy, pcl, pcl_in
            
            real*8 hsol,dhsolp,dhsolt
            real*8 hcl
c gaz 052522 define local variables
c iflg = 0, allocate memory
c iflg = 1, loop for property evaluation
c iflg = -1, initialize metrics  
c iflg = -2, finish metrics   
c istart = node number to start property calculation
c iend = node number to finish  property calculation 
c prop = property to calculate (den,vis,enth,all) for ngas always 'all'
c phase =  phase, ngas always is 'gas'         
c gaz 070321            
            real*8 dum1,dumb,dumc,value(9),value_a(9)
            integer istate, ifail, isol_print    
c gaz 082821
            real*8 pcl0, dxvsct, dxvscpc, dxvscp, cpa, dcpat
            real*8 hcg_zero
            parameter (hcg_zero = 0.3992908d0)
c gaz 100621
            integer eval_test_air
            parameter(eval_test_air = 1)            
        if(iflg.eq.0) then
c initialize and allocate memory
         if(.not.allocated(den_ngas)) then
c gaz 110324
          if(ico2.lt.0) then
           n_allocate_iso = n0
           n_allocate = 1
          else
           n_allocate_iso = n0
           n_allocate = n0
          endif
          allocate(den_ngas(n_allocate_iso,6))
          allocate(enth_ngas(n_allocate,6))
          allocate(visc_ngas(n_allocate_iso,6))
          allocate(xnl_ngas(n_allocate,6))
          allocate(den_ngas_old(n_allocate_iso,6))
          allocate(enth_ngas_old(n_allocate,6))
          allocate(visc_ngas_old(n_allocate_iso,6))
          allocate(xnl_ngas_old(n_allocate,6))
         endif
         if(.not.allocated(phi_old)) then
          allocate(phi_old(n_allocate))
          allocate(t_old(n_allocate))
          allocate(pci_old(n_allocate))
          allocate(ieval_flag(n_allocate))
          ieval_flag = 0
          phi_old = 0.0d0
          t_old = 0.0d0
          pci_old = 0.0d0
         endif 
         return
        else if(iflg.eq.1) then
c             
c loop structure             
c 

        do ii = istart,iend
        
c         
c        ngas table data
c      

       if(igrav.ne.0) then
          p_energy = -grav*cord(ii,igrav)
       else
          p_energy = 0.0d0
       endif 
        pl = phi(ii)
        tl = t(ii)
        pcl = pci(ii)


c check for identical calls   (established in h2o calls)
c assume the same criteria for air props
c
       if(ieval_flag(ii).ne.0.or.eval_test_air.ne.0) then          
       ieosd = ieos(ii)        
c gaz 070820 notes
c solubility is henry's law (temperature dependent)type 
        if(ieosd.le.2) then
         call air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)
        else
         xnl = 0.0
         dxnlp = 0.0
         dxnlpc = 0.0
         dxnlt = 0.0
        endif   
        xnl_ngas(ii,1) =  xnl
        xnl_ngas(ii,2) =  dxnlp
        xnl_ngas(ii,3) =  dxnlt
        xnl_ngas(ii,4) =  dxnlpc
c enthalpy of air under liquid pressure or vapor pressure (no pressure dependance
c gaz enthapy of solution = 0
c gaz 082921 looks like need to add enthalpy at T= 0.0        
        call air_cp(tl, cpa, dcpat)
        hcg=1.e-6*cpa*tl + hcg_zero
        hcl=hcg+hsol
        dhcgp=0.0
        dhcgt = 1.e-6 * (cpa + dcpat * tl) 
        
       enth_ngas(ii,1) = hcg
       enth_ngas(ii,2) = dhcgp 
       enth_ngas(ii,3) = dhcgt 
c-----------------------------------------------------------------------            
c air density       
       pcl0=0.101325d0
       roc0=1.292864d0
**** density of air eqn (31)
       drocpc=roc0*(273./(tl+273.))/pcl0
       roc=drocpc*pcl
       droct=-roc/(tl+273.)

      den_ngas(ii,1)= roc
      den_ngas(ii,2)= drocpc
      den_ngas(ii,3)= droct
c
c air enthalpy under partial pressure for non table model the same as
c  enth_ngas(ii,1)      
c 
       enth_ngas(ii,4) = hcg
       enth_ngas(ii,5) = dhcgp 
       enth_ngas(ii,6) = dhcgt      
c      
c air viscosity     
c          
         xvisc=182.7e-07+0.41e-07*(tl-18.0)
         dxvisct=0.41e-07
         dxviscpc = 0.0
         dxviscp  = 0.0
c
      visc_ngas(ii,1)= xvisc
      visc_ngas(ii,2)= dxviscp
      visc_ngas(ii,3)= dxvisct
      
c save values for node ii for comparison with next update
      
      den_ngas_old(ii,1)= roc
      den_ngas_old(ii,2)= drocpc
      den_ngas_old(ii,3)= droct

      enth_ngas_old(ii,1)= hcg
      enth_ngas_old(ii,2)= dhcgp
      enth_ngas_old(ii,3)= dhcgt
      enth_ngas_old(ii,4)= hcg
      enth_ngas_old(ii,5)= dhcgp
      enth_ngas_old(ii,6)= dhcgt

      visc_ngas_old(ii,1)= xvisc
      visc_ngas_old(ii,2)= dxviscp
      visc_ngas_old(ii,3)= dxvisct
      
      xnl_ngas_old(ii,1)= xnl
      xnl_ngas_old(ii,2)= dxnlp
      xnl_ngas_old(ii,3)= dxnlt
      xnl_ngas_old(ii,4)= dxnlpc

c save values for node ii for comparison with ii+1
c if abs(tl-tl_tast) < tol_t and    abs(pl-pl_tast) < tol_p use previous
c densities,enthalpies,viscosities      
      	roc_last	  =     roc   
      	drocpc_last   =     drocpc	   
      	droct_last    =     droct	 
          
      	hcg_last	  =     hcg   
      	dhcgp_last    =     dhcgp	   
      	dhcgt_last    =     dhcgt	
          hcg_last4	  =     hcg4   
      	dhcgp_last5   =     dhcgp5	   
      	dhcgt_last6   =     dhcgt6
	   
      	xvisc_last    =     xvisc	   
      	dxviscp_last  =     dxviscp   
      	dxvisct_last  =     dxvisct	 
          
          xnl_last	  =     xnl   
      	dxnlp_last    =     dxnlp	   
      	dxnlt_last    =     dxnlt
          dxnlpc_last   =     dxnlpc

      elseif(var_old_flag) then  
c use old unchanged values  if p,t(node) = p,t(node-1) 
c last properties rol-dvisvt stored on com_prop_data
c or if p,t(node) = p,t(node,last iter)
      
        den_ngas(ii,1)  = den_ngas_old(ii,1)
        den_ngas(ii,2)  = den_ngas_old(ii,2)
        den_ngas(ii,3)  = den_ngas_old(ii,3)
        enth_ngas(ii,1) = enth_ngas_old(ii,1)
        enth_ngas(ii,2) = enth_ngas_old(ii,2)
        enth_ngas(ii,3) = enth_ngas_old(ii,3)
        enth_ngas(ii,4) = enth_ngas_old(ii,4)
        enth_ngas(ii,5) = enth_ngas_old(ii,5)
        enth_ngas(ii,6) = enth_ngas_old(ii,6)        
        visc_ngas(ii,1) = visc_ngas_old(ii,1)
        visc_ngas(ii,2) = visc_ngas_old(ii,2)
        visc_ngas(ii,3) = visc_ngas_old(ii,3) 
        xnl_ngas(ii,1) = xnl_ngas_old(ii,1)	   
        xnl_ngas(ii,2) = xnl_ngas_old(ii,2)
        xnl_ngas(ii,3) = xnl_ngas_old(ii,3)
        xnl_ngas(ii,4) = xnl_ngas_old(ii,2)
      elseif(var_last_flag) then
          
        den_ngas(ii,1)  = roc_last
        den_ngas(ii,2)  = drocpc_last
        den_ngas(ii,3)  = droct_last
        enth_ngas(ii,1) = hcg_last
        enth_ngas(ii,2) = dhcgp_last
        enth_ngas(ii,3) = dhcgt_last
        enth_ngas(ii,4) = hcg_last4
        enth_ngas(ii,5) = dhcgp_last5
        enth_ngas(ii,6) = dhcgt_last6
        visc_ngas(ii,1) = xvisc_last
        visc_ngas(ii,2) = dxviscp_last
        visc_ngas(ii,3) = dxvisct_last 
        xnl_ngas(ii,1) = xnl_last	   
        xnl_ngas(ii,2) = dxnlp_last
        xnl_ngas(ii,3) = dxnlt_last
        xnl_ngas(ii,4) = dxnlpc_last
      endif 
c      
      enddo
      else if(iflg.eq.-1) then
c start counting and other metrics  
          ic_eval = 0
          itot_calls = 0
      else if(iflg.eq.-2) then
c end counting and other metrics 
       ratio_eval_tot =  ic_eval/(itot_calls + 1.e-9)  
c       write(ierr,100) ic_eval, itot_calls, ratio_eval_tot
      if(iout.ne.0) write(iout,100) ic_eval,itot_calls,ratio_eval_tot
      if(iptty.ne.0) write(iptty,100) ic_eval,itot_calls,ratio_eval_tot
c gaz 032922 increase the integer format i9 to i12          
100    format(/,' number of unique property calcs',1x,i12,1x,/,
     &    ' total number of prossible calls',1x,i12,1x,/,
     &    ' fraction of calcs to total calls',1x,f8.3) 
c      stop 
      endif
      return
      end   
      subroutine air_props_table(iflg, istart, iend, prop, phase)
c      
c subroutine to evaluate air properties using tables
c gaz 070921 initial coding
c air property calls assume water property calls are made earlier
c call optimization probably workd only for nodal sequence (all water and air props for m=1,2,..)  
c note that the derivatives wrt t and f(t,p) dfdt, dfdp (value_a(1:3) wheeas the storage for
c den_ngas(:,:) is f(t,p) dfdp, dfdt   
c      
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, iout, p_tol, t_tol,
     &     pc_tol, iptty
          use combi, only: cord
          use comdi, only: phi, t, ieos
          use comdti, only: n0
          use comfi, only: pci
          use comii, only: pmin_air_tabl

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate, n_allocate_iso
          character*9 prop, phase

          integer mid,mi,ieosd,kq

      	  real*8 pl,tl,psatl,dtsatp,dpsats,xv
            real*8 vis,xvisv,xvisl
            real*8 xrl,xrv
            real*8 rov,rol,ros
            real*8 enw,enl,env,ens
            real*8 dvisvt,dvisvp,dvislt,dvislp
            real*8 drovt,drovp,drovpc,dvvpc
            real*8 drvp,drv
            real*8 drolt,drolp
            real*8 drlp,drl
            real*8 dhvt,dhvp,dhvpc,dhlt,dhlp   
            real*8 pv, dhlpc,drolpc,dvlpc
            real*8 p_energy, pcl, pcl_in
            
            real*8 hsol,dhsolp,dhsolt
            real*8 hcl

c gaz 070321            
            real*8 dum1,dumb,dumc,value(9),value_a(9)
            integer istate, ifail, isol_print   

c gaz 092821 
            integer eval_test
            parameter(eval_test = 1)
        if(iflg.eq.0) then
c initialize and allocate memory
         if(.not.allocated(den_ngas)) then
          if(ico2.lt.0) then
           n_allocate_iso = n0
           n_allocate = 1
          else
           n_allocate_iso = n0
           n_allocate = n0
          endif
          allocate(den_ngas(n_allocate_iso,6))
          allocate(enth_ngas(n_allocate,6))
          allocate(visc_ngas(n_allocate_iso,6))
          allocate(xnl_ngas(n_allocate,6))
          allocate(den_ngas_old(n_allocate_iso,6))
          allocate(enth_ngas_old(n_allocate,6))
          allocate(visc_ngas_old(n_allocate_iso,6))
          allocate(xnl_ngas_old(n_allocate,6))
         endif
         if(.not.allocated(phi_old)) then
          allocate(phi_old(n_allocate_iso))
          allocate(t_old(n_allocate))
          allocate(pci_old(n_allocate_iso))
          allocate(ieval_flag(n_allocate))
          ieval_flag = 0
          phi_old = 0.0d0
          t_old = 0.0d0
          pci_old = 0.0d0
         endif 
         return
        else if(iflg.eq.1) then
c             
c loop structure             
c 

        do ii = istart,iend
        
c          
c        liquid phase table data
c      

       if(igrav.ne.0) then
          p_energy = -grav*cord(ii,igrav)
       else
          p_energy = 0.0d0
       endif 
        pl = phi(ii)
        tl = t(ii)
        pcl = pci(ii)


c check for identical calls   (established in h2o calls)
c assume the same criteria for air props
c
       if(ieval_flag(ii).ne.0.or.eval_test.eq.1) then          
       ieosd = ieos(ii) 
c ngas solubility
       
c gaz 070820 notes
c solubility is henry's law (temperature dependent)type 
        if(ieosd.le.2) then
         call air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)
        else
         xnl = 0.0
         dxnlp = 0.0
         dxnlpc = 0.0
         dxnlt = 0.0
        endif   
        xnl_ngas(ii,1) =  xnl
        xnl_ngas(ii,2) =  dxnlp
        xnl_ngas(ii,3) =  dxnlt
        xnl_ngas(ii,4) =  dxnlpc
c ngas table properties  
c enthalpy of air under liquid pressure
c gaz enthapy of solution = 0
        hsol=0.0
        dhsolt=0.0
        dhsolp=0.0   
        value_a = 0.0
        call air_properties_new(4,3,pl,tl,dum1,istate,
     &                 dumb,value_a,dumc)  

       hcg = value_a(4) + p_energy
       dhcgt = value_a(5)
       dhcgp = value_a(6)   
       enth_ngas(ii,1) = hcg
       enth_ngas(ii,2) = dhcgp 
       enth_ngas(ii,3) = dhcgt 
c-----------------------------------------------------------------------       
     
c air density     
c note pcl, not pl, is used
       value_a = 0.0          

       if(pcl.le.pmin_air_tabl(1)) then
        pcl_in = pmin_air_tabl(1)
       else
        pcl_in = pcl
       endif
       call air_properties_new(4,3,pcl_in,tl,dum1,istate,
     &                 dumb,value_a,dumc)            
      roc = value_a(1)
      droct = value_a(2)
      drocpc = value_a(3) 
      den_ngas(ii,1)= roc
      den_ngas(ii,2)= drocpc
      den_ngas(ii,3)= droct
c
c air enthalpy under partial pressure (note 4,5,6)
c 
       hcg4   = value_a(4)
       dhcgt6 = value_a(5)
       dhcgp5 = value_a(6)
       enth_ngas(ii,4) = hcg4
       enth_ngas(ii,5) = dhcgp5 
       enth_ngas(ii,6) = dhcgt6      
c      
c air viscosity     
c          
      xvisc = value_a(7)
      dxvisct = value_a(8)
      dxviscp = value_a(9)
c
      visc_ngas(ii,1)= xvisc
      visc_ngas(ii,2)= dxviscp
      visc_ngas(ii,3)= dxvisct
      
c save values for node ii for comparison with next update
      
      den_ngas_old(ii,1)= roc
      den_ngas_old(ii,2)= drocpc
      den_ngas_old(ii,3)= droct

      enth_ngas_old(ii,1)= hcg
      enth_ngas_old(ii,2)= dhcgp
      enth_ngas_old(ii,3)= dhcgt
      enth_ngas_old(ii,4)= hcg4
      enth_ngas_old(ii,5)= dhcgp5
      enth_ngas_old(ii,6)= dhcgt6

      visc_ngas_old(ii,1)= xvisc
      visc_ngas_old(ii,2)= dxviscp
      visc_ngas_old(ii,3)= dxvisct
      
      xnl_ngas_old(ii,1)= xnl
      xnl_ngas_old(ii,2)= dxnlp
      xnl_ngas_old(ii,3)= dxnlt
      xnl_ngas_old(ii,4)= dxnlpc

c save values for node ii for comparison with ii+1
c if abs(tl-tl_tast) < tol_t and    abs(pl-pl_tast) < tol_p use previous
c densities,enthalpies,viscosities      
      	roc_last	  =     roc   
      	drocpc_last   =     drocpc	   
      	droct_last    =     droct	 
          
      	hcg_last	  =     hcg   
      	dhcgp_last    =     dhcgp	   
      	dhcgt_last    =     dhcgt	
          hcg_last4	  =     hcg4   
      	dhcgp_last5   =     dhcgp5	   
      	dhcgt_last6   =     dhcgt6
	   
      	xvisc_last    =     xvisc	   
      	dxviscp_last  =     dxviscp   
      	dxvisct_last  =     dxvisct	 
          
          xnl_last	  =     xnl   
      	dxnlp_last    =     dxnlp	   
      	dxnlt_last    =     dxnlt
          dxnlpc_last   =     dxnlpc

      elseif(var_old_flag) then  
c use old unchanged values  if p,t(node) = p,t(node-1) 
c last properties rol-dvisvt stored on com_prop_data
c or if p,t(node) = p,t(node,last iter)
      
        den_ngas(ii,1)  = den_ngas_old(ii,1)
        den_ngas(ii,2)  = den_ngas_old(ii,2)
        den_ngas(ii,3)  = den_ngas_old(ii,3)
        enth_ngas(ii,1) = enth_ngas_old(ii,1)
        enth_ngas(ii,2) = enth_ngas_old(ii,2)
        enth_ngas(ii,3) = enth_ngas_old(ii,3)
        enth_ngas(ii,4) = enth_ngas_old(ii,4)
        enth_ngas(ii,5) = enth_ngas_old(ii,5)
        enth_ngas(ii,6) = enth_ngas_old(ii,6)        
        visc_ngas(ii,1) = visc_ngas_old(ii,1)
        visc_ngas(ii,2) = visc_ngas_old(ii,2)
        visc_ngas(ii,3) = visc_ngas_old(ii,3) 
        xnl_ngas(ii,1) = xnl_ngas_old(ii,1)	   
        xnl_ngas(ii,2) = xnl_ngas_old(ii,2)
        xnl_ngas(ii,3) = xnl_ngas_old(ii,3)
        xnl_ngas(ii,4) = xnl_ngas_old(ii,2)
      elseif(var_last_flag) then
          
        den_ngas(ii,1)  = roc_last
        den_ngas(ii,2)  = drocpc_last
        den_ngas(ii,3)  = droct_last
        enth_ngas(ii,1) = hcg_last
        enth_ngas(ii,2) = dhcgp_last
        enth_ngas(ii,3) = dhcgt_last
        enth_ngas(ii,4) = hcg_last4
        enth_ngas(ii,5) = dhcgp_last5
        enth_ngas(ii,6) = dhcgt_last6
        visc_ngas(ii,1) = xvisc_last
        visc_ngas(ii,2) = dxviscp_last
        visc_ngas(ii,3) = dxvisct_last 
        xnl_ngas(ii,1) = xnl_last	   
        xnl_ngas(ii,2) = dxnlp_last
        xnl_ngas(ii,3) = dxnlt_last
        xnl_ngas(ii,4) = dxnlpc_last
      endif 
c      
      enddo
      else if(iflg.eq.-1) then
c start counting and other metrics  
          ic_eval = 0
          itot_calls = 0
      else if(iflg.eq.-2) then
c end counting and other metrics 
       ratio_eval_tot =  ic_eval/(itot_calls + 1.e-9)  
c       write(ierr,100) ic_eval, itot_calls, ratio_eval_tot
      if(iout.ne.0) write(iout,100) ic_eval,itot_calls,ratio_eval_tot
      if(iptty.ne.0) write(iptty,100) ic_eval,itot_calls,ratio_eval_tot
c gaz 032922 increase the integer format i9 to i12          
100    format(/,' number of unique property calcs',1x,i12,1x,/,
     &    ' total number of prossible calls',1x,i12,1x,/,
     &    ' fraction of calcs to total calls',1x,f8.3) 
c      stop 
      endif
      return
      end            
c gaz 081421 added co2 capability 
        subroutine co2wh_props_table(iflg, istart, iend, prop, phase)
c      
c subroutine to evaluate co2 properties using tables
c gaz 081421  initial coding
c co2 uses same structure and variables as air
c air property calls assume water property calls are made earlier
c call optimization probably workd only for nodal sequence (all water and air props for m=1,2,..)  
c note that the derivatives wrt t and f(t,p) dfdt, dfdp (value_a(1:3) whereas the storage for
c den_ngas(:,:) is f(t,p) dfdp, dfdt   
c      
          use com_prop_data
          use comai, only: ico2, ierr, igrav, grav, iout, p_tol, t_tol,
     &     pc_tol,tcrit_h2o,iptty
          use combi, only: cord
          use comdi, only: phi, t, ieos, an, to
          use comdti, only: n0
          use comfi, only: pci,cnlf,cnvf
          use comii, only: pmin_co2wh_tabl

          implicit none
             
          integer iflg, istart, iend, ii, n_allocate
          character*9 prop, phase

          integer mid,mi,ieosd,kq

      	  real*8 pl,tl,psatl,dtsatp,dpsats,xv
            real*8 vis,xvisv,xvisl
            real*8 xrl,xrv
            real*8 rov,rol,ros
            real*8 enw,enl,env,ens
            real*8 dvisvt,dvisvp,dvislt,dvislp
            real*8 drovt,drovp,drovpc,dvvpc
            real*8 drvp,drv
            real*8 drolt,drolp
            real*8 drlp,drl
            real*8 dhvt,dhvp,dhvpc,dhlt,dhlp   
            real*8 pv, dhlpc,drolpc,dvlpc
            real*8 p_energy, pcl, pcl_in
            
            real*8 hsol,dhsolp,dhsolt
            real*8 hcl

c gaz 070321            
            real*8 dum1,dumb,dumc,value(9),value_a(9)
c gaz 092721  
            integer eval_test
            parameter(eval_test=1)
            real*8 henry_coeff, dpsatt,tratio,dtratiot
c gaz 093021 for taylor series version of solubility based on mao et al (2013)
            real*8 xnl0,dxnl0t,dxnl0p,tlm0,plm0
            
            integer istate, ifail, isol_print    
        if(iflg.eq.0) then
c initialize and allocate memory
         if(.not.allocated(den_ngas)) then
          if(ico2.lt.0) then
           n_allocate = 1
          else
           n_allocate = n0
          endif
          allocate(den_ngas(n_allocate,6))
          allocate(enth_ngas(n_allocate,6))
          allocate(visc_ngas(n_allocate,6))
          allocate(xnl_ngas(n_allocate,6))
          allocate(den_ngas_old(n_allocate,6))
          allocate(enth_ngas_old(n_allocate,6))
          allocate(visc_ngas_old(n_allocate,6))
          allocate(xnl_ngas_old(n_allocate,6))
         endif
         if(.not.allocated(phi_old)) then
          allocate(phi_old(n_allocate))
          allocate(t_old(n_allocate))
          allocate(pci_old(n_allocate))
          allocate(ieval_flag(n_allocate))
          ieval_flag = 0
          phi_old = 0.0d0
          t_old = 0.0d0
          pci_old = 0.0d0
         endif 
         return
        else if(iflg.eq.1) then
c             
c loop structure             
c 

        do ii = istart,iend
        
c          
c        liquid phase table data
c      

       if(igrav.ne.0) then
          p_energy = -grav*cord(ii,igrav)
       else
          p_energy = 0.0d0
       endif 
        pl = phi(ii)
        tl = t(ii)
        pcl = pci(ii)


c check for identical calls   (established in h2o calls)
c assume the same criteria for air props
c
       if(ieval_flag(ii).ne.0.or.eval_test.eq.1) then          
       ieosd = ieos(ii) 
c ngas solubility
       
c gaz 081421 this needs to be changed for co2
c solubility is henry's law (temperature dependent)type 
c derived from mao 2013 (programmed by M Osullivan)
c see com_prop_data for parameter description
        if(ieosd.le.2) then
          tlm0 = tlm(1)  
          plm0 = plm(1)
          xnl0 = xnlm(imaox(1,1))
          dxnl0t = dxnlmt(imaodxt(1,1))
          dxnl0p = dxnlmp(imaodxp(1,1))
c          dxnl0t = 0.0
c          xnl0  = 0.0
c          plm0 = 0.0
c          tlm0 = 0.0          
c          xnl = xnl0 + dxnl0t*(to(ii)-tlm0) + dxnl0p*(pcl-plm0) 
          xnl = xnl0 + dxnl0t*(to(ii)-tlm0) + dxnl0p*(pcl-plm0) 
c gaz 100121 test   dxnl0p       
          if(xnl.gt.0.d0) then
           dxnlpc = dxnl0p
           dxnlp=0.0
c           dxnlt = dxnl0t
          else
           xnl = 0.0
           dxnlp = 0.0
           dxnlpc = dxnl0p
           dxnlt = 0.0
          endif     
        else
         xnl = 0.0
         dxnlp = 0.0
         dxnlpc = 0.0
         dxnlt = 0.0
        endif   
        xnl_ngas(ii,1) =  xnl
        xnl_ngas(ii,2) =  dxnlp
        xnl_ngas(ii,3) =  dxnlt
        xnl_ngas(ii,4) =  dxnlpc
c ngas table properties  
c enthalpy of air under liquid pressure
c gaz enthapy of solution = 0
        hsol=0.0
        dhsolt=0.0
        dhsolp=0.0   
        value_a = 0.0
        call co2wh_properties_new(4,3,pl,tl,dum1,istate,
     &                 dumb,value_a,dumc)  

       hcg = value_a(4) + p_energy
       dhcgt = value_a(5)
       dhcgp = value_a(6)   
       enth_ngas(ii,1) = hcg
       enth_ngas(ii,2) = dhcgp 
       enth_ngas(ii,3) = dhcgt 
c-----------------------------------------------------------------------       
     
c co2 density     
c note pcl, not pl, is used
       value_a = 0.0          

       if(pcl.le.pmin_co2wh_tabl(1)) then
        pcl_in = pmin_co2wh_tabl(1)
       else
        pcl_in = pcl
       endif
       call co2wh_properties_new(4,3,pcl_in,tl,dum1,istate,
     &                 dumb,value_a,dumc)            
      roc = value_a(1)
      droct = value_a(2)
      drocpc = value_a(3) 
      den_ngas(ii,1)= roc
      den_ngas(ii,2)= drocpc
      den_ngas(ii,3)= droct
c
c co2 enthalpy under partial pressure (note 4,5,6)
c 
       hcg4   = value_a(4)
       dhcgt6 = value_a(5)
       dhcgp5 = value_a(6)
       enth_ngas(ii,4) = hcg4
       enth_ngas(ii,5) = dhcgp5 
       enth_ngas(ii,6) = dhcgt6      
c      
c co2 viscosity     
c          
      xvisc = value_a(7)
      dxvisct = value_a(8)
      dxviscp = value_a(9)
c
      visc_ngas(ii,1)= xvisc
      visc_ngas(ii,2)= dxviscp
      visc_ngas(ii,3)= dxvisct
      
c save values for node ii for comparison with next update
      
      den_ngas_old(ii,1)= roc
      den_ngas_old(ii,2)= drocpc
      den_ngas_old(ii,3)= droct

      enth_ngas_old(ii,1)= hcg
      enth_ngas_old(ii,2)= dhcgp
      enth_ngas_old(ii,3)= dhcgt
      enth_ngas_old(ii,4)= hcg4
      enth_ngas_old(ii,5)= dhcgp5
      enth_ngas_old(ii,6)= dhcgt6

      visc_ngas_old(ii,1)= xvisc
      visc_ngas_old(ii,2)= dxviscp
      visc_ngas_old(ii,3)= dxvisct
      
      xnl_ngas_old(ii,1)= xnl
      xnl_ngas_old(ii,2)= dxnlp
      xnl_ngas_old(ii,3)= dxnlt
      xnl_ngas_old(ii,4)= dxnlpc

c save values for node ii for comparison with ii+1
c if abs(tl-tl_tast) < tol_t and    abs(pl-pl_tast) < tol_p use previous
c densities,enthalpies,viscosities      
      	roc_last	  =     roc   
      	drocpc_last   =     drocpc	   
      	droct_last    =     droct	 
          
      	hcg_last	  =     hcg   
      	dhcgp_last    =     dhcgp	   
      	dhcgt_last    =     dhcgt	
          hcg_last4	  =     hcg4   
      	dhcgp_last5   =     dhcgp5	   
      	dhcgt_last6   =     dhcgt6
	   
      	xvisc_last    =     xvisc	   
      	dxviscp_last  =     dxviscp   
      	dxvisct_last  =     dxvisct	 
          
          xnl_last	  =     xnl   
      	dxnlp_last    =     dxnlp	   
      	dxnlt_last    =     dxnlt
          dxnlpc_last   =     dxnlpc

      elseif(var_old_flag) then  
c use old unchanged values  if p,t(node) = p,t(node-1) 
c last properties rol-dvisvt stored on com_prop_data
c or if p,t(node) = p,t(node,last iter)
      
        den_ngas(ii,1)  = den_ngas_old(ii,1)
        den_ngas(ii,2)  = den_ngas_old(ii,2)
        den_ngas(ii,3)  = den_ngas_old(ii,3)
        enth_ngas(ii,1) = enth_ngas_old(ii,1)
        enth_ngas(ii,2) = enth_ngas_old(ii,2)
        enth_ngas(ii,3) = enth_ngas_old(ii,3)
        enth_ngas(ii,4) = enth_ngas_old(ii,4)
        enth_ngas(ii,5) = enth_ngas_old(ii,5)
        enth_ngas(ii,6) = enth_ngas_old(ii,6)        
        visc_ngas(ii,1) = visc_ngas_old(ii,1)
        visc_ngas(ii,2) = visc_ngas_old(ii,2)
        visc_ngas(ii,3) = visc_ngas_old(ii,3) 
        xnl_ngas(ii,1) = xnl_ngas_old(ii,1)	   
        xnl_ngas(ii,2) = xnl_ngas_old(ii,2)
        xnl_ngas(ii,3) = xnl_ngas_old(ii,3)
        xnl_ngas(ii,4) = xnl_ngas_old(ii,2)
      elseif(var_last_flag) then
          
        den_ngas(ii,1)  = roc_last
        den_ngas(ii,2)  = drocpc_last
        den_ngas(ii,3)  = droct_last
        enth_ngas(ii,1) = hcg_last
        enth_ngas(ii,2) = dhcgp_last
        enth_ngas(ii,3) = dhcgt_last
        enth_ngas(ii,4) = hcg_last4
        enth_ngas(ii,5) = dhcgp_last5
        enth_ngas(ii,6) = dhcgt_last6
        visc_ngas(ii,1) = xvisc_last
        visc_ngas(ii,2) = dxviscp_last
        visc_ngas(ii,3) = dxvisct_last 
        xnl_ngas(ii,1) = xnl_last	   
        xnl_ngas(ii,2) = dxnlp_last
        xnl_ngas(ii,3) = dxnlt_last
        xnl_ngas(ii,4) = dxnlpc_last
      endif 
c      
      enddo
      else if(iflg.eq.-1) then
c start counting and other metrics  
          ic_eval = 0
          itot_calls = 0
      else if(iflg.eq.-2) then
c end counting and other metrics 
       ratio_eval_tot =  ic_eval/(itot_calls + 1.e-9)  
c       write(ierr,100) ic_eval, itot_calls, ratio_eval_tot
      if(iout.ne.0) write(iout,100) ic_eval,itot_calls,ratio_eval_tot
      if(iptty.ne.0) write(iptty,100) ic_eval,itot_calls,ratio_eval_tot
c gaz 032922 increase the integer format i9 to i12      
100    format(/,' number of unique property calcs',1x,i12,1x,/,
     &    ' total number of prossible calls',1x,i12,1x,/,
     &    ' fraction of calcs to total calls',1x,f8.3) 
c      stop 
      endif
      return
      end 


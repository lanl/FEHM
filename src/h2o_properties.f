      subroutine h2o_properties(iflg,iphase,var1,var2,var3,
     &     var4,prop,der1,der2,der3,der4)
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
!D1 To calculate water properties and derivatives (new subroutine form).
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/methane_properties.f_a  $
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 ?
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

      use comco2
      use comai, only : itsat 
      use comrxni, only : cden_flag

      implicit none

      integer iflg,iphase
      real*8 var1,var2,var3,var4
      real*8 prop,der1,der2,der3,der4
      real*8 tm0,pm0
      real*8 cpms,cpml,cpmv
      real*8 lhmsl,lhmlv
      real*8 denms,ddenmsp,ddenmst
      real*8 denml,ddenmlp,ddenmlt
      real*8 denmv,prefmv,pterm,tterm   
c     sold phase is immobile (set visms very large)
      real*8 visms  
      real*8 visml,pvisml1,pvisml2,tvisml
      real*8 vismv,dvismvp,dvismvt    
      real*8 tphasesl,tphaselv        
      real*8 dtphaseslp,dtphaselvp        
      real*8 pphasesl,pphaselv        
      real*8 cap_max                  
      real*8 fracw_min,fracw_max            
c     remove following declarations later
      real*8 rolref, comw, rcomd,pref,ela0,elpa1
      real*8 elpa2,elpa3,elta1,elta2,elta3,elpta,elp2ta,elpt2a
      real*8 elb0,elpb1,elpb2,elpb3,eltb1,eltb2,eltb3,elptb
      real*8 elp2tb,elpt2b,dla0,dlpa1,dlpa2,dlpa3,dlta1,dlta2,dlta3
      real*8 dlpta,dlp2ta,dlpt2a,dlb0,dlpb1,dlpb2,dlpb3,dltb1,dltb2
      real*8 dltb3,dlptb,dlp2tb,dlpt2b,vla0,vlpa1,vlpa2,vlpa3,vlta1
      real*8 enwn1,enwn2,enwn3,enwn,enwd1,enwd2,enwd3
      real*8 enwd,dhwpn1,dhwpn2,dhwpn,dhwpd,dhwp,dhwtn1
      real*8 dhwtn2,dhwtn,dhwtd,dhwt,rnwn1,rnwn2,rnwn3
      real*8 rnwd1,rnwd2,rnwd3,rnwn,rnwd,rnw,drlpn1,drlpn2
      real*8 drlpn,drolpd,drlen1,drlen2,drlen,droled
      real*8 vlta2,vlta3
      real*8 vlpta,vlp2ta,vlpt2a
      real*8 vlb0,vlpb1,vlpb2,vlpb3,vltb1,vltb2,vltb3
      real*8 vlptb,vlp2tb,vlpt2b
      real*8 x,x2,x3,tl,tl2,tl3,tref,tlx,tl2x,tlx2
      real*8 viln1,viln2,viln3,viln,vil
      real*8 vild1,vild2,vild3,vild
      real*8 dvlpn1,dvlpn2,dvlpn,dvilpd,dvlen1,dvlen2,dvlen,dviled
      real*8 eva0,evpa1,evpa2
      real*8 evpa3,evta1,evta2,evta3,evpta,evp2ta,evpt2a,evb0,evpb1
      real*8 evpb2,evpb3,evtb1,evtb2,evtb3,evptb,evp2tb,evpt2b
      real*8 ensn,ensn1,ensn2,ensn3,ensd1,ensd2,ensd3
      real*8 ensd,ens,dhvp1,dhvp2,dhvpn,dhvpd,dhvt1,dhvt2
      real*8 dhvtn,dhvtd,dva0,rnsn1,rnsn2,rnsn3,rnsd1,rnsd2
      real*8 rnsd3,rnsn,rnsd,rns,drspn1,drspn2,drspn,drospd
      real*8 drsen1,drsen2,drsen,drostd,visn1,visn2,visn3,visn,visd1
      real*8 visd2,visd3,visd,vis,xvisv,dvspn1,dvspn2,dvspn,dvispd
      real*8 dvsen1,dvsen2,dvsen,dvised
      real*8 dvpa1,dvpa2,dvpa3,dvta1,dvta2,dvta3,dvpta,dvp2ta,dvpt2a
      real*8 dvb0,dvpb1,dvpb2,dvpb3,dvtb1,dvtb2,dvtb3,dvptb,dvp2tb
      real*8 dvpt2b,vva0,vvpa1,vvpa2,vvpa3,vvta1,vvta2,vvta3,vvpta
      real*8 vvp2ta,vvpt2a,vvb0,vvpb1,vvpb2,vvpb3,vvtb1,vvtb2,vvtb3
      real*8 vvptb,vvp2tb,vvpt2b
      real*8 tsa0,tspa1,tspa2,tspa3,tspa4,tsb0,tspb1,tspb2,tspb3
      real*8 tspb4,tfunn,tfund,tfun,dtpsn,dtpsd      
      real*8 del_h, mol, derdel_h,x1,der,prop1,h_salt,dh_salt,prop2
      real*8 a(0:3,0:2)
      real*8 y,x4,y1,y2,y3,y4,a1,a2,a3,b1,b2,b3,b4,b5,guess,bsl,asl,f,sl
      real*8 df,dx,t1

      real*8 vo, phistar, kappa, vc, phi, num1, den1, dvodp, dphistardp
      real*8 dkappadp, dphidp, dnumdp, ddendp, dvodt, dphistardt
      real*8 dkappadt, dphidt, dnumdt, var5, ddendt
      integer i, j

c     cpms - heat capacity for methane solid
c     cpml - heat capacity for methane liquid
c     cpmv - heat capacity for methane vapor 
c     tm0 - reference temperature for methane
c     lhmsl - latent heat for methane solid to liquid phase change
c     lhmlv - latent heat for methane liquid to vapor phase change
      parameter(cpms=4.0d-3,cpml=4.0d-3,cpmv=0.002,pm0=0.0,tm0=0.0)
      parameter(lhmsl=0.1,lhmlv=2.0)                      
c     denms - reference density for methane solid
c     ddenmsp - derivative of solid density wrt pressure 
c     ddenmst - derivative of solid density wrt temperature
      parameter(denms=1000.,ddenmsp=0.5,ddenmst=-0.5)
c     denml - reference density for methane liquid
c     ddenmlp - derivative of liquid density wrt pressure 
c     ddenmlt - derivative of liquid density wrt temperature
      parameter(denml=1000.,ddenmlp=0.5,ddenmlt=-1.0)
c     denmv - reference density for methane gas    
c     prefmv - reference pressure for perfect gas law   
c     pterm,tterm - intermediate calculations
      parameter(denmv=1.0,prefmv=0.1)       
c     visms  reference viscosity for methane solid  
      parameter(visms=1.d10)                                     
c     visml - reference viscosity for methane liquid
c     pvisml1 - parameter used with liquid viscosity 
c     pvisml2 - parameter used with liquid viscosity 
c     tvisml - cutoff temperature for constant liquid viscosity 
c     parameter(visml=0.005,pvisml1=5.0)
c     parameter(pvisml2=0.1,tvisml=0.0001)
      parameter(visml=0.005,pvisml1=0.0017389)
      parameter(pvisml2=-0.0000457,tvisml=0.0001)
c     vismv- reference  viscosity for methane vapor  
c     dvismvp - derivative of vapor viscosity wrt pressure 
c     dvismvt - derivative of vapor viscosity wrt temperature
      parameter(vismv=1.d-5,dvismvp=-8.e-7,dvismvt=4.d-8)
c     tphasesl - phase change temperature for solid to liquid 
c     tphaselv - phase change temperature for liquid to vapor 
c     dtphaseslp - derivative phase change temperature (s-l) wrt pressure
c     dtphaselvp - derivative phase change temperature (l-v) wrt pressure 
c     parameter(tphasesl=0.0,tphaselv=100.0)                   
      parameter(tphasesl=-100.0,tphaselv=800.0)                   
c     parameter(dtphaseslp=1.0,dtphaselvp=250.0)                   
      parameter(dtphaseslp=1.0e-5,dtphaselvp=250.0e-5)                   
c     pphasesl - phase change pressure for solid to liquid 
c     pphaselv - phase change pressure for liquid to vapor 
      parameter(pphasesl=0.1,pphaselv=0.1)                   
c     cap_max - maximum capillary pressure wrt water_frac 
c     parameter(cap_max=1e-5)                  
c     parameter(cap_max=2.0)                  
      parameter(cap_max=0.1)                  
c     fracw_min - minimum water fraction for methane production
      parameter(fracw_min=0.01)               
c     fracw_max - maximum water fraction for methane production
      parameter(fracw_max=1.01)               

c     input ----------
c     iflag - designator for type of calculation
c     iphase - phase state                                   
c     var1 - variable 1(pressure)
c     var2 - variable 2(temperature)
c     var3 - variable 3
c     var4 - variable 4
c     output ----------
c     prop - property of material (methane)
c     der1 - derivative of property wrt variable 1
c     der2 - derivative of property wrt variable 2
c     der3 - derivative of property wrt variable 3
c     der4 - derivative of property wrt variable 4
c     -----------------

c     iflg=1, enthalpy and derivative wrt pressure temperature
c     iflg=2, density and derivative wrt pressure temperature
c     iflg=3, viscosity and derivative wrt pressure temperature
c     iflg=4, saturation (or melting) temperature and derivative wrt pressure
c     iflg=5, saturation (or melting) pressure and derivative wrt temperature
c     iflg=10, capillary pres (water frac) and derivative wrt water frac

      if(iflg.eq.1.and.iphase.eq.1) then
c     solid enthalpy and derivative wrt pressure and temperature
         prop=cpms*(var2-tm0)
         der1 = 0.0
         der2 = cpms
      else if(iflg.eq.1.and.iphase.eq.2) then
c     liquid enthalpy and derivative wrt pressure and temperature
       if(itsat.le.10) then 
         x=var1
         x2=x*x
         x3=x2*x
         x4=x3*x
         tl=var2
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
         ela0=0.25623465d-03
         elpa1=0.10184405d-02
         elpa2=0.22554970d-04
         elpa3=0.34836663d-07
         elta1=0.41769866d-02
         elta2=-0.21244879d-04
         elta3=0.25493516d-07
         elpta=0.89557885d-04
         elp2ta=0.10855046d-06
         elpt2a=-0.21720560d-06
c     denomenator coefficients
         elb0=0.10000000d+01
         elpb1=0.23513278d-01
         elpb2=0.48716386d-04
         elpb3=-0.19935046d-08
         eltb1=-0.50770309d-02
         eltb2=0.57780287d-05
         eltb3=0.90972916d-09
         elptb=-0.58981537d-04
         elp2tb=-0.12990752d-07
         elpt2b=0.45872518d-08
c     liquid enthalpy
         enwn1=ela0+elpa1*x+elpa2*x2+elpa3*x3
         enwn2=elta1*tl+elta2*tl2+elta3*tl3
         enwn3=elpta*tlx+elpt2a*tl2x+elp2ta*tlx2
         enwn=enwn1+enwn2+enwn3
         enwd1=elb0+elpb1*x+elpb2*x2+elpb3*x3
         enwd2=eltb1*tl+eltb2*tl2+eltb3*tl3
         enwd3=elptb*tlx+elpt2b*tl2x+elp2tb*tlx2
         enwd=enwd1+enwd2+enwd3
         enw=enwn/enwd
         prop=enw

c     derivatives of enthalpy
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
         der2=dhwt
         der1=dhwp
       else if(itsat.gt.10) then
c   call  pseudo-vap eos  
          tl = var2
          x = var1      
          call eos_aux(itsat,tl,x,1,1,prop,der2,der1)
       endif         
c     prop = 1037.7d0-4.75d0*(240.d0-var2)-0.1d0*(5.d0-var1)
c     der1 = 0.1d0
c     der2 = 4.75d0
c     prop=prop*1.d-3
c     der1 = der1*1.d-3
c     der2 = der2*1.d-3
c     expression from Michaelides
         if (ibrine.ne.0) then
            a(0,0) = -9633.6d0
            a(0,1) = -4080.d0
            a(0,2) = 286.49d0
            a(1,0) = 166.58d0
            a(1,1) = 68.577d0
            a(1,2) = -4.6856d0
            a(2,0) = -0.90963d0
            a(2,1) = -0.36524d0
            a(2,2) = 0.249667d-1
            a(3,0) = 0.17965d-2
            a(3,1) = 0.71924d-3
            a(3,2) = -0.49d-4
            if (cden_flag .eq. 2) then
c     Sum concentrations to get total moles/kg-water
               mol = var3
            else
c     Convert ppm of salt to moles/kg-water    
               mol = var3/58.44d3
            end if
c     'mol' in below equation is the molality
            del_h=0.d0
            derdel_h=0.d0
            do i = 0, 3
               do j = 0, 2
                  del_h=del_h+(a(i,j)*(var2**i)*(mol**j))
               enddo	
            enddo	
            do i = 1, 3
               do j = 0, 2
                  derdel_h=derdel_h+(a(i,j)*i*(var2**(i-1))*(mol**j))
               enddo
            enddo
            del_h=-1d-3*(4.184d0/(1000.d0+(58.44d0*mol)))*del_h
            derdel_h=1d-3*(4.184d0/(1000.d0+(58.44d0*mol)))*derdel_h
            x1=1000.d0/(1000.d0+58.44*mol)
            x2=58.44*mol/(1000.d0+58.44*mol)
            h_salt=41.293d0*(var2+273.15d0)+
     &           (3.3607d-2*((var2+273.15d0)**2.d0)/2.d0)
     &           -(1.3927d-5*((var2+273.15d0)**3.d0)/3.d0)
            dh_salt=41.293d0+
     &           (3.3607d-2*(var2+273.15d0))
     &           -(1.3927d-5*((var2+273.15d0)**2.d0))
c     convert units from J/Mole to MJ/Kg
            h_salt=h_salt*1d-3/58.44d0
            dh_salt=dh_salt*1d-3/58.44d0
            prop=x1*prop+x2*h_salt+mol*del_h
            der1=x1*der1
            der2=x1*der2+x2*dh_salt+mol*derdel_h
         endif
      else if(iflg.eq.1.and.iphase.eq.3) then
c     vapor enthalpy and derivative wrt pressure and temperature
         x=var1
         x2=x*x
         x3=x2*x
         x4=x3*x
         tl=var2
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
c     vapor enthalpy
         eva0=0.31290881d+00
         evpa1= -0.10000000d+01
         evpa2=  0.25748596d-01
         evpa3=  0.38846142d-03
         evta1=  0.11319298d-01
         evta2=  0.20966376d-04
         evta3=  0.74228083d-08
         evpta=  0.19206133d-02
         evp2ta= -0.10372453d-03
         evpt2a=  0.59104245d-07
         evb0=  0.12511319d+00
         evpb1= -0.36061317d+00
         evpb2=  0.58668929d-02
         evpb3=  0.99059715d-04
         evtb1=  0.44331611d-02
         evtb2=  0.50902084d-05
         evtb3= -0.10812602d-08
         evptb=  0.90918809d-03
         evp2tb= -0.26960555d-04
         evpt2b= -0.36454880d-06
         ensn1=eva0+evpa1*x+evpa2*x2+evpa3*x3
         ensn2=evta1*tl+evta2*tl2+evta3*tl3
         ensn3=evpta*tlx+evpt2a*tl2x+evp2ta*tlx2
         ensn=ensn1+ensn2+ensn3
         ensd1=evb0+evpb1*x+evpb2*x2+evpb3*x3
         ensd2=evtb1*tl+evtb2*tl2+evtb3*tl3
         ensd3=evptb*tlx+evpt2b*tl2x+evp2tb*tlx2
         ensd=ensd1+ensd2+ensd3
         ens=ensn/ensd
         prop=ens

c     derivatives of vapor enthalpy
         dhvp1=evpa1+2*evpa2*x+3*evpa3*x2+evpta*tl
         dhvp1=ensd*(dhvp1+2*evp2ta*tlx+evpt2a*tl2)
         dhvp2=evpb1+2*evpb2*x+3*evpb3*x2+evptb*tl
         dhvp2=ensn*(dhvp2+2*evp2tb*tlx+evpt2b*tl2)
         dhvpn=dhvp1-dhvp2
         dhvpd=ensd**2
         der1=dhvpn/dhvpd
         dhvt1=evta1+2*evta2*tl+3*evta3*tl2+evpta*x
         dhvt1=ensd*(dhvt1+evp2ta*x2+2*evpt2a*tlx)
         dhvt2=evtb1+2*evtb2*tl+3*evtb3*tl2+evptb*x
         dhvt2=ensn*(dhvt2+evp2tb*x2+2*evpt2b*tlx)
         dhvtn=dhvt1-dhvt2
         dhvtd=ensd**2
         der2=dhvtn/dhvtd
      else if(iflg.eq.1.and.iphase.eq.22) then

c     prop = 1045.77225136398d0-80.764d0*(3.5d0-var1)
         prop = 1049.8d0-80.764d0*(3.5d0-var1)
         der1 = 80.764d0
         der2 = 0.d0
         prop=prop*1.d-3
         der1 = der1*1.d-3
         der2 = der2*1.d-3
      else if(iflg.eq.1.and.iphase.eq.-22) then
         prop = 2802.6d0+3.063d0*(3.5d0-var1)
c     prop = 2799.8686959315d0+3.063d0*(3.5d0-var1)
         der1 = -3.063d0
         der2 = 0.d0
         prop=prop*1.d-3
         der1 = der1*1.d-3
         der2 = der2*1.d-3
      else if(iflg.eq.2.and.iphase.eq.1) then
c     solid density and derivative wrt pressure and temperature
         prop=denms + ddenmsp*(var1-pm0) + ddenmst*(var2-tm0)
         der1 = ddenmsp
         der2 = ddenmst
      else if(iflg.eq.2.and.iphase.eq.2) then
c     liquid density and derivative wrt pressure and temperature
c     liquid density
c     numerator coefficients
       if(itsat.le.10) then
         x=var1
         x2=x*x
         x3=x2*x
         x4=x3*x
         tl=var2
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
         dla0=  0.10000000d+01
         dlpa1=  0.17472599d-01
         dlpa2= -0.20443098d-04
         dlpa3= -0.17442012d-06
         dlta1=  0.49564109d-02
         dlta2= -0.40757664d-04
         dlta3=  0.50676664d-07
         dlpta=  0.50330978d-04  
         dlp2ta=  0.33914814d-06
         dlpt2a= -0.18383009d-06
c     denomenator coefficients        
         dlb0=  0.10009476d-02
         dlpb1=  0.16812589d-04
         dlpb2= -0.24582622d-07
         dlpb3= -0.17014984d-09
         dltb1=  0.48841156d-05
         dltb2= -0.32967985d-07
         dltb3=  0.28619380d-10
         dlptb=  0.53249055d-07
         dlp2tb=  0.30456698d-09
         dlpt2b= -0.12221899d-09
         rnwn1=dla0+dlpa1*x+dlpa2*x2+dlpa3*x3
         rnwn2=dlta1*tl+dlta2*tl2+dlta3*tl3
         rnwn3=dlpta*tlx+dlpt2a*tl2x+dlp2ta*tlx2
         rnwn=rnwn1+rnwn2+rnwn3
         rnwd1=dlb0+dlpb1*x+dlpb2*x2+dlpb3*x3
         rnwd2=dltb1*tl+dltb2*tl2+dltb3*tl3
         rnwd3=dlptb*tlx+dlpt2b*tl2x+dlp2tb*tlx2
         rnwd=rnwd1+rnwd2+rnwd3
         rnw=rnwn/rnwd
         prop=rnw
c     derivatives of density
         drlpn1=dlpa1+2*dlpa2*x+3*dlpa3*x2+dlpta*tl
         drlpn1=rnwd*(drlpn1+2*dlp2ta*tlx+dlpt2a*tl2)
         drlpn2=dlpb1+2*dlpb2*x+3*dlpb3*x2+dlptb*tl
         drlpn2=rnwn*(drlpn2+2*dlp2tb*tlx+dlpt2b*tl2)
         drlpn=drlpn1-drlpn2
         drolpd=rnwd**2
         der1=drlpn/drolpd
         drlen1=dlta1+2*dlta2*tl+3*dlta3*tl2+dlpta*x
         drlen1=rnwd*(drlen1+dlp2ta*x2+2*dlpt2a*tlx)
         drlen2=dltb1+2*dltb2*tl+3*dltb3*tl2+dlptb*x
         drlen2=rnwn*(drlen2+dlp2tb*x2+2*dlpt2b*tlx)
         drlen=drlen1-drlen2
         droled=rnwd**2
         der2=drlen/droled
       else if(itsat.gt.10) then
c   call  pseudo-vap eos   
          tl = var2
          x = var1     
          call eos_aux(itsat,tl,x,0,0,prop,der2,der1)
          continue
       endif
c     prop = 815.1d0+1.445d0*(240.d0-var2)-1.04d0*(5.d0-var1)
c     der1 = 1.04d0
c     der2 = -1.445d0
         
c     Below expression is from Haas JL (1976) 
         if (ibrine.ne.0) then
c     change water density (kg/m3) to specific volume in cm3/gm
c     change the salt conc. mol/kg from ppm, assume salt is mainly NaCl 
            if (cden_flag .eq. 2) then
c     Concentrations were summed to get total moles/kg-water
               mol = var3
            else
c     Convert ppm of salt to moles/kg-water    
               mol = var3/58.44d3
            end if
c            var3 = var3/58.44d3
            vo = 1.d3/prop
            vc = 3.1975d0
            phistar = -167.29d0+(448.55d0*vo)+(-261.07*vo*vo)
            kappa = (-13.644d0+(13.97*vo))*((vo/(vc-vo))**2.d0)
            phi = phistar + (kappa*dsqrt(mol))
            num1 = 1000.d0+(mol*58.445d0)
            den1 = (1000.d0*vo)+(mol*phi)
            prop = num1/den1

            dvodp = -der1*vo*vo*1.d-3
            dphistardp = (448.55d0-(2.d0*261.07d0*vo))*dvodp
            dkappadp = dvodp*(vo/((vc-vo)**2.d0))*(-(2.d0*13.644d0*vc)+
     &           (3.d0*13.97d0*vo*vc)-(13.97d0*vo*vo))
            dphidp = dphistardp + dsqrt(mol)*dkappadp
            dnumdp = 0.d0
            ddendp = 1000.d0*dvodp + mol*dphidp
            der1 = ((den1*dnumdp)-(num1*ddendp))/(den1*den1)

            dvodt = -der2*vo*vo*1.d-3
            dphistardt = (448.55d0-(2.d0*261.07d0*vo))*dvodt
            dkappadt = dvodt*(vo/((vc-vo)**2.d0))*(-(2.d0*13.644d0*vc)+
     &           (3.d0*13.97d0*vo*vc)-(13.97d0*vo*vo))
            dphidt = dphistardt + dsqrt(mol)*dkappadt
            dnumdt = 0.d0
            ddendt = 1000.d0*dvodt + var3*dphidt
            der2 = ((den1*dnumdt)-(num1*ddendt))/(den1*den1)

            prop = prop*1000.d0
            der1 = der1*1000.d0
            der2 = der2*1000.d0

c     convert the unit of salt conc. back to ppm.
c            var3 = var3*58.44d3
         endif
      else if(iflg.eq.2.and.iphase.eq.22) then

c     prop = 810.729360708d0+23.43088d0*(3.5d0-var1)
         prop = 809.74d0+23.43088d0*(3.5d0-var1)
         der1 = -23.43088d0
         der2 = 0.d0
      else if(iflg.eq.2.and.iphase.eq.-22) then
c     prop = 17.54606423d0-5.07889d0*(3.5d0-var1)
         prop = 17.526d0-5.07889d0*(3.5d0-var1)
         der1 = 5.07889d0
         der2 = 0.d0
      else if(iflg.eq.2.and.iphase.eq.3) then
c     vapor density and derivative wrt pressure and temperature
c     vapor density
         x=var1
         x2=x*x
         x3=x2*x
         x4=x3*x
         tl=var2
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
         dva0  =  0.15089524d-05
         dvpa1 =  0.10000000d+01
         dvpa2 = -0.10000000d+01
         dvpa3 = -0.16676705d-02
         dvta1 =  0.40111210d-07
         dvta2 =  0.25625316d-10
         dvta3 = -0.40479650d-12
         dvpta =  0.43379623d-01
         dvp2ta=  0.24991800d-02
         dvpt2a= -0.94755043d-04
         dvb0  =  0.12636224d+00
         dvpb1 = -0.30463489d+00
         dvpb2 =  0.27981880d-02
         dvpb3 =  0.51132337d-05
         dvtb1 =  0.59318010d-02
         dvtb2 =  0.80972509d-05
         dvtb3 = -0.43798358d-07
         dvptb =  0.53046787d-03
         dvp2tb= -0.84916607d-05
         dvpt2b=  0.48444919d-06
         rnsn1=dva0+dvpa1*x+dvpa2*x2+dvpa3*x3
         rnsn2=dvta1*tl+dvta2*tl2+dvta3*tl3
         rnsn3=dvpta*tlx+dvpt2a*tl2x+dvp2ta*tlx2
         rnsn=rnsn1+rnsn2+rnsn3
         rnsd1=dvb0+dvpb1*x+dvpb2*x2+dvpb3*x3
         rnsd2=dvtb1*tl+dvtb2*tl2+dvtb3*tl3
         rnsd3=dvptb*tlx+dvpt2b*tl2x+dvp2tb*tlx2
         rnsd=rnsd1+rnsd2+rnsd3
         rns=rnsn/rnsd
         prop=rns

c     derivatives of vapor density
         drspn1=dvpa1+2*dvpa2*x+3*dvpa3*x2+dvpta*tl
         drspn1=rnsd*(drspn1+2*dvp2ta*tlx+dvpt2a*tl2)
         drspn2=dvpb1+2*dvpb2*x+3*dvpb3*x2+dvptb*tl
         drspn2=rnsn*(drspn2+2*dvp2tb*tlx+dvpt2b*tl2)
         drspn=drspn1-drspn2
         drospd=rnsd**2
         der1=drspn/drospd
         drsen1=dvta1+2*dvta2*tl+3*dvta3*tl2+dvpta*x
         drsen1=rnsd*(drsen1+dvp2ta*x2+2*dvpt2a*tlx)
         drsen2=dvtb1+2*dvtb2*tl+3*dvtb3*tl2+dvptb*x
         drsen2=rnsn*(drsen2+dvp2tb*x2+2*dvpt2b*tlx)
         drsen=drsen1-drsen2
         drostd=rnsd**2
         der2=drsen/drostd

      else if(iflg.eq.3.and.iphase.eq.1) then
c     solid viscosity and derivative wrt pressure and temperature
c     solid is immobile
         prop=visms
         der1 = 0.0
         der2 = 0.0
      else if(iflg.eq.3.and.iphase.eq.2) then
c     liquid viscosity and derivative wrt pressure and temperature
       if(itsat.le.10) then
         vla0=  0.17409149d-02
         vlpa1=  0.18894882d-04
         vlpa2= -0.66439332d-07
         vlpa3= -0.23122388d-09
         vlta1= -0.31534914d-05
         vlta2=  0.11120716d-07
         vlta3= -0.48576020d-10
         vlpta=  0.28006861d-07
         vlp2ta=  0.23225035d-09
         vlpt2a=  0.47180171d-10
c     denomenator coefficients
         vlb0=  0.10000000d+01
         vlpb1=  0.10523153d-01
         vlpb2= -0.22658391d-05
         vlpb3= -0.31796607d-06
         vltb1=  0.29869141d-01
         vltb2=  0.21844248d-03
         vltb3= -0.87658855d-06
         vlptb=  0.41690362d-03
         vlp2tb= -0.25147022d-05
         vlpt2b=  0.22144660d-05
         x = var1
         x2 = x*x
         x3 = x2*x
         tl = var2
         tl2 = tl*tl
         tl3 = tl2*tl
         tlx = tl*x
         tl2x = tl2*x
         tlx2 = tl*x2
         viln1=vla0+vlpa1*x+vlpa2*x2+vlpa3*x3
         viln2=vlta1*tl+vlta2*tl2+vlta3*tl3
         viln3=vlpta*tlx+vlpt2a*tl2x+vlp2ta*tlx2
         viln=viln1+viln2+viln3
         vild1=vlb0+vlpb1*x+vlpb2*x2+vlpb3*x3
         vild2=vltb1*tl+vltb2*tl2+vltb3*tl3
         vild3=vlptb*tlx+vlpt2b*tl2x+vlp2tb*tlx2
         vild=vild1+vild2+vild3
         vil=viln/vild
         prop=vil
c     derivatives of liquid viscosity
         dvlpn1=vlpa1+2*vlpa2*x+3*vlpa3*x2+vlpta*tl
         dvlpn1=vild*(dvlpn1+2*vlp2ta*tlx+vlpt2a*tl2)
         dvlpn2=vlpb1+2*vlpb2*x+3*vlpb3*x2+vlptb*tl
         dvlpn2=viln*(dvlpn2+2*vlp2tb*tlx+vlpt2b*tl2)
         dvlpn=dvlpn1-dvlpn2
         dvilpd=vild**2
         der1=dvlpn/dvilpd
         dvlen1=vlta1+2*vlta2*tl+3*vlta3*tl2+vlpta*x
         dvlen1=vild*(dvlen1+vlp2ta*x2+2*vlpt2a*tlx)
         dvlen2=vltb1+2*vltb2*tl+3*vltb3*tl2+vlptb*x
         dvlen2=viln*(dvlen2+vlp2tb*x2+2*vlpt2b*tlx)
         dvlen=dvlen1-dvlen2
         dviled=vild**2
         der2=dvlen/dviled
       else if(itsat.gt.10) then
c   call  pseudo-vap eos  
          tl = var2
          x = var1         
          call eos_aux(itsat,tl,x,2,1,prop,der2,der1)
       endif               
         if (ibrine.ne.0) then
c     expression from Phillips et al. pg. 5. change salt conc. from
c     ppm to g moles/ kg H2O
            if (cden_flag .eq. 2) then
c     Sum concentrations to get total moles/kg-water
               mol = var3
            else
c     Convert ppm of salt to moles/kg-water    
               mol = var3/58.44d3
            end if
c            var3=var3/58.44d3
            prop1=prop
            prop2=1.d0+(0.0816d0*mol)+(0.012d0*mol*mol)+(0.000128d0*
     &           mol*mol*mol)+(0.000629d0*var2*(1.d0-
     &           dexp(-0.7d0*mol)))
            prop=prop1*prop2
            der1=0.d0
            der2=(der2*prop2)+(prop1*(0.000629d0*(1.d0-
     &           dexp(-0.7d0*mol))))
c            var3 = var3*58.44d3
         endif
      else if(iflg.eq.3.and.iphase.eq.3) then
c     vapor viscosity and derivative wrt pressure and temperature
c     vapor viscosity
         x = var1
         x2 = x*x
         x3 = x2*x
         tl = var2
         tl2 = tl*tl
         tl3 = tl2*tl
         tlx = tl*x
         tl2x = tl2*x
         tlx2 = tl*x2

         vva0= -0.13920783d-03
         vvpa1=  0.98434337d-02
         vvpa2= -0.51504232d-03
         vvpa3=  0.62554603d-04
         vvta1=  0.27105772d-04
         vvta2=  0.84981906d-05
         vvta3=  0.34539757d-07
         vvpta= -0.25524682d-03
         vvp2ta= 0.d0
         vvpt2a=  0.12316788d-05
         vvb0=  0.10000000d+01
         vvpb1=  0.10000000d+01
         vvpb2= -0.10000000d+01
         vvpb3= -0.10000000d+01
         vvtb1=  0.10000000d+01
         vvtb2=  0.10000000d+01
         vvtb3= -0.22934622d-03
         vvptb=  0.10000000d+01
         vvp2tb= 0.d0
         vvpt2b=0.25834551d-01

         visn1=vva0+vvpa1*x+vvpa2*x2+vvpa3*x3
         visn2=vvta1*tl+vvta2*tl2+vvta3*tl3
         visn3=vvpta*tlx+vvpt2a*tl2x+vvp2ta*tlx2
         visn=visn1+visn2+visn3
         visd1=vvb0+vvpb1*x+vvpb2*x2+vvpb3*x3
         visd2=vvtb1*tl+vvtb2*tl2+vvtb3*tl3
         visd3=vvptb*tlx+vvpt2b*tl2x+vvp2tb*tlx2
         visd=visd1+visd2+visd3
         vis=visn/visd
         prop=vis

c     derivatives of vapor viscosity
         dvspn1=vvpa1+2*vvpa2*x+3*vvpa3*x2+vvpta*tl
         dvspn1=visd*(dvspn1+2*vvp2ta*tlx+vvpt2a*tl2)
         dvspn2=vvpb1+2*vvpb2*x+3*vvpb3*x2+vvptb*tl
         dvspn2=visn*(dvspn2+2*vvp2tb*tlx+vvpt2b*tl2)
         dvspn=dvspn1-dvspn2
         dvispd=visd**2
         der1=dvspn/dvispd
         dvsen1=vvta1+2*vvta2*tl+3*vvta3*tl2+vvpta*x
         dvsen1=visd*(dvsen1+vvp2ta*x2+2*vvpt2a*tlx)
         dvsen2=vvtb1+2*vvtb2*tl+3*vvtb3*tl2+vvptb*x
         dvsen2=visn*(dvsen2+vvp2tb*x2+2*vvpt2b*tlx)
         dvsen=dvsen1-dvsen2
         dvised=visd**2
         der2=dvsen/dvised
      else if(iflg.eq.4.and.iphase.eq.1) then
c     phase-change(s-l) temperature and derivative wrt pressure
         prop = tphasesl + dtphaseslp*(var1-pphasesl)
         der1 = dtphaseslp
         der2 = 0.0
      else if(iflg.eq.4.and.iphase.eq.2) then
c     phase-change(l-v) temperature and derivative wrt pressure
c     prop = tphaselv + dtphaselvp*(var1-pphaselv)
c     der1 = dtphaselvp
c     der2 = 0.0
c     changed by RJP 4/12/04
c     introduced new equation valid from 0.01 deg C to 300 deg C	
c     the equation is a curve fit generated from the data that was
c     calculated using Wagner & Prub formulation given by NIST
c     J. Phy. Chem. Ref. Data 2002, Vol 31, No. 2
c     x = log(var1)
c     x2 = x*x
c     x3 = x2*x
c     x4 = x3*x
c     prop = 	0.0191d0*x4+0.4225d0*x3+4.6982d0*x2+43.623d0*x+179.8
c     der1 = (4.d0*0.0191d0*x3+3.d0*0.4225d0*x2+2.d0*4.6982d0*x+43.623d0)
c     &	/var1
c     der2 = 0.d0
c     RJP 12/23/04 remove below to above

         tsa0   = -0.25048121d-05
         tspa1  =  0.45249584d-02
         tspa2  =  0.33551528d+00
         tspa3  =  0.10000000d+01
         tspa4  =  0.12254786d+00
         tsb0   =  0.20889841d-06
         tspb1  =  0.11587544d-03
         tspb2  =  0.31934455d-02
         tspb3  =  0.45538151d-02
         tspb4  =  0.23756593d-03
         x=var1
         x2 = x*x
         x3 = x2*x
         x4 = x3*x
         tfunn=tsa0+tspa1*x+tspa2*x2+tspa3*x3+tspa4*x4
         tfund=tsb0+tspb1*x+tspb2*x2+tspb3*x3+tspb4*x4
         tfun=tfunn/tfund
         prop=tfun
         dtpsn=((tspa1+2.*tspa2*x+3.*tspa3*x2+4.*tspa4*x3)*tfund)-
     &        (tfunn*(tspb1+2.*tspb2*x+3.*tspb3*x2+4.*tspb4*x3))
         dtpsd=tfund**2
         der1=dtpsn/dtpsd
         der2=0.d0
         if (ibrine.ne.0) then
c     use Haas's formula.  Iterate on temperature.
            guess=prop
            y = var3
            a1 = 5.93582d-6
            a2 = -5.19386d-5
            a3 = 1.23156d-5
            b1 = 1.1542d-6
            b2 = 1.41254d-7
            b3 = -1.92476d-8
            b4 = -1.70717d-9
            b5 = 1.0539d-10
            asl = 1.d0+a1*y+a2*y*y+a3*y*y*y
            bsl = b1*y+b2*y*y+b3*y*y*y+b4*y*y*y*y+b5*y*y*y*y*y
            do i = 1, 1000
               f=dlog(guess)-(bsl*guess*dlog(prop))-(asl*dlog(prop))
               df = (1.d0/guess)-(bsl*dlog(prop))
               dx = f/df
               guess=guess-dx
               if (dabs(dx).lt.1d-6) goto 101
            enddo
 101        continue
            sl = asl+bsl*guess
            der1 = (sl*sl*guess)*der1/(prop*(sl-bsl*dlog(guess)))
            prop=guess
            der2 = 0.d0
         endif
      else if(iflg.eq.5.and.iphase.eq.1) then
c     phase-change(s-l) pressure and derivative wrt temperature
         prop = pphasesl + 1.0/dtphaseslp*(var2-tphasesl)
         der1 = 0.0
         der2 = 1.0/dtphaseslp
      else if(iflg.eq.5.and.iphase.eq.2) then
c     phase-change(l-v) pressure and derivative wrt temperature
c     prop = pphaselv + 1.0/dtphaselvp*(var2-tphaselv)
c     der1 = 0.0
c     der2 = 1.0/dtphaselvp
c     changed by RJP 4/12/04
c     introduced new equation valid from 0.01 deg C to 300 deg C	
c     the equation is a curve fit generated from the data that was
c     calculated using Wagner & Prub formulation given by NIST
c     J. Phy. Chem. Ref. Data 2002, Vol 31, No. 2
         x = var2
         x2 = x*x
         x3 = x2*x
         x4 = x3*x
         prop=(-2.d-10*x4)+(2.d-7*x3)-(1.d-4*x2)+(0.0307d0*x)-3.2095d0
         der1 = 0.d0
         der2 = (-(8.d-10*x3)+(6.d-7*x2)-(2.d-5*x)+0.0307d0)
         prop=10.d0**prop
         der2=prop*der2
         if(ibrine.ne.0) then
            y = var3
            a1 = 5.93582d-6
            a2 = -5.19386d-5
            a3 = 1.23156d-5
            b1 = 1.1542d-6
            b2 = 1.41254d-7
            b3 = -1.92476d-8
            b4 = -1.70717d-9
            b5 = 1.0539d-10
            asl = 1.d0+a1*y+a2*y*y+a3*y*y*y
            bsl = b1*y+b2*y*y+b3*y*y*y+b4*y*y*y*y+b5*y*y*y*y*y
            t1 = dexp(dlog(var2)/(asl+bsl*var2))
            x = t1
            x2 = x*x
            x3 = x2*x
            x4 = x3*x
            prop=(1.d-9*x4)-(2.d-7*x3)+(1.d-5*x2)-(3.d-4*x)+0.0031
            der1 = 0.d0
            sl = asl+bsl*t1
            der2 = t1*(((sl/var2)-(dlog(var2)*bsl))/(sl*sl))
         endif
      else if(iflg.eq.6.and.iphase.eq.1) then
c     liquid-vapor capillary pressure and derivative wrt liq. saturation
         prop = 0.0       
         der1 = 0.0
         der2 = 0.0
      else if(iflg.eq.6.and.iphase.eq.2) then
c     solid-liquid capillary pressure and derivative wrt liq. saturation
         prop = 0.0       
         der1 = 0.0
         der2 = 0.0
      else if(iflg.eq.7.and.iphase.eq.1) then
c     vapor relative permeability and derivative wrt liq. saturation
         prop = 1.0-var3      
         der1 = 0.0
         der2 = -1.0
      else if(iflg.eq.7.and.iphase.eq.2) then
c     liquid relative permeability and derivative wrt liq. saturation
         prop = var3      
         der1 = 0.0
         der2 = 1.0
      else if(iflg.eq.8.and.iphase.eq.1) then
c     solid relative permeability and derivative wrt liq. saturation
         prop = 0.0       
         der1 = 0.0
         der2 = 0.0
      else if(iflg.eq.8.and.iphase.eq.2) then
c     liquid relative permeability and derivative wrt liq. saturation
         prop = var3      
         der1 = 0.0
         der2 = 1.0
      else if(iflg.eq.9) then
c     source enthalpy returned when a temperature is specified
         x=var1
         x2=x*x
         x3=x2*x
         x4=x3*x
         tl=var2
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
         ela0=0.25623465d-03
         elpa1=0.10184405d-02
         elpa2=0.22554970d-04
         elpa3=0.34836663d-07
         elta1=0.41769866d-02
         elta2=-0.21244879d-04
         elta3=0.25493516d-07
         elpta=0.89557885d-04
         elp2ta=0.10855046d-06
         elpt2a=-0.21720560d-06
c     denomenator coefficients
         elb0=0.10000000d+01
         elpb1=0.23513278d-01
         elpb2=0.48716386d-04
         elpb3=-0.19935046d-08
         eltb1=-0.50770309d-02
         eltb2=0.57780287d-05
         eltb3=0.90972916d-09
         elptb=-0.58981537d-04
         elp2tb=-0.12990752d-07
         elpt2b=0.45872518d-08
c     liquid enthalpy
         enwn1=ela0+elpa1*x+elpa2*x2+elpa3*x3
         enwn2=elta1*tl+elta2*tl2+elta3*tl3
         enwn3=elpta*tlx+elpt2a*tl2x+elp2ta*tlx2
         enwn=enwn1+enwn2+enwn3
         enwd1=elb0+elpb1*x+elpb2*x2+elpb3*x3
         enwd2=eltb1*tl+eltb2*tl2+eltb3*tl3
         enwd3=elptb*tlx+elpt2b*tl2x+elp2tb*tlx2
         enwd=enwd1+enwd2+enwd3
         enw=enwn/enwd
         prop=enw

c     derivatives of enthalpy
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
         der2=dhwt
         der1=dhwp
c     expression from Michaelides
         if (ibrine.ne.0) then
            a(0,0) = -9633.6d0
            a(0,1) = -4080.d0
            a(0,2) = 286.49d0
            a(1,0) = 166.58d0
            a(1,1) = 68.577d0
            a(1,2) = -4.6856d0
            a(2,0) = -0.90963d0
            a(2,1) = -0.36524d0
            a(2,2) = 0.249667d-1
            a(3,0) = 0.17965d-2
            a(3,1) = 0.71924d-3
            a(3,2) = -0.49d-4
            if (cden_flag .eq. 2) then
c     Sum concentrations to get total moles/kg-water
               mol = var3
            else
c     Convert ppm of salt to moles/kg-water    
               mol = var3/58.44d3
            end if
c     'mol' in below equation is the molality
            del_h=0.d0
            derdel_h=0.d0
            do i = 0, 3
               do j = 0, 2
                  del_h=del_h+(a(i,j)*(var2**i)*(mol**j))
               enddo	
            enddo	
            do i = 1, 3
               do j = 0, 2
                  derdel_h=derdel_h+(a(i,j)*i*(var2**(i-1))*(mol**j))
               enddo
            enddo
            del_h=-1d-3*(4.184d0/(1000.d0+(58.44d0*mol)))*del_h
            derdel_h=1d-3*(4.184d0/(1000.d0+(58.44d0*mol)))*derdel_h
            x1=1000.d0/(1000.d0+58.44*mol)
            x2=58.44*mol/(1000.d0+58.44*mol)
            h_salt=41.293d0*(var2+273.15d0)+
     &           (3.3607d-2*((var2+273.15d0)**2.d0)/2.d0)
     &           -(1.3927d-5*((var2+273.15d0)**3.d0)/3.d0)
            dh_salt=41.293d0+
     &           (3.3607d-2*(var2+273.15d0))
     &           -(1.3927d-5*((var2+273.15d0)**2.d0))
c     convert units from J/Mole to MJ/Kg
            h_salt=h_salt*1d-3/58.44d0
            dh_salt=dh_salt*1d-3/58.44d0
            prop=x1*prop+x2*h_salt+mol*del_h
            der1=der1
            der2=x1*der2+x2*dh_salt+mol*derdel_h
         endif
      else if(iflg.eq.10.and.iphase.eq.1) then
c     liquid-vapor capillary pressure and derivative wrt water fraction
         prop = 0.0                   
         der1 = 0.0       
         der2 = 0.0
c     updated by Rajesh to take care of negative saturations 11/22/2003
c     prop = cap_max*(1.0-var1)
c     der1 = -cap_max
c     der2 = 0.0
         prop = 1.0-var1
         if (prop.le.0.0) then
            prop = 0.0
            der1 = 0.0
            der2 = 0.0
         else
            prop = cap_max*prop
            der1 = -cap_max
            der2 = 0.0
         endif
      else if(iflg.eq.11) then
c     water-methane relative permeability based on fractions 
c     for liquid water
c     indepedent of fluid state 
         if(var3.lt.fracw_max) then
            prop = var3                    
            der1 = 0.0        
            der2 = 0.0
            der3 = 1.0
         else
            prop = fracw_max               
            der1 = 0.0        
            der2 = 0.0
            der3 = 0.0
         endif
c     Modified by Rajesh 11/25/03 
c     if(var3.lt.fracw_max) then
c     prop = var3**4                    
c     der1 = 0.0        
c     der2 = 0.0
c     der3 = 4*var3**3
c     else
c     prop = fracw_max               
c     der1 = 0.0        
c     der2 = 0.0
c     der3 = 0.0
c     endif

      else if(iflg.eq.12) then
c     water-methane production factor 
c     for liquid water
c     indepedent of fluid state 
         if(var3.gt.fracw_min) then
            prop = (1.0-var3-var4)/var3             
            der1 = 0.0        
            der3 = (-1.0)/var3-(1.0-var3-var4)/(var3*var3)
            der4 = -1.0/var3
         else
            prop = 0.0                              
            der1 = 0.0        
            der3 = 0.0  
            der4 = 0.0 
         endif
      else if(iflg.eq.13) then
c     water-CO2 production factor 
c     for liquid water
c     independent of fluid state 
         if(var3.gt.fracw_min) then
            prop = (1.0-var3)/var3             
            der1 = 0.0        
            der3 = (-1.0)/(var3*var3)
            der4 = 1.0/(var3*var3)
         else
            prop = 0.0                              
            der1 = 0.0        
            der3 = 0.0  
            der4 = 0.0 
         endif

      endif
      return 
      end

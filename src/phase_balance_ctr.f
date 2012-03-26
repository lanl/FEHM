       subroutine phase_balance_ctr(iflg,i,phi0,tl0,sl0)
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
!!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 insure energy balance during phase changes
!D1
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 29-Nov-09, Programmer: George Zyvoloski
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 3.0
!D4 
!**********************************************************************

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comci
      use combi
      use comdti
      use comki  
      use comai

      implicit none

      integer, allocatable :: kq_dum(:), icount(:)
      integer i,j,ii,ij,jj,kb,iflg,iieosd
      integer max_iter
      character*4 vard
      real*8 pl,tl,sl,phi0,tl0, eps_phase
      real*8 energy_0,energy,d_energy,sl0,resid,delp
      parameter (max_iter = 10, eps_phase = 1.d-6)
      
c      if(iphasebal.eq.0) return
      if(iflg.eq.0) then
c
c input (no input yet)
c
            
      else if(iflg.eq.1) then
c 
c  iteration control (liquid to 2-phase )
c   
        pl = phi0
        tl = tl0
        sl = sl0
        iieosd = iieos(i)
        call energy_calc(pl,tl,sl,1,iieosd,energy_0,d_energy)
               sl = sl0
               do  ii=1,max_iter
                  call energy_calc(pl,tl,sl,2,iieosd,energy,d_energy)
                  resid=energy-energy_0
                  delp=-resid/d_energy
                  sl = sl + delp
                  if (abs(resid).le.eps_phase) go to 1000
               enddo
c check for  convergence
 1000          continue
               if(ii .gt. max_iter) then
                  if(iout.ne.0) 
     &                 write(iout,*) 'failed in phase_balance_ctr'
                  if(iptty.ne.0)
     *                 write(iptty,*) 'failed in phase_balance_ctr'
               end if
       
      endif   
      return 
      end
      subroutine energy_calc(pl,tl,sl,ieosd,iieosd,energy,daeh)
c      
c     c     adjust coefficients for thermo fits
c
      use comrxni
      use comii
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer ndummy,iieosl,mid,mi,ieosd,iieosd,kq
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
      real*8 drlpn,drolpd,drolp,drlen1,drlen2,drlen,droled,drolt,viln1
      real*8 viln2,viln3,viln,vild1,vild2,vild3,vild,vil,xvisl
      real*8 dvlpn1,dvlpn2,dvlpn,dvilpd,dvislp,dvlen1,dvlen2,dvlen
      real*8 dviled,dvislt,ensn1,ensn2,ensn3,ensn,ensd1,ensd2,ensd3
      real*8 ensd,ens,env,dhvp1,dhvp2,dhvpn,dhvpd,dhvt1,dhvt2
      real*8 dhvtn,dhvtd,dhvt,dhvp,rnsn1,rnsn2,rnsn3,rnsd1,rnsd2
      real*8 rnsd3,rnsn,rnsd,rns,rov,drspn1,drspn2,drspn,drospd
      real*8 drsen1,drsen2,drsen,drostd,visn1,visn2,visn3,visn,visd1
      real*8 visd2,visd3,visd,vis,xvisv,dvspn1,dvspn2,dvspn,dvispd
      real*8 dvisvp,dvsen1,dvsen2,dvsen,dvised,dvisvt,sl,qdis
      real*8 cp,por,vol,pldif,permsd,eskd,eqdum,rag,sig
      real*8 rop,rop1,daep1,dragp,dsigp,dsige
      real*8 den1,roe,dql,hprod,dhprdp,dhprde,htc,tbound,hflux
      real*8 dhflxp,sbound,cprd,edif,denrd,vfd,rl,dvfp,tfun
      real*8 tfunn,tfund,dtpsn,dtpsd,dpldt,psat,vfcal,rop2,daep2
      real*8 dqv,dhflxe,drovp,drovt
      real*8 xa,xa2,xa3,xa4
      real*8 p_energy
      real*8 dtps,dtd,damp,daep,damh,daeh,energy,d_energy
      real*8 cden_correction, cden_cor
      p_energy = 0.0

c     liquid phase coefficients
c     
c     liquid enthalpy
c     numerator coefficients
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
c     denomenator coefficients
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
c     liquid density
c     numerator coefficients
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
c     denomenator coefficients
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
c     liquid viscosity
c     numerator coefficients
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
c     denomenator coefficients
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
c     vapor phase coefficients
c     
c     vapor enthalpy
c     numerator coefficients
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
c     denomenator coefficients
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
c     vapor density
c     numerator coefficients
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
c     denomenator coefficients
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
c     vapor viscosity
c     numerator coefficients
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
c     denomenator coefficients
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


c     evaluate thermo functions and derivatives
         x=pl - phi_inc
         x2=x*x
         x3=x2*x
         x4=x3*x
         xa=pl
         xa2=xa*xa
         xa3=xa2*xa
         xa4=xa3*xa 
         if(ieosd.eq.2) then

c     two phase conditions

c     calculate temperature and dt/dp
            tfunn=tsa0+tspa1*xa+tspa2*xa2+tspa3*xa3+tspa4*xa4
            tfund=tsb0+tspb1*xa+tspb2*xa2+tspb3*xa3+tspb4*xa4
            tfun=tfunn/tfund
            tl=tfun
            dtpsn=((tspa1+2.*tspa2*xa+3.*tspa3*xa2+4.*tspa4*xa3)*tfund)-
     &           (tfunn*(tspb1+2.*tspb2*xa+3.*tspb3*xa2+4.*tspb4*xa3))
            dtpsd=tfund**2
            dtps=dtpsn/dtpsd
         else
c            tl=tl0
         endif
         tl2=tl*tl
         tl3=tl2*tl
         tlx=x*tl
         tl2x=tl2*x
         tlx2=tl*x2
         if(ieosd.ne.3) then

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
            enl=enw + p_energy

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
            dhlt=dhwt
            dhlp=dhwp

c     liquid density
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
            if(cden) then
c     Add correction for liquid species
               cden_cor = cden_correction(mi)
               rol = rol + cden_cor
            end if


c     derivatives of density
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

c     liquid viscosity
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

c     derivatives of liquid viscosity
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

         endif

         if(ieosd.ne.1) then

c     vapor enthalpy
            ensn1=eva0+evpa1*x+evpa2*x2+evpa3*x3
            ensn2=evta1*tl+evta2*tl2+evta3*tl3
            ensn3=evpta*tlx+evpt2a*tl2x+evp2ta*tlx2
            ensn=ensn1+ensn2+ensn3
            ensd1=evb0+evpb1*x+evpb2*x2+evpb3*x3
            ensd2=evtb1*tl+evtb2*tl2+evtb3*tl3
            ensd3=evptb*tlx+evpt2b*tl2x+evp2tb*tlx2
            ensd=ensd1+ensd2+ensd3
            ens=ensn/ensd
            env=ens + p_energy

c     derivatives of vapor enthalpy
            dhvp1=evpa1+2*evpa2*x+3*evpa3*x2+evpta*tl
            dhvp1=ensd*(dhvp1+2*evp2ta*tlx+evpt2a*tl2)
            dhvp2=evpb1+2*evpb2*x+3*evpb3*x2+evptb*tl
            dhvp2=ensn*(dhvp2+2*evp2tb*tlx+evpt2b*tl2)
            dhvpn=dhvp1-dhvp2
            dhvpd=ensd**2
            dhvp=dhvpn/dhvpd
            dhvt1=evta1+2*evta2*tl+3*evta3*tl2+evpta*x
            dhvt1=ensd*(dhvt1+evp2ta*x2+2*evpt2a*tlx)
            dhvt2=evtb1+2*evtb2*tl+3*evtb3*tl2+evptb*x
            dhvt2=ensn*(dhvt2+evp2tb*x2+2*evpt2b*tlx)
            dhvtn=dhvt1-dhvt2
            dhvtd=ensd**2
            dhvt=dhvtn/dhvtd
            
c     vapor density
            rnsn1=dva0+dvpa1*x+dvpa2*x2+dvpa3*x3
            rnsn2=dvta1*tl+dvta2*tl2+dvta3*tl3
            rnsn3=dvpta*tlx+dvpt2a*tl2x+dvp2ta*tlx2
            rnsn=rnsn1+rnsn2+rnsn3
            rnsd1=dvb0+dvpb1*x+dvpb2*x2+dvpb3*x3
            rnsd2=dvtb1*tl+dvtb2*tl2+dvtb3*tl3
            rnsd3=dvptb*tlx+dvpt2b*tl2x+dvp2tb*tlx2
            rnsd=rnsd1+rnsd2+rnsd3
            rns=rnsn/rnsd
            rov=rns

c     derivatives of vapor density
            drspn1=dvpa1+2*dvpa2*x+3*dvpa3*x2+dvpta*tl
            drspn1=rnsd*(drspn1+2*dvp2ta*tlx+dvpt2a*tl2)
            drspn2=dvpb1+2*dvpb2*x+3*dvpb3*x2+dvptb*tl
            drspn2=rnsn*(drspn2+2*dvp2tb*tlx+dvpt2b*tl2)
            drspn=drspn1-drspn2
            drospd=rnsd**2
            drovp=drspn/drospd
            drsen1=dvta1+2*dvta2*tl+3*dvta3*tl2+dvpta*x
            drsen1=rnsd*(drsen1+dvp2ta*x2+2*dvpt2a*tlx)
            drsen2=dvtb1+2*dvtb2*tl+3*dvtb3*tl2+dvptb*x
            drsen2=rnsn*(drsen2+dvp2tb*x2+2*dvpt2b*tlx)
            drsen=drsen1-drsen2
            drostd=rnsd**2
            drovt=drsen/drostd

c     vapor viscosity
            visn1=vva0+vvpa1*x+vvpa2*x2+vvpa3*x3
            visn2=vvta1*tl+vvta2*tl2+vvta3*tl3
            visn3=vvpta*tlx+vvpt2a*tl2x+vvp2ta*tlx2
            visn=visn1+visn2+visn3
            visd1=vvb0+vvpb1*x+vvpb2*x2+vvpb3*x3
            visd2=vvtb1*tl+vvtb2*tl2+vvtb3*tl3
            visd3=vvptb*tlx+vvpt2b*tl2x+vvp2tb*tlx2
            visd=visd1+visd2+visd3
            vis=visn/visd
            xvisv=vis

c     derivatives of vapor viscosity
            dvspn1=vvpa1+2*vvpa2*x+3*vvpa3*x2+vvpta*tl
            dvspn1=visd*(dvspn1+2*vvp2ta*tlx+vvpt2a*tl2)
            dvspn2=vvpb1+2*vvpb2*x+3*vvpb3*x2+vvptb*tl
            dvspn2=visn*(dvspn2+2*vvp2tb*tlx+vvpt2b*tl2)
            dvspn=dvspn1-dvspn2
            dvispd=visd**2
            dvisvp=dvspn/dvispd
            dvsen1=vvta1+2*vvta2*tl+3*vvta3*tl2+vvpta*x
            dvsen1=visd*(dvsen1+vvp2ta*x2+2*vvpt2a*tlx)
            dvsen2=vvtb1+2*vvtb2*tl+3*vvtb3*tl2+vvptb*x
            dvsen2=visn*(dvsen2+vvp2tb*x2+2*vvpt2b*tlx)
            dvsen=dvsen1-dvsen2
            dvised=visd**2
            dvisvt=dvsen/dvised
         endif

c     modify derivatives for 2-phase
         if(ieosd.eq.2) then
            drolp=drolp+drolt*dtps
            drovp=drovp+drovt*dtps
            dhlp=dhlp+dhlt*dtps
            dhvp=dhvp+dhvt*dtps
            dvisvp=dvisvp+dvisvt*dtps
            dvislp=dvislp+dvislt*dtps
            drolt=0.0
            drovt=0.0
            dhlt=0.0
            dhvt=0.0
            dvisvt=0.0
            dvislt=0.0
            dtd=0.0
         endif
         
c     compressed liquid
         if(ieosd.eq.1) then
            dtps=0.0
            sl=1.d0
            sv=0.0
            sig=0.d0
            rov=0.d0
            env=0.d0
            drovp=0.0
            drovt=0.0
            dhvp=0.0
            dhvt=0.0
            dvisvp=0.0
            dvisvt=0.0
            dtd=1.0
            xvisv=1.d20

c     accumulation terms
            den1=rol
            den=den1*por
            dene=(1.d0-por)*cp*tl+den*enl-por*pl

c     derivatives of accumulation terms
            rop=por*drolp
            rop2=rol*dporpl
            rop1=rop
            damp=(rop1+rop2)*dtin
            daep1=(rop*enl+por*rol*dhlp-por)
            daep2=dporpl*(-cp*tl+rol*enl-pl)
            daep=(daep1+daep2)*dtin
            roe=por*drolt
            damh=(rol*dportl+drolt*por)*dtin
            daeh=(-dportl*cp*tl+(1.-por)*cp+
     &           dportl*(enl*rol-pl)+por*(rol*dhlt+enl*drolt))*dtin

         end if

c     superheated vapour
         if(ieosd.eq.3) then
            dtps=0.0
            sl=0.d0
            sv=1.d0
            sig=1.d0
            rol=0.d0
            enl=0.d0
            drolp=0.0
            drolt=0.0
            dhlp=0.0
            dhlt=0.0
            dvislp=0.0
            dvislt=0.0
            dtd=1.0
            xvisl=1.d20

c     accumulation terms
            den1=rov
            den=den1*por
            dene=(1.d0-por)*cp*tl+den*env-por*pl

c     derivatives of accumulation terms
            rop=por*drovp
            damp =rop*dtin
            daep=(env*rop-por)*dtin
            roe=por*drovt
            damh =roe*dtin
            daeh =((1.d0-por)*cp+por*rov*dhvt+env*roe)*dtin
         end if         
c     two phase conditions
         if(ieosd.eq.2) then
c     saturation and relative permeabilities
            sv=1.d0-sl
c     accumulation terms
            den=por*(sl*rol+sv*rov)
            eqdum=sl*rol*enl+sv*rov*env-pl
            dene=((1.-por)*cp*tl+por*eqdum)
c     derivatives of accumulation terms
            rop=por*(sv*drovp+sl*drolp)
            damp =rop*dtin
            daep =((1.d0-por)*cp*dtps+por*(sv*drovp*env+sv*rov*dhvp+
     &           sl*drolp*enl+sl*rol*dhlp)-por)*dtin
            damh =por*(rol-rov)*dtin
            daeh =por*(rol*enl-rov*env)*dtin
         end if      
       energy = dene
      end
      

      subroutine ther_meth_h2o(iflg,ndummy)
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
!D1 To calculate the equation coeffients and derivatives for a
!D1 h2o-methane simulation.
!D1 solid-liquid-gas system.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20 
!D2 
!D2 Initial implementation: Date 01-June-02, Programmer: George Zyvoloski
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
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
      use commeth
      implicit none

      integer ndummy,iieosl,mid,mi,ieosd,kq
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
      real*8 dvisvp,dvsen1,dvsen2,dvsen,dvised,dvisvt,dtd,sl,qdis
      real*8 cp,por,vol,pldif,permsd,eskd,eqdum,rag,sig,dtrdp,dtrdt
      real*8 rop,rop1,damp,daep,daep1,damh,daeh,dragp,dsigp,dsige
      real*8 den1,roe,dql,hprod,dhprdp,dhprde,htc,tbound,hflux
      real*8 dhflxp,sbound,cprd,edif,denrd,vfd,rl,dvfp,tfun
      real*8 tfunn,tfund,dtpsn,dtpsd,dpldt,dpld3,psat,vfcal,rop2,daep2
      real*8 dqv,dhflxe,dtps,drovp,drovt
      real*8 dhflxem,dhflxpm
      real*8, allocatable :: sto1(:)


      real*8 psatd
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 frac_meth,frac_w,dfrac_w,dfrac_methm
      real*8 damm,daem,dsigm,drvm,drlm,dhprdm                        
      real*8 damw,daew,dsigw,drvw,drlw,dhprdw                        
      real*8 dhflxew,frac_g,dfrac_gm,dfrac_gw                         
      real*8 frac_hyd,dfrac_hyd                                       
      real*8 dskhydp,dskhydt,dskhydw,dskhydm           
      real*8 dqhhydp,dqhhydt,dqhhydw,dqhhydm 
      real*8 fracmb,methterm
      real*8 cden_correction, cden_cor
c     
c     gaz took it out of loop (can be specified source as well)
c     will over-ride above methane production
c     
      
      real*8 rlwm,drlwm3,drlwm4,drlwmp,drlwmt
      real*8 rgwm,drgwm3,drgwm4,drgwmp,drgwmt
      real*8 prf,dprf4                                                
      real*8 prodfac,dprodfac3,dprodfac4,fracwmin

      parameter (fracwmin=0.1)
      integer idum,iflg

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
      dtin=1.0/dtot

c     ****************************************************************
      if(iflg.eq.0) then
c     calculations that are common to methane,methane hydrate, and water
c     calculate pressure dependant porosity and derivatives

         if(iporos.ne.0) call porosi(1)
c     new cappilary pressure routine
         call cappr_hyd(1,0)

         allocate(rl_meth(n0))
         allocate(drl_meth3(n0))
         allocate(drl_meth4(n0))
         allocate(rl_h2o(n0))
         allocate(drl_h2o3(n0))
         allocate(drl_h2o4(n0))
         rl_h2o = 0.0
         rl_meth = 0.0
         
         drl_meth3 = 0.0
         drl_meth4 = 0.0
         drl_h2o3 = 0.0
         drl_h2o4 = 0.0

         do mid=1,neq
            mi=mid+ndummy
            call rlperm_hyd(ndummy,0,mi,rl_h2o(mi),0.0,0.0,
     &           drl_h2o3(mi),drl_h2o4(mi),rl_meth(mi),
     &           0.0,0.0,drl_meth3(mi),drl_meth4(mi))

         enddo 


c     
c     ****************************************************************
      else if(iflg.eq.-1) then
c     release memory

         deallocate(rl_meth)
         deallocate(drl_meth3)
         deallocate(drl_meth4)
         deallocate(rl_h2o)
         deallocate(drl_h2o3)
         deallocate(drl_h2o4)

c     ****************************************************************
      else if(iflg.eq.1) then
c     calculations for water    

         iieosl=0
         drvw = 0.0
         drlw = 0.0
         drvm = 0.0
         do mid=1,neq
            mi=mid+ndummy
c     phase information now contained in array ices
            ieosd=ieos(mi)

            pl=phi(mi)
            sl=s(mi)

c     
c     term to account for space taken up with water   
c     
            frac_w = fracw(mi)
            frac_hyd  = frachyd(mi)
            dfrac_w = 1.0
            dfrac_hyd = -1.0
            if(ieosd.eq.2) then
               
c     two phase (gas/liquid) conditions
c     calculate relative permeabilities of vapor-liquid

               call h2o_properties(6,1,dum1,dum2,sl,dum3,
     &              pcp(mi),dum4,dpcef(mi),dum5,dum6)

c     calculate relative permeabilities of vapor-liquid

               call h2o_properties(7,1,dum1,dum2,sl,dum3,
     &              rvf(mi),drvpf(mi),drvef(mi),dum5,dum6)
               xrv=rvf(mi) 
               drv=drvef(mi)
               drvp=drvpf(mi)

               call h2o_properties(7,2,dum1,dum2,sl,dum3,
     &              rlf(mi),drlpf(mi),drlef(mi),dum5,dum6)
               xrl=rlf(mi) 
               drl=drlef(mi)
               drlp=drlpf(mi)
            endif

            if(ieosd.eq.-2) then
               
c     two phase (liquid/solid) conditions
c     calculate relative permeabilities of vapor-liquid
c     rvf now containes the solid relative permeability

               call h2o_properties(6,2,dum1,dum2,sl,dum3,
     &              pcp(mi),dum4,dpcef(mi),dum5,dum6)

c     calculate relative permeabilities of vapor-liquid

               call h2o_properties(8,1,dum1,dum2,sl,dum3,
     &              rvf(mi),drvpf(mi),drvef(mi),dum5,dum6)
               xrv=rvf(mi) 
               drv=drvef(mi)
               drvp=drvpf(mi)

               call h2o_properties(8,2,dum1,dum2,sl,dum3,
     &              rlf(mi),drlpf(mi),drlef(mi),dum5,dum6)
               xrl=rlf(mi) 
               drl=drlef(mi)
               drlp=drlpf(mi)
            endif

            if(ieosd.eq.2) then

c     two phase (gas/liquid) conditions

c     calculate phase-change temperature and dt/dp
               call h2o_properties(4,2,pl,dum1,dum2,dum3,
     &              tl,dtps,dum4,dum5,dum6)
c     calculate phase-change pressure and dp/dt
               call h2o_properties(5,2,tl,dum1,dum2,dum3,
     &              psatd,dum4,dpsatt,dum5,dum6)
c     
            else if(ieosd.eq.-2) then

c     two phase (liquid/solid) conditions

c     calculate temperature and dt/dp
               call h2o_properties(4,1,pl,dum1,dum2,dum3,
     &              tl,dtps,dum4,dum5,dum6)
c     calculate pressure  and dp/dt
               call h2o_properties(5,1,tl,dum1,dum2,dum3,
     &              psatd,dum4,dpsatt,dum5,dum6)
            else
               tl=t(mi)
            endif
            if(abs(ieosd).ne.3) then

c     liquid enthalpy and derivatives

               call h2o_properties(1,2,pl,tl,dum2,dum3,
     &              enl,dhlp,dhlt,dum5,dum6)

c     liquid density and derivatives

               call h2o_properties(2,2,pl,tl,dum2,dum3,
     &              rol,drolp,drolt,dum5,dum6)

c     modify density(explicit update) if it changes with concentration
               if(cden) then
c     Add correction for liquid species
                  cden_cor = cden_correction(mi)
                  rol = rol + cden_cor
               end if

c     liquid viscosity and derivatives

               call h2o_properties(3,2,pl,tl,dum2,dum3,
     &              xvisl,dvislp,dvislt,dum5,dum6)

            endif

            if(ieosd.ge.2) then

c     vapor enthalpy and derivatives

               call h2o_properties(1,3,pl,tl,dum2,dum3,
     &              env,dhvp,dhvt,dum5,dum6)

c     vapor density and derivatives

               call h2o_properties(2,3,pl,tl,dum2,dum3,
     &              rov,drovp,drovt,dum5,dum6)

c     vapor viscosity and derivatives

               call h2o_properties(3,3,pl,tl,dum2,dum3,
     &              xvisv,dvisvp,dvisvt,dum5,dum6)

            endif
            if(ieosd.le.-2) then

c     properties of solid phase
c     properties of solid phase will be put into variables for vapor
c     it is assumed that solid and vapor can not exist simultaneously at the same node

c     calculate solid enthalpy and derivative wrt p and t

               call h2o_properties(1,1,pl,tl,dum2,dum3,
     &              env,dhvp,dhvt,dum5,dum6)

c     calculate solid density and derivative wrt p and t

               call h2o_properties(2,1,pl,tl,dum2,dum3,
     &              rov,drovp,drovt,dum5,dum6)

c     calculate solid viscosity (should be very large!) and derivative wrt p and t
c     should be zero

               call h2o_properties(3,1,pl,tl,dum2,dum3,
     &              xvisv,dvisvp,dvisvt,dum5,dum6)

            endif
c     calculate special capillary pressure as a function of water frac
c     gaz 1-12-04 implemented cappr_hyd
c     call h2o_properties(10,1,frac_w,tl,dum2,dum3,
c     &                              pcpmeth(mi),dpcefmeth(mi),
c     &                              dum4,dum5,dum6)
c     gaz 12-24-03
c     calculate special relative permeability for liquid water-methane gas
c     
c     call h2o_properties(11,1,pl,tl,frac_w,dum2,    
c     &                              rlwm,drlwmp,drlwmt,drlwm3,dum6)
            rlwm = rl_h2o(mi)
            drlwm3 = drl_h2o3(mi)
            drlwm4 = drl_h2o4(mi)
            drlwmp = 0.0
            drlwmt = 0.0

c     calculate solid gas generation term and derivative wrt p and t and frachyd
c     kinetic approach (note call to hydrate properties)
c     porh is pass-through variable ps(mi) to hydrate_properties

            porh = ps(mi)
            volh = sx1(mi)

            call hydrate_properties(3,1,pl,tl,fracw(mi),frachyd(mi),
     &           skhyd(mi),dskhydp,dskhydt,dskhydw,dskhydm)
c     
c     save rate and derivatives for enthalp calculation (passed through commeth)
c     

            source_hyd_temp = skhyd(mi)
            ds1_tmp = dskhydp
            ds2_tmp = dskhydt
            ds3_tmp = dskhydw
            ds4_tmp = dskhydm

c     calculate solid gas heat generation term and derivative wrt p and t and frachyd
c     kinetic approach
c     porh is pass-through variable ps(mi) to hydrate_properties

            porh = ps(mi)
            volh = sx1(mi)

            call hydrate_properties(4,1,pl,tl,fracw(mi),frachyd(mi),
     &           qhhyd(mi),dqhhydp,dqhhydt,dqhhydw,dqhhydm)

c     form dissociation flow term

            skwhyd(mi) = -ratw*skhyd(mi)

c     load water generation term derivatives 

	    dskhyd1(mi) = -ratw*dskhydp
	    dskhyd2(mi) = -ratw*dskhydt
	    dskhyd3(mi) = -ratw*dskhydw
	    dskhyd4(mi) = -ratw*dskhydm

c     modify derivatives for 2-phase (should now depend on p, not t)
            if(abs(ieosd).eq.2) then
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
               dtd=1.0
            endif

            dq(mi)=0.0
            qh(mi)=0.0
            dqh(mi)=0.0
            deqh(mi)=0.0
            deqw(mi)=0.0
            deqpw(mi)=0.0
            dqt(mi)=0.0
            qdis=sk(mi)
            cp=denr(mi)*cpr(mi)
            por=ps(mi)
            dporpl=dporp(mi)
            dportl=dport(mi)
            kq=ka(mi)
            vol=volume(mi)

            if(kq.lt.0) then
c     form pressure dependent flow term
               qh(mi)=0.
               sk(mi)=0.0
               if(pflow(mi).gt.0.0) then
                  pldif=pl-pflow(mi)
                  dpldt=0.0
               else
                  pldif=pl-(psatd-pflow(mi))
                  dpldt=-dpsatt
               endif
               permsd=wellim(mi)
               if(pldif.le.0.0d00.and.kq.eq.-2) permsd=0.
c     changed by rajesh to account for water fraction
c     qdis=permsd*(pldif)
c     dq(mi)=permsd
c     sk(mi)=qdis
c     
               qdis=permsd*rlwm*pldif
               sk(mi)=qdis
               dq(mi)=permsd*rlwm+drlwmp*permsd*pldif
               dq3(mi)=permsd*(pldif*drlwm3)
               dq4(mi)=permsd*(pldif*drlwm4)
               if(qdis.le.0.) then
                  eskd=eflow(mi)
                  qh(mi)=eskd*qdis
                  sk(mi)=qdis
                  dqh(mi)=eskd*dq(mi)
               endif
               dqt(mi)=permsd*(dpldt*rlwm+pldif*drlwmt)
            endif
c     
c     two phase conditions
            if(abs(ieosd).eq.2) then
c     saturation and relative permeabilities
               sv=1.d0-sl
c     accumulation terms
c     note new terms for fraction of water
               den=por*(sl*rol+sv*rov)*frac_w
               eqdum=(sl*rol*enl+sv*rov*env-pl)*frac_w
               dene=((1.-por)*cp*tl+por*eqdum)
c     production of steam
c     note xrv and xrl now also depend on frac_w 
c     xrv and xrl dependence on frac_w not implemented yet
               rag=rol*xvisv/rov/xvisl
               sig=xrv/(xrv+rag*xrl)
c     derivatives of accumulation terms
               rop=por*(sv*drovp+sl*drolp)*frac_w
               damp =rop*dtin
               daep =((1.d0-por)*cp*dtps+por*frac_w*(sv*drovp*env+
     &              sv*rov*dhvp+sl*drolp*enl+sl*rol*dhlp)-por)*dtin
               damh =por*frac_w*(rol-rov)*dtin
               daeh =por*frac_w*(rol*enl-rov*env)*dtin
c     add derivative wrt methane fraction
               damw=por*(sl*rol+sv*rov)*dfrac_w
               daew=por*(sl*rol*enl+sv*rov*env-pl)*dfrac_w
               damm=por*(sl*rol+sv*rov)*dfrac_hyd
               daem=por*(sl*rol*enl+sv*rov*env-pl)*dfrac_hyd
c     derivatives of sink terms
               if(qdis.gt.0.0) then
                  dragp=(rov*xvisl*(drolp*xvisv+rol*dvisvp)-rol*
     2                 xvisv*(drovp
     3                 *xvisl+rov*dvislp))/(xvisl*rov)**2
                  if(sig.ne.0.d0.and.sig.ne.1.d0) then
                     dsigp=-xrv*dragp*xrl/(xrv+rag*xrl)**2
                     dsige=((xrv+rag*xrl)*drv-xrv*(drv+rag*drl))
     3                    /(xrv+rag*xrl)**2
                     dsigw=((xrv+rag*xrl)*drvw-xrv*(drvm+rag*drlw))
     3                    /(xrv+rag*xrl)**2
                  else
                     dsigp=0.d0
                     dsige=0.d0
                     dsigw=0.d0
                  end if
               end if
            end if

c     compressed liquid
            if(ieosd.eq.1) then
               xrl=1.0
               xrv=0.0
               drl=0.0
               drv=0.0
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
               den=rol*por*frac_w
               dene=(1.d0-por)*cp*tl+por*(rol*enl-pl)*frac_w

c     derivatives of accumulation terms
               damp=(rol*dporpl+por*drolp)*dtin*frac_w
               damw=rol*por*dtin*dfrac_w
c     damm=rol*por*dtin*dfrac_hyd
               damm=0.0                     
               daep1=(dporpl*rol*enl+por*drolp*enl+por*rol*dhlp
     &              -dporpl*pl-por)*frac_w
c     gaz corrected (rajesh found) 12-17-03
c     daep2=-dporpl*cp*tl*frac_w
               daep2=-dporpl*cp*tl
               daep=(daep1+daep2)*dtin
               damh=(dportl*rol+por*drolt)*dtin*frac_w
               daeh=(-dportl*cp*tl+(1.-por)*cp+(dportl*(enl*rol-pl)+
     &              por*(rol*dhlt+enl*drolt))*frac_w)*dtin
               daew=por*(rol*enl-pl)*dtin*dfrac_w
c     daem=por*(rol*enl-pl)*dtin*dfrac_hyd
               daem=0.0                              
            end if

c     superheated vapour
            if(abs(ieosd).eq.3) then
               dtps=0.0
               xrl=0.0
               xrv=1.0
               drl=0.0
               drv=0.0
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
               den=rov*por*frac_w
               dene=(1.d0-por)*cp*tl+por*(rov*env-pl)*frac_w

c     derivatives of accumulation terms
               damp=(dporpl*rov+por*drovp)*dtin*frac_w
               damw=rov*por*dtin*dfrac_w
               damm=rov*por*dtin*dfrac_hyd
               daep1=(dporpl*rov*env+por*drovp*env+por*rov*dhvp
     &              -dporpl*pl-por)*frac_w
               daep2=-dporpl*cp*tl
               daep=(daep1+daep2)*dtin
               damh=(dportl*rov+por*drovt)*dtin*frac_w
               daeh=(-dportl*cp*tl+(1.-por)*cp+(dportl*(env*rov-pl)+
     &              por*(rov*dhvt+env*drovt))*frac_w)*dtin
               daew=por*(rov*env-pl)*dtin*dfrac_w
               daem=por*(rov*env-pl)*dtin*dfrac_hyd
            end if
            if(ps(mi).ne.0.0) then

c     store derivatives of accumulation terms
               dmef(mi)=damh
               dmpf(mi)=damp
               dmmf(mi)=damm
               dmwf(mi)=damw
               depf(mi)=daep
               deef(mi)=daeh
               dewf(mi)=daew
               demf(mi)=daem
               dtpa(mi)=dtps
               dtpae(mi)=dtd
               t(mi)=tl
               s(mi)=sl

c     save accumulation terms for possible volume changes
               sto1(mi)=den
               sto1(mi+neq)=dene
               dstm(mi)=por*rov*sv*vol
               dil(mi)=0.0
               dilp(mi)=0.0
               dile(mi)=0.0
               dilw(mi)=0.0
               dilm(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divw(mi)=0.0
               divm(mi)=0.0

c     permeability reduction factor

               prf = permhyd(mi)
               dprf4 = dpermhyd4(mi)

               if(abs(ieosd).ne.3) then

c     modify flow terms for new upwind scheme
c     note rlwm used for both phases of methane 
                  dql=prf*rlwm*rol*xrl/xvisl
                  dil(mi)=dql
                  dilp(mi)=prf*rlwm*(drolp*xrl/xvisl-
     &                 rol*xrl/xvisl**2*dvislp+rol*drlp/xvisl)
                  dile(mi)=prf*rlwm*(drolt*xrl/xvisl+rol*drl/xvisl
     &                 -rol*xrl/xvisl**2*dvislt)
                  dilw(mi)=prf*rlwm*rol*drlw/xvisl+
     &                 prf*drlwm3*rol*xrl/xvisl
                  dilm(mi)=dprf4*rol*xrl/xvisl
               endif
               dglp(mi)=drolp
               dgle(mi)=drolt
               enlf(mi)=enl
               delf(mi)=dhlp
               delef(mi)=dhlt

c     modify flow terms for new upwind scheme
               if(ieosd.ne.1) then
                  dqv=prf*rlwm*rov*xrv/xvisv
                  div(mi)=dqv
                  divp(mi)=prf*rlwm*(drovp*xrv/xvisv-
     &                 rov*xrv/xvisv**2*dvisvp+rov*drvp/xvisv)
                  dive(mi)=prf*rlwm*(drovt*xrv/xvisv+rov*drv/xvisv
     &                 -rov*xrv/xvisv**2*dvisvt)
                  divw(mi)=prf*rlwm*rov*drvw/xvisv+
     &                 prf*drlwm3*rov*xrv/xvisv
                  divm(mi)=dprf4*rov*xrv/xvisv
               endif
               dgvp(mi)=drovp
               dgve(mi)=drovt
               envf(mi)=env
               devf(mi)=dhvp
               devef(mi)=dhvt
               rolf(mi)=rol
               rovf(mi)=rov

               if(qdis.gt.0.0) then

c     organize source terms and derivatives
                  if(abs(ieosd).eq.2) then
                     hprod=sig*env+(1.0-sig)*enl
                     dhprdp=dsigp*env+sig*dhvp-dsigp*enl+(1.0-sig)*dhlp
                     dhprde=dsige*(env-enl)
                     dhprdw=dsigw*(env-enl)
                  endif
                  if(ieosd.eq.1) then
                     hprod=enl
                     dhprdp=dhlp
                     dhprde=dhlt
                  endif
                  if(abs(ieosd).eq.3) then
                     hprod=env
                     dhprdp=dhvp
                     dhprde=dhvt
                  endif
                  qh(mi)=hprod*qdis
                  dqh(mi)=dhprdp*qdis+hprod*dq(mi)
                  deqh(mi)=dhprde*qdis+hprod*dqt(mi)
                  deqw(mi)=dhprdw*qdis+hprod*dq3(mi)
               endif
               if(qdis.lt.0.) then
                  qh(mi)=qdis*eflow(mi)
                  dqh(mi)=dq(mi)*eflow(mi)
               endif
            endif
c     
c     add intercomponent heat flux
c     this term transfers heat from water to methane (and back)
c     also note this term only contains derivative wrt water   
c     derivative with respect to methane can also be derived 
c     
            if(qhflxmeth(mi).ne.0.0) then
               htc = qhflxmeth(mi)
               hflux=htc*(t(mi)-tmeth(mi))
               if(abs(ieosd).ne.2) then
                  dhflxp=0.0
                  dhflxe=htc
               else
                  dhflxp=htc*dtpa(mi)
                  dhflxe=0.0
               endif
c     this next if block contains derivative wrt methane component
               if(abs(ices(mi)).ne.2) then
                  dhflxew=-htc
               else
                  dhflxpm=-htc*dtpameth(mi)
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
                  if(abs(ieosd).ne.2) then
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
                  if(abs(ieosd).ne.2) then
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
                  edif=cprd*tl-eskd
                  permsd=wellim(mi)
                  qh(mi)=permsd*edif
                  deqh(mi)=permsd*cprd
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
            deni(mi)=(sto1(mi)  -denh(mi))*dtin
            denei(mi)=(sto1(mi+neq)-deneh(mi))*dtin
         enddo    

c     ****************************************************************
      else if(iflg.eq.2) then
c     
c     calculations for fluid methane
c     
c     note: heat capacity of rock thermal term removed
c     now included only with water conservation of energy equation
c     
         drvw = 0.0
         drlw = 0.0
         drlm = 0.0
         drvm = 0.0
         iieosl=0
         do mid=1,neq
            mi=mid+ndummy
c     phase information now contained in array ices
            ieosd=ices(mi)

            pl=phimeth(mi)
            sl=smeth(mi)
c     water and methane are independent
            frac_hyd  = frachyd(mi)
            frac_w = fracw(mi)
            frac_g = 1.0-fracw(mi)-frachyd(mi)
            dfrac_w = 1.0d00
            dfrac_gw = -1.0d00
            dfrac_gm = -1.0d00

c     
            if(ieosd.eq.2) then
               
c     two phase (gas/liquid) conditions
c     calculate relative permeabilities of vapor-liquid

               call methane_properties(6,1,dum1,dum2,sl,dum3,
     &              pcp(mi),dum4,dpcef(mi),dum5,dum6)

c     calculate relative permeabilities of vapor-liquid

               call methane_properties(7,1,dum1,dum2,sl,dum3,
     &              rvf(mi),drvpf(mi),drvef(mi),dum5,dum6)
               xrv=rvf(mi) 
               drv=drvef(mi)
               drvp=drvpf(mi)

               call methane_properties(7,2,dum1,dum2,sl,dum3,
     &              rlf(mi),drlpf(mi),drlef(mi),dum5,dum6)
               xrl=rlf(mi) 
               drl=drlef(mi)
               drlp=drlpf(mi)
            endif

            if(ieosd.eq.-2) then
               
c     two phase (liquid/solid) conditions
c     calculate relative permeabilities of vapor-liquid
c     rvf now containes the solid relative permeability

               call methane_properties(6,2,dum1,dum2,sl,dum3,
     &              pcp(mi),dum4,dpcef(mi),dum5,dum6)

c     calculate relative permeabilities of vapor-liquid

               call methane_properties(8,1,dum1,dum2,sl,dum3,
     &              rvf(mi),drvpf(mi),drvef(mi),dum5,dum6)
               xrv=rvf(mi) 
               drv=drvef(mi)
               drvp=drvpf(mi)

               call methane_properties(8,2,dum1,dum2,sl,dum3,
     &              rlf(mi),drlpf(mi),drlef(mi),dum5,dum6)
               xrl=rlf(mi) 
               drl=drlef(mi)
               drlp=drlpf(mi)
            endif

            if(ieosd.eq.2) then

c     two phase (gas/liquid) conditions

c     calculate phase-change temperature and dt/dp
               call methane_properties(4,2,pl,dum1,dum2,dum3,
     &              tl,dtps,dum4,dum5,dum6)
c     calculate phase-change pressure and dp/dt
               call methane_properties(5,2,tl,dum1,dum2,dum3,
     &              psatd,dum4,dpsatt,dum5,dum6)
c     
            else if(ieosd.eq.-2) then

c     two phase (liquid/solid) conditions

c     calculate temperature and dt/dp
               call methane_properties(4,1,pl,dum1,dum2,dum3,
     &              tl,dtps,dum4,dum5,dum6)
c     calculate pressure  and dp/dt
               call methane_properties(5,1,tl,dum1,dum2,dum3,
     &              psatd,dum4,dpsatt,dum5,dum6)
            else
               tl=tmeth(mi)
            endif
            if(abs(ieosd).ne.3) then

c     liquid enthalpy and derivatives

               call methane_properties(1,2,pl,tl,dum2,dum3,
     &              enl,dhlp,dhlt,dum5,dum6)

c     liquid density and derivatives

               call methane_properties(2,2,pl,tl,dum2,dum3,
     &              rol,drolp,drolt,dum5,dum6)

c     modify density(explicit update) if it changes with concentration
               if(cden) then
c     Add correction for liquid species
                  cden_cor = cden_correction(mi)
                  rol = rol + cden_cor
               end if

c     liquid viscosity and derivatives

               call methane_properties(3,2,pl,tl,dum2,dum3,
     &              xvisl,dvislp,dvislt,dum5,dum6)

            endif

            if(ieosd.ge.2) then

c     vapor enthalpy and derivatives

               call methane_properties(1,3,pl,tl,dum2,dum3,
     &              env,dhvp,dhvt,dum5,dum6)

c     vapor density and derivatives

               call methane_properties(2,3,pl,tl,dum2,dum3,
     &              rov,drovp,drovt,dum5,dum6)

c     vapor viscosity and derivatives

               call methane_properties(3,3,pl,tl,dum2,dum3,
     &              xvisv,dvisvp,dvisvt,dum5,dum6)

            endif
            if(ieosd.le.-2) then

c     properties of solid phase
c     properties of solid phase will be put into variables for vapor
c     it is assumed that solid and vapor can not exist simultaneously at the same node

c     calculate solid enthalpy and derivative wrt p and t

               call methane_properties(1,1,pl,tl,dum2,dum3,
     &              env,dhvp,dhvt,dum5,dum6)

c     calculate solid density and derivative wrt p and t

               call methane_properties(2,1,pl,tl,dum2,dum3,
     &              rov,drovp,drovt,dum5,dum6)

c     calculate solid viscosity (should be very large!) and derivative wrt p and t

               call methane_properties(3,1,pl,tl,dum2,dum3,
     &              xvisv,dvisvp,dvisvt,dum5,dum6)

            endif

c     calculate special relative permeability for liquid water-methane gas

c     gaz 12-24-03
c     call methane_properties(11,1,pl,tl,frac_w,frac_hyd,
c     &                              rgwm,drgwmp,drgwmt,drgwm3,drgwm4)
            rgwm = rl_meth(mi)
            drgwmp = 0.0
            drgwmt = 0.0
            drgwm3 = drl_meth3(mi)
            drgwm4 = drl_meth4(mi)

c     added liquid relative permeability rajesh 11/20/2003
            call h2o_properties(11,1,pl,tl,frac_w,dum2,    
     &           rlwm,drlwmp,drlwmt,drlwm3,dum6)
c     porh is pass-through variable ps(mi) to hydrate_properties

            porh = ps(mi)
            volh = sx1(mi)

            call hydrate_properties(3,1,pl,tl,fracw(mi),frachyd(mi),
     &           skhyd(mi),dskhydp,dskhydt,dskhydw,dskhydm)

c     
c     save rate and derivatives for enthalp calculation (passed through commeth)
c     

            source_hyd_temp = skhyd(mi)
            ds1_tmp = dskhydp
            ds2_tmp = dskhydt
            ds3_tmp = dskhydw
            ds4_tmp = dskhydm

c     calculate solid gas heat generation term and derivative wrt p and t and frachyd
c     kinetic approach
c     porh is pass-through variable ps(mi) to hydrate_properties

            porh = ps(mi)
            volh = sx1(mi)

            call hydrate_properties(4,1,pl,tl,fracw(mi),frachyd(mi),
     &           qhhyd(mi),dqhhydp,dqhhydt,dqhhydw,dqhhydm)

c     form dissociation flow term

            skmhyd(mi) = -ratgas*skhyd(mi)

c     load gas generation term derivatives 

	    dskhyd1(mi) = -ratgas*dskhydp
	    dskhyd2(mi) = -ratgas*dskhydt
	    dskhyd3(mi) = -ratgas*dskhydw
	    dskhyd4(mi) = -ratgas*dskhydm

c     modify derivatives for 2-phase (should now depend on p, not t)
            if(abs(ieosd).eq.2) then
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
               dtd=1.0
            endif

c     
c     gaz need to carry dq from water calcs
c     dq(mi)=0.0
            qhmeth(mi)=0.0
            dqh(mi)=0.0
            deqh(mi)=0.0
            dqt(mi)=0.0
            kq=kameth(mi)
            if(kq.ne.1) then
               skmeth(mi) =0.0
               dskmethw1(mi)= 0.0
               dskmethw3(mi)= 0.0
               dskmethw4(mi)= 0.0
            else
               dskmethw1(mi)= 0.0
               dskmethw3(mi)= 0.0
               dskmethw4(mi)= 0.0
            endif

            por=ps(mi) 
            dporpl=dporp(mi)
            dportl=dport(mi)

            vol=volume(mi)

c     gaz skmethw should be removed (plus derivatives)
c     enough storage already  (deleted 1-12-04 gaz)                           
c     skmethw(mi) =0.0
            if(sk(mi).gt.0.0.and.kq.eq.0) then
c     calculate production factor for liquid water-methane gas
               call h2o_properties(12,1,pl,tl,frac_w,frac_hyd,
     &              prodfac,dum3,dum4,
     &              dprodfac3,dprodfac4)
               skmeth(mi)=sk(mi)*prodfac
               dskmethw1(mi)=dq(mi)*prodfac
               dskmethw3(mi)=sk(mi)*dprodfac3+dq3(mi)*prodfac
               dskmethw4(mi)=sk(mi)*dprodfac4   
               dq(mi)=0.0
            else if(kq.lt.0) then
c     form pressure dependent flow term
c     gaz 2-19-04 added gas relative perms
               qhmeth(mi)=0.
               skmeth(mi)=0.0
c     tenma 6-6-2005
c     note this forces the methane pressure to equal the water pressure at wells
c     ie pressure equilibrium
c     
               pflowmeth(mi) = pflow(mi)
               wellmeth(mi) = wellim(mi)
               if(pflowmeth(mi).gt.0.0) then
                  pldif=pl-pflowmeth(mi)
                  dpldt=0.0
                  permsd=wellmeth(mi)
                  if(pldif.le.0.0d00.and.kq.eq.-2) permsd=0.
                  qdis=permsd*rgwm*(pldif)
                  skmeth(mi)=qdis
                  dskmethw1(mi) =permsd*rgwm
                  dskmethw3(mi) =permsd*drgwm3*pldif
                  dskmethw4(mi) =permsd*drgwm4*pldif
               else
c     this case for specified methane fraction
                  fracmb=-pflowmeth(mi)
                  pldif=frac_g-fracmb                     
                  permsd=wellmeth(mi)
                  if(pldif.le.0.0d00.and.kq.eq.-2) permsd=0.
                  qdis=permsd*(pldif)
                  skmeth(mi)=qdis
                  dq(mi) =0.0
                  dskmethw1(mi)=0.0                      
                  dskmethw3(mi)=permsd*dfrac_gw   
                  dskmethw4(mi)=permsd*dfrac_gm   
               endif
               if(qdis.le.0.) then
                  eskd=eflowmeth(mi)
                  qhmeth(mi)=eskd*qdis
                  skmeth(mi)=qdis
                  dqh(mi)=eskd*dq(mi)
               endif
            endif
c     identify flux term
            qdis=skmeth(mi)


c     two phase conditions
            if(abs(ieosd).eq.2) then
c     saturation and relative permeabilities
               sv=1.d0-sl
c     accumulation terms
c     note new terms for fraction of water
               den=por*(sl*rol+sv*rov)*frac_g
               eqdum=(sl*rol*enl+sv*rov*env-pl)*frac_g
c     dene=((1.-por)*cp*tl+por*eqdum)
               dene=por*eqdum
c     production of steam
c     note xrv and xrl now also depend on frac_g (and thus frac_g)
               rag=rol*xvisv/rov/xvisl
               sig=xrv/(xrv+rag*xrl)
c     derivatives of accumulation terms
               rop=por*(sv*drovp+sl*drolp)*frac_g
               damp =rop*dtin
c     daep =((1.d0-por)*cp*dtps+por*frac_g*(sv*drovp*env+
               daep =(por*frac_g*(sv*drovp*env+
     &              sv*rov*dhvp+sl*drolp*enl+sl*rol*dhlp)-por)*dtin
               damh =por*frac_g*(rol-rov)*dtin
               daeh =por*frac_g*(rol*enl-rov*env)*dtin
c     add derivative wrt methane fraction
               damm=por*(sl*rol+sv*rov)*dfrac_gm
               daem=por*(sl*rol*enl+sv*rov*env-pl)*dfrac_gm
c     add derivative wrt water fraction
               damw=por*(sl*rol+sv*rov)*dfrac_gw
               daew=por*(sl*rol*enl+sv*rov*env-pl)*dfrac_gw
c     derivatives of sink terms
               if(qdis.gt.0.0) then
                  dragp=(rov*xvisl*(drolp*xvisv+rol*dvisvp)-rol*
     2                 xvisv*(drovp
     3                 *xvisl+rov*dvislp))/(xvisl*rov)**2
                  if(sig.ne.0.d0.and.sig.ne.1.d0) then
                     dsigp=-xrv*dragp*xrl/(xrv+rag*xrl)**2
                     dsige=((xrv+rag*xrl)*drv-xrv*(drv+rag*drl))
     3                    /(xrv+rag*xrl)**2
                     dsigm=(-xrv*(drvm+rag*drlm))
     3                    /(xrv+rag*xrl)**2
                     dsigw=((xrv+rag*xrl)*drvw-xrv*(rag*drlw))
     3                    /(xrv+rag*xrl)**2
                  else
                     dsigp=0.d0
                     dsige=0.d0
                     dsigm=0.d0
                     dsigw=0.d0
                  end if
               end if
            end if

c     compressed liquid
            if(ieosd.eq.1) then
               xrl=1.0
               xrv=0.0
               drl=0.0
               drv=0.0
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
               den=rol*por*frac_g
c     dene=(1.d0-por)*cp*tl+por*(rol*enl-pl)*frac_g
               dene=por*(rol*enl-pl)*frac_g
c     Date 11/25/03 RJP added to compensate for -ve frac_g remove later
c     dene=dene+2.
c     derivatives of accumulation terms
               damp=(rol*dporpl+por*drolp)*dtin*frac_g
               damm=rol*por*dtin*dfrac_gm
               damw=rol*por*dtin*dfrac_gw
               daep1=(dporpl*rol*enl+por*drolp*enl+por*rol*dhlp
     &              -dporpl*pl-por)*frac_g
c     daep2=-dporpl*cp*tl
               daep2=0.0d00             
               daep=(daep1+daep2)*dtin
               damh=(dportl*rol+por*drolt)*dtin*frac_g
c     daeh=(-dportl*cp*tl+(1.-por)*cp+(dportl*(enl*rol-pl)+
               daeh=((dportl*(enl*rol-pl)+
     &              por*(rol*dhlt+enl*drolt))*frac_g)*dtin
               daem=por*(rol*enl-pl)*dtin*dfrac_gm
               daew=por*(rol*enl-pl)*dtin*dfrac_gw

            end if

c     superheated vapour
            if(abs(ieosd).eq.3) then
               dtps=0.0
               xrl=0.0
               xrv=1.0
               drl=0.0
               drv=0.0
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
               den=rov*por*frac_g
c     dene=(1.d0-por)*cp*tl+por*(rov*env-pl)*frac_g
               dene=por*(rov*env-pl)*frac_g

c     derivatives of accumulation terms
               damp=(dporpl*rov+por*drovp)*dtin*frac_g
               damm=rov*por*dtin*dfrac_gm
               damw=rov*por*dtin*dfrac_gw
               daep1=(dporpl*rov*env+por*drovp*env+por*rov*dhvp
     &              -dporpl*pl-por)*frac_g
c     daep2=-dporpl*cp*tl
               daep2=0.0d00
               daep=(daep1+daep2)*dtin
               damh=(dportl*rov+por*drovt)*dtin*frac_g
c     daeh=(-dportl*cp*tl+(1.-por)*cp+(dportl*(env*rov-pl)+
               daeh=((dportl*(env*rov-pl)+
     &              por*(rov*dhvt+env*drovt))*frac_g)*dtin
               daem=por*(rov*env-pl)*dtin*dfrac_gm
               daew=por*(rov*env-pl)*dtin*dfrac_gw

            end if
            if(ps(mi).ne.0.0) then

c     store derivatives of accumulation terms
               dmef(mi)=damh
               dmpf(mi)=damp
               dmmf(mi)=damm
               dmwf(mi)=damw
               depf(mi)=daep
               deef(mi)=daeh
               demf(mi)=daem
               dewf(mi)=daew
               dtpameth(mi)=dtps
               dtpaemeth(mi)=dtd
               tmeth(mi)=tl
               smeth(mi)=sl

c     save accumulation terms for possible volume changes
               sto1(mi)=den
               sto1(mi+neq)=dene
               dstm(mi)=por*rov*sv*vol
               dil(mi)=0.0
               dilp(mi)=0.0
               dile(mi)=0.0
               dilm(mi)=0.0
               dilw(mi)=0.0
               div(mi)=0.0
               divp(mi)=0.0
               dive(mi)=0.0
               divm(mi)=0.0
               divw(mi)=0.0

c     permeability reduction factor

               prf = permhyd(mi)
               dprf4 = dpermhyd4(mi)

               if(abs(ieosd).ne.3) then

c     modify flow terms for new upwind scheme
c     note rgwm used for both phases of methane 
                  dql=prf*rgwm*rol*xrl/xvisl
                  dil(mi)=dql
                  dilp(mi)=prf*rgwm*(drolp*xrl/xvisl-
     &                 rol*xrl/xvisl**2*dvislp+rol*drlp/xvisl)
                  dile(mi)=prf*rgwm*(drolt*xrl/xvisl+rol*drl/xvisl
     &                 -rol*xrl/xvisl**2*dvislt)
                  dilw(mi)=prf*rgwm*rol*drlw/xvisl+
     &                 prf*drgwm3*rol*xrl/xvisl
                  dilm(mi)=dprf4*rgwm*rol*xrl/xvisl+
     &                 prf*drgwm4*rol*xrl/xvisl
               endif
               dglp(mi)=drolp
               dgle(mi)=drolt
               enlfmeth(mi)=enl
               delf(mi)=dhlp
               delef(mi)=dhlt

c     modify flow terms for new upwind scheme
               if(ieosd.ne.1) then
                  dqv=prf*rgwm*rov*xrv/xvisv
                  div(mi)=dqv
                  divp(mi)=prf*rgwm*(drovp*xrv/xvisv-
     &                 rov*xrv/xvisv**2*dvisvp+rov*drvp/xvisv)
                  dive(mi)=prf*rgwm*(drovt*xrv/xvisv+rov*drv/xvisv
     &                 -rov*xrv/xvisv**2*dvisvt)
                  divw(mi)=prf*rgwm*rov*drvw/xvisv+
     &                 prf*drgwm3*rov*xrv/xvisv
c     divm(mi)=dprf4*rov*xrv/xvisv
c     divw(mi)=prf*rgwm*rov*drvw/xvisv+prf*drgwm3*rov*xrv/xvisv
                  divm(mi)=dprf4*rgwm*rov*xrv/xvisv+
     &                 prf*drgwm4*rov*xrv/xvisv
               endif
               dgvp(mi)=drovp
               dgve(mi)=drovt
               envfmeth(mi)=env
               devf(mi)=dhvp
               devef(mi)=dhvt
               rolfmeth(mi)=rol
               rovfmeth(mi)=rov

               if(qdis.ge.0.0) then
                  
c     organize source terms and derivatives
                  if(abs(ieosd).eq.2) then
                     hprod=sig*env+(1.0-sig)*enl
                     dhprdp=dsigp*env+sig*dhvp-dsigp*enl+(1.0-sig)*dhlp
                     dhprde=dsige*(env-enl)
                     dhprdm=dsigm*(env-enl)
                     dhprdw=dsigw*(env-enl)
                  endif
                  if(ieosd.eq.1) then
                     hprod=enl
                     dhprdp=dhlp
                     dhprde=dhlt
                     dhprdm=0.d0
                     dhprdw=0.d0
                  endif
                  if(abs(ieosd).eq.3) then
                     hprod=env
                     dhprdp=dhvp
                     dhprde=dhvt
                     dhprdm=0.d0
                     dhprdw=0.d0
                  endif
c     gaz 1-12-04 removed skmethw 
                  qhmeth(mi)=hprod*(skmeth(mi))
                  dqmh(mi)=dhprdp*(skmeth(mi))
     &                 +hprod*dskmethw1(mi)
c     gaz 1-12-04 gaz
c     &                      +hprod*dq(mi)
                  deqmh(mi)=dhprde*(skmeth(mi))
                  deqw(mi)=dhprdw*(skmeth(mi))
     &                 +hprod*dskmethw3(mi)
                  deqm(mi)=dhprdm*(skmeth(mi))
     &                 +hprod*dskmethw4(mi)
               endif
               if(qdis.lt.0.) then
                  qhmeth(mi)=qdis*eflowmeth(mi)
                  dqmh(mi)=dskmethw1(mi)*eflowmeth(mi)
               endif
            endif
         enddo

         do     mid=1,neq
            mi=mid+ndummy
            denmethi(mi)=(sto1(mi)  -denmethh(mi))*dtin
            denemethi(mi)=(sto1(mi+neq)-denemethh(mi))*dtin
         enddo    

c     ****************************************************************
      else if(iflg.eq.3) then
c     
c     calculations for methane hydrate
c     
         iieosl=0
         do mid=1,neq
            mi=mid+ndummy
c     phase information now contained in array ices
            ieosd=ices(mi)

            pl=phimeth(mi)
            tl=tmeth(mi)
c     water and methane are independent
c     gaz 2-25-2003
            frac_hyd  = frachyd(mi)
            dfrac_hyd = 1.0d00
            dfrac_w = -1.0d00
c     first zero out sources and derivatives
            qhhyd(mi)=0.0
            dq(mi)=0.0
            dqt(mi)=0.0
            dqh(mi)=0.0
            deqh(mi)=0.0
            deqpm(mi)=0.0
            deqm(mi)=0.0
            por = ps(mi)
            dporpl=dporp(mi)
            dportl=dport(mi)


c     properties of solid phase
c     it is assumed that solid and vapor can not exist simultaneously at the same node

c     calculate solid enthalpy and derivative wrt p and t

            call hydrate_properties(1,1,pl,tl,dum2,dum3,
     &           enl,dhlp,dhlt,dum5,dum6)

c     calculate solid density and derivative wrt p and t

            call hydrate_properties(2,1,pl,tl,dum2,dum3,
     &           rol,drolp,drolt,dum5,dum6)

c     porh is pass-through variable ps(mi) to hydrate_properties

            porh = ps(mi)
            volh = sx1(mi)

            call hydrate_properties(3,1,pl,tl,fracw(mi),frachyd(mi),
     &           skhyd(mi),dskhydp,dskhydt,dskhydw,dskhydm)
c     
c     save rate and derivatives for enthalp calculation (passed through commeth)
c     

            source_hyd_temp = skhyd(mi)
            ds1_tmp = dskhydp
            ds2_tmp = dskhydt
            ds3_tmp = dskhydw
            ds4_tmp = dskhydm
c     
c     calculate solid gas heat generation term and derivative wrt p and t and frachyd
c     kinetic approach
c     porh is pass-through variable ps(mi) to hydrate_properties

            porh = ps(mi)
            volh = sx1(mi)

            call hydrate_properties(4,1,pl,tl,fracw(mi),frachyd(mi),
     &           qhhyd(mi),dqhhydp,dqhhydt,dqhhydw,dqhhydm)

c     form dissociation flow term
            skhyd(mi) = rathyd*skhyd(mi)
c     form heat dissociation flow term (minus means add heat to system)
            qhhyd(mi) = rathyd*qhhyd(mi)

c     load gas generation term derivatives 

	    dskhyd1(mi) = rathyd*dskhydp
	    dskhyd2(mi) = rathyd*dskhydt
	    dskhyd3(mi) = rathyd*dskhydw
	    dskhyd4(mi) = rathyd*dskhydm
	    dqhhyd1(mi) = rathyd*dqhhydp
	    dqhhyd2(mi) = rathyd*dqhhydt
	    dqhhyd3(mi) = rathyd*dqhhydw
	    dqhhyd4(mi) = rathyd*dqhhydm

c     compressed liquid
            xvisl=1.d20
            xvisv=1.d20

c     accumulation terms
            den=rol*por*frac_hyd
            dene=por*(rol*enl-pl)*frac_hyd
            sto1(mi)=den
            sto1(mi+neq)=dene

c     derivatives of accumulation terms
            damp=(rol*dporpl+por*drolp)*dtin*frac_hyd
            damm=rol*por*dtin*dfrac_hyd
c     damw=rol*por*dtin*dfrac_w
            damw=0.0
            daep1=(dporpl*rol*enl+por*drolp*enl+por*rol*dhlp
     &           -dporpl*pl-por)*frac_hyd
c     daep2=-dporpl*cp*tl
            daep2=0.0d00             
            daep=(daep1+daep2)*dtin
            damh=(dportl*rol+por*drolt)*dtin*frac_hyd
c     daeh=(-dportl*cp*tl+(1.-por)*cp+(dportl*(enl*rol-pl)+
            daeh=((dportl*(enl*rol-pl)+
     &           por*(rol*dhlt+enl*drolt))*frac_hyd)*dtin
            daem=por*(rol*enl-pl)*dtin*dfrac_hyd
c     daew=por*(rol*enl-pl)*dtin*dfrac_w
            daew=0.0                              
            if(ps(mi).ne.0.0) then

c     store derivatives of accumulation terms
               dmef(mi)=damh
               dmpf(mi)=damp
               dmmf(mi)=damm
               dmwf(mi)=damw
               depf(mi)=daep
               deef(mi)=daeh
               demf(mi)=daem
               dewf(mi)=daew
               dtpameth(mi)=dtpa(mi)
               dtpaemeth(mi)=dtpae(mi)
               tmeth(mi)=tl
               smeth(mi)=s(mi)

            endif
         enddo

         do     mid=1,neq
            mi=mid+ndummy
            denhydi(mi)=(sto1(mi)  -denhydh(mi))*dtin
            denehydi(mi)=(sto1(mi+neq)-denehydh(mi))*dtin
         enddo    

c     end if block on different fluids 
      endif
      
      deallocate(sto1)

      return
      end

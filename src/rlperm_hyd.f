      subroutine rlperm_hyd(ndummy,iz,i,prop1,dprop11,dprop12,
     &           dprop13,dprop14,prop2,dprop21,dprop22,dprop23,dprop24)
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
!D1 To compute the relative permeabilities for vapor and liquid.
!D1 IN the presence of hydrate 
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 12-23-2003, Programmer: G. Zyvoloski
!D2                                     
!D2 $Log:   /pvcs.config/fehm90/src/rlperm_hyd.f_a  $
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
c
c calculates relative permiabilities and derivatives
c
      use comhi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comki
      use comki
      use commeth
      implicit none

      real*8 scutm,hmin,darcyf,tol_l,tol_u,su_cut,rlp_tolh
      real*8 dsh, denom_tol
      parameter(scutm = 1.d-03)
      parameter(hmin = 1.d-18)
      parameter(darcyf = 1.d12)
      parameter(tol_l  = 1.d-5)
      parameter(tol_u  = 1.d-5)
      parameter(su_cut = 0.99d00)
      parameter(rlp_tolh = 0.00)
      parameter(denom_tol = 1.d-08)


      integer iz,ndummy,i,irlpd,mi,ieosd,it,ir,j,num_models,ireg
      real*8 alpha,beta,alamda,alpi,smcut,slcut,fac,ds,dhp,rl,rv
      real*8 drls,drvs,rp1,rp2,rp3,rp4,denom,star,hp,rl1,rv1
      real*8 drls1,drvs1,akf,akm,porf,permb,sl
      real*8 smcutm,smcutf,alpham,alamdam,facf
      real*8 rpa1, rpa2, rpa3, rpa4, rpa5      
      logical null1,ex
      real*8 prop1,dprop11,dprop12,dprop13,dprop14
      real*8 prop2,dprop21,dprop22,dprop23,dprop24
      real*8 starv,drlh,drvh,drp1h,drp2h,rp5,rp6,fhyd
	real*8 sakamoto_a1, sakamoto_b1

      if(iz.eq.0) then
c     
c     calculate relative perm and derivatives
c     
c     note: must be called for each node i
            mi = i+ndummy
c           ieosd = ieos(mi)
            ieosd = 2        
            it = irlp(mi)
            if(it.eq.0) then
               irpd=0
            else
               irpd = irlpt(it)
            endif
            if(irpd.ge.3.and.irpd.le.8) then 
              ieosd = 2
            endif
            if(ieosd.eq.1) then
c     liquid only
               rl = 1.0
               rv = 0.0
               drls = 0.0
               drvs = 0.0
            elseif(ieosd.eq.3) then
c     superheated vapor
               rl = 0.0
               rv = 1.0
               drls = 0.0
               drvs = 0.0
            elseif(ieosd.eq.4) then
c     no permeability, no porosity
               rl = 0.0
               rv = 0.0
               drls = 0.0
               drvs = 0.0
            elseif(ieosd.eq.2) then
               rp1 = rp1f(it)
               rp2 = rp2f(it)
               rp3 = rp3f(it)
               rp4 = rp4f(it)
               sl = fracw(mi)
               sv = 1.0-sl-frachyd(mi)
               if(irpd.eq.1) then
c     two-phase mixture
c     linear function of saturations
                  rl = 0.0
                  rv = 0.0
                  drls = 0.0
                  drvs = 0.0
                  if(sl.gt.rp1) then
                     rl = 1.0
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=(sl-rp1)/denom
                        drls = 1.0/denom
                     endif
                  endif
                  if(sv.gt.rp2) then
                     rv = 1.0
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=(sv-rp2)/denom
                        drvs=-1.0/denom
                     endif
                  endif
	            if (rv.le.0.0) then
	              rv = 0.0
	              drvs = 0.0
	            endif
               else if(irpd.eq.10) then
c     two-phase mixture
c     linear function of saturations
c     minimum rlperms
                  rl = rp5f(it)
                  rv = rp6f(it)
                  drls = 0.0
                  drvs = 0.0
                  if(sl.ge.rp1) then
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=(sl-rp1)/denom+rp5f(it)
                        drls = 1.0/denom
                     endif
                     if(sl.gt.rp3.or.rl.gt.1.0) then
                      rl=1.0
                      drls=0.0
                     endif
                  endif
                  if(sv.ge.rp2) then
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=(sv-rp2)/denom+rp6f(it)
                        drvs=-1.0/denom
                     endif
                     if(sv.gt.rp4.or.rv.gt.1.0) then
                      rv=1.0
                      drvs=0.0
                     endif
                  endif
               elseif(irpd.eq.-1) then
c constant rl perms
                  rl = rp1
                  rv = rp2
                  drls = 0.0
                  drvs = 0.0
               elseif(irpd.eq.2) then
c     corey relationships(geothermal)
                  rl = 0.0
                  rv = 0.0
                  drls = 0.0
                  drvs = 0.0
                  denom = 1.0-rp1-rp2
                  star = (sl-rp1)/denom
                  starv = (sl+frachyd(mi)-rp1)/denom
                  ds = 1.0/denom
                  if(starv.le.0.0) then
                     rv = 1.0
                  else if(star.ge.1.0) then
                     rl = 1.0
                  else
                     rl = star**4
                     drls = 4.0*star**3*ds
                     rv=(1.0-starv)**2*(1.0-starv**2)
                     drvs=(-2.0*(1.0-starv)*(1.0-starv**2)+
     2                    (1.0-starv)**2*(-2.0*starv))*ds
                  endif
               elseif(irpd.eq.3) then
c rlp(h)
                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it),su_cut ,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                 dhp, rl, drls, rv, drvs,0)
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.8) then
c rlp(h) and vapor rlp
                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it),su_cut ,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                 dhp, rl, drls, rv, drvs,1)
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.9) then
c rlp(h) and vapor rlp
                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it),su_cut ,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  rpa1 = rp13f(it)
                  rpa2 = rp14f(it)
                  rpa3 = rp15f(it)
                  rpa4 = rp16f(it)
                  rpa5 = rp17f(it)
                  call vgrlpa(sl, star, rp3, rp4, hmin, hp,
     2                 dhp, rpa1, rpa2, rpa3, rpa4, rpa5,
     3                 rp1f(it), rp2f(it), rl, drls, rv, drvs)
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.5) then
c rlp(S)
                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it), su_cut,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  call  vgrlps(2, sl, star, rp3, rp4, rp1, rp2, 
     2                  tol_l, tol_u, rl, drls, rv, drvs )
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.4) then
c     
c     akf-fracture saturated permeability(rp15)
c     akm-matrix saturated permeability(rp16)
c     porf-fracture fraction(rp17)
c     
                  akf = rp17f(it)
                  akm = rp18f(it)
                  porf = rp19f(it)
                  
                  if( idpdp .ne. 0 .or. idualp .ne. 0 ) then
                     
                     if( mi .le. neq ) then
                      call vg_regions(2,ireg,mi,su_cut)
                        call vgcap( sl, rp11f(it), rp12f(it), rp13f(it),
     2                       rp14f(it), rp20f(it), rp21f(it), rp22f(it),
     3                       rp23f(it), rp16f(it), su_cut,    
     4                       cp3f(it),cp4f(it),hp, dhp, ireg    )
                        star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                        call vgrlp(sl, star, rp13f(it), rp14f(it), hmin,
     2                       hp, dhp, rl, drls, rv, drvs,0)
                        permb = akf * porf
                        
                     else
                      call vg_regions(1,ireg,mi,su_cut)
                        call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                       rp7f(it), rp8f(it), rp9f(it),
     3                       rp10f(it), rp6f(it), su_cut,      
     4                       cp1f(it),cp2f(it),hp, dhp, ireg    )
                        star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                        call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                       dhp, rl, drls, rv, drvs,0)
                        permb = akm * (1. - porf)
                     end if
                     
                  else
                     call vg_regions(1,ireg,mi,su_cut)
                     call vgcap( sl, rp1, rp2, rp3, rp4,  
     2                    rp7f(it), rp8f(it), rp9f(it),
     3                    rp10f(it), rp6f(it), su_cut,      
     4                    cp1f(it),cp2f(it),hp, dhp, ireg    )
                     star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                     call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                    dhp, rl, drls, rv, drvs,0)
                     star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                     call vgrlp(sl, star, rp13f(it), rp14f(it), hmin,
     2                    hp, dhp, rl1, drls1, rv1, drvs1,0)
                     permb = akf * porf + akm * (1. - porf)
                     rl=(akf*rl1*porf+akm*rl*(1.0-porf))/permb
                     drls=(akf*drls1*porf+akm*drls*(1.0-porf))/permb
                     rv = 1. - rl
                     drvs = -drls
                  endif
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
c     set capillary pressures
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp

               elseif(irpd.eq.6.or.irpd.eq.7) then
c     
c     akf-fracture saturated permeability(rp15)
c     akm-matrix saturated permeability(rp16)
c     porf-fracture fraction(rp17)
c     rlp(S)     
c     
                  akf = rp17f(it)
                  akm = rp18f(it)
                  porf = rp19f(it)
                  
                  if( idpdp .ne. 0 .or. idualp .ne. 0 ) then
                     
                     if( mi .le. neq )  then
                      call vg_regions(2,ireg,mi,su_cut)
                        call vgcap( sl, rp11f(it), rp12f(it), rp13f(it),
     2                       rp14f(it), rp20f(it), rp21f(it), rp22f(it),
     3                       rp23f(it), rp16f(it), su_cut,    
     4                       cp3f(it),cp4f(it),hp, dhp, ireg    )
                        star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                        call  vgrlps(2, sl, star, rp13f(it),
     2                  rp14f(it), rp11f(it), rp12f(it),
     3                  tol_l, tol_u, rl, drls, rv, drvs)
                        permb = akf * porf
                      if(irpd.eq.7) then
c calculate fracture term(from Sandia) if necessary
                       call rlp_frac(1,mi,sl,star,rp11f(it),
     &                    rp12f(it), rl,drls,1.-rl,-drls,rp24f(it))
                      endif
                        
                     else
                     call vg_regions(1,ireg,mi,su_cut)
                        call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                       rp7f(it), rp8f(it), rp9f(it),
     3                       rp10f(it), rp6f(it), su_cut,      
     4                       cp1f(it),cp2f(it),hp, dhp, ireg    )
                        star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                        call  vgrlps(2, sl, star, rp3, rp4, rp1, rp2, 
     2                  tol_l, tol_u, rl, drls, rv, drvs )
                        permb = akm * (1. - porf)
                      if(irpd.eq.7) then
c calculate fracture term(from Sandia) if necessary
                       call rlp_frac(1,mi,sl,star,rp1,rp2,
     &                    rl,drls,1.-rl,-drls,rp24f(it))
                      endif
                     end if
                     
                  else
                     call vg_regions(1,ireg,mi,su_cut)
                     call vgcap( sl, rp1, rp2, rp3, rp4,  
     2                    rp7f(it), rp8f(it), rp9f(it),
     3                    rp10f(it), rp6f(it), su_cut,      
     4                    cp1f(it),cp2f(it),hp, dhp, ireg    )
                     star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                     call  vgrlps(2, sl, star, rp3, rp4, rp1, rp2, 
     2                  tol_l, tol_u, rl, drls, rv, drvs )
                     star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                        call  vgrlps(2, sl, star, rp13f(it),
     2                  rp14f(it), rp11f(it), rp12f(it),
     3                  tol_l, tol_u, rl1, drls1, rv1, drvs1)
                     permb = akf * porf + akm * (1. - porf)
                     rl=(akf*rl1*porf+akm*rl*(1.0-porf))/permb
                     drls=(akf*drls1*porf+akm*drls*(1.0-porf))/permb
                     rv = 1. - rl
                     drvs = -drls
                  endif
c     upstream weighting of intrinsic permeability

                  if(iupk.gt.0) then
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
c     set capillary pressures
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp


            elseif(irpd.eq.11) then
c Brooks-Corey
            call bcrlp(sl,star,hp,dhp,rl,drls,rv,drvs)
           elseif(irpd.eq.12) then
c Function to fit Japanese data
c     two-phase mixture
c     linear function of saturations
c     residual saturations are functions 
c rp1f is constant term for irreducible water saturation
c rp2f is constant term for irreducible gas saturation
c rp5f is linear term for irreducible water saturation
c rp6f is linear term for irreducible gas saturation
               rl = 0.0
               rv = 0.0
               drls = 0.0
               drlh = 0.0
               drvs = 0.0
               drvh = 0.0
               rpa1 = rp1f(it)
               rpa2 = rp2f(it)
               rp3 = rp3f(it)
               rp4 = rp4f(it)
               rp5 = rp5f(it)
               rp6 = rp6f(it)
               sl = fracw(mi)
               sv = 1.0-sl
               fhyd =frachyd(mi)
               rp1 = rpa1 + rp5*fhyd
               rp2 = rpa2 + rp6*fhyd
               drp1h = rp5
               drp2h = rp6
                  if(sl.gt.rp1) then
                     rl = 1.0
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=(sl-rp1)/denom
                        drls = 1.0/denom
                        drlh = -drp1h/denom - (sl-rp1)/denom**2*drp1h
                     endif
                  endif
                  if(sv.gt.rp2) then
                     rv = 1.0
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=(sv-rp2)/denom
                        drvs= -1.0/denom
                        drvh = -drp2h/denom - (sv-rp2)/denom**2*drp2h
                     endif
                  endif
	            if (rv.le.0.0) then
	              rv = 0.0
	              drvs = 0.0
                    drvh = 0.0
	            endif
 
           elseif(irpd.eq.13) then
c Function to fit Japanese data
c     two-phase mixture
c     linear function of saturations
c     residual saturations are functions 
c rp1f is constant term for irreducible water saturation
c rp2f is constant term for irreducible gas saturation
c rp3f is maximum liquid relative permeability
c rp4f is maximum gas relative permeability
c rp5f is linear term for irreducible water saturation
c rp6f is linear term for irreducible gas saturation
c new model for fitting experiments march 9 2004
c introduced two new parameters
               rl = 0.0
               rv = 0.0
               drls = 0.0
               drlh = 0.0
               drvs = 0.0
               drvh = 0.0
               rpa1 = rp1f(it)
               rpa2 = rp2f(it)
               rp3 = rp3f(it)
               rp4 = rp4f(it)
               rp5 = rp5f(it)
               rp6 = rp6f(it)
               sl = fracw(mi)
               fhyd =frachyd(mi)
               rp1 = rpa1 + rp5*fhyd
               rp2 = rpa2 + rp6*fhyd
               drp1h = rp5
               drp2h = rp6
               denom = 1. - rp2 - rp1 - fhyd + denom_tol
               star = (sl - rp1)/denom
               ds = 1./denom
               dsh = (-drp1h)/denom - 
     &               (sl-rp1)*(-drp1h-drp2h-1.)/denom**2
c changed by RJP 6/1/04 to give option to specify the following
c parameters
			sakamoto_a1 = rp7f(it)
			sakamoto_b1 = rp8f(it)
c					sakamoto_a1=.4
c					sakamoto_b1=.4
					
					rl= sakamoto_a1*star
c changed the following expressions RJP 6/1/04
c					drls = ds
c					drlh = dsh
					drls = sakamoto_a1*ds
					drlh = sakamoto_a1*dsh
										
					rv = sakamoto_b1*(1.-star)
c changed the following expressions RJP 6/1/04
c					drvs = -ds
c					drvh = -dsh   
					drvs = -sakamoto_b1*ds
					drvh = -sakamoto_b1*dsh
c                       rl= star
c                       drls = ds
c                       drlh = dsh
      
c                       rv = 1. - star
c                       drvs = -ds
c                       drvh = -dsh
 
	            if (rv.le.0.0) then
	              rv = 0.0
	              drvs = 0.0
                    drvh = 0.0
	            endif
	            if (rl.le.0.0) then
	              rl = 0.0
	              drls = 0.0
                    drlh = 0.0
	            endif
	            if (rv.ge.1.0) then
	              rv = 1.0
	              drvs = 0.0
                    drvh = 0.0
	            endif
	            if (rl.ge.1.0) then
	              rl = 1.0
	              drls = 0.0
                    drlh = 0.0
	            endif
c              endif
           elseif(irpd.eq.14) then
c Rel perms with max rel perms
c     two-phase mixture
c     linear function of saturations
c     residual saturations are functions 
c rp1f is constant term for irreducible water saturation
c rp2f is constant term for irreducible gas saturation
c rp3f is maximum water rl perm
c rp4f is maximum gas rv perm
               rl = 0.0
               rv = 0.0
               drls = 0.0
               drlh = 0.0
               drvs = 0.0
               drvh = 0.0
               rp1 = rp1f(it)
               rp2 = rp2f(it)
               rp3 = 1.0
               rp4 = 1.0
               sl = fracw(mi)
               fhyd =frachyd(mi)
               sv = 1.0-sl-fhyd
                   if(sl.gt.rp1) then
                     rl = 1.0
c
c modified for max rl perm (rp3f)
c
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=rp3f(it)*(sl-rp1)/denom
                        drls = rp3f(it)*1.0/denom
                     endif
                  endif
c
c modified for max rv perm (rp4f)
c
                  if(sv.gt.rp2) then
                     rv = 1.0
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=rp4f(it)*(sv-rp2)/denom
                        drvs= -rp4f(it)*1.0/denom
                     endif
                  endif
               else if(irpd.eq.15) then
c     two-phase mixture for methane hydrate
c     linear function of saturations
               rl = 0.0
               rv = 0.0
               drls = 0.0
               drlh = 0.0
               drvs = 0.0
               drvh = 0.0
               rp1 = rp1f(it)
               rp2 = rp2f(it)
               rp3 = 1.0
               rp4 = 1.0
               sl = fracw(mi)
               fhyd =frachyd(mi)

                     rl = sl - fhyd - rp1 + rlp_tolh
                     drls = 1.0
                     drlh = -1.0
 
                     rv = 1.- sl - fhyd - rp2
                     drvs = -1.
                     drvh = -1.

	            if (rl.ge.1.0 - rp1 + rlp_tolh) then
	              rl = 1.0 - rp1 + rlp_tolh
	              drls = 0.0
                    drlh = 0.0
	            endif
	            if (rv.ge.1.0 - rp2) then
	              rv = 1.0 - rp2 
	              drvs = 0.0
                    drvh = 0.0
	            endif
	            if (rl.le.0.0) then
	              rl = 0.0
	              drls = 0.0
                    drlh = 0.0
	            endif
	            if (rv.le.0.0) then
	              rv = 0.0
	              drvs = 0.0
                    drvh = 0.0
	            endif
            endif
           endif
             if(irpd.eq.12) then
               prop1 = rl
               prop2 = rv
               dprop13 = drls
               dprop14 = drlh
               dprop23 = drvs
               dprop24 = drvh
            else if(irpd.eq.13) then
               prop1 = rl
               prop2 = rv
               dprop13 = drls
               dprop14 = drlh
               dprop23 = drvs
               dprop24 = drvh
            else if(irpd.eq.15) then
               prop1 = rl
               prop2 = rv
               dprop13 = drls
               dprop14 = drlh
               dprop23 = drvs
               dprop24 = drvh
              else
               prop1 = rl
               prop2 = rv
               dprop13 = drls
               dprop14 = 0.0
               dprop23 = drvs
               dprop24 = drvs
              endif

      endif
      return
      end

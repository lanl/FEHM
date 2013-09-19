      subroutine saltctr(iflg,ndummy,dum_salt)
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
!D1 Manage simulations with salt iteractions on porosity, permeability, 
!D1 diffusivity, density etc. 
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 3.1
!D2 
!D2 Initial implementation: George Zyvoloski 071713 using preexisting 
!D2 code written by Phil Stauffer and Dylan Harp
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
      use comchem
      use comsplitts 
      use comrxni, only : cden_flag, cden_sp, cden_spnam, mw_speci

      implicit none

      integer iflg,ndummy,i,nr1,nr2,nr3,nr4,nr5,mi      
      integer i1,i2,i3,i4,i5,i6,ilev,mlev,il,md,k,j, ja, kb, im, iphase
      integer ii,ij,icedc,iced,ieosd,icesd,idco2,duma,neqp1,ncont
      integer totvnum
   
      integer imped_ex_0
      real*8, allocatable :: aiped(:)
      real*8, allocatable :: sktmp(:)
      real*8, allocatable :: esktmp(:)
      integer, allocatable :: iflg_flowmac(:)
      real*8 tliquid,dummyreal
      real*8 tw,pc,tc,ensrc,delsrc,phase_frac,skmd,qhmd
      real*8 skmd1,skmd10,skmd2,skmd3,amco2w,amco2w1,amco21
      real*8 amco2w0,skmd21,skmd31,amco2f,amco2f0,amco2f1
      real*8 dtps,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real*8 pwv, permsb
      real*8 strd_co2,strd1
      real*8 eostol,ps_tol
      real*8 eosmg
      real*8 eosml
      real*8 deltc, delpc, tem(9)
      real*8 dum_salt
      character*80 input_msg, dummy_line
      character*32 cmsg(4), cdum
      character*5 cden_type
      parameter(deltc=31.1d0,delpc=3.78d0, ps_tol = 1.d-03)
      
      real*8 amaxflx, dums1,tolw, zero_e

      real*8 pcrit, tcrit,dumb
      real*8 psatd, tl1,tl2,dtps1,dtps2,dpsatt
      real*8 sl,sg,sw, pl, tl, tsolid
      real*8 denc,dencp,denct,enc,encp,enct,visc,viscp,visct
      real*8 denw,enwp,enwt,viswp,viswt
      real*8 ycmax, xwp, xcp, xcpb, mwc, mww, xc_prime, tmp1
      real*8 psatl,pv,pcl, dpsats
c
c salt related local variables
c
      real*8 spec1,fracc,ms,xf,af,bf,cf,df,ef,dumf,tltemp
      real*8 strd_salt, por_salt_min_0, pnx_new, ps_new
      real*8 xmsg(1), mw
      integer cnum,iieosd,inptorig,kk,msg(20),nwds,imsg(20)
      integer ico2d, ico2dc
      character*8 macro1
      parameter(eostol=0.0001d00)
      parameter(eosmg=1.0001d00)
      parameter(eosml=0.9999d0)
      parameter(strd1=1.0d0)
      parameter(tolw=1.d-90, zero_e=1.e-10)
      parameter(mwc=44.d0,mww=18.d0)
      parameter(imped_ex_0 = 3)
      parameter(por_salt_min_0=1.0d-8)
      integer  open_file, myfile_salt
      
      logical it_is_open
   
      save myfile_salt
c gaz debug 090113
      i = l
      if(isalt.eq.0) return  

       if(nr_stop.eq.0) then
        strd_salt=0.7d0
       else
        strd_salt = strd_iter
       endif

         if(iflg.eq.0) then

c     
c     allocate memory for salt arrays
c     
            if(.not.allocated(ps_old)) allocate(ps_old(n0))
            if(.not.allocated(pnx_old)) allocate(pnx_old(n0))
c     
c     zero out arrays for salt simulations, initialize salt parameters        
c     
            isalt_ps = 0
            isalt_pnx = 0
            ivaprsalt = 0
            ps_old = 0.0d00
            pnx_old = 0.0d00
            por_salt_min = por_salt_min_0

            macro = 'salt'
c     
c     read in "sub macros" for salt
c     

 100        continue
            read (inpt, '(a80)') wdd1
c Changed end check for compatibility with macro "off" option
            if (wdd1(1:7) .eq. 'saltend' .or. wdd1(1:7) .eq. 'endsalt'
     &           .or. wdd1(1:8) .eq. 'end salt') go to 200 
            if (wdd1(1:1) .eq. '#') go to 100 
            read (wdd1, '(a8)') macro1
            if (iout .ne. 0) write(iout, 50) macro1
            if (iptty .gt. 0) write(iptty, 50) macro1 
 50         format(3x, '**** salt sub macro : ', a8,' **** ' ) 
            if(macro1.eq.'saltvcon') then
c
c salt thermal conductivity
c
c set to arbitrary max of 20 models
c
      if(.not.allocated(vc1f)) then
       totvnum = 20
       allocate(vc1f(totvnum),vc2f(totvnum),vc3f(totvnum),vc4f(totvnum),
     &    vc5f(totvnum),vc6f(totvnum),vc7f(totvnum),vc8f(totvnum))   
      else
       totvnum = 20
       deallocate(vc1f,vc2f,vc3f,vc4f,vc5f,vc6f,vc7f,vc8f)  
       allocate(vc1f(totvnum),vc2f(totvnum),vc3f(totvnum),vc4f(totvnum),
     &    vc5f(totvnum),vc6f(totvnum),vc7f(totvnum),vc8f(totvnum))  
      endif
             ivcond = 1 
             call vcon (0, 0)
            elseif(macro1.eq.'saltppor') then             
c     
c  call porosi with model 6 or 7 for for salt (otherwise)
c  print out an error statement
c  iporos changed in call to porosi(0)  
c    
           iporos = 1
c
c allocate memory for models 6 or 7
c
c        parameters for temperature-dependent porosity
           if(.not.allocated(dporp)) allocate(dporp(n8))
           if(.not.allocated(dport)) allocate(dport(n8))
           if(.not.allocated(porTemp1)) then
            allocate(porTemp1(n8),porTemp2(n8),porTemp3(n8),
     &        porTemp4(n8))
           endif
             call porosi(0)    
c write out warning  for non-salt iporos
             if(iporos.ne.6.or.iporos.ne.7) then
              write(ierr,901) iporos
              if(iout.ne.0) write(iout,901) iporos
              if(iptty.ne.0) write(iptty,901) iporos
901    format('warning:porosity model ',i3,' used in salt simulation')
             endif   
            elseif (macro1.eq.'saltadif') then
c**** air-water vapor diffusion with salt (adif = 333 or 666)****
             iadif = 1
             read(inpt,*) tort
c write out warning for non-salt tort
             if(tort.ne.333.0.or.tort.ne.666.0) then
              write(ierr,9022) tort
              if(iout.ne.0) write(iout,9022) tort
              if(iptty.ne.0) write(iptty,9022) tort
9022    format('warning: adif model ',f7.2,' used in salt simulation')
             endif   
            elseif (macro1.eq.'saltvapr') then
c manage vapor pressure lowering with salt 
c Sparrow (2003) Desalination
c Sparrow has no capillary vapor pressure lowering
c ivaprsalt-vapor pressure lowering model
c ivaprsalt = 0, standard (traditional FEHM h2o vapor pressure fit)
c                no lowering (even with capillary pressury)
c ivaprsalt = 1, Sparrow vapor pressure model with no salt
c ivaprsalt = 2, Sparrow vapor pressure model with with salt
c ivaprsalt = 3, Sparrow vapor pressure model with with salt
c                and capillary pressure vapor pressure lowering (not operational)
c ivaprsalt = 4, traditional FEHM vapor pressure model without 
c                 capillary pressure vapor pressure lowering (no salt)
c ivaprsalt = 5, traditional FEHM vapor pressure model with 
c                 capillary pressure vapor pressure lowering (no salt
c ivaprsalt = 6, traditional FEHM pressure model with with salt
c                and NO capillary pressure vapor pressure lowering
c ivaprsalt = 7, traditional FEHM pressure model with with salt
c                and capillary pressure vapor pressure lowering

              ivapl = 0
              read(inpt,*) ivaprsalt
              if(ivaprsalt.eq.3.or.ivaprsalt.eq.5.or.
     &         ivaprsalt.eq.7) ivapl = 1 
              if(ivaprsalt.eq.2.or.ivaprsalt.eq.3.or.
     &          ivaprsalt.eq.6.or.ivaprsalt.eq.7) then
c check for errors (need tracer for salt)
               if(iccen.eq.0) then
                write(ierr,902) 
                if(iout.ne.0) write(iout,902)
                if(iptty.ne.0) write(iptty,902)
902      format('vapor pressure model needs salt tracer:stopping')
               stop       
               endif
              endif

            else if(macro1.eq.'saltden ') then              
c     
c manage salt density  
c 
c new form cden FEHM ver 3.2
         backspace (inpt)
         read (inpt, '(a80)') input_msg
         call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
         if (nwds .gt. 1) then
            cden_type = cmsg(2)
         else
            cden_type = 'default'
         end if
         select case (cden_type(1:5))
         case ('multi', 'table')
            cden_flag = 1
         case ('co2')
            cden_flag = 2
            if (nwds .gt. 2 .and. msg(3) .eq. 3) then
               cden_spnam = cmsg(3)
            else
               cden_spnam = 'Na'
            end if
         case default
            cden_flag = 0
            read(inpt,*) ispcden
            read(inpt,*) factcden
            if (nspeci .gt. 0) then
               if(ispcden.le.nspeci) then
                  cden = .true.
               else
                  if (iout .ne. 0) write (iout, 1005)
                  if( iptty .gt. 0) write (iptty, 1005)
               end if 
            end if        
         end select
         if(nspeci .eq. 0) then
            write (ierr, 1006)
            if (iout .ne. 0) write (iout, 1006)
            if (iptty .gt. 0) write (iptty, 1006)
         end if
 1005    format ('No solute transport, cden macro ignored')
 1006    format ('ispcden > nspeci, cden macro ignored')


            else if(macro1.eq.'saltnum ') then
c
c manage the implicit and explicit updates for salt parameters
c
c update averaging for perms 
c update averaging for porosity
c
301          read(inpt,'(a80)') wdd1(1:80)
             do i = 1,70
              if(wdd1(i:i+6).eq.'permavg') then
               read(wdd1(i+7:80),*) permavg_salt
               isalt_pnx = 1
               go to 301
              endif
              if(wdd1(i:i+5).eq.'poravg') then
               read(wdd1(i+6:80),*) poravg_salt
               isalt_ps = 1
               go to 301
              endif
              if(wdd1(i:i+5).eq.'pormin') then
               read(wdd1(i+6:80),*) por_salt_min
               go to 301
              endif
             enddo
c check for blank line
             do i = 1,80
              if(wdd1(i:i).ne.' ') then
                write(ierr,*) 'error in salt (saltnum), stopping' 
               if(iout.ne.0) write(iout,*) 
     &          'error in salt (saltnum), stopping'
               if(iptty.ne.0) write(iptty,*) 
     &          'error in salt (saltnum), stopping'
                stop
              endif
             enddo
            else if(macro1.eq.'saltflow') then


               
            else
               if (iout .ne. 0) write(iout,*) 
     &              'ERROR IN SALT INPUT(STOPPING)'
               write(ierr,*) 'ERROR IN SALT INPUT(STOPPING)'
            end if
 40         continue
            go to 100
 200        continue
                       
            macroread(8) = .TRUE.
            
        else if(iflg.eq.1) then
c     
c Sparrow has no capillary vapor pressure lowering
c ivaprsalt-vapor pressure lowering model
c ivaprsalt = 0, standard (traditional FEHM h2o vapor pressure fit)
c                no lowering (even with capillary pressury)
c ivaprsalt = 1, Sparrow vapor pressure model with no salt
c ivaprsalt = 2, Sparrow vapor pressure model with with salt
c ivaprsalt = 3, Sparrow vapor pressure model with with salt
c                and capillary pressure vapor pressure lowering
c ivaprsalt = 4, traditional FEHM vapor pressure model without 
c                 capillary pressure vapor pressure lowering
c ivaprsalt = 5, traditional FEHM vapor pressure model with 
c                 capillary pressure vapor pressure lowering
c ivaprsalt = 6, traditional FEHM pressure model with with salt
c                and NO capillary pressure vapor pressure lowering
c ivaprsalt = 7, traditional FEHM pressure model with with salt
c                and capillary pressure vapor pressure lowering

      dpsatt=0.0
c
      mi = ndummy
      tl = t(mi)
      pl = phi(mi)
      if(ieos(mi).eq.2) then
         if(ivaprsalt .le. 3) then
    

c PHS  11/28/12  Adding term for chemical lowering of vapor pressure.
c     only on species #1,  anl(i = 1,NEQ) 
c     Use simple linear fraction of mass calculation to reduce vapor pressure.
c      assume salt 58.55 g/mole   an is moles per kg water
c     Neeper suggests   Pwv = 0.76 at 1 mole NaCL/liter water
c  Originally done 6/5/2006 for species #2
c     Now in terms of mol fraction (moles/kg / 55.55 moles/kg water)

c            spec1 = an(mi)/55.55
c            fracc = -16.801*spec1*spec1 - 3.081*spec1 + 1.00
c
c            pcl=pl-fracc*psatl(tl,pcp(mi),dpcef(mi),dpsatt,dpsats,0)
c              pv=(pl-pcl)

c DRH  12/03/12 Added term for chemical lowering of vapor pressure.
c      Relation from B.S. Sparrow (2003) Desalination 159, 161-170
c		Determine mass of salt (ms), multiply by 58.55 g/mole * 1kg/1000g
c		assuming mass of water is 1kg

c		Determine mass fraction (xf) assuming mass of water (mw) is 1kg
c		following assumption in calculation of ms above, xf = ms / (ms + mw)

c        1         2         3         4         5         6        7 2       8    

        if(ivaprsalt .eq. 2 .or. ivaprsalt .eq. 7) then
			 ms = an(mi) * 58.55 / 1000.
	   xf = ms / ( ms + 1.0 )
        else
         xf = 0.0
        endif
		
			if(tl.ge.0.0.and.tl.le.150) then
		af = ( 0.9083 - 0.5690 * xf + 0.1945 * xf**2 - 3.7360
     & 		    * xf**3 + 2.820 * xf**4) * 1.e-3
		bf = (-0.0669 + 0.0582 * xf - 0.1668 * xf**2 + 0.6761
     &          * xf**3 - 2.091 * xf**4) * 1.e-3
		cf = ( 7.5410 - 5.1430 * xf + 6.4820 * xf**2 - 52.620
     &          * xf**3 + 115.7 * xf**4) * 1.e-6
		df = (-0.0922 + 0.0649 * xf - 0.1313 * xf**2 + 0.8024
     &          * xf**3 - 1.986 * xf**4) * 1.e-6
		ef = ( 1.2370 - 0.7530 * xf + 0.1448 * xf**2 - 6.9640
     &          * xf**3 + 14.61 * xf**4) * 1.e-9
			else if(tl.gt.150.and.tl.le.300.0) then
		af = (-3.2480 + 7.0810 * xf - 49.930 * xf**2 + 219.60
     &          * xf**3 - 308.5 * xf**4)
		bf = ( 0.0610 - 0.1185 * xf + 0.7916 * xf**2 - 3.4740
     &          * xf**3 + 4.882 * xf**4)
		cf = (-0.4109 + 0.6789 * xf - 4.1550 * xf**2 + 18.340
     &          * xf**3 - 25.89 * xf**4) * 1.e-3
		df = ( 1.1300 - 1.4320 * xf + 7.1690 * xf**2 - 33.170
     &          * xf**3 + 47.45 * xf**4) * 1.e-6
				ef = 0.0
			endif
c

c           Truncate function at 300C  and dump a statement to the screen
            if(tl.gt.300) then
                tltemp = 300.
c gaz debug 080613 note removed derivative above 300 C (function is constant)
                dpsatt = 0.0
                dpsats = 0.0
                write(ierr,903) tltemp,tl
                if(iout.ne.0) write(iout,903) tltemp,tl
                if(iptty.ne.0) write(iptty,903) tltemp,tl
903       format('Exceeded temp range ',f10.3,' for salt  T= ',f10.3)
            else
                tltemp = tl
            endif

c           Calculate vapor pressure using eq. 6 from Sparrow 2003, see ref above
              pv = af + bf*tltemp + cf*tltemp**2. + df*tltemp**3. +
     &                      ef*tltemp**4.

c           write(iout,*) 'temp pvap an(mi) xf',tltemp, pv,an(mi), xf
c           write(iout,*) 'af bf cf df ef', af,bf,cf,df,ef

c        1         2         3         4         5         6        7         8  

c           Calculate gas pressure

             dpsatt = bf + 2.0*cf*tltemp +  3.0*df*tltemp**2 +
     &                4.0*ef*tltemp**3

			     pcl = pl - pv
             pci(mi) = pcl
             dum_salt = dpsatt

c DRH 12/03/12 end.
         else if(ivaprsalt .ge. 4) then
c
c original FEHM vaopr pressure
c vapor pressure lowering in psatl enabled via ivapl parameter in comai
c
            pcl=pl-psatl(tl,pcp(mi),dpcef(mi),dpsatt,dpsats,0)
            pv=pl-pcl
            pci(mi) = pcl
            dum_salt = dpsatt
         end if

         endif
        else if(iflg.eq.2) then
c     
c   mineral reaction porosity change (called from csolve)
c
c calculate porosity change due to mineral reactions
         ps_delta_rxn = 0.
c DRH and PHS 12/4/12  So Close.  Just needed (1-porosity) 
c PHS and DHR 7/31/13  Changed to expression:
c      delta por = delta(mol/kg rock)*rock density(kg rock/m3 rock)*
c           *(1-porosity) (m3 rock/m3 tot)
c             *molecular weight(kg salt/mol)/ salt grain density (kg salt/m3 salt)
c      = delat m3 salt /m3 total.
c
c PHS 8/1/13  Also added call to porosi and reset an to anlo to maintain
c             constant mol/kg of salt.
c PHS 8/2/13  Added a throttle to stop the porosity loss below 1e-3 when porosity   
c             is decreasing (ps_delta_rxn -VE)

         do im = 1, nimm
            nsp = pimm(im)
            npn = npt(nsp)

            if(mw_mineral(im).ne.0)then
               do i = 1, n0
                  ja = i + npn
                  ps_delta_rxn(i) = -((an(ja)-anlo(ja))*denr(i)* 
     &                 (1-ps(i))*mw_mineral(im)/rho_mineral(im))

                  if(ps_delta_rxn(i).GT.0) then
                    an(ja) = anlo(ja)
                  else  
                    if(ps(i).GT.ps_tol) an(ja) = anlo(ja)
                    if(ps(i).LE.ps_tol) ps_delta_rxn(i) = 0.0
                 end if       

               enddo
            endif

            ps_delta_rxn_s = ps_delta_rxn
         enddo
c
c save porosity after flow simulation timestep
c
         ps_old(1:n0) = ps(1:n0)
         pnx_old(1:n0) = pnx(1:n0)
	
            call porosi(1)
c
c average porosity and permeability
c new porosity is available after tracer run
c should be called after tracer run is completed
c
       if(isalt_pnx.ne.0) then          
          do i = 1, neq
c
c average new and old permeabilities
c           
             pnx_new = pnx(i)
             pnx(i)  = pnx_old(i)*(1.0d0-permavg_salt) + 
     &                 pnx_new*permavg_salt
             pny(i)  = pnx(i)
             pnz(i)  = pnx(i)         
          enddo
        endif
c
c average new and old permeabilities
c   
       if(isalt_ps.ne.0) then          
          do i = 1, neq
           ps_new = ps(i)
           ps(i)  = ps_old(i)*(1.0d0-poravg_salt) + ps_new*poravg_salt
          enddo
       endif
       continue
         else if(iflg.eq.3) then
c
c save porosity after flow simulation timestep
c
c
         else if(iflg.eq.4) then
c
c average porosity and permeability
c new porosity is available after tracer run
c should be called after tracer run is completed
c

         else if(iflg.eq.-10) then   
c          
c   testing only src 
c   modify vapor pressure lowering for salt
c     
c
c two phase conditions test
c
       deallocate(an,pcp,dpcef)
       allocate(an(100),pcp(100),dpcef(100))
       write(ierr,*) 'pl,tl,dumf,pv,pv/dumf,xf, ms,an'
       an(1:100) = 30.0
       pcp = 0.0
       dpcef = 0.0
        n = 60
        dum1 = (340.-40.)/n
        dum2 = 30./n
        pl = 0.1
        tl = 30.
        do mi = 1, n
         tl = tl + dum1

         if( ivapl .ne. 0 ) then

c PHS  11/28/12  Adding term for chemical lowering of vapor pressure.
c     only on species #1,  anl(i = 1,NEQ) 
c     Use simple linear fraction of mass calculation to reduce vapor pressure.
c      assume salt 58.55 g/mole   an is moles per kg water
c     Neeper suggests   Pwv = 0.76 at 1 mole NaCL/liter water
c  Originally done 6/5/2006 for species #2
c     Now in terms of mol fraction (moles/kg / 55.55 moles/kg water)

c            spec1 = an(mi)/55.55
c            fracc = -16.801*spec1*spec1 - 3.081*spec1 + 1.00
c
c            pcl=pl-fracc*psatl(tl,pcp(mi),dpcef(mi),dpsatt,dpsats,0)
c              pv=(pl-pcl)

c DRH  12/03/12 Added term for chemical lowering of vapor pressure.
c      Relation from B.S. Sparrow (2003) Desalination 159, 161-170
c		Determine mass of salt (ms), multiply by 58.55 g/mole * 1kg/1000g
c		assuming mass of water is 1kg
			ms = an(mi) * 58.55 / 1000
c		Determine mass fraction (xf) assuming mass of water (mw) is 1kg
c		following assumption in calculation of ms above, xf = ms / (ms + mw)

c        1         2         3         4         5         6        7 2       8    

			xf = ms / ( ms + 1.0 )
			if(tl.ge.0.0.and.tl.le.150) then
		af = ( 0.9083 - 0.5690 * xf + 0.1945 * xf**2 - 3.7360
     & 		    * xf**3 + 2.820 * xf**4) * 1.e-3
		bf = (-0.0669 + 0.0582 * xf - 0.1668 * xf**2 + 0.6761
     &          * xf**3 - 2.091 * xf**4) * 1.e-3
		cf = ( 7.5410 - 5.1430 * xf + 6.4820 * xf**2 - 52.620
     &          * xf**3 + 115.7 * xf**4) * 1.e-6
		df = (-0.0922 + 0.0649 * xf - 0.1313 * xf**2 + 0.8024
     &          * xf**3 - 1.986 * xf**4) * 1.e-6
		ef = ( 1.2370 - 0.7530 * xf + 0.1448 * xf**2 - 6.9640
     &          * xf**3 + 14.61 * xf**4) * 1.e-9
			else if(tl.gt.150.and.tl.le.300.0) then
		af = (-3.2480 + 7.0810 * xf - 49.930 * xf**2 + 219.60
     &          * xf**3 - 308.5 * xf**4)
		bf = ( 0.0610 - 0.1185 * xf + 0.7916 * xf**2 - 3.4740
     &          * xf**3 + 4.882 * xf**4)
		cf = (-0.4109 + 0.6789 * xf - 4.1550 * xf**2 + 18.340
     &          * xf**3 - 25.89 * xf**4) * 1.e-3
		df = ( 1.1300 - 1.4320 * xf + 7.1690 * xf**2 - 33.170
     &          * xf**3 + 47.45 * xf**4) * 1.e-6
				ef = 0.0
			endif
c           Call psatl just to set dpsatt and dpsats derivatives
            dumf = psatl(tl,pcp(mi),dpcef(mi),dpsatt,dpsats,0)

c           Truncate function at 300C  and dump a statement to the screen
            if(tl.gt.300) then
                tltemp = 300.
                write(6,*) 'Exceeded temp range for salt  T:', tl,tltemp
            else 
                tltemp = tl
            endif

c           Calculate vaport pressure using eq. 6 from Sparrow 2003, see ref above
              pv = af + bf*tltemp + cf*tltemp**2. + df*tltemp**3. +
     &                      ef*tltemp**4.

c           write(iout,*) 'temp pvap an(mi) xf',tltemp, pv,an(mi), xf
c           write(iout,*) 'af bf cf df ef', af,bf,cf,df,ef

c        1         2         3         4         5         6        7         8  

c           Calculate gas pressure
			 pcl = pl - pv

         write(ierr,'(1p,8(1x,g14.5))') 
     &   pl,tl,dumf,pv,pv/dumf,ms, xf,an(mi)
c DRH 12/03/12 end.
         else
            pcl=pl-psatl(tl,pcp(mi),dpcef(mi),dpsatt,dpsats,0)
            pv=pl-pcl
         end if          
        enddo
c
        write(*,*) 'stopping in saltctr'
        stop
c

         else if(iflg.eq.5) then

         elseif(iflg.eq.6) then
               
         elseif(iflg.eq.7) then  
c
c additional output for salt 
c 
                  if (iout .ne. 0) write(iout,6015) 
                  if (iatty.gt.0) write(iatty,6015) 
                  if (iout .ne. 0) write(iout,6016) 
                  if (iatty.gt.0) write(iatty,6016) 
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
                           if (ieos(md) .eq. 2) then
                              iphase = 2
                           else
                              iphase = 1
                           end if
                        permsb =  pnx(md)*1.0e-6
                        pwv    =  phi(md) - pci(md)

                      if (iout .ne. 0)     
     &                  write(iout,6017) md,permsb,ps(md),thx(md)*1.e6,
     &                     pwv,dvas(md),ps_delta_rxn_s(md)

                      if (iatty .ne. 0)
     &                  write(iatty,6017) md,permsb,ps(md),thx(md)*1.e6,
     &                     pwv,dvas(md),ps_delta_rxn_s(md)
                     enddo
                  enddo

 600           format(2x,'Matrix Level = ',i1)
 6015        format(' - - - - - - - - - - - - - - - - - - - - - - - - ',
     &       '- - - - - - - - - - - - - - - - -')
 6016             format(1x,'  Node', 5x, 'perm (m2)', 6x,
     &                 'porosity', 7x, 'Kx W/(m K)', 5x,'Pwv (MPa)',
     &                6x, 'D*wv (m2/s)',2x,' ps_delta_rxn')
 6017             format(1x,i6,1x,6(3x,g12.5))

  
         elseif(iflg.eq.8) then
                  if (iout .ne. 0) write(iout,6015) 
                  if (iatty.gt.0) write(iatty,6015) 
                  if (iout .ne. 0) write(iout,6018) 
                  if (iatty.gt.0) write(iatty,6018) 
c
 6018             format(1x,'  Node', 5x, 'density (kg/m3)')
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
                           if (ieos(md) .eq. 2) then
                              iphase = 2
                           else
                              iphase = 1
                           end if

                      if (iout .ne. 0)     
     &                  write(iout,6017) md,rolf(md)

                      if (iatty .ne. 0)
     &                  write(iatty,6017) md,rolf(md)
                     enddo
                  enddo

  
         elseif(iflg.eq.9) then
         elseif(iflg.eq.10) then
            
c     advance variable in transient simulation
c     gaz done elsewhere for now
c
            
         elseif(iflg.eq.11) then
            
c     reset variables
         
           ps_old(1:n0) = ps(1:n0)     
            
         elseif(iflg.eq.-20) then

         elseif(iflg.eq.20) then

         elseif(iflg.eq.22) then


         elseif(iflg.eq.21) then
            
c     write surfer or tecplot contour files
            
            write(iscon,790)
 790        format('         ','X',12x,'Y',12x,'Z',6x,'   pres',12x,'t',
     &           8x,'an',8x,'ps',7x,'pci',7x,'s')
            j = ndummy
            do k=1,neq_primary
               if(izone_surf_nodes(k).eq.j) then
                  
                  write(iscon,855) 
     &                 cord(k,1),cord(k,2),cord(k,3),phi(k),
     &                 t(k),an(k),ps(k),pci(k),s(k)
               endif
            enddo
 855        format(1x,10(1x,g12.5))
         endif
      
      return
      end

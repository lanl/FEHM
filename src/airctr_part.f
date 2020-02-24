      subroutine airctr_part(iflg,ndummy,zone)
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
!D1 To manage the isothermal air-water calculations.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/airctr_part.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.2 Heat- and mass-transfer equations
!D3 2.3.3 Noncondensible gas flow equations
!D3 2.5.1 Implement time-step mechanism
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
      use comxi
      use davidi
      use comflow 
      use comsplitts 
      use com_part
      implicit none

      integer zone
      integer iflg,ndummy,ico2d,ndum,ieoss,iieoss,iwelbs,i,mid
      integer mi,i1,i2,ilev,mlev,il,md,irlpsv,irlptsv,nr1
      integer nr2,irdofsv
      integer ja
c gaz 110819 removed tref, pref (now global)       
      real*8 pssv,ssv,phisv,dmpfd,dmefd,dqd,rqd,qcd
      real*8  strd_part
      real*8 inflow_thstime,inen_thstime,denht,deneht
      
      real*8 dels,delp,dfdum11,dfdum12,dfdum21,dfdum22
      real*8 dfdum11i,dfdum12i,dfdum21i,dfdum22i,detdf
      real*8 fdum01,fdum02,sx1d
c  gaz 110819 tref, pref (now global)      
c      save tref,pref

C ? ico2d not passed in, not initialized
      ico2d = 0
c     return if no air present
      if(ico2.lt.0) then
c     iflg is used to tell if call is for initialization
         if(iflg.eq.0) then
            qtc=0.0
            qtotc=0.0
            amc=0.0
c     
c     read in reference pressure and temperature
c     
            read(inpt,*) tref,pref
c     
c     set max and min values of p and t
c     
            pmin(1)=-1.e05
            pmax(1)=1.e05
            tmin(1)=-1.e05
            tmax(1)=1.e05

c	Read in NAPL properties for SZ NAPL case

            if(ico2.eq.-3) then
               read(inpt,*) dennapl, viscnapl
            end if
         elseif(iflg.eq.-1) then

c     
c     calculate density,compressibility,and viscosity of water
c     
            ndum=neq
            irdofsv=irdof
            irdof=0
            neq=1
            dtot=1.0
            pssv=ps(1)
            ps(1)=1.0
            ieoss=ieos(1)
            iwelbs=iwelb
            iwelb=0
            ssv=s(1)
            irlpsv=irlp(1)
            irlptsv=irlpt(1)
            s(1)=1.0
            iieoss=iieos(1)
            ieos(1)=1
            iieos(1)=1
            phisv=phi(1)
            phi(1)=pref
            t(1)=tref
            irlp(1)=1
            irlpt(1)=1
            dmpfd=dmpf(1)
            dmefd=dmef(1)
            dqd=dq(1)
            call thermw(0)
            irdof=irdofsv
            neq=ndum
            ieos(1)=ieoss
            iieos(1)=iieoss
            phi(1)=phisv
            ps(1)=pssv
            s(1)=ssv
            iwelb=iwelbs
            irlp(1)=irlpsv
            irlpt(1)=irlptsv
c     reference density is in crl(1,1)
c     reference liquid viscosity is in crl(2,1)
c     reference compressibility is in crl(3,1)
c     reference pressure in crl(4,1)
c     reference air viscosity is in crl(5,1)
            crl(1,1)=rolf(1)
            crl(2,1)=1.0/(dil(1)/rolf(1))
            crl(3,1)=dmpf(1)/rolf(1)
c gaz 110819 pref, tref (global) read in scanin              
c            crl(4,1)=pref
            crl(5,1)=182.e-07
c            crl(6,1)=tref
            if(ico2d.eq.-1) then
               crl(7,1)=1.0
            else
               crl(7,1)=0.0
            endif

            dmpf(1)=dmpfd
            dmef(1)=dmefd
            dq(1)=dqd
c     
c     make thernmal conductivities small
c     
c GAZ 10-23-98 eliminate thermal conductivities for isothermal
c           do i=1,n0
c              thx(i)=1.e-30
c              thy(i)=1.e-30
c              thz(i)=1.e-30
c           enddo
c     
c     initialize 2-phase regions             
c     
            do i=1,n0
               ieos(i)=2
            enddo

         elseif(iflg.eq.1) then
c     
c     determine phase state
c
           if(irdof.ne.13) then
c pnx used to be set in the following loops -- currently removed     
            if (iad.eq.0) strd = 1.0        
               do  mid=1,neq
                  mi=mid+ndummy
c
                  if(s(mi).lt.1.0.and.ieos(mi).eq.1) then
                     if(iad.eq.1) strd = strd_iter
                     strd = strd_iter
                     ieos(mi)=2
c                    call phase_timcrl(days,day,dtot,dtotdm)
c
                  else if(s(mi).gt.0.0.and.ieos(mi).eq.3) then
                     if(iad.eq.1) strd = strd_iter
                     strd = strd_iter
                     ieos(mi)=2
c                    call phase_timcrl(days,day,dtot,dtotdm)
c
                  else if(s(mi).gt.1.0.and.ieos(mi).eq.2) then
                     if(iad.eq.1) strd = strd_iter
                     ieos(mi)=1
c                    call phase_timcrl(days,day,dtot,dtotdm)
c
                  else if(s(mi).le.0.0.and.ieos(mi).eq.2) then
                     if(iad.eq.1) strd = strd_iter
                     s(mi)=0.0
                     strd = strd_iter
                     ieos(mi)=3
c                    call phase_timcrl(days,day,dtot,dtotdm)
c
                  endif
               enddo
               do  mid=1,neq
                  mi=mid+ndummy
 
                  if(phi(mi).lt.1.d-2) then
                     if(iad.eq.1) strd = strd_iter
                     phi(mi)=1.d-2
                     if (so(mi).ge.1.0) then
                      s(mi)= 0.999
                      phi(mi)=5.d-2
                      strd = strd_iter
                     endif
                  endif
c
                  if(s(mi).lt.0.0) then
                    if(iad.eq.1) strd = strd_iter
                     s(mi)=0.0
c                    strd = strd_iter
                  endif
c
                  if(s(mi).gt.1.0) then
                     if(iad.eq.1) strd = strd_iter
c                    s(mi)=1.0
c                    strd = strd_iter
                  endif
c
               enddo
           endif
c gaz taken out 4-11-2001 
c              if(ihead.ne.0) then
c               do  mid=1,neq
c                 mi=mid+ndummy
c                 if(ka(mi).lt.0.and.phi(mi).lt.crl(4,1)) then
c                    pflow(mi)=crl(4,1)
c                    wellim(mi)=1.e-2
c                    esk(mi)=-1.0             
c                 endif
c               enddo
c              endif
         elseif(iflg.eq.2) then
          if(abs(irdof).ne.14) then
c     
c           update solution
c     
c       called externally for each zone (gaz 11-20-2003)
c          do zone = 1, 2

            nr1=nrhs_part(zone,1)
            nr2=nrhs_part(zone,2)
	      bp = 0.0

c           strd is passed through common
            do i=1,neq_part(zone)
               i1=i+nr1
               i2=i+nr2
               mi = index_part(zone,i)
               strd_part= repeat_part(mi)
               phi(mi)= phi(mi)-bp_part(zone,i1)*strd*
     &             strd_part
               s(mi)= s(mi)-bp_part(zone,i2)*strd*
     &             strd_part

               bp(index_part(zone,i)+nrhs(1)) = 
     &          + bp(index_part(zone,i)+nrhs(1))
     &          + bp_part(zone,i1)*strd*strd_part
               bp(index_part(zone,i)+nrhs(2)) = 
     &          + bp(index_part(zone,i)+nrhs(2))
     &          + bp_part(zone,i2)*strd*strd_part
            enddo
c
c turn off initial parent BC if necessary
c

c          enddo
c
c call cascade redistribution if requested
c
           if(iflux_ts.ne.0) then
            call cascade_sat(1)
           endif
c
          else
c     
c           update solution(with exact mass balance)
c     
            nr1=nrhs(1)
            nr2=nrhs(2)
c           strd is passed through common
            do i=1,neq
               i1=i+nr1
               i2=i+nr2
               phi(i)=phi(i)-bp(i1)*strd
            enddo
c
c           call wellrate(0,0) 
c           call thrair(0)
            call interblock_iso(0) 
          
          endif
         elseif(iflg.eq.3) then
            if(ico2.eq.-3) then
c     call isothermal SZ NAPL thermodynamics
               call thrsznapl(ndummy)
            else
c     call isothermal air water thermodynamics
               if(ifree.ne.0) then
                 call wtsictr(2)
                 call thrair(ndummy)
                 call wtsictr(4)
               else
                 call thrair(ndummy)
               endif
            end if
         elseif(iflg.eq.4) then
c     call equation generation and load a array

c    gensl2 must be passed the zone

            call gensl2_part(zone)

        elseif(iflg.eq.-4) then
c solve explicit equations
         call gensl2_explicit       
c
         elseif(iflg.eq.5) then
            if(ntty.eq.2) then
c     
c     output for air
c     
               if (iout .ne. 0) write(iout,803)
               if (iatty .ne. 0) write(iatty,803)
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
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
                     if (iout .ne. 0) write(iout,702) il
                     if (iatty .gt. 0) write(iatty,702) il
 702                 format(2x,'Matrix Level = ',i1)
                  endif
                  do i=1,mlev
                     md=  nskw(i+(il-1)*mlev)
                     rqd= sk(md)
                     qcd=qh(md)
                     if(ihead.ne.0) then
                        if (iout .ne. 0) write(iout,804) 
     &                       md,max(phi(md)-phi_inc,0.1d00),
     &                       0.d00,max(phi(md)-phi_inc,0.1d00),qcd
                        if(iatty.ne.0)
     &                       write(iatty,804) 
     &                       md,max(phi(md)-phi_inc,0.1d00),
     &                       0.d00,max(phi(md)-phi_inc,0.1d00),qcd
                     else
                        if (iout .ne. 0) write(iout,804) 
     &                       md,phi(md),pcp(md),
     &                       phi(md)-pcp(md),qcd
                        if(iatty.ne.0)
     &                       write(iatty,804) md,phi(md),pcp(md),
     &                       phi(md)-pcp(md),qcd
                     endif
                  enddo
               enddo
 803           format(/,20x,'Nodal Information (Vapor)',/,2x,'Node',3x,
     *              'air(vp) pres',3x,'cap pres ',3x,'liquid pres '
     *              ,3x,'source/sink(kg/s)')
 804           format(i6,5x,f8.3,1x,f10.3,4x,f10.3,5x,g12.4)
c     calculate global mass and energy flows
               if (iout .ne. 0) then
                  write(iout,703) 
                  write(iout,704) qtotei
                  write(iout,705) qtote
                  write(iout,706) qte
                  write(iout,707) dife
               end if
               if(iatty.ne.0) then
                  write(iatty,703) 
                  write(iatty,704) qtotei
                  write(iatty,705) qtote
                  write(iatty,706) qte
                  write(iatty,707) dife
               endif
 703           format(/,20x,'Global Mass Balances (Vapor)')
 704           format(1x,'Vapor discharge this time step: ',e14.6,' kg')
 705           format(1x,'Total vapor discharge: ',9x,e14.6,' kg')
 706           format(/,1x,'Net kg vapor discharge (total out-total ',
     &              'in): ',e14.6)
 707           format(1x,'Conservation Error: ',25x,e14.6)
            endif
         elseif(iflg.eq.6) then
c     
c     store sk in qc
c     initialize t,to,tini,iieos
c
c gaz 110919 tref now global        
c            tref=crl(6,1)
            do i=1,n
               qc(i)=sk(i)
               if(to(i).eq.0.0d00) then
                 t(i)=tref
                 to(i)=tref
                 tini(i)=tref
               endif
               iieos(i)=1
               if(ieos(i).eq.1) then
                  ieos(i)=2
                  s(i)=1.0
                  so(i)=1.0
               endif
             if(abs(irdof).eq.14) then
              denj(i)=1.0-s(i)
             endif
            enddo
c
c check also for free surface calcs
c
           if(ifree.ne.0) then
            call wtsictr(1)
c
           endif
         elseif(iflg.eq.7) then
c
c convert from head to pressure
c
          call headctr(3,0,0.0,0.0)

         elseif(iflg.eq.8) then
c
c convert from pressure to head
c
          call headctr(2,0,0.0,0.0)
         elseif(iflg.eq.9) then
c
c convert from head to pressure boundary value
c
          call headctr(6,0,0.0,0.0)

         elseif(iflg.eq.10) then
c
c convert from head to pressure boundary and initial 
c conditions (based on pair = pref at max height)
c correct for negative pressures
c
          call head_2phase(0)           

         endif
      endif

      return
      end
  


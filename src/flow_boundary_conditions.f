       subroutine flow_boundary_conditions(iflg)
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
C***********************************************************************
CD1 
CD1 PURPOSE
CD1 
CD1 To read and apply boundary and source conditions for FEHM   
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2
CD2 Initial implementation: NOV-96, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/flow_boundary_conditions.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:04:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
C***********************************************************************
CD3 
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 2.3.7 Sources and sinks
CD3 2.6   Provide Input/Output Data Files
CD3 3.0   INPUT AND OUTPUT REQUIREMENTS
CD3 
C***********************************************************************
CD4 
CD4 SPECIAL COMMENTS AND REFERENCES
CD4 
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4 
C***********************************************************************
c
c flow boundary condition input including time dependent
c and distributed (by volume) sources
c started george zyvoloski june 94
c finally finishing in nov 96
c       

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comii
      use comki
      use comxi
      use davidi

      implicit none

      integer iflg, i, mmodel_old
      integer iimodel
      integer, allocatable :: idum(:)
 
      real*8 days0,qair,qwat
C       real*8 sdum,pdum,denhdum,denehdum
      real*8 denhold,denehold,diffmass,diffener
c
      if(iboun.eq.0) return
c
      if(iflg.eq.0) then
c
c zero out memory and set memory as appropriate
c
         maxtimes=maxtimes+1
c
         if (.not. allocated (node_model)) allocate(node_model(n0))
         node_model = 0
         if (.not. allocated (min_model)) 
     &		allocate(min_model(maxmodel))
          min_model = 0
         if (.not. allocated (vtotw)) allocate(vtotw(maxmodel))
         if (.not. allocated (vtota)) allocate(vtota(maxmodel))
         if (.not. allocated (vtote)) allocate(vtote(maxmodel))
         if (.not. allocated (atotd)) allocate(atotd(maxmodel))
         if (.not. allocated (time)) allocate(time(maxtimes,maxmodel))
         if (.not. allocated (time_cycle)) 
     .        allocate(time_cycle(maxmodel))
         if (.not. allocated (time_type ))allocate(time_type(maxmodel))
	   if (.not. allocated (tunit_type ))
     &	    allocate(tunit_type(maxmodel))
         if (.not. allocated (time_interpolate ))
     &        allocate(time_interpolate(maxmodel))
         vtotw = 0.0d00
         vtota = 0.0d00
         vtote = 0.0d00
         time = 0.0d00
         time_cycle = 0.0d00
         time_type        = 0
         time_interpolate = 0
         tunit_type = 3
         fac_sec_days = 1./86400.
         fac_min_days = 1./1440.
         fac_year_days = 365.25
         if(iqa.ne.0) then
            if (.not. allocated (qa)) allocate(qa(n0))                  
            if (.not. allocated (sourcea)) 
     .           allocate(sourcea(maxtimes,maxmodel))   
            if (.not. allocated (sourcea_type)) 
     .           allocate(sourcea_type(maxmodel))   
            qa = 0.0d00
            sourcea = 0.0d00
            sourcea_type = 0
         endif
         if(iqf.ne.0) then
            if (.not. allocated (qw)) allocate(qw(n0))                  
            if (.not. allocated (qw0)) allocate(qw0(n0)) 
            if (.not. allocated (sourcef)) 
     2           allocate(sourcef(maxtimes,maxmodel))   
            if (.not. allocated (sourcef_type)) 
     2           allocate(sourcef_type(maxmodel))   
            qw = 0.0d00
            qw0 = 0.0d00
            sourcef = 0.0d00
            sourcef_type = 0
         endif
         if(ifd.ne.0) then
            if (.not. allocated (drain)) allocate(drain(n0))                  
            if (.not. allocated (drainar)) 
     .           allocate(drainar(maxtimes,maxmodel))   
            if (.not. allocated (drainar_type)) 
     .           allocate(drainar_type(maxmodel))   
            drain = 0.0d00
            seepfac = 0.0d00
            seepfac_type = 0
         endif
         if(isf.ne.0) then
            if (.not. allocated (sp)) allocate(sp(n0))                  
            if (.not. allocated (seepfac)) 
     .           allocate(seepfac(maxtimes,maxmodel))   
            if (.not. allocated (seepfac_type)) 
     .           allocate(seepfac_type(maxmodel))   
            sp = 0.0d00
            seepfac = 0.0d00
            seepfac_type = 0
         endif
         if(iqw.ne.0) then
            if (.not. allocated (qw)) allocate(qw(n0))                  
            if (.not. allocated (sourcew)) 
     .           allocate(sourcew(maxtimes,maxmodel))   
            if (.not. allocated (sourcew_type)) 
     .           allocate(sourcew_type(maxmodel))   
            qw = 0.0d00
            sourcew = 0.0d00
            sourcew_type = 0
         endif
         if(iqenth.ne.0) then
            if (.not. allocated (qenth)) allocate(qenth(n0))
            if (.not. allocated (sourcee)) 
     .           allocate(sourcee(maxtimes,maxmodel))   
            if (.not. allocated (sourcee_type)) 
     .           allocate(sourcee_type(maxmodel))   
            qenth = 0.0d00
            sourcee = 0.0d00
            sourcee_type = 0
         endif
         if(ienth.ne.0) then
            if (.not. allocated (enth)) allocate(enth(n0))
            if (.not. allocated (enthalpy)) 
     .           allocate(enthalpy(maxtimes,maxmodel))   
            if (.not. allocated (enthalpy_type)) 
     .           allocate(enthalpy_type(maxmodel))   
            enth = 0.0d00
            enthalpy = 0.0d00
            enthalpy_type = 0
         endif
         if(isatb.ne.0) then
            if (.not. allocated (satb)) allocate(satb(n0))
            if (.not. allocated (saturation)) 
     .           allocate(saturation(maxtimes,maxmodel))   
            if (.not. allocated (saturation_type)) 
     .           allocate(saturation_type(maxmodel))   
            satb = 0.0d00
            saturation = 0.0d00
            saturation_type = 0
         endif
         if(itempb.ne.0) then
            if (.not. allocated (tempb)) allocate(tempb(n0))
            if (.not. allocated (temperature)) 
     .           allocate(temperature(maxtimes,maxmodel))   
            if (.not. allocated (temperature_type)) 
     .           allocate(temperature_type(maxmodel))   
            tempb = 0.0d00
            temperature = 0.0d00
            temperature_type = 0
         endif
         if(ipresa.ne.0) then
            if (.not. allocated (presa)) allocate(presa(n0))
            if (.not. allocated (pressurea)) 
     .           allocate(pressurea(maxtimes,maxmodel))   
            if (.not. allocated (pressurea_type)) 
     .           allocate(pressurea_type(maxmodel))   
            presa = 0.0d00
            pressurea = 0.0d00
            pressurea_type = 0
         endif
         if(ipresw.ne.0) then
            if (.not. allocated (presw)) allocate(presw(n0))
            if (.not. allocated (pressurew)) 
     .           allocate(pressurew(maxtimes,maxmodel))   
            if (.not. allocated (pressurew_type)) 
     .           allocate(pressurew_type(maxmodel))   
            presw = 0.0d00
            pressurew = 0.0d00
            pressurew_type = 0
         endif
         if(imped.ne.0) then
            if (.not. allocated (impedance)) 
     .           allocate(impedance(maxtimes,maxmodel))
            if (.not. allocated (impedance_type)) 
     .           allocate(impedance_type(maxmodel))
         endif
c new initial value capability 12/03/98 GAZ
         if(ipresw_ini.ne.0) then
            if (.not. allocated (presw_ini)) allocate(presw_ini(n0))
            if (.not. allocated (pressurew_ini)) 
     .           allocate(pressurew_ini(maxtimes,maxmodel))   
            if (.not. allocated (pressurew_ini_type)) 
     .           allocate(pressurew_ini_type(maxmodel))   
            presw_ini = 0.0d00
            pressurew_ini = 0.0d00
            pressurew_ini_type = 0
         endif
         if(ipresa_ini.ne.0.or.isatb_ini.ne.0) then
            if (.not. allocated (presa_ini)) allocate(presa_ini(n0))
            if (.not. allocated (pressurea_ini)) 
     .           allocate(pressurea_ini(maxtimes,maxmodel))   
            if (.not. allocated (pressurea_ini_type)) 
     .           allocate(pressurea_ini_type(maxmodel))   
            presa_ini = 0.0d00
            pressurea_ini = 0.0d00
            pressurea_ini_type = 0
         endif
         if(ipresa_ini.ne.0.or.isatb_ini.ne.0) then
            if (.not. allocated (satb)) allocate(satb_ini(n0))
            if (.not. allocated (saturation_ini)) 
     .           allocate(saturation_ini(maxtimes,maxmodel))   
            if (.not. allocated (saturation_ini_type)) 
     .           allocate(saturation_ini_type(maxmodel))   
            satb_ini = 0.0d00
            saturation_ini = 0.0d00
            saturation_ini_type = 0
         endif
         if(itempb_ini.ne.0) then
            if (.not. allocated (tempb_ini)) allocate(tempb_ini(n0)) 
            if (.not. allocated (temperature_ini)) 
     .           allocate(temperature_ini(maxtimes,maxmodel))   
            if (.not. allocated (temperature_ini_type)) 
     .           allocate(temperature_ini_type(maxmodel))   
            tempb_ini = 0.0d00
            temperature_ini = 0.0d00
            temperature_ini_type = 0
         endif
         if(ixperm.ne.0) then
            if (.not. allocated (permx)) 
     &           allocate(permx(maxtimes,maxmodel))
            if (.not. allocated (permx_type))
     &           allocate(permx_type(maxmodel))
            permx = 0.0d00
            permx_type = 0
         endif
         if(iyperm.ne.0) then
            if (.not. allocated (permy)) 
     &           allocate(permy(maxtimes,maxmodel))
            if (.not. allocated (permy_type))
     &           allocate(permy_type(maxmodel))
            permy = 0.0d00
            permy_type = 0
         endif
         if(izperm.ne.0) then
            if (.not. allocated (permz)) 
     &           allocate(permz(maxtimes,maxmodel))
            if (.not. allocated (permz_type))
     &           allocate(permz_type(maxmodel))
            permz = 0.0d00
            permz_type = 0
         endif
         if(iwght.ne.0) then
            if (.not. allocated (weight_type)) 
     .           allocate(weight_type(maxmodel))   
            weight_type = 0
         endif
         if(isty.ne.0) then
            if (.not. allocated (steady_type)) 
     .           allocate(steady_type(maxmodel))   
         else
            if (.not. allocated (steady_type)) 
     .           allocate(steady_type(1))
         endif
         steady_type = 1
c always allocate space for timestep
         its = 1
         if (.not. allocated (timestep))
     .        allocate(timestep(maxtimes,maxmodel))   
         if (.not. allocated (timestep_type))
     .        allocate(timestep_type(maxmodel))   
         timestep_type = 0
         timestep = 0.0d00
         if(icm.ne.0) then
            if (.not. allocated (node_ch_model)) 
     .           allocate(node_ch_model(maxtimes,maxmodel))   
            if (.not. allocated (node_ch_model_type)) 
     .           allocate(node_ch_model_type(maxmodel))   
            node_ch_model = 0
            node_ch_model_type = 0
         endif
c
c set model number to zero
         mmodel=0
         
      else if(iflg.eq.1) then
c
c read boundary conditions
c        
c see if called before by checking node_model
c
         mmodel_old=mmodel
         call flow_boun(0
     &        ,n,ico2,idof
     &        ,days0,days,day,daynew,sx1
     &        ,inpt,iptty,iout,ierr,l,igrav,ihead)
c
c 
c read in model type for each node
c
         allocate(idum(n0))
         idum = 0
         narrays = 1
         itype(1) = 4
         default(1) = 0
         macro = "boun"
         igroup = 1
         call initdata2( inpt, ischk, n0, narrays,
     2        itype, default, macroread(7), macro, igroup, ireturn,
     3        i4_1=idum(1:n0) )
         macroread(7) = .TRUE.
         do i=1,n0
            if(idum(i).ne.0) then
               node_model(i)=idum(i)+mmodel_old
            endif
         enddo
         deallocate(idum)
         
      else if(iflg.eq.2) then
c
c consistency checks
c      
         call flow_boun(1
     &        ,n,ico2,idof
     &        ,days0,days,day,daynew,sx1
     &        ,inpt,iptty,iout,ierr,l,igrav,ihead)
c
      else if(iflg.eq.3) then
c
c change boundary conditions with time
c      
c change BC
c
         days0=max(days-day,0.0d00)
         allocate(idum(n0))
c
c find which nodes are associated with boun macro
c
         idum = 0
         do i = 1,neq
            iimodel = node_model(i)  
            if(iimodel.ne.0) then
               if(icm.ne.0) then
                  if(node_ch_model_type(iimodel).gt.0) then
                     idum(i)=node_ch_model(abs(time_type(iimodel)),
     &                    iimodel)
                  else
                     idum(i) = iimodel
                  endif
               else
                  idum(i) = iimodel           
               endif
            else
               idum(i) = iimodel           
            endif
         enddo
c
c load correct array of physics module
         if(idof.eq.1) then
c heat conduction
            call flow_boun(3
     &           ,n,ico2,idof
     &           ,days0,days,day,daynew,sx1
     &           ,inpt,iptty,iout,ierr,l,igrav,ihead)
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqenth.ne.0) then
                     if(qenth(i).ne.0.0) then
                        qflux(i)=qenth(i)
                        qflxm(i)= 0.0
                     endif
                  else if(itempb.ne.0) then
                     if(tempb(i).gt.0.0) then
                        qflux(i)=tempb(i)
                        qflxm(i)= sx1(i)
                        if(wellim(i).ne.0.0) qflxm(i)=wellim(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(itempb_ini.ne.0) then
                     if(tempb_ini(i).gt.0.0d00) then
                        t(i) = tempb_ini(i)
                        to(i) = t(i)
                        tempb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
         else if(ico2.eq.0) then
c pure water/heat
            call flow_boun(3
     &           ,n,ico2,idof
     &           ,days0,days,day,daynew,sx1
     &           ,inpt,iptty,iout,ierr,l,igrav,ihead)
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqenth.ne.0) then
                     if(qenth(i).ne.0.0) then
                        qflux(i)=qenth(i)
                        qflxm(i)= 0.0
                     endif
                  else if(itempb.ne.0) then
                     if(tempb(i).gt.0.0) then
                        qflux(i)=tempb(i)
                        qflxm(i)= sx1(i)
                        if(wellim(i).ne.0.0) qflxm(i)=wellim(i)
                     endif
                  else if(isatb.ne.0) then
                     if(satb(i).gt.0.0) then
                        qflux(i)=satb(i)
                        qflxm(i)= -sx1(i)
                        if(wellim(i).ne.0.0) qflxm(i)=-wellim(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw_ini.ne.0) then
                     if(presw_ini(i).gt.0.0d00) then
                        phi(i) = presw_ini(i)
                        pho(i) = phi(i)
                        presw_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(itempb_ini.ne.0) then
                     if(tempb_ini(i).gt.0.0) then
                        t(i) = tempb_ini(i)
                        to(i) = t(i)
                        tempb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(isatb_ini.ne.0) then
                     if(satb_ini(i).gt.0.0) then
                        s(i) = satb_ini(i)
                        so(i) = s(i)
                        ieos(i) = 2
                        if(s(i).le.0.0d00) ieos(i) = 1
                        if(s(i).ge.1.0d00) ieos(i) = 3
                        satb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw.ne.0) then
                     if(presw(i).ne.0.0) then
                        pflow(i)=presw(i)
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        ka(i)=-1
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresa.ne.0) then
                     if(presa(i).ne.0.0) then
                        pflow(i)=presa(i)
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        ka(i)=-1
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ienth.ne.0) then
                     if(enth(i).ne.0.0) then
                        esk(i)=enth(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqw.ne.0.or.iqf.ne.0) then
                     if(qw(i).ne.0.0) then
                        sk(i)=qw(i)
                        ka(i)=1
                     endif
                  endif
               endif
            enddo
         else if(ico2.gt.0) then
c air/water/heat
            call flow_boun(3
     &           ,n,ico2,idof
     &           ,days0,days,day,daynew,sx1
     &           ,inpt,iptty,iout,ierr,l,igrav,ihead)
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqenth.ne.0) then
                     if(qenth(i).ne.0.0) then
                        qflux(i)=qenth(i)
                        qflxm(i)= 0.0
                     endif
                  else if(itempb.ne.0) then
                     if(tempb(i).gt.0.0) then
                        qflux(i)=tempb(i)
                        qflxm(i)= sx1(i)
                        if(wellim(i).ne.0.0) qflxm(i)=wellim(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw.ne.0) then
                     if(presw(i).ne.0.0) then
                        pflow(i)=presw(i)
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(isatb.ne.0.0) then
                     if(satb(i).gt.0.0) then
                        pflow(i)=-satb(i)
                        ka(i)=-1
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                     else if(satb(i).lt.0.0) then
                        eskc(i)=-satb(i)
                        ka(i)=-1
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresa.ne.0.0) then
                     if(presa(i).ne.0.0) then
                        eskc(i)=-presa(i)
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ienth.ne.0.0) then
                     if(enth(i).ne.0.0) then
                        esk(i)=enth(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw_ini.ne.0) then
                     if(presw_ini(i).gt.0.0) then
                        phi(i) = presw_ini(i)
                        pho(i) = phi(i)
                        presw_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresa_ini.ne.0) then
                     if(presa_ini(i).gt.0.0) then
                        pci(i) = presa_ini(i)
                        pcio(i) = pci(i)
                        presa_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(itempb_ini.ne.0) then
                     if(tempb_ini(i).gt.0.0) then
                        t(i) = tempb_ini(i)
                        to(i) = t(i)
                        tempb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(isatb_ini.ne.0) then
                     if(satb_ini(i).gt.0.0) then
                        s(i) = satb_ini(i)
                        so(i) = s(i)
                        ieos(i) = 2
                        if(s(i).le.0.0d00) ieos(i) = 1
                        if(s(i).ge.1.0d00) ieos(i) = 3
                        satb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqw.ne.0.or.iqf.ne.0) then
                     if(qw(i).ne.0.0) then
                        sk(i)=qw(i)
                        ka(i)=1
                     endif
                  endif
               endif
            enddo
c       else if (ico2 .lt. 0) then
c 04/06/2005 MH
 	 else if (ico2 .lt. 0 .and. ice .ne. 2) then
c 04/06/2005 MH
c isothermal air/water
            call flow_boun(3
     &           ,n,ico2,idof
     &           ,days0,days,day,daynew,sx1
     &           ,inpt,iptty,iout,ierr,l,igrav,ihead)
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresa.ne.0.0) then
                     if(presa(i).ne.0.0) then
c zero pres indicated by -9999999.0
                        if(presa(i).eq.-9999999.0) then
                           pflow(i)=0.
                        else
                           pflow(i)=presa(i)
                        endif
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        if(satb(i).eq.1.0) esk(i)=1.0
                        if(satb(i).eq.0.0) esk(i)=-1.000
                        ka(i)=-1
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw.ne.0.0) then
                     if(presw(i).ne.0.0) then
c zero pres indicated by -9999999.0
                        if(presw(i).eq.-9999999.0) then
                           pflow(i)=0.
                        else
                           pflow(i)=presw(i)
                        endif
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        esk(i)=1.0    
                        ka(i)=-1
                     endif
                  endif
               endif
            enddo
c set seepage face conditions
            do i=1,n
               if(idum(i).ne.0) then
                  if(isf.ne.0.0) then
                     if(sp(i).ne.0.0) then
c zero pres indicated by -9999999.0
                        if(sp(i).eq.-9999999.0) then
                           pflow(i)=0.
                        else if(sp(i).lt.0.0) then
                           call headctr(5,i,pflow(i),abs(sp(i)))
                        else
                           pflow(i)=sp(i)
                        endif
                        if(sp(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) then 
                           wellim(i)=-sx1(i)*1.d06
                        else
                           wellim(i)=-abs(wellim(i))
                        endif
c           esk(i)=1.0    
                        ka(i)=-3
                     endif
                  endif
               endif
            enddo
c set free drainage conditions
            do i=1,n
               if(idum(i).ne.0) then
                  if(ifd.ne.0.0) then
                     if(drain(i).ne.0.0) then
c zero pres indicated by -9999999.0
                        if(drain(i).eq.-9999999.0) then
                           pflow(i)=0.
                        else
                           pflow(i)=drain(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        esk(i)=1.0    
                        ka(i)=-2
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(itempb.ne.0.0) then
                     if(tempb(i).ne.0.0) then
c convert head to pressure(head resides in the tempb array)
c zero head indicated by -9999999.0
                        if(tempb(i).eq.-9999999.0) then
                           pflow(i)=0.
                        else
                           pflow(i)=tempb(i)
                        endif
                        call headctr(5,i,pflow(i),pflow(i))
                        if(irdof.ne.13) then
                           if(pflow(i).lt.crl(4,1)) then
c GAZ(4-8-2001)            pflow(i)=crl(4,1)
c GAZ(4-8-2001)            satb(i)=1.e-6
                           else
                              satb(i)=1.0
                           endif
                        endif
                        ka(i)=-1
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(node_model(i).gt.0) then
                     if(ipresa_ini.ne.0) then
                        if(presa_ini(i).gt.0.0) then
                           phi(i) = presa_ini(i)
                           pho(i) = presa_ini(i)
                        endif
                     endif
                     if(isatb_ini.ne.0) then
                        if(satb_ini(i).gt.0.0) then
                           s(i) = satb_ini(i)
                           so(i) = satb_ini(i)
                        endif
                     endif
                     if(isatb_ini.ne.0.or.ipresa_ini.ne.0) then
                        if(satb_ini(i).gt.0.0.or.presa_ini(i).gt.0.0) 
     &                       then
                           denhold = denh(i)
                           denehold = deneh(i)
                           call accum_term(i,0)
                           diffmass = (denh (i)-denhold)*volume(i)
                           diffener = (deneh (i)-denehold)*volume(i)
                           am0=am0+diffmass
                           ame=ame+diffener
                           node_model(i)=-abs(node_model(i))
                        endif
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  qair=0.0
                  qwat=0.0
                  if(iqa.ne.0) then
                     if(qa(i).ne.0.0) then
                        qair =qa(i)
                        ka(i)=1
                     endif
                  endif
                  if(iqw.ne.0.or.iqf.ne.0) then
                     if(qw(i).ne.0.0) then
                        sk(i)=qw(i)
                        qwat =qw(i)
                        ka(i)=1
c     set specified air pres when flowing
                        if(ipresa.ne.0) then
                           if(presa(i).gt.0.0) esk(i)=-presa(i)
                           if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        endif
                     endif
                  endif
                  if(isatb.gt.0.0) then
                     if(satb(i).gt.0.0) then
                        esk(i)=satb(i)
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        ka(i)=-1
                     endif
                  endif
                  if(qair+qwat.ne.0) then
                     if(esk(i).ge.0.0) esk(i)=qwat/(qwat+qair)
                     qc(i)=qair+qwat
                  endif
               endif
            enddo
            
	 else if(ico2 .lt. 0 .and. ice .eq. 2 ) then
c     meth tenma 4/6/2005
            call flow_boun(3
c     & ,n,ico2,idof
     &           ,n,0,idof
     &           ,days0,days,day,daynew,sx1
     &           ,inpt,iptty,iout,ierr,l,igrav,ihead)
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqenth.ne.0) then
                     if(qenth(i).ne.0.0) then
                        qflux(i)=qenth(i)
                        qflxm(i)= 0.0
                     endif
                  else if(itempb.ne.0) then
                     if(tempb(i).gt.0.0) then
                        qflux(i)=tempb(i)
                        qflxm(i)= sx1(i)
                        if(wellim(i).ne.0.0) qflxm(i)=wellim(i)
                     endif
                  else if(isatb.ne.0) then
                     if(satb(i).gt.0.0) then
                        qflux(i)=satb(i)
                        qflxm(i)= -sx1(i)
                        if(wellim(i).ne.0.0) qflxm(i)=-wellim(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw_ini.ne.0) then
                     if(presw_ini(i).gt.0.0d00) then
                        phi(i) = presw_ini(i)
                        pho(i) = phi(i)
                        presw_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(itempb_ini.ne.0) then
                     if(tempb_ini(i).gt.0.0) then
                        t(i) = tempb_ini(i)
                        to(i) = t(i)
                        tempb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(isatb_ini.ne.0) then
                     if(satb_ini(i).gt.0.0) then
                        s(i) = satb_ini(i)
                        so(i) = s(i)
                        ieos(i) = 2
                        if(s(i).le.0.0d00) ieos(i) = 1
                        if(s(i).ge.1.0d00) ieos(i) = 3
                        satb_ini(i) = 0.0d00
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresw.ne.0) then
                     if(presw(i).ne.0.0) then
                        pflow(i)=presw(i)
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        ka(i)=-1
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ipresa.ne.0) then
                     if(presa(i).ne.0.0) then
                        pflow(i)=presa(i)
                        if(pflow(i).eq.0.0) then
                           pflow(i)=phini(i)
                        endif
                        if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
                        ka(i)=-1
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(ienth.ne.0) then
                     if(enth(i).ne.0.0) then
                        esk(i)=enth(i)
                     endif
                  endif
               endif
            enddo
            do i=1,n
               if(idum(i).ne.0) then
                  if(iqw.ne.0.or.iqf.ne.0) then
                     if(qw(i).ne.0.0) then
                        sk(i)=qw(i)
                        ka(i)=1
                     endif
                  endif
               endif
            enddo
c     meth 4/6/2005
         endif           
c     added call to adjust boundary heads if necessary(GAZ 4-8-2001)
         if(ihead.ne.0) call airctr(10,0)
         deallocate(idum)
      elseif(iflg.eq.4) then
c     
c     set models for transient simulation
c     
         if(isty.ne.0) then
            lchange = 0
            do i=1,mmodel
               steady_type(i) = abs(steady_type(i))
            enddo
         endif
      endif      
c     
      return
      end

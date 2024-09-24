      subroutine hydrate_equil(iflg,itype)
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
!D1 To calculate  hydrate equilibrium properties and derivatives.
!D1 Modify exixting terms for equilibrium
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 01-august-04, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/hydrate_equilibrium.f_a  $
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

      use comai
      use comei
      use davidi
      use commeth
      implicit none
      integer iflg,iphase,itype,i,j,k,ihyddc
      real*8 prop,der1,der2,der3,der4,w_frac,fracg,frach,tdis1
      real*8 pl,tl
      real*8 dum5,dum6,dum7,dum8
      real*8 gas_const,temp_conv,v_tol,vard
      parameter (gas_const = 8.314, temp_conv = 273.15, v_tol = 1.e-6)

c input ----------
c iflag - designator for type of calculation
c iphase - phase state                                   
c var1 - variable 1(pressure)
c var2 - variable 2(temperature)
c var3 - variable 3(water fraction)
c var4 - variable 4(hydrate fraction)
c output ----------
c prop - property 
c der1 - derivative of property wrt variable 1
c der2 - derivative of property wrt variable 2
c der3 - derivative of property wrt variable 3
c der4 - derivative of property wrt variable 4
c -----------------

c iflg=1, phase state calculations
c iflg=2, modify accumulation terms
c iflg=3, variable switching 
c iflg=4, remove rate terms (set = 0)             

      if(iflg.eq.1) then
c phase state and independent variable calculation
c assumes icectr has been called for hydrate phase (icectr(-6,0)
         do i = 1,n
            if(ihyd(i).lt.0) then
c temperature is less than the hydrate disociation temperature 
c independent variables are p,t,fracw or fracg, Sh(function other variables)
c    calculate hydrate fraction
            
               fracg = fracgas(i)
               w_frac = fracw(i)
               tl = tmeth(i)
               pl= phimeth(i)
               ihyddc = ihyd(i) 
               call hydrate_properties(6,0,pl,0.0d00,0.0d00,0.0d00,
     &              tdis1,der1,0.0d00,0.0d00,0.0d00)
               call hydrate_properties(8,ihyddc,pl,tl,w_frac,fracg,
     &              frach,dum5,dum6,der3,der4)
           
               if(tl.le.tdis1) then
                  if(ihyddc.eq.-1) then
                     a(nmat(25)+i) = der3
                  else
		     a(nmat(25)+i) = der4
                  endif
               endif
            else if(ihyd(i).eq.0) then
c pressure is at the hydrate forming pressure  
c independent variables are p,fracw,Sh
c t is calculated 
               pl = phimeth(i)
               tl = tmeth(i)
               call hydrate_properties(6,0,pl,0.0d00,0.0d00,0.0d00,
     &          tl,der1,0.0d00,0.0d00,0.0d00)
               a(nmat(25)+i) = der1
            else
c pressure is greater than the hydrate forming pressure 
c pressure is less than the hydrate forming pressure  
c independent variables are p,t,fracw,frachyd = const. 
            endif
         enddo
      else if(iflg.eq.2) then
c identify hydrate portion of fluid fractions
         allocate (fracw_tot(n))
         allocate (fracg_tot(n))
         do i = 1,n
            if(frachyd(i).gt.0.0d0) then
c calculate mobile fluid phases
               w_frac = fracw(i) - fracl*frachyd(i)
               fracg = fracgas(i) - fracv*frachyd(i)
               fracw_tot(i) = fracw(i)
               fracg_tot(i) = fracgas(i)
               fracw(i) = w_frac
               fracgas(i) = fracg
            endif
         enddo

      else if(iflg.eq.3) then
         do i = 1,n
            if(frachyd(i).gt.0) then
c calculate fluid phases
               fracw(i) = fracw_tot(i) 
               fracgas(i) = fracg_tot(i)
            endif
         enddo
         deallocate (fracw_tot)
         deallocate (fracg_tot)
     
      else if(iflg.eq.4) then
c    zero out rate terms 
         if(itype.eq.1) then
c    water
            do i = 1,n
               skwhyd(i) = 0.0
               dskhyd1(i) = 0.0
               dskhyd2(i) = 0.0
               dskhyd3(i) = 0.0
               dskhyd4(i) = 0.0
            enddo
         else if(itype.eq.2) then
c    methane
            do i = 1,n
               skwhyd(i) = 0.0
               dskhyd1(i) = 0.0
               dskhyd2(i) = 0.0
               dskhyd3(i) = 0.0
               dskhyd4(i) = 0.0
            enddo
         else if(itype.eq.3) then
c    hydrate
            do i = 1,n
               skwhyd(i) = 0.0
               dskhyd1(i) = 0.0
               dskhyd2(i) = 0.0
               dskhyd3(i) = 0.0
               dskhyd4(i) = 0.0
               qhhyd(i) = 0.0
               dqhhyd1(i) = 0.0
               dqhhyd2(i) = 0.0
               dqhhyd3(i) = 0.0
               dqhhyd4(i) = 0.0
            enddo

         endif

      endif
      return 
      end

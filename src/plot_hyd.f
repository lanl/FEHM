      subroutine plot_hyd(igf)
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Write to time history plot tape.   
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/plot_hyd.f_a  $
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

      integer i,igf,mi
      real*8 relhum,pwatersat,psat,dpdummy
      real*8 qhdum, skdum, tdum, phidum, pcidum
      real*8 pcpdum, sdum, headdum
      real*8 dumconv, dumconv1, rolconv
      real*8 pdum,fracwdum,frachdum, fracgdum

      if (ishis .gt. 0) then

         if (igf .eq. 0)  then
            ishyd = 102
            open(ishyd,file = 'methhydrate.his',status='unknown')
            write(ishyd,*)
            write(ishyd,*)        
            write(ishyd, '(a4)')  'hyd ' 
            write(ishyd,*) 
            write(ishyd,*) 
c**** write number of plot nodes, node numbers and coordinates ****

            write (ishyd, *)  m
            do i = 1, m
               mi = nskw(i)
               if (mi .gt. neq*2) then
                  mi = mi - neq*2
               else if (mi .gt. neq) then
                  mi = mi - neq
               end if
               write(ishyd, 6010) nskw(i), corz(mi,1), corz(mi,2),
     &              corz(mi,3)
 6010          format(i8,1p,3e16.7)
            end do
            
            write(ishyd, 6016) 
            

 6016       format(/,1x,'Headings for Methane ',/,1x,'node ,', 
     &           'flow enthalpy(Mj/kg), gas flow(kg/s), ',
     &           'gas total(kg),', ' temperature(deg C), ',
     &           'pressure(Mpa), gas fraction, hydrate fraction')
c     initialize accumulation
            allocate(skmeth_tot(m))
            
            skmeth_tot = 0.0    
            
         end if

         write(ishyd, *)  days, '  days'
         
         do i = 1, m
            mi = nskw(i)
            qhdum=max(qhmeth(mi),1.d-20)
            if(abs(skmeth(mi)).lt.1.d-20) then
               skdum = 0.0
               qhdum = 0.0
            else
               skdum=skmeth(mi)
               qhdum = qhdum/skdum
            endif
c     RJP 7/18/04 Put a conditional so that for the initial time step
c     dtotdm is not getting multiplied.
            if (igf.ne.0) then
               
               skmeth_tot(i) = skmeth_tot(i) + skdum*dtotdm
            endif
            tdum=max(tmeth(mi),1.d-20)
            phidum=max(phimeth(mi),1.d-20)
            if(idof_meth.ne.7) then
               fracwdum=max(fracw(mi),1.d-20)
               frachdum=max(frachyd(mi),1.d-20)
               fracgdum =1. - fracwdum- frachdum
            else
               fracwdum=max(fracw(mi)-fracl*frachyd(mi),1.d-20)
               frachdum=max(frachyd(mi),1.d-20)
               fracgdum=max(fracgas(mi)-fracv*frachyd(mi),1.d-20)
            endif
            write(ishyd, 6020)  mi, qhdum , skdum , skmeth_tot(i),
     &           tdum , phidum , fracgdum, frachdum

         end do

 6020    format(i8, 7(1x, g12.5))

      end if

      end

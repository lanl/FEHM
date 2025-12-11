      subroutine plot_co2(igf)
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
      use comco2
      implicit none

      integer i,igf,mi,isco2
      real*8 relhum,pwatersat,psat,dpdummy
      real*8 qhdum, skdum, tdum, phidum, pcidum
      real*8 pcpdum, sdum, headdum
      real*8 dumconv, dumconv1, rolconv
      real*8 pdum,fracwdum,fracldum, fracgdum

      if (ishis .gt. 0) then

         isco2 = 103
         if (igf .eq. 0)  then
            open(isco2,file = 'co2.his',status='unknown')
            write(isco2,*)
            write(isco2,*)        
            write(isco2, '(a4)')  'co2 ' 
            write(isco2,*) 
            write(isco2,*) 
c**** write number of plot nodes, node numbers and coordinates ****

            write (isco2, *)  m
            do i = 1, m
               mi = nskw(i)
               if (mi .gt. neq*2) then
                  mi = mi - neq*2
               else if (mi .gt. neq) then
                  mi = mi - neq
               end if
               write(isco2, 6010) nskw(i), corz(mi,1), corz(mi,2),
     &              corz(mi,3)
 6010          format(i8,1p,3e16.7)
            end do
            
            write(isco2, 6016) 
            

 6016       format(/,1x,'Headings for CO2 ',/,1x,'node ,', 
     &           'flow enthalpy(Mj/kg), gas flow(kg/s), ',
     &           'gas total(kg),', ' temperature(deg C), ',
     &           'pressure(Mpa), gas fraction, hydrate fraction')
c     initialize accumulation
            allocate(skco2_tot(m))
            
            skco2_tot = 0.0    
            
         end if

         write(isco2, *)  days, '  days'
         
         do i = 1, m
            mi = nskw(i)
            qhdum=max(qhco2(mi),1.d-20)
            if(abs(skco2(mi)).lt.1.d-20) then
               skdum = 0.0
               qhdum = 0.0
            else
               skdum=skco2(mi)
               qhdum = qhdum/skdum
            endif
c     RJP 7/18/04 Put a conditional so that for the initial time step
c     dtotdm is not getting multiplied.
            if (igf.ne.0) then
               
               skco2_tot(i) = skco2_tot(i) + skdum*dtotdm
            endif
            tdum=max(tco2(mi),1.d-20)
            phidum=max(phico2(mi),1.d-20)
            fracwdum=max(fw(mi),1.d-20)
            fracldum=max(fl(mi),1.d-20)
            fracgdum=max(fg(mi),1.d-20)
            write(isco2, 6020)  mi, qhdum , skdum , skco2_tot(i),
     &           tdum , phidum , fracwdum, fracldum, fracgdum

         end do

 6020    format(i8, 7(1x, g12.5))

      end if

      end

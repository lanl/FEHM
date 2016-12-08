      subroutine initchem
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
C**********************************************************************
CD1 
CD1 PURPOSE
CD1
CD1 To initialize the chemistry calculations
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2 
CD2 Initial implementation: ?, Programmer: Hari Viswanathan
CD2 new reaction method using chemical components
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/initchem.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:12   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:44   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:32 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2
C**********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.4 Solute-transport equations
CD3  2.4.6 Multiple, interacting solutes
CD3
C**********************************************************************
CD4 
CD4 SPECIAL COMMENTS AND REFERENCES
CD4 
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

      use comrxni
      use comchem
      use comdi
      use comdti

	use comco2, only: icarb, carbon_tracer

      implicit none

c local variables
      integer ic,in,idh,mi,mi2,im
      real*8 newimm



		carbon_tracer=0
      if(icarb.eq.1.and.co2_couple.eq.1)then

c loop over aqueus components

                 do ic = 1,ncpnt

		          if (cpntnam(ic).eq.'C'.or.cpntnam(ic).eq.'H2CO3')then

					carbon_tracer=ic

					exit

				  endif

				end do

	endif


c
c.....Set the initial guesses for the speciation of aqueous 
c.....component concentrations.
      do in = 1,n0
         t(in) = 25
      enddo

      do ic=1,ncpnt
         If(ifxconc(ic).eq.2) then
            do in=1,n0
               cpntsv(ic,in)=1.d-7
            enddo
         Else
            do in=1,n0
               cpntsv(ic,in)=1.d-9
            enddo
         Endif
      enddo

c.....If pH is variable and the pH is entered instead of the TOTAQ 
c.....of H+, calculate TOTH for the initial boundary conditions.

      do idh=1,ncpnt
         if(ifxconc(idh).eq.2) then
            do in=1,n0
               mi = in+(pcpnt(idh)-1)*n0
               If(an(mi).ne.0.d0) then

                  do ic=1,ncpnt
                     mi2 = in+(pcpnt(ic)-1)*n0
                     totaq(ic)=an(mi2)
                  enddo
                  call varph(in)
c.....If negative total component concentrations are permitted 
c.....for other species, update their concentations along with TOTH
c.....at this time
                  do ic=1,ncpnt
                     if(ifxconc(ic).eq.2) then
                       mi2 = in+(pcpnt(ic)-1)*n0
                       an(mi2) = totaq(ic)
                     endif
                  enddo
                  do ic = 1,ncpnt
                     cpntsv(ic,in)=cpnt(ic)
                  enddo
               endif
            enddo            
            exit
         endif
      enddo

c.....Check the necessity of DXCT calculation for each component
      call chckderiv
c.....If pH is variable and the pH is entered instead of the TOTAQ 
c.....of H+, calculate TOTH for the initial boundary conditions.

      do idh=1,ncpnt
         if(ifxconc(idh).eq.2) then
            do in=1,n0
               mi = in+(pcpnt(idh)-1)*n0
               If(cnsk(mi).ne.0.d0) then
                  do ic=1,ncpnt
                     mi2 = in+(pcpnt(ic)-1)*n0
                     totaq(ic)=abs(cnsk(mi2))
                  enddo
                  call varph(in)
c.....If negative total component concentrations are permitted 
c.....for other species, update their concentations along with TOTH
c.....at this time
                  do ic=1,ncpnt
                     if(ifxconc(ic).eq.2) then
                        mi2 = in+(pcpnt(ic)-1)*n0
                        if(cnsk(mi2).gt.0)then
                           cnsk(mi2)=totaq(ic)
                        else
                           cnsk(mi2)=-totaq(ic)
                        endif
                     endif
                  enddo
               endif
            enddo            
c.....This exit is necessary
            exit
         endif
      enddo
      do im = 1,nimm
         do in = 1,n0
            mi = in+(pimm(im)-1)*n0
            newimm = an(mi)
            if(newimm.le.conc_min)then
               an(mi)=0
               pd_flag(im,in)=1
            else
               pd_flag(im,in)=0
            endif
         enddo
      enddo
      end

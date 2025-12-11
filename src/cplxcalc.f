      subroutine cplxcalc(node,ix,complex_conc)
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
!D1 To compute the concentration of the aqueous complexes
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 Initial implementation: ?, Programmer: Hari Viswanathan
!D2 new reaction method using chemical components
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/cplxcalc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2
!**********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.4 Solute-transport equations
!D3 2.4.6 Multiple, interacting solutes
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

      use comchem
      use comdti
      use comdi
      use comrxni
      implicit none

      integer node
      integer ix
      integer ic
      integer in
      real*8 prod,tkeq(200)
      real*8 complex_conc

      prod=1
      do ic = 1, ncpnt
         if(spstoic(ix,ic).ne.0) 
     2        prod = prod*cpntsv(ic,node)**spstoic(ix,ic)
      enddo
      if(temp_model(ix).eq.'l')then
         tkeq(ix)=heq(ix,1)+heq(ix,2)*t(node)+heq(ix,3)*t(node)**2
         tkeq(ix)=10**tkeq(ix)
	elseif (temp_model(ix) .eq. 't') then
		tkeq(ix) = heq(ix, 1) + heq(ix, 2) * t(in) + heq(ix, 3) / t(in)
     2		+ heq(ix, 4) * log10(t(in)) + heq(ix, 5) / t(in) ** 2
		tkeq(ix) = 10 ** tkeq(ix)
      else
         tkeq(ix)=ckeq(ix)*exp((heq(ix,1)/gas_const*
     2        (1/(25+temp_conv)-1/(t(node)+temp_conv))))
      endif
      complex_conc = tkeq(ix)*prod
      complex_conc = max(complex_conc, 1.d-90)
      end
      

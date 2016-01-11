      subroutine co2h2o_combine(iflg,idofm)
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

C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine enforces variable equivalence between fluid
CD1  components and combines conservation equations
CD1  full jacobian (size depends on number of active variables)
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 
CD2    Rev 1.0   06/20/02 10:24:20   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.5.2 Solve nonlinear equation set at each time step
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 RJP 12/02/04 Major revision, modified to do only CO2 no hydrate
CD4
CD4
CD4
C***********************************************************************

      use comai
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comhi
      use davidi
      use comdti
      use comco2
      implicit none

      integer iflg, idofm, nmatd, index(6)
      integer i, icoupl_iter, id, idl 
      integer neqp1, nrhs1, nrhs2, nrhs3, nsizea, nsizea1 
      integer i1,i2,j
      real*8, allocatable :: dumz(:)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto5(:)   
      real*8  facr, fdum2, tollr, tolls
      
      if(idof_co2.eq.0) then
c     gaz 2-26-2003
         idofm = 6   
         return
      elseif(iprtype.eq.1) then
c     water-only problem, 2 equations originally formulated w/ 5 dof
c     reduce it to 2 dof
         if(iflg.eq.1) then
c     
c     reduce jacobian array from 5 to 2 
c     coefficients are calculated for 5 dof
c     1   2   3   4   5
c     6   7   8   9  10     
c     
c     2 dof
c     1   2
c     3   4
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
            do i = nmat(3)+1, nmat(4)
               a(i) = a(i+3*nmat(2))
            enddo

            do i = nmat(4)+1, nmat(5)
               a(i) = a(i+3*nmat(2))
            enddo
c     a(nmat(3)+1:nmat(4))  = a(nmat(6)+1:nmat(7))
c     a(nmat(4)+1:nmat(5))  = a(nmat(7)+1:nmat(8)) 
c     
            idofm = 2
c     
         else if(iflg.eq.2) then
c     
         endif                  

      elseif(iprtype.eq.2) then
c     co2-only problem, 2 equations originally formulated w/ 5 dof
c     reduce it to 2 dof
         if(iflg.eq.1) then
c     
c     reduce jacobian array from 5 to 2 
c     coefficients are calculated for 5 dof
c     11   12   13   14   15
c     16   17   18   19   20     
c     
c     2 dof
c     1   2
c     3   4
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1

            a(nmat(1)+1:nmat(2))  = a(nmat(11)+1:nmat(12))
            a(nmat(2)+1:nmat(3))  = a(nmat(12)+1:nmat(13))
            a(nmat(3)+1:nmat(4))  = a(nmat(16)+1:nmat(17))
            a(nmat(4)+1:nmat(5))  = a(nmat(17)+1:nmat(18)) 
c     
            bp(nrhs(1)+1:nrhs(2)) = bp(nrhs(3)+1:nrhs(4)) 
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(4)+1:nrhs(5)) 
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(3)+1:nrhs(4)) = 0.0d00                  
            bp(nrhs(4)+1:nrhs(5)) = 0.0d00                  
c     
c     note that the degree of freedom has been changed from 4 to 3
c     
            idofm = 2
c     
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     
            bp(nrhs(3)+1:nrhs(4)) = bp(nrhs(1)+1:nrhs(2))
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(2)+1:nrhs(3))
c     
c     
c     
            idofm = 2
c     
         endif                  

      elseif(iprtype.eq.3) then
c     
c     CO2-water problem without dissolution degrees of freedom 3
c     
         if(iflg.eq.1) then
c     
c     3 dof
c     1   2   3   4   5
c     6   7   8   9  10
c     11  12  13  14  15
c     16  17  18  19  20
c     
c     3 dof
c     1   2   3
c     4   5   6
c     7   8   9
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     enforce thermal equilibrium
c     
            do i = nmat(6)+1, nmat(7)
               a(i) = a(i)+a(i+10*nmat(2))
            enddo

            do i = nmat(7)+1, nmat(8)
               a(i) = a(i)+a(i+10*nmat(2))
            enddo

            do i = nmat(8)+1, nmat(9)
               a(i) = a(i)+a(i+10*nmat(2))
            enddo

            do i = nmat(4)+1, nmat(5)
               a(i) = a(i+2*nmat(2))
            enddo

            do i = nmat(5)+1, nmat(6)
               a(i) = a(i+2*nmat(2))
            enddo

            do i = nmat(6)+1, nmat(7)
               a(i) = a(i+2*nmat(2))
            enddo

            do i = nmat(7)+1, nmat(8)
               a(i) = a(i+4*nmat(2))
            enddo

            do i = nmat(8)+1, nmat(9)
               a(i) = a(i+4*nmat(2))
            enddo

            do i = nmat(9)+1, nmat(10)
               a(i) = a(i+4*nmat(2))
            enddo

c     a(nmat(6)+1:nmat(7))  = a(nmat(6)+1:nmat(7)) + 
c     &                         a(nmat(16)+1:nmat(17))
c     a(nmat(7)+1:nmat(8))  = a(nmat(7)+1:nmat(8)) + 
c     &                         a(nmat(17)+1:nmat(18))
c     a(nmat(8)+1:nmat(9))  = a(nmat(8)+1:nmat(9)) + 
c     &                         a(nmat(18)+1:nmat(19)) 

c     a(nmat(4)+1:nmat(5))  = a(nmat(6)+1:nmat(7))
c     a(nmat(5)+1:nmat(6))  = a(nmat(7)+1:nmat(8))
c     a(nmat(6)+1:nmat(7))  = a(nmat(8)+1:nmat(9)) 
c     a(nmat(7)+1:nmat(8))  = a(nmat(11)+1:nmat(12)) 
c     a(nmat(8)+1:nmat(9))  = a(nmat(12)+1:nmat(13)) 
c     a(nmat(9)+1:nmat(10)) = a(nmat(13)+1:nmat(14)) 
c     
c     bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3)) +
c     &                         bp(nrhs(4)+1:nrhs(5)) 
c     RJP 01/08/08 changed following
c     RJP 01/06/16 changed the following to avoid stack overflow for large problems
            do i = nrhs(2)+1, nrhs(3)
                bp(i) = bp(i)+bp(i+2*nrhs(2))
            enddo
c            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3))+
c     &           bp(nrhs(4)+1:nrhs(5))

c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(4)+1:nrhs(5)) = 0.0d00                  
c     
c     note that the degree of freedom has been changed from 4 to 3
c     
            idofm = 3
c     
         else if(iflg.eq.2) then
c     
         endif                  

      elseif(iprtype.eq.-3) then
c     isothermal co2-h2o problem without dissolution
         if(iflg.eq.1) then
c     
c     3 dof
c     1   2   3   4   5
c     6   7   8   9  10
c     11  12  13  14  15
c     16  17  18  19  20
c     
c     3 dof
c     1   2
c     3   4
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     a(nmat(2)+1:nmat(3))  = a(nmat(3)+1:nmat(4))
c     a(nmat(3)+1:nmat(4))  = a(nmat(11)+1:nmat(12)) 
c     a(nmat(4)+1:nmat(5))  = a(nmat(13)+1:nmat(14)) 
            do i = nmat(2)+1, nmat(3)
               a(i)=a(i+nmat(2))
            enddo
            do i = nmat(3)+1, nmat(4)
               a(i)=a(i+8*nmat(2))
            enddo
            do i = nmat(4)+1, nmat(5)
               a(i)=a(i+9*nmat(2))
            enddo
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(3)+1:nrhs(4))

c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(3)+1:nrhs(4)) = 0.0d00                  
            bp(nrhs(4)+1:nrhs(5)) = 0.0d00                  
c     
c     note that the degree of freedom has been changed from 3 to 2
c     
            idofm = 2
c     
         else if(iflg.eq.2) then
c     
            idofm = 3
         endif                  

      elseif(iprtype.eq.4) then
c     
c     water-CO2 problem, no air, with co2 dissolution in water                             
c     C of M water (row 1)
c     C of E water (row 2)
c     C of M co2 (row 3)
c     C of E co2 (row 4)
c     variables:P,T,Fracw/yco2. variables yco2 and fracw switch depending
c     upon whether max co2 solubility is achieved
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     3 independent variables are assumed here
c     condense 4 equations to 3 with thermal equilibrium
c     3 dof (4 equations)
c     1   2   3   
c     6   7   8   
c     11  12  13  
c     16  17  18
c     condensed to   
c     
c     3 dof (3 equations)
c     1   2   3
c     4   5   6
c     7   8   9
c     
c     3 dof and 3 equations as follows
c     dof          P		   T	  fw/yc									
c     CofMw         1         2         3   
c     Energy    (6+16)    (7+17)    (8+18)  
c     CofMc        11        12        13   
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     add thermal equations
		  do i = nmat(6)+1, nmat(7)
			a(i) = a(i)+a(i+10*nmat(2))
		  enddo
		  do i = nmat(7)+1, nmat(8)
			a(i) = a(i)+a(i+10*nmat(2))
		  enddo
		  do i = nmat(8)+1, nmat(9)
			a(i) = a(i)+a(i+10*nmat(2))
		  enddo

c            a(nmat(6)+1:nmat(7))   = a(nmat(6)+1:nmat(7)) 
c     &           + a(nmat(16)+1:nmat(17))
c            a(nmat(7)+1:nmat(8))   = a(nmat(7)+1:nmat(8)) 
c     &           + a(nmat(17)+1:nmat(18))
c            a(nmat(8)+1:nmat(9))   = a(nmat(8)+1:nmat(9)) 
c     &           + a(nmat(18)+1:nmat(19))
c     combine for thermal equilibrium
c     

c            bp(nrhs(2)+1:nrhs(3))  = bp(nrhs(2)+1:nrhs(3))   
c     &           + bp(nrhs(4)+1:nrhs(5))
c     RJP 01/06/16 Changed the following to avoid stack overflow crashes for large problems

            do i = nrhs(2)+1, nrhs(3)
                bp(i) = bp(i)+bp(i+2*nrhs(2))
            enddo

		  do i = nmat(4)+1, nmat(5)
			a(i) = a(i+2*nmat(2))
	      enddo
		  do i = nmat(5)+1, nmat(6)
			a(i) = a(i+2*nmat(2))
	      enddo
		  do i = nmat(6)+1, nmat(7)
			a(i) = a(i+2*nmat(2))
	      enddo
		  do i = nmat(7)+1, nmat(8)
			a(i) = a(i+4*nmat(2))
	      enddo
		  do i = nmat(8)+1, nmat(9)
			a(i) = a(i+4*nmat(2))
	      enddo
		  do i = nmat(9)+1, nmat(10)
			a(i) = a(i+4*nmat(2))
	      enddo
c            a(nmat(4)+1:nmat(5))   = a(nmat(6)+1:nmat(7))
c            a(nmat(5)+1:nmat(6))   = a(nmat(7)+1:nmat(8))
c            a(nmat(6)+1:nmat(7))   = a(nmat(8)+1:nmat(9)) 
c            a(nmat(7)+1:nmat(8))   = a(nmat(11)+1:nmat(12)) 
c            a(nmat(8)+1:nmat(9))   = a(nmat(12)+1:nmat(13)) 
c            a(nmat(9)+1:nmat(10))  = a(nmat(13)+1:nmat(14)) 
c     
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
c     bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3))+ 
c     &						 bp(nrhs(4)+1:nrhs(5)) 
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(4)+1:nrhs(5)) = 0.0d00                  
c     
            idofm = 3
c     
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays don't need to re-arrange bp 
c     (taken care of in N-R update code) change back the degree of 
c     freedom to 4
c     
            idofm = 5
c     
         endif                  


      elseif(idof_co2.eq.5) then
c     
c     full problem                              
c     C of M water (row 1)
c     C of E water (row 2)
c     C of M co2 (row 3)
c     C of E co2 (row 4)
c     Thermodynamic Eqbm. of co2 (row 5)
c     C of M air (row 6)
c     C of E air (row 7)
c     variables:P,T,Fracw,yco2, yair
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     5 independent variables are assumed here
c     condense 7 equations to 5 with thermal equilibrium
c     5 dof (7 equations)
c     1   2   3   4   5
c     6   7   8   9  10  
c     11  12  13  14  15
c     16  17  18  19  20
c     21  22  23  24  25
c     26  27  28  29  30
c     31  32  33  34  35

c     condensed to   
c     
c     5 dof (5 equations)
c     1   2   3   4   5
c     6   7   8   9  10
c     11  12  13  14  15
c     16  17  18  19  20
c     21  22  23  24  25
c     
c     The equations and variables are as follows
c     P         T        fw        yc         ya
c     CMW           1         2         3        4           5
c     E     (6+16+31) (7+17+32) (8+18+33) (9+19+34) (10+20+35)
c     CMC          11        12        13       14          15
c     TC           21        22        23       24          25
c     CMA          26        27        28       29          30
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
            a(nmat(6)+1:nmat(7))   = a(nmat(6)+1:nmat(7)) 
     &           + a(nmat(16)+1:nmat(17))
     &           + a(nmat(31)+1:nmat(32))
            a(nmat(7)+1:nmat(8))   = a(nmat(7)+1:nmat(8)) 
     &           + a(nmat(17)+1:nmat(18))
     &           + a(nmat(32)+1:nmat(33))
            a(nmat(8)+1:nmat(9))   = a(nmat(8)+1:nmat(9)) 
     &           + a(nmat(18)+1:nmat(19))
     &           + a(nmat(33)+1:nmat(34))
            a(nmat(9)+1:nmat(10))  = a(nmat(9)+1:nmat(10)) 
     &           + a(nmat(19)+1:nmat(20))
     &           + a(nmat(34)+1:nmat(35))
            a(nmat(10)+1:nmat(11)) = a(nmat(10)+1:nmat(11)) 
     &           + a(nmat(20)+1:nmat(21))
     &           + a(nmat(35)+1:nmat(36))
            a(nmat(11)+1:nmat(12))  = a(nmat(11)+1:nmat(12)) 
            a(nmat(12)+1:nmat(13))  = a(nmat(12)+1:nmat(13)) 
            a(nmat(13)+1:nmat(14))  = a(nmat(13)+1:nmat(14)) 
            a(nmat(14)+1:nmat(15))  = a(nmat(14)+1:nmat(15)) 
            a(nmat(15)+1:nmat(16))  = a(nmat(15)+1:nmat(16)) 
            a(nmat(16)+1:nmat(17))  = a(nmat(21)+1:nmat(22)) 
            a(nmat(17)+1:nmat(18))  = a(nmat(22)+1:nmat(23)) 
            a(nmat(18)+1:nmat(19))  = a(nmat(23)+1:nmat(24)) 
            a(nmat(19)+1:nmat(20))  = a(nmat(24)+1:nmat(25)) 
            a(nmat(20)+1:nmat(21))  = a(nmat(25)+1:nmat(26)) 
            a(nmat(21)+1:nmat(22))  = a(nmat(26)+1:nmat(27)) 
            a(nmat(22)+1:nmat(23))  = a(nmat(27)+1:nmat(28)) 
            a(nmat(23)+1:nmat(24))  = a(nmat(28)+1:nmat(29)) 
            a(nmat(24)+1:nmat(25))  = a(nmat(29)+1:nmat(30)) 
            a(nmat(25)+1:nmat(26))  = a(nmat(30)+1:nmat(31)) 
c     
c     combine for thermal equilibrium
c     
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3))   
     &           + bp(nrhs(4)+1:nrhs(5))
     &           + bp(nrhs(7)+1:nrhs(8))
c     just move up other equations
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(5)+1:nrhs(6))
            bp(nrhs(5)+1:nrhs(6)) = bp(nrhs(6)+1:nrhs(7))   
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
            bp(nrhs(6)+1:nrhs(8)) = 0.d0
c     
            idofm = 5
c     
c     call routine to further reduce the degrees of freedom from 5 to 4
c     
            call co2_array(1,idofm)
c     
         else if(iflg.eq.2) then
c     
            call co2_array(2,idofm)
c     
c     replace variables in correct arrays
c     
c     don't need to re-arrange bp (taken care of in N-R update code)
c     
c     change back the degree of freedom from 4 to 6
c     
            idofm = 5
c     
         endif                  


      elseif(idof_co2.eq.-6) then
c     
c     RJP 12/03/04 This is for CO2-water problem without hydrate
c     Salt precipitation is ignored.  The equations are
c     Conservation of Mass of Water (Row 1)
c     Conservation of Energy of Water (Row 2)
c     Conservation of Mass of CO2 (Row 3)
c     Conservation of Energy of CO2 (Row 4)
c     Assume thermal equilibrium
c     Variables are P, T, Fracw, FracCO2-gas
c     Degrees of freedom are reduced from 4 to 3.
c     If CO2 is single phase: the variables are P, T and Fracw
c     If CO2 is double phase: the variables are P, Fracw  and FracCO2g
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     4 dof
c     1   2   3   4
c     5   6   7   8
c     9  10  11  12
c     13  14  15  16
c     
c     Thermal equilibrium
c     3 equations
c     1      2      3      4
c     (5+13) (6+14) (7+15) (8+16)
c     9      10     11     12
c     3 dof (final) based on phase state
c     if ieosd is 2 then variables are P, FracCO2g, Fracw
c     if ieosd is 1 then variables are P, T, Fracw
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     First do thermal equilibrium
            a(nmat(5)+1:nmat(6))  = a(nmat(5)+1:nmat(6)) 
     &           + a(nmat(13)+1:nmat(14))
            a(nmat(6)+1:nmat(7))  = a(nmat(6)+1:nmat(7)) 
     &           + a(nmat(14)+1:nmat(15))
            a(nmat(7)+1:nmat(8))  = a(nmat(7)+1:nmat(8)) 
     &           + a(nmat(15)+1:nmat(16))
            a(nmat(8)+1:nmat(9))  = a(nmat(8)+1:nmat(9)) 
     &           + a(nmat(16)+1:nmat(17))
c     
c     combine for thermal equilibrium
c     
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3))   
     &           + bp(nrhs(4)+1:nrhs(5))
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(4)+1:nrhs(5)) = 0.d0
c     note that the degree of freedom has been changed from 4 to 3
c     
            idofm = 3
c     Modify the variables based on Phase change. If two phase
c     CO2, substitute for temperature with fraction of CO2 gas
c     Pressure derivative is already modified in ther_co2_h2o
c     
c     RJP 12/23/04 change below later
            do i = 1, neq
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do j = i1, i2
                  if (ices(i).eq.2) then
                     a(j+nmat(2)-(neq+1)) = a(j+nmat(4)-(neq+1))
                     a(j+nmat(6)-(neq+1)) = a(j+nmat(8)-(neq+1))
                     a(j+nmat(10)-(neq+1)) = a(j+nmat(12)-(neq+1))
                  endif
               enddo
            enddo
c     RJP 12/23/04 change above later
c     RJP 12/20/04 Modify array 'a' to reflect the change of degrees
c     of freedom. Current 'a'
c     1  2  3  4
c     5  6  7  8
c     9 10 11 12
c     should be modified as follows
c     4 = 5, 5 = 6, 6 = 7, 7 = 9, 8 = 10, 9 = 11
            a(nmat(4)+1:nmat(5)) = a(nmat(5)+1:nmat(6))
            a(nmat(5)+1:nmat(6)) = a(nmat(6)+1:nmat(7))
            a(nmat(6)+1:nmat(7)) = a(nmat(7)+1:nmat(8))
            a(nmat(7)+1:nmat(8)) = a(nmat(9)+1:nmat(10))
            a(nmat(8)+1:nmat(9)) = a(nmat(10)+1:nmat(11))
            a(nmat(9)+1:nmat(10)) = a(nmat(11)+1:nmat(12)) 
            
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     don't need to re-arrange bp (taken care of in N-R update code)
c     change back the degree of freedom from 3 to 4
c     
            idofm = 4
c     
         endif                  
      elseif(idof_co2.eq.-7) then
c     
c     RJP 04/16/05 This is for CO2 only problem. The equations for
c     CO2/water problem are further reduced to CO2 only equations. 
c     The original equations are 
c     Conservation of Mass of Water (Row 1)
c     Conservation of Energy of Water (Row 2)
c     Conservation of Mass of CO2 (Row 3)
c     Conservation of Energy of CO2 (Row 4)
c     Assume thermal equilibrium
c     Variables are P, T, Fracw, FracCO2-gas
c     Degrees of freedom are reduced from 4 to 2.
c     If CO2 is single phase: the variables are P, T
c     If CO2 is double phase: the variables are P, FracCO2g
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     4 dof
c     1   2   3   4
c     5   6   7   8
c     9  10  11  12
c     13  14  15  16
c     
c     3 dof (final) based on phase state
c     if ieosd is 2 then variables are P, FracCO2g
c     if ieosd is 1 then variables are P, T
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1

c     Move CO2 equations to replace water equations.

            bp(nrhs(1)+1:nrhs(2)) = bp(nrhs(3)+1:nrhs(4))
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(4)+1:nrhs(5))
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(3)+1:nrhs(4)) = 0.d0     
            bp(nrhs(4)+1:nrhs(5)) = 0.d0

c     note that the degree of freedom has been changed from 4 to 2
c     
            idofm = 2

c     First move the variables for CO2 in place of water
            a(nmat(1)+1:nmat(2))  = a(nmat(9)+1:nmat(10)) 
            a(nmat(2)+1:nmat(3))  = a(nmat(10)+1:nmat(11)) 
            a(nmat(3)+1:nmat(4))  = a(nmat(11)+1:nmat(12)) 
            a(nmat(4)+1:nmat(5))  = a(nmat(12)+1:nmat(13)) 
            a(nmat(5)+1:nmat(6))  = a(nmat(13)+1:nmat(14)) 
            a(nmat(6)+1:nmat(7))  = a(nmat(14)+1:nmat(15)) 
            a(nmat(7)+1:nmat(8))  = a(nmat(15)+1:nmat(16)) 
            a(nmat(8)+1:nmat(9))  = a(nmat(16)+1:nmat(17)) 
c     get rid of the water variable dependent terms
            a(nmat(3)+1:nmat(4))  = a(nmat(4)+1:nmat(5)) 
            a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6))
            a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7)) 
            a(nmat(6)+1:nmat(7))  = a(nmat(8)+1:nmat(9)) 

c     Modify the variables based on Phase change. If two phase
c     CO2, substitute for temperature with fraction of CO2 gas
c     Pressure derivative is already modified in ther_co2_h2o
c     
c     RJP 12/23/04 change below later
            do i = 1, neq
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do j = i1, i2
                  if (ices(i).ne.2) then
                     a(j+nmat(3)-(neq+1)) = a(j+nmat(4)-(neq+1))
                     a(j+nmat(4)-(neq+1)) = a(j+nmat(5)-(neq+1))
                  else
                     a(j+nmat(2)-(neq+1)) = a(j+nmat(3)-(neq+1))
                     a(j+nmat(3)-(neq+1)) = a(j+nmat(4)-(neq+1))
                     a(j+nmat(4)-(neq+1)) = a(j+nmat(6)-(neq+1))
                  endif
               enddo
            enddo
            
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     don't need to re-arrange bp (taken care of in N-R update code)
c     change back the degree of freedom from 2 to 3
c     
            idofm = 3
c     
         endif                  


      elseif(idof_co2.eq.1) then
c     
c     degrees of freedom reduced from 4 to 3
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     4 independent variables are assumed here(P1,P2,T1,S2)
c     condense dof from 4 to 3  (Pressure equilibrium)
c     combine energy equation
c     operations below assume arrays 3,4,7,8,9,10,13,14 = 0
c     
c     4 dof
c     1   2   3   4
c     5   6   7   8
c     9  10  11  12
c     13  14  15  16
c     
c     3 dof
c     1   2   3
c     4   5   6
c     7   8   9
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     enforce thermal equilibrium
c     
            a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6)) + 
     &           a(nmat(15)+1:nmat(16))
            a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7)) + 
     &           a(nmat(14)+1:nmat(15))
            a(nmat(6)+1:nmat(7))  = a(nmat(16)+1:nmat(17)) 
            a(nmat(7)+1:nmat(8))  = a(nmat(11)+1:nmat(12)) 
            a(nmat(8)+1:nmat(9))  = a(nmat(10)+1:nmat(11)) 
            a(nmat(9)+1:nmat(10)) = a(nmat(12)+1:nmat(13)) 
c     
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3)) +
     &           bp(nrhs(4)+1:nrhs(5)) 
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
            bp(nrhs(4)+1:nrhs(5)) = 0.0d00                  
c     
c     note that the degree of freedom has been changed from 4 to 3
c     
            idofm = 3
c     
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     
c     
c     put S2 in ccorrect array           
c     enforce P2=P1                      
c     
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(3)+1:nrhs(4))
            bp(nrhs(3)+1:nrhs(4)) = bp(nrhs(1)+1:nrhs(2))
c     
c     change back the degree of freedom from 3 to 4
c     
            idofm = 4
c     
         endif                  

      endif                  
      return
      end

       subroutine methh2o_combine(iflg,idofm)
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
!D1  PURPOSE
!D1
!D1  This subroutine enforces variable equivalence between fluid
!D1  components and combines conservation equations
!D1  full jacobian (size depends on number of active variables)
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/methh2o_combine.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  !D3  2.5.2 Solve nonlinear equation set at each time step
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

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
      use commeth
      implicit none

      integer iflg, idofm, nmatd, index(6)
      integer i, icoupl_iter, id, idl, i1, i2, ii, kb
      integer neqp1, nrhs1, nrhs2, nrhs3, nsizea, nsizea1 
      real*8, allocatable :: dumz(:)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto5(:)   
      real*8  facr, fdum2, tollr, tolls, dtp, dshsw
      
      if(idof_meth.eq.0) then
c     gaz 2-26-2003
         idofm = 6   
         return

      elseif(idof_meth.eq.2) then
c     
c     degrees of freedom reduced from 4 to 2
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     reduce jacobian array from 4 to 2 
c     assume thermal and pressure equilibrium
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
c     2 dof
c     1   2
c     3   4
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     enforce thermal equilibrium 4 to 3 dof
c     note f90 array equivalencing
c     also this assume arrays 3,4,7,8 and 9,10,13,14
c     are empty (only if distinct fluids)
c     
            a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6)) 
            a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7)) + 
     &           a(nmat(16)+1:nmat(17))
            a(nmat(6)+1:nmat(7))  = a(nmat(15)+1:nmat(16)) 
            a(nmat(8)+1:nmat(9))  = a(nmat(12)+1:nmat(13)) 
            a(nmat(9)+1:nmat(10)) = a(nmat(11)+1:nmat(12)) 
c     
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3)) +
     &           bp(nrhs(4)+1:nrhs(5)) 
c     
c     enforce thermal equilibrium 3 to 2 dof
c     
            a(nmat(1)+1:nmat(2))  = a(nmat(1)+1:nmat(2)) + 
     &           a(nmat(7)+1:nmat(8))
            a(nmat(2)+1:nmat(3))  = a(nmat(2)+1:nmat(3)) + 
     &           a(nmat(8)+1:nmat(9))
            a(nmat(3)+1:nmat(4))  = a(nmat(3)+1:nmat(4)) + 
     &           a(nmat(9)+1:nmat(10))

            bp(nrhs(1)+1:nrhs(2)) = bp(nrhs(1)+1:nrhs(2)) +
     &           bp(nrhs(3)+1:nrhs(4)) 

c     now assume form is 2 DOF
            a(nmat(1)+1:nmat(2))  = a(nmat(1)+1:nmat(2)) +
     &           a(nmat(3)+1:nmat(4))
            a(nmat(3)+1:nmat(4))  = a(nmat(4)+1:nmat(5)) + 
     &           a(nmat(6)+1:nmat(7))
            a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6)) 
c     
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
c     
c     enforce thermal equilibrium
c     
            bp(nrhs(3)+1:nrhs(4)) = bp(nrhs(1)+1:nrhs(2))
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(2)+1:nrhs(3))
c     
c     change back the degree of freedom from 2 to 4
c     
            idofm = 4
c     
         endif                  

      elseif(idof_meth.eq.3) then
c     
c     degrees of freedom reduced from 4 to 3
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     4 independent variables are assumed here
c     condense dof from 4 to 3  (thermal equilibrium)
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
            a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6)) 
            a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7)) + 
     &           a(nmat(16)+1:nmat(17))
            a(nmat(6)+1:nmat(7))  = a(nmat(15)+1:nmat(16)) 
            a(nmat(8)+1:nmat(9))  = a(nmat(12)+1:nmat(13)) 
            a(nmat(9)+1:nmat(10)) = a(nmat(11)+1:nmat(12)) 
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
c     enforce thermal equilibrium
c     
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(2)+1:nrhs(3))
c     
c     change back the degree of freedom from 3 to 4
c     
            idofm = 4
c     
         endif                  

      elseif(idof_meth.eq.4) then
c     
c     degrees of freedom reduced from 5 to 4
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     5 independent variables are assumed here
c     condense dof from 5 to 4  (with thermal equilibrium)
c     5 dof (but stored like 4 dof)
c     1   2   3   4   5
c     6   7   8   9  10
c     11  12  13  14  15
c     16  17  18  19  20
c     21  22  23  24  25

c     derivatives wrt frac_meth in 17,18,19,20

c     other terms for derivatives wrt frac_meth
c     
c     4 dof
c     1   2   3   4
c     5   6   7   8
c     9  10  11  12
c     13  14  15  16
c     
c     4 dof (final)
c     1   2   0  17
c     5   6   0  18
c     0  12  11  19
c     0  16  15  20
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     enforce temperature (methane) = temperature(water)
c     
            a(nmat(10)+1:nmat(11))  = a(nmat(12)+1:nmat(13)) 
            a(nmat(14)+1:nmat(15))  = a(nmat(16)+1:nmat(17)) 
c     D(WM)/D(FM)
            a(nmat(4)+1:nmat(5))  = a(nmat(17)+1:nmat(18)) 
c     D(WE)/D(FM)
            a(nmat(8)+1:nmat(9))  = a(nmat(18)+1:nmat(19)) 
c     D(MM)/D(FM)
            a(nmat(12)+1:nmat(13))  = a(nmat(19)+1:nmat(20)) 
c     D(ME)/D(FM)
            a(nmat(16)+1:nmat(17))  = a(nmat(20)+1:nmat(21)) 
c     
c     
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
c     note that the degree of freedom has been changed from 4 to 3
c     
            idofm = 4
c     
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     
c     
c     enforce thermal equilibrium
c     
c     meth fraction update
            bp(nrhs(5)+1:nrhs(6)) = bp(nrhs(4)+1:nrhs(5))   
c     tmeth = t (water)      
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(2)+1:nrhs(3))   
c     
c     
c     change back the degree of freedom from 3 to 4
c     
            idofm = 4
c     
         endif                  

      elseif(idof_meth.eq.5) then
c     
c     degrees of freedom reduced from 6 to 4
c     then reduced from 4 to 3 in call to hydrate_array
c     full problem                              
c     C of M water (row 1)
c     C of E water (row 2)
c     C of M methane (row 3)
c     C of E methane (row 4)
c     C of M hydrate (row 5)
c     C of E hydrate (row 6)
c     
c     variables:Pw,Tw,Pmeth,Tmeth,Fracw,Fracmeth
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     6 independent variables are assumed here
c     condense dof from 6 to 4  (with thermal,pressure equilibrium)
c     6 dof (but stored like 4 dof)
c     1   2   3   4   5   6
c     7   8   9  10  11  12
c     13  14  15  16  17  18
c     19  20  21  22  23  24
c     25  26  27  28  29  30
c     31  32  33  34  35  36

c     condensed to   
c     
c     4 dof
c     1   2   3   4
c     5   6   7   8
c     9  10  11  12
c     13  14  15  16
c     
c     4 dof (final) based on original 6 dof numbers
c     (1+3)             (2+4)              5          6
c     (7+19+31+9+21+33) (8+20+32+10+22+34) (11+23+35) (12+24+36)
c     (13+15)           (14+16)            17         18            
c     (25+27)           (26+28)            29         30
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
            a(nmat(1)+1:nmat(2))  = a(nmat(1)+1:nmat(2)) 
     &           + a(nmat(3)+1:nmat(4))
            a(nmat(2)+1:nmat(3))  = a(nmat(2)+1:nmat(3)) 
     &           + a(nmat(4)+1:nmat(5))
            a(nmat(3)+1:nmat(4))  = a(nmat(5)+1:nmat(6)) 
            a(nmat(4)+1:nmat(5))  = a(nmat(6)+1:nmat(7)) 
            a(nmat(5)+1:nmat(6))  = a(nmat(7)+1:nmat(8)) 
     &           + a(nmat(19)+1:nmat(20))
     &           + a(nmat(31)+1:nmat(32))
     &           + a(nmat(9)+1:nmat(10))
     &           + a(nmat(21)+1:nmat(22))
     &           + a(nmat(33)+1:nmat(34))
            a(nmat(6)+1:nmat(7))  = a(nmat(8)+1:nmat(9)) 
     &           + a(nmat(20)+1:nmat(21))
     &           + a(nmat(32)+1:nmat(33))
     &           + a(nmat(10)+1:nmat(11))
     &           + a(nmat(22)+1:nmat(23))
     &           + a(nmat(34)+1:nmat(35))
            a(nmat(7)+1:nmat(8))  = a(nmat(11)+1:nmat(12)) 
     &           + a(nmat(23)+1:nmat(24))
     &           + a(nmat(35)+1:nmat(36))
            a(nmat(8)+1:nmat(9))  = a(nmat(12)+1:nmat(13)) 
     &           + a(nmat(24)+1:nmat(25))
     &           + a(nmat(36)+1:nmat(37))
            a(nmat(9)+1:nmat(10))  = a(nmat(13)+1:nmat(14)) 
     &           + a(nmat(15)+1:nmat(16))
            a(nmat(10)+1:nmat(11))  = a(nmat(14)+1:nmat(15)) 
     &           + a(nmat(16)+1:nmat(17))
            a(nmat(11)+1:nmat(12))  = a(nmat(17)+1:nmat(18)) 
            a(nmat(12)+1:nmat(13))  = a(nmat(18)+1:nmat(19)) 
            a(nmat(13)+1:nmat(14))  = a(nmat(25)+1:nmat(26)) 
     &           + a(nmat(27)+1:nmat(28))
            a(nmat(14)+1:nmat(15))  = a(nmat(26)+1:nmat(27)) 
     &           + a(nmat(28)+1:nmat(29))
            a(nmat(15)+1:nmat(16))  = a(nmat(29)+1:nmat(30)) 
            a(nmat(16)+1:nmat(17))  = a(nmat(30)+1:nmat(31)) 
c     
c     combine for thermal equilibrium
c     
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3))   
     &           + bp(nrhs(4)+1:nrhs(5))
     &           + bp(nrhs(6)+1:nrhs(7))
c     just move up to 4th equation
            bp(nrhs(4)+1:nrhs(5)) = bp(nrhs(5)+1:nrhs(6))   
c     
c     zero out residual that has been combined so it does not interfere
c     with convergence
c     
c     note that the degree of freedom has been changed from 6 to 4
c     
            idofm = 4
c     
c     call routine to further reduce the degrees of freedom from 4 to 3
c     
            call hydrate_array(1,idofm)
c     
         else if(iflg.eq.2) then
c     
            call hydrate_array(2,idofm)
c     
c     replace variables in correct arrays
c     
c     don't need to re-arrange bp (taken care of in N-R update code)
c     
c     change back the degree of freedom from 4 to 6
c     
            idofm = 6
c     
         endif                  

      elseif(idof_meth.eq.1) then
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
      elseif(idof_meth.eq.7) then
c     
c     degrees of freedom reduced from 6 to 4 (equilibrium conditions)
c     then reduced from 4 to 3 in call to hydrate_array
c     full problem                              
c     C of M water (row 1)
c     C of E water (row 2)
c     C of M methane (row 3)
c     C of E methane (row 4)
c     C of M hydrate (row 5)
c     C of E hydrate (row 6)
c     
c     variables:Pw,Tw,Pmeth,Tmeth,Fracw,Frach
c     
         if(iflg.eq.1) then
c     
c     adjust storage
c     
c     6 independent variables are assumed here
c     condense dof from 6 to 4  (with thermal,pressure equilibrium)
c     6 dof (but stored like 4 dof)
c     1   2   3   4   5   6
c     7   8   9  10  11  12
c     13  14  15  16  17  18
c     19  20  21  22  23  24
c     25  26  27  28  29  30
c     31  32  33  34  35  36
c     in equilibrium calcs COM and COE for hydrate is not formed
c     condensed to   
c     
c     3 dof
c     1   2   3
c     4   5   6  
c     7   8   9    
c     
c     4  by 3 dof (final) based on original 6 dof numbers
c     (1+3)               (2+4)              5          6
c     (7+19+9+21)  (8+20+10+22)         (11+23)   (12+24)
c     (13+15)           (14+16)            17         18            
c     in equilibrium calcs, either t dependence or frach dependence
c     but not both
c     reduce to 

c     
            idofm = 3
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
            a(nmat(1)+1:nmat(2))  = a(nmat(1)+1:nmat(2)) 
     &           + a(nmat(3)+1:nmat(4))
            a(nmat(2)+1:nmat(3))  = a(nmat(2)+1:nmat(3)) 
     &           + a(nmat(4)+1:nmat(5))
            a(nmat(3)+1:nmat(4))  = a(nmat(5)+1:nmat(6)) 
            a(nmat(4)+1:nmat(5))  = a(nmat(6)+1:nmat(7)) 
            a(nmat(5)+1:nmat(6))  = a(nmat(7)+1:nmat(8)) 
     &           + a(nmat(19)+1:nmat(20))
     &           + a(nmat(9)+1:nmat(10))
     &           + a(nmat(21)+1:nmat(22))
            a(nmat(6)+1:nmat(7))  = a(nmat(8)+1:nmat(9)) 
     &           + a(nmat(20)+1:nmat(21))
     &           + a(nmat(10)+1:nmat(11))
     &           + a(nmat(22)+1:nmat(23))
            a(nmat(7)+1:nmat(8))  = a(nmat(11)+1:nmat(12)) 
     &           + a(nmat(23)+1:nmat(24))
            a(nmat(8)+1:nmat(9))  = a(nmat(12)+1:nmat(13)) 
     &           + a(nmat(24)+1:nmat(25))
            a(nmat(9)+1:nmat(10))  = a(nmat(13)+1:nmat(14)) 
     &           + a(nmat(15)+1:nmat(16))
            a(nmat(10)+1:nmat(11))  = a(nmat(14)+1:nmat(15)) 
     &           + a(nmat(16)+1:nmat(17))
            a(nmat(11)+1:nmat(12))  = a(nmat(17)+1:nmat(18)) 
            a(nmat(12)+1:nmat(13))  = a(nmat(18)+1:nmat(19)) 

c     
c     combine for thermal equilibrium
c     added + bp(nrhs(6)+1:nrhs(7)) here gaz 1-17-05
            bp(nrhs(2)+1:nrhs(3)) = bp(nrhs(2)+1:nrhs(3))   
     &           + bp(nrhs(4)+1:nrhs(5))
c     zero out residuals that have been combined 
            bp(nrhs(4)+1:nrhs(5)) =0.0
            bp(nrhs(5)+1:nrhs(6)) =0.0
            bp(nrhs(6)+1:nrhs(7)) =0.0
c     
c     condense either temperature or hydrate fraction
c     
            do i = 1, n
c     changed  5-18-05 gaz  
c     pressure is at the hydrate forming pressure 
c     independent variables are p,fracw,frac
c     dtp has been calculated and stored in a(nmat(25))
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do ii = i1, i2
                  kb = nelm(ii)
                  if(ihyd(kb).eq.0) then
                     dtp = a(nmat(25)+kb)
                     a(nmat(1)+ii-neqp1) = a(nmat(1)+ii-neqp1) +
     &                    a(nmat(2)+ii-neqp1)*dtp
                     a(nmat(5)+ii-neqp1) = a(nmat(5)+ii-neqp1) +
     &                    a(nmat(6)+ii-neqp1)*dtp
                     a(nmat(9)+ii-neqp1) = a(nmat(9)+ii-neqp1) +
     &                    a(nmat(10)+ii-neqp1)*dtp
                     a(nmat(4)+ii-neqp1) = a(nmat(4)+ii-neqp1) -
     &                    a(nmat(3)+ii-neqp1)*fracl
                     a(nmat(8)+ii-neqp1) = a(nmat(8)+ii-neqp1) -
     &                    a(nmat(7)+ii-neqp1)*fracl
                     a(nmat(12)+ii-neqp1) = a(nmat(12)+ii-neqp1) -
     &                    a(nmat(11)+ii-neqp1)*fracl
                     
                  else if(ihyd(kb).eq.-1) then
c     should cover both gas and water fractions (swf=sw-fracl*sh)
                     dshsw = a(nmat(25)+kb)
                     a(nmat(3)+ii-neqp1) = a(nmat(3)+ii-neqp1) -
     &                    a(nmat(3)+ii-neqp1)*dshsw*fracl +
     &                    a(nmat(4)+ii-neqp1)*dshsw
                     a(nmat(7)+ii-neqp1) = a(nmat(7)+ii-neqp1) -
     &                    a(nmat(7)+ii-neqp1)*dshsw*fracl +
     &                    a(nmat(8)+ii-neqp1)*dshsw
                     a(nmat(11)+ii-neqp1) = a(nmat(11)+ii-neqp1) -
     &                    a(nmat(11)+ii-neqp1)*dshsw*fracl +
     &                    a(nmat(12)+ii-neqp1)*dshsw
                  else if(ihyd(kb).eq.-2) then
c     remember come out of eos as dCM/dswf and dCM/dsh
c     first get hydrate derivatives (swf=sw-fracl*sh)
                     a(nmat(4)+ii-neqp1) = a(nmat(4)+ii-neqp1) -
     &                    a(nmat(3)+ii-neqp1)*fracl
                     a(nmat(8)+ii-neqp1) = a(nmat(8)+ii-neqp1) -
     &                    a(nmat(7)+ii-neqp1)*fracl
                     a(nmat(12)+ii-neqp1) = a(nmat(12)+ii-neqp1) -
     &                    a(nmat(11)+ii-neqp1)*fracl
c     should cover both gas and water fractions (use sg = 1-sw)
                     dshsw = a(nmat(25)+kb)
                     a(nmat(3)+ii-neqp1) = -a(nmat(3)+ii-neqp1) +
     &                    a(nmat(4)+ii-neqp1)*dshsw*fracl
                     a(nmat(7)+ii-neqp1) = -a(nmat(7)+ii-neqp1) +
     &                    a(nmat(8)+ii-neqp1)*dshsw*fracl
                     a(nmat(11)+ii-neqp1) = -a(nmat(11)+ii-neqp1) +
     &                    a(nmat(12)+ii-neqp1)*dshsw*fracl
                  endif
               enddo
            enddo
            do i = 1, n

               i1 = nelm(i)+1
               i2 = nelm(i+1) 
               do ii = i1, i2
                  kb = nelm(ii)
                  if(ihyd(kb).eq.0) then
c     put hydrate in second slot
                     a(nmat(2)+ii-neqp1) = a(nmat(4)+ii-neqp1)
                     a(nmat(6)+ii-neqp1) = a(nmat(8)+ii-neqp1)
                     a(nmat(10)+ii-neqp1) = a(nmat(12)+ii-neqp1)
                  endif
               enddo
            enddo
c     
c     now do last condesation (to complete 3 by 3)
c     
            a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6)) 
            a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7))   
            a(nmat(6)+1:nmat(7))  = a(nmat(7)+1:nmat(8)) 
            a(nmat(7)+1:nmat(8))  = a(nmat(9)+1:nmat(10)) 
            a(nmat(8)+1:nmat(9))  = a(nmat(10)+1:nmat(11)) 
            a(nmat(9)+1:nmat(10))  = a(nmat(11)+1:nmat(12))  
c     
         else if(iflg.eq.2) then
c     

         endif                  

      endif                  
      return
      end

      subroutine hydrate_array(iflg,idofm)
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
!D1 To simplify arrays in methane hydrate calculations
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/hydrate_array.f_a  $
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.2 Heat- and mass-transfer equations
!D3 2.3.3 Noncondensible gas flow equations
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
c
c called from methh20_combine
c
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

      integer iflg, idofm, nmatd, i, id, ii
      integer neqp1
      integer nh1(9), nh2(9), nh3(9), nh4(3)   

      data nh1  /1, 2, 3, 5, 6, 7, 9, 10, 11/
      data nh2  /4, 4, 4, 8, 8, 8, 12, 12, 12/
      data nh3  /13, 14, 15, 13, 14, 15, 13, 14, 15/
      data nh4  /4, 8, 12/
     
c
c      return 
c
      neqp1=neq+1
      if(iflg.eq.1) then  

c      nmatd is size of one subarray
       nmatd=nelm(neqp1)-neqp1
c
c combine parts of a array and bp array
c
      if(.not.allocated(a44i)) allocate(a44i(neq))
c a(nmat(16)+1:nmat(17)) is a44
       do i = 1,neq
        id = nelmdg(i)-neqp1
        a44i(i) = 1.0/a(id+nmat(16))
       enddo
c
       do ii = 1,9
        do i = 1,neq
         id = nelmdg(i)-neqp1
         a(id+nmat(nh1(ii))) = a(id+nmat(nh1(ii))) - 
     &    a(id+nmat(nh2(ii)))*a44i(i)*a(id+nmat(nh3(ii)))
        enddo
       enddo
       do ii = 1,3
        do i = 1,neq
         id = nelmdg(i)-neqp1
         bp(i+nrhs(ii)) = bp(i+nrhs(ii))- 
     &     a(id+nmat(nh4(ii)))*a44i(i)*bp(i+nrhs(4))
        enddo
       enddo
c 
c now arrange arrays
c
       a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6))
       a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7))
       a(nmat(6)+1:nmat(7))  = a(nmat(7)+1:nmat(8))
       a(nmat(7)+1:nmat(8))  = a(nmat(9)+1:nmat(10))
       a(nmat(8)+1:nmat(9))  = a(nmat(10)+1:nmat(11))
       a(nmat(9)+1:nmat(10))  = a(nmat(11)+1:nmat(12))
c
c      note that the degree of freedom has been changed from 4 to 3
c
       idofm = 3
c
      else  if(iflg.eq.2) then
c extract full solution
        do i = 1,neq
         id = nelmdg(i)-neqp1
         bp(i+nrhs(4)) = a44i(i)*(bp(i+nrhs(4))- 
     &     a(id+nmat(13))*bp(i+nrhs(1)) -
     &     a(id+nmat(14))*bp(i+nrhs(2)) -
     &     a(id+nmat(15))*bp(i+nrhs(3)))
        enddo
c
c      note that the degree of freedom has been changed from 3 to 4
c
       idofm = 4
c
      endif
      return
      end


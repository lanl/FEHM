      subroutine co2_array(iflg,idofm)
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
c called from co2h20_combine
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
      integer neqp1,nmati
      integer nh11(16),nh12(16),nh13(16),nh14(4)
      integer nh21(9),nh22(9),nh23(9),nh24(3)
      real*8 tmp

      data nh11/1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 21, 22, 23, 25/
      data nh12/4, 4, 4, 4, 9, 9, 9,  9, 14, 14, 14, 14, 24, 24, 24, 24/
      data nh13/16,17,18,20,16,17,18,20, 16, 17, 18, 20, 16, 17, 18, 20/
      data nh14 /4, 9, 14, 24/
      data nh21  / 1,  2,  3,  5,  6,  7,  9, 10, 11 /
      data nh22  / 4,  4,  4,  8,  8,  8, 12, 12, 12 /
      data nh23  /13, 14, 15, 13, 14, 15, 13, 14, 15 /
      data nh24  /4, 8, 12 /	

c
      neqp1=neq+1
      if(iflg.eq.1) then  

c      nmatd is size of one subarray
       nmatd=nelm(neqp1)-neqp1
c
c combine parts of a array and bp array
c
       if(.not.allocated(a44i)) allocate(a44i(neq))
c     Inverse of a(nmat(24)+1:nmat(25)) is a44
       do i = 1,neq
          id = nelmdg(i)-neqp1
          if(idofm.eq.5) then
             nmati = 19
          else
             nmati = 16
          endif
          a44i(i) = 1.0/a(id+nmat(nmati))
       enddo
c     
       if(idofm.eq.5) then
          do ii = 1,(idofm-1)*(idofm-1)
             do i = 1,neq
                id = nelmdg(i)-neqp1
                a(id+nmat(nh11(ii))) = a(id+nmat(nh11(ii))) - 
     &               a(id+nmat(nh12(ii)))*a44i(i)*a(id+nmat(nh13(ii)))
             enddo
          enddo
       else
          do ii = 1,(idofm-1)*(idofm-1)
             do i = 1,neq
                id = nelmdg(i)-neqp1
                a(id+nmat(nh21(ii))) = a(id+nmat(nh21(ii))) - 
     &               a(id+nmat(nh22(ii)))*a44i(i)*a(id+nmat(nh23(ii)))
             enddo
          enddo
       endif

       if(idofm.eq.4) then
          do ii = 1,idofm-1
             do i = 1,neq
                id = nelmdg(i)-neqp1
                bp(i+nrhs(ii)) = bp(i+nrhs(ii))- 
     &               a(id+nmat(nh24(ii)))*a44i(i)*bp(i+nrhs(4))
             enddo
          enddo
       elseif(idofm.eq.5) then
          do ii = 1,3
             do i = 1,neq
                id = nelmdg(i)-neqp1
                bp(i+nrhs(ii)) = bp(i+nrhs(ii))- 
     &               a(id+nmat(nh14(ii)))*a44i(i)*bp(i+nrhs(4))
             enddo
          enddo
          do i = 1,neq
             id = nelmdg(i)-neqp1
             tmp = bp(i+nrhs(4))
             bp(i+nrhs(4)) = bp(i+nrhs(5))- 
     &            a(id+nmat(nh14(ii)))*a44i(i)*bp(i+nrhs(4))
             bp(i+nrhs(5)) = tmp
          enddo
       endif
c     
c     now arrange arrays
c     
       if(idofm.eq.5) then
          a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6))
          a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7))
          a(nmat(6)+1:nmat(7))  = a(nmat(7)+1:nmat(8))
          a(nmat(7)+1:nmat(8))  = a(nmat(8)+1:nmat(9))
          a(nmat(8)+1:nmat(9))  = a(nmat(10)+1:nmat(11))
          a(nmat(9)+1:nmat(10))  = a(nmat(11)+1:nmat(12))
          a(nmat(10)+1:nmat(11))  = a(nmat(12)+1:nmat(13))
          a(nmat(11)+1:nmat(12))  = a(nmat(13)+1:nmat(14))
          a(nmat(12)+1:nmat(13))  = a(nmat(15)+1:nmat(16))
          a(nmat(13)+1:nmat(14))  = a(nmat(16)+1:nmat(17))
          a(nmat(14)+1:nmat(15))  = a(nmat(17)+1:nmat(18))
          a(nmat(15)+1:nmat(16))  = a(nmat(18)+1:nmat(19))
          a(nmat(16)+1:nmat(17))  = a(nmat(20)+1:nmat(21))
c     change the degrees of freedom from 5 to 4
          idofm = 4
       elseif(idofm.eq.4) then
          a(nmat(4)+1:nmat(5))  = a(nmat(5)+1:nmat(6))
          a(nmat(5)+1:nmat(6))  = a(nmat(6)+1:nmat(7))
          a(nmat(6)+1:nmat(7))  = a(nmat(7)+1:nmat(8))
          a(nmat(7)+1:nmat(8))  = a(nmat(9)+1:nmat(10))
          a(nmat(8)+1:nmat(9))  = a(nmat(10)+1:nmat(11))
          a(nmat(9)+1:nmat(10))  = a(nmat(11)+1:nmat(12))
c     change the degrees of freedom from 4 to 3
          idofm = 3
       endif	 
c     

      else  if(iflg.eq.2) then
c     extract full solution
         if(idofm.eq.4) then
            do i = 1,neq
               tmp = bp(i+nrhs(4))
               id = nelmdg(i)-neqp1
               bp(i+nrhs(4)) = a44i(i)*(bp(i+nrhs(5))- 
     &              a(id+nmat(16))*bp(i+nrhs(1)) -
     &              a(id+nmat(17))*bp(i+nrhs(2)) -
     &              a(id+nmat(18))*bp(i+nrhs(3)) -
     &              a(id+nmat(20))*bp(i+nrhs(4)))
               bp(i+nrhs(5)) = tmp
            enddo
c     change the degree of freedom from 4 to 5
            idofm = 5
         elseif(idofm.eq.3) then
            do i = 1,neq
               id = nelmdg(i)-neqp1
               bp(i+nrhs(4)) = a44i(i)*(bp(i+nrhs(4))- 
     &              a(id+nmat(13))*bp(i+nrhs(1)) -
     &              a(id+nmat(14))*bp(i+nrhs(2)) -
     &              a(id+nmat(15))*bp(i+nrhs(3)))
            enddo
c     change the degree of freedom from 3 to 4
            idofm = 4
         endif
      endif
      return
      end


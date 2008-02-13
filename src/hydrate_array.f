      subroutine hydrate_array(iflg,idofm)
c
c simplifies arrays in methane hydrate calculations
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
      integer neqp1, j, i1, i2, jj, kb
      integer nh1(9), nh2(9), nh3(9), nh4(3)
      real*8, allocatable ::  bpmod(:)   

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
c allocate(a44i(neq))
c a(nmat(16)+1:nmat(17)) is a44
      if(.not.allocated(a44i)) then
       allocate(a44i(neq))
      endif
       do i = 1,neq
        id = nelmdg(i)-neqp1
        a44i(i) = 1.0/a(id+nmat(16))
       enddo

c first modify rhs
       do ii = 1,3
        do i = 1,neq
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         do j = i1,i2
          id = nelm(j)
          jj = j-neqp1
          bp(i+nrhs(ii)) = bp(i+nrhs(ii))- 
     &     a(jj+nmat(nh4(ii)))*a44i(id)*bp(id+nrhs(4))
         enddo
        enddo
       enddo

c next combine (forward substitution) appropriate matrices
       do ii = 1,9
        do i = 1,neq
         i1 = nelm(i)+1
         i2 = nelm(i+1)
         do j = i1,i2
          id = j-neqp1
          kb = nelm(j)
          jj = nelmdg(kb)-neqp1
          a(id+nmat(nh1(ii))) = a(id+nmat(nh1(ii))) - 
     &     a(id+nmat(nh2(ii)))*a44i(kb)*a(jj+nmat(nh3(ii)))
         enddo
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
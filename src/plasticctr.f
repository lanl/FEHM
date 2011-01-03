      subroutine plasticctr(iflg)

      use comai
      use comsi

      implicit none

      integer iflg
      integer i,j,k,iterm,jterm
      logical null1
      character*10 macro1
      double precision, dimension(1:3, 1:3) :: tmp
      
      if(iflg.eq.initPlastic) then

         if(icnl.eq.0) then
           allocate(row(1:3,1:3,1:9))
           allocate(col(1:3,1:3,1:9))
           tmp(1,1) = 1; tmp(1,2) = 4; tmp(1,3) = 6;
           tmp(2,1) = 4; tmp(2,2) = 2; tmp(2,3) = 5;
           tmp(3,1) = 6; tmp(3,2) = 5; tmp(3,3) = 3;
           do i=1,3
             do j=1,3
               do k=1,9
                 iterm = floor((k-0.5)/3) + 1
                 jterm = k - (iterm-1)*3
                 row(i,j,k) = tmp(i,jterm)
                 col(i,j,k) = tmp(iterm,j)
!                 write(iout,*) i,j,k,row(i,j,k), col(i,j,k)
               enddo
             enddo
           enddo
         else
             allocate(row(1:2,1:2,1:4))
             allocate(col(1:2,1:2,1:4))
         endif

      elseif (iflg.eq.assemblePlastic) then
!        if(icnl.eq.0) then
!          call geneq_plastic_3D()
!        endif
      endif

      end subroutine plasticctr

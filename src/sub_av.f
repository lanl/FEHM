      subroutine solve(neq,b,nop,irb,iirb,npvt,piv,nrhs1,nrhs2,
     *    idof,x1,x2)
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

!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/sub_av.f_a  $ 
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:18:22   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!***********************************************************************

      implicit none
      real*8 x1(*),x2(*)
      real*8 b(*),piv(*)
      integer nop(*),irb(*),iirb(*),npvt(*)
      integer neq, idof
      integer nrhs1(*),nrhs2(*)

      call renumber_array(x2,x1,iirb,neq,idof,nrhs2,nrhs1)
      call sub_bksub1(neq,b,x2,nop,npvt,piv)
      call renumber_array(x1,x2,irb,neq,idof,nrhs1,nrhs2)

      return
      end

      subroutine solven(neq,b,nb,nop,irb,iirb,npvt,piv,nrhs1,nrhs2,
     *    idof,x1,x2)

      implicit none
      real*8 x1(*),x2(*)
      real*8 b(*),piv(*)
      integer nb(*),nop(*),irb(*),iirb(*),npvt(*)
      integer neq, idof
      integer nrhs1(*),nrhs2(*)

      call renumber_array(x2,x1,iirb,neq,idof,nrhs2,nrhs1)
      if(idof.eq.1) then
        call sub_bksub1(neq,b,x2,nop,npvt,piv)
      else if(idof.eq.2) then
        call sub_bksub2(neq,b,x2,nb,nrhs2,nop,npvt,piv)
      else if(idof.eq.3) then
        call sub_bksub3(neq,b,x2,nb,nrhs2,nop,npvt,piv)
      else if(idof.eq.4) then
        call sub_bksub4(neq,b,x2,nb,nrhs2,nop,npvt,piv)
      else if(idof.eq.5) then
        call sub_bksub5(neq,b,x2,nb,nrhs2,nop,npvt,piv)
      else if(idof.eq.6) then
        call sub_bksub6(neq,b,x2,nb,nrhs2,nop,npvt,piv)
      else
        call sub_bksubn(neq,idof,b,x2,nb,nrhs2,nop,npvt,piv)
      endif
      call renumber_array(x1,x2,irb,neq,idof,nrhs1,nrhs2)

      return
      end


      subroutine mv(neq,a,ncon,x,y)

      implicit none
      real*8 a(*),x(*),y(*)
      integer ncon(*)
      integer neq

      integer i,j,j1,j2,neqp1
      real*8 sum

      neqp1 = neq + 1

      do i=1,neq
        sum = 0.0d00
        j1 = ncon(i) + 1
        j2 = ncon(i+1)
        do j = j1,j2
          sum = sum + a(j - neqp1) * x(ncon(j))
        enddo
        y(i) = sum
      enddo

      return
      end

      subroutine mvn(neq,idof,a,na,ncon,x,y,nrhs_x,nrhs_y)

      implicit none
      real*8 a(*),x(*),y(*)
      integer na(*),ncon(*),nrhs_x(*),nrhs_y(*)
      integer idof,neq

      integer i,j,j1,j2,js,jd1,jd2,kb,neqp1
      real*8 sum(idof),xs(idof),s

      neqp1 = neq + 1

      do i=1,neq
         do j=1,idof
            sum(j)=0.0d0
         end do
         j1=ncon(i)+1
         j2=ncon(i+1)
         do j=j1,j2
            kb=ncon(j)
            if( idof.eq.1 ) then
               sum(1)=sum(1)+a(j - neqp1)*x(kb)
            else
               do jd1=1,idof
                 xs(jd1) = x(kb + nrhs_x(jd1))
               end do
               js=j-neqp1
               kb=0
               do jd1=1,idof
                  s=0.0d0
                  do jd2=1,idof
                     kb=kb+1
                     s=s+a(js+na(kb))*xs(jd2)
                  end do
                  sum(jd1)=sum(jd1)+s
               end do
            end if
         end do
         if( idof.eq.1 ) then
            y(i) = sum(1)
         else
            do j=1,idof
               y(i+nrhs_y(j)) = sum(j)
            end do
         end if
      end do

      return
      end

      subroutine stressperm_91(i)               
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

      
      use comdi
      use comsi
      implicit none
      
      integer i,iispmd,itable, flag_stars
      real*8 biot,e1i,e2i,e3i,erat,efac,epi,shpi, shpi0
      real*8 eff_x,eff_y,eff_z,permfac_x,permfac_y,permfac_z
      real*8 fac_x,fac_y,fac_z,eff_x0,eff_y0,eff_z0, eff_strs
      
c hardwiring this for now, for compatibility with stars
      flag_stars = 1

      biot=bulk(i)
      e1i = e1(i)
      e2i = e2(i)
      e3i = e3(i)
      erat=e2i/e1i
      efac=3.d0*e2i+2.d0*e3i
      epi=efac*biot
c     dpd=phi(i)-phini(i)
      shpi0=epi*phini(i)
      shpi=epi*phi(i)
c      eff_x0=str_x0(i)-shpi0
c      eff_y0=str_y0(i)-shpi0
c      eff_z0=str_z0(i)-shpi0
      eff_x=str_x(i)-(shpi-shpi0)
      eff_y=str_y(i)-(shpi-shpi0)
      eff_z=str_z(i)-(shpi-shpi0)
c      eff_x=eff_x-eff_x0
c      eff_y=eff_y-eff_y0
c      eff_z=eff_z-eff_z0

      eff_strs=(eff_x+eff_y+eff_z)/3.

      permfac_x=1.0
      permfac_y=1.0
      permfac_z=1.0

      if(flag_stars.eq.0) then
         if(eff_x.le.k_strs91(1,1)) then
            permfac_x=k_strs91(1,2)
         else
            do itable=2,nentries
               if(eff_x.lt.k_strs91(itable,1)) then
                  fac_x=(k_strs91(itable,2)-k_strs91(itable-1,2))
     &                 /(k_strs91(itable,1)-k_strs91(itable-1,1))
                  permfac_x=(eff_x-k_strs91(itable-1,1))*fac_x
     &                 +k_strs91(itable-1,2)
                  goto 9193
               endif
            enddo
            permfac_x=k_strs91(nentries,2)
 9193       continue
         endif
         if(eff_y.le.k_strs91(1,1)) then
            permfac_Y=k_strs91(1,3)
         else
            do itable=2,nentries
               if(eff_y.lt.k_strs91(itable,1)) then
                  fac_y=(k_strs91(itable,3)-k_strs91(itable-1,3))
     &                 /(k_strs91(itable,1)-k_strs91(itable-1,1))
                  permfac_y=(eff_y-k_strs91(itable-1,1))*fac_y
     &                 +k_strs91(itable-1,3)
                  goto 9194
               endif
            enddo
            permfac_Y=k_strs91(nentries,3)
 9194       continue
         endif
         if(eff_z.le.k_strs91(1,1)) then
            permfac_z=k_strs91(1,4)
         else
            do itable=2,nentries
               if(eff_z.lt.k_strs91(itable,1)) then
                  fac_z=(k_strs91(itable,4)-k_strs91(itable-1,4))
     &                 /(k_strs91(itable,1)-k_strs91(itable-1,1))
                  permfac_z=(eff_z-k_strs91(itable-1,1))*fac_z
     &                 +k_strs91(itable-1,4)
                  goto 9195
               endif
            enddo
            permfac_z=k_strs91(nentries,4)
 9195       continue
         endif
      elseif(flag_stars.eq.1) then
         if(eff_strs.le.k_strs91(1,1)) then
            permfac_x=k_strs91(1,2)
            permfac_y=k_strs91(1,3)
            permfac_z=k_strs91(1,4)
         else
            do itable=2,nentries
               if(eff_strs.lt.k_strs91(itable,1)) then
                  fac_x=(k_strs91(itable,2)-k_strs91(itable-1,2))
     &                 /(k_strs91(itable,1)-k_strs91(itable-1,1))
                  fac_y=(k_strs91(itable,3)-k_strs91(itable-1,3))
     &                 /(k_strs91(itable,1)-k_strs91(itable-1,1))
                  fac_z=(k_strs91(itable,4)-k_strs91(itable-1,4))
     &                 /(k_strs91(itable,1)-k_strs91(itable-1,1))
                  permfac_x=(eff_strs-k_strs91(itable-1,1))*fac_x
     &                 +k_strs91(itable-1,2)
                  permfac_y=(eff_strs-k_strs91(itable-1,1))*fac_y
     &                 +k_strs91(itable-1,3)
                  permfac_z=(eff_strs-k_strs91(itable-1,1))*fac_z
     &                 +k_strs91(itable-1,4)
                  goto 9197
               endif
            enddo
            permfac_x=k_strs91(nentries,2)
            permfac_Y=k_strs91(nentries,3)
            permfac_z=k_strs91(nentries,4)
 9197       continue
            
         endif
      endif
         
      pnx(i) = pnx0(i)*permfac_x
      pny(i) = pny0(i)*permfac_y
      pnz(i) = pnz0(i)*permfac_z

      if(permfac_x.ne.1.0.or.permfac_y.ne.1.0.or.
     &     permfac_z.ne.1.0) then
c      write(97,9791)i,eff_x0,eff_y0,eff_z0,
      write(97,9791)i,eff_x,eff_y,eff_z,permfac_x,permfac_y,permfac_z
 9791 format(i8,9(2x,f8.3))
      endif

      return
      end

c.....................................................................

      subroutine velocity_derivatives_lsq(inp1,
     $        dvd,dvxd,dvyd,dvzd,v,vv)
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1 
!D1 PURPOSE
!D1 
!D1 Compute first derivatives of velocity on a general grid using least
!D1 squares.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 07-Mar-05, Programmer: S Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/velocity_derivatives_lsq.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************
c s kelkar March 7 05.First derivatives of velocity on a general grid
c using least squares. Going up to second neighbours of the node i
c points weighted by inverse of distance.

      use comai, only : ierr, iptty, neq
      use combi, only : cord, nelm, istrw, isox,isoy,isoz, sx
      use comci, only : rolf
      use comdi, only : ps
      use comflow, only : a_axy
      use comsptr
      use comsk

      implicit none

      integer inp1,i,j,i1,i2,j1,j2,ipos,jj,jjj,ndone,kdone(1000)
      integer ij,ikb,j3,kb,iwsk,indx(8),ludflag, icount

      real*8 axy,uij,dludsk,epsilon
      real*8 dvd(3),dvxd(3),dvyd(3),dvzd(3),v(3),vv
      real*8 aarea,area,distij,sxkb,xn(3),term(8,8),rhs(8),fac(8)
      real*8 dx(3),fad(8),distikb,distij2,distkb2, rolfup

      epsilon = 1.e-22

      ndone=1
      kdone=0
      term=0.
      rhs=0.
      fac=0.
      fad=0.
      icount = 0

      kdone(1)=inp1
      
      i1=nelm(inp1)+1
      i2=nelm(inp1+1)

c loop over 1-st neighbours
      do j=i1,i2
         kb=nelm(j)
         if(kb.ne.inp1.and.ps(kb).gt.0.) then
            icount=icount+1
            ndone=ndone+1
            kdone(ndone)=kb
            ipos=j-neq-1
            axy=a_axy(ipos)
            iwsk=istrw(ipos)
            sxkb=sx(iwsk,isox)+sx(iwsk,isoy)+sx(iwsk,isoz)
            if(abs(sxkb).gt.epsilon) then
c     incrimental to the midpoint of i,j; and
c     distance between inp1 and kb, needed for calculating area
               distij2=0.
               do i=1,3
                  dx(i)=0.5*(cord(kb,i)-cord(inp1,i))
                  distij2=distij2+4.*dx(i)*dx(i)
               enddo
               distij=sqrt(distij2)
               aarea=-sxkb*distij
               area=abs(aarea)
c     direction cosins of the line joining inp1 and kb
               xn=2.*dx/distij
c     velocity from FEHM-mass flux accross the face
c     get upwinded density
               rolfup=rolf(kb)
               if(axy.gt.0) rolfup=rolf(inp1)
               uij=axy/area/rolfup
c     
               fac(1)=xn(1)
               fac(2)=xn(2)
               fac(3)=xn(3)
               fac(4)=xn(1)*dx(1)-xn(3)*dx(3)
               fac(5)=xn(1)*dx(2)+xn(2)*dx(1)
               fac(6)=xn(1)*dx(3)+xn(3)*dx(1)
               fac(7)=xn(2)*dx(2)-xn(3)*dx(3)
               fac(8)=xn(2)*dx(3)+xn(3)*dx(2)
c     distance weighting
c
               fad=fac/distij2
c               fad=fac
c     coefficients for  equations
               do jj=1,8
                  rhs(jj)=rhs(jj)+fad(jj)*uij
                  do jjj=1,8
                     term(jj,jjj)=term(jj,jjj)+fad(jj)*fac(jjj)
                  enddo
               enddo
            endif
         endif
      enddo

c     loop over 2-nd neighbours
      do ij=i1,i2
         ikb=nelm(ij)
         if(ikb.ne.inp1.and.ps(ikb).gt.0.) then
            j1=nelm(ikb)+1
            j2=nelm(ikb+1)
            do j=j1,j2
               kb=nelm(j)
               if(ps(kb).gt.0.) then
                  icount=icount+1
                  do j3=1,ndone
                     if(kb.eq.kdone(j3)) goto 99999
                  enddo
                  ndone=ndone+1
                  kdone(ndone)=kb
                  ipos=j-neq-1
c     note: this is the normal flux between the nodes ikb and kb
                  axy=a_axy(ipos)
                  iwsk=istrw(ipos)
                  sxkb=sx(iwsk,isox)+sx(iwsk,isoy)+sx(iwsk,isoz)
                  if(abs(sxkb).gt.epsilon) then
c     distance between ikb and kb, needed for calculating area
                     distikb=0.
                     do i=1,3
                        distikb=distikb+(cord(ikb,i)-cord(kb,i))**2.
                     enddo
                     distikb=sqrt(distikb)
                     aarea=-sxkb*distikb
                     area=abs(aarea)
c     velocity from FEHM-mass flux accross the face
c     get upwinded density
                     rolfup=rolf(kb)
                     if(axy.gt.0) rolfup=rolf(inp1)
                     uij=axy/area/rolfup
c     incrimental to the midpoint of ikb,kb, wrt inp1, and
c     distance between inp1 and the midpoint of ijk-kb
                     distij2=0.
                     do i=1,3
                        dx(i)=0.5*(cord(kb,i)+cord(ikb,i))-cord(inp1,i)
                        distij2=distij2+dx(i)*dx(i)
                     enddo
                     distij=sqrt(distij2)
c     direction cosins of the normal between ikb and kb
c     note that the flux given in a-axy is the projection along this 
c     direction, and not along the line inp1-mid(ikb-kb)
                     do i=1,3
                        xn(i)=(cord(kb,i)-cord(ikb,i))/distikb
                     enddo
c     
                     fac(1)=xn(1)
                     fac(2)=xn(2)
                     fac(3)=xn(3)
                     fac(4)=xn(1)*dx(1)-xn(3)*dx(3)
                     fac(5)=xn(1)*dx(2)+xn(2)*dx(1)
                     fac(6)=xn(1)*dx(3)+xn(3)*dx(1)
                     fac(7)=xn(2)*dx(2)-xn(3)*dx(3)
                     fac(8)=xn(2)*dx(3)+xn(3)*dx(2)
c     distance weighting
                     fad=fac/distij2
c     fad=fac
c     coefficients for  equations
                     do jj=1,8
                        rhs(jj)=rhs(jj)+fad(jj)*uij
                        do jjj=1,8
                           term(jj,jjj)=term(jj,jjj)+fad(jj)*fac(jjj)
                        enddo
                     enddo
                  endif
99999             continue
               endif
            enddo
         endif
      enddo

c may 20 2009 porosity check
      if(icount.gt.0) then
c     LU decomposition on the matrix term(8,8)
         call ludsk(term,8,8,indx,dludsk,ludflag)
c     stop if singular matrix in ludsk
         if(ludflag.gt.0 ) then
            write(ierr,*) 'stopvelocity_derivatives_lsq.'  
            write(ierr,*)'singular matrix in ludsk'
            write(ierr,*)'node #, row # = ',i,ludflag
            write(ierr,*) 'irray=', (irray(i,i1), i1=-3,3)
            stop
         endif
c     back substitute to solve the equation term*Vcoef=rhs for vcoef
c     where vcoef(1-8) has the coefficient for interpolated velocity
c     and derivatives contained in rhs upon return
         call lubsk(term,8,8,indx,rhs)
         
c     form velocity and its derivatives
         
         v(1)=rhs(1)
         v(2)=rhs(2)
         v(3)=rhs(3)
         
         dvxd(1)=rhs(4)
         dvxd(2)=rhs(5)
         dvxd(3)=rhs(6)
         
         dvyd(2)=rhs(7)
         dvyd(3)=rhs(8)
         
         dvyd(1)=dvxd(2)
         
         dvzd(1)=dvxd(3)
         dvzd(2)=dvyd(3)
         dvzd(3)=-(dvxd(1)+dvyd(2))
         
         vv=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         
         do i=1,3
            dvd(i)=(v(1)*dvxd(i)+v(2)*dvyd(i)+v(3)*dvzd(i))
         enddo
         dvd=dvd/vv
      else
         write(ierr,*)'STOP.velocity_derivatives_lsq. node inp1=',inp1
         write(ierr,*)' has no positive porosity neighbours. STOP.'
         if (iptty .ne. 0) then
            write(iptty,*)'STOP.velocity_derivatives_lsq. node inp1=',
     &           inp1
            write(iptty,*)' has no positive porosity neighbours. STOP.'
         end if
         stop
      endif

      return

      end

c.....................................................................

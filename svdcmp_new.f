      SUBROUTINE svdcmp_new(ierr,a,m,n,mp,np,w,v)
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
!D1 To use singular value decomposition to determine weights of nearest
!D1 curves for delay time interpolation.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 
!D2 Initial implementation: 15-DEC-02, Programmer: Bruce Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/svdcmp_new.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
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

      implicit none

      integer ierr
      INTEGER m,mp,n,np,NMAX
      REAL a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500)
      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
      g=0.0
      scale=0.0
      anorm=0.0
      do i=1,n
         l=i+1
         rv1(i)=scale*g
         g=0.0
         s=0.0
         scale=0.0
         if(i.le.m)then
            do k=i,m
               scale=scale+abs(a(k,i))
            enddo
            if(scale.ne.0.0)then
               do k=i,m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
               enddo
               f=a(i,i)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,i)=f-g
               do j=l,n
                  s=0.0
                  do k=i,m
                     s=s+a(k,i)*a(k,j)
                  enddo
                  f=s/h
                  do k=i,m
                     a(k,j)=a(k,j)+f*a(k,i)
                  enddo
               enddo
               do k=i,m
                  a(k,i)=scale*a(k,i)
               enddo
            endif
         endif
         w(i)=scale *g
         g=0.0
         s=0.0
         scale=0.0
         if((i.le.m).and.(i.ne.n))then
            do k=l,n
               scale=scale+abs(a(i,k))
            enddo
            if(scale.ne.0.0)then
               do k=l,n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
               enddo
               f=a(i,l)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,l)=f-g
               do k=l,n
                  rv1(k)=a(i,k)/h
               enddo
               do j=l,m
                  s=0.0
                  do k=l,n
                     s=s+a(j,k)*a(i,k)
                  enddo
                  do k=l,n
                     a(j,k)=a(j,k)+s*rv1(k)
                  enddo
               enddo
               do k=l,n
                  a(i,k)=scale*a(i,k)
               enddo
            endif
         endif
         anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      enddo
      do i=n,1,-1
         if(i.lt.n)then
            if(g.ne.0.0)then
               do j=l,n
                  v(j,i)=(a(i,j)/a(i,l))/g
               enddo
               do j=l,n
                  s=0.0
                  do k=l,n
                     s=s+a(i,k)*v(k,j)
                  enddo
                  do k=l,n
                     v(k,j)=v(k,j)+s*v(k,i)
                  enddo
               enddo
            endif
            do j=l,n
               v(i,j)=0.0
               v(j,i)=0.0
            enddo
         endif
         v(i,i)=1.0
         g=rv1(i)
         l=i
      enddo
      do i=min(m,n),1,-1
         l=i+1
         g=w(i)
         do j=l,n
            a(i,j)=0.0
         enddo
         if(g.ne.0.0)then
            g=1.0/g
            do j=l,n
               s=0.0
               do k=l,m
                  s=s+a(k,i)*a(k,j)
               enddo
               f=(s/a(i,i))*g
               do k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
               enddo
            enddo
            do j=i,m
               a(j,i)=a(j,i)*g
            enddo
         else
            do j= i,m
               a(j,i)=0.0
            enddo
         endif
         a(i,i)=a(i,i)+1.0
      enddo
      do k=n,1,-1
         do its=1,30
            do l=k,1,-1
               nm=l-1
               if((abs(rv1(l))+anorm).eq.anorm) goto 2
               if((abs(w(nm))+anorm).eq.anorm) goto 1
            enddo
 1          c=0.0
            s=1.0
            do i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               if((abs(f)+anorm).eq.anorm) goto 2
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=1.0/h
               c = (g*h)
               s=-(f*h)
               do j=1,m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
               enddo
            enddo
 2          z=w(k)
            if(l.eq.k)then
               if(z.lt.0.0)then
                  w(k)=-z
                  do j=1,n
                     v(j,k)=-v(j,k)
                  enddo
               endif
               goto 3
            endif
            if(its.eq.30) then
               write(ierr,*)'Stopping in svdcmp_new'
               write(ierr,*)'Fatal error in transfer function'
               stop
            end if
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=pythag(f,1.0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c =1.0
            s=1.0
            do j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=pythag(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f= (x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               do jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)= (x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
               enddo
               z=pythag(f,h)
               w(j)=z
               if(z.ne.0.0)then
                  z=1.0/z
                  c=f*z
                  s=h*z
               endif
               f= (c*g)+(s*y)
               x=-(s*g)+(c*y)
               do jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)= (y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
               enddo
            enddo
            rv1(l)=0.0
            rv1(k)=f
            w(k)=x
         enddo
 3       continue
      enddo
      return
      END



      FUNCTION pythag(a,b)
      REAL a,b,pythag
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
         pythag=absa*sqrt(1.+(absb/absa)**2)
      else
         if(absb.eq.0.)then
            pythag=0.
         else
            pythag=absb*sqrt(1.+(absa/absb)**2)
         endif
      endif
      return
      END


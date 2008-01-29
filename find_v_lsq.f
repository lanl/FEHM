      subroutine find_v_lsq(i,xxx)
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
!D1 To perform velocity calculation.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 14-Jan-02, Programmer: S Kelkar
!D2
!D2 $Log:   /pvcs.config/fehm90/src/find_v_lsq.f_a  $
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.6 Streamline particle-tracking module
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4  Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************
c... s kelkar  march2 04
c modified this subroutine to do a least square on 6 velocity 
c parameters needed to get polack velocities in the approximate 
c brick shape by matching against the fehm normal velocities 
c components. 

      use comai
      use combi
      use comflow
      use comsptr
      use comsk

      implicit none

      integer i,j,k,ll,ikb,kb,ipos,iwsk,flag_dim,i1,i2
      integer ibou, iboui, i31, i32, i41, i42, i51, i52
      integer indx3(3), indx4(4), indx5(5), indx6(6)
      integer idelete(3), ireplace(3), ibou_temp(7)

      real*8 sxkb,axy,aarea
      real*8 area,uu,area2
      real*8 epsilon
c...  s kelkar 3d OMR    2/26/04
      integer ii1,ii2,ludflag
      real*8 rarea,xxx
      real*8 psierr, vdotc, dludsk
      real*8 a3(3,3),b3(3)
      real*8 a4(4,4),b4(4)
      real*8 a5(5,5),b5(5)
      real*8 a6(6,6),b6(6)
      real*8 xm,ym,zm,rx,ry,rz,distij,xn,yn,zn
      real*8 avx1,avx2,avy1,avy2,avz1,avz2
      real*8 epomr,axydg
      real*8 dmass_5lsq, dmass_6lsq
      real*8 epsilon_diag, a6_max
c..................................

c.... s kelkar  march 10, 04, 3D ORM stuff............
c     irray(i,0) = +i : regular interior node, not a source/sink
c     irray(i,0) = -i : regular node that is a sink/source
c     irray(i,0) = 0 :  OMR node not on boundary
c     irray(i,0) = -i-1000 : OMR node on a external boundary
c.........................................................

      epsilon = 1.e-10
      epomr=1.e-22
      epsilon_diag=1.e-6

c......................................................
c...  s kelkar  march2 04
c     modified this subroutine to do a least square on 6 velocity 
c     parameters needed to get polack velocities in the approximate 
c     brick shape by matching against the fehm normal velocities 
c     components.
c     convention used in a6 and b6 matrices is
c     index 1: Vx1, ie x-velocity on the i3=-1 (left) face 
c     index 2: Vx2, ie x-velocity on the i3=+1 (right)face 
c     index 3: Vy1, ie y-velocity on the i3=-2 (front)face 
c     index 4: Vy2, ie y-velocity on the i3=+2 (back) face 
c     index 5: Vz1, ie z-velocity on the i3=-3 (bottom)face 
c     index 6: Vz2, ie z-velocity on the i3=+3 (top) face 

      if (iboulist(i,7) .ge. 3) return
c     zero the coefficient matrices
      a6=0.
      b6=0.
c     
      ii1=nelm(i)+1
      ii2=nelm(i+1)
      do i1=ii1,ii2
         ipos=i1-(neq+1)
         axy=a_axy(ipos)
         kb=nelm(i1)
         if(kb.ne.i) then
c     istrw is a pointer array for sx, corrosponding to the connection
c     i-i1. Generally(but not always) it is a scalar, the negative of the 
c     magnitude of (cross sectional area0/3. divided by the inter-nodal 
c     distance
            iwsk=istrw(i1-neq-1)
c     define the midpoint between i-kb
            xm=0.5*(cord(i,1)+cord(kb,1))
            ym=0.5*(cord(i,2)+cord(kb,2))
            zm=0.5*(cord(i,3)+cord(kb,3))
c     local coordinates of the center point
            rx=(xm-corn(i,1))/ddx(i)
            ry=(ym-corn(i,2))/ddy(i)
            rz=(zm-corn(i,3))/ddz(i)
c     distance between i and kb, needed for calculating area
            distij=0.
            do j=1,3
               distij=distij+(cord(i,j)-cord(kb,j))**2.
            enddo
            distij=sqrt(distij)
            sxkb=sx(iwsk,isox)+sx(iwsk,isoy)+sx(iwsk,isoz)
            aarea=-sxkb*distij
            area=abs(aarea)
            rarea=area*xxx
c     direction cosins of the normal to the control voulme face
            xn=(cord(kb,1)-cord(i,1))/distij
            yn=(cord(kb,2)-cord(i,2))/distij
            zn=(cord(kb,3)-cord(i,3))/distij
c     
            avx1=rarea*xn*(1.-rx)
            avx2=rarea*xn*rx
            avy1=rarea*yn*(1.-ry)
            avy2=rarea*yn*ry
            avz1=rarea*zn*(1.-rz)
            avz2=rarea*zn*rz
c     
            b6(1)=b6(1)+avx1*axy
            b6(2)=b6(2)+avx2*axy
            b6(3)=b6(3)+avy1*axy
            b6(4)=b6(4)+avy2*axy
            b6(5)=b6(5)+avz1*axy
            b6(6)=b6(6)+avz2*axy
c     
            avx1=rarea*avx1
            avx2=rarea*avx2
            avy1=rarea*avy1
            avy2=rarea*avy2
            avz1=rarea*avz1
            avz2=rarea*avz2
c     least sq coeff in d(psi-mass)/dvx1=0.
            a6(1,1)=a6(1,1)+avx1*xn*(1.-rx)
            a6(1,2)=a6(1,2)+avx1*xn*rx
            a6(1,3)=a6(1,3)+avx1*yn*(1.-ry)
            a6(1,4)=a6(1,4)+avx1*yn*ry
            a6(1,5)=a6(1,5)+avx1*zn*(1.-rz)
            a6(1,6)=a6(1,6)+avx1*zn*rz
c     least sq coeff in d(psi-mass)/dvx2=0.
c           a6(2,1)=a6(1,2)
            a6(2,2)=a6(2,2)+avx2*xn*rx
            a6(2,3)=a6(2,3)+avx2*yn*(1.-ry)
            a6(2,4)=a6(2,4)+avx2*yn*ry
            a6(2,5)=a6(2,5)+avx2*zn*(1.-rz)
            a6(2,6)=a6(2,6)+avx2*zn*rz
c     least sq coeff in d(psi-mass)/dvy1=0.
c           a6(3,1)=a6(1,3)
c           a6(3,2)=a6(2,3)
            a6(3,3)=a6(3,3)+avy1*yn*(1.-ry)
            a6(3,4)=a6(3,4)+avy1*yn*ry
            a6(3,5)=a6(3,5)+avy1*zn*(1.-rz)
            a6(3,6)=a6(3,6)+avy1*zn*rz
c     least sq coeff in d(psi-mass)/dvy2=0.
c           a6(4,1)=a6(1,4)
c           a6(4,2)=a6(2,4)
c           a6(4,3)=a6(3,4)
            a6(4,4)=a6(4,4)+avy2*yn*ry
            a6(4,5)=a6(4,5)+avy2*zn*(1.-rz)
            a6(4,6)=a6(4,6)+avy2*zn*rz
c     least sq coeff in d(psi-mass)/dvz1=0.
c           a6(5,1)=a6(1,5)
c           a6(5,2)=a6(2,5)
c           a6(5,3)=a6(3,5)
c           a6(5,4)=a6(4,5)
            a6(5,5)=a6(5,5)+avz1*zn*(1.-rz)
            a6(5,6)=a6(5,6)+avz1*zn*rz
c     least sq coeff in d(psi-mass)/dvz2=0.
c           a6(6,1)=a6(1,6)
c           a6(6,2)=a6(2,6)
c           a6(6,3)=a6(3,6)
c           a6(6,4)=a6(4,6)
c           a6(6,5)=a6(5,6)
            a6(6,6)=a6(6,6)+avz2*zn*rz

         endif

      enddo

c     the coefficient matrix LSQ is symmetric
      do i1=2,6
         do kb=1,i1-1
            a6(i1,kb)=a6(kb,i1)
         enddo
      enddo
c...  s kelkar June 27 06...................................
c     put small numbers on diagonal if it has 0 value
      a6_max=0.d0
      do i1=1,6
         do kb=1,6
            if(a6_max.lt.abs(a6(i1,kb))) a6_max=abs(a6(i1,kb))
         enddo
      enddo
      a6_max=epsilon_diag*a6_max
      do i1=1,6
         if(abs(a6(i1,i1)).lt.a6_max) a6(i1,i1)=a6_max
      enddo
c..................................................................
      do j = 1, 7
         ibou_temp(j) = iboulist(i,j)
      end do

 100  iboui = ibou_temp(7)
      idelete = 0
      ireplace = 0
      do j = 1, iboui
         ibou=ibou_temp(j)
         select case (ibou)
         case (-1)
            idelete(j) = 1
            ireplace(j)= 2
         case (+1)
            idelete(j) = 2
            ireplace(j)= 1
         case (-2)
            idelete(j) = 3
            ireplace(j)= 4
         case (+2)
            idelete(j) = 4
            ireplace(j)= 3
         case (-3)
            idelete(j) = 5
            ireplace(j)= 6
         case (+3)
            idelete(j) = 6
            ireplace(j)= 5
         end select
      end do

      ludflag = 0
      select case (iboui)
      case (0)
         call ludsk(a6,6,6,indx6,dludsk,ludflag)
         if (ludflag .gt. 0) then
            call ludsk_error
            goto 100
         end if
         call lubsk(a6,6,6,indx6,b6)
      case (1)
         i51=0
         do i1=1,6
            if(i1.ne.idelete(1)) then
               i51=i51+1
               b5(i51)=b6(i1)
               i52=0
               do i2=1,6
                  if(i2.ne.idelete(1)) then
                     i52=i52+1
                     a5(i51,i52)=a6(i1,i2)
                  endif
               enddo
            endif
         enddo
         call ludsk(a5,5,5,indx5,dludsk,ludflag)
         if (ludflag .gt. 0) then
            call ludsk_error
            goto 100
         end if
         call lubsk(a5,5,5,indx5,b5)
c reassemble the b6 vector
         i51=0
         do i1=1,6
            if(i1.ne.idelete(1)) then
               i51=i51+1
               b6(i1)=b5(i51)
            endif
         enddo
         b6(idelete(1))=b6(ireplace(1))
      case (2) 
         i41=0
         do i1=1,6
            if(i1.ne.idelete(1).and.i1.ne.idelete(2)) then
               i41=i41+1
               b4(i41)=b6(i1)
               i42=0
               do i2=1,6
                  if(i2.ne.idelete(1).and.i2.ne.idelete(2)) then
                     i42=i42+1
                     a4(i41,i42)=a6(i1,i2)
                  endif
               enddo
            endif
         enddo
         call ludsk(a4,4,4,indx4,dludsk,ludflag)
         if (ludflag .gt. 0) then
            call ludsk_error
            goto 100            
         end if
         call lubsk(a4,4,4,indx4,b4)
c reassemble the b6 vector
         i41=0
         do i1=1,6
            if(i1.ne.idelete(1).and.i1.ne.idelete(2)) then
               i41=i41+1
               b6(i1)=b4(i41)
            endif
         enddo
         b6(idelete(1))=b6(ireplace(1))
         b6(idelete(2))=b6(ireplace(2))
      case (3)
         i31 = 0
         do i1=1,6
            if(i1.ne.idelete(1).and.i1.ne.idelete(2).and.
     &           i1.ne.idelete(3)) then
               i31=i31+1
               b3(i31)=b6(i1)
               i32=0
               do i2=1,6
                  if(i2.ne.idelete(1).and.i2.ne.idelete(2).and.
     &                 i2.ne.idelete(3)) then
                     i32=i32+1
                     a4(i31,i32)=a6(i1,i2)
                  endif
               enddo
            endif
         enddo
         call ludsk(a3,3,3,indx3,dludsk,ludflag)
         if (ludflag .gt. 0) goto 200
         call lubsk(a3,3,3,indx3,b3)
c reassemble the b6 vector
         i31=0
         do i1=1,6
            if(i1.ne.idelete(1).and.i1.ne.idelete(2).and.
     &           i1.ne.idelete(3)) then
               i31=i31+1
               b6(i1)=b3(i31)
            endif
         enddo
         b6(idelete(1))=b6(ireplace(1))
         b6(idelete(2))=b6(ireplace(2))
         b6(idelete(3))=b6(ireplace(3))
      end select

c put velocities into the ggg array the +- signs are to make it 
c consistant with the conventions in set_velocities_struct.
      ggg(i,-1)=-b6(1)
      ggg(i,+1)=+b6(2)
      ggg(i,-2)=-b6(3)
      ggg(i,+2)=+b6(4)
      ggg(i,-3)=-b6(5)
      ggg(i,+3)=+b6(6)

c correct ggg for source/sink nodes
c divide by iboui to partition the flux between the boundary faces
      iboui = iboulist(i,7)
      axydg=a_axy(nelmdg(i)-1-neq)
      do j = 1, iboui
c     if(abs(axydg).gt.epomr) then
         ibou = iboulist(i,j)
         call getarea(i,ibou,area)
         ggg(i,iboulist(i,j))=axydg/(xxx*area*iboui)
c     end if
      end do

      return

 200  if(ludflag.gt.0 ) then
         write(ierr,*) 'stop find_v_lsq singular matrix in ludsk'
         write(ierr,*)'node #, row # = ',i,ludflag
         write(ierr,*) 'irray=', (irray(i,i1), i1=-3,3)
         write(ierr,*) 'iboulist=', (iboulist(i,i1), i1=1,4)
         write(ierr,*) 'ibou_temp=', (ibou_temp(i1), i1=1,4)
         write(ierr,*) 'b6=', b6
         write(ierr,*) 'a6=', ((a6(i1,i2), i1=1,6), i2=1,6)
         if (iptty .ne. 0) then
            write(iptty,*) 'stop find_v_lsq singular matrix in ludsk'
            write(iptty,*)'node #, row # = ',i,ludflag
            write(iptty,*) 'irray=', (irray(i,i1), i1=-3,3)
            write(iptty,*) 'iboulist=', (iboulist(i,i1), i1=1,4)
            write(iptty,*) 'ibou_temp=', (ibou_temp(i1), i1=1,4)
            write(iptty,*) 'b6=', b6
            write(iptty,*) 'a6=', ((a6(i1,i2), i1=1,6), i2=1,6)
         end if
         stop
      endif

      contains

      subroutine ludsk_error

      implicit none 

      ibou_temp(7) = iboui + 1
      select case (iboui)
      case (0)
         if (ludflag .eq. 1) then
            ibou_temp(1) = -1
         else if (ludflag .eq. 2) then
            ibou_temp(1) = 1
         else if (ludflag .eq. 3) then
            ibou_temp(1) = -2
         else if (ludflag .eq. 4) then
            ibou_temp(1) = 2
         else if (ludflag .eq. 5) then
            ibou_temp(1) = -3
         else if (ludflag .eq. 6) then
            ibou_temp(1) = 3
         end if
      case (1)
         if (ludflag .eq. 1) then
            if (idelete(1) .eq. 1) then
               ibou_temp(2) = +1
            else
               ibou_temp(2) = -1
            end if
         else if (ludflag .eq. 2) then
            if (idelete(1) .le. 2) then
               ibou_temp(2) = -2
            else
               ibou_temp(2) = +1
            end if
         else if (ludflag .eq. 3) then
            if (idelete(1) .le. 3) then
               ibou_temp(2) = +2
            else
               ibou_temp(2) = -2
            end if
         else if (ludflag .eq. 4) then
            if (idelete(1) .le. 4) then
               ibou_temp(2) = -3
            else
               ibou_temp(2) = +2
            end if
         else if (ludflag .eq. 5) then
            if (idelete(1) .le. 5) then
               ibou_temp(2) = +3
            else
               ibou_temp(2) = -3
            end if
         end if
      case (2)
         if (idelete(1) .eq. 1 .or. idelete(2) .eq. 1) then
            if (ludflag .eq. 1) then
               if (idelete(1) .eq. 2 .or. idelete(2) .eq. 2) then
                  ibou_temp(3) = -2
               else
                  ibou_temp(3) = +1
               end if
            else if (ludflag .eq. 2) then
               if (idelete(1) .le. 3 .or. idelete(2) .le. 3) then
                  ibou_temp(3) = +2
               else
                  ibou_temp(3) = -2
               end if
            else if (ludflag .eq. 3) then
               if (idelete(1) .le. 4 .or. idelete(2) .le. 4) then
                  ibou_temp(3) = -3
               else
                  ibou_temp(3) = +2
               end if
            else if (ludflag .eq. 4) then
               if (idelete(1) .eq. 6 .or. idelete(2) .eq. 6) then
                  ibou_temp(3) = -3
               else
                  ibou_temp(3) = +3
               end if
            end if
         else if (idelete(1) .eq. 2 .or. idelete(2) .eq. 2) then
            if (ludflag .eq. 1) then
               ibou_temp(3) = -1
            else if (ludflag .eq. 2) then
               if (idelete(1) .eq. 3 .or. idelete(2) .eq. 3) then
                  ibou_temp(3) = +2
               else
                  ibou_temp(3) = -2
               end if
            else if (ludflag .eq. 3) then
               if (idelete(1) .le. 4 .or. idelete(2) .le. 4) then
                  ibou_temp(3) = -3
               else
                  ibou_temp(3) = +2
               end if
            else if (ludflag .eq. 4) then
               if (idelete(1) .eq. 6 .or. idelete(2) .eq. 6) then
                  ibou_temp(3) = -3
               else
                  ibou_temp(3) = +3
               end if
            end if
         else if (idelete(1) .eq. 3 .or. idelete(2) .eq. 3) then
            if (ludflag .eq. 1) then
               ibou_temp(3) = -1
            else if (ludflag .eq. 2) then
               ibou_temp(3) = +1
            else if (ludflag .eq. 3) then
               if (idelete(1) .eq. 4 .or. idelete(2) .eq. 4) then
                  ibou_temp(3) = -3
               else
                  ibou_temp(3) = +2
               end if
            else if (ludflag .eq. 4) then
               if (idelete(1) .eq. 6 .or. idelete(2) .eq. 6) then
                  ibou_temp(3) = -3
               else
                  ibou_temp(3) = +3
               end if
            end if
         else if (idelete(1) .ge. 4 .or. idelete(2) .ge. 4 ) then
            if (ludflag .eq. 1) then
               ibou_temp(3) = -1
            else if (ludflag .eq. 2) then
               ibou_temp(3) = +1
            else if (ludflag .eq. 3) then
               ibou_temp(3) = -2
            else if (ludflag .eq. 4) then
               if (idelete(1) .eq. 4 .or. idelete(2) .eq. 4 ) then
                  if (idelete(1) .eq. 6 .or. idelete(2) .eq. 6) then
                     ibou_temp(3) = -3
                  else
                     ibou_temp(3) = +3
                  end if
               else if (idelete(1) .eq. 5 .or. idelete(2) .eq. 5 ) then
                  ibou_temp(3) = +2
               end if
            end if
         end if
      end select

      end subroutine ludsk_error

      end subroutine find_v_lsq

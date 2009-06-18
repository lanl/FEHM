      subroutine porosity_gradient_log(i,i1,j1,k1)
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
!D1 Calculate  gradient of ln(porosity) assuming an exponential 
!D1 variation in the porosity between the nodes, ie, d(ln(phi))=constant
!D1 needed for evaluating the GRAD(phi)/(phi) term in random walk 
!D1 models.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 
!D2 Initial implementation: Sep 2002, Programmer: S. Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/porosity_gradient_log.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:40   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:12:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
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
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!D4 D.Grad(phi)/(phi) term from Tompson and Gelhar, WRR, Oct 90, eq. 8
!D4
!**********************************************************************

      use comai
      use combi
      use comdi
      use comsptr
      implicit none

      integer i,j,ix,iy,iz,i1,j1,k1

c.........for debugging s kelkar 7.31.01
      integer jj,nfence,nodesf(10000)
c............................

      real*8 dx,dy,dz
      real*8 porx1,porx2,pory1,pory2,porz1,porz2,pori
      real*8 length_factor, fac1
      parameter(length_factor=1.)
      
c zvd 13-May-08 change ps(*) to ps_trac(*)      
      porx1=ps_trac(i)
      pory1=ps_trac(i)
      porz1=ps_trac(i)

c recall neighbouring node numbers and
c account for boundary nodes, s kelkar, 12/19/02

      ix=irray(i,i1)
      if(ix.le.0) then
         ix=irray(i,-i1)
         if (ix.le.0) then
         end if
      end if
      iy =irray(i ,2*j1)
      if(iy.le.0) then
         iy = irray(i,-2*j1)
         if(iy.le.0) then
         end if
      end if
      iz =irray(i ,3*k1)
      if(iz.le.0) then
         iz = i
         if(iz.le.0) then
         end if
      end if
c
c     form distances between the neighbours of node i
c     recall porosities and form derivatives
      if (ix .gt. 0) then
         dx=(cord(i,1)-cord(ix,1))/length_factor
         porx2=ps_trac(ix)
         fac1=porx1/porx2
         if(fac1.gt.0.) then
            dpordx(1)=dlog(fac1)/dx
         else
            write(ierr,*)'STOP.ERROR in porosity_gradient.'
            write(ierr,*)'porx1/porx2 le 0.'
            write(ierr,*)'i,ix,iy,iz=',i,ix,iy,iz
            write(ierr,*)'porx1=', porx1, 'porx2=', porx2
            stop 
         endif
      end if
      if (iy .gt. 0) then
         dy=(cord(i,2)-cord(iy,2))/length_factor
         pory2=ps_trac(iy)
         fac1=pory1/pory2
         if(fac1.gt.0.) then
            dpordx(2)=dlog(fac1)/dy
         else
            write(ierr,*)'STOP.ERROR in porosity_gradient.'
            write(ierr,*)'pory1/pory2 le 0.'
            write(ierr,*)'i,ix,iy,iz=',i,ix,iy,iz
            write(ierr,*)'pory1=', pory1, 'pory2=', pory2
            stop
         endif
      end if
      if (iz .gt. 0) then
         dz=(cord(i,3)-cord(iz,3))/length_factor
         porz2=ps_trac(iz)
         fac1=porz1/porz2
         if(fac1.gt.0.) then
            dpordx(3)=dlog(fac1)/dz
         else
            write(ierr,*)'STOP.ERROR in porosity_gradient.'
            write(ierr,*)'porz1/porz2 le 0.'
            write(ierr,*)'i,ix,iy,iz=',i,ix,iy,iz
            write(ierr,*)'porz1=', porz1, 'porz2=', porz2
            stop
         endif
      end if
      
      
      return
      
      end


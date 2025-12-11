      subroutine porosity_gradient_omr(inp1)
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 On an OMR grid, using a least squares to interpolate,
!D1 calculate  gradient of ln(porosity) assuming an exponential 
!D1 variation in the porosity between the nodes, ie, d(ln(phi))=constant
!D1 needed for evaluating the GRAD(phi)/(phi) term in random walk 
!D1 models.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version
!D2 
!D2 Initial implementation: April 2005, Programmer: S. Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/porosity_gradient_log.f_a  $
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

      integer inp1, i,j,kb,ii,jj,i1,i2

      real*8 dx(3),dxw(3), cmat(3,3),rhs(3),cmat_inv(3,3)
      real*8 porkb,pori,pord,distij,distij2
      real*8 length_factor, fac1
      parameter(length_factor=1.)

      cmat=0.
      rhs=0.
      
      pori=ps_trac(inp1)

      i1=nelm(inp1)+1
      i2=nelm(inp1+1)
c loop over 1-st neighbours
      do j=i1,i2
         kb=nelm(j)
         if(kb.ne.inp1) then
            porkb=ps_trac(kb)
            pord=porkb-pori
c     displacement and distance between inp1 and kb
            distij2=0.
            do i=1,3
               dx(i)=(cord(kb,i)-cord(inp1,i))
               distij2=distij2+dx(i)*dx(i)
            enddo
            distij=sqrt(distij2)
c distance weighting
            dxw=dx/distij2
c calculate the coefficient matrix and rhs of the lsq equn
            do ii=1,3
               rhs(ii)=rhs(ii)+pord*dxw(ii)
               do jj=ii,3
                  cmat(ii,jj)=cmat(ii,jj)+dx(ii)*dxw(jj)
               enddo
            enddo
         endif
      enddo
c cmat is symmetric
      do ii=2,3
         do jj=1,ii-1
            cmat(ii,jj)=cmat(jj,ii)
         enddo
      enddo
c inverse of cmat
      call matinv3(cmat,cmat_inv)
c solve the equation
      do ii=1,3
         dpordx_omr(inp1,ii)=0.
         do jj=1,3
            dpordx_omr(inp1,ii)=dpordx_omr(inp1,ii)+
     $           cmat_inv(ii,jj)*rhs(jj)
         enddo
      enddo
      
      return
      
      end subroutine porosity_gradient_omr

c...............................................................

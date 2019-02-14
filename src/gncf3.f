      subroutine gncf3(nrq,nele,nga,neu,nsl,icsh,neumax,aj)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  To generate 3-D finite element coefficients.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gncf3.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:50   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Wed May 01 14:32:12 1996   gaz
CD2 correction for mixed elements
CD2 
CD2    Rev 1.6   Fri Apr 26 15:31:26 1996   gaz
CD2 error corrections for 2-d triangles in 3-d space
CD2 
CD2    Rev 1.5   Thu Jan 18 09:37:00 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.4   Wed Jan 10 12:46:24 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   Tue Jan 09 14:06:20 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.2   03/18/94 15:57:42   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   03/08/94 15:46:02   llt
CD2 moved the allocation and deallocation of array aj to the top and 
CD2 bottom of gencof (instead of inside gncf2 and gncf3).
CD2 
CD2    Rev 1.0   01/20/94 10:24:36   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  This subroutine calculates the element coefficient nrq for 
CD4  element nele in group neu with integration point nga.  If 
CD4  nrq is equal to 1 then the shape function and jacobian terms
CD4  are also calculated.
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5
CD5   aj              REAL*8   I
CD5   icsh            INT      I    Flag to indicate shape function and
CD5                                   jacobian terms should be calculated
CD5   nele            INT      I    Element
CD5   neu             INT      I    Element group
CD5   neumax          INT      I
CD5   nga             INT      I    Integration point
CD5   nrq             INT      O    Element coefficient position
CD5   nsl             INT
CD5
CD5 Interface Tables
CD5
CD5   None
CD5
CD5 Files
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6
CD6                            COMMON
CD6   Identifier      Type     Block  Description
CD6
CD6   bcoef           REAL*8   fff    Scratch storage array for coefficient
CD6                                     calculations
CD6   cord            REAL*8   fbs    Contains the coordinates of each node
CD6   dp              REAL*8   fbs    Contains weights for integration points
CD6                                     (prisms, triangles)
CD6   dr              REAL*8   fbs    Contains weights for integration
CD6                                     points (bricks, rectangles)
CD6   eta             REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   exci            REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   istrs           INT      faai   Parameter indicating if the stress
CD6                                     solution is enabled
CD6   ivf             INT      faai   ?
CD6   ldn             INT      param  Maximum array space needed for jacobian
CD6                                     array matrix
CD6   maxor           INT      param  Length of gmres working arrays
CD6   n0              INT      param  Maximum number of nodes allowed
CD6   n3              INT      param  3 * n0, storage parameter
CD6   n6              INT      param  Array storage for ice solution
CD6   n7              INT      param  Array storage for tracer solution
CD6   n8              INT      param  Array storage for variable  porosity
CD6                                     solution
CD6   nbd             INT      param  180 * n0 maximum array space for
CD6                                     incomplete lu decomposition matrix
CD6   nei             INT      faai   Total number of elements in the problem
CD6   nelm            INT      fbb    Initially information about nodes in each
CD6                                     element, later nodal connectivity
CD6                                     information
CD6   nelmd           INT      param  Maximum array space for element
CD6                                     connectivity array and (later) the nodal
CD6                                     connectivity array
CD6   nelucm          INT      param  nbd / 64
CD6   nnop            INT      param  Maximum array space for lu decomposition
CD6   nq              INT      param  Maximum array space for each finite
CD6                                     element coefficient array associated
CD6                                     with the stress solution
CD6   nr              INT      param  Maximum space allowed for each finite
CD6                                     element coefficient array
CD6   ns              INT      faai   Number of nodes per element
CD6   si              REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   w               REAL*8   fbs    Finite element shape functions
CD6   wr              REAL*8   fbs    Finite element shape functions
CD6                                     (rectangles)
CD6   wx              REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to x
CD6   wxp             REAL*8   fbs    Derivative of shape functions with
CD6                                      respect to x (prisms)
CD6   wxr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to x (rectangles)
CD6   wy              REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y
CD6   wyp             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y (prisms)
CD6   wyr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y (rectangles)
CD6   wz              REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to z
CD6   wzp             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to z (prisms)
CD6   wzr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to z (rectangles)
CD6   xt              REAL*8   fbs    Parameters needed in element calculations
CD6   yt              REAL*8   fbs    Parameters needed in element calculations
CD6   zt              REAL*8   fbs    Parameters needed in element calculations
CD6
CD6 Global Subprograms
CD6
CD6   Identifier      Type     Description
CD6
CD6   lubksub0        N/A      Solves the set of n linear equations A*x = b.
CD6   ludcmp0         N/A      Given an n by n matrix A, with physical
CD6                              dimension np, this routine replaces it
CD6                              by the LU decompostion of a row wise
CD6                              permutation of itself.
CD6   pebi3           N/A      Calculate the internodal volume (volume between
CD6                              nodes) using perpendicular bisectors.
CD6   shap3p          N/A      Evaluate 3-D prism element shape functions
CD6                              at quadrature points.
CD6   shap3r          N/A      Evaluate 3-D brick element shape functions
CD6                              at quadrature points.
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7   None
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   a1              REAL*8
CD7   a2              REAL*8
CD7   a3              REAL*8
CD7   a4              REAL*8
CD7   a5              REAL*8
CD7   a6              REAL*8
CD7   a7              REAL*8
CD7   a11             REAL*8
CD7   a12             REAL*8
CD7   a13             REAL*8
CD7   a21             REAL*8
CD7   a22             REAL*8
CD7   a23             REAL*8
CD7   a31             REAL*8
CD7   a32             REAL*8
CD7   a33             REAL*8
CD7   b1              REAL*8
CD7   b2              REAL*8
CD7   b3              REAL*8
CD7   b4              REAL*8
CD7   b5              REAL*8
CD7   b6              REAL*8
CD7   b7              REAL*8
CD7   c1              REAL*8
CD7   c2              REAL*8
CD7   c3              REAL*8
CD7   c4              REAL*8
CD7   c5              REAL*8
CD7   c6              REAL*8
CD7   c7              REAL*8
CD7   cord1           REAL*8
CD7   cord2           REAL*8   
CD7   cord3           REAL*8  
CD7   cpb             REAL*8
CD7   detja           REAL*8
CD7   dnga            REAL*8
CD7   dterm           REAL*8
CD7   ij              INT
CD7   iq              INT
CD7   ir              INT
CD7   jz              INT
CD7   k               INT
CD7   kb              INT
CD7   kjz             INT
CD7   knum            INT
CD7   sa11            REAL*8
CD7   sa12            REAL*8
CD7   sa13            REAL*8
CD7   sa21            REAL*8
CD7   sa22            REAL*8
CD7   sa23            REAL*8
CD7   sa31            REAL*8
CD7   sa32            REAL*8
CD7   sa33            REAL*8
CD7   vold            REAL*8
CD7   volt            REAL*8
CD7   vpebi           REAL*8
CD7   wxnga           REAL*8
CD7   wynga           REAL*8
CD7   wznga           REAL*8
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************

      use combi
      use comci
      use comdi
      use comei
      use comdti
      use comai
      implicit none
      
      integer icsh, nele, neu, neumax, nga, nrq, nsl, i, j
      integer ij, indx(4), iq, ir, jz, k, kb, kjz, knum
      integer ich_pebi_sv
      real*8  aj(neumax,*)
      real*8  cpb(4,4), vpebi(4)
      real*8  a1, a2, a3, a4, a5, a6, a7
      real*8  a11, a12, a13, a21, a22, a23, a31, a32, a33, area, area2
      real*8  b1, b2, b3, b4, b5, b6, b7
      real*8  c1, c2, c3, c4, c5, c6, c7
      real*8  cord1, cord2, cord3, detja, dnga, dterm
      real*8  sa11, sa12, sa13, sa21, sa22, sa23, sa31, sa32, sa33
      real*8  vold, vold4, volt, wxnga, wynga, wznga
      real*8 vol_tol
      parameter (vol_tol = -1.d-12)

      
c     if icsh=1 calculate jacobian information
c     
      if ( icsh.eq.1 )  then
c  
         if ( nsl.eq.8 )  then
c     procedure for 8 node brick elements
            call shap3r(nga)
            do jz=1,nsl
               kb=nelm((nele-1)*ns+jz)
               xt(jz)=cord(kb,1)
               yt(jz)=cord(kb,2)
               zt(jz)=cord(kb,3)
            enddo
            
            a5=0.125*(-xt(1)-xt(2)+xt(3)+xt(4)-xt(5)-xt(6)+xt(7)+xt(8))
            a2=0.125*(xt(1)-xt(2)+xt(3)-xt(4)+xt(5)-xt(6)+xt(7)-xt(8))
            a6=-0.125*(-xt(1)-xt(2)+xt(3)+xt(4)+xt(5)+xt(6)-xt(7)-xt(8))
            a4=-0.125*(xt(1)-xt(2)+xt(3)-xt(4)-xt(5)+xt(6)-xt(7)+xt(8))
            a1=0.125*(-xt(1)+xt(2)+xt(3)-xt(4)-xt(5)+xt(6)+xt(7)-xt(8))
            a3=-0.125*(-xt(1)+xt(2)+xt(3)-xt(4)+xt(5)-xt(6)-xt(7)+xt(8))
            a7=-0.125*(xt(1)+xt(2)+xt(3)+xt(4)-xt(5)-xt(6)-xt(7)-xt(8))
            b5=0.125*(-yt(1)-yt(2)+yt(3)+yt(4)-yt(5)-yt(6)+yt(7)+yt(8))
            b2=0.125*(yt(1)-yt(2)+yt(3)-yt(4)+yt(5)-yt(6)+yt(7)-yt(8))
            b6=-0.125*(-yt(1)-yt(2)+yt(3)+yt(4)+yt(5)+yt(6)-yt(7)-yt(8))
            b4=-0.125*(yt(1)-yt(2)+yt(3)-yt(4)-yt(5)+yt(6)-yt(7)+yt(8))
            b1=0.125*(-yt(1)+yt(2)+yt(3)-yt(4)-yt(5)+yt(6)+yt(7)-yt(8))
            b3=-0.125*(-yt(1)+yt(2)+yt(3)-yt(4)+yt(5)-yt(6)-yt(7)+yt(8))
            b7=-0.125*(yt(1)+yt(2)+yt(3)+yt(4)-yt(5)-yt(6)-yt(7)-yt(8))
            c5=0.125*(-zt(1)-zt(2)+zt(3)+zt(4)-zt(5)-zt(6)+zt(7)+zt(8))
            c2=0.125*(zt(1)-zt(2)+zt(3)-zt(4)+zt(5)-zt(6)+zt(7)-zt(8))
            c6=-0.125*(-zt(1)-zt(2)+zt(3)+zt(4)+zt(5)+zt(6)-zt(7)-zt(8))
            c4=-0.125*(zt(1)-zt(2)+zt(3)-zt(4)-zt(5)+zt(6)-zt(7)+zt(8))
            c1=0.125*(-zt(1)+zt(2)+zt(3)-zt(4)-zt(5)+zt(6)+zt(7)-zt(8))
            c3=-0.125*(-zt(1)+zt(2)+zt(3)-zt(4)+zt(5)-zt(6)-zt(7)+zt(8))
            c7=-0.125*(zt(1)+zt(2)+zt(3)+zt(4)-zt(5)-zt(6)-zt(7)-zt(8))
            
            sa11=a1+a2*eta(nga)+a3*exci(nga)+a4*eta(nga)*exci(nga)
            sa21=a5+a2*si(nga)+a6*exci(nga)+a4*si(nga)*exci(nga)
            sa31=a7+a3*si(nga)+a6*eta(nga)+a4*si(nga)*eta(nga)
            sa12=b1+b2*eta(nga)+b3*exci(nga)+b4*eta(nga)*exci(nga)
            sa22=b5+b2*si(nga)+b6*exci(nga)+b4*si(nga)*exci(nga)
            sa32=b7+b3*si(nga)+b6*eta(nga)+b4*si(nga)*eta(nga)
            sa13=c1+c2*eta(nga)+c3*exci(nga)+c4*eta(nga)*exci(nga)
            sa23=c5+c2*si(nga)+c6*exci(nga)+c4*si(nga)*exci(nga)
            sa33=c7+c3*si(nga)+c6*eta(nga)+c4*si(nga)*eta(nga)
       
            sa11=  -sa11
            sa21=  -sa21
            sa31=  -sa31
            sa12=  -sa12
            sa22=  -sa22
            sa32=  -sa32
            sa13=  -sa13
            sa23=  -sa23
            sa33=  -sa33

         endif
         if ( nsl.eq.6 )  then
c     procedure for 6 node prism
            call shap3p(nga)
c     
c     define geometric coefficients
c     
            sa11=0.
            sa12=0.
            sa22=0.
            sa21=0.
            sa13=0.
            sa23=0.
            sa31=0.
            sa32=0.
            sa33=0.
            do jz=1,nsl
               kb=nelm((nele-1)*ns+jz)
               cord1=cord(kb,1)
               cord2=cord(kb,2)
               cord3=cord(kb,3)
               wxnga=-wxp(nga,jz)
               wynga=-wyp(nga,jz)
               wznga=-wzp(nga,jz)
               sa11=sa11+wxnga*cord1
               sa12=sa12+wxnga*cord2
               sa21=sa21+wynga*cord1
               sa22=sa22+wynga*cord2
               sa13=sa13+wxnga*cord3
               sa23=sa23+wynga*cord3
               sa31=sa31+wznga*cord1
               sa32=sa32+wznga*cord2
               sa33=sa33+wznga*cord3
            enddo
            
         endif
      endif
      if ( nsl.eq.3.and.nga.le.1 )  then
c     2-d triangles in 3-d space
         if(nrq.eq.1) then
            do jz=1,3
               kb=nelm((nele-1)*ns+jz)
               xt(jz)=cord(kb,1)
               yt(jz)=cord(kb,2)
               zt(jz)=cord(kb,3)
            enddo
            call area2d_tri(1,xt(1),yt(1),zt(1),
     &           xt(2),yt(2),zt(2),xt(3),yt(3),zt(3),                 
     &           volt,vpebi(1),vpebi(2),vpebi(3),
     &           cpb(1,1),cpb(2,1),cpb(3,1))

            do jz=1,3
               kb=nelm((nele-1)*ns+jz)
               vpebi(jz)=vpebi(jz)*thic(kb)
               cpb(jz,1)=cpb(jz,1)*thic(kb)
            enddo
            bcoef(neu,1)=vpebi(1)
            bcoef(neu,2)=0.0
            bcoef(neu,3)=0.0
            bcoef(neu,4)=0.0
            bcoef(neu,5)=vpebi(2)
            bcoef(neu,6)=0.0
            bcoef(neu,7)=0.0
            bcoef(neu,8)=0.0
            bcoef(neu,9)=vpebi(3)
c     save cpb info
            bcoef(neu,9+1)= 0.0         
            bcoef(neu,9+2)= -cpb(1,1)
            bcoef(neu,9+3)= -cpb(3,1)
            bcoef(neu,9+4)= -cpb(1,1)         
            bcoef(neu,9+5)= 0.0              
            bcoef(neu,9+6)= -cpb(2,1)        
            bcoef(neu,9+7)= -cpb(3,1)        
            bcoef(neu,9+8)= -cpb(2,1)        
            bcoef(neu,9+9)= 0.0              
            go to 999
         else if(nrq.eq.2) then
c     calculate x area/dx
            bcoef(neu,1)=0.0
            bcoef(neu,2)=bcoef(neu,9+2)
            bcoef(neu,3)=bcoef(neu,9+3)
            bcoef(neu,4)=bcoef(neu,9+4)
            bcoef(neu,5)=0.0
            bcoef(neu,6)=bcoef(neu,9+6)
            bcoef(neu,7)=bcoef(neu,9+7)
            bcoef(neu,8)=bcoef(neu,9+8)
            bcoef(neu,9)=0.0
            goto 999
         else if(nrq.ge.3) then
c     calculate y area/dy
            bcoef(neu,1)=0.0
            bcoef(neu,2)=0.0
            bcoef(neu,3)=0.0
            bcoef(neu,4)=0.0
            bcoef(neu,5)=0.0
            bcoef(neu,6)=0.0
            bcoef(neu,7)=0.0
            bcoef(neu,8)=0.0
            bcoef(neu,9)=0.0
            goto 999
         endif
      else if(nsl.eq.3.and.nga.gt.1) then
         go to 999
      endif
      if ( nga.le.1 )  then
         if ( nsl.eq.4 )  then
            if(ivf.eq.0.or.nrq.gt.4) then
c     
c     form shape functions by matrix inversion
c     wx will contain the matrix
c     
               if(nrq.eq.1) then
                  do jz=1,4
                     kb=nelm((nele-1)*ns+jz)
                     wx(jz,1)=1.0
                     wx(jz,2)=cord(kb,1)
                     wx(jz,3)=cord(kb,2)
                     wx(jz,4)=cord(kb,3)
                  enddo
c     
c     solve for inverse and determinant
c     
                  call ludcmp0(wx,4,8,indx,dterm)
                  do iq=1,4
                     do ir=1,4
                        cpb(iq,ir)=0.0
                     enddo
                     cpb(iq,iq)=1.0
                  enddo
                  do ir=1,4
                     call lubksb0(wx,4,8,indx,cpb(1,ir))
                     dterm=dterm*wx(ir,ir)
                  enddo
                  volt=-dterm/6.0
c     vold4 is part apportioned to each node
                  vold4=volt/4.0
c     
c     partition volume equally
c     
                  bcoef(neu,1)=vold4
                  bcoef(neu,2)=0.0
                  bcoef(neu,3)=0.0
                  bcoef(neu,4)=0.0
                  bcoef(neu,5)=0.0
                  bcoef(neu,6)=vold4
                  bcoef(neu,7)=0.0
                  bcoef(neu,8)=0.0
                  bcoef(neu,9)=0.0
                  bcoef(neu,10)=0.0
                  bcoef(neu,11)=vold4
                  bcoef(neu,12)=0.0
                  bcoef(neu,13)=0.0
                  bcoef(neu,14)=0.0
                  bcoef(neu,15)=0.0
                  bcoef(neu,16)=vold4
c     save cpb info
                  bcoef(neu,16+1)= volt     
                  bcoef(neu,16+2)= volt
                  bcoef(neu,16+3)= volt
                  bcoef(neu,16+4)= volt
                  bcoef(neu,16+5)= cpb(2,1)         
                  bcoef(neu,16+6)= cpb(2,2)         
                  bcoef(neu,16+7)= cpb(2,3)         
                  bcoef(neu,16+8)= cpb(2,4)         
                  bcoef(neu,16+9)= cpb(3,1)         
                  bcoef(neu,16+10)= cpb(3,2)         
                  bcoef(neu,16+11)= cpb(3,3)         
                  bcoef(neu,16+12)= cpb(3,4)         
                  bcoef(neu,16+13)= cpb(4,1)         
                  bcoef(neu,16+14)= cpb(4,2)         
                  bcoef(neu,16+15)= cpb(4,3)         
                  bcoef(neu,16+16)= cpb(4,4)         
               else if(nrq.eq.2) then
c     calculate x area/dx
                  bcoef(neu,1)= 0.0
                  bcoef(neu,2)= bcoef(neu,16+5)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+5)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+5)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,2)
                  bcoef(neu,6)= 0.0
                  bcoef(neu,7)= bcoef(neu,16+6)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+6)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,3)
                  bcoef(neu,10)= bcoef(neu,7)
                  bcoef(neu,11)= 0.0
                  bcoef(neu,12)= bcoef(neu,16+7)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,4)
                  bcoef(neu,14)= bcoef(neu,8)
                  bcoef(neu,15)= bcoef(neu,12)
                  bcoef(neu,16)= 0.0
               else if(nrq.eq.3) then
c     calculate y area/dy
                  bcoef(neu,1)= 0.0
                  bcoef(neu,2)= bcoef(neu,16+9)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+9)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+9)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,2)
                  bcoef(neu,6)= 0.0
                  bcoef(neu,7)= bcoef(neu,16+10)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+10)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,3)
                  bcoef(neu,10)= bcoef(neu,7)
                  bcoef(neu,11)= 0.0
                  bcoef(neu,12)= bcoef(neu,16+11)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,4)
                  bcoef(neu,14)= bcoef(neu,8)
                  bcoef(neu,15)= bcoef(neu,12)
                  bcoef(neu,16)= 0.0
               else if(nrq.eq.4) then
c     calculate z area/dz
                  bcoef(neu,1)= 0.0
                  bcoef(neu,2)= bcoef(neu,16+13)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+13)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+13)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,2)
                  bcoef(neu,6)= 0.0
                  bcoef(neu,7)= bcoef(neu,16+14)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+14)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,3)
                  bcoef(neu,10)= bcoef(neu,7)
                  bcoef(neu,11)= 0.0
                  bcoef(neu,12)= bcoef(neu,16+15)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,4)
                  bcoef(neu,14)= bcoef(neu,8)
                  bcoef(neu,15)= bcoef(neu,12)
                  bcoef(neu,16)= 0.0
               else if(nrq.eq.5) then
c     ni*dnj/dx
                  bcoef(neu,1)=0.25*bcoef(neu,16+5)*bcoef(neu,16+1)
                  bcoef(neu,2)=0.25*bcoef(neu,16+6)*bcoef(neu,16+1)
                  bcoef(neu,3)=0.25*bcoef(neu,16+7)*bcoef(neu,16+1)
                  bcoef(neu,4)=0.25*bcoef(neu,16+8)*bcoef(neu,16+1)
                  bcoef(neu,5)=bcoef(neu,1)
                  bcoef(neu,6)=bcoef(neu,2)
                  bcoef(neu,7)=bcoef(neu,3)
                  bcoef(neu,8)=bcoef(neu,4)
                  bcoef(neu,9)=bcoef(neu,1)
                  bcoef(neu,10)=bcoef(neu,2)
                  bcoef(neu,11)=bcoef(neu,3)
                  bcoef(neu,12)=bcoef(neu,4)
                  bcoef(neu,13)=bcoef(neu,1)
                  bcoef(neu,14)=bcoef(neu,2)
                  bcoef(neu,15)=bcoef(neu,3)
                  bcoef(neu,16)=bcoef(neu,4)
               else if(nrq.eq.6) then
c     ni*dnj/dy
                  bcoef(neu,1)=0.25*bcoef(neu,16+9)*bcoef(neu,16+1)
                  bcoef(neu,2)=0.25*bcoef(neu,16+10)*bcoef(neu,16+1)
                  bcoef(neu,3)=0.25*bcoef(neu,16+11)*bcoef(neu,16+1)
                  bcoef(neu,4)=0.25*bcoef(neu,16+12)*bcoef(neu,16+1)
                  bcoef(neu,5)=bcoef(neu,1)
                  bcoef(neu,6)=bcoef(neu,2)
                  bcoef(neu,7)=bcoef(neu,3)
                  bcoef(neu,8)=bcoef(neu,4)
                  bcoef(neu,9)=bcoef(neu,1)
                  bcoef(neu,10)=bcoef(neu,2)
                  bcoef(neu,11)=bcoef(neu,3)
                  bcoef(neu,12)=bcoef(neu,4)
                  bcoef(neu,13)=bcoef(neu,1)
                  bcoef(neu,14)=bcoef(neu,2)
                  bcoef(neu,15)=bcoef(neu,3)
                  bcoef(neu,16)=bcoef(neu,4)
               else if(nrq.eq.7) then
c     ni*dnj/dz
                  bcoef(neu,1)=0.25*bcoef(neu,16+13)*bcoef(neu,16+1)
                  bcoef(neu,2)=0.25*bcoef(neu,16+14)*bcoef(neu,16+1)
                  bcoef(neu,3)=0.25*bcoef(neu,16+15)*bcoef(neu,16+1)
                  bcoef(neu,4)=0.25*bcoef(neu,16+16)*bcoef(neu,16+1)
                  bcoef(neu,5)=bcoef(neu,1)
                  bcoef(neu,6)=bcoef(neu,2)
                  bcoef(neu,7)=bcoef(neu,3)
                  bcoef(neu,8)=bcoef(neu,4)
                  bcoef(neu,9)=bcoef(neu,1)
                  bcoef(neu,10)=bcoef(neu,2)
                  bcoef(neu,11)=bcoef(neu,3)
                  bcoef(neu,12)=bcoef(neu,4)
                  bcoef(neu,13)=bcoef(neu,1)
                  bcoef(neu,14)=bcoef(neu,2)
                  bcoef(neu,15)=bcoef(neu,3)
                  bcoef(neu,16)=bcoef(neu,4)
c     code for stress
               else if(nrq.eq.11.and.istrs.ne.0) then
c     dni/dx*dnj/dy
                  bcoef(neu,1)= bcoef(neu,16+5)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,2)= bcoef(neu,16+5)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+5)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+5)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,16+6)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,6)= bcoef(neu,16+6)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,7)= bcoef(neu,16+6)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+6)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,16+7)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,10)= bcoef(neu,16+7)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,11)= bcoef(neu,16+7)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,12)= bcoef(neu,16+7)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,16+8)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,14)= bcoef(neu,16+8)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,15)= bcoef(neu,16+8)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,16)= bcoef(neu,16+8)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
               else if(nrq.eq.12.and.istrs.ne.0) then
c     dni/dy*dnj/dx
                  bcoef(neu,1)= bcoef(neu,16+9)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,2)= bcoef(neu,16+9)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+9)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+9)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,16+10)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,6)= bcoef(neu,16+10)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,7)= bcoef(neu,16+10)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+10)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,16+11)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,10)= bcoef(neu,16+11)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,11)= bcoef(neu,16+11)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,12)= bcoef(neu,16+11)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,16+12)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,14)= bcoef(neu,16+12)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,15)= bcoef(neu,16+12)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,16)= bcoef(neu,16+12)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
               else if(nrq.eq.13.and.istrs.ne.0) then
c     dni/dx*dnj/dz
                  bcoef(neu,1)= bcoef(neu,16+5)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,2)= bcoef(neu,16+5)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+5)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+5)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,16+6)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,6)= bcoef(neu,16+6)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,7)= bcoef(neu,16+6)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+6)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,16+7)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,10)= bcoef(neu,16+7)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,11)= bcoef(neu,16+7)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,12)= bcoef(neu,16+7)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,16+8)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,14)= bcoef(neu,16+8)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,15)= bcoef(neu,16+8)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,16)= bcoef(neu,16+8)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
               else if(nrq.eq.14.and.istrs.ne.0) then
c     dni/dz*dnj/dx
                  bcoef(neu,1)= bcoef(neu,16+13)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,2)= bcoef(neu,16+13)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+13)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+13)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,16+14)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,6)= bcoef(neu,16+14)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,7)= bcoef(neu,16+14)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+14)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,16+15)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,10)= bcoef(neu,16+15)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,11)= bcoef(neu,16+15)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,12)= bcoef(neu,16+15)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,16+16)*bcoef(neu,16+5)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,14)= bcoef(neu,16+16)*bcoef(neu,16+6)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,15)= bcoef(neu,16+16)*bcoef(neu,16+7)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,16)= bcoef(neu,16+16)*bcoef(neu,16+8)*
     &                 bcoef(neu,16+1)
               else if(nrq.eq.15.and.istrs.ne.0) then
c     dni/dy*dnj/dz
                  bcoef(neu,1)= bcoef(neu,16+9)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,2)= bcoef(neu,16+9)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+9)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+9)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,16+10)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,6)= bcoef(neu,16+10)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,7)= bcoef(neu,16+10)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+10)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,16+11)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,10)= bcoef(neu,16+11)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,11)= bcoef(neu,16+11)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,12)= bcoef(neu,16+11)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,16+12)*bcoef(neu,16+13)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,14)= bcoef(neu,16+12)*bcoef(neu,16+14)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,15)= bcoef(neu,16+12)*bcoef(neu,16+15)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,16)= bcoef(neu,16+12)*bcoef(neu,16+16)*
     &                 bcoef(neu,16+1)
               else if(nrq.eq.16.and.istrs.ne.0) then
c     dni/dz*dnj/dy
                  bcoef(neu,1)= bcoef(neu,16+13)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,2)= bcoef(neu,16+13)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,3)= bcoef(neu,16+13)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,4)= bcoef(neu,16+13)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,5)= bcoef(neu,16+14)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,6)= bcoef(neu,16+14)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,7)= bcoef(neu,16+14)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,8)= bcoef(neu,16+14)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,9)= bcoef(neu,16+15)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,10)= bcoef(neu,16+15)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,11)= bcoef(neu,16+15)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,12)= bcoef(neu,16+15)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,13)= bcoef(neu,16+16)*bcoef(neu,16+9)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,14)= bcoef(neu,16+16)*bcoef(neu,16+10)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,15)= bcoef(neu,16+16)*bcoef(neu,16+11)*
     &                 bcoef(neu,16+1)
                  bcoef(neu,16)= bcoef(neu,16+16)*bcoef(neu,16+12)*
     &                 bcoef(neu,16+1)
               endif
               go to 999
            endif
         endif
      endif
      if(ivf.eq.0.or.nrq.gt.4) then
c     
c     save jacobian information
c     
         a11=sa22*sa33-sa23*sa32
         a21=-sa21*sa33+sa23*sa31
         a31=sa21*sa32-sa22*sa31
         a12=-sa12*sa33+sa13*sa32
         a22=sa11*sa33-sa13*sa31
         a32=-sa11*sa32+sa12*sa31
         a13=sa12*sa23-sa13*sa22
         a23=-sa11*sa23+sa13*sa21
         a33=sa11*sa22-sa12*sa21
         detja=sa11*sa22*sa33+sa12*sa23*sa31+sa13*sa21*sa32
     *        -sa13*sa22*sa31-sa23*sa32*sa11-sa33*sa21*sa12
         aj(neu,1)=-a11
         aj(neu,2)=-a12
         aj(neu,3)=-a13
         aj(neu,4)=-a21
         aj(neu,5)=-a22
         aj(neu,6)=-a23
         aj(neu,7)=-a31
         aj(neu,8)=-a32
         aj(neu,9)=-a33
         aj(neu,10)=-detja
      else
c     
c     pebi calcs
c     
c     pebi tetrahedrals
         do jz=1,nsl
            kb=nelm((nele-1)*ns+jz)
            xt(jz)=cord(kb,1)
            yt(jz)=cord(kb,2)
            zt(jz)=cord(kb,3)
         enddo
         if(nrq.eq.1) then
            ich_pebi_sv = ich_pebi
            if(ireord.eq.10) then
               call area_vol_tet(xt,yt,zt,cpb,vpebi)
c check for all negative volumes
               if(vpebi(1).lt.vol_tol.and.vpebi(2).lt.vol_tol.and.
     &            vpebi(3).lt.vol_tol.and.vpebi(4).lt.vol_tol) then
                  ich_pebi_sv = ich_pebi 
                  ich_pebi = 1
c write element number to err file
              write(ierr,*) 'element ',nele,' has all neg node vols, ',
     &          'negative rule applied'        
               endif
            else 
               call pebi3(xt,yt,zt,cpb,vpebi)
c check for all negative volumes
               if(vpebi(1).lt.vol_tol.and.vpebi(2).lt.vol_tol.and.
     &            vpebi(3).lt.vol_tol.and.vpebi(4).lt.vol_tol) then
                  ich_pebi_sv = ich_pebi 
                  ich_pebi = 1  
c write element number to err file
               write(ierr,*) 'element ',nele,' has all neg node vols, ',
     &          'negative rule applied'        
               endif                  
            endif
            if(ich_pebi.ne.0) then
                do i = 1, 4
                  vpebi(i) = - vpebi(i)
                  do j = 1, 4  
                   cpb(i,j) = - cpb(i,j)
                  enddo
                enddo
            endif
            ich_pebi = ich_pebi_sv
c     calculate volumes
            bcoef(neu,1)=vpebi(1)
            bcoef(neu,2)=0.0
            bcoef(neu,3)=0.0
            bcoef(neu,4)=0.0
            bcoef(neu,5)=0.0
            bcoef(neu,6)=vpebi(2)
            bcoef(neu,7)=0.0
            bcoef(neu,8)=0.0
            bcoef(neu,9)=0.0
            bcoef(neu,10)=0.0
            bcoef(neu,11)=vpebi(3)
            bcoef(neu,12)=0.0
            bcoef(neu,13)=0.0
            bcoef(neu,14)=0.0
            bcoef(neu,15)=0.0
            bcoef(neu,16)=vpebi(4)
            bcoef(neu,16+1)= 0.0      
            bcoef(neu,16+2)= -cpb(1,2)
            bcoef(neu,16+3)= -cpb(1,3)
            bcoef(neu,16+4)= -cpb(1,4)
            bcoef(neu,16+5)= -cpb(2,1)
            bcoef(neu,16+6)= 0.0 
            bcoef(neu,16+7)= -cpb(2,3)
            bcoef(neu,16+8)= -cpb(2,4)
            bcoef(neu,16+9)= -cpb(3,1)
            bcoef(neu,16+10)= -cpb(3,2)
            bcoef(neu,16+11)= 0.0 
            bcoef(neu,16+12)= -cpb(3,4)
            bcoef(neu,16+13)= -cpb(4,1)
            bcoef(neu,16+14)= -cpb(4,2)
            bcoef(neu,16+15)= -cpb(4,3)
            bcoef(neu,16+16)= 0.0 
            go to 999
         else if(nrq.eq.2) then
c     calculate x area/dx
            bcoef(neu,1)= bcoef(neu,16+1)
            bcoef(neu,2)= bcoef(neu,16+2)
            bcoef(neu,3)= bcoef(neu,16+3)
            bcoef(neu,4)= bcoef(neu,16+4)
            bcoef(neu,5)= bcoef(neu,16+5)
            bcoef(neu,6)= bcoef(neu,16+6)
            bcoef(neu,7)= bcoef(neu,16+7)
            bcoef(neu,8)= bcoef(neu,16+8)
            bcoef(neu,9)= bcoef(neu,16+9)
            bcoef(neu,10)= bcoef(neu,16+10)
            bcoef(neu,11)= bcoef(neu,16+11)
            bcoef(neu,12)= bcoef(neu,16+12)
            bcoef(neu,13)= bcoef(neu,16+13)
            bcoef(neu,14)= bcoef(neu,16+14)
            bcoef(neu,15)= bcoef(neu,16+15)
            bcoef(neu,16)= bcoef(neu,16+16)
            go to 999
         else if(nrq.eq.3.or.nrq.eq.4) then
c     calculate y area/dy
            bcoef(neu,1)=0.0
            bcoef(neu,2)=0.0
            bcoef(neu,3)=0.0
            bcoef(neu,4)=0.0
            bcoef(neu,5)=0.0
            bcoef(neu,6)=0.0
            bcoef(neu,7)=0.0
            bcoef(neu,8)=0.0
            bcoef(neu,9)=0.0
            bcoef(neu,10)=0.0
            bcoef(neu,11)=0.0
            bcoef(neu,12)=0.0
            bcoef(neu,13)=0.0
            bcoef(neu,14)=0.0
            bcoef(neu,15)=0.0
            bcoef(neu,16)=0.0
            go to 999
         endif
      endif
c     
c     load shape function for brick or prism
c     
      if ( nsl.eq.8 )  then
c     bricks
         do  kjz=1,8
            w(nga,kjz)=wr(nga,kjz)
            wx(nga,kjz)=wxr(nga,kjz)
            wy(nga,kjz)=wyr(nga,kjz)
            wz(nga,kjz)=wzr(nga,kjz)
         enddo
         dnga=dr(nga)
      elseif(nsl.eq.6) then
         do kjz=1,6
            w(nga,kjz)=wp(nga,kjz)
            wx(nga,kjz)=wxp(nga,kjz)
            wy(nga,kjz)=wyp(nga,kjz)
            wz(nga,kjz)=wzp(nga,kjz)
         enddo
         dnga=dp(nga)
      elseif(nsl.eq.4) then
         do kjz=1,4
            w(nga,kjz)=wr(1,kjz)
            wx(nga,kjz)=wxr(1,kjz)
            wy(nga,kjz)=wyr(1,kjz)
            wz(nga,kjz)=wzr(1,kjz)
         enddo
         dnga=dr(nga)
      endif
c     
c     assemble global derivatives
c     
      if ( nrq.eq.1 )  then
c     nk*njz
         do k=1,nsl
            knum=(k-1)*nsl
c            do jz=1,nsl
c               ij=knum+jz
c               bcoef(neu,ij)=bcoef(neu,ij)+w(nga,k)*w(nga,jz)*dnga*
c     &              aj(neu,10)
c            enddo
            ij = knum + k
            bcoef(neu,ij)=bcoef(neu,ij)+w(nga,k)*dnga*aj(neu,10)
         enddo
      endif
      if ( nrq.eq.2.or.nrq.eq.20 )  then
c     (d(nk/dx)*(dnjz/dx)
c     recall jacobian information
         a11=aj(neu,1)
         a12=aj(neu,2)
         a13=aj(neu,3)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     *              (a11*wx(nga,k)+a12*wy(nga,k)+a13*wz(nga,k))*
     5              (a11*wx(nga,jz)+a12*wy(nga,jz)+a13*wz(nga,jz))*
     &              dnga/detja
            enddo
         enddo
      endif
      if ( nrq.eq.3.or.nrq.eq.21 )  then
c     (d(nk/dy)*(dnjz/dy)
c     recall jacobian information
         a21=aj(neu,4)
         a22=aj(neu,5)
         a23=aj(neu,6)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     *              (a21*wx(nga,k)+a22*wy(nga,k)+a23*wz(nga,k))*
     5              (a21*wx(nga,jz)+a22*wy(nga,jz)+a23*wz(nga,jz))*
     &              dnga/detja
            enddo
         enddo
      endif
      if ( nrq.eq.4.or.nrq.eq.22 )  then
c     (d(nk/dz)*(dnjz/dz)
c     recall jacobian information
         a31=aj(neu,7)
         a32=aj(neu,8)
         a33=aj(neu,9)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     *              (a31*wx(nga,k)+a32*wy(nga,k)+a33*wz(nga,k))*
     5              (a31*wx(nga,jz)+a32*wy(nga,jz)+a33*wz(nga,jz))*
     &              dnga/detja
            enddo
         enddo
      endif

c     c for 3-d stress
      if ( istrs.ne.0 )  then
         if ( nrq.eq.11 )  then
c     (dnk/dx)*(dnjz/dy)
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            a13=aj(neu,3)
            a21=aj(neu,4)
            a22=aj(neu,5)
            a23=aj(neu,6)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 (a11*wx(nga,k)+a12*wy(nga,k)+a13*wz(nga,k))*
     5                 (a21*wx(nga,jz)+a22*wy(nga,jz)+a23*wz(nga,jz))*
     &                 dnga/detja
               enddo
            enddo
         endif
         if ( nrq.eq.12 )  then
c     (dnk/dy)*(dnjz/dx)
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            a13=aj(neu,3)
            a21=aj(neu,4)
            a22=aj(neu,5)
            a23=aj(neu,6)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 (a21*wx(nga,k)+a22*wy(nga,k)+a23*wz(nga,k))*
     5                 (a11*wx(nga,jz)+a12*wy(nga,jz)+a13*wz(nga,jz))*
     &                 dnga/detja
               enddo
            enddo
         endif
         if ( nrq.eq.13 )  then
c     (dnk/dx)*(dnjz/dz)
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            a13=aj(neu,3)
            a31=aj(neu,7)
            a32=aj(neu,8)
            a33=aj(neu,9)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 (a11*wx(nga,k)+a12*wy(nga,k)+a13*wz(nga,k))*
     5                 (a31*wx(nga,jz)+a32*wy(nga,jz)+a33*wz(nga,jz))*
     &                 dnga/detja
               enddo
            enddo
         endif
         if ( nrq.eq.14 )  then
c     (dnk/dz)*(dnjz/dx)
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            a13=aj(neu,3)
            a31=aj(neu,7)
            a32=aj(neu,8)
            a33=aj(neu,9)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 (a31*wx(nga,k)+a32*wy(nga,k)+a33*wz(nga,k))*
     5                 (a11*wx(nga,jz)+a12*wy(nga,jz)+a13*wz(nga,jz))*
     &                 dnga/detja
               enddo
            enddo
         endif
         if ( nrq.eq.15 )  then
c     (dnk/dy)*(dnjz/dz)
c     recall jacobian information
            a21=aj(neu,4)
            a22=aj(neu,5)
            a23=aj(neu,6)
            a31=aj(neu,7)
            a32=aj(neu,8)
            a33=aj(neu,9)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 (a21*wx(nga,k)+a22*wy(nga,k)+a23*wz(nga,k))*
     5                 (a31*wx(nga,jz)+a32*wy(nga,jz)+a33*wz(nga,jz))*
     &                 dnga/detja
               enddo
            enddo
         endif
         if ( nrq.eq.16 )  then
c     (dnk/dz)*(dnjz/dy)
c     recall jacobian information
            a21=aj(neu,4)
            a22=aj(neu,5)
            a23=aj(neu,6)
            a31=aj(neu,7)
            a32=aj(neu,8)
            a33=aj(neu,9)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 (a31*wx(nga,k)+a32*wy(nga,k)+a33*wz(nga,k))*
     5                 (a21*wx(nga,jz)+a22*wy(nga,jz)+a23*wz(nga,jz))*
     &                 dnga/detja
               enddo
            enddo
         endif
      if ( nrq.eq.17 )  then
c     nk*(dnjz/dx)
c     recall jacobian information
         a11=aj(neu,1)
         a12=aj(neu,2)
         a13=aj(neu,3)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     &              w(nga,k)*(a11*wx(nga,jz)+a12*wy(nga,jz)+a13*
     &              wz(nga,jz))*dnga
            enddo
         enddo
      endif
      if ( nrq.eq.18 )  then
c     nk*(dnjz/dy)
c     recall jacobian information
         a21=aj(neu,4)
         a22=aj(neu,5)
         a23=aj(neu,6)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     &              w(nga,k)*(a21*wx(nga,jz)+a22*wy(nga,jz)+a23*
     &              wz(nga,jz))*dnga
            enddo
         enddo
      endif
      if ( nrq.eq.19 )  then
c     nk*(dnjz/dz)
c     recall jacobian information
         a31=aj(neu,7)
         a32=aj(neu,8)
         a33=aj(neu,9)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     &              w(nga,k)*(a31*wx(nga,jz)+a32*wy(nga,jz)+a33*
     &              wz(nga,jz))*dnga
            enddo
         enddo
      endif
      if ( nrq.eq.23 )  then
c     nk*(dnjz/dx)
c     recall jacobian information
         a11=aj(neu,1)
         a12=aj(neu,2)
         a13=aj(neu,3)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     &              w(nga,jz)*(a11*wx(nga,k)+a12*wy(nga,k)+a13*
     &              wz(nga,k))*dnga
            enddo
         enddo
      endif
      if ( nrq.eq.24 )  then
c     nk*(dnjz/dy)
c     recall jacobian information 
         a21=aj(neu,4)
         a22=aj(neu,5)
         a23=aj(neu,6)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     &              w(nga,jz)*(a21*wx(nga,k)+a22*wy(nga,k)+a23*
     &              wz(nga,k))*dnga
            enddo
         enddo
      endif
      if ( nrq.eq.25 )  then
c     nk*(dnjz/dz)
c     recall jacobian information
         a31=aj(neu,7)
         a32=aj(neu,8)
         a33=aj(neu,9)
         detja=aj(neu,10)
         do k=1,nsl
            knum=(k-1)*nsl
            do jz=1,nsl
               ij=knum+jz
               bcoef(neu,ij)=bcoef(neu,ij)+
     &              w(nga,jz)*(a31*wx(nga,k)+a32*wy(nga,k)+a33*
     &              wz(nga,k))*dnga
            enddo
         enddo
      endif      
      endif
 999  continue
      
      return
      end



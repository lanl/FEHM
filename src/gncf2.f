      subroutine gncf2(nrq,nele,nga,neu,nsl,icsh,neumax,aj)
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
CD1  To generate 2-D finite element coefficients. 
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gncf2.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:12   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:44   pvcs
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
CD2    Rev 1.7   Wed May 01 14:31:42 1996   gaz
CD2 correctioin for mixed elements
CD2 
CD2    Rev 1.6   Thu Jan 18 09:36:54 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.5   Wed Jan 10 12:47:44 1996   hend
CD2 Ammended Purpose in Prolog
CD2 
CD2    Rev 1.4   Wed Jan 10 12:37:14 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   Tue Jan 09 14:05:34 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.2   03/18/94 15:57:40   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   03/08/94 15:45:58   llt
CD2 moved the allocation and deallocation of array aj to the top and bottom of
CD2 gencof (instead of inside gncf2 and gncf3).
CD2 
CD2    Rev 1.0   01/20/94 10:24:34   pvcs
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
CD4  element nele in group neu with integration point nga.  If icsh
CD4  is equal to 1 then the shape function and jacobian terms are
CD4  also calculated.
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
CD6   bcoef           REAL*8   fff    
CD6   dr              REAL*8   fbs    Contains weights for integration
CD6                                     points (bricks, rectangles)
CD6   eta             REAL*8   fbs    Local coordinates in a finite element of
CD6                                     the numerical integration points
CD6   istrs           INT      faai   Parameter indicating if the stress
CD6                                     solution is enabled
CD6   ivf             INT      faai   
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
CD6                                     connectivity array
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
CD6   wxr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to x (rectangles)
CD6   wy              REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y
CD6   wyr             REAL*8   fbs    Derivative of shape functions with
CD6                                     respect to y (rectangles)
CD6   xt              REAL*8   fbs    Parameters needed in element calculations
CD6   yt              REAL*8   fbs    Parameters needed in element calculations
CD6
CD6
CD6 Global Subprograms
CD6
CD6   Identifier      Type     Description
CD6   
CD6   pebi            N/A      Calculate internodal areas using perpendicular 
CD6                              bisectors.
CD6                              
CD6   shap2r          N/A      Evaluate 2-D finite element shape functions  
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
CD7   a11             REAL*8   
CD7   a12             REAL*8   
CD7   a21             REAL*8   
CD7   a22             REAL*8   
CD7   area            REAL*8  
CD7   area2           REAL*8  
CD7   b1              REAL*8   
CD7   b2              REAL*8   
CD7   b3              REAL*8   
CD7   b4              REAL*8   
CD7   bt              REAL*8   
CD7   ct              REAL*8   
CD7   detja           REAL*8   
CD7   dnga            REAL*8   
CD7   ij              INT
CD7   jt              INT
CD7   jz              INT
CD7   k               INT
CD7   kb              INT
CD7   knum            INT
CD7   sa11            REAL*8   
CD7   sa12            REAL*8   
CD7   sa21            REAL*8   
CD7   sa22            REAL*8   
CD7   third           REAL*8  
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

      integer icsh, nele, neu, neumax, nga, nrq, nsl
      integer ij, jt(6), jz, k, kb, knum 
      real*8  aj(neumax,*), bt(3), ct(3)
      real*8  a1, a2, a3, a4, a11, a12, a21, a22, area, area2
      real*8  b1, b2, b3, b4, detja, dnga
      real*8  sa11, sa12, sa21, sa22, third
      data jt/1,2,3,1,2,3/
c     
c     calculate jacobian information once and save it
c     
      if(icsh.eq.1) then
c     
         if(nsl.eq.4) then
c     procedure for 4 node quad elements
            call shap2r(nga)
            do jz=1,nsl
               kb=nelm((nele-1)*ns+jz)
               xt(jz)=cord(kb,1)
               yt(jz)=cord(kb,2)
            enddo
            
            a1=0.250*(-xt(1)+xt(2)+xt(3)-xt(4))
            a2=0.250*(xt(1)-xt(2)+xt(3)-xt(4))
            a3=0.250*(-xt(1)-xt(2)+xt(3)+xt(4))
            a4=a2
            b1=0.250*(-yt(1)+yt(2)+yt(3)-yt(4))
            b2=0.250*(yt(1)-yt(2)+yt(3)-yt(4))
            b3=0.250*(-yt(1)-yt(2)+yt(3)+yt(4))
            b4=b2
            
            sa11=a1+a2*eta(nga)
            sa12=b1+b2*eta(nga)
            sa21=a3+a4*si(nga)
            sa22=b3+b4*si(nga)
            a11=sa22
            a22=sa11
            a21=-sa21
            a12=-sa12
            detja=sa11*sa22-sa12*sa21
            aj(neu,1)=a11
            aj(neu,2)=a12
            aj(neu,4)=a21
            aj(neu,5)=a22
            aj(neu,10)=detja
         endif
      endif
c     
c     code for triangles
c     
      if(nga.le.1) then
         if(nsl.eq.3) then
c     
c     procedure for 3 node triangle
c     
c     define geometric coefficients
c     
            third=1./3.d0
            do jz=1,nsl
               kb=nelm((nele-1)*ns+jz)
               xt(jz)=cord(kb,1)
               yt(jz)=cord(kb,2)
            enddo
            if(ivf.eq.0.or.nrq.gt.3) then
               do jz=1,nsl
                  bt(jz)=yt(jt(jz+1))-yt(jt(jz+2))
                  ct(jz)=xt(jt(jz+2))-xt(jt(jz+1))
               enddo
               area2=ct(3)*bt(2)-ct(2)*bt(3)
               area=area2*0.5
               do jz=1,nsl
                  wx(nga,jz)=bt(jz)/area2
                  wy(nga,jz)=ct(jz)/area2
                  w(nga,jz)=third
               enddo
            else
               call pebi(xt,yt,bt,ct)
c     bt are the areas/delx
c     ct are the volumes
               if(nrq.eq.1) then
c     calculate volumes
                  bcoef(neu,1)=ct(1)
                  bcoef(neu,2)=0.0
                  bcoef(neu,3)=0.0
                  bcoef(neu,4)=0.0
                  bcoef(neu,5)=ct(2)
                  bcoef(neu,6)=0.0
                  bcoef(neu,7)=0.0
                  bcoef(neu,8)=0.0
                  bcoef(neu,9)=ct(3)
                  goto 999
               else if(nrq.eq.2) then
c     calculate x area/dx
                  bcoef(neu,1)=0.0
                  bcoef(neu,2)=-bt(1)
                  bcoef(neu,3)=-bt(3)
                  bcoef(neu,4)=-bt(1)
                  bcoef(neu,5)=0.0
                  bcoef(neu,6)=-bt(2)
                  bcoef(neu,7)=-bt(3)
                  bcoef(neu,8)=-bt(2)
                  bcoef(neu,9)=0.0
                  goto 999
               else if(nrq.eq.3) then
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
            endif
         endif
      endif
c     
c     
c     code for quads
c     
      if(nsl.eq.4) then
c     
c     identify weighting factor for numerical integration
c     
         dnga=dr(nga)
c     
c     assemble global derivatives
c     
         if(nrq.eq.1) then
c     nk*njz
               do k=1,nsl
                  knum=(k-1)*nsl
c                  do jz=1,nsl
                     ij=knum+k
                     bcoef(neu,ij)=bcoef(neu,ij)+wr(nga,k)*
     &                 dnga*aj(neu,10)
c                    ij = knum + k
c                    bcoef(neu,ij)=bcoef(neu,ij)+w(nga,k)*dnga*aj(neu,10)
c                  enddo
               enddo
        endif        
       if(nrq.eq.2.or.nrq.eq.15) then     
c     (d(nk/dx)*(dnjz/dx)
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+(a11*wxr(nga,k)+a12*
     &                 wyr(nga,k))*
     5                 (a11*wxr(nga,jz)+a12*wyr(nga,jz))*dnga/detja
               enddo
            enddo
         endif
       if(nrq.eq.3.or.nrq.eq.16) then  
c     (d(nk/dy)*(dnjz/dy)
c     recall jacobian information
            a21=aj(neu,4)
            a22=aj(neu,5)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+(a21*wxr(nga,k)+a22*
     &                 wyr(nga,k))*
     5                 (a21*wxr(nga,jz)+a22*wyr(nga,jz))*dnga/detja
               enddo
            enddo
         endif

c     c for 2-d stress 
         if(istrs.ne.0) then
            if(nrq.eq.11) then
c     (dnk/dx)*(dnjz/dy)
c     recall jacobian information
               a11=aj(neu,1)
               a12=aj(neu,2)
               a21=aj(neu,4)
               a22=aj(neu,5)
               detja=aj(neu,10)
               do k=1,nsl
                  knum=(k-1)*nsl
                  do jz=1,nsl
                     ij=knum+jz
                     bcoef(neu,ij)=bcoef(neu,ij)+
     *                    (a11*wxr(nga,k)+a12*wyr(nga,k))*
     5                    (a21*wxr(nga,jz)+a22*wyr(nga,jz))*dnga/detja
                  enddo
               enddo
            endif
            if(nrq.eq.12) then
c     (dnk/dy)*(dnjz/dx)
c     recall jacobian information
               a11=aj(neu,1)
               a12=aj(neu,2)
               a21=aj(neu,4)
               a22=aj(neu,5)
               detja=aj(neu,10)
               do k=1,nsl
                  knum=(k-1)*nsl
                  do jz=1,nsl
                     ij=knum+jz
                     bcoef(neu,ij)=bcoef(neu,ij)+
     *                    (a21*wxr(nga,k)+a22*wyr(nga,k))*
     5                    (a11*wxr(nga,jz)+a12*wyr(nga,jz))*dnga/detja
                  enddo
               enddo
            endif
         if(nrq.eq.13) then
c                     
c     nk*(dnj/dx)
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 wr(nga,k)*(a11*wxr(nga,jz)+a12*wyr(nga,jz))*dnga
               enddo
            enddo
         endif         
         if(nrq.eq.14) then
c     
c     nk*(dnj/dy)
c     recall jacobian information
            a21=aj(neu,4)
            a22=aj(neu,5)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 wr(nga,k)*(a21*wxr(nga,jz)+a22*wyr(nga,jz))*dnga
               enddo
            enddo
         endif  
         if(nrq.eq.17) then
c     
c     (dnk/dx)*nj
c     recall jacobian information
            a11=aj(neu,1)
            a12=aj(neu,2)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 wr(nga,jz)*(a11*wxr(nga,k)+a12*wyr(nga,k))*dnga
               enddo
            enddo
         endif
         if(nrq.eq.18) then
c     
c     (dnk/dy)*nj
c     recall jacobian information
            a21=aj(neu,4)
            a22=aj(neu,5)
            detja=aj(neu,10)
            do k=1,nsl
               knum=(k-1)*nsl
               do jz=1,nsl
                  ij=knum+jz
                  bcoef(neu,ij)=bcoef(neu,ij)+
     *                 wr(nga,jz)*(a21*wxr(nga,k)+a22*wyr(nga,k))*dnga
               enddo
            enddo
         endif                
         endif
      endif
c     
c     code for triangles
c     
      if(nga.le.1) then

         if(nsl.eq.3) then
c     
c     assemble global derivatives
c     
            if(nrq.eq.1) then
c     nk*njz
               do k=1,nsl
                  knum=(k-1)*nsl
c                  do jz=1,nsl
c                     ij=knum+jz
c                     bcoef(neu,ij)=bcoef(neu,ij)+w(nga,k)*w(nga,jz)*area
                     ij = knum + k
                     bcoef(neu,ij)=bcoef(neu,ij)+w(nga,k)*area
c                  enddo
               enddo
            endif
            if(nrq.eq.2) then
c     (d(nk/dx)*(dnjz/dx)
               do k=1,nsl
                  knum=(k-1)*nsl
                  do jz=1,nsl
                     ij=knum+jz
                     bcoef(neu,ij)=bcoef(neu,ij)+wx(nga,k)*wx(nga,jz)*
     &                    area
                  enddo
               enddo
            endif
            if(nrq.eq.3) then
c     (d(nk/dy)*(dnjz/dy)
               do k=1,nsl
                  knum=(k-1)*nsl
                  do jz=1,nsl
                     ij=knum+jz
                     bcoef(neu,ij)=bcoef(neu,ij)+wy(nga,k)*wy(nga,jz)*
     &                    area
                  enddo
               enddo
            endif
            if(nrq.eq.5) then
c     nk*(dnjz/dx)
               do k=1,nsl
                  knum=(k-1)*nsl
                  do jz=1,nsl
                     ij=knum+jz
                     bcoef(neu,ij)=bcoef(neu,ij)+
     *                    w(nga,k)*wx(nga,jz)*area
                  enddo
               enddo
            endif
            if(nrq.eq.6) then
c     nk*(dnjz/dy)
               do k=1,nsl
                  knum=(k-1)*nsl
                  do jz=1,nsl
                     ij=knum+jz
                     bcoef(neu,ij)=bcoef(neu,ij)+
     *                    w(nga,k)*wy(nga,jz)*area
                  enddo
               enddo
            endif
c     c for 2-d stress
            if(istrs.ne.0) then
               if(nrq.eq.11) then
c     (dnk/dx)*(dnjz/dy)
                  do k=1,nsl
                     knum=(k-1)*nsl
                     do jz=1,nsl
                        ij=knum+jz
                        bcoef(neu,ij)=bcoef(neu,ij)+wx(nga,k)*
     &                       wy(nga,jz)*area
                     enddo
                  enddo
               endif
               if(nrq.eq.12) then
c     (dnk/dy)*(dnjz/dx)
                  do k=1,nsl
                     knum=(k-1)*nsl
                     do jz=1,nsl
                        ij=knum+jz
                        bcoef(neu,ij)=bcoef(neu,ij)+wy(nga,k)*
     &                       wx(nga,jz)*area
                     enddo
                  enddo
               endif
            endif
         endif
      endif
      
 999  continue
      
      return
      end

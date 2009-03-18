      subroutine split(iflg)
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
CD1  To split rectangles (bricks) into triangles (tetrahedrals) and
CD1  average so grid orientation is not present.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/split.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Thu Jan 18 10:50:36 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.6   Fri Jan 12 17:58:32 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.5   Thu Jan 11 11:39:28 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   Tue Jan 09 14:13:42 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.3   01/28/95 13:55:48   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.2   05/11/94 16:13:42   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.1   03/18/94 15:58:06   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:12   pvcs
CD2 original version in process of being certified
CD2 
c got rid of b matrix storage gaz
c added return for 6-noded prisms
c 12/21/94 gaz rearranged some memory calls
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
CD5   iflg            INT      I    Flag denoting whether elements are being
CD5                                   split or averaged.
CD5                                   iflg = 0, split elements
CD5                                   iflg = 1, average elements
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
CD6   cord            REAL*8   fbs    Contains the coordinates of each node
CD6   icnl            INT      faai   Problem dimension
CD6   ivf             INT      faai   Flag for finite volume calculations
CD6   n               INT      faai   Total number of nodes
CD6   nei             INT      faai   Total number of elements in the problem
CD6   nelm            INT      fbb    Information about nodes in each element
CD6   nelmd           INT      param  Maximum array space for element
CD6                                     connectivity array
CD6   ns              INT      faai   Number of nodes per element
CD6   sx              REAL*8   fbc    Contains finite element geometric
CD6                                     coefficients necessary for heat and mass
CD6                                     transfer simulation
CD6   sx1             REAL*8   fbc    Contains volume associated with each node
CD6
CD6
CD6 Global Subprograms
CD6
CD6   Identifier      Type     Description
CD6
CD6   mmgetblk        N/A      Allocate memory to an array
CD6   mmrelblk        N/A      Deallocates space for an array
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
CD7   dis13           REAL*8   
CD7   dis17           REAL*8   
CD7   dis24           REAL*8   
CD7   dis28           REAL*8   
CD7   dis35           REAL*8   
CD7   dis46           REAL*8   
CD7   dismin          REAL*8   
CD7   dist1x          REAL*8   
CD7   dist2x          REAL*8   
CD7   dist1y          REAL*8   
CD7   dist2y          REAL*8   
CD7   i               INT      
CD7   icode           INT      
CD7   ie              INT      
CD7   ipnelmt         POINTER   Pointer to array nelmt
CD7   istrwe          INT      
CD7   j               INT      
CD7   nelmt           INT      
CD7   nelt            INT      
CD7   neltt           INT      
CD7   nept            INT      
CD7   nloc            INT      
CD7   npoine          INT      
CD7   x1              REAL*8   
CD7   x2              REAL*8   
CD7   x3              REAL*8   
CD7   x4              REAL*8   
CD7   x5              REAL*8   
CD7   x6              REAL*8   
CD7   x7              REAL*8   
CD7   x8              REAL*8   
CD7   y1              REAL*8   
CD7   y2              REAL*8   
CD7   y3              REAL*8   
CD7   y4              REAL*8   
CD7   y5              REAL*8   
CD7   y6              REAL*8   
CD7   y7              REAL*8   
CD7   y8              REAL*8   
CD7   z1              REAL*8   
CD7   z2              REAL*8   
CD7   z3              REAL*8   
CD7   z4              REAL*8   
CD7   z5              REAL*8   
CD7   z6              REAL*8   
CD7   z7              REAL*8   
CD7   z8              REAL*8   
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************

      use combi
      use comdti
      use comai
      implicit none

      integer icode,ie,iflg,istrwe,i,j,nloc(8),nelt,neltt,nept,npoine
      integer, allocatable :: nelmt(:)
      real*8 dis13, dis17, dis24, dis28, dis35, dis46, dismin
      real*8 dist1x, dist2x, dist1y, dist2y
      real*8 x1, x2, x3, x4, x5, x6, x7, x8
      real*8 y1, y2, y3, y4, y5, y6, y7, y8
      real*8 z1, z2, z3, z4, z5, z6, z7, z8

      if(iflg.eq.0) then
c     split elements
c     return if 3-d and ns = 4 (can't divide tets further)
         if(icnl.eq.0.and.ns.eq.4) return
c     return if 3-d and ns = 6 (can't divide prisms  further)
         if(icnl.eq.0.and.ns.eq.6) return
c     return if 2-d and ns = 3 (can't divide triangles further)
         if(icnl.ne.0.and.ns.eq.3) return
c     return if ns <= 2 can't divide 1-d elements
         if (ns .le. 2) return
         if(icnl.ne.0) then
c     2-d elements
c     nelt=element count
            if(ivf.eq.0) then
               allocate(nelmt(nei*12))
            else if(ivf.gt.0) then
               allocate(nelmt(nei*6))
            endif
            nelt=0
            do 10 ie=1,nei
c     identify nodes in element ie
               npoine=(ie-1)*ns
               nloc(1)=nelm(npoine+1)
               nloc(2)=nelm(npoine+2)
               nloc(3)=nelm(npoine+3)
               nloc(4)=nelm(npoine+4)
               if(ivf.eq.0) then
c     check if triangle is already defined
                  if(nloc(4).ne.nloc(3).and.nloc(4).ne.0) then
c     divide rectangles into 4 triangles
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(3)
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(3)
                     nelmt(nept+3)=nloc(4)
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(4)
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(3)
                     nelmt(nept+3)=nloc(4)
                  else if(nloc(4).eq.0) then
c     already a triangle
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(3)
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(3)
                  endif
               else if(ivf.gt.0) then

c     for pebi,divide along small diagonal
                  if(nloc(4).ne.nloc(3).and.nloc(4).ne.0) then
                     dist1x=cord(nloc(1),1)-cord(nloc(3),1)
                     dist2x=cord(nloc(2),1)-cord(nloc(4),1)
                     dist1y=cord(nloc(1),2)-cord(nloc(3),2)
                     dist2y=cord(nloc(2),2)-cord(nloc(4),2)
                     dis13=sqrt(dist1x*dist1x+dist1y*dist1y)
                     dis24=sqrt(dist2x*dist2x+dist2y*dist2y)
                     if(dis13.le.dis24) then
c     divide rectangles into 2 triangles
                        nelt=nelt+1
                        nept=(nelt-1)*3
                        nelmt(nept+1)=nloc(1)
                        nelmt(nept+2)=nloc(2)
                        nelmt(nept+3)=nloc(3)
                        nelt=nelt+1
                        nept=(nelt-1)*3
                        nelmt(nept+1)=nloc(1)
                        nelmt(nept+2)=nloc(3)
                        nelmt(nept+3)=nloc(4)
                     else
                        nelt=nelt+1
                        nept=(nelt-1)*3
                        nelmt(nept+1)=nloc(1)
                        nelmt(nept+2)=nloc(2)
                        nelmt(nept+3)=nloc(4)
                        nelt=nelt+1
                        nept=(nelt-1)*3
                        nelmt(nept+1)=nloc(2)
                        nelmt(nept+2)=nloc(3)
                        nelmt(nept+3)=nloc(4)
                     endif
                  else
c     already a triangle
                     nelt=nelt+1
                     nept=(nelt-1)*3
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(3)
                  endif
               endif
 10         continue

c     set new ns and nei
            ns=3
            nei=nelt
            nemx=12
            neigh=12

c     set new elements
            neltt=nelt*3
            deallocate(nelm)
            allocate(nelm(neltt))
            do ie=1,neltt
               nelm(ie)=nelmt(ie)
            enddo

         else if(icnl.eq.0) then
c     3-d elements
            if(ivf.eq.-99) then
               allocate(nelmt(nei*96))
            else
               allocate(nelmt(nei*24))
            endif
            nelt=0
            do 50 ie=1,nei
c     identify local nodes
               npoine=(ie-1)*ns
c               nloc(1)=nelm(npoine+1)
c               nloc(2)=nelm(npoine+2)
c               nloc(3)=nelm(npoine+3)
c               nloc(4)=nelm(npoine+4)
c               nloc(5)=nelm(npoine+5)
c               nloc(6)=nelm(npoine+6)
c               nloc(7)=nelm(npoine+7)
c               nloc(8)=nelm(npoine+8)
c
	         nloc(5)=nelm(npoine+1)
               nloc(6)=nelm(npoine+2)
               nloc(7)=nelm(npoine+3)
               nloc(8)=nelm(npoine+4)
               nloc(1)=nelm(npoine+5)
               nloc(2)=nelm(npoine+6)
               nloc(3)=nelm(npoine+7)
               nloc(4)=nelm(npoine+8)
c     divide into tetrahedrals
               if(ivf.eq.-99) then
c     only doing this for brick elements
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(6)
                  nelmt(nept+3)=nloc(5)
                  nelmt(nept+4)=nloc(8)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(2)
                  nelmt(nept+3)=nloc(6)
                  nelmt(nept+4)=nloc(8)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(4)
                  nelmt(nept+3)=nloc(2)
                  nelmt(nept+4)=nloc(8)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(2)
                  nelmt(nept+4)=nloc(8)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(2)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(6)
                  nelmt(nept+4)=nloc(8)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(3)
                  nelmt(nept+2)=nloc(7)
                  nelmt(nept+3)=nloc(6)
                  nelmt(nept+4)=nloc(8)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(4)
                  nelmt(nept+3)=nloc(2)
                  nelmt(nept+4)=nloc(6)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(2)
                  nelmt(nept+4)=nloc(6)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(1)
                  nelmt(nept+3)=nloc(8)
                  nelmt(nept+4)=nloc(6)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(5)
                  nelmt(nept+3)=nloc(8)
                  nelmt(nept+4)=nloc(6)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(3)
                  nelmt(nept+2)=nloc(4)
                  nelmt(nept+3)=nloc(8)
                  nelmt(nept+4)=nloc(6)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(7)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(8)
                  nelmt(nept+4)=nloc(6)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(1)
                  nelmt(nept+4)=nloc(5)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(7)
                  nelmt(nept+3)=nloc(3)
                  nelmt(nept+4)=nloc(5)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(3)
                  nelmt(nept+2)=nloc(2)
                  nelmt(nept+3)=nloc(1)
                  nelmt(nept+4)=nloc(5)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(8)
                  nelmt(nept+3)=nloc(7)
                  nelmt(nept+4)=nloc(5)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(2)
                  nelmt(nept+2)=nloc(7)
                  nelmt(nept+3)=nloc(6)
                  nelmt(nept+4)=nloc(5)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(2)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(7)
                  nelmt(nept+4)=nloc(5)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(1)
                  nelmt(nept+4)=nloc(7)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(3)
                  nelmt(nept+3)=nloc(2)
                  nelmt(nept+4)=nloc(7)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(1)
                  nelmt(nept+2)=nloc(2)
                  nelmt(nept+3)=nloc(5)
                  nelmt(nept+4)=nloc(7)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(2)
                  nelmt(nept+2)=nloc(6)
                  nelmt(nept+3)=nloc(5)
                  nelmt(nept+4)=nloc(7)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(1)
                  nelmt(nept+3)=nloc(5)
                  nelmt(nept+4)=nloc(7)
                  nelt=nelt+1
                  nept=(nelt-1)*4
                  nelmt(nept+1)=nloc(4)
                  nelmt(nept+2)=nloc(5)
                  nelmt(nept+3)=nloc(8)
                  nelmt(nept+4)=nloc(7)
               else
c     pebi tetrahedrals
c     divide hexahedrons into tetrahedrons
c     form main diagonals,chech for smallest
c     
c     identify coordinates
                  x1=cord(nloc(1),1)
                  x2=cord(nloc(2),1)
                  x3=cord(nloc(3),1)
                  x4=cord(nloc(4),1)
                  x5=cord(nloc(5),1)
                  x6=cord(nloc(6),1)
                  x7=cord(nloc(7),1)
                  x8=cord(nloc(8),1)
                  y1=cord(nloc(1),2)
                  y2=cord(nloc(2),2)
                  y3=cord(nloc(3),2)
                  y4=cord(nloc(4),2)
                  y5=cord(nloc(5),2)
                  y6=cord(nloc(6),2)
                  y7=cord(nloc(7),2)
                  y8=cord(nloc(8),2)
                  z1=cord(nloc(1),3)
                  z2=cord(nloc(2),3)
                  z3=cord(nloc(3),3)
                  z4=cord(nloc(4),3)
                  z5=cord(nloc(5),3)
                  z6=cord(nloc(6),3)
                  z7=cord(nloc(7),3)
                  z8=cord(nloc(8),3)

c     check main diagonal
                  dis17=(x1-x7)**2 + (y1-y7)**2 + (z1-z7)**2
                  dis35=(x3-x5)**2 + (y3-y5)**2 + (z3-z5)**2
                  dis46=(x4-x6)**2 + (y4-y6)**2 + (z4-z6)**2
                  dis28=(x8-x2)**2 + (y8-y2)**2 + (z8-z2)**2
                  dismin=min(dis17,dis35,dis46,dis28)

                  if(dis17.eq.dismin) then
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(3)
                     nelmt(nept+4)=nloc(7)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(3)
                     nelmt(nept+3)=nloc(4)
                     nelmt(nept+4)=nloc(7)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(4)
                     nelmt(nept+3)=nloc(8)
                     nelmt(nept+4)=nloc(7)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(6)
                     nelmt(nept+3)=nloc(2)
                     nelmt(nept+4)=nloc(7)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(5)
                     nelmt(nept+3)=nloc(6)
                     nelmt(nept+4)=nloc(7)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(1)
                     nelmt(nept+2)=nloc(8)
                     nelmt(nept+3)=nloc(5)
                     nelmt(nept+4)=nloc(7)
                  else if(dis35.eq.dismin) then
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(3)
                     nelmt(nept+2)=nloc(4)
                     nelmt(nept+3)=nloc(1)
                     nelmt(nept+4)=nloc(5)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(3)
                     nelmt(nept+2)=nloc(1)
                     nelmt(nept+3)=nloc(2)
                     nelmt(nept+4)=nloc(5)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(3)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(6)
                     nelmt(nept+4)=nloc(5)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(3)
                     nelmt(nept+2)=nloc(6)
                     nelmt(nept+3)=nloc(7)
                     nelmt(nept+4)=nloc(5)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(3)
                     nelmt(nept+2)=nloc(4)
                     nelmt(nept+3)=nloc(2)
                     nelmt(nept+4)=nloc(5)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(3)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(7)
                     nelmt(nept+4)=nloc(5)
                  else if(dis46.eq.dismin) then
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(4)
                     nelmt(nept+2)=nloc(1)
                     nelmt(nept+3)=nloc(2)
                     nelmt(nept+4)=nloc(6)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(4)
                     nelmt(nept+2)=nloc(2)
                     nelmt(nept+3)=nloc(3)
                     nelmt(nept+4)=nloc(6)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(4)
                     nelmt(nept+2)=nloc(5)
                     nelmt(nept+3)=nloc(1)
                     nelmt(nept+4)=nloc(6)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(4)
                     nelmt(nept+2)=nloc(8)
                     nelmt(nept+3)=nloc(5)
                     nelmt(nept+4)=nloc(6)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(4)
                     nelmt(nept+2)=nloc(3)
                     nelmt(nept+3)=nloc(7)
                     nelmt(nept+4)=nloc(6)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(4)
                     nelmt(nept+2)=nloc(7)
                     nelmt(nept+3)=nloc(8)
                     nelmt(nept+4)=nloc(6)
                  else if(dis28.eq.dismin) then
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(4)
                     nelmt(nept+3)=nloc(1)
                     nelmt(nept+4)=nloc(8)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(1)
                     nelmt(nept+3)=nloc(5)
                     nelmt(nept+4)=nloc(8)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(3)
                     nelmt(nept+3)=nloc(4)
                     nelmt(nept+4)=nloc(8)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(5)
                     nelmt(nept+3)=nloc(6)
                     nelmt(nept+4)=nloc(8)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(6)
                     nelmt(nept+3)=nloc(7)
                     nelmt(nept+4)=nloc(8)
                     nelt=nelt+1
                     nept=(nelt-1)*4
                     nelmt(nept+1)=nloc(2)
                     nelmt(nept+2)=nloc(7)
                     nelmt(nept+3)=nloc(3)
                     nelmt(nept+4)=nloc(8)
                  endif
               endif
 50         continue

c     set new ns and nei
            ns=4
            nei=nelt
            nemx=96
            neigh=27

c     set new elements
            neltt=nelt*4
            deallocate(nelm)
            allocate(nelm(neltt))
            do ie=1,neltt
               nelm(ie)=nelmt(ie)
            enddo
         endif
         deallocate(nelmt)
      else if(iflg.eq.1) then

c     divide coefficients by 4 because 4 times the volume was used
c     to form the 24 tetrahedrals above
c     do averaging calcs only if ivf=0
         istrwe=iad
         if(icnl.ne.0.and.ivf.eq.0) then
            do i=1,2
               do j=1,istrwe
                  sx(j,i)=sx(j,i)*0.5
               enddo
            enddo
c     now reduce volumes
            do i=1,n
               sx1(i)=sx1(i)*0.5
            enddo
         else if(icnl.eq.0.and.ivf.eq.-99) then
            do i=1,3
               do j=1,istrwe
                  sx(j,i)=sx(j,i)/16.0
               enddo
            enddo
c     now reduce volumes
            do i=1,n
               sx1(i)=sx1(i)/16.0
            enddo
         endif
      endif
      
 999  continue
      
      return
      end

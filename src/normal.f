      subroutine  normal  ( iflg,fdum2 )
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine normalizes the matrix equations and returns the 
CD1  sum of the residuals.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/normal.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:28   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Wed Jan 10 14:52:24 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.3   11/16/95 10:29:18   llt
CD2 min function doesn't work inside a statement function on the ibm.
CD2 
CD2    Rev 1.2   09/29/95 16:11:10   llt
CD2 added small number on dividing -- was dividing by zerol
CD2 getting min of fadn & fbdn function -- was infinity
CD2 
CD2    Rev 1.1   03/18/94 15:59:40   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:58   pvcs
CD2 original version in process of being certified
CD2 
c version FEHM5.1J changes
c 16-july-91
c defered zeroing out dfmp(i) etc to geneq1 (took it out here
c 16-july-91
c set tollr=tmch*0.1 if appropriate
c 23-july-91
c zeroed out appropriate derivatives in gensl1(4) for zero porosity
c 19-august-91
c changed equivalence statements to reflect smaller common/fcc/
c 11-nov-91
c changed gensl1 to reflect symmetry in assembly(geneq1)
c 12-nov-91
c changed gensl4 to have changes made in gensl1
c 13-nov-91
c changed gensl3 to reflect changes in geneq3
c 15-nov-91
c changed nr3 to nrhs(3) in gensl4
c 16-nov-91
c changed nsizea=nmat(8).. to =nmat(9)..
c 29-nov-91
c changed to generic intrisic functions
c 23-dec-91
c cleaned up some coding in rdof routines
c added some new code to rd2dof
c 24-dec-91
c deleted do 40 loop(set sol1,sol2=0)
c 30-dec-91
c still playing with residual adjustment in rd2dof
c added other call to solve in rd2dof
c put in decision making ability in rd2dof
c 2-jan-92
c put changes in rd2dof ,do loops before 2nd call to solve
c 13-jan-92
c more changes to rd2dof :iteration on back subafter solve
c 16-jan-92
c started modifying rd3dof to be like rd2dof
c modified equivalence statement in gensl4
c 18-jan-92
c some mods of rd3dof
c put returns before call to solve if residuals are small(rd2dof,rd3dof)
c 26-feb-92
c modified coding near call to geneq3 to position equations correctly
c actually changed call to geneq3 to geneq1 in gensl1 and genslc
c 2-feb-93
c took out 1.0e-20 in fdum2 sum ,replace with tollah    
c 11-feb-93
c added subroutine normal_uc for uncoupled normalization
c 3-apr-93
c added variable switching in gensl4 (p,s,t) to (p,t,s)
c 7-apr-93
c took out bpc in gensl4
c 9-apr-93
c put in switch routines for equation ordering
c 12-apr-93
c put nmat stuff in normal ,got rid of nb
c 12-april-93
c delete code for iflg=1 from normal(irdof?)
c 27-apr-93
c added islord=-1 option
c changed call ti switchb,changed gensl1,gensl2
c 7-may-93 llt
c equivalences changed to pointers
c 7-may-93 llt
c replaced gotos embedded in code with a goto at end of subroutine
c 8-june-93
c fixed call to rd3dof
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  3.4.2     Solve Nonlinear Equation at Each Time Step
CD3
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
C**********************************************************************

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comhi
      use davidi
      use comdti
      use comai

      implicit none
      integer i, i1, i2, idgm, iflg, ikd, ilq 
      integer j, jm, neq2, neqp1
      real*8  a1, a2, a3, a4, a5, ali, ali3, alm, alm3
      real*8  a11, a12, a13, a21, a22, a23, a31, a32, a33 
      real*8  an11, an12, an21, an22
      real*8  b1, b2, b3, det
      real*8  fad, fadn, fbd, fbdn, fcd, fcdn, fdm2i, fdum2
      real*8  p1, p2, p3, p4, p5, p6, p7, p8, p9, p11, p12, p21, p22
      real*8  tollah
      parameter( tollah=1.d-40)
      
c**** linear algebra ****
      alm (a1,a2,b1,b2)       =       a1*b1+a2*b2
      ali (a1,a2,a3,a4,b1)    =  b1/( (a1*a4-a2*a3) + 1e-30 )
      ali3(a1,a2,a3,a4,a5)    =     ( a1*a2-a3*a4       )/a5
      alm3(a1,a2,a3,b1,b2,b3) =       a1*b1+a2*b2+a3*b3
      
      neqp1  =  neq+1
      if ( iflg .eq. 0 )  then

c**** heat and mass transfer solution ****
         fdum2  =  0.0
         do i=1,neq
            i1     =  nelm(i  )+1
            i2     =  nelm(i+1)
            idgm   =  nelmdg(i)-neqp1
            a11    =  a(idgm+nmat(1))
            a12    =  a(idgm+nmat(2))
            a21    =  a(idgm+nmat(3))
            a22    =  a(idgm+nmat(4))

            p11    =  ali( a11,a12,a21,a22,a22 )
            p12    = -ali( a11,a12,a21,a22,a12 )
            p21    = -ali( a11,a12,a21,a22,a21 )
            p22    =  ali( a11,a12,a21,a22,a11 )

            fad    =  bp(i+nrhs(1))
            fbd    =  bp(i+nrhs(2))
            fadn   =  alm( p11,p12,fad,fbd )
            fbdn   =  alm( p21,p22,fad,fbd )

            bp(i+nrhs(1)) =  fadn
            bp(i+nrhs(2)) =  fbdn

            do j=i1,i2
               jm =  j-neqp1
               if (jm .eq. idgm )  then
                  a(jm+nmat(1)) =  1.0
                  a(jm+nmat(2)) =  0.0
                  a(jm+nmat(3)) =  0.0
                  a(jm+nmat(4)) =  1.0
               else
                  a11    =  a(jm+nmat(1))
                  a12    =  a(jm+nmat(2))
                  a21    =  a(jm+nmat(3))
                  a22    =  a(jm+nmat(4))
                  an11   =  alm( p11,p12,a11,a21 )
                  an12   =  alm( p11,p12,a12,a22 )
                  an21   =  alm( p21,p22,a11,a21 )
                  an22   =  alm( p21,p22,a12,a22 )
                  a(jm+nmat(1)) =  an11
                  a(jm+nmat(2)) =  an12
                  a(jm+nmat(3)) =  an21
                  a(jm+nmat(4)) =  an22
               endif
            enddo
            fdum2  =  fdum2+fadn**2+fbdn**2
         enddo
      endif
      if ( iflg .ne .0 )  then
c**** heat , mass , and noncondensible present - full jacobian ****
         neq2   =  neq*2
         fdum2  =  0.0
         do i=1,neq
            i1     =  nelm(i  )+1
            i2     =  nelm(i+1)
            idgm   =  nelmdg(i)-neqp1
            a11    =  a(idgm+nmat(1))
            a12    =  a(idgm+nmat(2))
            a13    =  a(idgm+nmat(3))
            a21    =  a(idgm+nmat(4))
            a22    =  a(idgm+nmat(5))
            a23    =  a(idgm+nmat(6))
            a31    =  a(idgm+nmat(7))
            a32    =  a(idgm+nmat(8))
            a33    =  a(idgm+nmat(9))
            fad    =  bp(i+nrhs(1))
            fbd    =  bp(i+nrhs(2))
            fcd    =  bp(i+nrhs(3))
            det    =   a11*a22*a33+a12*a23*a31+a13*a21*a32
     *           -a13*a22*a31-a11*a23*a32-a12*a21*a33
            p1     =  ali3( a22,a33,a23,a32,det )
            p2     =  ali3( a13,a32,a12,a33,det )
            p3     =  ali3( a12,a23,a13,a22,det )
            p4     =  ali3( a23,a31,a21,a33,det )
            p5     =  ali3( a11,a33,a13,a31,det )
            p6     =  ali3( a13,a21,a11,a23,det )
            p7     =  ali3( a21,a32,a22,a31,det )
            p8     =  ali3( a12,a31,a11,a32,det )
            p9     =  ali3( a11,a22,a12,a21,det )
            fadn   =  alm3( p1,p2,p3,fad,fbd,fcd )
            fbdn   =  alm3( p4,p5,p6,fad,fbd,fcd )
            fcdn   =  alm3( p7,p8,p9,fad,fbd,fcd )
            fdm2i  =  fadn**2+fbdn**2+fcdn**2+tollah  
            bp(i+nrhs(1)) = fadn
            bp(i+nrhs(2)) = fbdn
            bp(i+nrhs(3)) = fcdn
            do ikd=i1,i2
               ilq    =  ikd-neqp1
               if (ilq .eq. idgm )  then
                  a(ilq+nmat(1)) =  1.0
                  a(ilq+nmat(2)) =  0.0
                  a(ilq+nmat(3)) =  0.0
                  a(ilq+nmat(4)) =  0.0
                  a(ilq+nmat(5)) =  1.0
                  a(ilq+nmat(6)) =  0.0
                  a(ilq+nmat(7)) =  0.0
                  a(ilq+nmat(8)) =  0.0
                  a(ilq+nmat(9)) =  1.0
               else
                  a11    =  a(ilq+nmat(1))
                  a12    =  a(ilq+nmat(2))
                  a13    =  a(ilq+nmat(3))
                  a21    =  a(ilq+nmat(4))
                  a22    =  a(ilq+nmat(5))
                  a23    =  a(ilq+nmat(6))
                  a31    =  a(ilq+nmat(7))
                  a32    =  a(ilq+nmat(8))
                  a33    =  a(ilq+nmat(9))
                  a(ilq+nmat(1)) =  alm3(p1,p2,p3,a11,a21,a31)
                  a(ilq+nmat(2)) =  alm3(p1,p2,p3,a12,a22,a32)
                  a(ilq+nmat(3)) =  alm3(p1,p2,p3,a13,a23,a33)
                  a(ilq+nmat(4)) =  alm3(p4,p5,p6,a11,a21,a31)
                  a(ilq+nmat(5)) =  alm3(p4,p5,p6,a12,a22,a32)
                  a(ilq+nmat(6)) =  alm3(p4,p5,p6,a13,a23,a33)
                  a(ilq+nmat(7)) =  alm3(p7,p8,p9,a11,a21,a31)
                  a(ilq+nmat(8)) =  alm3(p7,p8,p9,a12,a22,a32)
                  a(ilq+nmat(9)) =  alm3(p7,p8,p9,a13,a23,a33)
               endif
            enddo
            fdum2 =  fdum2+fdm2i
         enddo
      endif

      return
      end

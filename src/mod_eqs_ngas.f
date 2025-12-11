      subroutine mod_eqs_ngas(neq,nelm,nelmdg,a,nmat,r,nrhs,phi,pci,
     &                   pcp,dpcef,t,ieos)
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
CD1  This subroutine modifies the Newton-Raphson algebraic equations
CD1  for the case when no noncondensible gas is present. It prevents
CD1  the matrix from possibly going singular in that case.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/mod_eqs_ngas.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:26   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:24   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Mon Mar 31 08:39:56 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.1   Thu Jun 27 12:59:42 1996   gaz
CD2 added prolog
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.5.2 Solve nonlinear equation set at each time step
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
c
c modify equations if no gas is present

      implicit none

      real*8 a(*),r(*),pci(*),phi(*)
      real*8 psatl,pcp(*),dpcef(*),t(*)
      real*8 dpsatt,dpsats,toll3
      integer nelm(*),nelmdg(*),nmat(*),nrhs(*),ieos(*)
      integer neq,neqp1,i,ii,ic,ipiv  
      integer i1,i2,jj,i3,i4,kb,nrhs3
c
      toll3=1.d-13
      nrhs3=nrhs(3)
      neqp1=neq+1
c
      do i=1,neq
       if(abs(r(i+nrhs3)).le.toll3.and.pci(i).le.0.0
     &    .and.ieos(i).eq.2) then
        i1=nelm(i)+1
        i2=nelm(i+1)
         do jj=i1,i2
            ic=jj-neqp1
            a(nmat(3)+ic)=0.0
            a(nmat(6)+ic)=0.0
            a(nmat(7)+ic)=0.0
            a(nmat(8)+ic)=0.0 
            a(nmat(9)+ic)=0.0
          enddo
        ipiv=nelmdg(i)-neqp1
        a(nmat(3)+ipiv)=0.0
        a(nmat(6)+ipiv)=0.0
        r(i+nrhs3)= phi(i)-psatl(t(i),pcp(i),dpcef(i),dpsatt,dpsats,0)
        a(nmat(7)+ipiv)=1.0
        a(nmat(8)+ipiv)=-dpsats
        a(nmat(9)+ipiv)=-dpsatt
       endif
      enddo
c
      return
      end

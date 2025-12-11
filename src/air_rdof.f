      subroutine air_rdof(iflg,irdof,ndummy,nelm,nmat,nrhs,a,bp,sx1)
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
CD1  This subroutine reduces the degree of freedom solution
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/air_rdof.f_a  $
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Mon Apr 14 12:40:38 1997   gaz
CD2 improved formulation
CD2 
CD2    Rev 1.9   Wed Feb 14 10:19:00 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.8   Wed Jan 10 10:45:16 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.7   Tue Jan 09 15:23:24 1996   hend
CD2 Revised Prolog
CD2 
CD2    Rev 1.6   Tue Jan 09 15:19:18 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.5   11/15/95 09:31:44   gaz
CD2 removed phase change criteria for Richards equation
CD2 
CD2    Rev 1.4   08/18/95 09:48:36   llt
CD2 neq was already defined, removed for cray
CD2 
CD2    Rev 1.3   06/02/95 10:15:22   llt
CD2 added header and gaz changes
CD2
CD2    Rev 1.2   05/01/95 15:19:52   gaz
CD2 variable update and derivative control for 1dof solution (Richards) (gaz)
CD2
CD2    Rev 1.1   03/23/95 19:00:58   gaz
CD2 gaz changed to 1 dof solution
CD2
CD2    Rev 1.0   03/10/95 11:58:10   llt
CD2 new routine
CD2
C**********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.3.3 Noncondensible gas flow equations
CD3
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  Not implemented for dual continuum
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

      use comci
      use comdi
      use comii
      use comdti
      use comai

      implicit none
      integer nelm(*),nmat(*),nrhs(*)
      integer ndummy,irdof,iset 
      integer iadf(21)
      integer neqp1,i,jj,i1,i2,kb,ipos,iflg
      real*8 a(*),bp(*),sx1(*)
      real*8 a1,a2,a3,a4,b1,b2,b3,b4
      real*8 a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16
      data iadf/1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3/
c
      if(iflg.eq.0.and.ico2.lt.0.and.idpdp.eq.0) then
c     
c     modify solution matrix for 1dof solution       
c     
         if(irdof.eq.1) then
            neqp1=neq+1
c     first mark variables for the switch
            do i=1,neq
               if(so(i).lt.1.0) then
                  ieos(i) = 2
               else
                  ieos(i) = 1
               endif
            enddo
            do i=1,neq
c     rearrange the jacobian matrix
               i1=nelm(i)+1
               i2=nelm(i+1)
               iset = 0
               do jj=i1,i2
                  kb=nelm(jj)
                  if(ieos(kb).eq.2) then
                     ipos=jj-neqp1
                     a1=a(ipos+nmat(1))
                     a2=a(ipos+nmat(2))
                     a3=a(ipos+nmat(3))
                     a4=a(ipos+nmat(4))
                     a(ipos+nmat(1))=a2
                     a(ipos+nmat(2))=a1
                     a(ipos+nmat(3))=a4
                     a(ipos+nmat(4))=a3
                  else
                     iset = 1
                  endif
               enddo
               if(iset.eq.0) then
                  do jj=i1,i2
                     kb=nelm(jj)
                     ipos=jj-neqp1
                     a1=a(ipos+nmat(1))
                     a2=a(ipos+nmat(2))
                     a3=a(ipos+nmat(3))
                     a4=a(ipos+nmat(4))
                     a(ipos+nmat(1))=a3
                     a(ipos+nmat(2))=a4
                     a(ipos+nmat(3))=a1
                     a(ipos+nmat(4))=a2
                  enddo
                  b1 = bp(i+nrhs(1))
                  b2 = bp(i+nrhs(2))
                  bp(i+nrhs(1)) = b2
                  bp(i+nrhs(2)) = b1
               endif
            enddo
         endif
      else if(iflg.eq.0.and.ico2.gt.0) then
c     
c     modify solution matrix for 1dof solution       
c     
         neqp1=neq+1
         do i=1,neq
c     modify balance equations
            i1=nelm(i)+1
            i2=nelm(i+1)
            b1=bp(i+nrhs(1))
            b2=bp(i+nrhs(2))
            b3=bp(i+nrhs(3))
            if(iadf(iad).eq.1) then
               do jj=i1,i2
                  kb=nelm(jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat(1))
                  a2=a(ipos+nmat(2))
                  a3=a(ipos+nmat(3))
                  a4=a(ipos+nmat(4))
                  a5=a(ipos+nmat(5))
                  a6=a(ipos+nmat(6))
                  a7=a(ipos+nmat(7))
                  a8=a(ipos+nmat(8))
                  a9=a(ipos+nmat(9))
                  a(ipos+nmat(1))=a1
                  a(ipos+nmat(2))=a2
                  a(ipos+nmat(3))=a3
                  a(ipos+nmat(4))=a4
                  a(ipos+nmat(5))=a5
                  a(ipos+nmat(6))=a6
                  a(ipos+nmat(7))=a7
                  a(ipos+nmat(8))=a8
                  a(ipos+nmat(9))=a9
               enddo
               bp(i+nrhs(1))=b1              
               bp(i+nrhs(2))=b2              
               bp(i+nrhs(3))=b3              
            else if(iadf(iad).eq.2) then
               do jj=i1,i2
                  kb=nelm(jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat(1))
                  a2=a(ipos+nmat(2))
                  a3=a(ipos+nmat(3))
                  a4=a(ipos+nmat(4))
                  a5=a(ipos+nmat(5))
                  a6=a(ipos+nmat(6))
                  a7=a(ipos+nmat(7))
                  a8=a(ipos+nmat(8))
                  a9=a(ipos+nmat(9))
                  a(ipos+nmat(1))=a4
                  a(ipos+nmat(2))=a5
                  a(ipos+nmat(3))=a6
                  a(ipos+nmat(4))=a7
                  a(ipos+nmat(5))=a8
                  a(ipos+nmat(6))=a9
                  a(ipos+nmat(7))=a1
                  a(ipos+nmat(8))=a2
                  a(ipos+nmat(9))=a3
               enddo
               bp(i+nrhs(1))=b2              
               bp(i+nrhs(2))=b3              
               bp(i+nrhs(3))=b1              
            else if(iadf(iad).eq.3) then
               do jj=i1,i2
                  kb=nelm(jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat(1))
                  a2=a(ipos+nmat(2))
                  a3=a(ipos+nmat(3))
                  a4=a(ipos+nmat(4))
                  a5=a(ipos+nmat(5))
                  a6=a(ipos+nmat(6))
                  a7=a(ipos+nmat(7))
                  a8=a(ipos+nmat(8))
                  a9=a(ipos+nmat(9))
                  a(ipos+nmat(1))=a7
                  a(ipos+nmat(2))=a8
                  a(ipos+nmat(3))=a9
                  a(ipos+nmat(4))=a1
                  a(ipos+nmat(5))=a2
                  a(ipos+nmat(6))=a3
                  a(ipos+nmat(7))=a4
                  a(ipos+nmat(8))=a5
                  a(ipos+nmat(9))=a6
               enddo
               bp(i+nrhs(1))=b3
               bp(i+nrhs(2))=b1
               bp(i+nrhs(3))=b2
            endif
         enddo
      else if(iflg.eq.0.and.ico2.lt.0.and.
     &        idpdp.ne.0.and.irdof.ne.0) then
c     
c     modify solution matrix for 4dof solution       
c     
         neqp1=neq+1
         do i=1,neq
c     modify balance equations
            i1=nelm(i)+1
            i2=nelm(i+1)
            b1=bp(i+nrhs(1))
            b2=bp(i+nrhs(2))
            b3=bp(i+nrhs(3))
            b4=bp(i+nrhs(4))
            if(s(i).gt.0.0.and.s(i).lt.1.0.and.
     &           s(i+neq).gt.0.0.and.s(i+neq).lt.1.0) then
               do jj=i1,i2
                  kb=nelm(jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat(1))
                  a2=a(ipos+nmat(2))
                  a3=a(ipos+nmat(3))
                  a4=a(ipos+nmat(4))
                  a5=a(ipos+nmat(5))
                  a6=a(ipos+nmat(6))
                  a7=a(ipos+nmat(7))
                  a8=a(ipos+nmat(8))
                  a9=a(ipos+nmat(9))
                  a10=a(ipos+nmat(10))
                  a11=a(ipos+nmat(11))
                  a12=a(ipos+nmat(12))
                  a13=a(ipos+nmat(13))
                  a14=a(ipos+nmat(14))
                  a15=a(ipos+nmat(15))
                  a16=a(ipos+nmat(16))
                  a(ipos+nmat(1))=a2
                  a(ipos+nmat(2))=a4
                  a(ipos+nmat(3))=a1
                  a(ipos+nmat(4))=a3
                  a(ipos+nmat(5))=a10
                  a(ipos+nmat(6))=a12
                  a(ipos+nmat(7))=a9
                  a(ipos+nmat(8))=a11
                  a(ipos+nmat(9))=a6
                  a(ipos+nmat(10))=a8
                  a(ipos+nmat(11))=a5
                  a(ipos+nmat(12))=a7
                  a(ipos+nmat(13))=a14
                  a(ipos+nmat(14))=a16
                  a(ipos+nmat(15))=a13
                  a(ipos+nmat(16))=a15
               enddo
            endif 
            bp(i+nrhs(1))=b1              
            bp(i+nrhs(2))=b3              
            bp(i+nrhs(3))=b2 
            bp(i+nrhs(4))=b4             
         enddo
      else if(iflg.eq.1.and.ico2.lt.0.and.idpdp.eq.0) then
c     
c     extract solution for water and air
c     
         if(irdof.eq.1) then
            do i=1,neq

               if(ieos(i).eq.2) then
                  b1=bp(i+nrhs(1))
                  b2=bp(i+nrhs(2))
                  bp(i+nrhs(1))=b2              
                  bp(i+nrhs(2))=b1              
               endif
            enddo
         endif
      else if(iflg.eq.1.and.ico2.lt.0.and.
     &        idpdp.ne.0.and.irdof.ne.0) then
c     
c     extract solution for water-air dual permeability
c     
         do i=1,neq
            if(s(i).gt.0.0.and.s(i).lt.1.0.and.
     &           s(i+neq).gt.0.0.and.s(i+neq).lt.1.0) then
               b1=bp(i+nrhs(1))
               b2=bp(i+nrhs(2))
               b3=bp(i+nrhs(3))
               b4=bp(i+nrhs(4))
               bp(i+nrhs(1))=b3              
               bp(i+nrhs(2))=b1   
               bp(i+nrhs(3))=b4 
               bp(i+nrhs(4))=b2
            endif   
         enddo
      else if(iflg.eq.1.and.ico2.gt.0) then
c     
c     extract solution for heat-mass-air
c     
         do i=1,neq
            b1=bp(i+nrhs(1))
            b2=bp(i+nrhs(2))
            b3=bp(i+nrhs(3))
            bp(i+nrhs(1))=b1              
            bp(i+nrhs(2))=b2   
            bp(i+nrhs(3))=b3   
         enddo
      endif 
      return
      end     

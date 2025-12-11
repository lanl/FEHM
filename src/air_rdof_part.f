      subroutine air_rdof_part(iflg,irdof,ndummy,nelm,nmat,nrhs,a,
     &     bp,sx1,zone)
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
!***********************************************************************
!D1
!D1  PURPOSE
!D1
!D1  This subroutine reduces the degree of freedom solution
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY
!D2 
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/air_rdof_part.f_a  $
!D2 
!**********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  2.3.2 Heat- and mass-transfer equations
!D3  2.3.3 Noncondensible gas flow equations
!D3
!**********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4  Not implemented for dual continuum
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comci
      use comdi
      use comii
      use comdti
      use comai
      use com_part

      implicit none
      integer zone
      integer nelm(*),nmat(*),nrhs(*)
      integer ndummy,irdof 
      integer iadf(21)
      integer neqp1,i,jj,i1,i2,kb,ipos,iflg
      real*8 a(*),bp(*),sx1(*)
      real*8 a1,a2,a3,a4,b1,b2,b3,b4 
      real*8 a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16
      data iadf/1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3/
c
      if(iflg.eq.0.and.ico2.lt.0.and.idpdp.eq.0) then
c       
c modify solution matrix for 1dof solution       
c       
         if(irdof.eq.1) then
            neqp1=neq_part(zone)+1
            do i=1,neq_part(zone)
c rearrange the jacobian matrix
c         if(so(i).gt.0.0d00.and.so(i).lt.1.d00) then
               if(so(index_part(zone,i)).gt.0.0d00.and.
     &            so(index_part(zone,i)).le.2.d00) then
                  i1=nelm_part(zone,i)+1
                  i2=nelm_part(zone,i+1)
                  do jj=i1,i2
                     kb=nelm_part(zone,jj)
                     ipos=jj-neqp1
                     a1=a(ipos+nmat_part(zone,1))
                     a2=a(ipos+nmat_part(zone,2))
                     a3=a(ipos+nmat_part(zone,3))
                     a4=a(ipos+nmat_part(zone,4))
                     a(ipos+nmat_part(zone,1))=a2
                     a(ipos+nmat_part(zone,2))=a1
                     a(ipos+nmat_part(zone,3))=a4
                     a(ipos+nmat_part(zone,4))=a3
                  enddo
               else if(so(index_part(zone,i)).le.0.0d00) then
                  i1=nelm_part(zone,i)+1
                  i2=nelm_part(zone,i+1)
                  do jj=i1,i2
                     kb=nelm_part(zone,jj)
                     ipos=jj-neqp1
                     a1=a(ipos+nmat_part(zone,1))
                     a2=a(ipos+nmat_part(zone,2))
                     a3=a(ipos+nmat_part(zone,3))
                     a4=a(ipos+nmat_part(zone,4))
                     a(ipos+nmat_part(zone,1))=a3
                     a(ipos+nmat_part(zone,2))=a4
                     a(ipos+nmat_part(zone,3))=a1
                     a(ipos+nmat_part(zone,4))=a2
                  enddo
                  b1=bp_part(zone,i+nrhs_part(zone,1))
                  b2=bp_part(zone,i+nrhs_part(zone,2))
                  bp_part(zone,i+nrhs_part(zone,1)) = b2
                  bp_part(zone,i+nrhs_part(zone,2)) = b1
               endif
            enddo
         endif
      else if(iflg.eq.0.and.ico2.gt.0) then
c       
c modify solution matrix for 1dof solution       
c       
         neqp1=neq_part(zone)+1
         do i=1,neq_part(zone)
c modify balance equations
            i1=nelm_part(zone,i)+1
            i2=nelm_part(zone,i+1)
            b1=bp_part(zone,i+nrhs_part(zone,1))
            b2=bp_part(zone,i+nrhs_part(zone,2))
            b3=bp_part(zone,i+nrhs_part(zone,3))
            if(iadf(iad).eq.1) then
               do jj=i1,i2
                  kb=nelm_part(zone,jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat_part(zone,1))
                  a2=a(ipos+nmat_part(zone,2))
                  a3=a(ipos+nmat_part(zone,3))
                  a4=a(ipos+nmat_part(zone,4))
                  a5=a(ipos+nmat_part(zone,5))
                  a6=a(ipos+nmat_part(zone,6))
                  a7=a(ipos+nmat_part(zone,7))
                  a8=a(ipos+nmat_part(zone,8))
                  a9=a(ipos+nmat_part(zone,9))
                  a(ipos+nmat_part(zone,1))=a1
                  a(ipos+nmat_part(zone,2))=a2
                  a(ipos+nmat_part(zone,3))=a3
                  a(ipos+nmat_part(zone,4))=a4
                  a(ipos+nmat_part(zone,5))=a5
                  a(ipos+nmat_part(zone,6))=a6
                  a(ipos+nmat_part(zone,7))=a7
                  a(ipos+nmat_part(zone,8))=a8
                  a(ipos+nmat_part(zone,9))=a9
               enddo
               bp_part(zone,i+nrhs_part(zone,1))=b1              
               bp_part(zone,i+nrhs_part(zone,2))=b2              
               bp_part(zone,i+nrhs_part(zone,3))=b3              
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
                  kb=nelm_part(zone,jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat_part(zone,1))
                  a2=a(ipos+nmat_part(zone,2))
                  a3=a(ipos+nmat_part(zone,3))
                  a4=a(ipos+nmat_part(zone,4))
                  a5=a(ipos+nmat_part(zone,5))
                  a6=a(ipos+nmat_part(zone,6))
                  a7=a(ipos+nmat_part(zone,7))
                  a8=a(ipos+nmat_part(zone,8))
                  a9=a(ipos+nmat_part(zone,9))
                  a(ipos+nmat_part(zone,1))=a7
                  a(ipos+nmat_part(zone,2))=a8
                  a(ipos+nmat_part(zone,3))=a9
                  a(ipos+nmat_part(zone,4))=a1
                  a(ipos+nmat_part(zone,5))=a2
                  a(ipos+nmat_part(zone,6))=a3
                  a(ipos+nmat_part(zone,7))=a4
                  a(ipos+nmat_part(zone,8))=a5
                  a(ipos+nmat_part(zone,9))=a6
               enddo
               bp_part(zone,i+nrhs_part(zone,1))=b3
               bp_part(zone,i+nrhs_part(zone,2))=b1
               bp_part(zone,i+nrhs_part(zone,3))=b2
            endif
         enddo
      else if(iflg.eq.0.and.ico2.lt.0.and.
     &        idpdp.ne.0.and.irdof.ne.0) then
c       
c modify solution matrix for 4dof solution       
c       
         neqp1=neq_part(zone)+1
         do i=1,neq_part(zone)
c modify balance equations
            i1=nelm_part(zone,i)+1
            i2=nelm_part(zone,i+1)
            b1=bp_part(zone,i+nrhs_part(zone,1))
            b2=bp_part(zone,i+nrhs_part(zone,2))
            b3=bp_part(zone,i+nrhs_part(zone,3))
            if(s(index_part(zone,i)).gt.0.0.and.
     &         s(index_part(zone,i)).lt.1.0.and.
     &         s(index_part(zone,i+neq_part(zone))).gt.0.0.and.
     &         s(index_part(zone,i+neq_part(zone))).lt.1.0) then
               do jj=i1,i2
                  kb=nelm_part(zone,jj)
                  ipos=jj-neqp1
                  a1=a(ipos+nmat_part(zone,1))
                  a2=a(ipos+nmat_part(zone,2))
                  a3=a(ipos+nmat_part(zone,3))
                  a4=a(ipos+nmat_part(zone,4))
                  a5=a(ipos+nmat_part(zone,5))
                  a6=a(ipos+nmat_part(zone,6))
                  a7=a(ipos+nmat_part(zone,7))
                  a8=a(ipos+nmat_part(zone,8))
                  a9=a(ipos+nmat_part(zone,9))
                  a10=a(ipos+nmat_part(zone,10))
                  a11=a(ipos+nmat_part(zone,11))
                  a12=a(ipos+nmat_part(zone,12))
                  a13=a(ipos+nmat_part(zone,13))
                  a14=a(ipos+nmat_part(zone,14))
                  a15=a(ipos+nmat_part(zone,15))
                  a16=a(ipos+nmat_part(zone,16))
                  a(ipos+nmat_part(zone,1))=a2
                  a(ipos+nmat_part(zone,2))=a4
                  a(ipos+nmat_part(zone,3))=a1
                  a(ipos+nmat_part(zone,4))=a3
                  a(ipos+nmat_part(zone,5))=a10
                  a(ipos+nmat_part(zone,6))=a12
                  a(ipos+nmat_part(zone,7))=a9
                  a(ipos+nmat_part(zone,8))=a11
                  a(ipos+nmat_part(zone,9))=a6
                  a(ipos+nmat_part(zone,10))=a8
                  a(ipos+nmat_part(zone,11))=a5
                  a(ipos+nmat_part(zone,12))=a7
                  a(ipos+nmat_part(zone,13))=a14
                  a(ipos+nmat_part(zone,14))=a16
                  a(ipos+nmat_part(zone,15))=a13
                  a(ipos+nmat_part(zone,16))=a15
               enddo
            endif 
            bp_part(zone,i+nrhs_part(zone,1))=b1              
            bp_part(zone,i+nrhs_part(zone,2))=b3              
            bp_part(zone,i+nrhs_part(zone,3))=b2 
            bp_part(zone,i+nrhs_part(zone,4))=b4             
         enddo
      else if(iflg.eq.1.and.ico2.lt.0.and.idpdp.eq.0) then
c
c extract solution for water and air
c
         if(irdof.eq.1) then
            do i=1,neq_part(zone)
c        if(so(i)).gt.0.0d00.and.so(i).lt.1.d00) then
               if(so(index_part(zone,i)).gt.0.0d00.and.
     &            so(index_part(zone,i)).le.2.d00) then
                  b1=bp_part(zone,i+nrhs_part(zone,1))
                  b2=bp_part(zone,i+nrhs_part(zone,2))
                  bp_part(zone,i+nrhs_part(zone,1))=b2              
                  bp_part(zone,i+nrhs_part(zone,2))=b1              
               endif
            enddo
         endif
      else if(iflg.eq.1.and.ico2.lt.0.and.
     &        idpdp.ne.0.and.irdof.ne.0) then
c
c extract solution for water-air dual permeability
c
         do i=1,neq_part(zone)
            if(s(index_part(zone,i)).gt.0.0.and.
     &         s(index_part(zone,i)).lt.1.0.and.
     &         s(index_part(zone,i+neq_part(zone))).gt.0.0.and.
     &         s(index_part(zone,i+neq_part(zone))).lt.1.0) then
               b1=bp_part(zone,i+nrhs_part(zone,1))
               b2=bp_part(zone,i+nrhs_part(zone,2))
               b3=bp_part(zone,i+nrhs_part(zone,3))
               b4=bp_part(zone,i+nrhs_part(zone,4))
               bp_part(zone,i+nrhs_part(zone,1))=b3              
               bp_part(zone,i+nrhs_part(zone,2))=b1   
               bp_part(zone,i+nrhs_part(zone,3))=b4 
               bp_part(zone,i+nrhs_part(zone,4))=b2
            endif   
         enddo
      else if(iflg.eq.1.and.ico2.gt.0) then
c
c extract solution for heat-mass-air
c
         do i=1,neq_part(zone)
            b1=bp_part(zone,i+nrhs_part(zone,1))
            b2=bp_part(zone,i+nrhs_part(zone,2))
            b3=bp_part(zone,i+nrhs_part(zone,3))
            bp_part(zone,i+nrhs_part(zone,1))=b1              
            bp_part(zone,i+nrhs_part(zone,2))=b2   
            bp_part(zone,i+nrhs_part(zone,3))=b3   
         enddo
      endif 
      return
      end     


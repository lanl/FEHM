      subroutine explicit(nactive) 
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
!D1 PURPOSE
!D1 
!D1 To make the equations explicit (i.e., if timestep fails we will 
!D1 take an explicit time step)
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/explicit.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3  2.5.1 Implement time-step mechanism
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comji
      use davidi
      use comflow
      implicit none
c
      integer neqp1,i,j,jmia,nactive
      real*8 sx1d,tmche,tolerance
      real*8 bp1max,bp2max
c
      if(sssol(1:1).ne.'e') then
         return
      else if(sssol(2:2).eq.'1') then
         tolerance=1.d-1
      else if(sssol(2:2).eq.'2') then
         tolerance=1.d-2
      else if(sssol(2:2).eq.'3') then
         tolerance=1.d-3
      else if(sssol(2:2).eq.'4') then
         tolerance=1.d-4
      else if(sssol(2:2).eq.'5') then
         tolerance=1.d-5
      else if(sssol(2:2).eq.'6') then
         tolerance=1.d-6
      else if(sssol(2:2).eq.'7') then
         tolerance=1.d-7
      else if(sssol(2:2).eq.'8') then
         tolerance=1.d-8
      else if(sssol(2:2).eq.'9') then
         tolerance=1.d-9
      else if(sssol(2:2).eq.'0') then
         tolerance=1.d-10
      endif
      neqp1=neq+1
c
c only implemented for 2 dof problems
c       bp1max=0.0
c       bp2max=0.0
c
      do i=1,neq
         nopt(i)=1
c        bp1max=max(bp1max,bp(i+nrhs(1)))
c        bp2max=max(bp2max,bp(i+nrhs(2)))
      enddo
c       tmche=min(bp1max,bp2max,tmch)*tolerance
      tmche=tmch*tolerance
c check if any active nodes exist
 100  continue
      nactive=0
      do i=1,neq
         if(abs(bp(i+nrhs(1))).gt.tmche.or.
     &        abs(bp(i+nrhs(2))).gt.tmche) then
            nactive = nactive +1
         endif
      enddo
c if noactive nodes exist tighten tolerance
      if(nactive.eq.0) then
         tmche=tmche*0.1
         go to 100 
      endif
      do i=1,neq
         if(abs(bp(i+nrhs(1))).gt.tmche.or.
     &        abs(bp(i+nrhs(2))).gt.tmche) then
            do j=nelm(i)+1,nelm(i+1)
               nopt(nelm(j)) = 0
            enddo
         endif
      enddo
      nactive=0
      do i=1,neq
         if(nopt(i).eq.1) then
            sx1d=sx1(i)
            jmia=nelmdg(i)-neqp1
            do j=nelm(i)+1,nelm(i+1)
               a(j+nmat(1)-neqp1)=0.0                    
               a(j+nmat(2)-neqp1)=0.0                    
               a(j+nmat(3)-neqp1)=0.0                    
               a(j+nmat(4)-neqp1)=0.0                    
            enddo
            a(jmia+nmat(1))=1.0 
            a(jmia+nmat(2))=0.0
            a(jmia+nmat(3))=0.0
            a(jmia+nmat(4))=1.0 
            bp(i+nrhs(1))= 0.0
            bp(i+nrhs(2))= 0.0
         else
            nactive = nactive +1
         endif
      enddo
      return
      end

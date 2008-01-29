      subroutine varctr(iflg)
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
CD1  This subroutine manages the head input and changes pressures
CD1  into heads and vice versa           
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation: 6-FEB-97, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/varctr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:02 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS                  
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
**********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comrxni
      use comii
      use davidi
      implicit none


c
      integer iflg,icode
      integer mi,neqp1,i,i1,i2,j,jj,ja 
      real*8 a1,a2,b1,b2,aa1,aa2,aa3,aa4,bp1,bp2
c================================================================
      if(ivar.eq.0) return
c================================================================
      if(iflg.eq.0) then
         allocate(var1(neq),var2(neq),dvar11(neq),dvar12(neq))
         allocate(dvar21(neq),dvar22(neq))
c       call mmgetblk ("var1","comdi",ipvar1,neq,2, icode)
c       call mmgetblk ("var2","comdi",ipvar2,neq,2, icode)
c       call mmgetblk ("dvar11","comdi",ipdvar11,neq,2, icode)
c       call mmgetblk ("dvar12","comdi",ipdvar12,neq,2, icode)
c       call mmgetblk ("dvar21","comdi",ipdvar21,neq,2, icode)
c       call mmgetblk ("dvar22","comdi",ipdvar22,neq,2, icode)
c
c read input       
c ivar=1 : pressure and enthalpy are variables
c
       read(inpt,*) ivar
       if(ivar.ne.0) then
        ivar=1
       endif

      else if(iflg.eq.1) then
c
c modify initial values and BC's
c

      else if(iflg.eq.2) then
c
c convert equations to pressure enthalpy     
c
c 
c modify NR equations
c
c test equations of the form 
c v1=a1*phi + b1*t
c v2=a2*phi + b2*t
c a1=0.1, b1=10., a2=1.0, b2=100. 

       a1=0.1
       b1=10.0
       a2=1.0
       b2=100.
       neqp1=neq+1
       do i=1,neq
        i1=nelm(i)+1
        i2=nelm(i+1)
        do jj=i1,i2
         j=nelm(jj)
         ja=jj-neqp1
          aa1=a(nmat(1)+ja)
          aa2=a(nmat(2)+ja)
          aa3=a(nmat(3)+ja)
          aa4=a(nmat(4)+ja)
          a(nmat(1)+ja)=aa1/a1+aa2/b1       
          a(nmat(2)+ja)=aa1/a2+aa2/b2       
          a(nmat(3)+ja)=aa3/a1+aa4/b1       
          a(nmat(3)+ja)=aa3/a2+aa4/b2       
        enddo
       enddo

      else if(iflg.eq.3) then
c
c convert corrections to pressure temperature  
c
          do mi=1,neq
           bp1 = bp(mi+nrhs(1))
           bp2 = bp(mi+nrhs(2))
c          bp(mi+nrhs(1))= 
          enddo

      endif
c 
      return
      end                

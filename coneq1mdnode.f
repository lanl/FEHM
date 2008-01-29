      subroutine coneq1mdnode(matnum, spec_num)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To compute the Jacobian and residual terms of the concentration
CD1 equation associated with multiply defined node connections.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/coneq1mdnode.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:48   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.1   Thu Jun 06 15:11:02 1996   robinson
CD2 Fixed implementation for dpdp
CD2 
CD2    Rev 1.0   Fri May 24 09:53:20 1996   hend
CD2 Initial Implementation
CD2 
C**********************************************************************
CD3
CD3 SPECIAL COMMENTS AND REFERENCES
CD3 
CD3  Requirements from SDN: 10086-RD-2.20-00
CD3    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD3    FEHM Application Version 2.20
CD3
C**********************************************************************
CD4
CD4 REQUIREMENTS TRACEABILITY
CD4 
CD4 2.3.4 Solute-transport equations
CD4 
C**********************************************************************

      use comflow
      use comcouple
      use comrxni
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer i,j,ii,itr,icd,ii1,ii2,idg,iq,jmi,jml,jmia,jm,neqp1,ij
      integer ij1,ij2,iz,kb,iau,ial,kz
      real*8 anli,anvi,danli,danvi,anlri,anvri,danlri
      real*8 axy,dlaei,danvri,dlaekb,vxy,dvaei,dvaekb,toldil
      integer matnum,spec_num
      parameter(toldil = 1.d-20)

c NOTE---IMPLEMENTED FOR DPDP!!!!!! but only one species
c nmat_sol(1)=nmatb(4) and nrhs_sol(1)=nrhs(2)
      do i=1,n
         if(i.gt.neq.and.idualp.eq.0) then
            icd=neq
            nmat_sol(1)=nmatb(4)
            nrhs_sol(1)=nrhs(2)
         else
            icd=0
            nmat_sol(1)=0
            nrhs_sol(1)=0
         endif
         iz=i-icd
         do j=1,abs(mdnodes(iz))
            ii=mdnode(iz,j)
            if (ii.gt.iz) then

               itr=i+npn
               anli=anl(itr)
               anvi=anv(itr)
               danli=danl(itr)
               danvi=danv(itr)
               anlri=anli*rolf(i)
               anvri=anvi*rovf(i)
               danlri=danli*rolf(i)
               danvri=danvi*rovf(i)

               ii1=nelm(iz)+1
               ii2=nelm(iz+1)
               idg=nelmdg(iz)-ii1+1
               neqp1=neq+1
               jmi=nelmdg(iz)
               jmia=jmi-neqp1

               iq=1
               do jm=jmi+1,ii2
                  if(nelm(jm).eq.ii) then
                     it8(iq)=nelm(jm)+icd
                     it9(iq)=jm-ii1+1
                     it10(iq)=istrw(jm-neqp1)
                     it11(iq)=jm-neqp1
                     ij1=nelm(nelm(jm))+1
                     ij2=nelmdg(nelm(jm))-1
                     do  ij=ij1,ij2
                        if(nelm(ij).eq.iz) then
                           it12(iq)=ij-neqp1
                        endif
                     enddo
                  endif
               enddo

               if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2) then
c     liquid phase calculations
c     add dispersion term for liquid 
                  jm=1
                  kb=it8(jm)
                  kz=kb-icd
                  iau=it11(jm)
                  ial=it12(jm)
                  jml=nelmdg(kz)-neqp1
                  dlaei=max(abs(a(jmia+nmat_sol(matnum))),
     &                 abs(a(ial+nmat_sol(matnum))))*iyear*1e6
                  dlaekb=-(dlaei/danlri)*rolf(kb)*danl(kb+npn)
                  axy=-(dlaei/danlri)*(anl(kb+npn)*rolf(kb)-anlri)
                  bp(iz+nrhs_sol(spec_num))=
     &                 bp(iz+nrhs_sol(spec_num))+axy
                  bp(kz+nrhs_sol(spec_num))=
     &                 bp(kz+nrhs_sol(spec_num))-axy
                  a(jmia+nmat_sol(matnum))=
     &                 a(jmia+nmat_sol(matnum))+dlaei
                  a(ial+nmat_sol(matnum))=
     &                 a(ial+nmat_sol(matnum))-dlaei
                  a(iau+nmat_sol(matnum))=
     &                 a(iau+nmat_sol(matnum))+dlaekb
                  a(jml+nmat_sol(matnum))=
     &                 a(jml+nmat_sol(matnum))-dlaekb
               endif

               if (icns(nsp).eq.-1.or.abs(icns(nsp)).eq.2) then
c     vapor phase calculations
c     add dispersion term for vapor 
                  jm=1
                  kb=it8(jm)
                  kz=kb-icd
                  iau=it11(jm)
                  ial=it12(jm)
                  jml=nelmdg(kz)-neqp1
                  dvaei=max(abs(a(jmia+nmat_sol(matnum))),
     &                 abs(a(ial+nmat_sol(matnum))))*iyear*1e6
                  dvaekb=-(dvaei/danvri)*rovf(kb)*danv(kb+npn)
                  vxy=-(dvaei/danvri)*(anv(kb+npn)*rovf(kb)-anvri)
                  bp(iz+nrhs_sol(spec_num))=
     &                 bp(iz+nrhs_sol(spec_num))+vxy
                  bp(kz+nrhs_sol(spec_num))=
     &                 bp(kz+nrhs_sol(spec_num))-vxy
                  a(jmia+nmat_sol(matnum))=
     &                 a(jmia+nmat_sol(matnum))+dvaei
                  a(ial+nmat_sol(matnum))=
     &                 a(ial+nmat_sol(matnum))-dvaei
                  a(iau+nmat_sol(matnum))=
     &                 a(iau+nmat_sol(matnum))+dvaekb
                  a(jml+nmat_sol(matnum))=
     &                 a(jml+nmat_sol(matnum))-dvaekb
               endif
            endif
         enddo
      enddo

      return
      end
      
      

      subroutine setord
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
CD1  This subroutine sets up the order of the solution for the 
CD1  equations at each node.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/setord.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed May 29 08:19:46 1996   gaz
CD2 corrected re-number for 6dof
CD2 
CD2    Rev 1.5   Thu Jan 11 11:00:06 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   Tue Jan 09 14:10:50 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.3   11/15/95 15:21:30   gaz
CD2 added 4dof reordering
CD2 
CD2    Rev 1.2   01/28/95 13:55:34   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 15:59:52   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:28   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
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

      use comci
      use combi
      use davidi
      use comrxni
      use comdti
      use comcouple
      use comai
      implicit none

      integer i, idofdb, neqp1, nsizea, nsizeb, mdof
      integer neqp1_primary

c load nb and nmat matrix
      neqp1=neq+1
      nsizea=nelm(neqp1)-neqp1
c changed size of nb
      if(gdpm_flag.ne.0) then
         neqp1_primary=neq_primary+1
         if (igauss .gt. 1) then
            nsizeb=nop(neqp1_primary)-neqp1_primary
         else
            nsizeb = 1
         end if
      else
         if (igauss .gt. 1) then
            nsizeb=nop(neqp1)-neqp1
         else
            nsizeb = 1
         end if      
      endif
      nb(1)=0
      nmat(1)=0
      do i=2,36
         nb(i)=nb(i-1)+nsizeb
         nmat(i)=nmat(i-1)+nsizea
      enddo
      do i=1,6
         nrhs(i)=(i-1)*neq
      enddo
      idofdb=idof
      if(idofdb.le.1) then
         nmat(1)=0
         nrhs(1)=0
         nrhs(2)=neq
      else if(idofdb.eq.2) then
         if(islord.eq.0) then
            nmatb(1)=0
            nmatb(2)=nsizea
            nmatb(3)=nmatb(2)+nsizea
            nmatb(4)=nmatb(3)+nsizea
            nrhsb(1)=0
            nrhsb(2)=neq
            nb(1)=0
            nb(2)=nsizeb
            nb(3)=nb(2)+nsizeb
            nb(4)=nb(3)+nsizeb
         else if(islord.eq.1) then
            nmatb(2)=0
            nmatb(1)=nsizea
            nmatb(4)=nmatb(1)+nsizea
            nmatb(3)=nmatb(4)+nsizea
            nrhsb(2)=0
            nrhsb(1)=neq
            nb(2)=0
            nb(1)=nsizeb
            nb(4)=nb(1)+nsizeb
            nb(3)=nb(4)+nsizeb
         else if(islord.eq.-1) then
            nmatb(4)=0
            nmatb(3)=nsizea
            nmatb(2)=nmatb(3)+nsizea
            nmatb(1)=nmatb(2)+nsizea
            nrhsb(2)=0
            nrhsb(1)=neq
            nb(4)=0
            nb(3)=nsizeb
            nb(2)=nb(3)+nsizeb
            nb(1)=nb(2)+nsizeb
         endif
      else if(idofdb.eq.3) then
         if(islord.eq.0) then
            nmatb(1)=0
            nmatb(2)=nsizea
            nmatb(3)=nmatb(2)+nsizea
            nmatb(4)=nmatb(3)+nsizea
            nmatb(5)=nmatb(4)+nsizea
            nmatb(6)=nmatb(5)+nsizea
            nmatb(7)=nmatb(6)+nsizea
            nmatb(8)=nmatb(7)+nsizea
            nmatb(9)=nmatb(8)+nsizea
            nrhsb(1)=0
            nrhsb(2)=neq     
            nrhsb(3)=2*neq  
         else if(islord.eq.1) then
            nmatb(1)=0
            nmatb(3)=nsizea
            nmatb(2)=2*nsizea
            nmatb(4)=3*nsizea
            nmatb(6)=4*nsizea
            nmatb(5)=5*nsizea
            nmatb(7)=6*nsizea
            nmatb(9)=7*nsizea
            nmatb(8)=8*nsizea
            nrhsb(1)=0
            nrhsb(2)=2*neq   
            nrhsb(3)=neq    
         else if(islord.eq.2) then
            nmatb(2)=0
            nmatb(1)=nsizea
            nmatb(3)=2*nsizea
            nmatb(5)=3*nsizea
            nmatb(4)=4*nsizea
            nmatb(6)=5*nsizea
            nmatb(8)=6*nsizea
            nmatb(7)=7*nsizea
            nmatb(9)=8*nsizea
            nrhsb(1)=neq
            nrhsb(2)=0 
            nrhsb(3)=2*neq    
         else if(islord.eq.3) then
            nmatb(3)=0
            nmatb(2)=nsizea
            nmatb(1)=2*nsizea
            nmatb(6)=3*nsizea
            nmatb(5)=4*nsizea
            nmatb(4)=5*nsizea
            nmatb(9)=6*nsizea
            nmatb(8)=7*nsizea
            nmatb(7)=8*nsizea
            nrhsb(3)=0
            nrhsb(2)=neq
            nrhsb(1)=2*neq               
         endif
      else if(idofdb.eq.4) then
         if(islord.eq.0) then
            nmatb(1)=0
            nmatb(2)=nsizea
            nmatb(3)=nmatb(2)+nsizea
            nmatb(4)=nmatb(3)+nsizea
            nmatb(5)=nmatb(4)+nsizea
            nmatb(6)=nmatb(5)+nsizea
            nmatb(7)=nmatb(6)+nsizea
            nmatb(8)=nmatb(7)+nsizea
            nmatb(9)=nmatb(8)+nsizea
            nmatb(10)=nmatb(9)+nsizea
            nmatb(11)=nmatb(10)+nsizea
            nmatb(12)=nmatb(11)+nsizea
            nmatb(13)=nmatb(12)+nsizea
            nmatb(14)=nmatb(13)+nsizea
            nmatb(15)=nmatb(14)+nsizea
            nmatb(16)=nmatb(15)+nsizea
            nrhsb(1)=0
            nrhsb(2)=neq     
            nrhsb(3)=2*neq  
            nrhsb(4)=3*neq  
         else if(islord.eq.1) then
            nmatb(1)=0
            nmatb(02)=(03-1)*nsizea
            nmatb(03)=(02-1)*nsizea
            nmatb(04)=(04-1)*nsizea
            nmatb(05)=(09-1)*nsizea
            nmatb(06)=(11-1)*nsizea
            nmatb(07)=(10-1)*nsizea
            nmatb(08)=(12-1)*nsizea
            nmatb(09)=(05-1)*nsizea
            nmatb(10)=(07-1)*nsizea
            nmatb(11)=(06-1)*nsizea
            nmatb(12)=(08-1)*nsizea
            nmatb(13)=(13-1)*nsizea
            nmatb(14)=(15-1)*nsizea
            nmatb(15)=(14-1)*nsizea
            nmatb(16)=(16-1)*nsizea
            nrhsb(1)=0
            nrhsb(2)=(3-1)*neq   
            nrhsb(3)=(2-1)*neq    
            nrhsb(4)=(4-1)*neq    
         else if(islord.eq.2) then
            nmatb(01)=(02-1)*nsizea
            nmatb(02)=(04-1)*nsizea
            nmatb(03)=(01-1)*nsizea
            nmatb(04)=(03-1)*nsizea
            nmatb(05)=(10-1)*nsizea
            nmatb(06)=(12-1)*nsizea
            nmatb(07)=(09-1)*nsizea
            nmatb(08)=(11-1)*nsizea
            nmatb(09)=(06-1)*nsizea
            nmatb(10)=(08-1)*nsizea
            nmatb(11)=(05-1)*nsizea
            nmatb(12)=(07-1)*nsizea
            nmatb(13)=(14-1)*nsizea
            nmatb(14)=(16-1)*nsizea
            nmatb(15)=(13-1)*nsizea
            nmatb(16)=(15-1)*nsizea
            nrhsb(1)=0
            nrhsb(2)=(3-1)*neq   
            nrhsb(3)=(2-1)*neq    
            nrhsb(4)=(4-1)*neq    
         endif
      else if(idofdb.eq.6) then
         if(islord.ne.1) then
            do i=1,36
               nmatb(i)=nsizea*(i-1)
            enddo
            nrhsb(1)=0
            nrhsb(2)=neq
            nrhsb(3)=2*neq
            nrhsb(4)=3*neq
            nrhsb(5)=4*neq
            nrhsb(6)=5*neq
         else
            nmatb(1)=0
            nmatb(02)=(04-1)*nsizea
            nmatb(03)=(02-1)*nsizea
            nmatb(04)=(05-1)*nsizea
            nmatb(05)=(03-1)*nsizea
            nmatb(06)=(06-1)*nsizea
            nmatb(07)=(19-1)*nsizea
            nmatb(08)=(22-1)*nsizea
            nmatb(09)=(20-1)*nsizea
            nmatb(10)=(23-1)*nsizea
            nmatb(11)=(21-1)*nsizea
            nmatb(12)=(24-1)*nsizea
            nmatb(13)=(07-1)*nsizea
            nmatb(14)=(10-1)*nsizea
            nmatb(15)=(08-1)*nsizea
            nmatb(16)=(11-1)*nsizea
            nmatb(17)=(09-1)*nsizea
            nmatb(18)=(12-1)*nsizea
            nmatb(19)=(25-1)*nsizea
            nmatb(20)=(28-1)*nsizea
            nmatb(21)=(26-1)*nsizea
            nmatb(22)=(29-1)*nsizea
            nmatb(23)=(27-1)*nsizea
            nmatb(24)=(30-1)*nsizea
            nmatb(25)=(13-1)*nsizea
            nmatb(26)=(16-1)*nsizea
            nmatb(27)=(14-1)*nsizea
            nmatb(28)=(17-1)*nsizea
            nmatb(29)=(15-1)*nsizea
            nmatb(30)=(18-1)*nsizea
            nmatb(31)=(31-1)*nsizea
            nmatb(32)=(34-1)*nsizea
            nmatb(33)=(32-1)*nsizea
            nmatb(34)=(35-1)*nsizea
            nmatb(35)=(33-1)*nsizea
            nmatb(36)=(36-1)*nsizea
            nrhsb(1)=0
            nrhsb(2)=(4-1)*neq
            nrhsb(3)=(2-1)*neq
            nrhsb(4)=(5-1)*neq
            nrhsb(5)=(3-1)*neq
            nrhsb(6)=(6-1)*neq
         endif
      endif

      nb_sol(1)=0
      nmat_sol(1)=0
      do i = 2, mdof_sol**2
         nb_sol(i)=nb_sol(i-1)+nsizeb
         nmat_sol(i)=nmat_sol(i-1)+nsizea
      enddo
      do i = 1, mdof_sol
         nrhs_sol(i)=(i-1)*neq
      enddo

      return
      end

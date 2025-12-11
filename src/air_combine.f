       subroutine air_combine(iflg)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.      
!***********************************************************************
CD1
CD1
CD1  This subroutine changes variable in air-water problems
CD1  PC version uses do loop instead of F90 because os stack overflow 
CD1  large problems
CD1  gaz 112219 used only for Richard's Ed (jswitch ne 0)      
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 
CD2    Rev 1.0   02/18/08 10:24:20   GAZ
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD1  PURPOSE 
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  CD3  2.5.2 Solve nonlinear equation set at each time step
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

      use comai
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comhi
      use comii
      use davidi
      use comdti
      use commeth
      implicit none

      integer iflg, idofm, nmatd 
      integer i, icoupl_iter, id, idl, i1, i2, jj, j, kb
      integer icd,ii1,ii2
      integer neqp1, nrhs1, nrhs2, nrhs3, nsizea, nsizea1 
      real*8, allocatable :: dumz(:)
      real*8, allocatable :: dumn(:)
      real*8, allocatable :: sto5(:)   
      real*8  facr, fdum2, tollr, tolls, dtp, dshsw, bpsave
      real*8 cap_tol,bpsave2
      parameter (cap_tol= 1.e-2)

c     put memory management in

c     
c     degrees of freedom reduced from 2 to 1 for continuum
c     degrees of freedom reduced from 4 to 2 for dpdp
c     
      if(idpdp.ne.0) then
c     
c     dpdp code
c     

         if(iflg.eq.1) then
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     shift jacobean array for two phase region
c     s is variable
c     
c     first neq nodes    
            icd = 0
            do i = 1,neq 
               ii1=nelm(i-icd)+1
               ii2=nelm(i-icd+1)
               do j = ii1,ii2
                  kb = nelm(j) + icd
                  if(ieos(kb).ne.1) then
                     jj = j - neqp1
                     a(jj+nmat(1)) = a(jj+nmat(2))
                  endif
               enddo 
            enddo
c     
c     neq +1 to 2*neq nodes 
c     
            icd = neq
            do i = neq+1,2*neq 
               ii1=nelm(i-icd)+1
               ii2=nelm(i-icd+1)
               do j = ii1,ii2
                  kb = nelm(j) + icd
                  if(ieos(kb).ne.1) then
                     jj = j - neqp1
                     a(jj+nmat(11-joff)) = a(jj+nmat(12-joff))
                  endif
               enddo 
            enddo
c     
c     f-m interaction terms (1 to neq)
c     
            icd = 0
            do i = 1,neq 
               j = nelmdg(i-icd)
               if(ieos(i).ne.1) then
                  jj = j - neqp1
                  a(jj+nmat(9-joff)) = a(jj+nmat(10-joff))
               endif
            enddo
c     
c     f-m interaction terms (neq to 2*neq)
c     
            icd = neq
            do i = neq+1,2*neq 
               j = nelmdg(i-icd)
               if(ieos(i).ne.1) then
                  jj = j - neqp1
                  a(jj+nmat(3)) = a(jj+nmat(4))
               endif
            enddo

c     
c     condense arrays (first block is OK)
c     
c            a(nmat(2)+1:nmat(3)) = a(nmat(3)+1:nmat(4))
c            a(nmat(3)+1:nmat(4)) = a(nmat(9-joff)+1:nmat(10-joff))
c            a(nmat(4)+1:nmat(5)) = a(nmat(11-joff)+1:nmat(12-joff))
c gaz 040416 PC stack overflow correction
             do i = 1, nmat(2)
               a(nmat(2) + i) = a(nmat(3) + i)
             enddo
             do i = 1,nmat(2)
               a(nmat(3) + i) = a(nmat(9-joff) + i)
             enddo
             do i = 1,nmat(2)
               a(nmat(4) + i) = a(nmat(11-joff) + i)
             enddo
c     
c     arrange residual arrays
c     
            do i = 1,neq
               bp(i+nrhs(2)) = bp(i+nrhs(3))
            enddo

            idofm = 2
c     
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     
            do i = 1,neq 
               bpsave  = bp(i+nrhs(1))
               bpsave2 = bp(i+nrhs(2))
               if(ieos(i).ne.1) then
                  bp(i+nrhs(1)) = 0.0d0 
                  bp(i+nrhs(2)) = bpsave
c gaz 110819 pref, tref (global) read in scanin crl(4,1) repaced with pref                   
                  phi(i) = pref
               else 
c     bp(i+nrhs(1)) is correct       
                  bp(i+nrhs(2)) = 0.0d0	   	   
               endif	
               j = i +neq  
               if(ieos(j).ne.1) then
                  bp(i+nrhs(3)) = 0.0d0 
                  bp(i+nrhs(4)) = bpsave2
c gaz 110819 pref, tref (global) read in scanin crl(4,1) repaced with pref                   
                  phi(j) = pref 
               else 
                  bp(i+nrhs(3)) = bpsave2
                  bp(i+nrhs(4)) = 0.0d0	   	   
               endif	
            enddo

c     
            idofm = 4
c     
         endif                  

      else  
c     
c     continuum
c     

         if(iflg.eq.1) then
c     
            neqp1=neq+1
c     nmatd is size of one subarray
            nmatd=nelm(neqp1)-neqp1
c     
c     continuum code
c     
            do i = 1,n
               i1=nelm(i)+1
               i2=nelm(i+1)

               do j = i1,i2
                  kb = nelm(j)  
                  if(ieos(kb).ne.1) then 
                     jj = j - neqp1
                     a(jj+nmat(1)) = a(jj+nmat(2))
                  endif
               enddo 
            enddo
c     
            idofm = 1
c     
         else if(iflg.eq.2) then
c     
c     replace variables in correct arrays
c     
            do i = 1,n
               if(ieos(i).ne.1) then
                  bpsave = bp(i+nrhs(1))
                  bp(i+nrhs(1)) = 0.0d0 
                  bp(i+nrhs(2)) = bpsave 
c gaz 110819 pref, tref (global) read in scanin crl(4,1) repaced with pref 
                  phi(i) = pref
               else 
                  bp(i+nrhs(2)) = 0.0d0	   	   
               endif	  
            enddo
c     
            idofm = 2
c     
         endif                  
      endif
      return
      end

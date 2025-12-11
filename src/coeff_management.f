      subroutine coeff_management(iflg)   
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
CD1  To zero out(not eliminate) negative area connections)
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD4  gaz june 19 2012
CD4  Initial implementation for doubly defined nodes only.
CD4
C***********************************************************************

      use davidi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer iflg,i,j,ij,kb,neqp1,ii,jm, inoar,inoarp
      integer node,ipar,npar,isx, neqp1_old
      integer n_ncon,kbmin,kbmax, nr_old, jk
      integer num_md,k,jj,kk,i1,i2,icode,neq_old         
      real*8 sx_md, sx_nr, area_dx, area_tol, dis  
      real*8 area_lower_limit
      parameter (area_lower_limit = 1.d-5)        


c=======================================================================
c 
c iflg      - flag to regulate the activities of this subroutine  
c sx        - array of finite element flow coefficients
c=======================================================================
c 
c don't do for anisotropic,dpdp, or gdkm problems
c
      if(ianpe.ne.0.or.gdkm_flag.ne.0.or.idpdp.ne.0) return
c            
      if(iflg.eq.0) then
c     read instructions from input file
c    

      else if(iflg.eq.1) then
c zero out negative coefficients 
         neqp1 = neq+1
         jk = 0
         do i = 1,neq
          i1 = nelmdg(i) + 1
          i2 = nelm(i+1)
          do j = i1,i2
           iw = istrw(j-neqp1)
           area_tol = (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))
            if(area_tol.gt.area_lower_limit) then
             if(isoy.eq.1) then
              sx(iw,isox)= -area_lower_limit/3.
              sx(iw,isoy)= -area_lower_limit/3.
              sx(iw,isoz)= -area_lower_limit/3.
             else
              sx(iw,isox)= -area_lower_limit
              sx(iw,isoy)= 0.
              sx(iw,isoz)= 0.             
             endif
             jk = jk + 1
            endif
          enddo
         enddo
       if(iptty.ne.0.and.jk.gt.0)
     & write (iptty,*) '>>>> number of neg coefficients zeroed = ', jk
       if(iout.ne.0.and.jk.gt.0)
     & write (iout,*) '>>>> number of neg coefficients zeroed = ', jk
       numcoef_neg = jk
      else if(iflg.eq.2) then
c remove zero coefficients      
      endif
      return
      end  
      

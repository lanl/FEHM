      subroutine check_sx(neq,ianpe,ndim_sx,ncon,istrw,
     &     nelmdg,ischk,iptty,iout,ierr,min_nr,isox,isoy,isoz,intg,
     &     sx1,sx)
!***********************************************************************
!  Copyright, 1995, 2004,  The  Regents of the University of California.
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
CD1  To check volumes and finite element flow coefficients.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/check_sx.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Mon Mar 31 08:31:16 1997   gaz
CD2 changed in checking geometric coefs
CD2 
CD2    Rev 1.7   Wed Jan 17 11:14:36 1996   zvd
CD2 Passed ierr value to check_sx from datchk
CD2 
CD2    Rev 1.6   Wed Jan 17 10:10:34 1996   llt
CD2 added type for ierr
CD2 
CD2    Rev 1.5   Wed Jan 17 09:48:40 1996   zvd
CD2 Added write to ierr and added minor updates to prolog
CD2 
CD2    Rev 1.4   Wed Jan 10 10:56:26 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.3   Tue Jan 09 15:12:04 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   08/16/95 10:58:42   zvd
CD2 Corrected write to iatty when unassigned.
CD2 
CD2    Rev 1.1   04/25/95 09:14:20   llt
CD2 added log history information
CD2
CD2    Rev 1.0   01/28/95 13:54:02   llt
CD2 water balance equation was modified
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.2 Finite-Element Coefficient Generation
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

c-------------------------------------------------------------------
c this routine is called from subroutine datchk
c this routine does not call any other routines
c-------------------------------------------------------------------

      implicit none
      integer neq,ndim_sx,min_nr,ncon(*),istrw(*),nelmdg(*)
      integer i,j,i2,iout,iptty,num_neg,isx,ipiv,ischk,max_num
      integer ianpe,intg
      real*8 sx1(*)
      real*8 sx(ndim_sx,*),sumsx,zero_sx,vol_max,vol_min
      real*8 sum_area_coef
      character*5 pintg
      parameter (zero_sx=-1.d-10)
      integer ierr,isox,isoy,isoz

c===================================================================
c     neq        - number of nodes
c     ndim_sx    - first dimension of coefficient matrix
c     ncon       - connectivity array
c     istrw      - pointer array for coefficients
c     nelmdg     - position in connectivity array of diagonal term
c     i          - do loop index
c     j          - do loop index
c     kb         - neighbor node 
c     isx        - index for coefficient matrix
c     ipiv       - do loop index 
c     num_neg    - number of negative volumes or coefficients
c     iptty      - unit number for output file
c     ischk      - unit number for datacheck file
c     max_num    - maximun number of negative volume or coefficients
c                - before the program stops
c     min_nr     - minumum storage needed for fe coefficients
c     sx1        - volume array
c     sx         - flow coefficient array
c     sumsx      - intermediate sum of coefficients
c     zero_sx    - tolerance for flow coeficients
c     vol_max    - maximum volume or flow coefficient
c===================================================================
c     
      max_num=neq/20
c     
c     first check volumes
c     
      if(ianpe.ne.0) then
         write(ischk,*) 'not checking FE coeffs. (anisotropy)'
         return
      endif
      num_neg=0
      vol_max=-1.d20
      do i=1,neq
         vol_max=max(sx1(i),vol_max)
         if(abs(sx1(i)).lt.zero_sx) then
            sx1(i)=0.0
         endif
         if(sx1(i) .lt. 0.0) then
            num_neg=num_neg+1
            if (iptty .ne. 0) write(iptty,100) i,sx1(i)
            if (iout .ne. 0) write(iout,100) i,sx1(i)
            if (ischk .ne. 0) write(ischk,100) i,sx1(i)
 100        format(1x,'warning : negative volume at node ',i8,
     *           ' volume = ',g15.4)
            if(num_neg .ge. max_num) then         
               if (iptty .ne. 0) write(iptty, 110) 
               write(ierr, 110) 
               if (iout .ne. 0) write(iout, 110)
               if (ischk .ne. 0) write(ischk, 110) 
               stop
            endif
 110        format(1x,'too many negative volumes : stopping')
         endif
      enddo
      if(num_neg.ne.0) then
         if (iptty .ne. 0) write(iptty, 120) vol_max
         if (iout .ne. 0) write(iout, 120) vol_max
         if (ischk .ne. 0) write(ischk, 120) vol_max
      endif
 120  format(1x,'maximum volume = ', g15.4)
c     
c     next check coefficients(assume isotropic conditions)
c     
c     change bad coefficient maximum
      max_num=neq/20*27
c GAZ
      max_num=1000000000
      vol_max=-1.d20
      vol_min=+1.d20
      num_neg=0
      min_nr=0
      sum_area_coef = 0.0d00
      do i=1,neq
         ipiv=nelmdg(i)
         i2=ncon(i+1)
         do j=ipiv+1,i2
c GAZ 1-25-03 put abs( )
            isx=abs(istrw(j-(neq+1)))
            min_nr=max(isx,min_nr)
            sumsx=-(sx(isx,isox)+sx(isx,isoy)+sx(isx,isoz))
            sum_area_coef = sum_area_coef + sumsx
            vol_max=max(sumsx,vol_max)
            vol_min=min(sumsx,vol_min)
            if(sumsx .lt. zero_sx) then
               num_neg=num_neg+1
c              if (iptty .ne. 0) then
c                 write(iptty,200) i,ncon(j),sumsx
c              end if
c              write(iout,200) i,ncon(j),sumsx
               if (ischk .ne. 0) write(ischk,200) i, ncon(j), sumsx
 200           format(1x,'warning : negative fe coeffcient at node ',i7,
     &              ' neighbor node ',i7,' coeff sum = ',g15.4)
               if(num_neg.ge.max_num) then
                  if (iptty .ne. 0) write(iptty, 210)
                  write(ierr, 210) 
                  if (iout .ne. 0) write(iout, 210) 
                  if (ischk .ne. 0) write(ischk, 210)
                  stop
               endif
 210           format(1x, 'too many negative coefficients : stopping')
c            sx(isx,1)=0.0
c            sx(isx,2)=0.0
c            sx(isx,3)=0.0
            endif
         enddo
      enddo
      if(num_neg.ne.0) then
               
         if (iptty .ne. 0) write(iptty,300) num_neg
         if (iout .ne. 0) write(iout,300) num_neg
         if (ischk .ne. 0) write(ischk,300) num_neg
 300     format(1x, 'warning : ', i10,' neg area coefficients' 
     &        ,' warning')

         if (iptty .ne. 0) then
            write(iptty,*) 'maximum area divided by length = ',vol_max
            write(iptty,*) 'minimum area divided by length = ',vol_min
         end if
         if (iout .ne. 0) then
            write(iout,*) 'maximum area divided by length = ',vol_max
            write(iout,*) 'minimum area divided by length = ',vol_min
         end if
         if (ischk .ne. 0) then
            write(ischk,*) 'maximum area divided by length = ',vol_max
            write(ischk,*) 'minimum area divided by length = ',vol_min
         end if
      endif
c 
      if(intg.le.0) pintg = 'node '
      if(intg.gt.0) pintg = 'gauss'     
      if (iptty .ne. 0) then
         write(iptty,*) 
     &    'volumes and fe coefficients checked intg = ', pintg
         write(iptty,*) 'Sum of area/dis for control volumes = ',
     &        sum_area_coef
      end if
      if (ischk .ne. 0) then
         write(ischk,*)
     &     'volumes and fe coefficients checked intg = ', pintg
         write(ischk,*) 'Sum of area/dis for control volumes = ',
     &        sum_area_coef
      end if
c     
      return
      end                

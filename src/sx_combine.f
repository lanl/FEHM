      subroutine sx_combine(iflg)
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
CD1  This subroutine combines the dimensions of a variable, used in 
CD1  generating the finite element coefficients. It also breaks
CD1  connections when Boundary Conditions connect nodes.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/sx_combine.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:15:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:26   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Mon Mar 31 08:43:36 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.2   Fri Apr 26 16:33:44 1996   gaz
CD2 major changes including directional componets
CD2 
CD2    Rev 1.1   Thu Jan 11 11:53:20 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.0   11/16/95 09:37:50   llt
CD2 Initial revision.
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  3.1     Finite Element Coefficient Generation
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
C***********************************************************************

      use davidi
      use comii
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      real*8, allocatable :: sxtot(:)
      real*8, allocatable :: sum_itfc(:)
c
      real*8  xi,yi,zi,xkb,ykb,zkb,dis2
      real*8  area_dis, sx2c,sx2c_tol
      integer iflg,i,i1,i2
      integer jj,kb,neqp1,icon, iw_con
      integer icode
      parameter (sx2c_tol=1.d-10)
c
c  return if stress solution enabled
c
      if(istrs.ne.0) return
c      
      neqp1=neq+1
      if(iflg.eq.-3) then
c     isotropic minimun storage case
         if(icnl.eq.0) then
            do i=1,nr
               sx(i,1)=sx(i,1)/3.0d00
            enddo
         else
            do i=1,nr
               sx(i,1)=sx(i,1)/2.0d00
            enddo
         endif
      elseif(iflg.eq.-2) then
c     isotropic minimun storage case
c     first combine(from fehm,subroutine anonp)
         allocate(sxtot(nr))
         if(icnl.eq.0) then
            do i=1,nr
               sxtot(i)=(sx(i,1)+sx(i,2)+sx(i,3))/3.0
            enddo
         else
            do i=1,nr
               sxtot(i)=(sx(i,1)+sx(i,2)+sx(i,3))/2.0
            enddo
         endif
         deallocate(sx)
         allocate(sx(nr,1))
         do i=1,nr
            sx(i,1)=sxtot(i) 
         enddo
         deallocate(sxtot)
      elseif(iflg.eq.-1) then
c     first combine(from fehm,subroutine anonp)
         do i=1,nr
            sx(i,1)=sx(i,1)+sx(i,2)+sx(i,3)
            sx(i,2)=0.0
            sx(i,3)=0.0
         enddo
c     divide into directional componets
         do i=1,neq
            xi=cord(i,1)
            yi=cord(i,2)
            zi=cord(i,3)
            i1=nelmdg(i)+1
            i2=nelm(i+1)
            do jj=i1,i2
               kb=nelm(jj)
               iw=istrw(jj-neqp1)
               xkb=cord(kb,1)
               ykb=cord(kb,2)
               zkb=cord(kb,3)
               dis2=(xkb-xi)**2+(ykb-yi)**2+(zkb-zi)**2
               if(dis2.gt.0.0) then
                  area_dis=sx(iw,1)          
                  sx(iw,1)=area_dis*(xkb-xi)**2/dis2
                  sx(iw,2)=area_dis*(ykb-yi)**2/dis2
                  sx(iw,3)=area_dis*(zkb-zi)**2/dis2
               endif
            enddo      
         enddo      

      elseif(iflg.eq.0) then
c     
c     first combine(from storsx read file)
c     
         allocate(sxtot(nr))
c     call mmgetblk ("sxtot","sx_combine",ipsxtot,nr,
c     &        2, icode)
c     
         do i=1,nr
c            sxtot(i)=-sqrt(sx(i,1)**2+sx(i,2)**2+sx(i,3)**2)
            sxtot(i) = sx(i,1)+sx(i,2)+sx(i,3)
            sx(i,1)=0.0
            sx(i,2)=0.0
            sx(i,3)=0.0
         enddo
c     divide into directional components
         do i=1,neq
            xi=cord(i,1)
            yi=cord(i,2)
            zi=cord(i,3)
            i1=nelm(i)+1
            i2=nelm(i+1)
            do jj=i1,i2
               kb=nelm(jj)
               iw=istrw(jj-neqp1)
               xkb=cord(kb,1)
               ykb=cord(kb,2)
               if(icnl.eq.0) then
                zkb=cord(kb,3)
                dis2=(xkb-xi)**2+(ykb-yi)**2+(zkb-zi)**2
                if(dis2.gt.0.0) then
                  area_dis=sxtot(iw)          
                  sx(iw,1)=area_dis*(xkb-xi)**2/dis2
                  sx(iw,2)=area_dis*(ykb-yi)**2/dis2
                  sx(iw,3)=area_dis*(zkb-zi)**2/dis2
                endif                
               else
                dis2=(xkb-xi)**2+(ykb-yi)**2
                 if(dis2.gt.0.0) then
                  area_dis=sxtot(iw)          
                  sx(iw,1)=area_dis*(xkb-xi)**2/dis2
                  sx(iw,2)=area_dis*(ykb-yi)**2/dis2
                endif
               endif
            enddo      
         enddo      
c     
c     call mmrelblk ("sxtot","sx_combine",ipsxtot, icode)
         deallocate(sxtot)
c     
      else if(iflg.eq.1.and.ianpe.eq.0) then
c     
c     break connections for specified pressure nodes to other specified
c     pressure nodes
c     
c     check for previously defined mdnodes
c     don't do this for anisotropy            
c     
c     break connections if neccessary
c     ka(i) lt 0 is pressure BC,ieos(i) lt 0  is large rservoir
c     gaz 4-10-01 (only for ka(i)=-1)
c     
c     gaz 10-19-2001 can't all break connections with solid phase yet
C  zvd 09-09-2005 use istrw_itfc instead of istrw to store broken connections
         red_factor(nitfcpairs+1)=0.d0
         
         if(ice.eq.0) then
            icon=0
            do i=1,neq
               if(ka(i).lt.0.or.ieos(i).lt.0) then
                  i1=nelmdg(i)+1 
                  i2=nelm(i+1)
                  do jj=i1,i2
                     kb=nelm(jj)
                     if(ka(kb).lt.0) then
                        istrw_itfc(jj-neqp1)=nitfcpairs+1
                        icon=icon+1
                     endif
                  enddo
               endif
            enddo
            do i=1,neq
               if(ka(i).lt.0.or.ieos(i).lt.0) then
                  i1=nelmdg(i)+1 
                  i2=nelm(i+1)
                  do jj=i1,i2
                     kb=nelm(jj)
                     if(ka(kb).ne.0.and.ka(kb).ne.-3) then
                        istrw_itfc(jj-neqp1)=nitfcpairs+1
                        icon=icon+1
                     endif
                  enddo
               endif
            enddo
         else
            icon=0
            do i=1,neq
               if(ka(i).lt.0) then
                  i1=nelmdg(i)+1 
                  i2=nelm(i+1)
                  do jj=i1,i2
                     kb=nelm(jj)
                     if(ka(kb).eq.-1) then
                        istrw_itfc(jj-neqp1) = nitfcpairs+1
                        icon=icon+1
                     endif
                  enddo
               endif
            enddo
            do i=1,neq
               if(ka(i).eq.-1) then
                  i1=nelmdg(i)+1 
                  i2=nelm(i+1)
                  do jj=i1,i2
                     kb=nelm(jj)
                     if(ka(kb).ne.0.and.ka(kb).ne.-3) then
                        istrw_itfc(jj-neqp1) = nitfcpairs+1
                        icon=icon+1
                     endif
                  enddo
               endif
            enddo
         endif
c     
c     write out warning that connections have been broken
c     
         if(icon.gt.0) then
            if (iout .ne. 0) 
     &           write(iout,*) icon,' BC to BC connection(s) found ',
     &           '(now set=0.0)'
            if(iptty.ne.0)
     &           write(iptty,*) icon,' BC to BC connection(s) found ',
     &           '(now set=0.0)'
         endif
c         go to 678
            allocate(sum_itfc(neq))
            sum_itfc = 0.0
            do i=1,neq
                  i1=nelmdg(i)+1 
                  i2=nelm(i+1)
                  do jj=i1,i2                                          
                     iw = istrw(jj-neqp1)
                     sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
                     if(istrw_itfc(jj-neqp1).ne.nitfcpairs+1) then
                      kb = nelm(jj)
                      sum_itfc(i) = sum_itfc(i) + sx2c
                      sum_itfc(kb) = sum_itfc(kb) + sx2c                      
                     endif
                  enddo
            enddo 
            icon = 0
            do i = 1, neq
               if(abs(sum_itfc(i)).lt.sx2c_tol) then
                write(ierr,678) i,sum_itfc (i)   
                icon = icon + 1
               endif
            enddo     
 678   format('warning: all area connections to node ',i5,
     &    ' are near zero',', sum of area connections = ', 1pg13.4)
         if(icon.gt.0) then
            if (ierr .ne. 0) 
     &           write(ierr,*) icon,' isolated gridblocks(s) found ',
     &       '           connections restored for isolated gridblocks'      
            if (iout .ne. 0) 
     &           write(iout,*) icon,' isolated gridblocks(s) found ',
     &       '           connections restored for isolated gridblocks'             
            if(iptty.ne.0)
     &           write(iptty,*) icon,' isolated gridblocks(s) found ',
     &       '           connections restored for isolated gridblocks'             
         endif

            do i = 1, neq
                i1 = nelm(i)+1
                i2 = nelm(i+1)
                do jj = i1, i2
                 kb = nelm(jj)               
                 if(abs(sum_itfc(kb)).lt.sx2c_tol.or.
     &            abs(sum_itfc(i)).lt.sx2c_tol) then
                   istrw_itfc(jj-neqp1) = 0
                 endif
               enddo 
            enddo           
          
c     
      else if(iflg.eq.2) then
c     
c     set last connection = 0.0
c     
         sx(nr,1)=0.0
         sx(nr,2)=0.0
         sx(nr,3)=0.0
      else if(iflg.eq.3) then
c     
c     set last connection = 0.0
c     
         sx(nr,1)=0.0
      endif
      return
      end

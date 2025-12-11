      subroutine md_nodes(iflg,neighbors,ii)   
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
CD1  To manage multiply defined nodes.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2  $Log:   /pvcs.config/fehm90/src/mdnodes.f_a  $
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:22   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Mon Mar 31 08:39:44 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.10   Wed Jun 12 15:21:12 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.9   Mon Jun 03 11:18:14 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.8   Fri May 31 10:56:50 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.7   Fri Apr 26 15:49:12 1996   gaz
CD2 more mdnode changes
CD2 
CD2    Rev 1.6   Wed Feb 07 11:53:10 1996   gaz
CD2 corrections for carl's mdnode definition
CD2 only implimented for air-water isothermal
CD2 
CD2    Rev 1.5   Thu Jan 18 10:58:00 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.4   Fri Jan 12 17:52:00 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.3   Wed Jan 10 14:45:54 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   08/09/95 18:53:16   llt
CD2 moved implicit none up to top, for HP
CD2 added PVCS log history
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3  
CD3  2.2 Finite-Element Coefficient Generation
CD3  2.6 Provide Input/Output Data Files
CD3  3.0 INPUT AND OUTPUT REQUIREMENTS
CD3  
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
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
      integer neighbors(*),i_elim,iw_min, icoef         
      real*8 sx_md, sx_nr, area_dx, area_tol, dis
      parameter(area_tol = 1.d-30)
      parameter(inoarp=200)             

      integer, allocatable :: ncon_new(:)
      integer, allocatable :: istrw_new(:)
      integer, allocatable :: idum(:)
      integer, allocatable :: imd_con(:,:)
      real*8, allocatable :: area_con(:,:)
      real*8, allocatable :: sx_new(:,:)
      real*8, allocatable :: dum(:)

      save n_ncon,i_elim
      save imd_con
      save area_con, neq_old

c=======================================================================
c 
c iflg      - flag to regulate the activities of this subroutine
c n         - node investigated(iflg=1)
c           - total number of nodes(iflg=2)
c mdnodes   - list of connected nodes(iflg=2)
c neighbors - array of neighbors of node n(iflg=1)
c mdnodes   - number of multiply defined nodes associated with a given node   
c sx        - array of finite element flow coefficients
c=======================================================================
c       
      if(imdnode.eq.0) return
      if(iflg.eq.0) then
c     read as list
         read(inpt,*) num_md,mdmax,i_elim,sx_mult
c     
c     create storage for arrays
c     
         allocate(mdnode(neq_primary,mdmax))
         allocate(area_con(neq_primary,mdmax))
         allocate(mdnodes(neq_primary))
         area_con = 0.0d00
c     
c     node = node number
c     ipar=1 : parent node
c     ipar=0 : possible child
c     ipar=0 and npar ne node, node is child of parent npar
c     mdmax = max number of children
         do i=1,neq_primary
            mdnodes(i)=0
            do j=1,mdmax
               mdnode(i,j)=0
            enddo
         enddo
         do jj=1,num_md
            read(inpt,'(a80)') wdd1
            read(wdd1,*,end=20) node,ipar,npar,area_dx  
            go to 30
 20         continue
            area_dx = 0.0d00
 30         continue
            if(node.ne.npar) then
               mdnodes(npar)=mdnodes(npar)+1
               if(mdnodes(npar).gt.mdmax) then
                  write (ierr, 900) npar, mdnodes(npar), mdmax
                  if (iout .ne. 0) write(iout,900) npar, mdnodes(npar),
     &                 mdmax
                  if (iptty.ne.0) write(iptty,900) npar, mdnodes(npar),
     &                 mdmax
 900             format('>>> more new connections than mdmax: stopping',
     &                 /,'node ',i9,' new connections ',i6,' max ', i6)
                  stop
               endif
               mdnode(npar,mdnodes(npar))=node
c     use reciprocity to complete set      
               mdnodes(node)=mdnodes(node)+1
               if(mdnodes(node).gt.mdmax) then
                  write (ierr, 900) node, mdnodes(node), mdmax
                  if (iout .ne. 0) write(iout,900) node, mdnodes(node),
     &                 mdmax
                  if(iptty.ne.0) write(iptty,900) node, mdnodes(node),
     &                 mdmax
                  stop
               endif
c check for coincident nodes
               if(icnl.eq.0) then
                  dis=(cord(npar,1)-cord(node,1))**2
     &                 + (cord(npar,2)-cord(node,2))**2
     &                 + (cord(npar,3)-cord(node,3))**2
               else
                  dis=(cord(npar,1)-cord(node,1))**2
     &                 + (cord(npar,2)-cord(node,2))**2
               endif
c              if(dis.le.area_tol) then
               if(dis.le.-1.0    ) then
                  if (iout .ne. 0) write(iout,901) npar, node 
                  if(iptty.ne.0) write(iptty,901) npar, node 
 901              format('>>> mdnode msg: nodes ',2i9,' coincident',/,
     &                 'area set to zero, sx_mult takes precidence')
                  area_dx = 0.0    
               endif
               mdnode(node,mdnodes(node))=npar
               area_con(node,mdnodes(node)) = -abs(area_dx)
               area_con(npar,mdnodes(npar)) = -abs(area_dx)
            endif
         enddo
c
c eliminate parent nodes
c
         if(i_elim.gt.0) then
            do i=1,neq
               if(mdnodes(i).gt.1) then
c identify one child(kk)   
                  kk=mdnode(i,1)
c identify another child(ii)   
                  jj=mdnode(i,2)
c replace parent with child in child mdnode array
c assumes child mdnode array has only one parent
                  mdnode(kk,1)=jj
                  mdnode(jj,1)=kk
c eliminate parent node(below action is sufficient)
                  mdnodes(i)=0
               endif
            enddo
         endif               

      else if(iflg.eq.-1) then
c new ifblock to complement wellctr
c wellctr(0) must be called first in order to have
c correct memory allocated
c gaz 10-22-03
         read(inpt,*) num_md,mdmax,i_elim,sx_mult
c     
c     create storage for arrays
c     
         allocate(mdnode(neq_primary,mdmax))
         allocate(area_con(neq_primary,mdmax))
         allocate(mdnodes(neq_primary))
         area_con = 0.0d00
c     
c     node = node number
c     ipar=1 : parent node
c     ipar=0 : possible child
c     ipar=0 and npar ne node, node is child of parent npar
c     mdmax = max number of children
         do i=1,neq_primary
            mdnodes(i)=0
            do j=1,mdmax
               mdnode(i,j)=0
            enddo
         enddo
         do jj=1,num_md
            read(inpt,'(a80)') wdd1
            read(wdd1,*,end=21) node,ipar,npar,area_dx  
            go to 31
 21         continue
            area_dx = 0.0d00
 31         continue
            if(node.ne.npar) then
               mdnodes(npar)=mdnodes(npar)+1
               if(mdnodes(npar).gt.mdmax) then
                  write (ierr, 910) npar, mdnodes(npar), mdmax
                  if (iout .ne. 0) write(iout,910) npar, mdnodes(npar),
     &                 mdmax
                  if (iptty.ne.0) write(iptty,910) npar, mdnodes(npar),
     &                 mdmax
 910             format('>>> more new connections than mdmax: stopping',
     &                 /,'node ',i9,' new connections ',i6,' max ', i6)
                  stop
               endif
               mdnode(npar,mdnodes(npar))=node
c     use reciprocity to complete set      
               mdnodes(node)=mdnodes(node)+1
               if(mdnodes(node).gt.mdmax) then
                  write (ierr, 910) npar, mdnodes(node), mdmax
                  if (iout .ne. 0) write(iout,910) node, mdnodes(node), 
     &                 mdmax
                  if(iptty.ne.0)write(iptty,910) node, mdnodes(node), 
     &                 mdmax
                  stop
               endif
c check for coincident nodes
               if(icnl.eq.0) then
                  dis=(cord(npar,1)-cord(node,1))**2
     &                 + (cord(npar,2)-cord(node,2))**2
     &                 + (cord(npar,3)-cord(node,3))**2
               else
                  dis=(cord(npar,1)-cord(node,1))**2
     &                 + (cord(npar,2)-cord(node,2))**2
               endif
c              if(dis.le.area_tol) then
               if(dis.le.-1.0    ) then
                  if (iout .ne. 0) write(iout,911) npar, node 
                  if (iptty.ne.0) write(iptty,911) npar,node 
 911              format('>>> mdnode msg: nodes ',2i9,' coincident',/,
     &                 'area set to zero, sx_mult takes precidence')
                  area_dx = 0.0    
               endif
               mdnode(node,mdnodes(node))=npar
               area_con(node,mdnodes(node)) = -abs(area_dx)
               area_con(npar,mdnodes(npar)) = -abs(area_dx)
            endif
         enddo
c
c eliminate parent nodes
c
         if(i_elim.gt.0) then
            do i=1,neq
               if(mdnodes(i).gt.1) then
c identify one child(kk)   
                  kk=mdnode(i,1)
c identify another child(ii)   
                  jj=mdnode(i,2)
c replace parent with child in child mdnode array
c assumes child mdnode array has only one parent
                  mdnode(kk,1)=jj
                  mdnode(jj,1)=kk
c eliminate parent node(below action is sufficient)
                  mdnodes(i)=0
               endif
            enddo
         endif               

         
      else if(iflg.eq.1) then
c     modify connectivity array 
c     called from anonp when generating ncon     
         if(mdnodes(ii).ne.0) then
            neighbors(ii) = ii
            do i=1,mdnodes(ii)
               kb=mdnode(ii,i)
               neighbors(kb)=kb
            enddo   
         endif
      else if(iflg.eq.2) then
c     modify coefficient arrays
         neqp1=neq+1
         do i=1,neq
            if(mdnodes(i).ne.0) then
               do jj=1,mdnodes(i)
                  kk=mdnode(i,jj)
                  if(kk.gt.i) then  
                     i1=nelmdg(i)+1
                     i2=nelm(i+1)
                     do j=i1,i2
                        isx=istrw(j-neqp1)
                        if(kk.eq.nelm(j)) then
                           sx(isx,isox) =  0.0
                           sx(isx,isoy) =  0.0
                           sx(isx,isoz) =  0.0
                        endif
                     enddo        
                  endif             
               enddo     
            endif
         enddo
c
c eliminate parent nodes
c
         if(i_elim.gt.0) then
            do i=1,neq
               if(mdnodes(i).ge.2) then
c identify one child(kk)   
                  kk=mdnode(i,1)
c identify another child(ii)   
                  ii=mdnode(i,2)
c replace parent with child in child mdnode array
c assumes child mdnode array has only one parent
                  mdnode(kk,1)=ii
                  mdnode(ii,1)=kk
c relace parent with child in connectivity array
                  i1=nelm(kk)+1
                  i2=nelm(kk+1)
                  do j=i1,i2
                     if(i.eq.nelm(j)) then
                        nelm(j)=ii           
                        go to 110
                     endif
                  enddo        
 110              continue
                  i1=nelm(ii)+1
                  i2=nelm(ii+1)
                  do j=i1,i2
                     if(i.eq.nelm(j)) then
                        nelm(j)=kk           
                        go to 210
                     endif
                  enddo        
 210              continue
c eliminate parent node(below action is sufficient)
                  mdnodes(i)=0
               endif
            enddo
         endif               

      else if(iflg.eq.3.and.i_elim.eq.0) then
         call geneqmdnode
      else if(iflg.eq.4) then
c     
c     modify connectivity and coefficients when a stor file is used
c     think about adding nodes
c     just for read in stor files
c     
         n_ncon=ii            
         do i=1,neq_primary 
            if(mdnodes(i).ne.0) then
               n_ncon=n_ncon+mdnodes(i)
            endif
         enddo
c     
c     add space for extra nodes( old neq is abs(imdnode))
c     
         if(imdnode.lt.0) then
            neq_old = abs(imdnode)
c     at least add mdmax+2 spaces (node+connections+pointer)
            do i=neq_old+1,neq_primary
               n_ncon=n_ncon + mdmax+2
            enddo
         else
            neq_old=neq_primary
         endif
         ii=n_ncon
      else if(iflg.eq.5) then
c     
         neqp1_old=neq_old+1
         neqp1 = neq_primary +1

         allocate(istrw_new(n_ncon-neqp1))
         allocate(ncon_new(n_ncon))
         allocate(idum(2*neq_primary))
         allocate(dum(n_ncon-neqp1))
         j=neqp1
c        sx_nr = sx(nr,isox) + sx(nr,isoy) + sx(nr,isoz)
c        icoef = nr-1
c gaz 11-09-2001 need to check out implications
         icoef = nr
c        nr_old = nr      
c  gaz 1-23-03
         nr_old = nr-1      
         ncon_new(1)=neqp1
         do i=1,neq_primary
            idum(i)=0
            idum(neq_primary+i)=0
         enddo
         do i=1,neq_primary
            idum(i)=1
            kbmin=neq_primary
            kbmax=0
            kbmin=min(kbmin,i)
            kbmax=max(kbmax,i)
            if(i.le.neq_old) then
               i1=nelm(i)+1     
               i2=nelm(i+1) 
               do jj=i1,i2
                  kb=nelm(jj)
                  kbmin=min(kbmin,kb)
                  kbmax=max(kbmax,kb)
                  idum(kb)=1
                  idum(neq_primary+kb)=istrw(jj-neqp1_old)
                  if(idum(neq_primary+kb).eq.nr) 
     &                    idum(neq_primary+kb) = -99999999
               enddo
            endif 
            do k=1,mdnodes(i)
               kb=mdnode(i,k)
               kbmin=min(kbmin,kb)
               kbmax=max(kbmax,kb)
               idum(kb)=1
               if(kb.gt.i) then
                  icoef = icoef +1
                  if(abs(area_con(i,k)).gt.area_tol) then
                     idum(neq_primary+kb) = icoef
                  else
                     idum(neq_primary+kb) = -icoef
                  endif
                  dum(abs(icoef)) = area_con(i,k)
               else if(kb.lt.i) then
                  do jk = nelmdg(kb)+1, ncon_new(kb+1)
                     if(ncon_new(jk).eq.i) then
                        idum(neq_primary+kb) = istrw_new(jk-neqp1)
                     endif
                  enddo
               endif
            enddo
            do k=kbmin,kbmax
               if(idum(k).ne.0) then
                  j=j+1
                  ncon_new(j)=k
                  istrw_new(j-neqp1)=idum(k+neq_primary)
                  if(k.eq.i) nelmdg(i)=j
                  idum(k)=0
                  idum(k+neq_primary)=0
               endif
            enddo
            ncon_new(i+1)=j
         enddo
c
c deallocate and allocate some arrays
c
         deallocate(nelm,istrw)
         allocate (nelm(j))
         allocate (istrw(j-neqp1))
         do i=1,j
            nelm(i)=ncon_new(i)
         enddo
         do i=1,j-neqp1 
            if(istrw_new(i).eq.-99999999) istrw_new(i) = icoef+1
            istrw(i)=istrw_new(i)
         enddo
         nr = icoef+1

         if(isoy.ne.1) then
            allocate (sx_new(icoef+1,3))
         else
            allocate (sx_new(icoef+1,1))
         endif

         do i=1,nr_old
            sx_new(i,isox) = sx(i,isox)
            sx_new(i,isoy) = sx(i,isoy)
            sx_new(i,isoz) = sx(i,isoz)
         enddo
         
         deallocate(sx)
         if(isoy.ne.1) then
            allocate (sx(icoef+1,3))
         else
            allocate (sx(icoef+1,1))
         endif
         
         do i=1,nr_old
            sx(i,isox) = sx_new(i,isox)
            sx(i,isoy) = sx_new(i,isoy)
            sx(i,isoz) = sx_new(i,isoz)
         enddo

         do i = nr_old+1,icoef
            if(isoy.ne.1) then
               sx(i,1) = dum(i)             
               sx(i,2) = 0.0d00          
               sx(i,3) = 0.0d00          
            else
               sx(i,1) = dum(i)             
            endif
         enddo


         deallocate(idum,dum,istrw_new,ncon_new,sx_new)
      elseif(iflg.eq.6) then
c
c calculate flow coefficient for MD nodes
c
c gaz 6-1-2001 coding to let mdnode point to largest area coefficient
c initially set iw to minus so can tell mdnodes connections
         sx(nr,isox)=0.0  
         sx(nr,isoy)=0.0  
         sx(nr,isoz)=0.0  
c important neqp1 used only in iflg=6 for istrw needs
c to be based on neq,not neq_primary
         neqp1=neq+1
         do i=1,neq_primary
            if(mdnodes(i).eq.1) then
               if(abs(area_con(i,1)).lt.area_tol) then
                  sx_md=+1.e20
                  iw_min = 0 
                  i1=nelm(i)+1
                  i2=nelm(i+1)
                  do jm=i1,i2
                   kb=nelm(jm)
                   if(kb.ne.i) then
                     iw=istrw(jm-neqp1)
                     if(iw.gt.0) then
                       if(sx_md.gt.sx(iw,isox)+
     &                    sx(iw,isoy)+sx(iw,isoz)) then
                        sx_md = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
                        iw_min = iw
                       endif
                     endif
                   endif
                  enddo
                  do jm=i1,i2
                     kb=nelm(jm)
                     if(kb.eq.mdnode(i,1)) then
                        istrw(jm-neqp1) = -iw_min
                     endif
                  enddo
               endif
            else if(mdnodes(i).gt.1) then
               do jj=1,mdnodes(i)
                  if(abs(area_con(i,jj)).lt.area_tol) then
                     sx_md=+1.e20
                     iw_min = 0 
                     j=mdnode(i,jj)
                     ij=min(i,j)
                     i1=nelm(ij)+1
                     i2=nelm(ij+1)
                     do jm=i1,i2
                        iw=istrw(jm-neqp1)
                        kb=nelm(jm)
                       if(kb.ne.ij) then
                        if(iw.gt.0) then
                         if(sx_md.gt.sx(iw,isox)+
     &                       sx(iw,isoy)+sx(iw,isoz)) then
                           sx_md = sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
                           iw_min = iw
                         endif
                        endif
                       endif
                     enddo
                     i1=nelm(i)+1
                     i2=nelm(i+1)
                     do jm=i1,i2  
                        kb=nelm(jm)
                        if(kb.eq.j) then
                           istrw(jm-neqp1) = -iw_min
                        endif
                     enddo
                  endif                       
               enddo
            endif                       
         enddo

         deallocate(area_con)
         inoar = 0
         if (ischk .ne. 0) then
            write(ischk,*) 
            write(ischk,*) 'After connecting multiply defined nodes..' 
         end if
         do i= 1,neq_primary
            i2=nelm(i+1)   
            do jm=nelmdg(i)+1,i2
c        istrw(jm-neqp1) = abs(istrw(jm-neqp1))
c        iw = istrw(jm-neqp1)
               iw = abs(istrw(jm-neqp1))
               sx_md = sx(iw,isox) + sx(iw,isoy) + sx(iw,isoz)
               if(abs(sx_md).lt.area_tol) then 
                  inoar = inoar +1
                  if(inoar.le.inoarp) then
                     kb = nelm(jm)
                     if (ischk .ne. 0) write(ischk,*) 
     &                    'connection with no area,i = ',i,' j = ',kb
                  endif
               endif
            enddo
         enddo
         if (ischk .ne. 0) then
            write(ischk,*) 'total connections with no area = ',inoar 
            write(ischk,*)
         end if
      endif
      return
      end  
      

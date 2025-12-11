      subroutine gencof(nodeu,neluc,nsf,ineluf,inoduf,neluf,ncon,
     *     nodeuf,iplace,noodum,nopdum,ndimnoo,maxisx)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
CD1  To generate finite element coefficients (mixed elements) and 
CD1  perform numerical integration of the elements.
CD1  stress coefficients are assumed non symmetric
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/gencof.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:04   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:10   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Fri Apr 26 15:13:12 1996   gaz
CD2 corrections for 2-d elements in 3-d space
CD2 
CD2    Rev 1.10   Thu Jan 18 09:22:46 1996   zvd
CD2 Modifications to prolog
CD2 
CD2    Rev 1.9   Fri Jan 12 17:47:54 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.8   Wed Jan 10 11:40:06 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.7   Wed Jan 10 09:47:18 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.6   Tue Jan 09 14:02:54 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.5   11/15/95 16:15:50   gaz
CD2 changes for sx(iw,3) instead of sx(iw,6)
CD2 
CD2    Rev 1.4   01/28/95 13:54:46   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.3   06/03/94 15:39:04   gaz
CD2 Made icsh parameter change for prisms and quads. This fixed a bug in the inetra
CD2 gration of these elements.
CD2 
CD2    Rev 1.2   03/18/94 15:57:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   03/08/94 15:45:56   llt
CD2 moved the allocation and deallocation of array aj to the top and bottom of
CD2 gencof (instead of inside gncf2 and gncf3).
CD2 
CD2    Rev 1.0   01/20/94 10:24:02   pvcs
CD2 original version in process of being certified
CD2 
c  11/14/94 gaz added call to thickness
c  12/15/94 gaz changed to noodum,nopdum
c  12/22/94 gaz passed isx to zero our sx
c 1/3/95 gaz removed call to thickness
c 1/4/95 gaz removed line to zero out sumsx
c 1/4/95 gaz skipped calc for nrq=8,9,10
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
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5
CD5   ineluf          INT      
CD5   inoduf          INT      
CD5   iplace          INT
CD5   maxisx          INT
CD5   ncon            INT      
CD5   ndimnoo         INT
CD5   neluc           INT      
CD5   neluf           INT      
CD5   nodeu           INT      
CD5   nodeuf          INT       
CD5   noodum          INT      
CD5   nopdum          INT      
CD5   nsf             INT       
CD5
CD5 Interface Tables
CD5
CD5   None
CD5
CD5 Files
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6
CD6                            COMMON
CD6   Identifier      Type     Block  Description
CD6
CD6   bcoef           REAL*8   fff    Scratch storage array for coefficient 
CD6                                     calculations
CD6   icnl            INT      faai   Problem dimension
CD6   intg            INT      faai   Indicates integration type used
CD6   ipbcoef         POINTER  fff    Pointer to variable array bcoef
CD6
CD6   istrs           INT      faai   Parameter indicating if the stress
CD6                                     solution is enabled
CD6   istrw           INT      fbb    Starting positions in sx(nr,9) array of
CD6                                     finite element coefficients for each
CD6                                     node
CD6   lenreal         INT      param  Converts bits to words for allocating
CD6                                     memory
CD6   n0              INT      param  Maximum number of nodes allowed
CD6   nei             INT      faai   Total number of elements in the problem
CD6   nelm            INT      fbb    Initially information about nodes in each
CD6                                     element, later nodal connectivity
CD6                                     information
CD6   nelmdg          INT      fbb    Contains position of (i,i) element in
CD6                                     connectivity array
CD6   nelucm          INT      param  nbd / 64
CD6   neq             INT      faai   Number of nodes, not including dual
CD6                                     porosity nodes
CD6   ni              INT      faai   Number of integration points per element
CD6   nop             INT      fbb    Matrix sparsity structure for lu
CD6                                     decomposition
CD6   nr              INT      param  Maximum space allowed for each finite
CD6                                     element coefficient array
CD6   ns              INT      faai   Number of nodes per element
CD6   sx              REAL*8   fbc    Contains finite element geometric
CD6                                     coefficients necessary for heat and mass
CD6                                     transfer simulation
CD6   sx1             REAL*8   fbc    Contains volume associated with each node
CD6   sxs             REAL*8   fbc    Contains more finite element geometric
CD6                                     coefficients (ie., those necessary for
CD6                                     the stress module)
CD6
CD6 Global Subprograms
CD6
CD6   Identifier      Type     Description
CD6
CD6   gncf2           N/A      Generate 2-D finite element coefficients
CD6   gncf3           N/A      Generate 3-D finite element coefficients
CD6   mmgetblk        N/A      Allocates space for an array
CD6   mmrelblk        N/A      Deallocates space for an array
CD6   zeror_out       N/A      Assign value of zero to all elements of a
CD6                              real array
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7   None
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   aj              REAL*8
CD7   dumm            REAL*8
CD7   i               INT      Loop index
CD7   i0              INT      Loop index
CD7   i1              INT      
CD7   i2              INT      
CD7   i7              INT      
CD7   i8              INT      
CD7   icsh            INT
CD7   icode           INT      Error return code
CD7   ij              INT
CD7   ipaj            POINTER  Pointer to variable array aj
CD7   ipdumm          POINTER  Pointer to variable array dumm
CD7   iq              INT      Loop index      
CD7   intgo           INT
CD7   ipiv            INT
CD7   iw0             INT
CD7   isx             INT
CD7   j               INT      Loop index
CD7   j0              INT      Loop index
CD7   jset            INT      Loop index
CD7   k               INT
CD7   kb              INT
CD7   knum            INT      
CD7   kset            INT      Loop index
CD7   lz              INT
CD7   nele            INT
CD7   nelu            INT
CD7   neu             INT
CD7   nga             INT
CD7   noder           INT
CD7   nrq             INT      
CD7   nrqd            INT      Loop index
CD7   nrs             INT
CD7   nsl             INT
CD7   nterm           INT
CD7   numkset         INT
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************

      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer nodeu, neluc, ndimnoo, maxisx
      integer noodum(ndimnoo,*),nopdum(ndimnoo,*),iplace(*),ns2dum
      integer nsf(*),ncon(*),inoduf(*),ineluf(*),neluf(*),nodeuf(*)
      integer i, i1, i2, i7, i8, icsh, ij, intgo, ipiv, iq, isx, iw0
      integer j, jset, k, kb, knum, kset, lz, nele, nelu, neu, ns1, ns2
      integer nga, noder, nsl, nrq, nrqd, nrs, nterm, numkset
      integer icrem
      integer neqp1,jj
      real*8 tenth, arem
      real*8, allocatable ::  aj(:,:)
      real*8, allocatable :: dumm(:)
      real*8 sxx, syy, szz
      logical printc
      logical proc

      if(icnl.eq.0.and.ns.eq.4) then
         ns2dum=32
      else if(icnl.eq.0.and.ns.eq.3) then
         ns2dum=18
      else
         ns2dum=ns*ns
      endif 
      allocate(bcoef(neluc,ns2dum),dumm(n0),aj(neluc,10))
      do i=1,neluc
         do j=1,ns2dum
            bcoef(i,j)=0.0
         enddo
      enddo

c     note repeated element types
c     note repeated nodal geometries
c     
      if (istrs.ne.0) then
         if(icnl.eq.0) then
            nterm = 25         
            ns1 = 11
            ns2 = 25
         else
            nterm = 18
            ns1 = 11
            ns2 = 18
         endif
      else
         nterm=04
         ns1 = 0
         ns2 = 0
      endif
c     
c     determine maximum number of integration points
c     
      if (icnl.eq.0) then
         if(nsf(1).eq.4) then
            ni=1
         else
            ni=8
         endif
      else
         ni=4
      endif
c     save integration type
      intgo=intg
      tenth = max(neluc/10.d0,1.d0)
      icrem = 0
c     
c     start out with intg = -1  (close to finv)  
c     changed later as necessary (to intgo)
c     
      intg = -1
c     
c     loop on equation coefficient types
c     
      do nrqd=1,nterm
         nrq=nrqd
c     
c     don't loop for nrq = 8,9, or 10
c     
         if(nrq.lt.8.or.nrq.gt.10) then
            if ( nrq.eq.1) icsh=1
            if ( nrq.eq.5 )  then
               if ( intg.gt.0 )  then
                  intg=-1
                  icsh=1
               else
                  icsh=1
               endif
            endif
            if ( nrq.eq.7 )  then
               if (intgo.gt.0)  then
                  intg=intgo
                  icsh=1
               else
                  icsh=1
               endif
            endif
            if(nrq.eq.20) then
               icsh = 1
               intg = intgo
            endif
c     
c     zero out bcoef for next geometric type
c     
            do i=1,neluc
               nele=ineluf(i)
               nsl=nsf(nele)
               do j=1,nsl*nsl
                  bcoef(i,j)=0.0
               enddo
            enddo
c     
c     loop on integration points
c     
            printc = .false.
            if(nrq.eq.1)printc = .true. 
            if(nrq.eq.3.and.icnl.ne.0)printc = .true. 
            if(nrq.eq.4.and.icnl.eq.0)printc = .true. 
            if(nrq.eq.ns1.and.istrs.ne.0)printc = .true. 
            if(nrq.eq.ns2.and.istrs.ne.0)printc = .true. 
            do nga=1,ni
c     
c     loop on unique elements
c     
               icrem = 0
               do neu=1,neluc
                  arem = mod(dble(neu),tenth) 
                  if(arem.eq.0.and.nga.eq.ni.and.printc) then
                     icrem = icrem + 10
                     if (iout.ne.0) write(iout,1600) nrq, icrem
                     if (iptty.ne.0) write(iptty,1600) nrq, icrem
                  endif
1600   format(1x,'calcs for coef number = ',i5,1x,
     &   'percent complete = ',i5)               
                  nele=ineluf(neu)
c     calculate number of nodes
                  nsl=nsf(nele)
c     
c     shape function information calculated on nrq=1
c     jacobian information is calculated and stored on nrq=1
c     
c     three dimensional elements
                  if ( icnl.eq.0 )  then
                     if(nrq.eq.1) intg = 1
                     if(nrq.eq.23) intg = 1
                     if(nrq.eq.24) intg = 1
                     if(nrq.eq.25) intg = 1
                     if ( nga.le.nsl )  then
                        if(nsl.ge.6) icsh=1
                        call gncf3(nrq,nele,nga,neu,nsl,icsh,neluc,aj)
                     endif
                  endif
c reset integration type                  
                  intg = intgo
c     two dimensional elements
                  if ( icnl.ne.0 )  then
                     if(nrq.eq.1) intg = -1
                     if ( nga.le.nsl )  then
                        if(nsl.eq.4) icsh=1
                        call gncf2(nrq,nele,nga,neu,nsl,icsh,neluc,aj)
                     endif
                  endif
               enddo
            enddo
            
c     set integration type back to old one
c     should be taken care of above
c     intg=intgo

c     assemble sx(lz,nrq) from bcoef(nele,ij)
c     loop on unique node groups
c     
            isx=0
            dumm=0
            do noder=1,nodeu
c     find representative node of that group
               i=inoduf(noder)
c     perform calculations on node i
               i7=1
               i8=iplace(i)
               do j=i7,i8
                  nele=nopdum(i,j)
                  proc=.true.
                  if ( nele.eq.0) proc=.false.
                  if ( proc )  then
c     find set that nele belongs to
                     nelu=neluf(nele)
c     identify position of nodei in element nele
                     k=noodum(i,j)
                     nsl=nsf(nelu)
                     knum=(k-1)*nsl
c     add contribution from element nele
                     do lz=1,nsl
                        kb=nelm((nele-1)*ns+lz)
                        if ( kb.ne.0 )  then
                           ij=knum+lz
                           dumm(kb)=dumm(kb)+bcoef(nelu,ij)
                        endif
                     enddo
                  endif
               enddo
c     
c     load sx
c     
               i1=ncon(i)+1
               ipiv=nelmdg(i)
               i2=ncon(i+1)
               iw0=istrw(i)

               if ( nrq.eq.1 )  then
c     use mass lumping for capacitance term
                  do iq=i1,i2
                     kb=ncon(iq)
                     if(kb.eq.i) then
                      sx1(i)=sx1(i)+dumm(kb)
                      dumm(kb) = 0.0
                     endif
                  enddo
               endif
               if ( nrq.ne.1.and.nrq.le.4 )  then
c     rest of terms
                  do iq=ipiv+1,i2
                     kb=ncon(iq)
                     isx=isx+1
                     sx(isx,nrq-1)=sx(isx,nrq-1)+dumm(kb)
                     dumm(kb) = 0.0
                  enddo
                  
                  go to  777
                  if(nga.ge.ni) then
                     write(iout,*) 'nga ', nga, 'nrq ', nrq
                     write(iout,*) 'i,kb, iq, sz (kb) '
                     
                     i1 = ncon(i)+1
                     i2 = ncon(i+1)
                     do jj = i1,i2
                        kb = ncon(jj)
                        write(iout,*) i, kb,  nrq, dumm(kb) 
                     enddo

                  endif
 777              continue
               endif
               if(istrs.ne.0) then
               if ( nrq.ge.ns1.and.nrq.le.ns2-3 )  then
c     geometric coefficients (xy,xz,yz) for stress equations
                  nrs=nrq-10
                  do iq=i1,i2
                     kb=ncon(iq)
                     isx=isx+1
                     sxs(isx,nrs)=sxs(isx,nrs)+dumm(kb)
                     dumm(kb) = 0.0
                  enddo
               endif
               if ( nrq.ge.ns2-2.and.nrq.le.ns2 )  then
c     geometric coefficients (xy,xz,yz) for stress equations
                  nrs=nrq-10
                  do iq=i1,i2
                     kb=ncon(iq)
                     isx=isx+1
                     sxs(isx,nrs)=sxs(isx,nrs)+dumm(kb)
                     dumm(kb) = 0.0
                  enddo
               endif
               endif
            enddo
         endif
      enddo
c     
c     fill in other volumes from element groupings
c     
      do i=1,neq
         sx1(i)=sx1(inoduf(nodeuf(i)))
      enddo



      deallocate(bcoef,dumm,aj)
      
      return
      end

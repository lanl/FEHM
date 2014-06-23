      subroutine gdpm_corr (iflg)
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
CD1 PURPOSE
CD1
CD1 Manage gdpm connections when simulating rate-limited conditions.
CD1 IE make the connection value large for those nodes that serve only
CD1 to link gdpm nodes (diffusion only) to flow field
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 01-JAN-2007    G. Zyvoloski    Initial implementation.
CD2                        
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   iz              INT      I    Flag to indicate type of boundary
CD3                                   condition modification to be performed
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4   
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4
CD4 Global Subprograms
CD4
CD4   None
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   Identifier      Type     Description
CD5
CD5   sx1bc           REAL*8   /1.e12/
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   cqtout          REAL*8   Total produced tracer mass for each species at
CD5                              the boundary
CD5   cqtin           REAL*8   Total lost tracer mass for each species at
CD5                              the boundary
CD5   cqtrxn          REAL*8   Total produced tracer mass for each species at
CD5                              the boundary (reaction)
CD5   i               INT      Loop index
CD5   j               INT      Loop index
CD5   qtb             REAL*8   Total outflow for time step at the boundary
CD5   qtcb            REAL*8   Total mass injected at the boundary
CD5   qteb            REAL*8   Total energy outflow for time step at the
CD5                              boundary
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 ????? matrix pre-processing
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN bcon
CPS 
CPS 
CPS END bcon
CPS
C***********************************************************************
c
c adjust connections between some nodes and its neighbors 
c intended to be used with gdpm macro but might have general utility
c
      use combi
      use comdti      
      use comai 
      use comei
      use davidi
      use comki        
      implicit none

      integer ii, i, iflg, jj, i1, i2, i3, i4, kb, nmat_id, neqp1
      integer kk, ipiv, ipiv1, kc
      integer, allocatable ::  igdpm_rate_nodes_dum(:)      
      character*20 dummy
      logical null1, sheat_flag, stran_flag
      real*8  aijsave, aiidiff, corr_tol
      parameter (corr_tol = 1.e-18)
      
      if(igdpm_rate.le.0) return 
      
      if (iflg.eq.0) then
c     
c     read nodes and processes to be modified
c     
         allocate(igdpm_rate_nodes(n0))

         igdpm_rate = 1
         macro = 'cgdp'
	 
         
         read (inpt, '(a80)') wdd1
         select case (wdd1(1:4))
         case ('heat')
            read (wdd1, *) dummy,val_conh
            val_conh = max(val_conh*1.d-6,corr_tol)
            sheat_flag = .true.    
            igdpm_rate = 1
         case('tran')
            read (wdd1, *) dummy, val_conh
            val_conh = max(val_conh*1.d-6,corr_tol)
            stran_flag = .true.
            igdpm_rate = 2               
         end select
      
      
         igroup = 1
         narrays = 1
         itype(1) = 4
         default(1) = 0
      
c     
c     read in initial node  identifier here
c     
         call initdata2( inpt, ischk, n0, narrays,
     &        itype, default, macroread(8), macro, igroup, ireturn,
     &        i4_1 = igdpm_rate_nodes(1:n0))

      else if (iflg .eq. -1) then     
c     
c     expand the array size needed to store next nearest nodes  
c     just counting nodes now
c     
         do i= 1,n0
            if(igdpm_rate_nodes(i).gt.0) then
               i1 = nelm(i)+1
               i2 = nelm(i+1)
               do jj = i1,i2
                  kb = nelm(jj)
                  if(igdpm_rate_nodes(kb).eq.0) then
c     igdpm_rate_nodes(kb) = -1
                  endif
               enddo 
            endif
         enddo
         ngdpm_rate = 0
         kc = 0 
         do i= 1,n0
            if(igdpm_rate_nodes(i).gt.0) then
               ngdpm_rate = ngdpm_rate + 1
            else if(igdpm_rate_nodes(i).lt.0) then
               kc = kc+1
            endif
         enddo
         if (iout .ne. 0) write (iout,*)'gdpm primary nodes ',
     &        ngdpm_rate, 'next neighbor nodes ', kc 	 
      else if (iflg .eq. 1) then
c     
c     remove resistance from identified nodes 
c     
c     should be called just before array normalization and order switching
c     
         neqp1 = neq+1
         if(igdpm_rate.eq.1) then
c     heat conduction terms       
            
            do i = 1,n0
               if(igdpm_rate_nodes(i).gt.0) then        
                  call gdpm_geneqh(i)
               endif
            enddo 
         else if(igdpm_rate.eq.2) then
         endif
      end if

      end

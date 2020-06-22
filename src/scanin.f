      subroutine scanin
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
CD1 Scan input file for parameters needed prior to data input.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 13-JAN-94    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/scanin.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.12   Wed Jun 12 15:21:20 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.11   Mon Jun 03 11:18:34 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.10   Fri May 31 10:59:28 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.9   Wed Feb 07 12:13:28 1996   gaz
CD2 checked for bous macro
CD2 
CD2    Rev 1.8   Thu Feb 01 16:05:08 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.7   Tue Jan 09 14:08:20 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.6   12/13/95 10:30:24   robinson
CD2 Incorporated new zeolite hydration module
CD2 
CD2    Rev 1.5   12/13/95 08:47:00   gaz
CD2 changes to accodate variable rlp number of data
CD2 
CD2    Rev 1.4   04/25/95 10:11:20   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.3   01/31/95 16:24:42   llt
CD2 addition to the modification for the new particle tracking module
CD2 
CD2    Rev 1.2   01/28/95 14:03:58   llt
CD2 modified for new particle tracking module
CD2
CD2    Rev 1.1   03/18/94 16:03:54   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:26   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   inpt                     I    Main input data file.
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
c add grouping stuff 8/1/94
CD4 max_groups, ngroups, group, pos, n_couple_species, mdof_sol
CD4 num_particles
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   icnl            INT      faai   Problem dimension
CD4   idpdp           INT      faai   Parameter which indicates if the double
CD4                                     porosity / double permeability
CD4                                     solution is enabled
CD4   idualp          INT      faai   Parameter which indicates if the dual
CD4                                     porosity solution is enabled
CD4   inpt            INT      faai   Unit number for input file
CD4   neq             INT      faai   Number of nodes, not including dual
CD4   nspeci          INT      fdd1i  Number of species for tracer solution
CD4   wdd1            CHAR     faac   Alternate character input string
c add grouping stuff 8/1/94
CD4 ngroups         int     comrxni.h  The number of groups of 
CD4                                    coupled species
CD4 group           int     comrxni.h  Array indicating the coupling
CD4                                    between species
CD4 pos             int     comrxni.h  Array containing the species
CD4                                    numbers of the coupled species
CD4                                    in each group
CD4 n_couple_species int    comrxni.h  Array containg the number of
CD4                                    coupled species contained in
CD4                                    each grouping
CD4 mdof_sol     int     comrxni.h  The maximum number of degrees
CD4                                    of freedom necessary to solve
CD4                                    the tracer solution
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
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   adumm           REAL*8   Dummy real variable for input
CD5   cdumm           CHAR     Dummy character variable for input
CD5   idumm           INT      Dummy integer variable for input
CD5   ja              INT      Loop index input
CD5   jb              INT      Loop index input
CD5   jc              INT      Loop index input
CD5   macro           CHAR     Current macro being read
c changed to incorporate coupling 8/4/94
CD5   igrp          INT      Do loop index indicating the current
CD5                            group number
CD5   pos_index       INT      Index denoting the position of a 
CD5                            coupled species in a grouping
CD5   trac_flag       INT      Flag to search for group macro
CD5                            only if the rxn macro has been read
CD5   ispecies        INT      Species number for variable we are
CD5                            computing rate contribution for
CD5   
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
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9 2.8 Provide Multiple Realization Option
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
CPS BEGIN scanin
CPS 
CPS   rewind input tape
CPS      
CPS   REPEAT
CPS     read input line
CPS     EXIT IF macro read is stop
CPS     ELSE IF macro read is alti
CPS        read number of nodes
CPS     ELSE IF macro read is coor
CPS        read number of nodes
CPS     ELSE IF macro read is ctrl
CPS        read a line
CPS        REPEAT
CPS          read a line
CPS        UNTIL null line is found
CPS          read two lines
CPS          read geometry parameter
CPS     ELSE IF macro read is dual
CPS        set idualp to 1
CPS     ELSE IF macro read is dpdp
CPS        set idpdp to 1
CPS     ELSE IF macro read is trac
CPS        read 3 lines
CPS        read number of tracer species
CPS     ELSE IF macro read is rxn
CPS        read keyword
CPS        IF keyword is gr
CPS           read in the number of groupings
CPS           read in grouping indicator which indicates the coupling 
CPS           ... between species
CPS        FOR each group
CPS          FOR each species
CPS            IF the ith species is to be included in the ...
CPS            ... current group THEN
CPS              set up position array to indicate the species ...
CPS              ... number to be coupled
CPS              deterimine the number of coupled species in a group
CPS              determine the maximum # of degrees of freedom ...
CPS              ... needed to solve the tracer equations
CPS            ENDIF
CPS          ENDFOR each species
CPS        ENDFOR each group
c end change
CPS    
CPS     ELSE IF macro read is ptrk
CPS        set nspeci to 1
CPS        Read in the number of particles
CPS     END IF
CPS   UNTIL end of input file is found
CPS      
CPS   rewind input file
CPS   IF grup was not read then set variables to solve for each
CPS      species separately
CPS      set the total number of coupled groups to the total
CPS      number of species in the system
CPS      FOR each group
CPS        construct default grouping information
CPS      ENDFOR
CPS   ENDIF
CPS   
CPS END scanin
CPS
C***********************************************************************

      use comai
      use combi
      use comchem
      use comco2
      use comcouple
      use comdi
      use comdti
      use compart
      use comriv
      use comrlp, only: rlpnew, ntable, ntblines,nphases,rlp_phase,
     +  rlp_group
      use comrxni
      use comsi
      use comsk
      use comsplitts
      use comsptr
      use comwellphys
      use comwt
      use comxi, only : nmfil, cform
      use comzeoli
      use davidi
	use trxnvars

      implicit none

      logical blank
      integer ith, iscan, isorp

      logical null1, opened, done, ok, null_new
      integer idumm, ja, jb, jc, numtime
      character* 4 cdumm, macro, ctmp
      real*8 adumm, rdum1
      integer igrp,pos_index,trac_flag,ispeci
      integer ic2
      integer, allocatable :: group_mat(:,:)
      integer matnum,itemp,jj,ii,idum
      character*80 dumstring
      character*100 dummy_line, filename
      integer msg(20)
      integer nwds
      real*8 xmsg(20)
      integer imsg(20)
      character*32 cmsg(20)
      character*80 table_name
      integer kz,iz2,iz2p,i1,i2,ipivt,sehindexl,sehindexv,neqp1
      integer ireaddum, inptread, open_file
      integer locunitnum, kk, j
      integer ngdpm_models, nsize_layer
      integer, allocatable :: itemporary(:)
      integer idum1, idum2, ilines, i
      integer icount, tprp_num
      integer jjj, isimnum, realization_num,maxrp
      logical nulldum, found_end, intfile_ex

      real*8 rflag
        maxrp = 30
        if(.not.allocated (rlp_phase)) then
         allocate (rlp_phase(20, maxrp),rlp_group(20))
        endif
c**** read startup parameters ****
      rewind inpt
      iccen = 0
      ice = 0
      ico2 = 0
      icarb = 0
c gaz 112817      
      iwater_table = 0 
      icgts = 0
      idoff = 1
      ihead = 0
      irlp_fac = 0
      iflxc = 0 
      iflux_ts = 0
      istrs = 0
      isalt = 0
      iread_rcoll = 0           
      total_colloids = 0
      total_daughters = 0
      total_irrev = 0
      total_rev = 0
      maxprobdivs = 0
      ianpe = 0
      ivboun = 0
      nflxt = 0
      nspeci = 0
      numd = 0
      numsorp = 0
      numtime = 0
      ncplx = 0
      numrxn = 0
      numvcon = 0
      nriver = 0
      rlp_flag = 0
      rxn_flag = 0
      ipermstr1 = 0
      ipermstr2 = 0
      ipermstr3 = 0
      ipermstr4 = 0
      ipermstr5 = 0
      ipermstr6 = 0
      ipermstr8 = 0
      ipermstr11 = 0
      ipermstr31 = 0
      iactive = 0

c zvd 17-Aug-09 move boun flag initializations here
      iqa=0
      ixa = 0
      iha = 0
      ipha = 0
      itha = 0
      iqw=0
      iqf=0
      iqenth=0
      isatb=0
      ienth=0
      itempb=0
      ipresa=0
      ipresw=0
      imped=0
      its=0
      ifd=0
      isf=0
      ixperm = 0
      iyperm = 0
      izperm = 0
      i_init = 0
c new initial value stuff 12/3/98 GAZ
      itempb_ini=0
      ipresa_ini=0
      ipresw_ini=0
      isatb_ini=0
      icm=0
      iwght = 1
      isty = 0
      icoef_neg = 0
      nflxz = 0
      sv_hex_tet = .false.
      ipr_tets = 0
c gaz 100118 set nrlp here to enable mulpiple rlp macros      
      nrlp = 0
c zvd 17-Aug-09 move boun flag initializations here

      if(allocated(izone_free_nodes)) izone_free_nodes = 0
      if(allocated(move_type)) move_type=0
      compute_flow = .true.
      reverse_flow = .false.
      pod_flag = .false.
      omr_flag = .false.
      save_omr = .false.
      sptr_flag = .false.
      fperm_flag = .false.
c gaz 070619 variables used in newer itfc      
      nitf_use = .false.
      ncol_use = .false.

 10   continue
      filename = ''
      read (inpt, '(a80)', END = 50) dumstring
      if (dumstring(1:1) .eq. '#') go to 10
      read (dumstring, '(a4)') macro
      call parse_string2(dumstring,imsg,msg,xmsg,cmsg,nwds)
      if (nwds .gt. 1) then
c Check for "OFF" keyword for skipping macro
         found_end = .false.
         do i = 2, nwds
            if (msg(i) .eq. 3) then
               if (cmsg(i) .eq. 'off' .or. cmsg(i) .eq. 'OFF') then
                  call skip_macro (macro, inpt, found_end)
                  exit
               end if
            end if
         end do
         if (found_end) goto 10
      end if
      if (macro.eq.'stop') then
         go to 40
      else if(macro.eq.'file') then
! Read filename to ensure we skip this line and don't interpret start
! of filename as a macro
         read (inpt, '(a100)') filename
      else if (macro .eq. 'cden') then
         if (nwds .gt. 1) then
            if (cmsg(2) .eq. 'multi' .or.  cmsg(2) .eq. 'table') then
               cden_flag = 1
            else if (cmsg(2) .eq. 'co2') then
               cden_flag = 2
               if (nwds .gt. 2 .and. msg(3) .eq. 3) then
                  cden_spnam = cmsg(3)
               else
                  cden_spnam = 'Na'
               end if
            end if
         end if
      else if (macro .eq. 'cont') then
         call start_macro(inpt, locunitnum, macro)
! Read over keywords if avs, tecplot or surfer
         read (locunitnum, '(a3)', END = 50) macro
         if(macro(1:3) .eq. 'avs' .or. macro(1:3) .eq. 'tec' .or.
     &        macro(1:3) .eq. 'sur' .or. macro .eq. 'csv') then
            ok = .false.
            do
               read (locunitnum, *, END = 60) dumstring
               if (dumstring(1:1) .eq. 'e' .or. dumstring(1:1) .eq. 
     &              'E') then
                  ok = .true.
                  exit
! Check for other macros that may use 'end' keyword
               else if (dumstring(1:4) .eq. 'boun' .or. dumstring(1:4)
     &                 .eq. 'hist' .or. dumstring(1:4) .eq. 'rest' 
     &                 .or. dumstring(1:4) .eq. 'stea' .or. 
     &                 dumstring(1:4) .eq. 'rlpm' .or. dumstring(1:4) 
     &                 .eq. 'trxn') then
                  exit
               end if
            end do
 60         if (.not. ok) then
               write (ierr, 55) 'ENDAVS'
               write (ierr, 56) 'cont'
               if (iptty .ne. 0) then 
                  write (iptty, 55) 'ENDAVS'
                  write (iptty, 56) 'cont'
               end if
               stop
            end if
         end if
         call done_macro(locunitnum)
      else if(macro .eq. 'hist') then
         call start_macro(inpt, locunitnum, macro)
! Read through keywords so no conflict with actual macro names
         ok = .false.
         do
            read (locunitnum, '(a3)', END = 70) dumstring
            if (dumstring(1:3) .eq. 'end' .or. dumstring(1:3) .eq. 
     &           'END') then
               ok = .true.
               exit
! Check for other macros that may use 'end' keyword
            else if (dumstring(1:4) .eq. 'boun' .or. dumstring(1:4)
     &              .eq. 'cont' .or. dumstring(1:4) .eq. 'rest' .or.
     &              dumstring(1:4) .eq. 'stea' .or. dumstring(1:4) 
     &              .eq. 'rlpm' .or. dumstring(1:4) .eq. 'trxn') then
               exit
            end if
         end do
 70      if (.not. ok) then
            write (ierr, 55) 'END'
            write (ierr, 56) 'hist'
            if (iptty .ne. 0) then 
               write (iptty, 55) 'END'
               write (iptty, 56) 'hist'
            end if
            stop
         end if
         call done_macro(locunitnum)
      else if(macro.eq.'zone') then
         call start_macro(inpt, locunitnum, macro)
         call done_macro(locunitnum)
      else if(macro.eq.'zonn') then
         call start_macro(inpt, locunitnum, macro)
         call done_macro(locunitnum)
      else if(macro.eq.'rflo') then
         compute_flow = .false.
         backspace inpt
         read (inpt, *) dumstring
         if (dumstring(1:5) .eq. 'rflor') reverse_flow = .true.
      else if (macro.eq.'hcon') then
         idoff = -1
      else if (macro .eq. 'flxz') then
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) nflxz
c Read zones to name flxz output files in inhist if necessary
         if(.not.allocated(iflxz)) allocate(iflxz(max(1,nflxz)))
         read (locunitnum, *)(iflxz(i),i=1,nflxz)
         call done_macro(locunitnum)
      else if (macro .eq. 'iter') then
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum,*)g1, g2, g3, tmch, overf
         read (locunitnum, *) irdof, islord, iback, icoupl, rnmax
         if (ihead .eq. 1 .and. irdof .ne. 13) then
            irdof = 13
         end if
          if (jswitch .ne. 0 .and. irdof .eq. 13) then
            irdof = 0
         end if
         call done_macro(locunitnum)
      else if (macro.eq.'sol ') then
         call start_macro(inpt, locunitnum, macro)
         read(locunitnum,*) idoff, intg
         if (idoff .le. 0) idoff = -1
         call done_macro(locunitnum)
      else if (macro .eq. 'time') then
         numtime = numtime + 1
         if (numtime .gt. 1) then
            write (ierr,*) 'More than one time macro in input file'
            write (ierr,*) 'STOPPING'
            if (iptty .ne. 0) then
               write (iptty,*) 
     &              'More than one time macro in input file'
               write (iptty,*) 'STOPPING'
            end if
            stop
         end if
         call start_macro(inpt, locunitnum, macro)
         read(locunitnum, *)
 5       read(locunitnum, '(a)') dummy_line
         if (.not. null1(dummy_line)) then
            icgts = icgts + 1
            goto 5
         end if
         call done_macro(locunitnum)
      else if (macro .eq. 'phys'.or. macro .eq. 'ndph') then
c     find number of wellphysics or non darcy models models 
c     eg drift flux
         iwellp_flag = 1
         call start_macro(inpt, locunitnum, macro)

         nwellphy = 0
 111     continue
         read(locunitnum,'(a80)') wdd1
         if(.not. null1(wdd1)) then
            backspace locunitnum
            read(locunitnum,*) idumm
            backspace locunitnum
            nwellphy = nwellphy + 1
            if(idumm.eq.1.or.idumm.eq.2) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
            elseif(idumm.eq.-1) then

            end if
         else
            go to 211
         endif
         go to 111
 211     continue

         call done_macro(locunitnum)              
      else if (macro.eq.'rlp ') then
c     old model
c     find number of relative perm models 
c     use nrlp (comai.h) to transfer size of rlp arrays
c     no longer use variable ichng so value can be saved
c     check for read from other file
         rlp_flag = 1
         call start_macro(inpt, locunitnum, macro)
c gaz 100118 will set nrlp = 0 at start on routine so multiple rlp macros can be used
c         nrlp = 0
 11      continue
         read(locunitnum,'(a80)') wdd1
         if(.not. null1(wdd1)) then
            backspace locunitnum
            read(locunitnum,*) idumm
            backspace locunitnum
            ichng=ichng+1
            nrlp = nrlp + 1
            if(idumm.eq.1) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
            elseif(idumm.eq.-1) then
               read(locunitnum,*) idumm,(adumm,ja=1,4)
            elseif(idumm.eq.2) then
               read(locunitnum,*) idumm,(adumm,ja=1,4)
            elseif(idumm.eq.3) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
            elseif(idumm.eq.-4) then
               read(locunitnum,*) idumm,(adumm,ja=1,15)
            elseif(idumm.eq.4) then
               read(locunitnum,*) idumm,(adumm,ja=1,15)
            elseif(idumm.eq.5) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
            elseif(idumm.eq.6) then
               read(locunitnum,*) idumm,(adumm,ja=1,15)
            elseif(idumm.eq.7) then
               read(locunitnum,*) idumm,(adumm,ja=1,16)
            elseif(idumm.eq.8) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
            elseif(idumm.eq.9) then
               read(locunitnum,*) idumm,(adumm,ja=1,10)
            elseif(idumm.eq.10) then
               read(locunitnum,*) idumm,(adumm,ja=1,8)
            elseif(idumm.eq.11) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
            elseif(idumm.eq.12) then
               read(locunitnum,*) idumm,(adumm,ja=1,10)
            elseif(idumm.eq.13) then
               read(locunitnum,*) idumm,(adumm,ja=1,10)
c            elseif(idumm.eq.14) then
c               read(locunitnum,*) idumm,(adumm,ja=1,6)
c            elseif(idumm.eq.15) then
c               read(locunitnum,*) idumm,(adumm,ja=1,6)
c    temma add 2005/11/07
            elseif(idumm.eq.14) then
                read(locunitnum,*) idumm,(adumm,ja=1,12)
            elseif(idumm.eq.15) then
                read(locunitnum,*) idumm,(adumm,ja=1,12)
            elseif(idumm.eq.16) then
               read(locunitnum,*) idumm,(adumm,ja=1,12)
            elseif(idumm.eq.17) then
               read(locunitnum,*) idumm,(adumm,ja=1,14)
        nphases(nrlp)=3
        rlp_phase(nrlp,1)=20
        rlp_phase(nrlp,2)=21
        rlp_phase(nrlp,3)=22
        rlp_group(nrlp)=nrlp

            elseif(idumm.eq.18) then
                read(locunitnum,*) idumm,(adumm,ja=1,14)
            elseif(idumm.eq.19) then
               read(locunitnum,*) idumm,(adumm,ja=1,7)
            elseif(idumm.eq.20) then
               read(locunitnum,*) idumm,(adumm,ja=1,16)
            elseif(idumm.eq.21) then
               read(locunitnum,*) idumm,(adumm,ja=1,6)
c    temma add 2005/11/07
c            elseif(idumm.eq.14) then
c                read(locunitnum,*) idumm,(adumm,ja=1,12)
c            elseif(idumm.eq.15) then
c                read(locunitnum,*) idumm,(adumm,ja=1,12)
c            elseif(idumm.eq.18) then
c               read(locunitnum,*) idumm,(adumm,ja=1,14)
c            elseif(idumm.eq.19) then
c               read(locunitnum,*) idumm,(adumm,ja=1,14)
c            elseif(idumm.eq.20) then
c               read(locunitnum,*) idumm,(adumm,ja=1,16)
            else
c Undefined rlp model, stop
               write (ierr, 12) nrlp, idumm
               if (iout .ne. 0) write (iout, 12) nrlp, idumm
               if (iptty .ne. 0) write (iptty, 12) nrlp, idumm
               stop
            end if
         else
            go to 21
         endif
         go to 11
 12      format ('STOP: Error in rlp macro at entry', i3, 
     &        'Invalid model specified: ', i3)
 21      continue
         call done_macro(locunitnum)

      else if (macro .eq. 'rlpm') then
         call start_macro(inpt, locunitnum, macro)
         if (.not. rlpnew) nrlp = 0
         rlp_flag = 1
         rlpnew = .true.
         idum2 = 0
         do
            read (locunitnum, '(a80)') dumstring
            if (null1(dumstring) .or. dumstring(1:3) .eq. 'end' .or. 
     &           dumstring(1:3) .eq. 'END') then
               exit
            else if (dumstring(1:5) .eq. 'group' .or. dumstring(1:5)
     &              .eq. 'GROUP') then
               nrlp = nrlp + 1
            else if (dumstring(1:3) .eq. 'tab' .or. dumstring(1:3)
     &              .eq. 'TAB') then
c Count number of tables that will be read, 
c     and number of entries in each table
               ntable = ntable + 1
               read (locunitnum, '(a80)') dumstring
               if (dumstring(1:4) .eq. 'file') then
                  read (locunitnum, '(a80)') filename
                  idum = open_file (filename, 'old')
                  do
                     read (idum, '(a80)', end = 24) dumstring
                     if (null_new(dumstring) .or.  dumstring(1:3) .eq. 
     &                    'end' .or. dumstring(1:3) .eq. 'END') exit
                     call parse_string2(dumstring, imsg, msg, xmsg, 
     &                    cmsg, nwds)
c Don't count header lines (header lines should start with a character)
                     if (msg(1) .ne. 3) ntblines = ntblines + 1
                  end do
 24               close (idum)
               else
                  backspace (locunitnum)
                  do
                     read (locunitnum, '(a80)') dumstring
                     if (null_new(dumstring) .or.  dumstring(1:3) .eq. 
     &                    'end' .or. dumstring(1:3) .eq. 'END') exit
                     ntblines = ntblines + 1
                  end do
               end if
            end if
         end do
         call done_macro(locunitnum)

      else if (macro.eq.'boun') then
c     find number of boun models 
c     
c     check for read from other file
         call start_macro(inpt, locunitnum, macro)

 31      continue
c         read(locunitnum,'(a80)') wdd1
c         if(.not. null1(wdd1)) then
c            backspace locunitnum
         read(locunitnum,'(a80)') wdd1(1:80)
         if(wdd1(1:4).eq.'mode') then
            maxmodel=maxmodel+1
         else if(wdd1(1:2).eq.'cy') then
            read(locunitnum,*) idumm
            maxtimes=max(maxtimes,idumm)
         else if(wdd1(1:2).eq.'ti') then
            read(locunitnum,*) idumm, rdum1
            if (rdum1 .ne. 0.0) then
               maxtimes=max(maxtimes,idumm+1)
            else
               maxtimes=max(maxtimes,idumm)
            end if
         else if(wdd1(1:3).eq.'huf') then
            iha=1 
         else if(wdd1(1:2).eq.'hu') then
            iha=1             
         else if(wdd1(1:2).eq.'ph') then
            ipha=1 
         else if(wdd1(1:2).eq.'th') then
            itha=1 
         else if(wdd1(1:3).eq.'chm') then
            icm=1
         else if(wdd1(1:2).eq.'sa') then
            iqa=1
         else if(wdd1(1:3).eq.'fxa') then
            ixa=1
         else if(wdd1(1:3).eq.'swf') then
            iqf=1
         else if(wdd1(1:2).eq.'sw') then
            iqw=1
         else if(wdd1(1:4).eq.'sco2') then
            iqco2=1
         else if(wdd1(1:3).eq.'sfh') then
            isf=1
         else if(wdd1(1:2).eq.'sf') then
            isf=1
         else if(wdd1(1:2).eq.'fd') then
            ifd=1
         else if(wdd1(1:3).eq.'dfd') then
            ifd=1
         else if(wdd1(1:2).eq.'se') then
            iqenth=1
         else if(wdd1(1:3).eq.'dsa') then
            iqa=1
         else if(wdd1(1:3).eq.'dsw') then
            iqw=1
         else if(wdd1(1:5).eq.'dsco2') then
            iqco2=1
         else if(wdd1(1:3).eq.'dse') then
            iqenth=1
         else if(wdd1(1:2).eq.'si') then
            isatb_ini=1
         else if(wdd1(1:2).eq.'s ') then
            isatb=1
         else if(wdd1(1:3).eq.'pwi') then
            ipresw_ini=1
         else if(wdd1(1:3).eq.'pwo') then
            ipresw=1
            isatb=1
         else if(wdd1(1:2).eq.'pw') then
            ipresw=1
            isatb=1
         else if(wdd1(1:3).eq.'pai') then
            ipresa_ini=1
         else if(wdd1(1:3).eq.'pao') then
            ipresa=1
            isatb=1
         else if(wdd1(1:2).eq.'pa') then
            ipresa=1
            isatb=1
         else if(wdd1(1:4).eq.'tran') then
            isty=1
         else if(wdd1(1:3).eq.'tmi') then
            itempb_ini=1
         else if(wdd1(1:2).eq.'tc') then
            itempb2=1
         else if(wdd1(1:2).eq.'ts') then
            its=1
         else if(wdd1(1:2).eq.'t') then
            itempb=1
         else if(wdd1(1:2).eq.'hd') then
            itempb=1
         else if(wdd1(1:3).eq.'hdo') then
            itempb=1
         else if(wdd1(1:1).eq.'h') then
            isatb=1
         else if(wdd1(1:2).eq.'hu') then
            iha=1 
         else if(wdd1(1:2).eq.'ph') then
            ipha=1 
         else if(wdd1(1:2).eq.'th') then
            itha=1 
         else if(wdd1(1:4).eq.'wgtp') then
            iwght=2
         else if(wdd1(1:4).eq.'wgtu') then
            iwght=3
         else if(wdd1(1:3).eq.'wgt') then
            iwght=1
         else if(wdd1(1:2).eq.'if') then
            imped=1
         else if(wdd1(1:2).eq.'en' .and. wdd1(1:3).ne.'end') then
            ienth=1
         else if(wdd1(1:2).eq.'ft') then
            ienth=1
          else if(wdd1(1:3).eq.'fen') then
            ienth=1 
         else if(wdd1(1:2).eq.'kx') then
            ixperm=1
         else if(wdd1(1:2).eq.'ky') then
            iyperm=1
         else if(wdd1(1:2).eq.'kz') then
            izperm=1
         else if(wdd1(1:20).eq.'                    ') then
            if(null1(wdd1)) go to 41
         else if(wdd1(1:3).eq.'end') then
            go to 41
         endif
         go to 31
c         end if
 41      continue

         call done_macro(locunitnum)

      else if (macro .eq. 'alti')  then
c**** We need to know number of nodes if alternate geometry data is used ****
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) cdumm, neq
         call done_macro(locunitnum)

      else if (macro(1:3) .eq. 'nap' .or. macro .eq. 'szna')  then
c     need to know if napl-water  is envoked
         ico2 = -3  
c gaz 112817         
      else if (macro(1:3)  .eq.  'eos')  then
c     need to know if a table with water props is read
         call start_macro(inpt, locunitnum, macro)
          read(locunitnum,'(a80)') wdd1(1:80)
          if(wdd1(1:5).eq.'table') iwater_table = 1   
          table_name = trim(wdd1(7:80))
          inquire(file=table_name, exist=intfile_ex)
          if(intfile_ex) iwater_table = 2    
          if(iwater_table.eq.2) nmfil(31) = table_name
         call done_macro(locunitnum)
         
      else if (macro(1:3)  .eq.  'air')  then
c     need to know if air-water  is envoked
         ico2 = -2 
c gaz 110819 reading tref, pref here (these variables are now global          
         call start_macro(inpt, locunitnum, macro)
          read(locunitnum,'(a80)') wdd1(1:80)
          read(locunitnum,*) tref, pref
         call done_macro(locunitnum)           
      else if (macro .eq. 'ice ' .or. macro .eq. 'meth') then
         ice = 1
         
c RJP added carbon flag
      else if (macro .eq. 'carb') then
         icarb = 1
         
      else if (macro  .eq.  'boun')  then
c     need to know if boussinesq approx is envoked
         iboun = 1           
         
      else if (macro  .eq.  'head')  then
c     need to know if head input instead of pressures
c airwater isothermal
c bous automatically enabled
         ihead = 1           
         icons = 100000
         irdof = 13

      else if (macro  .eq.  'wtsi')  then
c need to know if simplified water table is invoked
         ifree = 1
      else if (macro .eq. 'salt') then
c     DRH 1/2/2013: need to know if modeling salt
         isalt = 1
      else if (macro  .eq.  'coor')  then
c**** We need to know number of nodes if coor is in inpt ****
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) neq
         call done_macro(locunitnum)
         
      else if (macro  .eq.  'elem')  then
c**** We need to know number of elements  if elem is in inpt ****
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) idumm,nei 
         call done_macro(locunitnum)
         
      else if (macro  .eq.  'ctrl')  then
         call start_macro(inpt,locunitnum, macro)
         read (locunitnum, *) idumm, rdum1, north
         maxor = north + 1
         igauss = 1
         
 20      continue
         read (locunitnum, '(a80)') wdd1
         if (null1(wdd1)) go to 30
         backspace locunitnum
         read (locunitnum, *) ja, jb, jc, idumm
         igauss = max (igauss, idumm)
         if(ja  .eq.  0) go to 30
         go  to  20
         
 30      read (locunitnum, *)  adumm
         read (locunitnum, *)  idumm
c**** We need to know problem geometry ****
         read (locunitnum, *)  icnl
         
         call done_macro(locunitnum)
c**** We need to know if dual porosity or 
c**** dual porosity/dual permeability problem ****
      else if(macro .eq. 'dual') then
         idualp = 1
         
      else if (macro .eq. 'fper') then
c need to know if fperms are being used
         fperm_flag = .true.

      else if (macro .eq. 'frlp') then
c need to know if rel perm factor is used
         irlp_fac = 1
      else if (macro .eq. 'flxo') then
         iflxc = iflxc +1
         call start_macro(inpt, locunitnum, macro)
         read(locunitnum,*) nflx
          ivelo = 1
         nflxt = nflxt + nflx
         call done_macro(locunitnum)
      else if (macro .eq. 'dvel') then
         iflxc = iflxc +1
         call start_macro(inpt, locunitnum, macro)
         read(locunitnum,*) nflx
          ivelo = -1
         nflxt = nflxt + nflx
         call done_macro(locunitnum)
      else if (macro .eq. 'dpdp') then
         idpdp = 1

c	Section for determining Generalized Dual Porosity Model (GDPM)
c	and Generalized Dual Permeability Model (GDKM) parameters
         
      else if(macro .eq. 'gdpm' .or. macro .eq. 'gdkm') then
         if(macro .eq. 'gdkm') then
          gdkm_flag = 1
          backspace locunitnum
          read(locunitnum,'(a80)') dumstring
          gdkm_new = .false.
c----------------------------------------
c    Shaoping add, 10/23/2017
c         do i = 1, 80
          do i = 1, 78
c----------------------------------------
           if(dumstring(i:i+2).eq.'new') then
             gdkm_new =.true.
             go to 690
           endif
          enddo
         endif
690     call start_macro(inpt, locunitnum, macro)
         if(.not.gdkm_new) then
           read (locunitnum, *) gdpm_flag, ngdpmnodes
         else
            ngdpmnodes = -999 
         endif
         
         if(gdkm_flag.eq.1) then
          gdpm_flag = 1
         endif
         maxgdpmlayers = 0
         ngdpm_models = 0

c	Code scans the file to determine the number of models and the
c	maximum number of node points in the matrix so that array
c	sizes can be set and arrays allocated in allocmem
         
 1000    continue
         read(locunitnum,'(a80)') dumstring
         if (.not.null1(dumstring)) then
            backspace locunitnum
            ngdpm_models = ngdpm_models + 1
c gaz 091516            
            if(gdkm_flag.eq.0) then
             read(locunitnum,*) nsize_layer,adumm,(adumm,i=1,
     &      nsize_layer)
            else
             read(locunitnum,*) nsize_layer,adumm
             nsize_layer = 1 
            endif
            maxgdpmlayers = max(nsize_layer,maxgdpmlayers)
            goto 1000
         end if

c	At this point, we have the number of models in the gdpm
c	formulation (ngdpm_models), and the maximum number of layers
c       in any of the models

c     Allocate arrays needed for gdpm

         ngdpm_models = max(1,ngdpm_models)
         maxgdpmlayers = max(1,maxgdpmlayers)
         if(.not.allocated(ngdpm_layers)) then
            allocate(ngdpm_layers(0:ngdpm_models))
            allocate(vfrac_primary(0:ngdpm_models))
            allocate(gdpm_x(0:ngdpm_models,maxgdpmlayers))
            allocate(gdpm_vol(0:ngdpm_models,maxgdpmlayers))
            allocate(wgt_length(ngdpm_models))
            ngdpm_layers = 0
            vfrac_primary = 0.
            gdpm_x = 0.
         end if
         if(gdkm_flag.ne.0) then
          if(.not.allocated(gdkm_dir)) then  
           allocate(gdkm_dir(0:ngdpm_models))
           gdkm_dir = 0
          endif
         endif

         call done_macro(locunitnum)
         
      else if (macro .eq. 'rare') then
c**** over write area/d terms and volumes ****
c**** can be used with gdkm and gdpm models *****
c**** need to count entries
         call start_macro(inpt, locunitnum, macro)         
         icoef_replace = 1
         call coef_replace_ctr(-1)   
         call done_macro(locunitnum)         
c     Section for determining Element enrichment
      else if(macro .eq. 'enri') then
  
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) enri_flag, nenrinodes
         maxenrichlayers = 0
         ienrich_models = 0

c	Code scans the file to determine the number of models and the
c	maximum number of node points in the matrix so that array
c	sizes can be set and arrays allocated in allocmem
         
 1010    continue
         read(locunitnum,'(a80)') dumstring
         if (.not.null1(dumstring)) then
            backspace locunitnum
            ienrich_models = ienrich_models + 1
            read(locunitnum,*) nsize_layer,adumm,(adumm,i=1,
     &      nsize_layer)
            maxenrichlayers = max(nsize_layer,maxenrichlayers)
            goto 1010
         end if
c
c     Allocate arrays needed for enri
c

         call done_macro(locunitnum) 
         
c     RJP 12/14/06 added following for river/wellbore	      
c     Section for determining implicit river or well model (riv )
c     parameters
         
      else if(macro .eq. 'rive'.or. macro .eq. 'well') then
         
         call start_macro(inpt, locunitnum, macro)
 785     read (locunitnum, '(a9)') dumstring(1:9)
         if (dumstring(1:1) .eq. '#') go to 785
         if (dumstring(1:9) .eq. 'wellmodel') then
            read (locunitnum, *) nriver,iriver

            if(iriver.eq.1) then
               npoint_riv = 0
               do i = 1,nriver
                  read (locunitnum, *) idum1,ii,jj		      
                  idum = ii 
                  do kk = 1,jj
                     npoint_riv = npoint_riv + 1
                     do j = 1,idum
                        npoint_riv = npoint_riv + 1
                     enddo
                  enddo
                  do ii = 1,jj+1
                     read(locunitnum,'(a80)') wdd1
                     if(null1(wdd1)) goto 441
                     backspace locunitnum
                     read (locunitnum, *) idum1
                     if (idum1 .lt. 0) goto 441
                  enddo
 441              continue
               enddo
            else if(iriver.eq.2) then
c	 
               allocate(coor_dum(nriver,3))
               do i = 1,nriver
                  read (locunitnum, *) ii, 
     &                 coor_dum(ii,1),coor_dum(ii,2),coor_dum(ii,3)
               enddo
               read(locunitnum,*) n_well_seg
	       allocate(iwell_seg(n_well_seg,2))
	       allocate(well_rad(n_well_seg))
	       allocate(well_dz(n_well_seg))
               do i = 1,n_well_seg
                  read(locunitnum,*) j,iwell_seg(j,1),iwell_seg(j,2),
     &                 well_rad(j),well_dz(j)	       
               enddo
               npoint_riv = 0
               max_seg_div = 0
               do j = 1,n_well_seg
                  ii = iwell_seg(j,1)
                  kk = iwell_seg(j,2)
                  adumm = (coor_dum(ii,1)-coor_dum(kk,1))**2 +
     &                 (coor_dum(ii,2)-coor_dum(kk,2))**2 + 
     &                 (coor_dum(ii,3)-coor_dum(kk,3))**2 
                  adumm = sqrt(adumm)
                  jj = adumm/well_dz(j) + 1 
                  max_seg_div = max(max_seg_div,jj)
                  npoint_riv = npoint_riv + jj
                  nnelm_riv = nnelm_riv + (jj-1)
               enddo
               deallocate(coor_dum,iwell_seg,well_rad,well_dz)
            endif
         endif
         call done_macro(locunitnum) 	        
      else if (macro .eq. 'init') then
        i_init = 1  
      else if (macro .eq. 'itfc') then
         call start_macro(inpt, locunitnum, macro)
         opened = .false.
         interface_flag = 1
         nitfcpairs = 0
         ncoldpairs = 0
         nitfcitfc = 0
         nitfcsizes = 0
c     Flow interface portion of input
         kk = 0
 6000    continue
         read(locunitnum,'(a80)') dummy_line
         if(null1(dummy_line)) then
            goto 6001
         else
            nitfcpairs = nitfcpairs + 1
         end if
         goto 6000
c     Transport part of interface input
 6001    continue
         read(locunitnum,'(a80)') dummy_line
         if(null1(dummy_line)) goto 6002
 6003    continue
         read(locunitnum,'(a80)') dummy_line
         if(null1(dummy_line)) then
            goto 6002
         else
            kk = kk+1
            backspace locunitnum
            read(locunitnum,*) idum1, idum2, rflag
            if(rflag.lt.0) then
               nitfcitfc = nitfcitfc + 1
               ilines = 0
               read(locunitnum,'(a80)') dummy_line
               if(dummy_line(1:4).eq.'file') then
                  read(locunitnum,'(a100)') filename
                  inptread = open_file(filename,'old')
                  opened = .true.
               else
                  backspace locunitnum
                  inptread = locunitnum
                  opened = .false.
               end if
 7001          continue
               read(inptread,'(a80)') dummy_line
               if(null1(dummy_line)) then
                  goto 7002
               else
                  ilines = ilines + 1
                  goto 7001
               end if
            end if
 7002       continue
            if(opened) then
               close (inptread)
               opened = .false.
            end if
            nitfcsizes = max(ilines,nitfcsizes)
         end if
         goto 6003
 6002    continue
         if(nitfcpairs.gt.0) nitf_use = .true.
         if(kk.gt.0) ncol_use = .true.
         nitfcpairs = max(1,nitfcpairs)
         ncoldpairs = max(1,kk)
         call done_macro(locunitnum)

      else if (macro .eq. 'ngas') then
c     need to know if noncondensible present
         ico2 = 3
         
      else if (macro .eq. 'trac') then
         iccen = 1
         call start_macro(inpt, locunitnum, macro)
         
c---------  Hari added  6/10/04
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:3).eq.'rip') then
            read(locunitnum,'(a80)') dummy_line
         else
            backspace(locunitnum)
         end if
c---------------------- Hari         
 
         read (locunitnum, '(a5)') dumstring
         if (dumstring(1:5) .eq. 'userc') then
c check for file keyword, and read past filename if present
            read (locunitnum, '(a4)') dumstring
            if (dumstring(1:4) .eq. 'file')   then
               read (locunitnum, *)
            else
               backspace locunitnum
            end if
         end if
         read (locunitnum, *) 
         read (locunitnum, *) 
         read (locunitnum, '(a4)')dumstring
         if(dumstring(1:4).eq.'tpor')then
            tpor_flag = .TRUE.
            do
               read(locunitnum,'(a80)') dumstring
               if (null1(dumstring)) exit
            enddo
         else
            tpor_flag = .FALSE.
            backspace locunitnum
         endif
         read (locunitnum, *) nspeci
          if(.not.allocated(species)) allocate (species(nspeci))
          species=' '
         trac_flag=1
         read(locunitnum,'(a4)') dumstring
         if (dumstring(1:4).eq.'ldsp') then
! longitudinal and transverse dispersion coefficients keyword was read
! read another line
            read(locunitnum,'(a4)') dumstring
         endif
         if ((dumstring(1:4).ne.'dspl').and.(dumstring(1:4)
     &        .ne.'dspb').and.(dumstring(1:4).ne.'dspb')) then
! Same dispersivity and diffusivity are not used for all liquid / vapor 
! calculations
            backspace locunitnum
         else
            do
               read(locunitnum,'(a80)') dumstring
               if (null1(dumstring)) exit
! Count number of regions with same dispersivity and diffusivity
               numd = numd + 1
            enddo
            do
               read(locunitnum,'(a80)') dumstring
               if (null1(dumstring)) exit
            enddo
         endif
         ncpnt = 0
         nimm = 0
         nvap = 0
         do ii = 1,nspeci
            read(locunitnum,*) ispeci
            if(ispeci.eq.1.or.abs(ispeci).eq.2)then
               ncpnt = ncpnt + 1
            elseif(ispeci.eq.0)then
               nimm = nimm + 1
            else
               nvap = nvap + 1
            endif
            isorp = 0
            if (ispeci.eq.0) then
! Solid, do nothing     
               continue
            else   
               do
                  read(locunitnum,'(a80)') dumstring
                  if(null1(dumstring)) exit
! Count number of sorption models for this species
                  isorp = isorp + 1
! Liquid or vapor, abs(ispeci).eq.1, do nothing else
                  if (abs(ispeci).eq.2) then
! Henry's law species read additional input line if data is on two lines
                     call parse_string(dumstring, imsg, msg, xmsg, 
     &                    cmsg, nwds)
                     if ((numd .eq. 0 .and. nwds .lt. 16) .or. 
     &                    (numd .gt. 0 .and. nwds .lt. 8)) then
                        read(locunitnum,*)
                     end if
                  endif
               end do
               do
                  read(locunitnum,'(a80)') dumstring
                  if (null1(dumstring)) exit
               enddo
               if (abs(ispeci).eq.2) then
! Henry's law species, model parameters
                  read (locunitnum,*) idum
                  if (idum.eq.1 .or. idum.eq.2 .or. idum.eq.3) then
                     continue
                  else
                     write(ierr,*)' ** Using Old Input '
                     write(ierr,*)' Enter Temperature Dependency '
                     write(ierr,*)' Model Number: 1 - Van Hoff '
                     write(ierr,*)' 2 - awwa model, see manual'
                     write(ierr,*)' for details **'
                     if (iptty .gt. 0) then
                        write(iptty,*)' Using Old Input '
                        write(iptty,*)' Enter Temperature Dependency '
                        write(iptty,*)' Model Number: 1 - Van Hoff '
                        write(iptty,*)' 2 - awwa model, see manual'
                        write(iptty,*)' for details **'
                     end if
                     stop
                  endif
               endif
            endif
            do
               read(locunitnum,'(a80)') dumstring
               if (null1(dumstring)) exit
            enddo
            do
               read(locunitnum,'(a80)') dumstring
               if (null1(dumstring)) exit
            enddo
! maximum number of sorption models for any species
            numsorp = max (isorp, numsorp)
         enddo
          
         call done_macro(locunitnum)

      else if(macro(1:3) .eq. 'rxn') then
         call start_macro( inpt, locunitnum, macro)
         rxn_flag = 1
         read(locunitnum,*)
         read(locunitnum,*)ncplx,numrxn
         read (locunitnum, *)
         read (locunitnum, *)ngroups
         allocate(group_mat(ncpnt,ncpnt))
         group_mat = 0
         if(irun.eq.1) then
            if (ngroups .ne. 0) then
               allocate(group(ngroups,ncpnt),pos(ngroups,ncpnt))
               allocate(n_couple_species(ngroups),fzero(ngroups))
            end if
         end if
         group = 0
         pos = 0
         n_couple_species = 0
         do igrp = 1, ngroups
            read(locunitnum, *)(group(igrp,ic),ic = 1,ncpnt)
         enddo
c     add grouping stuff 8/1/94
         mdof_sol = 0
         do igrp = 1, ngroups
            pos_index = 0
            do ic = 1, ncpnt
               if(group(igrp,ic).ne.0)then
                  do ic2 = 1, ncpnt
                     if(group(igrp,ic2).ne.0)then
                        group_mat(ic,ic2)=1
                     else
                        group_mat(ic,ic2)=0
                     endif
                  enddo
                  pos_index = pos_index + 1
                  pos(igrp,pos_index)= ic
                  n_couple_species(igrp)=n_couple_species(igrp)+1
                  mdof_sol = max(n_couple_species(igrp),
     2                 mdof_sol)
               endif
            enddo
         enddo
         dimdrc = 0
         matnum = 0
         do ic = 1, ncpnt
            itemp = 0
            do ic2 = 1,ncpnt
               matnum = matnum + 1
               if(group_mat(ic,ic2).ne.0)then
                  dimdrc=dimdrc+1
                  matpos(matnum)=dimdrc
                  itemp = itemp + 1
                  drcpos(ic,itemp)=ic2
               endif
            enddo
            nderivs(ic)=itemp
         enddo
         call done_macro(locunitnum)
         
      else if (macro .eq. 'trxn') then
         iccen = 1
         call start_macro(inpt, locunitnum, 'trxn')
	inpttmp = inpt
	inpt = locunitnum
         call trxninit
	inpt = inpttmp
	 call done_macro(locunitnum)
         if(rxn_flag .eq. 0) then
            if (iout .ne. 0) write (iout, 8000)
            if (iptty .ne. 0) write(iptty, 8000)
         endif
 8000    format ('Reactions disabled for this simulation.')
         if(rxn_flag .eq. 1) then
            allocate(group_mat(ncpnt, ncpnt))
		group_mat = 0
            do igrp = 1, ngroups
               pos_index = 0
               do ic = 1, ncpnt
                  if(group(igrp,ic).ne.0)then
                     do ic2 = 1, ncpnt
                        if(group(igrp,ic2).ne.0)then
                           group_mat(ic,ic2)=1
                        else
                           group_mat(ic,ic2)=0
                        endif
                     enddo
                     pos_index = pos_index + 1
                     pos(igrp,pos_index)= ic
                     n_couple_species(igrp)=n_couple_species(igrp)+1
                     mdof_sol = max(n_couple_species(igrp),
     2                    mdof_sol)
                  endif
               enddo
            enddo
            dimdrc = 0
            matnum = 0
            do ic = 1, ncpnt
               itemp = 0
               do ic2 = 1,ncpnt
                  matnum = matnum + 1
                  if(group_mat(ic,ic2).ne.0)then
                     dimdrc=dimdrc+1
                     matpos(matnum)=dimdrc
                     itemp = itemp + 1
                     drcpos(ic,itemp)=ic2
                  endif
               enddo
               nderivs(ic)=itemp
            enddo
         endif
         
      else if (macro .eq. 'mptr') then
         ptrak=.true.
         ctmp=macro
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) nspeci,maxlayers,max_particles
c     Line with pout
         read(locunitnum,*) idum1
c     Read past group with tcurve
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:6).eq.'tcurve') then
            read(locunitnum,'(a80)') dummy_line
            read(locunitnum,'(a80)') dummy_line
         else
            backspace(locunitnum)
         endif  
  
c     Read past group with zptr
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).eq.'zptr') then
            read(locunitnum,*) ipzone
            allocate(itemporary(ipzone))
            read(locunitnum,*) (itemporary(i),i=1,ipzone)
            deallocate(itemporary)
         else
            backspace(locunitnum)
         end if
c     read rseed line
         read(locunitnum,*) idum1

c*** water table rise modification
c     Read past group with wtri
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).eq.'wtri') then
            read(locunitnum,*)water_table
            wtrise_flag = .true.
         else
            wtrise_flag = .false.
            water_table = -1.d+10
            backspace(locunitnum)
         end if
c*** water table rise modification

c     read daycs line
         read(locunitnum,*) rdum1
c     read file name if this option is used
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).eq.'file') then
            read(locunitnum,'(a80)') dummy_line
         else
            backspace(locunitnum)
         end if
c     Read afm if it exists
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:3).ne.'afm') then
            backspace(locunitnum)
         end if

c zvd 02/28/07 Add for free water diffusion coefficient and tortuosity
c     Read dfree keyword if it exists
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).ne.'dfre') then
            backspace(locunitnum)
         end if

c     read until next black line to get past dispersivities
         blank = .false.
         do while(.not.blank)
            read(locunitnum,'(a80)') dummy_line
            if(null1(dummy_line)) blank = .true.
         end do
c     read until next blank line to get past itrc lines
         blank = .false.
         do while(.not.blank)
            read(locunitnum,'(a80)') dummy_line
            if(null1(dummy_line)) blank = .true.
         end do
c     Loop through each species, read ith, determine if the
c     particles have a size distribution
         ncolspec = 0
         ncolsizes = 0
         do i = 1, nspeci
            read(locunitnum,*) ith
c     Ready to look for size distribution, count the number
c     of species with size distributions, count the number
c     of entries, record the maximum number of entries

            read(locunitnum,'(a80)') dummy_line
            if(dummy_line(1:4).eq.'size') then
               ncolspec = ncolspec + 1
               jj = 0
 2991          continue
               read(locunitnum,'(a80)') dummy_line
               if (null1(dummy_line)) then
                  goto 2992
               else
                  jj = jj + 1
                  goto 2991
               end if
 2992          continue

               ncolsizes = max(jj,ncolsizes)
            else
               backspace(locunitnum)
            end if

CHari 01-Nov-06 read in parameters for the colloid diversity model
c find max_probdivs

         read(locunitnum,'(a80)') wdd1
         if(wdd1(1:4).eq.'dive') then
! Flag for non-colloid daughter products
            read(locunitnum,*) idum
            if (idum.eq.0) then
               total_colloids = total_colloids + 1
            else
               total_daughters = total_daughters + 1
               goto 18933
            endif
            read(locunitnum,'(a80)') dummy_line
            if (dummy_line(1:4) .eq. 'file') then
               do icount = 1, 80
                  filename(icount:icount) = ' '
               enddo
               read(locunitnum,'(a80)') filename
               if(total_colloids.eq. 1)then
                  iread_rcoll = open_file(filename,'old')
               endif
            else
               backspace locunitnum
            end if
            read(locunitnum,'(a80)') dummy_line
            call parse_string(dummy_line,imsg,msg,xmsg,cmsg,nwds)
            tprp_num = imsg(1)
! Use multiple simulation number - irun for GoldSim or msim run
            isimnum = int(irun+0.01)           
            if(isimnum<=0) isimnum = 1
            if (ripfehm .eq. 0 .and. nwds .ge. 2) isimnum = imsg(2)
!            if(tprp_num.ge.11.and.tprp_num.le.14.and.
!     2           total_colloids.eq.1) then
            if(tprp_num.ge.11.and.tprp_num.le.14) then
               if(tprp_num.lt.13) then
                  if(iread_rcoll .eq. 0) then
                     write (ierr, 18929)
                     if (iout .ne. 0) write (iptty, 18929)
                     if (iptty .ne. 0) write (iptty, 18929)
                     stop
                  end if
18929             format ('STOP - CDF table data must be input using ',
     &                    'an external data file -- check mptr macro.')
                  read(iread_rcoll,'(a80)',end=18928) dummy_line
                  read(iread_rcoll,*,end=18928)idum
                  do while (idum.ne.isimnum)
                     nulldum=.false.
                     do while(.not.nulldum)
                        read(iread_rcoll,'(a80)',end=18928)dummy_line
                        nulldum=null1(dummy_line)
                     enddo
                     read(iread_rcoll,*,end=18928)idum
                  enddo
18928             if (idum .ne. isimnum) then
                     write(ierr,18930) tprp_num, isimnum,
     &                    trim(filename)
                     if (iout .ne. 0) write(iout, 18930) tprp_num,
     &                    isimnum, trim(filename)
                     if (iptty .ne. 0) write(iptty, 18930) tprp_num,
     &                    isimnum, trim(filename)
                     stop
                  end if
18930             format ('STOP - Error in inmptr for tprpflag = ', i2,
     &                 ' realization_num ', i4, ' not found in ', a)
                  jjj=0
                  read(iread_rcoll,'(a80)',end=18931) dummy_line
                  nulldum=null1(dummy_line)
                  do while(.not.nulldum)
                     jjj=jjj+1
                     read(iread_rcoll,'(a80)',end=18931) dummy_line
                     nulldum=null1(dummy_line)
                  enddo
18931             if (jjj .lt. 2) then

                     write(ierr,18932) tprp_num, isimnum,
     &                    trim(filename)
                     if (iout .ne. 0) write(iout, 18932) tprp_num,
     &                    isimnum, trim(filename)
                     if (iptty .ne. 0) write(iptty, 18932) tprp_num,
     &                    isimnum, trim(filename)
                     stop

                  end if
18932             format ('STOP - Error in inmptr for tprpflag = ',
     &                 i2, ' realization_num ', i4,
     &                 ' data not found in ', a)

                  maxprobdivs = max(jjj,maxprobdivs)
                  rewind(iread_rcoll)
               else
                  if(iread_rcoll .eq. 0) read(locunitnum,'(a80)') wdd1
               endif
            endif

            read(locunitnum,'(a80)') wdd1
cHari setup pointer arrays for reversible and irreversible colloids
            if(wdd1(1:4).eq.'irre')then
               total_irrev = total_irrev +1
            elseif(wdd1(1:4).eq.'reve')then
               total_rev = total_rev+1
            else
               write(ierr,*)'must specify irre or reve in mptr'
               if(iout.ne.0)then
                  write(iout,*)'must specify irre or reve in mptr'
               endif
               if(iptty.ne.0)then
                  write(iptty,*)'must specify irre or reve in mptr'
               endif

            endif
         else
            backspace(locunitnum)
         endif



c	read past transport parameter information
18933    continue
         read(locunitnum, *) idum1
         do iscan = 1, idum1
            read(locunitnum,'(a80)') dummy_line
         end do
c	read past source information
         read(locunitnum, *) idum1
         do iscan = 1, idum1
c     read until next black line to get past current source
            blank = .false.
            do while(.not.blank)
               read(locunitnum,'(a80)') dummy_line
               if(null1(dummy_line)) blank = .true.
            end do
         end do
c	Ready to loop back to the next species
            
         end do

         call done_macro(locunitnum)

      else if (macro .eq. 'ptrk') then
	   ptrak=.true.
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum, *) max_particles
         nspeci = 1
c         maxlayers = n0
c     Determine actual number of transport maylayers there are

c*** water table rise modification
c     Read past group with wtri
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).eq.'wtri') then
            read(locunitnum,*)water_table
         else
            backspace(locunitnum)
         end if
c*** water table rise modification

c    Read 'rip' string  
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:3).eq.'rip') then
            read(locunitnum,'(a80)') dummy_line
         end if
         read(locunitnum,'(a80)') dummy_line

         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:6).eq.'tcurve') then
            read(locunitnum,'(a80)') dummy_line
            read(locunitnum,'(a80)') dummy_line
         else
            backspace(locunitnum)
         endif  
  
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).eq.'zptr') then
            read(locunitnum,*) ipzone
            allocate(itemporary(ipzone))
            read(locunitnum,*) (itemporary(i),i=1,ipzone)
            deallocate(itemporary)
         else
            backspace(locunitnum)
         end if
         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).eq.'file') then
            read(locunitnum,'(a80)') dummy_line
         else
            backspace(locunitnum)
         end if

c     Get past afm keyword if it exists

         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:3).ne.'afm') then
            backspace(locunitnum)
         end if

c zvd 02/28/07 Add for free water diffusion coefficient and tortuosity
c     Get past dfree keyword if it exists

         read(locunitnum,'(a80)') dummy_line
         if(dummy_line(1:4).ne.'dfre') then
            backspace(locunitnum)
         end if

c     Look for colloid size information

        read(locunitnum,'(a80)') dummy_line
        if(dummy_line(1:4).eq.'size') then
           ncolspec = 1
           jj = 0
1991       continue
           read(locunitnum,'(a80)') dummy_line
           if (null1(dummy_line)) then
              goto 1992
           else
              jj = jj + 1
              goto 1991
           end if
1992       continue

           ncolsizes = jj
        else
           backspace(locunitnum) 
        end if

CHari 01-Nov-06 read in parameters for the colloid diversity model
c find max_probdivs

         read(locunitnum,'(a80)') wdd1
         if(wdd1(1:4).eq.'dive') then
            total_colloids = 1

            read(locunitnum,'(a80)') dummy_line
            if (dummy_line(1:4) .eq. 'file') then
               do icount = 1, 80
                  filename(icount:icount) = ' '
               enddo
               read(locunitnum,'(a80)') filename
               if(total_colloids.eq. 1)then
                  iread_rcoll = open_file(filename,'old')
               endif
            else
               backspace locunitnum
            end if
            read(locunitnum,'(a80)') dummy_line
            call parse_string(dummy_line,imsg,msg,xmsg,cmsg,nwds)
            tprp_num = imsg(1)
! Use multiple simulation number - irun for GoldSim or msim run
            isimnum = int(irun+0.01)
            if(isimnum<=0) isimnum = 1
            if (ripfehm .eq. 0 .and. nwds .ge. 2) isimnum = imsg(2)
            if(tprp_num.ge.11.and.tprp_num.le.14.and.
     2           total_colloids.eq.1) then
               maxprobdivs = 1
               if(tprp_num.lt.13) then
                  if(iread_rcoll .eq. 0) then
                     write (ierr, 17929)
                     if (iout .ne. 0) write (iptty, 17929)
                     if (iptty .ne. 0) write (iptty, 17929)
                     stop
                  end if
17929             format ('STOP - CDF table data must be input using ',
     &                    'an external data file -- check ptrk macro.')
                  read(iread_rcoll,'(a80)',end=17928) dummy_line
                  read(iread_rcoll,*,end=17928)idum
                  do while (idum.ne.isimnum)
                     nulldum=.false.
                     do while(.not.nulldum)
                        read(iread_rcoll,'(a80)',end=17928)dummy_line
                        nulldum=null1(dummy_line)
                     enddo
                     read(iread_rcoll,*,end=17928)idum
                  enddo
17928             if (idum .ne. isimnum) then
                     write(ierr,17930) tprp_num, isimnum,
     &                    trim(filename)
                     if (iout .ne. 0) write(iout, 17930) tprp_num,
     &                    isimnum, trim(filename)
                     if (iptty .ne. 0) write(iptty, 17930) tprp_num,
     &                    isimnum, trim(filename)
                     stop
                  end if
17930             format ('STOP - Error in inptrk for tprpflag = ', i2,
     &                 ' realization_num ', i4, ' not found in ', a)
                  jjj=0
                  read(iread_rcoll,'(a80)',end=17931) dummy_line
                  nulldum=null1(dummy_line)
                  do while(.not.nulldum)
                     jjj=jjj+1
                     read(iread_rcoll,'(a80)',end=17931) dummy_line
                     nulldum=null1(dummy_line)
                  enddo
17931             if (jjj .lt. 2) then
                     write(ierr,17932) tprp_num, isimnum,
     &                    trim(filename)
                     if (iout .ne. 0) write(iout, 17932) tprp_num,
     &                    isimnum, trim(filename)
                     if (iptty .ne. 0) write(iptty, 17932) tprp_num,
     &                    isimnum, trim(filename)
                     stop
                  end if
17932             format ('STOP - Error in inptrk for tprpflag = ', i2,
     &                 ' realization_num ', i4, ' data not found in ',
     &                 a)

                  maxprobdivs = jjj
                  rewind(iread_rcoll)
               else
                  if(iread_rcoll .eq. 0) read(locunitnum,'(a80)') wdd1
               endif
            endif
            read(locunitnum,'(a80)') wdd1
cHari setup pointer arrays for reversible and irreversible colloids
            if(wdd1(1:4).eq.'irre')then
               total_irrev = 1
            elseif(wdd1(1:4).eq.'reve')then
               total_rev = 1
            else
               write(ierr,*)'must specify irre or reve in ptrk'
               if(iout.ne.0)then
                  write(iout,*)'must specify irre or reve in ptrk'
               endif
               if(iptty.ne.0)then
                  write(iptty,*)'must specify irre or reve in ptrk'
               endif
            endif
         else
            backspace(locunitnum)
         endif


c     We are now at the point where we can count the
c     number of maxlayers

         maxlayers = 0
 991     continue
         read(locunitnum,'(a80)') dummy_line
         if(null1(dummy_line)) then
            goto 99
         else
            maxlayers = maxlayers + 1
         end if
         goto 991
 99      continue
         maxlayers = max(1,maxlayers)
         call done_macro(locunitnum)
c***********sptr***groemer 5/1/98***************************************
      else if (macro .eq. 'sptr') then
         nspeci = 1
         sptrak=.true.
         nzbtc = 0
         nsize_tprp = 0
         nplum = 0
         itensor = -999
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum,'(a80)') dummy_line

c     Section to read past keywords

         done=.false.
         do while(.not.done)
            done=.true.
            read(locunitnum,'(a80)') dummy_line
            
            read_flags:select case(dummy_line(1:2))
            case('tc', 'TC')
! If tcurve data is to be input read two more lines, zvd 18-Nov-02
               done = .false.
               read(locunitnum,'(a80)') dummy_line
               read(locunitnum,'(a80)') dummy_line
! If omr grid is being used
            case ('om', 'OM')
               done = .false.
               omr_flag = .true.
               read(locunitnum,'(i8)') omr_nodes
               read (locunitnum,'(a80)') dummy_line
               if (dummy_line(1:4).eq.'file') then
                  save_omr = .TRUE.
                  read (locunitnum,'(a100)') nmfil(23)
                  read (locunitnum,'(a11)') cform(23)
               else
                  save_omr = .FALSE.
                  backspace (locunitnum)
               end if
            case('po', 'PO')
               done = .false.
! If POD basis funtion derivatives are needed, zvd 16-Aug-04
               if (dummy_line(3:3) .eq. 'd' .or. 
     &              dummy_line(3:3) .eq. 'D') then
                  pod_flag = .true.
! Else po represents porosity
               end if
            case ('tp', 'TP')
               done = .false.
               if (dummy_line(3:3) .eq. 'o' .or. 
     &              dummy_line(3:3) .eq. 'O') then
! We are reading transport porosities 'tpor'
                  tpor_flag = .TRUE.
                  read(locunitnum,'(a80)') dumstring
                  if (dumstring(1:4) .ne. 'file') then
                     do
                        if (null1(dumstring)) exit
                        read(locunitnum,'(a80)') dumstring
                     enddo
                  else
                     read(locunitnum,'(a80)') dumstring
                  end if
               else
! Else we are reading the transport property 'tprp' models
                  if (itensor .eq. 0) itensor = 5
                  do
                     read(locunitnum,'(a80)') dummy_line
                     if(null1(dummy_line)) exit
                     nsize_tprp = nsize_tprp + 1
                     read(dummy_line, *) tprp_num
                     if(tprp_num.ge.11.and.tprp_num.le.14) then
                        div_flag=.true.
                        read(locunitnum,'(a80)') dummy_line
                        if (dummy_line(1:4) .eq. 'file') then
                           backspace(locunitnum)
                           backspace(locunitnum)
                           read(locunitnum,*) tprp_num,realization_num
                           read(locunitnum,'(a80)') dummy_line
                        else
c     We are reading eqn parameters from tprp line
                           backspace(locunitnum)
                        end if
                        if (dummy_line(1:4) .eq. 'file') then
                           do icount = 1, 80
                              filename(icount:icount) = ' '
                           enddo
                           read(locunitnum,'(a80)') filename  
                           iread_rcoll = open_file(filename,'old')
                           maxprobdivs = 1
                        end if
                        if(tprp_num.lt.13) then
                           if(iread_rcoll .eq. 0) then
                              write (ierr, 19929)
                              if (iout .ne. 0) write (iptty, 19929)
                              if (iptty .ne. 0) write (iptty, 19929)
                              stop
                           end if
19929                      format('STOP - CDF table data must be ',
     &                          'input using an external data file -- ',
     &                          'check sptr macro.')
                     
                           read(iread_rcoll,'(a80)', end=19928) 
     &                          dummy_line
                           read(iread_rcoll,*,end=19928)idum
                           do while (idum.ne.realization_num)
                              nulldum=.false.
                              do while(.not.nulldum)
                                 read(iread_rcoll,'(a80)',end=19928)
     &                                dummy_line
                                 nulldum=null1(dummy_line)
                              enddo
                              read(iread_rcoll,*,end=19928)idum
                           enddo
19928                      if (idum .ne. realization_num) then
                              write(ierr,19930) tprp_num, 
     &                             realization_num, trim(filename)
                              if (iout .ne. 0) write(iout, 19930) 
     &                             tprp_num, realization_num, 
     &                             trim(filename)
                              if (iptty .ne. 0) write(iptty, 19930) 
     &                             tprp_num, realization_num, 
     &                             trim(filename)
                              stop                
                           end if
19930                      format ('STOP - Error in insptr for ',
     &                          'tprpflag = ', i2, ' realization_num ', 
     &                          i4, ' not found in ', a)
                           jjj=0
                           read(iread_rcoll,'(a80)', end=19931) 
     &                          dummy_line
                           nulldum=null1(dummy_line)
                           do while(.not.nulldum)
                              jjj=jjj+1
                              read(iread_rcoll,'(a80)',end=19931) 
     &                             dummy_line
                              nulldum=null1(dummy_line)
                           enddo
19931                      if (jjj .lt. 2) then
                              write(ierr,19932) tprp_num, 
     &                             realization_num, trim(filename)
                              if (iout .ne. 0) write(iout, 19932) 
     &                             tprp_num, realization_num, 
     &                             trim(filename)
                              if (iptty .ne. 0) write(iptty, 19932)  
     &                             tprp_num, realization_num, 
     &                             trim(filename)
                              stop
                           end if
19932                      format ('STOP - Error in insptr for ', 
     &                          'tprpflag = ', i2, ' realization_num ', 
     &                          i4, ' data not found in ', a)

                           maxprobdivs = jjj
                           rewind(iread_rcoll)
c                           exit
                        endif
                        exit
                     end if
                  end do

c Read past the itrc group (model assignment)
                  do
                     read(locunitnum,'(a80)') dummy_line
                     if(null1(dummy_line)) exit
                  end do
               end if
            case('sa', 'SA')
               done = .false.

            case('pe', 'PE')
               done = .false.
            case('de', 'DE')
               done = .false.
            case('pr', 'PR')
               done = .false.
            case('te', 'TE')
               done = .false.
            case('zo', 'ZO')
               done = .false.
            case('id', 'ID')
               done = .false.
            case('wt', 'WT')
               done = .false.
            case('vo','VO')
               done = .false.
            case('xy','XY')
               done = .false.
            case ('ip', 'IP') 
               done = .false.
            case ('tr', 'TR') 
               done = .false.
            case ('in', 'IN')
               done = .false.
            case ('ex', 'EX')
               done = .false.
            case ('ou', 'OU')
               done = .false.
            case('zb','ZB')
               done = .false.
               read(locunitnum,*) ireaddum
               if(ireaddum.gt.0) then
                  read(locunitnum,'(a80)') dummy_line
                  nzbtc = ireaddum
               else
                  nzbtc = 0
               end if
            case('pl','PL')
c march 14, 02 s kelkar, from /scratch/ymp/kelkar/fehmv2.12.......
c     Read past plum input
               done = .false.
               read(locunitnum,*) ireaddum
               if(ireaddum.gt.0) then
                  read(locunitnum,'(a80)') dummy_line
                  nplum = ireaddum
               else
                  nplum = 0
               end if
            case('cl','CL')
c Jan 5, 06 S kelkar read past keyword cliff
               done = .false.
               if (dummy_line(1:6).eq.'cliffg') then
                  read(locunitnum,'(a80)') dummy_line
                  read(locunitnum,'(a80)') dummy_line
                  read(locunitnum,'(a80)') dummy_line
               else if (dummy_line(1:6).eq.'cliffr') then
                  read(locunitnum,'(a80)') dummy_line
               end if
            case('co','CO')
c Jan 27, 06 S kelkar read past keyword corner
               done = .false.
            case('ca','CA')
c Jan 13, 05 S kelkar read past capture node list
               done = .false.
               do
                  read(locunitnum,'(a80)') wdd1
                  if (wdd1(1:4) .eq. 'file') then
                     read (locunitnum,'(a80)') wdd1
                     exit
                  else if (null1(wdd1)) then
                     exit
                  end if
               end do
            case('sp','SP')
c Jan 27, 05 S kelkar read past "spring" node list
               done = .false.
               do
                  read(locunitnum,'(a80)') wdd1
                  if (wdd1(1:4) .eq. 'file') then
                     read (locunitnum,'(a80)') wdd1
                     exit
                  else if (null1(wdd1)) then
                     exit
                  end if
               end do
            case default
! The first non-character line read should contain courant factor, etc.
               if (itensor .eq. -999) then
                  done = .false.
                  call parse_string(dummy_line, imsg, msg, xmsg, cmsg,
     &                 nwds)
                  if(nwds.ge.3) then
                     itensor = imsg(3)
                  else
                     itensor = 0
                  end if
               else if (itensor .eq. 0) then
! We are done with keywords and haven't ready any transport properties
! This is really old style input, so we needed to read past the
! transport property line
                  nsize_tprp = max(1,nsize_tprp)
               else
                  backspace (locunitnum)
               end if
                  
            end select read_flags
         end do

c...................................................................
 
         read (locunitnum, '(a80)') dumstring
         call parse_string(dumstring, imsg, msg, xmsg, cmsg, nwds)
         if (msg(2).eq.1) then
            ist=imsg(2)
         else
            ist=xmsg(2)
         end if
         if (nwds .gt. 2) then
            do i = 3, nwds
               if (msg(i) .eq. 3) then
                  if (cmsg(i)(1:4) .eq. 'save') sptr_flag = .true.
               end if
            end do
         end if   
         read (locunitnum, *) nx,ny,nz
         read (locunitnum,'(a80)') dummy_line
         read (locunitnum,'(a80)') dummy_line
         read (locunitnum,'(a80)') dummy_line         

         num_part = 0
         if(ist.eq.2) then
            if(.not. null1(dummy_line)) then
! There shouldn't be any more data for sptr macro
               write (ierr, 110) ist
               if (iptty .ne. 0) write(iptty, 110) ist
              stop
            else
               if(icnl.ne.0) nz = 1
               num_part=nx*ny*nz
            end if
         else if (ist .eq. 0 .or. abs(ist) .eq. 1 .or. ist .eq. 3) then
            if(null1(dummy_line)) then
! There should be more data for sptr macro
               write (ierr, 120) ist
               if (iptty .ne. 0) write(iptty, 120) ist
               stop
            else
               if (dummy_line(1:4) .eq. 'file') then
                  read(locunitnum,'(a100)') filename
                  inptread = open_file(filename,'old')
                  opened = .true.
                  read (inptread,'(a80)') dummy_line
                  if (dummy_line(1:4) .eq. 'FEHM') then
! This file was generated by FEHM and the title, last seed and 
! ending time were written to the file
                     read(inptread,'(a80)') dummy_line
                     read(inptread,'(a80)') dummy_line
                  else
                     backspace inptread
                  end if
               else
                  backspace locunitnum
                  inptread = locunitnum
                  opened = .false.
               end if
               do
                  read (inptread,'(a80)') dummy_line
                  if (dummy_line(1:4) .eq. 'zbtc' .or.
     &                 null1(dummy_line)) exit
                  num_part = num_part + 1
               end do
               if (opened) then
                  close (inptread)
                  opened = .false.
               end if
            end if
         else
            write(ierr,130) ist
            if (iptty .ne. 0) write(iptty, 130) ist
            stop            
         end if
 110     format ('Error at end of sptr macro, input not terminated',
     &        ' with blank line, ist = ', i2)
 120     format ('Error at end of sptr macro, no particle position',
     &        ' input for ist = ', i2)
 130     format ('Error in sptr macro, illegal value for ist:', i2)

         call done_macro(locunitnum)
c***********sptr***groemer 5/1/98***************************************

      else if (macro .eq. 'ppor') then
c need porosity model
         call start_macro(inpt, locunitnum, macro)
         read (locunitnum,*) iporos
         call done_macro(locunitnum)
c gaz 041620 (need two phases to work with)
      else if(macro .eq. 'rich') then
         jswitch = 1
         if(irdof.eq.13) irdof = 0
      else if (macro .eq. 'vcon') then
! need number of variable conductivity models
         call start_macro(inpt, locunitnum, macro)

         do
            read(locunitnum,'(a80)') dumstring
            if (null1(dumstring)) exit
            numvcon = numvcon + 1
         enddo

         call done_macro(locunitnum)
      else if (macro .eq. 'vbou') then
         ivboun = 1
      else if (macro .eq. 'zeol') then
         izeolites = 1
      else if (macro .eq. 'strs') then
         call start_macro(inpt, locunitnum, macro)
         read(locunitnum,*) istrs, ihms
         i = 0
         do
            read(locunitnum,'(a80)') dumstring
            if(dumstring(1:9).eq.'stressend' .or. dumstring(1:8) .eq.
     &         'end strs' .or. dumstring(1:7) .eq. 'endstrs') exit
            if(dumstring(1:9).eq.'permmodel') then
               do 
                  read(locunitnum,'(a80)') dumstring
                  if (null1(dumstring)) exit
                  read(dumstring,*) idumm
	          backspace locunitnum
	          i = i + 1
	          if(idumm.eq.1) then
 	           read(locunitnum,*) idumm
c	          else if (idumm .eq. 2 .or. idumm .eq. 4) then
	          else if (idumm .eq. 2) then
                     read(locunitnum,*) idumm, (adumm, ja = 1, 9) 
	          else if (idumm .eq. 3) then
                     read(locunitnum,*) idumm, (adumm, ja = 1, 9) 
	          else if (idumm .eq. 4) then
                     read(locunitnum,*) idumm, (adumm, ja = 1, 14)   
	          else if(idumm.eq.5) then
                     read(locunitnum,*) idumm,(adumm,ja=1,6)
	          else if(idumm.eq.6) then
                     read(locunitnum,*) idumm,(adumm,ja=1,11)
	          else if(idumm.eq.7) then
                     read(locunitnum,*) idumm,(adumm,ja=1,6)
	          else if(idumm.eq.8) then
                     read(locunitnum,*) idumm,(adumm,ja=1,9)         
                  else if(idumm.eq.11) then
                     read(locunitnum,*) idumm,(adumm,ja=1,3)
                  else if(idumm.eq.21) then
                     read(locunitnum,*) idumm,(adumm,ja=1,13)
                  else if(idumm.eq.22) then
                     read(locunitnum,*) idumm,(adumm,ja=1,10)
                  else if(idumm.eq.222) then
                     read(locunitnum,*) idumm,(adumm,ja=1,9)
                     read(locunitnum,*) (adumm,ja=1,10)
                  else if(idumm.eq.31) then
                     read(locunitnum,*) idumm,(adumm,ja=1,2)      
                  else if(idumm.eq.91) then
c                     open(unit=97,file='debug_permmodel_91.dat')
                     read(locunitnum,*) 
                     read(locunitnum,*)
                  else
                     read(locunitnum,*) idumm             
                  endif
               enddo
            endif
         enddo
         if(i.ge.1) then
          allocate(ispmt(i))
          allocate(spm1f(i),spm2f(i),spm3f(i),spm4f(i),spm5f(i))
          allocate(spm6f(i),spm7f(i),spm8f(i),spm9f(i),spm10f(i))
          allocate(spm11f(i),spm12f(i),spm13f(i),spm14f(i))
          allocate(spm1f222(i),spm2f222(i),spm3f222(i),spm4f222(i))
          allocate(spm5f222(i),spm6f222(i),spm7f222(i),spm8f222(i))
          allocate(spm9f222(i),spm10f222(i))
          ispmt = 0
          spm1f = 0.
          spm2f = 0.
          spm3f = 0.
          spm4f = 0.
          spm5f = 0.
          spm6f = 0.
          spm7f = 0.
          spm8f = 0.
          spm9f = 0.
          spm10f = 0.
          spm11f = 0.
          spm12f = 0.
          spm13f = 0.
          spm14f = 0.
          spm1f222 = 0.
          spm2f222 = 0.
          spm3f222 = 0.
          spm4f222 = 0.
          spm5f222 = 0.
          spm6f222 = 0.
          spm7f222 = 0.
          spm8f222 = 0.
          spm9f222 = 0.
          spm10f222 = 0.
         endif
         call done_macro(locunitnum)
         
      end if
      
      go to 10
      
 40   rewind locunitnum
      if(ngroups.eq.0)then
         ngroups=ncpnt
         if(irun.eq.1) then
            allocate(fzero(max(ngroups,1)))
            allocate(group(max(ngroups,1),ncpnt))
            allocate(pos(max(ngroups,1),ncpnt))
            allocate(n_couple_species(max(ngroups,1)))
         end if
         mdof_sol = 1
         n_couple_species = 0
         do igrp = 1, ngroups
            pos(igrp,1)= igrp
            n_couple_species(igrp)=1
         enddo
      endif
      if(idpdp.ne.0)mdof_sol=mdof_sol*2
      if(allocated(group_mat)) deallocate(group_mat)

      if ((iccen.eq.1).and.ptrak) then
         write(ierr,*) 'ERROR: Can not have both particle tracking' 
         write(ierr,*) '(ptrk) and tracer input (trac).'
         write(ierr,*) 'Code Aborted in scanin'
         if (iatty.ne.0) then
            write(iatty,*) 'ERROR: Can not have both particle tracking'
            write(iatty,*) '(ptrk) and tracer input (trac).'
            write(iatty,*) 'Code Aborted in scanin'
         endif
         stop
      endif

c gaz added 112014    
      
       if (iccen.eq.1 .and. rxn_flag .eq. 0) then
c     Check that no solid species are used without reaction
         if (nimm .ne. 0) then
            write(ierr,*) 'ERROR: Can not using solid species in trac',
     &           ' without rxn'
            write(ierr,*) 'Code Aborted in scanin'
            if (iatty.ne.0) then
               write(iatty,*) 'ERROR: Can not using solid species in ',
     &           'trac without rxn'
               write(iatty,*) 'Code Aborted in scanin'
            end if
            stop
        end if
      end if 
      
      return

 50   write (ierr, 55) 'STOP'
      if (iptty .ne. 0) 
     &     write (iptty, 55) 'STOP'
      stop
 55   format ("Missing '", a, "' statement in input. STOPPING")
 56   format ("Check macro: ", a)
      end subroutine scanin

      subroutine zone(cnum, infile)
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
CD1 Create FEHM zones using geometry or node lists
CD1 Enter properties using geometric description.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-JAN-94    Z. Dash        22      Add prolog/major cleanup.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/zone.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:25:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Wed Jun 12 15:21:26 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.8   Mon Jun 03 11:18:46 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.7   Fri May 31 10:55:24 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.6   Fri Feb 16 13:59:46 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.5   Fri Feb 02 14:34:12 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   Fri Jan 12 18:04:06 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.3   08/22/95 13:47:30   llt
CD2 quotes around internal read name doesn't work on IBM.
CD2 
CD2    Rev 1.2   08/18/95 10:32:26   llt
CD2 needed quotes around ltest for cray to read
CD2 
CD2    Rev 1.1   03/18/94 16:04:00   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:42   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   cnum            INT      I    Number of times zone has been called
CD3   infile          INT      I    Unit number of file to be read.
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
CD3   inzone                   I    Zone input file.
CD3   iout                     O    General output file.
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
CD4   icnl            INT      faai   Problem dimension
CD4   idpdp           INT      faai   Parameter which indicates if the double
CD4                                     porosity / double permeability
CD4                                     solution is enabled
CD4   idualp          INT      faai   Parameter which indicates if the dual
CD4                                     porosity solution is enabled
CD4   inpt            INT      faai   Unit number for input file
CD4   inzone          INT      faai   Unit number for zone input file
CD4   izonef          INT      fbb    Zone in which each node is located
CD4   lenintg         INT      param  Converts bits to words for allocating
CD4                                     memory
CD4   macroread(18)   LOGICAL  macro  Flag denoting if macro zone has been read
CD4   wdd1            CHAR     faac   Alternate character input string
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   mmgetblk                 Allocate memory to an array
CD4   mmrelblk                 Deallocate array memory
CD4   null1            LOGICAL  Check for null lines or 0's in lines
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
CD5   i               INT      Loop index
CD5   icode           INT      Return value from mmgetblk, mmrelblk routines
CD5   ipncord         POINTER  Pointer to variable array ncord
CD5   izone           INT      Identifier for zone being defined
CD5   izonel          INT      Identifier for last zone defined
CD5   ltest           CHAR     Variable for reading character input     
CD5   macro           CHAR     Current macro being read
CD5   ncord           REAL*8   Nodes found within a zone
CD5   nin             INT      Number of nodes assigned to zone
CD5   nodez           INT      Number of node nearest given coordinates
CD5   nsl             INT      Number of coordinates in element
CD5   xg              REAL*8   X coordinate defining node
CD5   xz              REAL*8   X coordinates defining zone
CD5   yg              REAL*8   Y coordinate defining node
CD5   yz              REAL*8   Y coordinates defining zone
CD5   zg              REAL*8   Z coordinate defining node
CD5   zz              REAL*8   Z coordinates defining zone
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
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
CD8 Zones are a convenient way to group nodes and specify properties,
CD8 but they are not an essential function of the code.
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
CPS BEGIN  zone
CPS 
CPS   backup a line in file to be read
CPS      
CPS   FOR each node
CPS       set zone to zero
CPS   END FOR
CPS      
CPS   END IF
CPS   
CPS   IF this is a 3D problem
CPS      set nsl to 8
CPS   ELSE
CPS      set nsl to 4
CPS   END IF
CPS   
CPS   REPEAT
CPS     read input line
CPS     IF macro read is zone
CPS        REPEAT
CPS          read input line
CPS        EXIT IF null line is read
CPS          reread the input line to get zone identification number
CPS          IF the zone id is > 0
CPS             read input line
CPS             IF macros list and nnum are not read
CPS                IF this zone isn't the same as the last zone
CPS                   set node count to zero
CPS                END IF
CPS                read zone X coordinates
CPS                read zone Y coordinates
CPS                IF this is a 3D problem
CPS                   read zone Z coordinates
CPS                END IF
CPS                call setzone to find nodes in zone
CPS             ELSE IF macro read is list
CPS                set node count to zero
CPS                REPEAT
CPS                  read input line
CPS                  IF this is a 3D problem
CPS                     read node X, Y, and Z coordinates from input line
CPS                  ELSE IF this is a 2D problem
CPS                     read node X and Y coordinates from input line, . . .
CPS                     . . . set Z coordinate to 0
CPS                  END IF
CPS                  call near3 to determine number of node closest to . . .
CPS                  . . . coordinates read
CPS                  set zone id number for identified node
CPS                UNTIL null line is read
CPS             ELSE IF macro read is nnum
CPS                read number of nodes belonging to zone and the node numbers
CPS                FOR each node read
CPS                    set zone id number for the node
CPS                END FOR
CPS             END IF
CPS             write zone number and each node in zone to data check file
CPS          END IF
CPS        UNTIL null line is read
CPS     END IF
CPS     
CPS     set node count to zero
CPS     FOR each node in the problem
CPS         IF the node is not in a defined zone
CPS            save its number and increment count of nodes not in a zone
CPS         END IF
CPS     END FOR
CPS     
CPS     IF there were nodes that were not in a zone
CPS        write warning and node numbers of unassigned nodes to the . . .
CPS        . . . data check file
CPS     END IF
CPS     
CPS     IF this is a dual porosity problem
CPS        set number of dual porosity nodes
CPS        FOR each node
CPS            set zone for first and second matrix level dual porosity node
CPS        END FOR
CPS     ELSE IF this a dual porosity/dual permeability problem
CPS        set number of dual porosity/dual permeability nodes
CPS        FOR each node
CPS            set zone for first matrix level dual porosity/dual . . .
CPS            . . . permeability node
CPS        END FOR
CPS     ELSE
CPS        set number of nodes
CPS     END IF
CPS     
CPS   UNTIL zone macro has been read
CPS   
CPS END  zone
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      use trxnvars
      implicit none

      logical null1, null_new, cdum
      integer cnum, i, infile, izone, izonel, nin, nodez, nsl
      integer nxy, icnl_old, nin_old, i_old, ja, jb, jc, izonn
      integer izunit, open_file
      character* 4 macro, cmacro
      character(20), allocatable :: znametmp(:)
      character*80 ltest
      real*8 xg, xz(8), yg, yz(8), zg, zz(8)
      real*8 tol_zone, zxy_min, zxy_max 
      integer imodel, j, n_n_n, zmaxtmp
c gaz 022818
      integer neq_t
      integer zone_dpadd, i3d_2d, i3d_rad, num_zones, lsize
      integer, allocatable :: znumtmp(:)
c ncord now in combi
c      integer, allocatable :: ncord(:)
      integer, allocatable :: izonef_old(:)
      integer, allocatable :: zone_list(:), tmp_list(:)
      character*20 zonetmp
      character*30 zonesavename
      integer k, curzone
      logical ex, zone_check

      save zmaxtmp
      allocate(ncord(n0))

c     Dual perm or dual porosity value to add to get zones
      zonesavename(1:4) = 'zone'
      izone_save = 0
      zone_dpadd = 100
      lsize = 100
      allocate (zone_list(lsize))
      num_zones = 0

      izonel = 0
      nin = 0
      backspace infile

      if (icnl .eq. 0)  then
         nsl = 8
      else
         nsl = 4
      end if

 50   continue
      i3d_2d = 0
      i3d_rad = 0
      read (infile, '(a80)')  wdd1
      read (wdd1,*) macro
      do i = 5,77
         if(wdd1(i:i+3).eq.'conv') i3d_2d = 1
         if(wdd1(i:i+2).eq.'rad') i3d_rad = 1
         if(wdd1(i:i+3).eq.'save') izone_save = 1
      enddo
      if(i3d_2d.eq.1.and.icnl.eq.0) then
         write(ierr,*) 'i3d_2d parameter ignored for 3d problem'
         if(iout.ne.0) 
     &        write(iout,*) 'i3d_2d parameter ignored for 3d problem'
         if(iptty.ne.0) 
     &        write(iptty,*) 'i3d_2d parameter ignored for 3d problem'
      endif
      if (.not. allocated(zonenames)) then
         zmaxtmp = 100
         allocate (zonenames(zmaxtmp), zonenums(zmaxtmp))
      end if
      if (macro .eq. 'zone'.or.macro .eq. 'zonn')  then
         cmacro = macro
         if(macro .ne. 'zonn') then
            izonef = 0
            izonn = 0
	    numzones = 0
	    zonenames = '*'
	    zonenums = 0
	    zonemax = 0            
         else
            allocate(izonef_old(n0))
            izonef_old = izonef
            izonn = 1
         endif
      endif
 60   continue
      read (infile, '(a80)') wdd1
      if (null1(wdd1)) go to 90
      backspace infile
      !read(infile, *) izone
      read(infile, *) zonetmp
      do k = 1, numzones
         if (zonenames(k) .eq. zonetmp) then
            curzone = k
            goto 63
         endif
      enddo
      curzone = numzones + 1
      if (zonetmp(1:1) .ne. '-') numzones = numzones + 1
 63   if (curzone .gt. zmaxtmp) then
         allocate (znametmp(zmaxtmp),znumtmp(zmaxtmp))
         znametmp = zonenames
         znumtmp = zonenums 
         deallocate (zonenames, zonenums)
         allocate (zonenames(zmaxtmp*10), zonenums(zmaxtmp*10))
         zonenames(1:zmaxtmp) = znametmp
         zonenums(1:zmaxtmp) = znumtmp
         deallocate (znametmp, znumtmp)
         zmaxtmp = zmaxtmp*10
      end if
      zonenames(curzone) = zonetmp
      read(zonenames(curzone), *, err=61) zonenums(curzone)
      goto 62
 61   zonenums(curzone) = zonemax + 1
 62   zonemax = max(zonemax + 1, zonenums(curzone))
      izone = zonenums(curzone)

c zvd 01/04/2012 Keep track of zones that are defined so auto generated double permeability or porosity nodes can be reported
      num_zones = num_zones + 1
      if (num_zones .gt. lsize) then
         allocate (tmp_list(lsize))
         tmp_list = zone_list
         deallocate (zone_list)
         allocate (zone_list(lsize*2))
         zone_list(1:lsize) = tmp_list
         lsize = lsize*2
         deallocate (tmp_list)
      end if
      zone_list(num_zones) = izone
             
c     Determine if zone_dpadd needs to be increased to 1000
      if(izone.gt.99) zone_dpadd = 1000

      if (izone .gt. 0) then
c     first check if saved already  
c     check if list or nnum occurs          
        call check_save_zone(1,izone,zone_check)     
         read(infile, '(a4)') macro
         if (macro .ne. 'list' .and. macro .ne. 'nnum' .and.
     &        macro .ne. 'xyli'.and. macro.ne. 'all ' .and.
     &        macro .ne. 'jajb') then
            backspace  infile
            if (izone .ne. izonel) nin = 0
            if(i3d_rad.eq.1.and.icnl.eq.0) then
c     radial input for 3D problems (uses 2D input)
               read (infile, *) (xz(i), i = 1, 4)
               read (infile, *) (yz(i), i = 1, 4)  
c               call setzone(izone, nin, ncord, 4, xz, yz, zz, 1)
               call setzone(izone, nin, 4, xz, yz, zz, 1)
            else if(i3d_2d.eq.0.or.icnl.eq.0) then
               read (infile, *) (xz(i), i = 1, nsl)
               read (infile, *) (yz(i), i = 1, nsl)
               if (icnl .eq. 0) then
c**** 3-d calculation ****
                  read (infile, *)  (zz(i), i = 1, nsl)
               end if
c               call setzone(izone, nin, ncord, nsl, xz, yz, zz, 0)
               call setzone(izone, nin, nsl, xz, yz, zz, 0)
            else
c     3-d zones in 2D model (for consistency when extracting slices in 3d)               
               read (infile, *) (xz(i), i = 1, 8)
               xz(3) = xz(2)
               xz(4) = xz(1)                 
               xz(1) = xz(5)
               xz(2) = xz(6)
               read (infile, *) (yz(i), i = 1, 8)
c     note we overwrite yz with zz                 
               read (infile, *) (zz(i), i = 1, 8)
               yz(1) = zz(5)
               yz(2) = zz(6)
               yz(3) = zz(2)
               yz(4) = zz(1)
c               call setzone(izone, nin, ncord, nsl, xz, yz, zz, 0) 
               call setzone(izone, nin,  nsl, xz, yz, zz, 0)  
            endif
         else if(macro .eq. 'xyli') then
c     read in nodes in zone from xy list
            nxy = 0
            nin = 0
            i = 0
            read(infile,*) tol_zone, zxy_min, zxy_max
 71         read(infile, '(a80)') ltest
            if(.not.null1(ltest)) then
               read(ltest, *, end = 81, err = 81) xg, yg
               i_old = i
               icnl_old=icnl
               icnl=1
               call near3(xg,yg,0.0,i,0)
               icnl=icnl_old
               if(i_old.eq.i) go to 71
               xg=cord(i,1)
               yg=cord(i,2)
               nxy = nxy + 1
            else
               goto 81
            end if
            xz(1)=xg-tol_zone
            xz(2)=xg+tol_zone
            xz(3)=xg+tol_zone
            xz(4)=xg-tol_zone
            xz(5)=xg-tol_zone
            xz(6)=xg+tol_zone
            xz(7)=xg+tol_zone
            xz(8)=xg-tol_zone
            yz(1)=yg-tol_zone
            yz(2)=yg-tol_zone
            yz(3)=yg+tol_zone
            yz(4)=yg+tol_zone
            yz(5)=yg-tol_zone
            yz(6)=yg-tol_zone
            yz(7)=yg+tol_zone
            yz(8)=yg+tol_zone
            zz(1)=zxy_max
            zz(2)=zxy_max
            zz(3)=zxy_max
            zz(4)=zxy_max
            zz(5)=zxy_min
            zz(6)=zxy_min
            zz(7)=zxy_min
            zz(8)=zxy_min
            nin_old = 0
c            call setzone(izone, nin_old, ncord(nin+1:neq),            
            call setzone(izone, nin_old, 
     &           nsl, xz, yz, zz, 0)
            nin=nin+nin_old
            go to 71
 81         continue
         else if(macro .eq. 'list') then
c     read in coordinates for nodes in zone
            nin = 0
 70         read(infile, '(a80)') ltest
            if(.not.null_new(ltest)) then
               if(icnl .eq. 0) then
                  read(ltest, *, end = 80, err = 80) xg, yg, zg
               else
                  if(i3d_2d.eq.1) then
                     read(ltest, *, end = 80, err = 80) xg, yg, zg
                     yg = zg
                     zg = 0.0
                  else
                     read(ltest, *, end = 80, err = 80) xg, yg
                     zg = 0.0
                  endif
               end if
            else
               goto 80
            endif
            nin = nin + 1
            call near3(xg, yg, zg, nodez, 0)
            ncord(nin) = nodez
            izonef(nodez) = izone
            go to 70
 80         continue
         else if(macro .eq. 'nnum') then
c     read in nodes belonging to zone
            read(infile, *) nin, (ncord (i), i = 1, nin)
            do i = 1, nin
               if (ncord(i) .gt. n0) then
                  write (ierr, 6008) cmacro
                  write (ierr, 6009) ncord(i), n0
                  if (iout .ne. 0) write (iout, 6009) ncord(i), n0
                  if (iptty .gt. 0) write (iptty, 6009) ncord(i), n0
                  stop
               else
                  izonef(ncord(i)) = izone
               end if
            end do
         else if(macro .eq. 'all ') then     
c check for multiple porosity models
          if(gdpm_flag.ne.0.or.gdkm_flag.ne.0.or.
     &       idualp.ne.0.or.idpdp.ne.0) then 
              if(gdpm_flag.ne.0) neq_t = neq_primary
              if(gdkm_flag.ne.0) neq_t = neq_primary
              if(idualp.ne.0) neq_t = neq
              if(idpdp.ne.0) neq_t = neq
              do i = 1, neq_t
               izonef(i) = izone
              enddo
          else    
            do i = 1, n0
               izonef(i) = izone
            enddo 
          endif
         else if(macro .eq. 'jajb') then
            do
               read (infile, '(a80)') ltest
               if(null_new(ltest)) exit
               read (ltest, *) ja, jb, jc
               do i = ja, jb, jc
                  izonef(i) = izone
               end do
            end do
         endif
 6008    format (' **** Invalid input: macro ', a4, ' ****')
 6009    format(' **** Invalid node specified, ', i8, 
     .        ' is greater than ', 'n0 (', i8, ' ): stopping ****')
c**** print out nodes in izone ****
         nin = 0
c     Change to n0 (used to be neq) - BAR 12-15-99
         do i = 1, n0
            if(izonef(i) .eq. izone) then
               nin = nin + 1
               ncord(nin) = i
            endif
         end do
         if (ischk .ne. 0) then
            write(ischk, 6010) nin, izone
 6010       format(/, 1x, i8,' nodes contained in zone = ',
     &           i10, /)
            write(ischk, 6011) (ncord(i), i = 1, nin)
 6011       format (10i8)
         end if
c if zonesave ne 0, then save it to a file
c first check if saved already
c         call check_save_zone(1,izone,zone_check)
         if(izone_save.ne.0) then            
c create a file to save the zone information          
          write(zonesavename(5:9),'(i5)') izone+10000
          zonesavename(5:5) = '0'
          zonesavename(10:14) = '.save'
          ex = .false.
c          inquire (file = zonesavename, exist = ex)
          if(ex) then
c abort zone save because zonefile exists
          if (ischk .ne. 0) then
           write(ischk, 6014)  izone, zonesavename(1:14)
6014      format(1x, 'zone 'i5,' file creation aborted, ',a14,' exists')
          endif
          if (iout .ne. 0) then
           write(iout, 6014)  izone, zonesavename(1:14)   
          endif     
          if (iptty .ne. 0) then 
           write(iptty, 6014)  izone, zonesavename(1:14)   
          endif  
          go to 6015
          endif
          izunit=open_file(zonesavename,'unknown')
c write zone in nnum style
          if (ischk .ne. 0) then
           write(ischk, 6013)  izone, zonesavename(1:14)
6013      format(1x, 'zone 'i5,' saved in file  ', a14)          
          endif
          if (iout .ne. 0) then
           write(iout, 6013)  izone, zonesavename(1:14)   
          endif     
          if (iptty .ne. 0) then
           write(iptty, 6013)  izone, zonesavename(1:14)   
          endif 
          write(izunit,*) izone
          write(izunit,*) 'nnum'
c   write out nodes belonging to zone
          nin = 0
          do i = 1, n0
            if(izonef(i) .eq. izone) then
               nin = nin + 1
               ncord(nin) = i
            endif
         end do
          write(izunit,*) nin, (ncord(i), i = 1,nin)
c add elements to saved zone        
          call zone_elem(1,izunit,zonesavename,izone,nin)
          call zone_elem(2,izunit,zonesavename,izone,nin)
          close(izunit)
6015      continue  
       endif
c end of saved zones                  
         izonel = izone
         go to 60
      endif
c     check which nodes don't belong to a zone
c gaz 022818
90    nin = 0
c     Changed to neq_primary (used to be neq) BAR - 12-15-99
      do i = 1, neq_primary
         if(izonef(i) .eq. 0) then
            nin = nin + 1
            ncord(nin) = i
         endif
      end do
      if (nin .ne. 0) then
         if (ischk .ne. 0) write(ischk, 6012) nin, cnum
c     write(ischk, 6011) (ncord(i), i = 1, nin)
      end if
 6012 format(/, 1x, i8,
     &     ' nodes not assigned to a zone in call # ', i10)


c     Assign zones for GDPM nodes

      if(gdpm_flag.ne.0.and.gdkm_flag.eq.0) then

c     Set zones for GDPM nodes for the case in which zone has
c     already been called
c     Convention: all GDPM nodes are assigned a zone that
c     is 100 + the zone number of the primary node
c     unless the zone numbers declared are greater than 100,
c     then we use 1000 (zone_dpadd is the variable)

         n_n_n = neq_primary
         do i = 1, neq_primary
            
c     Loop over all GDPM nodes for primary node i
c     ngdpm_layers(imodel) = 0 for imodel = 0 (i.e. no GDPM nodes)
            imodel = igdpm(i)
            do j = 1, ngdpm_layers(imodel)
               n_n_n = n_n_n + 1
               if(izonn.eq.1) then
                  if(izonef(i).ne.izonef_old(i)) then              
c     Only assign the zone number this way if
c     it hasn't already been assigned a non-zero value
c     for example, in a zone with the nnum option
                     izonef(n_n_n) = izonef(i) + zone_dpadd
                  endif
               else
c gaz new 061817 to make consistent with dpdp
                 if(izonef(i).eq.0) then  
                   izonef(n_n_n) = izonef(i) + zone_dpadd
                 endif
               endif       
            end do 
         end do

      end if



c     check for dual porosity or dpdp solution
      if(idualp .eq. 1) then
c     dual porosity solution
         n = neq+neq+neq
         if (ischk .ne. 0) then
            write(ischk, *) 'dual porosity solution'
            write(ischk, *) 'first matrix level zone = ',
     .           'fracture level zone + ',zone_dpadd
            write(ischk, *) 'second matrix level zone = ', 
     .           'fracture level zone + ',zone_dpadd*2
         end if

c     This loop changed to set zones to their value plus 100
c     only if the nodes have not been explicitly set in the
c     zone definition. This allows the user to set the matrix
c     nodes in the zone macro and not have the code default
c     to zone number plus 100 (or 200).
c     Loop changed to accomodate the new zone_dpadd variable

         do  i = 1, neq
            if(izonn.eq.1)then
               if(izonef(i).ne.izonef_old(i)) then
                  izonef(i + neq) = izonef(i) + zone_dpadd
               endif
            else
               izonef(i + neq) = izonef(i) + zone_dpadd
            end if
            if(izonn.eq.1)then
               if(izonef(i).ne.izonef_old(i)) then
                  izonef(i + neq + neq) = izonef(i) + 2*zone_dpadd
               endif
            else
               izonef(i + neq + neq) = izonef(i) + 2*zone_dpadd
            end if
         end do
      else if(idpdp .ne. 0. or. gdkm_flag .ne. 0) then
c     dpdp solution or gdkm solution
        if(idpdp.ne.0) then
         n = neq+neq
         neq_t = neq
        else
c     gdkm : neq is total nodes            
         n = neq
         neq_t = neq_primary
        endif
         if (ischk .ne. 0) then
            write(ischk, *) 'dual porosity/dual permeability ',
     &           'solution'
            write(ischk, *) 'first matrix level zone = ',
     .           'fracture level zone + ',zone_dpadd
         end if
c     This loop changed to set zones to their value plus 100
c     only if the nodes have not been explicitly set in the
c     zone definition. This allows the user to set the matrix
c     nodes in the zone macro and not have the code default
c     to zone number plus 100.
c     Loop changed to accomodate the new zone_dpadd variable
         k = neq_t
         do i = 1, neq_t
c first check for gdkm primary node
          if(gdkm_flag.ne.0) then
           if(igdpm(i).ne.0) then
            k = k+1
            if(izonn.eq.1)then
               if(izonef(i).ne.izonef_old(i)) then
                  izonef(k) = izonef(i) + zone_dpadd
               endif
            else
               if(izonef(k).eq.0.and.izonef(i).ne.0) then
                  izonef(k) = izonef(i) + zone_dpadd
               endif
            end if
           endif
          else
            if(izonn.eq.1)then
               if(izonef(i).ne.izonef_old(i)) then
                  izonef(i+neq_t) = izonef(i) + zone_dpadd
               endif
            else
               if(izonef(i+neq_t).eq.0.and.izonef(i).ne.0) then
                  izonef(i+neq_t) = izonef(i) + zone_dpadd
               endif
            end if              
          endif          
         end do
         if (ischk .ne. 0) then
            do i = 1, num_zones
               nin = 0
               k = neq_t
               do j = 1, neq_t 
                if(gdkm_flag.ne.0) then
                 if(igdpm(j).gt.0) then
                  k = k +1
                  if(izonef(k) .eq. zone_list(i) + zone_dpadd) then  
                     nin = nin + 1
                     ncord(nin) = k
                  endif
                 endif
                else
                 if(izonef(j+neq_t) .eq. zone_list(i) + zone_dpadd) then
                     nin = nin + 1
                     ncord(nin) = j+neq_t   
                  endif
                endif
               end do
               if (nin .ne. 0) then
                  write(ischk, 6010) nin, zone_list(i) + zone_dpadd
                  write(ischk, 6011) (ncord(j), j = 1, nin)
               end if
            end do   
         end if      
      else
         n = neq
      endif
      go to 100
c     end if
      go  to  50
 100  continue
      deallocate(zone_list)  
      macroread(18) = .TRUE.
      deallocate(ncord)
      if(allocated(izonef_old)) deallocate(izonef_old)
      return
      end
      subroutine zone_elem(iflg,izunit,zonesavename,izone,nin)
c add element information to saved zones
c
      use combi
      use comdti
      use comai
            
      implicit none
      integer iflg,izone,izunit,i_elem_cover
      integer i,ii,j,ie,nin,ic,ic1,k,ib,id
c      ncord now in combi
c      integer ncord(*)
      integer, allocatable :: nopdum(:,:)
      integer, allocatable :: noodum(:,:)
      integer, allocatable :: iplace(:) 
      integer, allocatable :: iplace1(:) 
      integer, allocatable :: ielem_used(:)  
      integer, allocatable :: nei_list(:)
      integer, allocatable :: ncord_new(:)
      character*30 zonesavename
      save nopdum,noodum,iplace,ielem_used,iplace1
      save nei_list,ncord_new
      parameter(i_elem_cover = 0)
      if(iflg.eq.1) then
c          
c find elements associated with each node 
c
c calculate # connections for each node
c
      if(.not.allocated(iplace)) then
       allocate(iplace(n0))
      else
c only need to do this section once          
c this idicates the arrays are ready to use for iflg = 2         
       return
      endif
      iplace = 0
      do ie=1,nei
         do j=1,ns
            i=nelm((ie-1)*ns+j)
            if (i .ne. 0) then
               iplace(i)=iplace(i)+1
            endif
         enddo
      enddo 
      nemx = 0
      do i = 1, n0
       nemx = max(nemx,iplace(i))
      enddo
      allocate(nopdum(n0,nemx))
      allocate(ielem_used(nei))
      iplace=0
      do ie=1,nei
         do j=1,ns
            i=nelm((ie-1)*ns+j)
            if (i .ne. 0) then
               iplace(i)=iplace(i)+1
               nopdum(i,iplace(i))=ie
            endif
         enddo
      enddo
      else if(iflg.eq.2) then
c just use element - node relationship
       if(.not.allocated(iplace1)) allocate (iplace1(n0))
       iplace1 = 0
       ielem_used = 0
       do ii = 1, nin
        i = ncord(ii)
        iplace1(i) = 1
        do j = 1, iplace(i)
          ie = nopdum(i,j)
          ielem_used(ie) = ielem_used(ie) + 1
        enddo
       enddo
c iplace is the element list connected to node i in the current zone       
c iplace1 is list of nodes in current, (izone) zone       
c look for elements that contain only nodes in current zone
       if(.not.allocated(nei_list)) allocate(nei_list(nei))
       nei_list = 0
       ic = 0
       ib = 0
       id = 0
       do ie = 1, nei
        if(ielem_used(ie).gt.0) then
         if(ielem_used(ie).eq.ns) id = id +1
          ib = ib + 1
          do j= 1, ns
            i=nelm((ie-1)*ns+j)  
            if(iplace1(i).eq.0.and.i_elem_cover.eq.0) go to 201
          enddo
c decide on elemnt coverage option (i_elem_cover)          
          ic = ic + 1
          nei_list(ic) = ie
          go to 202
201       continue  
          write(ierr,*) 'zone = ',izone,' elem no = ',ie,' node = ',i,
     &      ' mixed zones'
202       continue          
        endif
       enddo
          write(ierr,*) 'number of elements that contain at least one '
     &     ,'node in zone ', izone, ' num elements ', ib
          write(ierr,*) 'number of elements that contain only  '
     &     ,'nodes in zone ', izone, ' num elements ', id
c now add element information
       write(izunit,*) 'elem'
       write(izunit,*) ic, ns
       if(.not.allocated(ncord_new)) allocate (ncord_new(n0))
c check for expanded node 
       ncord_new(1:nin) = ncord(1:nin)
       ic1 = nin
       do j = 1, ic
        ie = nei_list(j)
        write(izunit,*) ie, (nelm((ie-1)*ns+k),k = 1,ns)
        do k = 1,ns
          i = nelm((ie-1)*ns+k) 
          if(iplace1(i).eq.0) then
            ic1 = ic1 + 1
            iplace1(i) = -1
            ncord_new(ic1) = i
          endif
         enddo
       enddo
c check consistency     
c       do j = 1, ic
c        ie = nei_list(j)
c       enddo
c correct for expanded node list
      if (ic1.gt.nin) then
       write(izunit,'(a13)') 'nnum expanded'
       write(izunit,*) ic1-nin
       write (izunit,*) (ncord_new(i),i = nin + 1, ic1)
       write(izunit,*) '    '
      endif
c close file
      if(.not.allocated(zonesavenames)) 
     &       allocate(zonesavenames(maxsvzone))
       izone_sv_cnt = izone_sv_cnt + 1
       zonesavenames(izone_sv_cnt) = zonesavename
       close(izunit)
      else if(iflg.eq.-1) then
c deallocate memory
c gz 103018 check allocate status (some arrays may not have been allocated)
       if(allocated(iplace)) deallocate(iplace)
       if(allocated(iplace1)) deallocate(iplace1)
       if(allocated(nopdum)) deallocate(nopdum)
       if(allocated(ielem_used)) deallocate(ielem_used)
       if(allocated(nei_list)) deallocate(nei_list)
       if(allocated(ncord_new)) deallocate(ncord_new)
      endif
      return 
      end
      subroutine check_save_zone(iflg,izone,zone_check)
c
c check for saved zones
c 
      use comai
      use combi, only : izonesavenum, izone_save, maxsvzone 
      
      implicit none
      
      integer iflg,izone,i 
      logical zone_check
      if(iflg.eq.1) then
        if(.not.allocated(izonesavenum)) 
     &       allocate(izonesavenum(maxsvzone))           
c check for previous zone usage for a saved zone           
          do i = 1, num_sv_zones
           if(izone.eq.izonesavenum(i)) then
c write error messages if a preveous saved zone exists     
            write(ierr,10) izone
            if(iout.ne.0) write(iout,10) izone
            if(iptty.ne.0) write(iptty,10) izone
            stop
           endif
           enddo
c increment the number of saved zones if this is a sved zone        
       if(izone_save.ne.0) then
          num_sv_zones = num_sv_zones + 1
          izonesavenum(num_sv_zones) = izone           
       endif
      else
      endif
10    format('saved zone',1x,i6,1x,'exists, cannot define again',
     &       ' stopping')      
      return
      end
      subroutine zone_saved_manage(iflg,izunit,idz,nin,
     &  n_elem,zone_saved)
c
c  manage saved zone files
c
      use comai
      use combi
      use avsio, only : ioconcentration
      use comdi, only : nspeci
      
      implicit none
      
      integer iflg,i,izunit,idz,nin,n_elem,ns_sv,j,ie,k,k1,k2,k3,ij
      integer izunit1,izunit2,izunit3, max_lines,var_count,n_n,n_e,il
      integer n_n_c,n_e_c, var_node, var_count2, max_line_char,length
      integer open_file
      integer node_temp, n_xyz_node, n_zone, n_zone_max, nzone_cnt
      integer n_zone_last, iop_conc
      character*30 zonesavename, char_temp
      character*200 file_flux, file_scalar
      character*1100 line_temp1, line_temp2
      character*200 string_temp
      logical ex,op,zone_saved
      parameter (max_lines = 100000000, max_line_char = 1000)
      real*8 dum_v1(20), dum_v2(20), dum_var, time_temp
      real*8 time_temp_last
      real*4 caz(2)
      real*8 tajj,tyming
      integer idum_v1(20)
      save tajj
c gaz 080417 might want zone saved in general    
c      if(.not.sv_combine) return
c       
      if(iflg.eq.0) then
      if(.not.sv_combine) return
c call timing function
        if(iout.ne.0) then
         write(iout,100) 
        endif
        if(iptty.ne.0) then
         write(iptty,100) 
        endif
100     format(1x,/,'>> Combining files for SoilVision application  <<')        
        tajj = tyming(caz) 
      elseif(iflg.eq.3) then
        if(.not.sv_combine) return
c write out cpu time for combing files
        tajj = tyming(caz) - tajj
        if(iout.ne.0) then
         write(iout,101) tajj
        endif
        if(iptty.ne.0) then
         write(iptty,101) tajj
        endif
101     format(1x,'CPU time for combining SV files = ',1p,g12.4,' sec')
      else if(iflg.eq.1) then
c open and read saved zone file if they exist
            zonesavename(1:14) = 'zone00000.save'
            write(zonesavename(5:9),'(i5)') idz+10000
            zonesavename(5:5) = '0'
c inquire here and if saved file exists then use it 
c otehrwise revert to old form 
            zone_saved =.false.
            inquire (file = zonesavename, exist = zone_saved)
            if(zone_saved) then
             izunit=open_file(zonesavename,'old')
             read(izunit,*) 
             read(izunit,*)
             if(.not.allocated(ncord)) allocate(ncord(neq))
             if(.not.allocated(ncord_inv)) allocate(ncord_inv(neq))
             read(izunit,*) nin, (ncord(i), i =1, nin)
             ncord_inv = 0
             do i = 1, nin
              ncord_inv(ncord(i)) = i
             enddo
             read(izunit,*)
             read(izunit,*) n_elem, ns_sv
             allocate(elem_temp(n_elem,ns_sv))
             do j = 1, n_elem
              read(izunit,*) ie, (elem_temp(j,i), i = 1, ns_sv)
             enddo
             read(izunit,'(a13)',end = 402) char_temp
             if(char_temp(1:13).eq.'nnum expanded') then
              read(izunit,*) j
              read(izunit,*) (ncord(i), i = nin+1, nin+j)
              do i = nin+1, nin+j
               ncord_inv(ncord(i)) = i
              enddo
              nin = nin +j
             endif
402             close(izunit) 
            endif
      else if(iflg.eq.2) then
c combine SV flux files with scalar files
        if (sv_combine) then
        izunit2 = open_file('soil_vision_beta.dat','unknown')
        do i = 1, icflux
         file_flux = contour_flux_files(i)
          izunit = open_file(file_flux,'old')
          do j = 1, max_line_char
            if(file_flux(j:j+4).eq.'_vec_') then
              file_scalar = file_flux
              file_scalar(j:j+4) = '_sca_'
              ex = .false.
              inquire(file = file_scalar, exist = ex)
              if(ex) then
c now we can combine files  
               izunit1 = open_file(file_scalar,'old')
c modify 2nd line variables 
               if(icnl.eq.0) then                                
                line_temp2(1:37) = '"Liquid X Volume Flux (m3/[m2 s])"'
                line_temp2(38:74) = '"Liquid Y Volume Flux (m3/[m2 s])"'
                line_temp2(75:111) ='"Liquid Z Volume Flux (m3/[m2 s])"'
               else
                line_temp2(1:37) = '"Liquid X Volume Flux (m3/[m2 s])"'
                line_temp2(38:74) = '"Liquid Y Volume Flux (m3/[m2 s])"'
                line_temp2(75:111) ='                                  '
               endif     
c the first countour file contains header information
               if(i.eq.1) then
                read(izunit1,'(a)', end = 501) line_temp1
                write(izunit2,'(a)') line_temp1
                read(izunit1,'(a)', end = 501) line_temp1
                read(izunit,*)
                read(izunit,*)
                var_count = 0
                do k = 1, max_line_char
                 if(line_temp1(k:k).eq.'"') var_count = var_count+1
                 if(line_temp1(k:k+5).eq.'"node"') var_node =var_count
                 if(line_temp1(k:k+3).eq.'"   ') then
                  line_temp1(k+2:k+111) = line_temp2(1:111)
                  go to 500
                 endif
                enddo 
 500            continue
                var_node = var_node/2 + 1
                var_count = var_count/2 
c calculate number of flux variables                
                if(icnl.eq.0) then
                 var_count2 = var_node + 3
                else
                 var_count2 = var_node + 2
                endif
                write(izunit2,'(a)') line_temp1
c  other files all set             
               else
                   
               endif
                do k = 1, max_lines
                  read(izunit1,'(a)', end = 501) line_temp1
                  read(izunit,'(a)', end = 501) line_temp2
                  if(line_temp1(1:4).eq.'ZONE') then
c find node and element numbers 
                   do k1 = 1, max_line_char
                      if(line_temp1(k1:k1+2).eq.'N =') then
                       read(line_temp1(k1+3:k1+11),'(i9)') n_n
                       n_n_c = 1
                      elseif(line_temp1(k1:k1+2).eq.'E =') then
                       read(line_temp1(k1+3:k1+11),'(i9)') n_e
                       n_e_c = 1
                       go to 502
                      endif
                   enddo
502                continue                   
                   write(izunit2,'(a)') line_temp1
                  else
                   backspace izunit1
                   backspace izunit
                   if(n_n_c.le.n_n) then
                    read(izunit1,*, end = 501) 
     &               (dum_v1(k1), k1 = 1, var_count)
                    write(line_temp2,504) 
     &               (dum_v1(k1), k1 = 1, var_count)
                    k2 = 14*(var_node-1)+1
                    ij = dum_v1(var_node)
                    write(line_temp2(k2:k2+13),'(i14)') ij
                    read(izunit,*, end = 501)                     
     &               (dum_v2(k1), k1 = 1, var_count2)
                    k2 = var_count*14
                    if(icnl.eq.0)then
                     k3 = (var_count+3)*14
                     il = var_count2 -2
                    else
                     k3 = (var_count+2)*14
                     il = var_count2 -1
                    endif
                    write(line_temp2(k2+1:k3),503) 
     &                 (dum_v2(k1),k1 =il, var_count2)
                    write(izunit2,'(a)') line_temp2(1:k3)
c                    write(izunit2,503)
c     &              (dum_v1(k1), k1 = 1, var_count),
c     &              (dum_v2(k1), k1 = 4, 5) 
                    n_n_c = n_n_c + 1
503    format(40(1x,1p,g13.6))  
504    format(40(1x,1p,g13.6))                  
                   else
c read and print element information                       
                    read(izunit1,*, end = 501) 
     &               (idum_v1(k1), k1 = 1, ns_in) 
                    read(izunit,*, end = 501) 
     &               (idum_v1(k1), k1 = 1, ns_in)                     
                    write(izunit2,'(50(1x,i10))')
     &              (idum_v1(k1), k1 = 1, ns_in) 
                    n_e_c = n_e_c + 1
                   endif
                  endif
                enddo
              endif  
            endif
          enddo   
501     continue
        enddo
c
       else  
c no Soil Vision files           
        return  
       endif
      else if(iflg.eq.4) then
c add transport output to SV file
       if (sv_combine.and.ioconcentration.eq.1) then
        if(icnl.eq.0) then
         n_xyz_node = 4
        else
         n_xyz_node = 3  
        endif
        izunit2 = open_file('soil_vision_beta.dat','unknown')   
        izunit3 = open_file('soil_vision_beta_conc.dat','unknown')
         read(izunit2,'(a)') line_temp1
         write(izunit3,'(a)') line_temp1
         read(izunit2,'(a)') line_temp1   
         line_temp2 = trim(line_temp1)
c adjust header in SV file         
         write(string_temp,'(a)') ' "Aqueous_Species_000"'
         do i = 1, nspeci
          write(string_temp(19:21),'(i3)') i+100
          string_temp(19:19) = '0'
          line_temp2 = trim(line_temp2) // trim(string_temp)
         enddo
         write(izunit3,'(a)') line_temp2  
c i is time count (different conc files-outer loop)         
         i = 0
         n_zone_max = 0
c nzone_cnt is zone count same conc file - inner loop)          
         nzone_cnt = 0
         time_temp_last = -1.0
         iop_conc = 0
599      continue   
c gaz use soil_vision_beta.dat(izunit2) as the template         
         read(izunit2,'(a)') line_temp1
         write(izunit3,'(a)') trim(line_temp1) 
         read(line_temp1,'(a4,i6,a21,f17.1,a10,i10,a5,i9)')
     &      string_temp(1:4), 
     &      n_zone, string_temp(11:32),time_temp,string_temp (1:10),        
     &      node_temp,string_temp(1:5), n_elem
          if(time_temp.eq.time_temp_last) then
           nzone_cnt = nzone_cnt + 1
          else
            n_zone_max = max(nzone_cnt,n_zone_max)
            time_temp_last = time_temp
            nzone_cnt = 1
            i = i + 1
            iop_conc = 0
          endif          
c gaz 070418 need zone number (first zone,second zone,..)          
          if(n_zone.ne.0) then

            file_scalar = contour_conc_files(i) 
            ex = .false.
            inquire(file = file_scalar, exist = ex)
            if(ex.and.iop_conc.eq.0) then
c now we can combine files  
              izunit1 = open_file(file_scalar,'old')
              iop_conc = 1
            endif
            if(i.eq.1.and.nzone_cnt.eq.1) then
             read(izunit1,'(a)') line_temp1
             read(izunit1,'(a)') line_temp1
            endif
          endif
          read(izunit1,'(a)') line_temp1
            do j = 1, node_temp
             read(izunit1,*,end = 601) (dum_var, k = 1, n_xyz_node),
     &          (dum_v1(k), k = 1, nspeci)
             read(izunit2,'(a)') line_temp1
             length = len_trim(line_temp1) 
             write(line_temp1(length+2:length+202),*) 
     &         (dum_v1(k), k = 1, nspeci)
             write(izunit3,'(a)') trim(line_temp1)
            enddo
c read element list from both files, write element list
            do j = 1, n_elem
              read(izunit1,'(a)') line_temp1
              read(izunit2,'(a)') line_temp1
              write(izunit3,'(a)') trim(line_temp1)
            enddo
c icconc is the number of printout times (in use module comai)  
c gaz 070718            
c            if(i.eq.icconc.and.nzone_cnt.eq.n_zone_max) go to 602
c            if(i.eq.icconc*n_zone_max) go to 602
c gaz debug 070818set i max  = 15            
             if(i.eq.icconc.and.nzone_cnt.eq.n_zone_max) go to 602
          go to 599
c

c
       go to 602
601    continue
       if(iptty.ne.0) write(iptty,*) 'read error in '
     &  ,'subroutine zone_saved_manage. stopping'
       if(iout.ne.0) write(iout,*) 'read error in '
     &  ,'subroutine zone_saved_manage. stopping'  
       stop
602    continue  
c rename file
       close(izunit3)
       call rename('soil_vision_beta_conc.dat','soil_vision.dat')
       continue
c 
      endif
      else if(iflg.eq.-1) then
c delete saved zone files
       do i = 1, izone_sv_cnt
        zonesavename = zonesavenames(i)
        ex = .false.
        op = .false.
        inquire (file = zonesavename, exist = ex)
         if(ex) then
          inquire (file = zonesavename, opened = op)   
          if(.not.op) izunit=open_file(zonesavename,'unknown')
          close(izunit, status = 'delete')
         endif
       enddo
      endif
      return
      end
      
      
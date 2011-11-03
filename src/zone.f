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
CD1
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
      implicit none

      logical null1, null_new, cdum
      integer cnum, i, infile, izone, izonel, nin, nodez, nsl
      integer nxy, icnl_old, nin_old, i_old, ja, jb, jc, izonn
      character* 4 macro, cmacro
      character*80 ltest
      real*8 xg, xz(8), yg, yz(8), zg, zz(8)
      real*8 tol_zone, zxy_min, zxy_max 
      integer imodel, j, n_n_n
      integer zone_dpadd, i3d_2d, i3d_rad

      integer, allocatable :: ncord(:)
      integer, allocatable :: izonef_old(:)

      allocate(ncord(n0))

c     Dual perm or dual porosity value to add to get zones
      zone_dpadd = 100

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
      enddo
      if(i3d_2d.eq.1.and.icnl.eq.0) then
       write(ierr,*) 'i3d_2d parameter ignored for 3d problem'
       if(iout.ne.0) 
     &   write(iout,*) 'i3d_2d parameter ignored for 3d problem'
       if(iptty.ne.0) 
     &   write(iptty,*) 'i3d_2d parameter ignored for 3d problem'     
      endif
      if (macro .eq. 'zone'.or.macro .eq. 'zonn')  then
         cmacro = macro
         if(macro .ne. 'zonn') then
            izonef = 0
            izonn = 0
         else
           allocate(izonef_old(neq))
           izonef_old = izonef
           izonn = 1
         endif
         endif
 60      continue
         read (infile, '(a80)') wdd1
         if (null1(wdd1)) go to 90
         backspace infile
         read(infile, *) izone
c     Determine if zone_dpadd needs to be increased to 1000
         if(izone.gt.99) zone_dpadd = 1000

         if (izone .gt. 0) then
c check if list or nnum occurs
            read(infile, '(a4)') macro
            if (macro .ne. 'list' .and. macro .ne. 'nnum' .and.
     &           macro .ne. 'xyli'.and. macro.ne. 'all ' .and.
     &           macro .ne. 'jajb') then
               backspace  infile
               if (izone .ne. izonel) nin = 0
               if(i3d_rad.eq.1.and.icnl.eq.0) then
c radial input for 3D problems (uses 2D input)
                read (infile, *) (xz(i), i = 1, 4)
                read (infile, *) (yz(i), i = 1, 4)  
                call setzone(izone, nin, ncord, 4, xz, yz, zz, 1)         
               else if(i3d_2d.eq.0.or.icnl.eq.0) then
                read (infile, *) (xz(i), i = 1, nsl)
                read (infile, *) (yz(i), i = 1, nsl)
                if (icnl .eq. 0) then
c**** 3-d calculation ****
                   read (infile, *)  (zz(i), i = 1, nsl)
                end if
                call setzone(izone, nin, ncord, nsl, xz, yz, zz, 0)
               else
c 3-d zones in 2D model (for consistency when extracting slices in 3d)               
                read (infile, *) (xz(i), i = 1, 8)
                xz(3) = xz(2)
                xz(4) = xz(1)                 
                xz(1) = xz(5)
                xz(2) = xz(6)
                read (infile, *) (yz(i), i = 1, 8)
c note we overwrite yz with zz                 
                read (infile, *) (zz(i), i = 1, 8)
                yz(1) = zz(5)
                yz(2) = zz(6)
                yz(3) = zz(2)
                yz(4) = zz(1)
                call setzone(izone, nin, ncord, nsl, xz, yz, zz, 0)  
               endif
            else if(macro .eq. 'xyli') then
c read in nodes in zone from xy list
               nxy = 0
               nin = 0
               i = 0
               read(infile,*) tol_zone, zxy_min, zxy_max
 71            read(infile, '(a80)') ltest
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
                 call setzone(izone, nin_old, ncord(nin+1:neq),
     &                nsl, xz, yz, zz, 0)
                 nin=nin+nin_old
               go to 71
 81            continue
            else if(macro .eq. 'list') then
c read in coordinates for nodes in zone
               nin = 0
 70            read(infile, '(a80)') ltest
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
 80            continue
            else if(macro .eq. 'nnum') then
c read in nodes belonging to zone
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
              do i = 1, n0
               izonef(i) = izone
              enddo         
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
 6008       format (' **** Invalid input: macro ', a4, ' ****')
 6009       format(' **** Invalid node specified, ', i8, 
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
 6010          format(/, 1x, i8,' nodes contained in zone = ',
     &              i10, /)
               write(ischk, 6011) (ncord(i), i = 1, nin)
 6011          format (10i8)
            end if
            izonel = izone
            go to 60
         endif
c check which nodes don't belong to a zone
 90      nin = 0
c     Changed to neq_primary (used to be neq) BAR - 12-15-99
         do i = 1, neq_primary
            if(izonef(i) .eq. 0) then
               nin = nin + 1
               ncord(nin) = i
            endif
         end do
         if (nin .ne. 0) then
            if (ischk .ne. 0) write(ischk, 6012) nin, cnum
c            write(ischk, 6011) (ncord(i), i = 1, nin)
         end if
 6012    format(/, 1x, i8,
     &        ' nodes not assigned to a zone in call # ', i10)


c     Assign zones for GDPM nodes

         if(gdpm_flag.ne.0) then

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
                   endif       
               end do 
            end do
         end if



c check for dual porosity or dpdp solution
         if(idualp .eq. 1) then
c dual porosity solution
            n = neq+neq+neq
            if (ischk .ne. 0) then
               write(ischk, *) 'dual porosity solution'
               write(ischk, *) 'first matrix level zone = ',
     .              'fracture level zone + ',zone_dpadd
               write(ischk, *) 'second matrix level zone = ', 
     .              'fracture level zone + ',zone_dpadd*2
            end if

c	This loop changed to set zones to their value plus 100
c	only if the nodes have not been explicitly set in the
c	zone definition. This allows the user to set the matrix
c	nodes in the zone macro and not have the code default
c	to zone number plus 100 (or 200).
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
         else if(idpdp .ne. 0) then
c dpdp solution
            n = neq+neq
            if (ischk .ne. 0) then
               write(ischk, *) 'dual porosity/dual permeability ',
     &              'solution'
               write(ischk, *) 'first matrix level zone = ',
     .              'fracture level zone + ',zone_dpadd
            end if
c	This loop changed to set zones to their value plus 100
c	only if the nodes have not been explicitly set in the
c	zone definition. This allows the user to set the matrix
c	nodes in the zone macro and not have the code default
c	to zone number plus 100.
c     Loop changed to accomodate the new zone_dpadd variable

            do i = 1, neq
               if(izonn.eq.1)then
                if(izonef(i).ne.izonef_old(i)) then
                  izonef(i + neq) = izonef(i) + zone_dpadd
                 endif
               else
                if(izonef(i+neq).eq.0) then
                 izonef(i + neq) = izonef(i) + zone_dpadd
                endif
               end if
            end do
         else
            n = neq
         endif
         go to 100
c      end if
      go  to  50
 100  continue

      macroread(18) = .TRUE.
      deallocate(ncord)
      if(macro .eq. 'zonn') deallocate(izonef_old)
      return
      end
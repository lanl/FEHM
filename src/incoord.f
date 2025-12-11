      subroutine incoord
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
CD1 Control reading of input coordinate data.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 20-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/incoord.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.2   Tue Jan 30 10:37:14 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:03:22   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:24:50   pvcs
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
CD3   incoor                   I    Coordinate input file
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
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   incoor          INT      faai   Unit number for coordinate input file
CD4   mlz             INT      faai   Out of bounds node
CD4   nei             INT      faai   Total number of elements in the problem
CD4   nelm            INT      fbb    Initially information about nodes in each
CD4                                     element, later nodal connectivity
CD4                                     information
CD4   neq             INT      faai   Number of nodes
CD4   ns              INT      faai   Number of nodes per element
CD4   wdd1            CHAR     faac   Character input string
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   gendat                   Generate mesh through interpolation of input
CD4                              coordinates
CD4   geoin                    Read in alternate geometric data
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
CD5   mb              INT      Node number
CD5   mc              INT      Node/element count for interpolation of
CD5                              coordinate data by gendat
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
CPS BEGIN incoord
CPS 
CPS   REPEAT
CPS   
CPS     read line from input file
CPS     
CPS     IF the coor macro is read
CPS        write macro identifier and input file unit to output file 
CPS        and tty if enabled
CPS        set initial node count to 1
CPS        read in number of nodes
CPS        
CPS        REPEAT
CPS          read a line from the input file
CPS          EXIT IF the line is null (no more coordinate data)
CPS          reread the input line with coordinate format
CPS          IF the node number is negative
CPS             call gendat to generate coordinate data through 
CPS             interpolation
CPS          END IF
CPS          set node count to absolute value of current node number
CPS        UNTIL there is no more coordinate data to be read
CPS        
CPS     ELSE IF the elem macro is read
CPS        write macro identifier and input file unit to output file 
CPS        and tty if enabled
CPS     
CPS        set initial element count to 1
CPS        read in number of nodes per element and total number of nodes
CPS        IF the number of nodes per element is negative
CPS           set the out of bounds flag to 1
CPS           take the absolute value of the number of nodes per element
CPS        END IF
CPS        
CPS        REPEAT
CPS          read a line from the input file
CPS          EXIT IF the line is null (no more element data)
CPS          reread the input line with element format
CPS          IF the element number is negative
CPS             call gendat to generate element data through 
CPS             interpolation
CPS          END IF
CPS          set element count to absolute value of current node number
CPS        UNTIL there is no more element data to be read
CPS        
CPS     ELSE IF the alti macro is read
CPS        write macro identifier and input file unit to output file 
CPS        and tty if enabled
CPS        call geoin to read in geometric data with alternate file 
CPS        format
CPS     ELSE IF stop macro is read
CPS        write macro identifier and input file unit to output file 
CPS        and tty if enabled
CPS        EXIT repeat loop
CPS     END IF
CPS     
CPS   UNTIL stop is found
CPS   
CPS END incoord
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      use comxi
      implicit none

      logical null1
      integer i, j, nodi, mb, mc, istat, ielemorder 
      integer ic, nc1, nc2, nc3, nc4
      character*4 macro
      character*4 elem_order
      character*80 dum_char
c parser variables
      character*80 input_msg
      integer msg(2)
      integer nwds,sehtemp
      real*8 xmsg(2)
      integer imsg(2)
      character*32 cmsg(2)

! Is it an unformatted file?
         
      if (cform(3)(1:3).eq.'UNF' .or. cform(3)(1:3).eq.'unf') then
! Read unformatted coordinate file
!coordinate data
         read (incoor, ERR=999, IOSTAT=istat) macro
         if (iout .ne. 0) write(iout, 6010) macro, 'incoor', incoor
         if (iptty .gt. 0) write(iptty, 6010) macro, 
     *        'incoor', incoor
         read (incoor, ERR=999, IOSTAT=istat) neq_primary
         read (incoor, ERR=999, IOSTAT=istat) (nodi, cord(i,1), 
     *        cord(i,2), cord(i,3), i=1,neq_primary)
!blank line separator
         read (incoor, ERR=999, IOSTAT=istat) macro
!element data
         read (incoor, ERR=999, IOSTAT=istat) macro
         if (iout .ne. 0) write(iout, 6010) macro, 'incoor', incoor
         if (iptty .gt. 0) write(iptty, 6010) macro, 'incoor', 
     *        incoor
         read (incoor) ns, nei
         ns_in = ns
         nei_in = nei
         do j = 1, nei
            read (incoor, ERR=999, IOSTAT=istat) mb, (nelm(
     *           (iabs(mb) - 1) * ns + i), i = 1, ns)
         end do
!blank line separator
         read (incoor, ERR=999, IOSTAT=istat) macro
 999     if (istat .ne. 0) then
            write (ierr, *) 'Error reading unformatted ',
     *           'coordinate file'
            if (iptty .ne. 0) write (ierr, *) 'Error reading ',
     *           'unformatted coordinate file'
            close (incoor)
            stop
         end if
      else
! Read formatted coordinate file

 10      continue

         read (incoor, '(a4)', end = 40) macro

         if (macro .eq. 'coor') then
c**** node coordinate data ****
            if (iout .ne. 0) write(iout, 6010) macro, 'incoor', incoor
            if (iptty .gt. 0) write(iptty, 6010) macro, 
     *           'incoor', incoor
 6010       format(1x, '**** input title : ', a4, ' **** ', a6, 
     *           ' = ', i3, ' ****')

            mc = 1
            read (incoor, *) neq_primary

 20         continue
            read (incoor, '(a80)') wdd1
            if (null1(wdd1)) go to 25 
            backspace incoor
            read (incoor, *) mb, cord(iabs(mb), 1),  
     *           cord(iabs(mb), 2), cord(iabs(mb), 3)
            if (mb .lt. 0) call gendat(iabs(mb), mc, 1)
            mc = iabs(mb)
            go  to  20
 25         continue
         else if (macro .eq. 'elem') then
c**** element node data ****

            backspace incoor
            read (incoor,'(a80)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
            macro = cmsg(1)
            if (nwds .lt. 2) then
               elem_order = '    '
            else
               if (msg(2) .ne. 3) then
                  elem_order = '    '
               else
                  elem_order = cmsg(2)
               end if
            end if

            if (iout .ne. 0) write(iout, 6010) macro, 'incoor', incoor
            if (iptty .gt. 0) write(iptty, 6010) macro, 'incoor', 
     *           incoor

            mc = 1

            read (incoor, *) ns, nei
            ns_in = ns
            nei_in = nei

            if(ns .lt. 0) then
               mlz = 1
               ns = iabs(ns)
            else
c gaz debug 032219                
c               mlz = 0 gaz changed back "no change of mlz"
            end if

            if(elem_order.eq.'trad'.or.ns.ne.8) then
             ielemorder = 0
            else
             ielemorder = 1
            endif
            if(elem_order.eq.'trad'.and.ns.eq.4.and.icnl.eq.0) then
             ich_pebi = 1
            else
             ich_pebi = 0
            endif
            
 30         continue

            read (incoor, '(a80)') wdd1
            if (null1(wdd1)) go to 35 

            backspace incoor
            read(incoor, *) mb, 
     *           (nelm((iabs(mb) - 1) * ns + i), i = 1, ns)
            if (mb .lt. 0) call gendat (iabs(mb), mc, 2)
            mc = iabs(mb)
            go  to  30
 35         continue
c reverse order for non traditional hexes
c exclude mixed in prisms (nelm(ic+7)=0)
            if(ielemorder.eq.1) then
             do j = 1,nei              
               ic = (j-1)*ns 
               if(nelm(ic+7).ne.0) then
                nc1 = nelm(ic+1)
                nc2 = nelm(ic+2)
                nc3 = nelm(ic+3)
                nc4 = nelm(ic+4)
                nelm(ic+1) = nelm(ic+5)
                nelm(ic+2) = nelm(ic+6)
                nelm(ic+3) = nelm(ic+7)
                nelm(ic+4) = nelm(ic+8)
                nelm(ic+5) = nc1
                nelm(ic+6) = nc2
                nelm(ic+7) = nc3
                nelm(ic+8) = nc4
               endif 
             enddo
            endif         
         else if (macro .eq. 'alti') then
c**** alternate element input ****

            if (iout .ne. 0) write(iout, 6010) macro, 'incoor', incoor
            if (iptty .gt. 0) write(iptty, 6010) macro, 'incoor', 
     *           incoor

            call geoin

         else if (macro .eq. 'fdm ') then
c**** finite difference block-centered input

            call structured(0)

c RJP 12/13/06 modified below
c            go  to  40

c         else if (macro .eq. 'rive' .or. macro .eq. 'well') then
c**** river or well input
c            nriver = 1
c            call river_ctr(0)
c            call river_ctr(1)

         else if (macro .eq. 'stop') then
            
            if (iout .ne. 0) write(iout, 6010) macro, 'incoor', incoor
            if (iptty .gt. 0) write(iptty, 6010) macro, 'incoor', 
     *           incoor

            go  to  40

         end if

         go  to  10

 40      continue
      end if

      end

      subroutine innode (macro)
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
CD1 Read/find node numbers for output.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 22-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/innode.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:56   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Wed Jun 12 15:21:06 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.8   Mon Jun 03 11:18:00 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.7   Fri May 31 10:53:30 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.6   Tue Jan 30 13:06:52 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.5   08/22/95 14:58:58   llt
CD2 m2 was defined twice, removed for Cray
CD2 
CD2    Rev 1.4   08/22/95 13:48:40   llt
CD2 quotes around internal read name doesn't work on IBM
CD2 
CD2    Rev 1.3   08/18/95 10:28:08   llt
CD2 m was already defined, removed for cray
CD2 needed quotes around wddl for cray to read
CD2 
CD2    Rev 1.2   08/07/95 11:48:06   awolf
CD2 Now allows specifying node macro by blocks or zones
CD2 
CD2    Rev 1.1   03/18/94 16:03:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:08   pvcs
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
CD3   macro           CHAR     I    Macro data being read
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   inpt                     I    Input data file
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
CD4   inpt            INT      faai   Unit number for input file
CD4   nskw            INT      fddi   Contains nodes for print-out
CD4   nskw2           INT      fddi   Contains nodes for tty print-out
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   near3                    Determine the nearest node to a given set of
CD4                              coordinates (x, y, z)
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
CD5   m               INT      Number of nodes to print to output file
CD5   m2              INT      Number of nodes to print to tty
CD5   nodew           INT      Node to be printed
CD5   xc              REAL*8   X coordinate for node
CD5   yc              REAL*8   X coordinate for node
CD5   zc              REAL*8   X coordinate for node
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
CPS BEGIN innode
CPS 
CPS   set number of nodes for printout and tty output to 0
CPS   IF macro being read is node 
CPS      read number of nodes for printout
CPS   ELSE
CPS      read number of nodes for printout and number for tty output
CPS   END IF
CPS   
CPS   IF number of nodes for printout is > 0
CPS      read node numbers
CPS      FOR each node read
CPS          IF node number is negative
CPS             read node coordinates
CPS             call near3 to find node to printout
CPS          END IF
CPS      END FOR
CPS   END IF
CPS   
CPS   IF number of nodes for tty output is > 0
CPS      read node numbers
CPS      FOR each node read
CPS          IF node number is negative
CPS             read node coordinates
CPS             call near3 to find node to printout
CPS          END IF
CPS      END FOR
CPS   END IF
CPS   
CPS END innode
CPS
C***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      implicit none

      logical null1
      real*8 xmsg(4)
      integer i, cnt, cnt2, cnt3, ja, jb, jc, ja_tmp, ja_match
c gaz 090523
      integer i1, i2, ii, nin
      integer nflag, msg(4), imsg(4), nwds, ierr_flag
      character*4 macro
      character*5 myform
      character*32 cmsg(4)
      character*80 chdum

      ierr_flag = 0
c**** read node numbers for output ****
      if(macro .eq. 'node') then
         nflag = 1
         cnt = m
      else if (macro .eq. 'nod2') then
         nflag = 2
         cnt = m2
      else if (macro .eq. 'nod3') then
         nflag = 3
         cnt = m3
      end if

      read(inpt,'(a5)') myform
      select case (myform)
      case ('block', 'bloc')
         do
            read(inpt,'(a80)') chdum
            if (null1(chdum)) exit
            read(chdum,*) ja,jb,jc
            if (ja.lt.0) then
c we are dealing with zone format
                 ja=abs(ja)   
                 ja_tmp =ja
                 call zone_saved(2,'node',ja_tmp,1, nin) 
               if(nin.eq.0) then
c not a saved zone
                i1 = 1
                i2 = n0
               else 
c a saved zone
                i1 = 1
                i2 = nin
               endif
               do ii = i1, i2
c gaz 101223 use nin = 0 to indicate zone is not saved
                  if(nin.eq.0) then
                   i = ii
                   ja_match = izonef(i)
                  else
                   i = ncord(ii)
                   ja_match = ja
                  endif
                  if (ja_match.eq.ja) then
                     cnt = cnt + 1
                     select case (nflag)
                     case (1)
                        nskw(cnt)=i
                     case (2)
                        nskw2(cnt)=i
                     case (3)
                        nskw3(cnt)=i
                     end select
                  endif
               enddo
            else
c nodes are in block format
               if (jb .eq. 0 .and. jc .eq. 0) then
                  jb = n0
                  jc = 1
               end if
               do i=ja,jb,jc
                  cnt = cnt + 1
                  select case (nflag)
                  case (1)
                     nskw(cnt)=i
                  case (2)
                     nskw2(cnt)=i
                  case (3)
                     nskw3(cnt)=i
                  end select
               enddo
            endif
         end do
         select case (nflag)
         case (1)
            m = cnt
         case (2)
            m2 = cnt
         case (3)
            m3 = cnt
         end select
         return
      case ('azone','azon')
c**** read in history by zone output selections ****
         out_zones = .TRUE.
         call inhstz
         return
      case default
         backspace inpt
         read(inpt,'(a80)') chdum
         call parse_string(chdum,imsg,msg,xmsg,cmsg,nwds)
         cnt = imsg(1) + xmsg(1)
         select case (nflag)
         case (1)
            if (cnt .gt. 0) then
               call get_nodes (nskw, cnt, m, ierr_flag)
               m = m + cnt
            else
               ierr_flag = 1
            end if
         case (2)
            if (nwds .ge. 2) then
! Old input format
               cnt2 = imsg(2) + xmsg(2)      
               if (cnt .gt. 0) then
                  call get_nodes (nskw, cnt, m, ierr_flag)
                  m = m + cnt 
               end if
               if (cnt2 .gt. 0) then
                  call get_nodes (nskw2, cnt2, m2, ierr_flag)
                  m2 = m2 + cnt2
               else
                  ierr_flag = 1
               end if
            else
               cnt2 = cnt
               if (cnt2 .gt. 0) then
                  call get_nodes (nskw2, cnt2, m2, ierr_flag)
                  m2 = m2 + cnt
               else
                  ierr_flag = 1
               end if
            end if
         case (3)
            if(iporos.ne.-4) then
               if (iout .ne. 0) write(iout,200)
               if (iptty .gt. 0) write(iptty,200)
               if(.not.allocated(nskw3)) allocate(nskw3(n0))
            endif
            if (nwds .ge. 3) then
! Old input format
               cnt2 = imsg(2) + xmsg(2)
               cnt3 = imsg(3) + xmsg(3)
               if (cnt .gt. 0) then
                  call get_nodes (nskw, cnt, m, ierr_flag)
                  m = m + cnt
               end if
               if (cnt2 .gt. 0) then
                  call get_nodes (nskw2, cnt2, m2, ierr_flag)
                  m2 = m2 + cnt2
               end if
               if (cnt3 .gt. 0) then
                  call get_nodes (nskw3, cnt3, m3, ierr_flag)
                  m3 = m3 + cnt3
               else
                  ierr_flag = 1
               end if
            else
               cnt3 = cnt
               if (cnt3 .gt. 0) then
                  call get_nodes (nskw3, cnt3, m3, ierr_flag)
                  m3 = m3 + cnt3
               else
                  ierr_flag = 1
               end if
            end if
         end select
      end select

! Input error
      if (ierr_flag .eq. 1) then
         write (ierr, 500) macro
         if (iout .ne. 0) write (iout, 500) macro
         if (iptty .gt. 0) write (iptty, 500) macro
      else if (ierr_flag .eq. 2) then
         stop
      end if
      
 200  format(/
     &     ,'>>> Warning: printout for ppor model -4 requested',/
     &     ,'>>> But no model,Input for nod3 macro read but not used',/)
 500  format ('Warning: Invalid number of nodes specified for macro ', 
     &     a4)

      end subroutine innode

      subroutine get_nodes (node_array, cnt, m, ierr_flag)

      use comai, only : inpt, ierr, iout, iptty, gdkm_flag
      use comdti, only : n0
      use comki, only : macro
      implicit none
      integer i, nodew, cnt, m, ierr_flag
      integer node_array(*)
      real*8 xc, yc, zc

      read (inpt, *) (node_array(i), i = m + 1, m + cnt)
c**** read in coordinates if node_array(i) < 0 and find node number ****
!****  or Stop if node number is outside of problem domain ****
      do  i = m + 1, m + cnt
         if (node_array(i) .lt. 0) then
            read (inpt, *) xc, yc, zc
            call near3 (xc, yc, zc, nodew, 0)
            node_array(i) =  nodew
c gaz 022717            
         else if (node_array(i) .gt. n0. and . gdkm_flag.eq.0) then
            ierr_flag = 2
            write (ierr, 400) macro
            write (ierr, 300) node_array(i), n0
            if (iout .ne. 0) write (iout, 400) macro
            if (iout .ne. 0) write (iout, 300) node_array(i), n0
            if (iptty .gt. 0) write (iptty, 400) macro
            if (iptty .gt. 0) write (iptty, 300) node_array(i), n0
         end if
      end do

 400  format (' **** Invalid input: macro ', a4, ' ****')
 300  format(' **** Invalid node specified,  ', i8, 
     .        ' is greater than ', 'n0 (', i8, ' ): stopping ****')
      end subroutine get_nodes
      

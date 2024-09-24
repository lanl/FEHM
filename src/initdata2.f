      subroutine initdata2( in_number,
     2     out_number,
     3     npoints,
     4     narrays,
     6     itype,
     7     default,
     8     readflag,
     9     macro,
     t     igroup,
     1     ireturn, r8_1, r8_2, r8_3, r8_4, r8_5,i4_1, 
     2     i4_2, i4_3, i4_4, i4_5 )
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To read in an arbitrary number of lines of data and set parameter
CD1 values at given nodes.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-14-93     B. Robinson    00022   Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/initdata2.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Tue Jan 30 12:49:24 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.7   Fri Jan 12 17:50:32 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.6   09/07/95 08:58:14   gaz
CD2 added format statement(9004) for integer outpu
CD2 
CD2    Rev 1.5   11/28/94 14:20:56   llt
CD2 Changed .eq. to .eqv. in logical statement, so would run on ibm.
CD2 
CD2    Rev 1.4   06/20/94 11:09:00   zvd
CD2  
CD2 
CD2    Rev 1.3   03/18/94 15:55:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.2   02/14/94 11:12:20   zvd
CD2  
CD2 
CD2    Rev 1.1   02/14/94 11:04:48   zvd
CD2  
CD2 
CD2    Rev 1.0   01/20/94 10:25:04   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier        Type     Use  Description
CD3 
CD3 in_number         int       I   Tape number of input file
CD3 out_number        int       I   Tape number for writing error
CD3                                     checking information
CD3 npoints           int       I   Number of nodes at which values
CD3                                    are set
CD3 narrays           int       I   Number of values to read (i.e.
CD3                                    number of arrays being
CD3                                    initialized)
CD3 itype             int       I   Flag denoting whether array is
CD3                                     integer or real*8
CD3 default           int       I   Default value to set for all
CD3                                     positions in array
CD3 readflag          logical   I   Flag denoting if data
CD3                                     initialization is being
CD3                                     performed
CD3 macro             char      I   Macro being read
CD3 igroup            int       I   Group number within macro
CD3 ireturn           int       O   Return flag
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 Identifier    Use  Description
CD3
CD3 in_number      I   Main input file containing the data
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4 
CD4 Identifier   Type   Description
CD4 
CD4 welbor       N/A    Compute wellbore solution
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8 
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN initdata
CPS 
CPS allocate storage in temporary arrays
CPS initialize values in arrays to 0
CPS 
CPS FOR each array
CPS 
CPS   IF this array is a real*8
CPS     Set pointer for real*8 array
CPS     
CPS     IF initialization is to be performed
CPS       FOR each node
CPS           Set parameter values to their default values
CPS       ENDFOR
CPS     ENDIF
CPS   
CPS   ELSEIF this array is an integer
CPS     Set pointer for integer array
CPS     
CPS     IF initialization is to be performed
CPS       FOR each node
CPS         Set parameter values to their default values
CPS       ENDFOR
CPS     ENDIF
CPS   
CPS   ELSE there is an error
CPS     Set error flag to indicate unsuccessful read
CPS     ERROREXIT
CPS   ENDIF
CPS   
CPS ENDFOR each array
CPS   
CPS Initialize count of number of lines to 0
CPS 
CPS LOOP to read in a group of data
CPS     
CPS   null1 - determine if there is more data to read
CPS       
CPS EXITIF there was no data to read
CPS     
CPS   Read next line of data
CPS   
CPS   EXITIF the first loop index is 0 (no more data, old input format)
CPS   
CPS   EXITIF there is a read error after setting error flag
CPS   
CPS   Add one to the running total of the number of lines read
CPS   
CPS   IF input is by zone
CPS     FOR each node
CPS       IF this node is in this zone
CPS         Set flag to indicate parameter value has been...
CPS         ... explicitly set
CPS         FOR each array
CPS 
CPS           IF this array is a real*8
CPS           
CPS             Set pointer for real*8 array
CPS             Set parameter value
CPS           
CPS           ELSE this array is an integer
CPS             Set pointer for integer array
CPS           
CPS             Set pointer for integer array
CPS             Set parameter value
CPS           
CPS           ENDIF
CPS   
CPS         ENDFOR
CPS         
CPS       ENDIF
CPS     ENDFOR
CPS   ELSE input is by node
CPS   
CPS     IF the values are to be set at every node
CPS       Reset do loop indexes
CPS     ENDIF
CPS     
CPS     FOR each specified node
CPS     
CPS       IF this node number is within the range of possible node...
CPS       ... numbers
CPS         Set flag to indicate parameter value has been...
CPS         ... explicitly set
CPS         FOR each array
CPS 
CPS           IF this array is a real*8
CPS             Set pointer for real*8 array
CPS             Set parameter value
CPS             Set flag to indicate parameter value has been...
CPS             ... explicitly set 
CPS           ELSE this array is an integer
CPS             Set pointer for real*8 array
CPS             Set parameter value
CPS             Set flag to indicate parameter value has been...
CPS             ... explicitly set 
CPS           ENDIF
CPS         
CPS         ENDFOR
CPS         
CPS       ENDIF
CPS         
CPS     ENDFOR
CPS   ENDIF
CPS   
CPS ENDLOOP
CPS 
CPS Initialize node counters
CPS 
CPS FOR each node
CPS   IF the value was not explicitly set at this node
CPS     Store node number
CPS     Increase counter by 1
CPS   ELSE the value was explicitly set
CPS     Store node number
CPS     Increase counter by 1
CPS   ENDIF
CPS ENDFOR
CPS 
CPS IF all values were explicitly set
CPS   Write message
CPS ELSEIF no values were explicitly set
CPS   Write message
CPS ELSEIF more values were set than not set
CPS   Write node numbers at which values were not explicitly set
CPS ELSE more values were not set than set
CPS   Write node numbers at which values were explicitly set
CPS ENDIF
CPS   
CPS ERRORSEGMENT
CPS   deallocate space in temporary arrays
CPS   IF there was an error in specifying the data types or reading data
CPS      Write error message to error output
CPS      IF tty output is being used
CPS         Write error message to tty
CPS      END IF
CPS      IF the error was in specifying data type
CPS         Write error message to error output
CPS         IF tty output is being used
CPS            Write error message to tty
CPS         END IF
CPS      ELSE IF there was a read error
CPS         Write error message to error output
CPS         IF tty output is being used
CPS            Write error message to tty
CPS         END IF
CPS      END IF
CPS      terminate program
CPS   ENDIF
CPS ENDSEGMENT
CPS 
CPS END initdata
CPS 
C**********************************************************************

      use combi
      use comdti
      use comai
      implicit none

      logical null1,readflag
      integer inode,max_arrays
      integer i,ii,izunit,nin
      integer open_file
      parameter(max_arrays = 10 )
      integer in_number,out_number,npoints,narrays,ireturn
      integer itype(*),iarray,ipoint,inumber,ja,jb,jc,icode,nfound
      integer nnotfound,igroup
      real*8 default(max_arrays),values(max_arrays)
      integer n_realcount, n_intcount
      character*80 strtot
      character*4 macro
      character*30 zonesavename
      integer, allocatable :: ifind(:)
      integer, allocatable :: isset(:)
      integer, allocatable :: notset(:)
      real*8, optional :: r8_1(:)
      real*8, optional :: r8_2(:)
      real*8, optional :: r8_3(:)
      real*8, optional :: r8_4(:)
      real*8, optional :: r8_5(:)
      integer, optional :: i4_1(:)
      integer, optional :: i4_2(:)
      integer, optional :: i4_3(:)
      integer, optional :: i4_4(:)
      integer, optional :: i4_5(:)
      logical ex,op

      allocate(ifind(npoints),isset(npoints),notset(npoints))
      ifind=0
      isset=0
      notset=0
      
c     initialize integer and real counters

      n_realcount = 0
      n_intcount = 0

c     Set values to default values for all positions in each array

      do iarray = 1, narrays
         if( itype(iarray) .eq. 8 ) then
            n_realcount = n_realcount + 1
            if( .not. readflag ) then
               if(n_realcount.eq.1) then
                  r8_1(1:npoints) = default(iarray)
               elseif(n_realcount.eq.2) then
                  r8_2(1:npoints) = default(iarray)
               elseif(n_realcount.eq.3) then
                  r8_3(1:npoints) = default(iarray)
               elseif(n_realcount.eq.4) then
                  r8_4(1:npoints) = default(iarray)
               elseif(n_realcount.eq.5) then
                  r8_5(1:npoints) = default(iarray)
               else
                  write (ierr, *) 'Fatal error, too many real inputs ',
     &                 'to initdata2'
                  if (out_number .ne. 0) then
                     write(out_number,*) 'Fatal error, too many'
                     write(out_number,*) 'real inputs to initdata2'
                  end if
                  stop
               end if
            end if
         elseif( itype(iarray) .eq. 4 ) then
            n_intcount = n_intcount + 1
            if( .not. readflag ) then
               if(n_intcount.eq.1) then
                  i4_1(1:npoints) = nint(default(iarray))+initdata_pad
               elseif(n_intcount.eq.2) then
                  i4_2(1:npoints) = nint(default(iarray))
               elseif(n_intcount.eq.3) then
                  i4_3(1:npoints) = nint(default(iarray))
               elseif(n_intcount.eq.4) then
                  i4_4(1:npoints) = nint(default(iarray))
               elseif(n_intcount.eq.5) then
                  i4_5(1:npoints) = nint(default(iarray))
               else
                  write (ierr, *) 'Fatal error, too many integer ',
     &                 'inputs to initdata2'
                  if (out_number .ne. 0) then
                     write(out_number,*) 'Fatal error, too many'
                     write(out_number,*) 'integer inputs to initdata2'
                  end if
                  stop
               end if
            end if
         else
            ireturn = -1
            goto 9000
         end if
      end do

      ireturn = 0

c     Loop over all lines of data to read
 1000 continue
         read(in_number, '(a80)') strtot
         if( null1(strtot) .eqv. .TRUE. ) goto 2000
         backspace in_number
         read(in_number, *, ERR = 3000) ja, jb, jc, 
     2        (values(iarray), iarray = 1, narrays )
         if (ja .eq. 0) goto 2000
         goto 4000
 3000    continue
         inumber = ireturn + 1
         ireturn = -2
         goto 9000
 4000    continue
         ireturn = ireturn + 1
c     Input by zones is first, or else input is by node
         if( ja .lt. 0 ) then
c gaz 060617         
c check for saved zonefile             
          zonesavename(1:14) = 'zone00000.save'
          write(zonesavename(5:9),'(i5)') abs(ja)+10000
          zonesavename(5:5) = '0'
          ex = .false.
          op = .false.
          inquire (file = zonesavename, exist = ex)
          if(ex) then
           inquire (file = zonesavename, opened = op)   
           if(.not.op) then
            izunit=open_file(zonesavename,'unknown')
           else
            inquire(file = zonesavename, number=izunit)
            close(izunit)
            izunit=open_file(zonesavename,'unknown')
           endif
           read(izunit,*)  ja
           read(izunit,*)  wdd
           read(izunit,*) nin
           if(allocated(ncord)) deallocate(ncord)
           allocate(ncord(nin))
           backspace izunit
           read(izunit,*) nin, (ncord(i), i =1, nin)
           close(izunit)
          endif  
          if(.not.ex) then
            do inode = 1, npoints
               if( izonef(inode) .eq. abs(ja) ) then
                  ifind(inode) = 1
                  n_realcount = 0
                  n_intcount = 0
                  call zone_char_fill(1,inode,abs(ja))
                  do iarray = 1, narrays
                     if( itype(iarray) .eq. 8 ) then
                        n_realcount = n_realcount + 1
                        if(n_realcount.eq.1) then
                           r8_1(inode) = values(iarray)
                        elseif(n_realcount.eq.2) then
                           r8_2(inode) = values(iarray)
                        elseif(n_realcount.eq.3) then
                           r8_3(inode) = values(iarray)
                        elseif(n_realcount.eq.4) then
                           r8_4(inode) = values(iarray)
                        elseif(n_realcount.eq.5) then
                           r8_5(inode) = values(iarray)
                        end if
                     else
                        n_intcount = n_intcount + 1
                        if(n_intcount.eq.1) then
                           i4_1(inode) = nint(values(iarray))
     &                      +initdata_pad
                        elseif(n_intcount.eq.2) then
                           i4_2(inode) = nint(values(iarray))
                        elseif(n_intcount.eq.3) then
                           i4_3(inode) = nint(values(iarray))
                        elseif(n_intcount.eq.4) then
                           i4_4(inode) = nint(values(iarray))
                        elseif(n_intcount.eq.5) then
                           i4_5(inode) = nint(values(iarray))
                        end if
                     end if
                  end do
               end if
c gaz 103118 set to             
            end do
           else
c gaz 111716  
            do ii = 1, nin
               inode = ncord(ii)
                  ifind(inode) = 1
                  n_realcount = 0
                  n_intcount = 0
                  do iarray = 1, narrays
                     if( itype(iarray) .eq. 8 ) then
                        n_realcount = n_realcount + 1
                        if(n_realcount.eq.1) then
                           r8_1(inode) = values(iarray)
                        elseif(n_realcount.eq.2) then
                           r8_2(inode) = values(iarray)
                        elseif(n_realcount.eq.3) then
                           r8_3(inode) = values(iarray)
                        elseif(n_realcount.eq.4) then
                           r8_4(inode) = values(iarray)
                        elseif(n_realcount.eq.5) then
                           r8_5(inode) = values(iarray)
                        end if
                     else
                        n_intcount = n_intcount + 1
                        if(n_intcount.eq.1) then
                           i4_1(inode) = nint(values(iarray))
     &                      +initdata_pad
                        elseif(n_intcount.eq.2) then
                           i4_2(inode) = nint(values(iarray))
                        elseif(n_intcount.eq.3) then
                           i4_3(inode) = nint(values(iarray))
                        elseif(n_intcount.eq.4) then
                           i4_4(inode) = nint(values(iarray))
                        elseif(n_intcount.eq.5) then
                           i4_5(inode) = nint(values(iarray))
                        end if
                     end if
                  end do
c                end if
            end do   
            deallocate(ncord)
          endif
         else
            if( ja .eq. 1 .and. jb .eq. 0 .and. jc .eq. 0 ) then
               ja = 1
               jb = npoints
               jc = 1
               call zone_char_fill(2,1,npoints)
            end if
            if( jb .gt. npoints ) jb = npoints
            ifind(ja:jb:jc) = 1
            n_realcount = 0
            n_intcount = 0
            
            do iarray = 1, narrays
               if( itype(iarray) .eq. 8 ) then
                  n_realcount = n_realcount + 1
                  if(n_realcount.eq.1) then
                     r8_1(ja:jb:jc) = values(iarray)
                  elseif(n_realcount.eq.2) then
                     r8_2(ja:jb:jc) = values(iarray)
                  elseif(n_realcount.eq.3) then
                     r8_3(ja:jb:jc) = values(iarray)
                  elseif(n_realcount.eq.4) then
                     r8_4(ja:jb:jc) = values(iarray)
                  elseif(n_realcount.eq.5) then
                     r8_5(ja:jb:jc) = values(iarray)
                  end if
               else
                  n_intcount = n_intcount + 1
                  if(n_intcount.eq.1) then
                     i4_1(ja:jb:jc) = nint(values(iarray))+initdata_pad
                  elseif(n_intcount.eq.2) then
                     i4_2(ja:jb:jc) = nint(values(iarray))
                  elseif(n_intcount.eq.3) then
                     i4_3(ja:jb:jc) = nint(values(iarray))
                  elseif(n_intcount.eq.4) then
                     i4_4(ja:jb:jc) = nint(values(iarray))
                  elseif(n_intcount.eq.5) then
                     i4_5(ja:jb:jc) = nint(values(iarray))
                  end if
               end if
            end do
         end if
         goto 1000
 2000 continue

c     Error checking for irregular input (missing nodes, etc.)
      nfound = 0
      nnotfound = 0
      do inode = 1, npoints
         if( ifind(inode) .eq. 0 ) then
            nnotfound = nnotfound + 1
            notset(nnotfound) = inode
         else
            nfound = nfound + 1
            isset(nfound) = inode
         end if
      end do
      if (out_number .ne. 0) then
         write(out_number,*)
         write(out_number,*) 'Macro: ', macro, '  Group:', igroup
         if( nfound .eq. npoints) then
            write(out_number,*) 'Values were explicitly set for each ',
     &           'node'
            write(out_number,*)
         elseif(nnotfound .eq. npoints) then
            write(out_number,*) 'Values were set to defaults for each ',
     &           'node'
            write(out_number,*) 'or were set in a previous invokation ',
     &           'of this macro'
            write(out_number,*)
         elseif(nfound .ge. nnotfound) then
            write(out_number,*) 'Data not explicitly set for the ',
     &           'following nodes:'
            write(out_number,9004) (notset(inode), inode = 1, nnotfound)
            write(out_number,*) 'Default values used for these nodes'
            write(out_number,*) 
            write(out_number,*) 'or were set in a previous invokation ',
     &           'of this macro'
         else
            write(out_number,*) 'Data explicitly set only for the ',
     &           'following nodes:'
            write(out_number,9004) (isset(inode), inode = 1, nfound)
            write(out_number,*) 'The remainder are set to the default ',
     &           'value'
            write(out_number,*) 'or were set in a previous invokation ',
     &           'of this macro'
            write(out_number,*)
         end if
      end if
 9000 continue
      deallocate(ifind,isset,notset)
      if(ireturn .eq. -1 .or. ireturn .eq. -2) then
         write (ierr, 9001) iarray, macro, igroup
         if (iptty .ne. 0) write (iptty, 9001) iarray, macro, igroup
         if (iout .ne. 0) write (iout, 9001) iarray, macro, igroup
         if (ireturn .eq. -1) then
            write (ierr, 9002)
            if (iptty .ne. 0) write (iptty, 9002)
            if (iout .ne. 0) write (iout, 9002)
         else if ( ireturn .eq. -2 ) then
            write (ierr, 9003) inumber
            if (iptty .ne. 0) write (iptty, 9003) inumber
            if (iout .ne. 0) write (iout, 9003) inumber
         end if
         stop
      end if
 9001 format (/, 'Fatal error - for array number ', i8,
     .     /, 'macro - ', a4, /, 'Group number - ', i8)
 9002 format ('Something other than a real or integer has been',
     .     ' specified')
 9003 format ('Line number - ', i8, /, 'Bad input, check this line')
 9004 format (1x,10i8) 
      end
         subroutine zone_char_fill(iflg,inode,nzone_used)
c 
c gaz 101323 intial coding (tract zones used be a given node)
c
        use comai
        use combi
        use comdti, only : n0
        implicit none
        integer  jb, jc, iflg, nzone_used, inode
        character*5 dumzone

c allocate memory if needed
        if(.not.allocated(zones_char)) allocate(zones_char(n0))
        if(iflg.eq.1) then 
c gaz 062723  write nodes to zone_char  for no saved zone
          write(dumzone(1:5),'(i5)') nzone_used    
           jb = inode      
           do jc = 1,25           
            if(zones_char(jb)(jc:jc+4).eq.dumzone(1:5)) then
             go to 8001
            else if(zones_char(jb)(jc:jc+4).eq.'    ') then
             write( zones_char(jb)(jc:jc+4),'(i5)') nzone_used
             go to 8001
            endif
           enddo
8001      continue                           
       
        else if(iflg.eq.2) then
c gaz 101523  write nodes to zone_char  for ja,jb,jc =' 1 0 0'
          write(dumzone(1:5),'(a5)') ' all ' 
c gaz 101523  inode = 1, nzone_use = n0
          do jb  = inode, nzone_used
           do jc = 1,25           
            if(zones_char(jb)(jc:jc+4).eq.dumzone(1:5)) then
             go to 8002
            else if(zones_char(jb)(jc:jc+4).eq.'    ') then
             write( zones_char(jb)(jc:jc+4),'(a5)') ' all '
             go to 8002
            endif
           enddo
8002      continue           
          enddo                          

        endif
        return
        end
      subroutine zone_saved(iflg,macr,nzone_saved,icall_sv,input_unit)
c subroutine created gaz 060523
c gaz 060617 algoithm created        
c manage saved zones 
      use combi
      use comdti, only : n0
      use comai
      use comdi, only : izone_conv_nodes
      implicit none
      integer i,ii,izunit,nin,j
      integer iflg,nzone_saved,icall_sv
c gaz 061223
      integer input_unit, maxlines
      integer open_file
      character*30 zonesavename
c gaz 110123
      character*4 macr
      logical ex,op
      parameter (maxlines = 10000)
        if(iflg.eq.0) then
c initialize if necessary   
        else if(iflg.eq.1) then 
c identify and count saved zones
c gaz 061223                   
         icall_sv = 0
         do ii = 1, maxlines
         read(input_unit,'(a80)') wdd(1:80)
         if(wdd(1:4).eq.'stop') go to 600
120      if(wdd(1:4).eq.'zone') then
         read(input_unit,'(a80)') wdd(1:80)
121      do i = 1, 76
           if(wdd(i:i+3).eq.'save') then
            icall_sv = icall_sv +1
            read(wdd(1:i-1),*) izonesavenum(icall_sv)
            go to 122
           endif
          enddo
122       read(input_unit,'(a80)') wdd(1:80)
          if(wdd(1:4).eq.'zone') then
           go to 120
          else if(wdd(1:4).eq.'stop') then
           go to 600
          else
           go to 121
          endif
         endif
         enddo
600      num_sv_zones = icall_sv
         return
        else if(iflg.eq.2) then   
          zonesavename(1:14) = 'zone00000.save'
          write(zonesavename(5:9),'(i5)') abs(nzone_saved)+10000
          zonesavename(5:5) = '0'
          ex = .false.
          op = .false.
          inquire (file = zonesavename, exist = ex)
c gaz 061523 
          if(ex) then
           inquire (file = zonesavename, opened = op)   
           if(.not.op) then
            izunit=open_file(zonesavename,'unknown')
           else
            inquire(file = zonesavename, number=izunit)
            close(izunit)
            izunit=open_file(zonesavename,'unknown')
           endif
           read(izunit,*)  nzone_saved
           read(izunit,*)  wdd
           read(izunit,*) nin
           if(allocated(ncord)) deallocate(ncord)
           allocate(ncord(nin))
           backspace izunit
           read(izunit,*) nin, (ncord(i), i =1, nin)
           close(izunit)
c gaz 090523 don't need conv info here
c send info to innode to read zones in ncord
           input_unit = nin
           if(macr.eq.'conv') then
            do i = 1, nin
             ii = ncord(i)
             izone_conv_nodes(ii+n0*(icall_sv-1)) = nzone_saved
            enddo
           else if(macr.eq.'node') then
c gaz node is singe call
c zone info used from ncord
            continue
           else
           endif
          else
c gaz 101023 input_unit = 0 means not a saved zone 
           input_unit = 0
          endif
          else if(iflg.eq.3) then
c gaz 062723  close files
           do i = 1, num_sv_zones 
              ii = izonesavenum(i)
              zonesavename(1:14) = 'zone00000.save'
              write(zonesavename(5:9),'(i5)') abs(ii)+10000
              zonesavename(5:5) = '0'
              ex = .false.
              op = .false.
              inquire (file = zonesavename, exist = ex)
              if(ex) then
               inquire (file = zonesavename, opened = op)
               if(op) then
                inquire(file = zonesavename, number = izunit)
                close (izunit,status='delete')
               else
                izunit=open_file(zonesavename,'old')
                close (izunit,status='delete')
               endif
              endif
           enddo
          continue
         else
          nzone_saved = 0
         endif
         return
         end


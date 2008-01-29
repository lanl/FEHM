      subroutine initfc
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Read interface reduction factor data. 
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: ?, Programmer: ?
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/initfc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.12 Mass transport at interfaces
!D3 2.6    Provide Input/Output Data Files
!D3 3.0    INPUT AND OUTPUT REQUIREMENTS
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comai
      use comdti
      use combi
      use comki
      use comdi
      implicit none

      character*80 dummy_line, filename
      logical null1, opened
      integer jj, kk, ll, isizes, ilines
      integer inptread, open_file

      macro = 'itfc'
      
c     Set current zone info in array for later use in figuring
c     interface factors for connections

      izonef_itfc = izonef
      opened = .false.



c**** read interface reduction factor data ****

      jj = 0
      kk = 0
      isizes = 0
      red_factor(0) = 1.
      ftn_factor(0) = 1.
 6000       continue
            read(inpt,'(a80)') dummy_line
            if(null1(dummy_line)) then
               goto 6001
            else
               backspace inpt
               jj = jj + 1
               read(inpt,*) zone_pair(jj,1),zone_pair(jj,2),
     2              red_factor(jj)
               red_factor(jj) = max(1.d-30,red_factor(jj))
            end if
            goto 6000
6001       continue
            read(inpt,'(a80)') dummy_line
            if(null1(dummy_line)) then
               goto 6002
            else
               backspace inpt
               kk = kk + 1
               if(kk.eq.1) read(inpt,*) (filter_flag(ll),ll=1,nspeci)
               read(inpt,*) zonec_pair(kk,1),zonec_pair(kk,2),
     2              ftn_factor(kk)
               if(ftn_factor(kk).lt.0) then
c     Option to read from file
                  read(inpt,'(a80)') dummy_line
                  if(dummy_line(1:4).eq.'file') then
                     read(inpt,'(a80)') filename
                     inptread = open_file(filename,'old')
                     opened = .true.
                  else
                     backspace inpt
                     inptread = inpt
                     opened = .false.
                  end if
                  isizes = isizes + 1
                  itfcsize(kk) = isizes
                  ilines = 0
7001              continue
                  read(inptread,'(a80)') dummy_line
                  if(null1(dummy_line)) then
                    goto 7002
                  else
                    ilines = ilines + 1
                    backspace inptread
                    read(inptread,*) itfcporsize(ilines,isizes),
     2                           itfcprobsize(ilines,isizes)
                    goto 7001
                  end if
               end if
7002           continue
               if(opened) close (inptread)
            end if
            goto 6001
6002       continue
      

      
      end

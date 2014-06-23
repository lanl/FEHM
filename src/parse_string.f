      subroutine parse_string(line,imsg,msg,xmsg,cmsg,nwds)
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
CD1 Parse character string into tokens.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/parse_string.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:34   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
C**********************************************************************
CD3
CD3 REQUIREMENTS TRACEABILITY
CD3
CD3 N/A
CD3
C**********************************************************************
CD4
CD4 SPECIAL COMMENTS
CD4
CD4 Character strings containing only d,D,e,E and integers will
CD4 not be accepted.
CD4
C***********************************************************************

      implicit none

      integer max_entries, line_length
      parameter(max_entries=20)
      character(*) line
      integer imsg(max_entries)
      integer msg(max_entries)
      real*8 xmsg(max_entries)
      character*32 cmsg(max_entries)
      integer nwds
      logical finished

      integer ndex(max_entries,2),i,begin,entrynum,isinteger,isreal,i2

      line_length = len(line)
      entrynum=1
      begin=1
      do i=1,line_length
         if (((line(i:i).eq.' ').or.(line(i:i).eq.achar(9)))
     &        .and.(begin.eq.0)) then
            ndex(entrynum,2)=i-1
            begin=1
            entrynum=entrynum+1
         else if ((line(i:i).eq.' ').or.(line(i:i).eq.achar(9))) then
            continue
         else if (begin.eq.1) then
            ndex(entrynum,1)=i
            begin=0
         else
            continue
         endif
      enddo
      if (begin.eq.1) entrynum=entrynum-1
      nwds=entrynum

      do i=1,nwds
         isinteger=1
         isreal=1
         do i2=ndex(i,1),ndex(i,2)
            if (ndex(i,1) .eq. ndex(i,2)) then
c check if just have +, -, e, E, d, or D
               if ( (line(i2:i2).eq.'+') .or. (line(i2:i2).eq.'-') .or.
     &              (line(i2:i2).eq.'e') .or. (line(i2:i2).eq.'E') .or.
     &              (line(i2:i2).eq.'d') .or. (line(i2:i2).eq.'D') ) 
     &              then
                  isreal = 0
                  isinteger = 0
                  exit
               end if
            end if
            if ((line(i2:i2).ne.'+').and.(line(i2:i2).ne.'-').and.
     &           ((iachar(line(i2:i2)).gt.57).or.
     &           (iachar(line(i2:i2)).lt.48))) isinteger=0
            if ((line(i2:i2).ne.'+').and.(line(i2:i2).ne.'-').and.
     &           ((iachar(line(i2:i2)).gt.57).or.
     &           (iachar(line(i2:i2)).lt.48)).and.
     &           (line(i2:i2).ne.'.').and.
     &           (line(i2:i2).ne.'E').and.
     &           (line(i2:i2).ne.'e').and.
     &           (line(i2:i2).ne.'D').and.
     &           (line(i2:i2).ne.'d')) isreal=0
         enddo
         imsg(i) = 0
         xmsg(i) = 0.
         cmsg(i) = ''
         if ((isreal.eq.1).and.(isinteger.eq.0)) then
            msg(i)=2
            read(line(ndex(i,1):ndex(i,2)),*) xmsg(i)
         else if (isinteger.eq.1) then
            msg(i)=1
            read(line(ndex(i,1):ndex(i,2)),*) imsg(i)
         else
            msg(i)=3
            cmsg(i)=line(ndex(i,1):ndex(i,2))
         endif
      enddo

      end

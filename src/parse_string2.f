      subroutine parse_string2(line,imsg,msg,xmsg,cmsg,nwds)
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
      real*8 rnum
      character*32 cmsg(max_entries), word
      integer nwds, ifdebug, inum
      logical finished, isint, isrl

      integer ndex(max_entries,2),i,begin,entrynum

      imsg = 0
      xmsg = 0.
      cmsg = ''
      msg = 0
      nwds = 0
	
      ifdebug = 0
      line_length = len(line)
      entrynum=1
      begin=1
      do i=1,line_length
         if (((line(i:i).eq.' ').or.(line(i:i).eq.achar(9)))
     &        .and.(begin.eq.0)) then
            ndex(entrynum,2)=i-1
            begin=1
            entrynum=entrynum+1
            if (entrynum .gt. max_entries) then
                entrynum = max_entries
                exit
            end if
         else if ((line(i:i).eq.' ').or.(line(i:i).eq.achar(9))) then
            continue
         else if (begin.eq.1) then
            ndex(entrynum,1)=i
            begin=0
         else if (i .eq. line_length) then
            ndex(entrynum,2)=i
         else
            continue
         endif
      enddo
      if (begin.eq.1) entrynum=entrynum-1
      nwds=entrynum

      do i=1,nwds
         isint = .false.
         isrl = .false.
         word = line(ndex(i,1):ndex(i,2))

         if (ifdebug .eq. 1) then
            write (6, *) 'ndex ', ndex(i,1), ndex(i,2)
            write (6, *) 'string ', line(ndex(i,1):ndex(i,2))
         end if

         call isinteger (word, inum, isint)

         if (.not. isint) then
            imsg(i) = 0
            call isreal (word, rnum, isrl)
            if (.not. isrl) then
               cmsg(i) = word
               msg(i) = 3
            else
               msg(i) = 2
               xmsg(i) = rnum
            end if
         else
            msg(i) = 1
            imsg(i) = inum
         end if
         
      enddo

      end

      subroutine isinteger (word, int, isintgr)

      implicit none
      logical isintgr
      integer len,int,i,j,iflag,m,itop,ibot
      character*(*) word

C     integer is +- numbers
      len = len_trim(word)
      iflag=1
      int=0
      j=1
      if(word(1:1).eq.'-') then
         iflag=-1
         j=2
      elseif(word(1:1).eq.'+') then
         j=2
      endif
      do i=j,len
         m=ichar(word(i:i))
         itop=ichar('9')
         ibot=ichar('0')
         if(m.ge.ibot.and.m.le.itop) then
            int=int*10 + m-ibot
         else
            isintgr=.false.
            return
         endif
      enddo
      isintgr=.true.
      int=iflag*int

      end subroutine isinteger

      subroutine isreal (word, xnum, isrl)

      implicit none
      integer len,int,j,iflag,ii,m,ibot,itop
      character*(*) word
      logical isrl
      real*8 xnum,power,ten

      len = len_trim(word)
      isrl=.false.
      ten=10.0d0
      xnum=0.0d0
      iflag=1
      j=1
      if(word(1:1).eq.'-') then
         iflag=-1
         j=j+1
      elseif (word(1:1).eq.'+') then
         j=j+1
      endif
C     look for base
      int=0
      ii=j
      itop=ichar('9')
      ibot=ichar('0')
      m=ichar(word(ii:ii))
      if ((m.gt.itop.or.m.lt.ibot).and.word(ii:ii).ne.'.') then
         isrl=.false.
         return
      endif
      do while (m.ge.ibot.and.m.le.itop)
         int=int*10+m-ichar('0')
         ii=ii+1
         m=ichar(word(ii:ii))
      enddo
      if(word(ii:ii).eq.'.') then
C     get fraction part of base
         ii=ii+1
         if(ii.gt.len) then
            isrl=.true.
            xnum=int*iflag
            go to 9999
         endif
         xnum=int
         power=.1d0
         m=ichar(word(ii:ii))
         do while (m.ge.ibot.and.m.le.itop)
            xnum=xnum+(m-ibot)*power
            power=power*.1d0
            ii=ii+1
            m=ichar(word(ii:ii))
         enddo
         xnum=xnum*iflag
         if(ii.gt.len) then
            isrl=.true.
            go to 9999
         endif
      endif

      if(word(ii:ii).eq.'e'.or.word(ii:ii).eq.'E'.or.
     *     word(ii:ii).eq.'+'.or.word(ii:ii).eq.'-'.or.
     *     word(ii:ii).eq.'d'.or.word(ii:ii).eq.'D') then
         if(ii.eq.len) go to 9999
         if(xnum.eq.0.0d0.and.int.eq.0.and.ii.eq.1) go to 9999
         if(xnum.eq.0d0) xnum=int*iflag
         iflag=1
         if(word(ii:ii).eq.'-') iflag=-1
C     get exponent
         ii=ii+1
         int=0
         if(word(ii:ii).eq.'+') ii=ii+1
         if(word(ii:ii).eq.'-') then
            ii=ii+1
            iflag=-1
         endif
         m=ichar(word(ii:ii))
         do while
     *        (m.ge.ibot.and.m.le.itop.and.
     *        ii.le.len)
            int=int*10+m-ichar('0')
            ii=ii+1
            m=ichar(word(ii:ii))
         enddo
         if(ii.le.len) go to 9999
         if(iflag.eq.1) then
            xnum=xnum*ten**int
         else
            xnum=xnum/ten**int
         endif
         isrl=.true.
      endif
 9999 continue
      return
      end subroutine isreal

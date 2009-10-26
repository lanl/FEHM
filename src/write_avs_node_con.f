      subroutine write_avs_node_con(icall,an,anv,npt,neq,nspeci,
     .     lu, ifdual)
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
CD1 Output concentration fields from FEHM
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-SEP-93    Carl Gable     22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/write_avs_node_con.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Fri Feb 02 14:19:38 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.8   09/11/95 17:30:46   awolf
CD2 Set so all concentration values less than 1e-20 are set to 1e-20
CD2 for use with avs
CD2 
CD2    Rev 1.7   08/07/95 11:49:36   awolf
CD2 Fixed for avs output of dpdp
CD2 
CD2    Rev 1.6   04/10/95 15:16:48   llt
CD2 changed max function for IBM compatibility
CD2 
CD2    Rev 1.5   04/10/95 14:44:04   robinson
CD2 Minimum output concentration is now 0 rather than 1.e-20
CD2 
CD2    Rev 1.4   01/28/95 16:14:06   llt
CD2 needed endif
CD2 
CD2    Rev 1.3   01/28/95 14:21:10   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.2   01/20/95 13:29:12   tam
CD2 Changed format for strings from * to a56, kept length to 80 so left justified
CD2 
CD2    Rev 1.1   11/29/94 13:42:32   llt
CD2 Removed protying, so could run under cc compile, instead of c compiler. 
CD2 Used #ifdef to determine format for routine name, depending on machine.
CD2 (Changes made by tam.)
CD2 
CD2    Rev 1.0   08/23/94 15:33:54   llt
CD2 Original version
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
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
CD4 None
CD4
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
CPS BEGIN 
CPS   
CPS END 
CPS 
C***********************************************************************

      use avsio
      use comai, only : altc, days, icnl, jdate, jtime, nei_in,
     &     ns_in, verno, wdd, neq_primary, ivf, ifdm_elem
      use combi, only : corz, izonef, nelm
      use comchem
      use comdi, only : nsurf, izone_surf, izone_surf_nodes, icns
      use comrxni
      use comdti
      implicit none

      integer add_dual, maxcon, iz, idz, iendz, il, open_file
      integer neq,nspeci,lu,ifdual,icall,length,i1,i2
      integer icord1, icord2, icord3, iaq, ivap, isolid
      integer npt(*)
      integer, allocatable ::  nelm2(:)
      parameter (maxcon = 100)
      real*8 an(n0,nspeci)
      real*8 anv(n0,nspeci)
      character*60, allocatable :: title(:)
      character*14 tailstring
      character*8 dual_char
      character*3 dls
      character*5 char_type
      character*60 fstring
      character*30 cordname(3)
      character*150 :: string = '', tecstring = '', sharestring = ''
      real*8 write_array(maxcon)
      integer i, ic, im, in, iv, ix, istep, j, n
      integer t1(maxcon),itotal2,write_total
      integer irxn_title
      real*8 complex_conc

      save tecstring, sharestring

      data cordname(1) / 'X (m)' /, 
     &     cordname(2) / 'Y (m)' /,
     &     cordname(3) / 'Z (m)' /  

      if(ifdual .eq. 0)then
         istep = 0
         add_dual=0
         dual_char = ''
         tailstring = '_con_node'
      else
         istep = maxcon
         add_dual=neq
         dual_char = 'Dual '
         tailstring = '_con_dual_node'
      endif

      if (iocord .ne. 0) then
         select case (icnl)
         case (1, 4)
            icord1 = 1
            icord2 = 2
            icord3 = 1
         case (2, 5)
            icord1 = 1
            icord2 = 3
            icord3 = 2
         case(3, 6)
            icord1 = 1
            icord2 = 3
            icord3 = 1
         case default
            icord1 = 1
            icord2 = 3
            icord3 = 1
         end select
      end if

      if (iozone .ne. 0 ) then
         iendz = nsurf
      else 
         iendz = 1
         idz = iozone
      end if

      if (altc(1:4) .eq. 'avsx') then
         dls = ' : '
      else if (altc(1:3) .eq. 'sur') then
         dls = ', '
      else
         dls = ' '
      end if
      
      if (altc(1:3) .ne. 'sur') then
         call namefile2(icall,lu,ioformat,tailstring,0)
! file will be opened in zone loop for surfer
      end if

      if (altc(1:3) .eq. 'tec') then
         if (icall .eq. 1 .and. iogeo .eq. 1) then
            iz = 1
            tecstring = ''
            select case (ns_in)
            case (5,6,8)
               write (tecstring, 135) neq, nei_in, 'FEBRICK'
            case (4)
               if (icnl .eq. 0) then
                  write (tecstring, 135) neq, nei_in, 
     &                 'FETETRAHEDRON'
               else
                  write (tecstring, 135) neq, nei_in, 
     &                 'FEQUADRILATERAL'
               end if
            case (3)
               write (tecstring, 135) neq, nei_in, 'FETRIANGLE'
            case (2)
               write (tecstring, 135) neq, nei_in, 'FELINESEG'
            case (0)
! fdm grid
               write (tecstring, '(a)') ''
               write (string, '(a)') ''
            end select
            if (icnl .eq. 0 .and. ns_in .ne. 0) then
               if (iozid .eq. 0) then
                  write (string, 140) '1-3', iz
               else
                  write (string, 140) '1-3, 5', iz
               end if
            else if (icnl .ge. 1 .and. ns_in .ne. 0) then
               if (iozid .eq. 0) then
                  write (string, 140) '1-2', iz
               else
                  write (string, 140) '1-2, 4', iz
               end if            
            end if
         else if (icall .eq. 1 .and. iocord .ne. 0) then
            if (icnl .eq. 0) then
               if (iozid .eq. 0) then
                  string = '1-3'
               else
                  string = '1-3, 5'
               end if
            else
               if (iozid .eq. 0) then
                  string = '1-2'
               else
                  string = '1-2, 4'
               end if
            end if
            if (iozone .ne. 0) sharestring = string
         else if (icall .eq. 1 .and. iozid .ne. 0) then
            string = '2'
            if (iozone .ne. 0) sharestring = string
         end if
      end if

      if(rxn_flag.eq.0)then
         if(nspeci .gt. maxcon)then
            write(lu,*)'--------------------------------------------'
            write(lu,*)'ERROR: WRITE_AVS_NODE_CON'
            write(lu,*)'nspeci = ',nspeci,
     2           ' is greater than maxcon = ',maxcon
            write(lu,*)'Subroutine only able to handle up to', maxcon,
     &           'tracers'
            write(lu,*)'--------------------------------------------'
            return
         endif

         
c     error
c     

c     Formatted write to accomodate the different way fortran 90
c     does unformatted writes (f77 way no longer worked because
c     it put commas between the numbers, avs didn't like that)

c----------------------------------------
c     PHS 4/27/00   Added new format for avsxpress!!
c---------------------------------------
         fstring = ''
         if (altc(1:3) .eq. 'tec') then
            itotal2 = nspeci+iocord+iozid+1
         else
            itotal2 = nspeci+iocord+iozid
         end if
         allocate (title(itotal2))
         do i = 1, iocord
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(i) = trim(cordname(i)) // ', (m)'
            else
               title(i) = trim(cordname(i))
            end if
         end do
         if (altc(1:3) .eq. 'tec') then
            title(i) = 'Node'
            i = i + 1
         end if
         if (iozid .eq. 1) then
            title(i) = 'Zone'
            i = i + 1
         end if
         j = i
         iaq = 0
         ivap = 0
         isolid = 0
         do i = 1, nspeci
            select case (icns(i))
            case (1, 2, -2)
               iaq = iaq + 1
               title(j) = trim(dual_char)//trim(cpntnam(iaq))
            case (0)
               isolid = isolid + 1
               title(j) = trim(dual_char)//trim(immnam(isolid))
            case (-1)
               ivap = ivap + 1
               title(j) = trim(dual_char)//trim(vapnam(ivap))
            end select
            j = j+ 1
         end do
      
         if(altc(1:4).eq.'avsx') then
            write(fstring, 664) itotal2
            write(lu,fstring) days,(trim(title(i)),i=1,itotal2)
         else if (altc(1:3) .eq. 'avs') then
            write(lu,'(i3, 100(1x, i1))') itotal2,(1,i=1,nspeci+iocord)
            write(lu, '(a)') (trim(title(i)),i=1,itotal2)
         else if (altc(1:3) .eq. 'tec') then
            if (icall .eq. 1) then
               write(lu, 98) verno, jdate, jtime, trim(wdd)
               write (fstring, 99) itotal2
               write(lu, fstring) 'VARIABLES = ', (trim(title(i)), 
     &              i=1,itotal2)
            end if
            if (iozone .ne. 0) write (lu, 97) trim(timec_string)
         end if

         do iz = 1, iendz
c     Zone loop
            if (iozone .ne. 0) then
               idz = izone_surf(iz)
            end if
            if (altc(1:3) .eq. 'tec') then
               if (iozone .eq. 0) then
                  write (lu, 94) trim(timec_string), trim(tecstring)
               else
                  if (icall .gt. 1 .and. iozone .ne. 0) then
                     write (tecstring, 125) trim(sharestring), iz
                  end if
                  write (lu, 95) idz, trim(tecstring)
               end if
               if (icall .eq. 1 .and. iz .eq. iendz) then
                  if (iogeo .eq. 1) then
                     tecstring = trim(tecstring) // trim(string)
                  else if (iocord .ne. 0 .or. iozid .ne. 0) then
                     write (tecstring, 125) trim(string), iz
                  end if   
               end if
            else if (altc(1:3) .eq. 'sur') then
               call namefile2(icall,lu,ioformat,tailstring,idz)
               fstring = ''
               write (fstring, 96) itotal2
               write(lu, fstring) (trim(title(i)), i=1,itotal2)
            end if
            fstring = ''
            if (iozid .eq. 0 .and. iocord .eq. 0) then
               write (fstring, 333) itotal2, dls
            else if (iozid .eq. 0 .and. iocord .ne. 0) then
               if (altc(1:3) .ne. 'tec') then
                  write (fstring, 333) itotal2, dls
               else
                  if (icall .eq. 1) then
                     write (fstring, 104) iocord, nspeci
                  else
                     write (fstring, 333) itotal2 - iocord, dls
                  end if
               end if
            else if (iozid .eq. 1 .and. iocord .eq. 0) then
               if (altc(1:3) .ne. 'tec' .or. icall .eq. 1) then
                  write (fstring, 334) dls, nspeci, dls
               else
                  write (fstring, 333) itotal2, dls
               end if
            else if (iozid .eq. 1 .and. iocord .ne. 0) then
               if (altc(1:3) .ne. 'tec') then
                  write (fstring, 335) iocord, dls, nspeci, dls,
     &                 nspeci, dls
               else
                  if (icall .eq. 1) then
                     write (fstring, 105) iocord, nspeci
                  else
                     write (fstring, 333) itotal2, dls
                  end if
               end if
            end if               
            
            do i = 1, neq
               if (iozone .ne. 0) then
                  if (izone_surf_nodes(i).ne.idz) goto 199
               end if
               if (altc(1:3) .eq. 'tec' .and. iocord .ne. 0) then
                  if (icall .eq. 1 .and. iozid .eq. 0) then
                     write(lu, fstring) (corz(i,j), j = icord1, icord2,
     +                    icord3), i, (min(1.0d+20, max(1.0d-20,
     +                    an(i+add_dual,n))), n=1,nspeci)
                  else if (icall .eq. 1 .and. iozid .eq. 1) then
                     write(lu, fstring) (corz(i,j), j = icord1, icord2,
     +                    icord3), i,izonef(i), (min(1.0d+20, 
     +                    max(1.0d-20,an(i+add_dual,n))), n=1,nspeci)
                  else
                     write(lu, fstring) i, (min(1.0d+20, max(1.0d-20,
     +                    an(i+add_dual,n))), n=1,nspeci)
                  end if

               else if (altc(1:3) .eq. 'tec' .and. iozid .ne. 0) then
                  if (icall .eq. 1) then
                     write(lu, fstring) i, izonef(i), (min(1.0d+20,
     +                    max(1.0d-20, an(i+add_dual,n))), n=1,nspeci)
                  else
                     write(lu, fstring) i, (min(1.0d+20, max(1.0d-20,
     +                    an(i+add_dual,n))), n=1,nspeci)
                  end if
                  
               else if (iocord .ne. 0) then
                  if (iozid .eq. 0) then
                     write(lu, fstring) i, (corz(i,j), j = icord1, 
     +                    icord2, icord3), (min(1.0d+20, max(1.0d-20,
     +                    an(i+add_dual,n))), n=1,nspeci)
                  else
                     write(lu, fstring) i, (corz(i,j), j = icord1,
     +                    icord2, icord3), izonef(i), (min(1.0d+20, 
     +                    max(1.0d-20, an(i+add_dual,n))), n=1,nspeci)
                  end if

               else
                  if (iozid .eq. 0) then
                     write(lu, fstring) i,(min(1.0d+20, max(1.0d-20,
     +                    an(i+add_dual,n))), n=1,nspeci)
                  else
                     write(lu, fstring) i, izonef(i), (min(1.0d+20,
     +                    max(1.0d-20, an(i+add_dual,n))), n=1,nspeci)
                  end if
               end if
 199        enddo
            if (altc(1:3) .eq. 'sur') close (lu)
         end do
         deallocate (title)

c---------------------------------------------------
c     PHS 4/28/2000    end of changes for xpress format
c     PHS 8/10/2000    more changes for xpress below. . .
c     for when rxn_flag NE zero
c---------------------------------------------------
      else                      ! IF rxn_flag  NE zero

         fstring = ''
         if(ncplx.ne.0)then
            itotal2 = ncpntprt*2+nimmprt+ncplxprt+nvapprt+iocord+iozid
         else
            itotal2 = ncpntprt+nimmprt+nvapprt+iocord+iozid
         endif

         allocate (title(itotal2))
         do i = 1,itotal2
            t1(i)=1
         enddo
c=======================================================
c     Titles for the rxn stuff  (like Np[aq] ,(Moles/kg H20))
c=======================================================
         do i = 1, iocord
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(i) = trim(cordname(i)) // ', (m)'
            else
               title(i) = trim(cordname(i))
            end if
         end do

         if (iozid .eq. 1) then
            title(i) = 'Zone'
         end if


         irxn_title = iocord + iozid         
         if(altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x')
     2        write(lu,'(i3, 100(1x, i1))')itotal2,(t1(i),i=1,itotal2)


         do i = 1,ncpntprt
            ic = cpntprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(cpntnam(ic))
               write(lu,190) trim(title(irxn_title))
            else
               title(irxn_title)= trim(dual_char)
     &              //trim(cpntnam(ic))//' (Moles/kg H20)'
            end if
         enddo

         do i = 1,nimmprt
            im = immprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(immnam(im))
               write(lu,200) trim(title(irxn_title))
            else
               title(irxn_title)= trim(dual_char)
     &              //trim(immnam(im))//' (Moles/kg Rock)'
            end if
         enddo

         do i = 1,nvapprt
            iv = vapprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(vapnam(iv))
               write(lu,210) trim(title(irxn_title))
            else
               title(irxn_title)= trim(dual_char)
     &              //trim(vapnam(iv))//' (Moles/kg Vapor)'
            end if
         enddo

         if(ncplx.ne.0)then
            do i = 1,ncpntprt
               irxn_title = irxn_title + 1
               ic = cpntprt(i)
               if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
                  title(irxn_title) = trim(dual_char)//'Free Ion '
     &                 //trim(cpntnam(ic))
                  write(lu,220) trim(title(irxn_title))
               else
                  title(irxn_title)=trim(dual_char) //'Free Ion '
     &                 //trim(cpntnam(ic))//' (Moles/kg H20)'
               end if
            enddo
         endif

         do i = 101, ncplxprt+100
            ix = cplxprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(cplxnam(ix))
               write(lu,190) trim(title(irxn_title))
            else
               title(irxn_title)=trim(dual_char)
     &              //trim(cplxnam(ix))//' (Moles/kg H20)'
            end if
         enddo
         
         fstring = ''
         if(altc(1:4).eq.'avsx') then
            write(fstring, 664) itotal2
            write(lu, fstring) days,(trim(title(i)),i=1,itotal2)

         else if (altc(1:3) .eq. 'tec') then
            if (icall .eq. 1) then
               write (fstring, 99) itotal2+1
               write(lu, 98) verno, jdate, jtime, trim(wdd)
               write(lu, fstring) 'VARIABLES = ', (trim(title(i)), 
     &              i=1,iocord), 'node', (trim(title(i)), 
     &              i=iocord+1,itotal2)
            end if
            if (iozone .ne. 0) write (lu, 97) trim(timec_string)
         end if

c=================================================
c     Write out node info to the  _con file
c=================================================

         do iz = 1, iendz
! Zone loop
            if (iozone .ne. 0) then
               idz = izone_surf(iz)
            end if
            if (altc(1:3) .eq. 'tec') then
               if (iozone .eq. 0) then
                  write (lu, 94) trim(timec_string), trim(tecstring)
               else
                  if (icall .gt. 1 .and. iozone .ne. 0) then
                     write (tecstring, 125) trim(sharestring), iz
                  end if
                  write (lu, 95) idz, trim(tecstring)
               end if
               if (icall .eq. 1 .and. iz .eq. iendz) then
                  if (iogeo .eq. 1) then
                     tecstring = trim(tecstring) // trim(string)
                  else if (iocord .ne. 0 .or. iozid .ne. 0) then
                     write (tecstring, 125) trim(string), iz
                  end if   
               end if
            else if (altc(1:3) .eq. 'sur') then
               call namefile2(icall,lu,ioformat,tailstring,idz)
               write (fstring, 96) itotal2
               write(lu, fstring) (trim(title(i)), i=1,itotal2)
            end if
            do in = 1,neq
               if (iozone .ne. 0) then
                  if (izone_surf_nodes(in).ne.idz) goto 299
               end if               
               j=0
               if (iocord .ne. 0) then
                  do i = icord1, icord2, icord3
                     j = j+1
                     write_array(j)=corz(in,i)
                  end do
               end if
               do i = 1,ncpntprt
                  ic = cpntprt(i)
                  j=j+1
                  write_array(j)=an(in+add_dual,pcpnt(ic))
               enddo
               do i = 1,nimmprt
                  im = immprt(i)
                  j=j+1
                  write_array(j)=min(1.0d+40,
     2                 max(1.0d-90,an(in+add_dual,pimm(im))))
               enddo
               do i = 1,nvapprt
                  iv = vapprt(i)
                  j=j+1            
                  write_array(j)=min(1.0d+40,
     2                 max(1.0d-90,an(in+add_dual,pvap(iv))))
               enddo
               if(ncplx.ne.0)then
                  do i = 1,ncpntprt
                     ic = cpntprt(i)
                     j=j+1
                     write_array(j)=min(1.0d+40,
     2                    max(1.0d-90,cpntsv(ic,in+add_dual)))
                  enddo
               endif
               do i = 101,ncplxprt+100
                  ix = cplxprt(i)
                  j=j+1
                  call cplxcalc(in+add_dual,ix,complex_conc)
                  write_array(j)=min(1.0d+40,
     2                 max(1.0d-90,complex_conc))
               enddo
               write_total = j
               fstring = ''
               if (iozid .eq. 0 .and. iocord .eq. 0) then
                  write (fstring, 333) write_total, dls
               else if (iozid .eq. 0 .and. iocord .ne. 0) then
                  if (altc(1:3) .ne. 'tec') then
                     write (fstring, 333) write_total, dls
                  else
                     if (icall .eq. 1) then
                        write (fstring, 104) iocord, write_total-iocord
                     else
                        write (fstring, 333) write_total - iocord, dls
                     end if
                  end if
               else if (iozid .eq. 1 .and. iocord .eq. 0) then
                  if (altc(1:3) .ne. 'tec' .or. icall .eq. 1) then
                     write (fstring, 334) dls, write_total, dls
                  else
                     write (fstring, 333) write_total, dls
                  end if
               else if (iozid .eq. 1 .and. iocord .ne. 0) then
                  if (altc(1:3) .ne. 'tec') then
                     write (fstring, 335) iocord, dls, dls, 
     &                    write_total - iocord, dls
                  else
                     if (icall .eq. 1) then
                        write (fstring, 105) iocord, write_total-iocord
                     else
                        write (fstring, 333) write_total - iocord, dls
                     end if
                  end if
               end if
       
               if (altc(1:3) .eq. 'tec' .and. iocord .ne. 0) then
                  if (icall .eq. 1 .and. iozid .eq. 0) then
                     write(lu,fstring) (write_array(i),i=1,iocord), 
     &                    in, (write_array(i),i=iocord+1,write_total)
                  else if (icall .eq. 1 .and. iozid .eq. 1) then
                     write(lu,fstring) (write_array(i),i=1,iocord), in,
     &                    izonef(in), (write_array(i),
     &                    i=iocord+1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i),
     &                    i=iocord+1,write_total)
                  end if

               else if (altc(1:3) .eq. 'tec' .and. iozid .ne. 0) then
                  if (icall .eq. 1) then
                     write(lu,fstring) in, izonef(in), 
     &                    (write_array(i),i=1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i),
     &                    i=iocord+1,write_total)
                  end if

               else if (iocord .ne. 0) then
                  if (iozid .eq. 0) then
                     write(lu,fstring) in, (write_array(i),
     &                    i=1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i), 
     &                    i=1,iocord), izonef(in), 
     &                    (write_array(i),i=iocord+1,write_total)
                  end if
               else
                  if (iozid .ne. 0) then
                     write(lu,fstring) in, izonef(in), 
     &                    (write_array(i),i=1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i), 
     &                    i=1,write_total)
                  end if
               end if
            enddo
            if(altc(1:3) .eq. 'sur') close (lu)
 299     end do
         deallocate (title)
      endif

      if (icall .eq. 1 .and. altc(1:3) .eq. 'tec' .and. iogeo .eq. 1)
     &     then
! Read the element connectivity and write to tec file
         il = open_file(geoname,'old')
! avsx geometry file has an initial line that starts with neq_primary
         allocate(nelm(ns_in))
         read(il,*) i
         if (i .ne. neq_primary) backspace il
         do i = 1, neq
            read(il,*)
         end do
         do i = 1, nei_in
            read (il,*) i1,i2,char_type,(nelm2(j), j=1,ns_in)
            write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
         end do
         deallocate(nelm2)
         close (il)
      end if
c gaz added element output  (hex only) for fdm generated grid   
      if (icall .eq. 1 .and. altc(1:3) .eq. 'tec' .and. ivf .eq. -1
     &     .and. ifdm_elem. eq. 1) then
c first generate elements      
         call structured(4)
         il = open_file('fdm_elem.macro','old')
         read(il,*) 
         read(il,*)  ns_in , nei_in
         allocate (nelm2(ns_in))
         read(il,*) nei_in, ns_in
         do i = 1, nei_in
            read (il,*) i1, (nelm2(j), j=1,ns_in)
            write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
         end do
         deallocate(nelm2)
         close (il)
      end if      
      if (altc(1:3) .ne. 'sur') close (lu)

c 94   format('ZONE T = "Simulation time ',1p,g16.9,' days"', a)
 94   format('ZONE T = ', a, a)
 95   format('ZONE T = "',i4.4,'"', a)
 96   format("('node', ", i3, "(', ', a))")
c 97   format('TEXT T = "Simulation time ',1p,g16.9,' days"')
 97   format('TEXT T = ', a)
 98   format ('TITLE = "', a30, 1x, a11, 1x, a8, 1x, a, '"')
 99   format ('(a11, ', i3, "('",' "',"', a, '",'"',"'))")
 104  format ('(1x, ', i3, '(g16.9, 1x), i10.10, ', i3, '(1x, g16.9))')
 105  format ('(1x, ', i3, '(g16.9, 1x), i10.10, 1x, i4, ', i3, 
     &     '(1x, g16.9))')
 125  format(', VARSHARELIST = ([', a,'] = ', i4, ')')
 135  format(', N = ', i8, ', E = ', i8, ', DATAPACKING = POINT',
     &     ', ZONETYPE = ', a)
 140  format(', VARSHARELIST = ([', a,'] = ', i4, '), ',
     &     'CONNECTIVITYSHAREZONE = 1')
 664  format("('nodes at ',e16.9,' days ',",i3,"(' : ', a))")
 190  format(a,', (Moles/kg H20)')
 200  format(a,', (Moles/kg Rock)')
 210  format(a,', (Moles/kg Vapor)')
 220  format('Free Ion ',a,', (Moles/kg H20)')
C     No zid
 333  format("(1x, i10.10, ", i3, "('", a, "', g16.9))")
C     Zid no coordinates
 334  format("(1x, i10.10, '", a, "', i4, ", i3, "('", a, "', g16.9))")
C     Coordinates and zid
 335  format("(1x, i10.10, ", i3, "('", a, "', g16.9), '", a, "', i4,",
     &     i3, "('", a, "', g16.9))")
      return
      end

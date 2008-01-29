      subroutine storsx
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
CD1 Store and retrieve element coefficients.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 14-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/storsx.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.18   Mon Mar 31 08:43:06 1997   gaz
CD2 new format from gable plus other changes
CD2 
CD2    Rev 1.17   Tue Jul 09 15:34:20 1996   zvd
CD2 Added format for writing problem title to store file.
CD2 
CD2    Rev 1.16   Tue May 14 14:32:22 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.15   Fri Apr 26 16:17:28 1996   gaz
CD2 lots of changes for mdnodes
CD2 
CD2    Rev 1.14   Thu Apr 18 13:31:34 1996   hend
CD2 Needed to allocate 3x for sx even if don't read all 3
CD2 
CD2    Rev 1.13   Thu Apr 18 13:12:12 1996   hend
CD2 Added parser for new stor file input
CD2 
CD2    Rev 1.12   Thu Apr 18 10:27:54 1996   gaz
CD2 changes to account for multiply defined nodes
CD2 
CD2    Rev 1.11   Fri Feb 16 11:36:10 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.10   Fri Feb 02 12:13:22 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.9   Fri Jan 12 17:59:08 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.8   11/15/95 16:17:40   gaz
CD2 changes for sx(iw,3) instead of sx(iw,6)
CD2 
CD2    Rev 1.7   06/07/95 09:07:08   llt
CD2 increased precision written on sx arrays
CD2 
CD2    Rev 1.6   04/25/95 10:14:30   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.5   01/28/95 13:56:00   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.4   11/29/94 18:22:36   llt
CD2 Changed length of jdate to 11 characters for ibm
CD2
CD2    Rev 1.3   06/20/94 11:09:08   zvd
CD2  
CD2 
CD2    Rev 1.2   03/18/94 15:55:56   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/28/94 11:52:36   zvd
CD2 Corrected problem of writing to coefficient storage file when it should be a
CD2 read only file.
CD2 
CD2    Rev 1.0   01/20/94 10:28:24   pvcs
CD2 original version in process of being certified
CD2 
c 12/22/94 gaz allocate storage space for coefficients
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
CD3   iatty                    O    File used for check information
CD4   iout                     O    File used for general code output
CD3   isstor                   I/O  File used for storing/reading element
CD3                                   coefficient values
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
CD4   iatty           INT      faai   Unit number for check file
CD4   iout            INT      faai   Unit number for output file
CD4   isstor          INT      faai   Unit number for element coefficient file
CD4   istrw           INT      fbb    Starting positions in sx array of finite
CD4                                     element coefficients for each node
CD4   iw              INT      faai   Number of storage locations needed to
CD4                                     store geometric input types
CD4   lda             INT      faai   Parameter which specifies if the geometric
CD4                                     coefficients are saved
CD4   nelm            INT      fbb    ?Initially information about nodes in each
CD4                                     element, later nodal connectivity
CD4                                     information
CD4   nelmdg          INT      fbb    Contains position of (i,i) element in
CD4                                     connectivity array
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4   nr              INT      param  Maximum space allowed for each finite
CD4                                     element coefficient array
CD4   sx              REAL*8   fbc    Contains finite element geometric
CD4                                     coefficients necessary for heat and mass
CD4                                     transfer simulation
CD4   sx1             REAL*8   fbc    Contains volume associated with each node
CD4   sxs             REAL*8   fbc    Contains more finite element geometric
CD4                                     coefficients (ie., those necessary for
CD4                                     the stress module)
CD4
CD4 Global Subprograms
CD4
CD4   read_sx
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
CD5   iwtotl          INT      Number of storage locations needed to
CD5                              store geometric input types
CD5   j               INT      Loop index
CD5   ncont           INT      Number of positions for which information needs
CD5                              to be stored
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
CD9 2.2 Finite-Element Coefficient Generation
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
CPS BEGIN  storsx
CPS 
CPS   IF coefficient storage file exists
CPS 
CPS      IF geometric coefficents should  be stored
CPS   
CPS         store the element coefficient information 
CPS      
CPS      ELSE if geometric coefficients should be read
CPS   
CPS         switch on format type of coefficient file
CPS
CPS           Allocate space needed for coefficients
CPS           call read_sx to retrieve the element
CPS           coefficient information
CPS           write storage requirements to output file 
CPS      
CPS           IF check file is being used 
CPS              write storage requirements to check file 
CPS           END IF
CPS      
CPS           IF space needed for coefficients exceeds maximum allocated
CPS              write termination message to output and error files
CPS              IF tty output is enabled
CPS                 write termination message to terminal
CPS              END IF
CPS              terminate program
CPS           END IF
CPS      
CPS      END IF
CPS      
CPS   ELSE
CPS      
CPS      write termination message to output and error files
CPS      IF tty output is enabled
CPS         write termination message to terminal
CPS      END IF
CPS      terminate program
CPS         
CPS   END IF
CPS   
CPS END  storsx
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      use comxi
      implicit none
c
      integer i, iwtotl, j, ncont,  neq_old,  ncoef
      integer ncont_new, neq_save, ncont_primary
c parser variables
      character*80 input_msg
      integer msg(10)
      integer nwds,sehtemp
      real*8 xmsg(10)
      integer imsg(10)
      character*32 cmsg(10)
c local
      logical opend
      logical exists
      integer ilen, rlen, flen
      integer ityp, max_con  
      character*100 filename, tail
      character*72 cline
      character*32 sxformat
      character*3 stat_var


c
c BEGIN
! Has a storage file been specified?
 9000 if (isstor .ne. 0) then

         if (lda .lt. 0) then
!Storing coefficients
            inquire(unit=isstor,name=filename,form=sxformat)
! If the file already contains data we don't want to overwrite it 
            stat_var='old'
            if(sxformat(1:9).eq.'FORMATTED') then
               read(isstor,'(a72)',end=900) cline
            else if(sxformat(1:11).eq.'UNFORMATTED') then
               read(isstor,end=900) cline
            endif
c This file already contains a stor file so a new one should be created 
            if (iout .ne. 0) write(iout,901)
            if(iptty.ne.0) write(iptty,901)
 901        format(1x,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     &         /,1x,'>>> changing name of new *.stor (old file exists)',
     &         /,1x,'>>> new file name is fehmn_temp.stor')
c Do not overwrite this one until we a sure
            inquire(file='fehmn_temp.stor',exist=exists)
            if(exists) then
               write (ierr, 902)
               if (iout .ne. 0) write(iout,902)
               if(iptty.ne.0) write(iptty,902)
 902           format(1x,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     &              /,1x,'>>> name fehmn_temp.stor is used : stopping',
     &              /,1x,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
               stop
            endif
            filename(16:72) =
     &      '                                                       '
            filename(1:15)='fehmn_temp.stor'
            stat_var='new'
 900        continue
            close(isstor)

            if (lda .eq.-1)  then
C**** store coefficients (formatted) ****
               open(isstor,file=filename,status=stat_var,
     &              form='formatted')

 1000          format(a30, 3x, a8, 3x, a8)

               if(gdpm_flag.ne.0) then
c gdpm nodes exist just stor primary nodes
                  iwtotl = iw
                  ncont_primary  = nelm(neq_primary + 1)

                  write(isstor, 1000)  verno, jdate, jtime
                  write(isstor, '(a80)') wdd
                  if(istrs.ne.0) then
                     ncoef=6
                  else
                     if(isoy.eq.1) then
                        ncoef=1
                     else if(isoy.eq.2) then
                        ncoef=3
                     endif
                  endif
                  write(isstor, '(5i10)'  ) iwtotl, neq_primary, 
     &                 ncont_primary, ncoef, max_con
                  write(isstor, '(5(1pe20.10))') 
     &                 (sx1(i), i = 1, neq_primary)
                  write(isstor, '(5i10)'  ) 
     &                 (nelm(i), i = 1, ncont_primary)
                  write(isstor, '(5i10)'  ) 
     &                 (istrw(i), i = 1, ncont_primary)
                  write(isstor, '(5i10)'  ) 
     &                 (nelmdg(i), i = 1, neq_primary)

                  if (istrs .eq. 0)  then

                     if(isoy.eq.2) then
                        write(isstor,'(5(1pe20.12))') 
     &                       (sx(j,1)+sx(j,2)+sx(j,3), j = 1, iwtotl)
                        write(isstor,'(5(1pe20.12))') (0.0d00 
     &                       , j = 1, iwtotl)
                        write(isstor,'(5(1pe20.12))') (0.0d00 
     &                       , j = 1, iwtotl)

                     else if(isoy.eq.1) then
                        if(icnl.eq.0) then
                           write(isstor,'(5(1pe20.12))')
     &                          (3.0d00*sx(j,1), j = 1, iwtotl)
                        else               
                           write(isstor,'(5(1pe20.12))')
     &                          (2.0d00*sx(j,1), j = 1, iwtotl)
                        endif
                        
                     endif
                  else
                     write(isstor,'(5g20.12)') (sx(j,1)+sx(j,2)+sx(j,3)
     &                    , j = 1, iwtotl)
                     do i = 2, 6
                        write(isstor,'(5(1pe20.12))') (0.0, 
     &                       j = 1, iwtotl)
                     end do
                     do i = 1, 6
                        write(isstor,'(5(1pe20.12))') (sxs(j,i),
     2                       j = 1, iwtotl)
                     end do
                  end if
               else
c no gdpm nodes
                  iwtotl = iw
                  ncont  = nelm(neq + 1)

                  write(isstor, 1000)  verno, jdate, jtime
                  write(isstor, '(a80)') wdd
                  if(istrs.ne.0) then
                     ncoef=6
                  else
                     if(isoy.eq.1) then
                        ncoef=1
                     else if(isoy.eq.2) then
                        ncoef=3
                     endif
                  endif
                  write(isstor, '(5i10)'  )
     &                 iwtotl, neq, ncont, ncoef, max_con
                  write(isstor, '(5(1pe20.12))') 
     &                 (sx1(i), i = 1, neq)
                  write(isstor, '(5i10)'  ) 
     &                 (nelm(i), i = 1, ncont)
                  write(isstor, '(5i10)'  ) 
     &                 (istrw(i), i = 1, ncont)
                  write(isstor, '(5i10)'  ) 
     &                 (nelmdg(i), i = 1, neq)

                  if (istrs .eq. 0)  then

                     if(isoy.eq.2) then
                        write(isstor,'(5(1pe20.12))') 
     &                       (sx(j,1)+sx(j,2)+sx(j,3), j = 1, iwtotl)
                        write(isstor,'(5(1pe20.12))') (0.0d00 
     &                       , j = 1, iwtotl)
                        write(isstor,'(5(1pe20.12))') (0.0d00 
     &                       , j = 1, iwtotl)

                     else if(isoy.eq.1) then
                        if(icnl.eq.0) then
                           write(isstor,'(5(1pe20.12))')
     &                          (3.0d00*sx(j,1), j = 1, iwtotl)
                        else               
                           write(isstor,'(5(1pe20.12))')
     &                          (2.0d00*sx(j,1), j = 1, iwtotl)
                        endif
                        
                     endif
                  else
                     write(isstor,'(5g20.12)') (sx(j,1)+sx(j,2)+sx(j,3)
     &                    , j = 1, iwtotl)
                     do i = 2, 6
                        write(isstor,'(5(1pe20.12))') (0.0, 
     &                       j = 1, iwtotl)
                     end do
                     do i = 1, 6
                        write(isstor,'(5(1pe20.12))') (sxs(j,i),
     2                       j = 1, iwtotl)
                     end do
                  end if
               endif

            else if(lda.eq.-2) then
C**** store coefficients (unformatted) ****
               open(isstor,file=filename,status=stat_var,
     &              form='unformatted')

               if(gdpm_flag.ne.0) then
c gdpm nodes exist just stor primary nodes
                  iwtotl = iw
                  ncont_primary  = nelm(neq_primary + 1)

                  write(isstor, 1000)  verno, jdate, jtime
                  write(isstor, '(a80)') wdd
                  if(istrs.ne.0) then
                     ncoef=6
                  else
                     if(isoy.eq.1) then
                        ncoef=1
                     else if(isoy.eq.2) then
                        ncoef=3
                     endif
                  endif
                  write(isstor) iwtotl, neq_primary,
     &                 ncont_primary, ncoef, max_con
                  write(isstor) 
     &                 (sx1(i), i = 1, neq_primary)
                  write(isstor) 
     &                 (nelm(i), i = 1, ncont_primary)
                  write(isstor) 
     &                 (istrw(i), i = 1, ncont_primary)
                  write(isstor) 
     &                 (nelmdg(i), i = 1, neq_primary)

                  if (istrs .eq. 0)  then

                     if(isoy.eq.2) then
                        write(isstor) (sx(j,1)+sx(j,2)+sx(j,3)
     &                       , j = 1, iwtotl)
                        write(isstor) (0.0d00 
     &                       , j = 1, iwtotl)
                        write(isstor) (0.0d00 
     &                       , j = 1, iwtotl)

                     else if(isoy.eq.1) then
                        if(icnl.eq.0) then
                           write(isstor)
     &                          (3.0d00*sx(j,1), j = 1, iwtotl)
                        else               
                           write(isstor)
     &                          (2.0d00*sx(j,1), j = 1, iwtotl)
                        endif

                     endif
                  else
                     write(isstor) (sx(j,1)+sx(j,2)+sx(j,3)
     &                    , j = 1, iwtotl)
                     do i = 2, 6
                        write(isstor) (0.0, j = 1, iwtotl)
                     end do
                     do i = 1, 6
                        write(isstor) (sxs(j,i),
     2                       j = 1, iwtotl)
                     end do
                  end if
               else
c no gdpm nodes

                  iwtotl = iw
                  ncont  = nelm(neq + 1)

                  cline = verno // ' ' // jdate // ' ' // jtime
                  write(isstor) cline
                  write(isstor) wdd(1:72)
                  max_con=0
                  if(istrs.ne.0) then
                     ncoef=6
                  else
                     if(isoy.eq.1) then
                        ncoef=1
                     else if(isoy.eq.2) then
                        ncoef=3
                     endif
                  endif
                  write(isstor) iwtotl, neq, ncont, ncoef, max_con
                  write(isstor) (sx1(i), i = 1, neq)
                  write(isstor) (nelm(i), i = 1, ncont)
                  write(isstor) (istrw(i), i = 1, ncont)
                  write(isstor) (nelmdg(i), i = 1, neq)

                  if (istrs .eq. 0)  then

                     if(isoy.eq.2) then
                        write(isstor) (sx(j,1)+sx(j,2)+sx(j,3)
     &                       , j = 1, iwtotl)
                        write(isstor) (0.0d00 
     &                       , j = 1, iwtotl)
                        write(isstor) (0.0d00 
     &                       , j = 1, iwtotl)

                     else if(isoy.eq.1) then
                        if(icnl.eq.0) then
                           write(isstor)                   
     &                          (3.0d00*sx(j,1), j = 1, iwtotl)
                        else               
                           write(isstor)                   
     &                          (2.0d00*sx(j,1), j = 1, iwtotl)
                        endif

                     endif
                  else
                     write(isstor) (sx(j,1)+sx(j,2)+sx(j,3)
     &                    , j = 1, iwtotl)
                     do i = 2, 6
                        write(isstor) (0.0, j = 1, iwtotl)
                     end do
                     do i = 1, 6
                        write(isstor) (sxs(j,i),
     2                       j = 1, iwtotl)
                     end do

                  end if
               end if
            endif
            if (iout .ne. 0) then
               write(iout,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!! '
               write(iout,*) 'Coefficients written to *.stor file '
               write(iout,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!! '
            end if
            if(iptty.gt.0) then
               write(iptty,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!! '
               write(iptty,*) 'Coefficients written to *.stor file '
               write(iptty,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!! '
            endif

         else if (lda .gt. 0) then

c**** retrieve coefficients ****

            if(lda.eq.1 .or. lda.eq.5 .or. lda.eq.6) then
               sxformat(1:5) = 'ascii'
               ityp = 2
            else if(lda.eq.2 .or. lda.eq.7 .or. lda.eq.8) then
               sxformat(1:5) = 'unfor'
               ityp = 3
            else if(lda.eq.3) then
               sxformat(1:5) = 'binar'
               ityp = 1
            else if(lda.eq.4) then
               sxformat(1:5) = 'unkno'
               ityp = 0
            endif

            inquire(unit=isstor,name=filename)
            inquire(unit=isstor,opened=opend)
            if (opend) close(isstor)

 1500       format (1x, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
 1501       format (1x, 'stor file has neq less than data file: quit')
 1502       format (1x, 'error in parsing beginning of stor file')
 1503       format (1x, 'stor file has unrecognized format: quit')
 1504       format (1x, 'Coefficients read from *.stor file ')

c-----------READ ASCII  format stor file ------------------------------
            if (sxformat(1:5).eq.'ascii') then

               if(iptty.ne.0) then
                 write(iptty,*)'WARNING: ASSUMING NEW ASCII STOR FORMAT'
               end if

               open(isstor,file=filename ,status='old',form='formatted')

               read (isstor,'(a72)', END=5000) cline
               read (isstor,'(a72)', END=5000) cline
               read (isstor,'(5i10)') iwtotl,neq_old,ncont,sehtemp,
     &              max_con

               if(neq_old.lt.neq_primary) then
c enable mdnodes
                  imdnode=-neq_old
               else if(neq_old.gt.neq_primary) then
c error occurs
                  write(ierr, 1500)
                  write(ierr, 1501)
                  write(ierr, 1500)
                  if (iout .ne. 0) then
                     write(iout, 1500)
                     write(iout, 1501)
                     write(iout, 1500)
                  end if
                  if(iptty.gt.0) then
                     write(iptty, 1500)
                     write(iptty, 1501)
                     write(iptty, 1500)
                  endif
                  stop
               endif
            
c     add 1 more conection to be used for BCtoBC and mdnode
               nr=iwtotl+1

c     call mdnode to arrange additional space 
c     for connectivity when we have mdnodes

               if(imdnode.ne.0) then
                  ncont_new=ncont
                  call md_nodes(4,0,ncont_new)
               else
                  ncont_new=ncont
               endif
c
               deallocate(nelm)
               allocate(nelm(ncont_new),istrw(ncont_new))
c gaz 11-09-2001
c allocation of istrw_itfc and istrw_cold done in startup
               nelmd=ncont_new

               read (isstor,*)  (sx1(i), i = 1, neq_old)
               read (isstor,*)  (nelm(i), i = 1, ncont)
               read (isstor,*)  (istrw(i), i = 1, ncont)
               read (isstor,*)  (nelmdg(i), i = 1, neq_old)

               if (istrs .eq. 0)  then
                  if (sehtemp.eq.1) then
                     isoy=1
                     isoz=1
                     allocate(sx(nr,1))
                     sx=0
                     call read_sx(nr,iwtotl,1,isstor,0,ityp)
                  else if (sehtemp.eq.3.or.sehtemp.eq.4) then
                     allocate(sx(nr,3))
                     call read_sx(nr,iwtotl,3,isstor,0,ityp)
                  else if (sehtemp.eq.6) then
                     allocate(sx(nr,6))
                     call read_sx(nr,iwtotl,6,isstor,0,ityp)
                  endif
               else
                  allocate(sx(nr,6))
                  call read_sx(nr,iwtotl,6,isstor,0,ityp)
                  allocate(sxs(nr,6))
                  call read_sx(nr,iwtotl,6,isstor,1,ityp)
               end if

c-----------READ UNFORMATTED  format stor file ------------------------
            elseif (sxformat(1:5).eq.'unfor') then
               if(iptty.ne.0) then
                  write(iptty,*)'WARNING: ASSUMING NEW UNFORMATTED',
     &                 ' STOR FORMAT'
               end if

               open(isstor,file=filename ,status='old',
     &              form='unformatted')

               read (isstor, END=5000) cline
               read (isstor, END=5000) cline
               read (isstor)  iwtotl, neq_old, ncont, sehtemp, max_con
c adding possibility that the *.stor file may have less nodes than
c what was defined in the data file. We assume that means additional
c nodes were defined for multiply connected situations. Further
c all the new nodes will be after neq
c
               if(neq_old.lt.neq_primary) then
c enable mdnodes
                  imdnode=-neq_old
               else if(neq_old.gt.neq_primary) then
c error occurs
                  write(ierr, 1500)
                  write(ierr, 1501)
                  write(ierr, 1500)
                  if (iout .ne. 0) then
                     write(iout, 1500)
                     write(iout, 1501)
                     write(iout, 1500)
                  end if
                  if(iptty.gt.0) then
                     write(iptty, 1500)
                     write(iptty, 1501)
                     write(iptty, 1500)
                  endif
                  stop
               endif
c
c always add one more conection to be used for BC to BC
c and mdnode connections
c
               nr=iwtotl+1
c
c call mdnode to arrange additional space for connectivity
c when we have mdnodes
c
               if(imdnode.ne.0) then
                  ncont_new=ncont
                  call md_nodes(4,0,ncont_new)
               else
                  ncont_new=ncont
               endif
c
               deallocate(nelm)
               allocate(nelm(ncont_new),istrw(ncont_new))
c gaz 11-09-2001
c allocation of istrw_itfc and istrw_cold done in startup
               nelmd=ncont_new

               read (isstor)  (sx1(i), i = 1, neq_old)
               read (isstor)  (nelm(i), i = 1, ncont)
               read (isstor)  (istrw(i), i = 1, ncont)
               read (isstor)  (nelmdg(i), i = 1, neq_old)
               if (istrs .eq. 0)  then
                  if (sehtemp.eq.1) then
                     isoy=1
                     isoz=1
                     allocate(sx(nr,1))
                     sx=0
                     call read_sx(nr,iwtotl,1,isstor,0,ityp)
                  else if (sehtemp.eq.3.or.sehtemp.eq.4) then
                     allocate(sx(nr,3))
                     call read_sx(nr,iwtotl,3,isstor,0,ityp)
                  else if (sehtemp.eq.6) then
                     allocate(sx(nr,6))
                     call read_sx(nr,iwtotl,6,isstor,0,ityp)
                  endif
               else
                  allocate(sx(nr,6))
                  call read_sx(nr,iwtotl,6,isstor,0,ityp)
                  allocate(sxs(nr,6))
                  call read_sx(nr,iwtotl,6,isstor,1,ityp)
               end if


c-----------READ OLD ASCII  format stor file --------------------------
c        if format not recognized, try old ascii format
            elseif (sxformat(1:7).eq.'unknown') then

               open(isstor,file=filename ,status='old',form='formatted')

               read (isstor, '(a80)', END=5000)     input_msg
               read (isstor, '(a80)', END=5000)     input_msg
               read (isstor, '(a80)') input_msg
               call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
               iwtotl=imsg(1)
               neq_old=imsg(2)
               ncont=imsg(3)
               if (nwds.eq.4) then
                  sehtemp=imsg(4)
               else if (nwds.eq.5) then
                  sehtemp=imsg(4)
               else if (nwds.eq.3) then
                  sehtemp=6
               else
                  write(ierr, 1500)
                  write(ierr, 1502)
                  write(ierr, 1500)
                  if (iout .ne. 0) then
                     write(iout, 1500)
                     write(iout, 1502)
                     write(iout, 1500)
                  end if
                  if(iptty.gt.0) then
                     write(iptty, 1500)
                     write(iptty, 1502)
                     write(iptty, 1500)
                  endif
                  stop
               endif
c
c adding possibility that the *.stor file may have less nodes than
c what was defined in the data file. We assume that means additional
c nodes were defined for multiply connected situations. Further
c all the new nodes will be after neq
c
               if(neq_old.lt.neq_primary) then
c enable mdnodes
                  imdnode=-neq_old
               else if(neq_old.gt.neq_primary) then
c error occurs
                  write(ierr, 1500)
                  write(ierr, 1501)
                  write(ierr, 1500)
                  if (iout .ne. 0) then
                     write(iout, 1500)
                     write(iout, 1501)
                     write(iout, 1500)
                  end if
                  if(iptty.gt.0) then
                     write(iptty, 1500)
                     write(iptty, 1501)
                     write(iptty, 1500)
                  endif
                  stop
               endif
c 
c always add one more conection to be used for BC to BC 
c and mdnode connections
c
               nr=iwtotl+1
c
c call mdnode to arrange additional space for connectivity
c when we have mdnodes
c        
               if(imdnode.ne.0) then
                  ncont_new=ncont
                  call md_nodes(4,0,ncont_new)
               else
                  ncont_new=ncont
               endif
c
               deallocate(nelm)
               allocate(nelm(ncont_new),istrw(ncont_new))
c gaz 11-09-2001
c allocation of istrw_itfc and istrw_cold done in startup
               nelmd=ncont_new

               read (isstor, *,err=990)  (sx1(i), i = 1, neq_old)
               read (isstor, *,err=990)  (nelm(i), i = 1, ncont)
               read (isstor, *,err=990)  (istrw(i), i = 1, ncont)
               read (isstor, *,err=990)  (nelmdg(i), i = 1, neq_old)
               if (istrs .eq. 0)  then
                  if (sehtemp.eq.1) then
                     isoy=1
                     isoz=1
                     allocate(sx(nr,1))
                     sx=0
                     call read_sx(nr,iwtotl,1,isstor,0,ityp)
                  else if (sehtemp.eq.3) then
                     allocate(sx(nr,3))
                     call read_sx(nr,iwtotl,3,isstor,0,ityp)
                  else if (sehtemp.eq.6) then
                     allocate(sx(nr,6))
                     call read_sx(nr,iwtotl,6,isstor,0,ityp)
                  endif
               else 
                  allocate(sx(nr,6))
                  call read_sx(nr,iwtotl,6,isstor,0,ityp)
                  allocate(sxs(nr,6))
                  call read_sx(nr,iwtotl,6,isstor,1,ityp)
               end if
       
c     made it without encountering read errors
               sxformat='oldascii'
 990           if (sxformat(1:7).eq.'unknown') then
                  write(ierr, 1500)
                  write(ierr, 1503)
                  write(ierr, 1500)
                  if (iout .ne. 0) then
                     write(iout, 1500)
                     write(iout, 1503)
                     write(iout, 1500)
                  end if
                  if(iptty.gt.0) then
                     write(iptty, 1500)
                     write(iptty, 1503)
                     write(iptty, 1500)
                  endif
                  stop
               endif

            endif
c-----END switch on reading stor formats--------------------------------
            if(iptty.ne.0) then
               write(iptty,*)'Coefficient file read with format ',
     &              sxformat(1:15)
            end if

c
c insure isotropy for non stress solutions
c
            if(isoy.eq.2) then
               if(gdpm_flag.ne.0) then
                  neq_save = neq
                  neq = neq_primary
                  call sx_combine(0)
                  neq = neq_save
               else
                  call sx_combine(0)
               end if
c done now in startup.f
c     call sx_combine(1)
               call sx_combine(2)
            else if(isoy.eq.1) then
               call sx_combine(-3)
c done now in startup.f
c     call sx_combine(1)
               call sx_combine(3)
            endif

c**** printout positions required for geometric coefficients ****
c
            if (iout .ne. 0) then
               write(iout, 1500)
               write(iout, 1504)
               write(iout, 1500)
            end if
            if(iptty.gt.0) then
               write(iptty, 1500)
               write(iptty, 1504)
               write(iptty, 1500)
            endif
c
            if (iout .ne. 0) write(iout, 2000) nr, nr
            if (iptty .gt. 0) write(iptty, 2000) nr,nr       
 2000       format(/,1x,'storage for geometric coefficients ', i10,
     *           ' in common(nr) ', i10)

            if (iwtotl .gt. nr)  then

               if (iout .ne. 0) write(iout, 3000)
               write(ierr, 3000)
               if (iptty .ne. 0) write(iptty, 3000)
 3000          format(/, 1x, 'program terminated because of ',
     *              'insufficient storage')

               stop

            end if
c
c
c calculate additional storage if necessary for multiply defined nodes
c
            if(imdnode.ne.0) then
               call md_nodes(5,0,ncont_new)
            endif

         end if

         close  (isstor)
         
      else
! A storage file was not specified 
         if (iout .ne. 0) write(iout, 4000)
         write(ierr, 4000)
         if (iptty .ne. 0) write(iptty, 4000)
 4000    format(/, 1x, 'program terminated because ',
     *        'coefficient storage file not specified')

         stop

      end if

! Write an unformatted file based on ascii data read in
      if (lda .ge. 5) then
         stat_var = '   '
         iw = iwtotl
         flen = len_trim (filename)
         do i = flen, 1, -1
            if (filename(i:i) .eq. '.') then
               rlen = i-1
               ilen = flen-rlen
               tail(1:ilen) = filename(i:flen)
               goto 9001
            end if
            if (i .eq. 1) rlen = flen
         end do
 9001    if ((flen + 3) .gt. 100) rlen = rlen-3
         filename(rlen+1:rlen+3) = 'UNF'
         filename(rlen+4:rlen+4+ilen) = tail(1:ilen)
         if (lda .eq. 5 .or. lda .eq. 7) then
            lda = -2
            open(isstor,file=filename,form='unformatted')
         else if (lda .eq. 6 .or. lda .eq. 8) then
            filename(rlen+1:rlen+3) = 'FOR'
            lda = -1
            open(isstor,file=filename,form='formatted')
         end if
         goto 9000
      end if
      return

! Data could not be read
 5000 continue
 5001 format (/, 1x, 'program terminated because ',
     *     'data could not be read from storage file')
      if (iout .ne. 0) write(iout, 5001)
      write(ierr, 5001)
      if (iptty .ne. 0) write(iptty, 5001)
      stop  
 
      end

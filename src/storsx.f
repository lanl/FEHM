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
      integer ityp
      integer :: max_con = 0
      character*100 filename, tail
      character*72 cline
      character*3 stat_var


c     
c     BEGIN
c     Has a storage file been specified?
      if (isstor .ne. 0) then

         if (lda .lt. 0) then
c     Storing coefficients

            call storsx_write

         else if (lda .gt. 0) then
c**** retrieve coefficients ****

            call storsx_read
            
            if (lda .ge. 5) then
c     Write file based on data read in
               stat_var = '   '
               iw = iwtotl
               filename = nmfil(12)
               flen = len_trim (filename)
               do i = flen, 1, -1
                  if (filename(i:i) .eq. '.') then
                     rlen = i-1
                     ilen = flen-rlen
                     tail(1:ilen) = filename(i:flen)
                     exit
                  end if
                  if (i .eq. 1) rlen = flen
               end do
               if ((flen + 3) .gt. 100) rlen = rlen-3
               filename(rlen+1:rlen+3) = 'UNF'
               filename(rlen+4:rlen+4+ilen) = tail(1:ilen)
               if (lda .eq. 5 .or. lda .eq. 7) then
c     Write an unformatted file based on data read in
                  lda = -2
                  open(isstor,file=filename,form='unformatted')
               else if (lda .eq. 6 .or. lda .eq. 8) then
c     Write a formatted file based on data read in
                  filename(rlen+1:rlen+3) = 'FOR'
                  lda = -1
                  open(isstor,file=filename,form='formatted')
               end if
               
               call storsx_write

            end if
         end if
      else
c     A storage file was not specified 

         if (iout .ne. 0) write(iout, 4000)
         write(ierr, 4000)
         if (iptty .ne. 0) write(iptty, 4000)
 4000    format(/, 1x, 'program terminated because ',
     *        'coefficient storage file not specified')

         stop

      end if

      end subroutine storsx

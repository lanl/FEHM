      subroutine write_avs_ucd_header
     1       (lu,verno,jdate,wdd,neq,nei,num_ndata,num_cdata,num_mdata)
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
CD1 Output AVS UCD header information.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 $Log:   /pvcs.config/fehm90/src/write_avs_ucd_header.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:46   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Fri Feb 02 14:24:52 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   01/20/95 12:47:30   tam
CD2 No change.
CD2 
CD2    Rev 1.1   11/29/94 13:42:42   llt
CD2 Removed protying, so could run under cc compile, instead of c compiler. 
CD2 Used #ifdef to determine format for routine name, depending on machine.
CD2 (Changes made by tam.)
CD2 
CD2    Rev 1.0   08/23/94 15:35:34   llt
CD2 Original version
CD2
CD2 10-SEP-93    Carl Gable     22      Initial implementation.
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

      implicit none

      integer lu,neq,nei,num_ndata,num_cdata,num_mdata
      character*30 verno
      character*11 jdate
      character*80 wdd

c     use formats to force 0 leading blanks in the avs file header
c     To left justify - aw must be <= length of string
      write(lu,300)verno,jdate
 300  format('# ',a30,3x,a11)
      write(lu,355)
     1'#  AVS UNSTRUCTURED CELL DATA (UCD) FROM FEHM                  '
      write(lu,310)wdd
  310 format('# ',a72)

c     directions for concatination of files
      write(lu,355)
     1'# ************************************************************ '
      write(lu,355)
     1'#  To prepare files for input to avs one must                  '
      write(lu,355)
     1'#  concatinate  header/geometry/node_value files.              '
      write(lu,355)
     1'#  For example, if your FEHM input file was fe.dat,            '
      write(lu,355)
     1'#  headers are fe10001_sca_head fe10001_vec_head, ...,         '
      write(lu,355)
     1'#  mesh geometry will be in fe10001_geo,                       '
      write(lu,355)
     1'#  field output will be in fe10001_sca_node,                   ' 
      write(lu,355)
     1'#            fe10001_vec_node, fe10001_con_dual_node           '
      write(lu,355)
     1'#                                                              '
      write(lu,355)
     1'#  A UCD input file can be produced using                      '
      write(lu,355)
     1'#   cat fe10001_sca_head fe10001_geo fe10001_sca_node >        '
      write(lu,355)
     1'#       fe10001_sca_node.inp                                   '
      write(lu,355)
     1'#                                                              '
      write(lu,355)
     1'#  The UNIX foreach command is useful for processing           '
      write(lu,355)
     1'#   multiple files. Also use the shell script fehm2avs         '
      write(lu,355)
     1'#   to perform automatic processing of all output.             '
      write(lu,355)
     1'# ************************************************************ '
      write(lu,320)neq,nei,num_ndata,num_cdata,num_mdata

      close (lu)

  320 format(i10.10,2x,i10,2x,i10,2x,i10,2x,i10)
  355 format(a55)

      return
      end

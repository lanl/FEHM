      module comxi
!    comxi
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Include file containing passed parameters used for input/output.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 09-13-93     B. Robinson    00022    Initial Implementation
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comxi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:02 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   06/20/94 11:02:30   zvd
!D2 Added name for error output file.
!D2 
!D2    Rev 1.1   03/18/94 16:23:28   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:22   pvcs
!D2 original version in process of being certified
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!***********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   None
!D4
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block  Description
!D4
!D4   ***** COMMON Block faax variables *****
!D4   blank           CHAR     faax   Empty character string for file names
!D4   cform           CHAR ARY faax   I/O file formats
!D4   cstats          CHAR ARY faax   I/O file status
!D4   cuser           CHAR     faax   String written for user subroutinenumber
!D4                                     use
!D4   iowork          CHAR ARY faax   File use/type designator
!D4   nmfil           CHAR ARY faax   I/O file names:
!D4                                   1 - Control file name
!D4                                   2 - Main input file name
!D4                                   3 - Coordinate input file name
!D4                                   4 - Zone input file name
!D4                                   5 - Main output file name
!D4                                   6 - Restart input file name
!D4                                   7 - Restart output file name
!D4                                   8 - Simulation history output file name
!D4                                   9 - Solute history output file name
!D4                                   10 - Contour plot output file name
!D4                                   11 - Dual porosity or dpdp contour plot
!D4                                          output file name
!D4                                   12 - Coefficient storage output file name
!D4                                   13 - Input check output file name
!D4                                   14 - Error output file name
!D4                                   15 - Pest output file name
!D4                                   16 - Pest auxilliary output file name
!D4                                   17 - Streamline particle tracking
!D4                                          output file 1 name
!D4                                   18 - Streamline particle tracking
!D4                                          output file 2 name
!D4                                   19 - Streamline particle tracking
!D4                                          output file 3 name
!D4                                   20 - Streamline particle tracking
!D4                                          output file 4 name
!D4                                   21 - Streamline particle tracking
!D4                                          output file 5 name
!D4                                   22 - Streamline particle tracking
!D4                                          output file 6 name
!D4                                   23 - Streamline particle tracking
!D4                                          output file 7 name
!D4                                   24 - Extracted flow model (flow 
!D4                                          macro) file name 
!D4                                   25 - Streamline particle tracking
!D4                                          output file 8 name
!D4                                   26 - Streamline particle tracking
!D4                                          output file 9 name
!D4                                   27 - Column node data file for wtsi
!D4                                   28 - Unformatted file with symbolic
!D4                                          factorization for incomplete
!D4                                          factorization for the solver
!D4                                          preconditioner (nop.temp)
!D4   nmfily          CHAR ARY faax   File messages
!D4   root_name       CHAR            Root name for optional output files
!D4   suffix          CHAR ARY faax   Default file suffixes
!D4   
!D4   ***** COMMON Block faay variables *****
!D4   ex              LOGICAL  faay   Logical for existence check (T/F)
!D4   isw             INT ARY  faay   Error flag for opening files,
!D4                                     isw = 1 -- No error
!D4                                     isw = 2 -- Error
!D4   nmmax           INT      faay   Number of I/O files
!D4   nmmaxa          INT      faay   Number of additional I/O files
!D4   nufilb          INT ARY  faay   File unit numbers
!D4   
!***********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!***********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!***********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 None
!D9
!***********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS 
!PS None
!PS
!***********************************************************************

      logical ex

      integer nmmax
! gaz 071020 increased size  of nmmaxa , increased to 15 081921 for co2wh     
      parameter(nmmax = 15)

      integer nmmaxa
      parameter(nmmaxa = 18)

      integer isw(nmmax+nmmaxa), nufilb(nmmax+nmmaxa)

      character* 6  suffix(nmmax+nmmaxa)
      character* 6  iowork(nmmax+nmmaxa)
      character* 7  cstats(nmmax+nmmaxa)
      character* 9  cuser
      character*11  cform(nmmax+nmmaxa)
      character*100 blank, nmfil(nmmax+nmmaxa), nmfily(3), root_name

      end module comxi

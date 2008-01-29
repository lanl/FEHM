      subroutine  rddpdp
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
CD1 To read the input data for a dpdp solution.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-08-93     B. Robinson    00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/rddpdp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jun 12 15:21:16 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.4   Mon Jun 03 11:18:22 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.3   Fri May 31 10:40:48 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.2   Thu Feb 01 15:33:20 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 15:47:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:58   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use      Description
CD3 
CD3 None
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier   Type        Description
CD4 
CD4 
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 Identifier   Type        Description
CD4 
CD4 inpt, idpdp, ipvolf1, volf1, ipapuv1, apuv1, ischk,narrays, pointer,
CD4 itype, default, igroup, ireturn,
CD4 macroread
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type       Description
CD4 
CD4 initdata  N/A        Reads in data and sets array values
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 narrays      int         Number of arrays being read in
CD5 pointer      int         Integer array of pointer values for
CD5                              variables being read
CD5 itype        int         Array of variable types being input
CD5 default      real*8      Array of default values for input arrays
CD5 ierrflag     int         Flag denoting whether input error
CD5                              information is to be written
CD5 macro        char        Current macro being read
CD5 igroup       int         Current group number in this macro
CD5 ireturn      int         Returned error flag from input subroutine
CD5 npoints      int         Number of points at which to set array
CD5                             values
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 N/A
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
CD9 2.4.9 Double-porosity/double-permeability formulation
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN rddpdp
CPS 
CPS Read in dpdp flag
CPS Set parameters for reading in fracture porosities
CPS       
CPS initdata - read and set fracture porosity values
CPS Set parameters for reading in fracture areas per unit volume
CPS initdata - read and set fracture areas per unit volume
CPS 
CPS Set flag to denote that the dpdp macro has been called
CPS       
CPS END rddpdp
CPS 
C**********************************************************************

      use combi
      use comhi
      use comdti
      use comai
      use comki
      implicit none
      
      integer npoints

c     read in dpdp parameter
      read(inpt,*) idpdp

c     read in fracture porosities
      narrays = 1
      itype(1) = 8
      default(1) = 1.
      npoints = neq
      macro = "dpdp"
      igroup = 2
      
      call initdata2(inpt, ischk, npoints, narrays,  itype,
     &     default, macroread(6), macro, igroup, ireturn,
     2     r8_1=volf1(1:npoints) )
      
c     Only reset those that change before reading apuv1
      default(1) = 10.
      igroup = 3
      
      call initdata2(inpt, ischk, npoints, narrays, itype,
     &     default, macroread(6), macro, igroup, ireturn,
     2     r8_1=apuv1(1:npoints) )
      
      macroread(6) = .TRUE.
      
      return
      end



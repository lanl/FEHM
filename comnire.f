      module  comnire
!    comnire
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
!D1 Include file containing passed parameters used for input/output.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comnire.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 08:57:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 N/A
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

      integer  i_yama, imod_ten, kkk 
      integer  imod3_yama, imod4_yama, imod5_yama
      integer  imod_yama(50)

      real*8  yama_ww
      real*8  yama_ff
      real*8  yama_phi
c
      character(100) ::  yama_filename1, yama_filename2
c
      end module comnire

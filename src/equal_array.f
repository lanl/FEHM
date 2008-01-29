      subroutine equal_array(a,b,neq,idof,nrhs1,nrhs2)
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
CD1 To set array a equal to array b.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 3-2-94       G. Zyvoloski   97      Initial implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/equal_array.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Thu Sep 12 08:25:10 1996   robinson
CD2 Prolog Changes
CD2 
CD2    Rev 1.5   Fri May 31 15:02:00 1996   gaz
CD2 added variable so each array can have a different nrhs
CD2 
CD2    Rev 1.4   05/12/94 10:04:36   llt
CD2 corrected pvcs log info
CD2 
CD2    Rev 1.3   05/11/94 16:18:38   llt
CD2 bug fixes - gaz
CD2 
CD2    Rev 1.2   03/23/94 14:43:54   robinson
CD2 Added prologs
CD2
CD2    Rev 1.1   03/18/94 16:05:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   03/02/94 08:46:26   pvcs
CD2 original version
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 a            real*8   O      Array being set to the value of array b
CD3 b            real*8   I      Array that array is being set to
CD3 neq          integer  I      Number of entries in the array for
CD3                                  each degree of freedom
CD3 idof         integer  I      Number of degrees of freedom
CD3 nrhs1        integer  I      Pointer integer array denoting the
CD3                                  first entry in the a array for
CD3                                  each degree of freedom
CD3 nrhs2        integer  I      Pointer integer array denoting the
CD3                                  first entry in the b array for
CD3                                  each degree of freedom
CD3 
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 NONE
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4 
CD4 NONE
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 NONE
CD4 
CD4 Global Subprograms
CD4 
CD4 NONE
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 NONE
CD5 
CD5 Local Types
CD5
CD5 NONE
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 j            int         Do loop index for each degree of freedom
CD5 i            int         Do loop index for each value in array
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 N/A
CD6
C**********************************************************************
CD7
CD7 SPECIAL COMMENTS
CD7 
CD7 No requirements are directly satisfied by this routine, but it is
CD7 one operation of the overall solution of the linear equations.
CD7
C**********************************************************************
CD8
CD8 REQUIREMENTS TRACEABILITY
CD8 
CD8 N/A
CD8
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN equal_array
CPS 
CPS FOR each degree of freedom
CPS   FOR each element in the array
CPS     Set value in array a to the value in array b
CPS   ENDFOR
CPS ENDFOR
CPS 
CPS END equal_array
CPS
C**********************************************************************
c
c set two arrays equal                      
c
      implicit none
      integer neq,i,j,idof,nrhs1(*),nrhs2(*),n1,n2
      real*8 a(*),b(*)   
      do j=1,idof
         n1 = nrhs1(j)
         n2 = nrhs2(j)
         if(n1.eq.n2) then
            n1=n2+1
            n2=n2+neq
            do i=n1,n2
              a(i)=b(i)
            enddo
         else
            do i=1,neq
              a(i+n1)=b(i+n2)
            enddo
         endif
      enddo
      return
      end

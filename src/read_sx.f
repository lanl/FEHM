      subroutine read_sx(nr_dim,nr,num2,isstor,sxflag,ityp)
!***********************************************************************
!  Copyright, 1995, 2004,  The  Regents of the University of California.
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
CD1  PURPOSE
CD1
CD1  To read finite element coefficients.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/read_sx.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:08   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Fri Apr 26 15:57:54 1996   gaz
CD2 passed array dimension through (different than nr)
CD2 
CD2    Rev 1.4   Wed Jan 17 09:58:40 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.3   Thu Jan 11 13:21:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.2   Thu Jan 11 12:49:38 1996   gaz
CD2 added end of file check
CD2 
CD2    Rev 1.1   04/25/95 10:08:02   llt
CD2 retrieved lost log history information
CD2
CD2    Rev 1.0   01/28/95 13:55:24   llt
CD2 water balance equation was modified
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.2 Finite-Element Coefficient Generation
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4  Requirements from SDN: 10086-RD-2.20-00
CD4    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD4    FEHM Application Version 2.20
CD4
C***********************************************************************

      use comai , only : ierr
      use combi
      implicit none
c
      integer nr_dim,nr,num2,isstor,ityp
c      real*8 sx(nr_dim,num2)
      integer i,j, sxflag
c
c BEGIN

c-----read unformatted file
      if(ityp.eq.3) then

        if(sxflag.eq.0) then
        do i = 1, num2
          read (isstor, end = 999)  (sx(j,i), j = 1, nr)
        end do
        else
        do i = 1, num2
          read (isstor, end = 999)  (sxs(j,i), j = 1, nr)
        end do
        end if

c-----read ascii file
      elseif(ityp.eq.2) then

        if(sxflag.eq.0) then
        do i = 1, num2
          read (isstor,*, end = 999)  (sx(j,i), j = 1, nr)
        end do
        else
        do i = 1, num2
          read (isstor,*, end = 999)  (sxs(j,i), j = 1, nr)
        end do
        end if

c-----read binary file
      elseif(ityp.eq.1) then
          write(ierr,*) 'readsx: Binary format not implemented.' 
          goto 999

c-----read old style ascii file
      else

        if(sxflag.eq.0) then
        do i = 1, num2
          read (isstor, *, end = 999)  (sx(j,i), j = 1, nr)
        end do
        else
        do i = 1, num2
          read (isstor, *, end = 999)  (sxs(j,i), j = 1, nr)
        end do
        end if

      endif

 999  continue
      return
      end

      subroutine storage_derivatives(iz,idum)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine allocates memory for derivatives.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/storage_derivatives.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.0   Wed May 29 16:51:30 1996   llt
CD2 changed file file name to match subroutine name
CD2 
CD2    Rev 1.13   Fri Feb 16 10:19:18 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.12   Tue Jan 30 15:53:58 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.11   Fri Jan 12 17:51:54 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.10   Wed Jan 10 14:40:10 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.9   08/07/95 13:13:16   awolf
CD2 fixed typo
CD2 
CD2    Rev 1.8   08/07/95 11:07:22   awolf
CD2 drc allocation changed fro n0 to n7a for multi-species
CD2 
CD2    Rev 1.7   04/25/95 09:49:24   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.6   01/28/95 13:55:14   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.5   07/28/94 10:07:42   robinson
CD2 Took handling of denci out of this routine, put in allocmem
CD2
CD2    Rev 1.4   03/24/94 08:02:40   llt
CD2 Added $Log to top
CD2
CD2    Rev 1.3   03/23/94 14:43:16   robinson
CD2 Additional cleanup of memory management
CD2
CD2    Rev 1.2   03/18/94 15:59:38   gaz
CD2 Added solve_new and cleaned up memory management.
CD2
CD2    Rev 1.1   03/03/94 08:34:54   llt
CD2 add $Log to top
CD2
CD2    Rev 1.0   03/03/94 08:32:18   pvcs
CD2 original version
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable.  See Special Comments
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  This is a general purpose utility routine called as needed.
CD4  
C***********************************************************************

      use comci
      use comdti
      use davidi 
      implicit none

      integer startpos,endpos
      integer n0_subst
      integer n7_subst
      integer n7a_subst
      integer kgmres_subst
      integer nbd_subst
      integer nstorepiv_subst
      integer n0_irdof 
      integer n0_vapor
      integer iz, idum
      if(iz.le.-1) then
         if(iz.eq.-1) then
            n0_irdof = n0
            if (irdof .eq. 13) then
               n0_vapor = 1
            else
               n0_vapor = n0_irdof
            end if
            n0_subst = n0
            n7_subst = n7
            n7a_subst = n7a
            kgmres_subst = kgmres
            nbd_subst = nbd
            nstorepiv_subst = nstorepiv
         else
c     small sized array used as placeholder for the case of
c     no fluid flow solution - see also startup.f
            n0_subst = 1
            n7_subst = 1
            n7a_subst = 1
            n0_irdof = 1 
            n0_vapor = 1
            kgmres_subst = 1
            nbd_subst = 1
            nstorepiv_subst = 1
         end if
c     Allocate big array
         if(idum.eq.1) allocate(bigblock(nbigblock))
         bigblock = 0
         
c     Set pointers to sections of the arrays for the solver storage
         startpos = 1
         endpos = kgmres_subst
         gmres => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + nbd_subst
         b => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + nstorepiv_subst
         piv => bigblock(startpos:endpos)
         
         
c     Set pointers to sections of the arrays for the derivatives
         startpos = 1
!         endpos = n7_subst
!         akc => bigblock(startpos:endpos)
         
!         startpos = endpos + 1
!         endpos = endpos + n7_subst
!         danl => bigblock(startpos:endpos)
         
!         startpos = endpos + 1
!         endpos = endpos + n7_subst
!         danv => bigblock(startpos:endpos)
         
!         startpos = endpos + 1
!         endpos = endpos + n0_irdof
         endpos = n0_vapor
         ddvac => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         ddvae => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         ddvap => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         deef => bigblock(startpos:endpos)
         rlf => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         delef => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         delf => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         denei => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_subst
         deni => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         denvac => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         denvae => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         denvap => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         depf => bigblock(startpos:endpos)
         drlef => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         deqh => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         devef => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + 2*n0_vapor
         devf => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + 2*n0_irdof
         dgle => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_subst
         dglp => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         dgve => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         dgvp => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + 2*n0_irdof
         dile => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_subst
         dilp => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         dive => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         divp => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         dmef => bigblock(startpos:endpos)
         rvf => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_subst
         dmpf => bigblock(startpos:endpos)
         drlpf => bigblock(startpos:endpos)
         
!         startpos = endpos + 1
!         endpos = endpos + n0_irdof
!         dpcef => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_subst
         dq => bigblock(startpos:endpos)
         drvef => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         dqh => bigblock(startpos:endpos)
         drvpf => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_subst
         dqt => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n7a_subst
         drc => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         dstm=> bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         dtpa => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         dtpae => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         dva => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_irdof
         enlf => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         enva => bigblock(startpos:endpos)
         
         startpos = endpos + 1
         endpos = endpos + n0_vapor
         envf => bigblock(startpos:endpos)
         
      end if
      
      if(iz.eq.0) then
         
c     Zero out derivative arrays
!         akc=0
!         danl=0
!         danv=0
         ddvac=0
         ddvae=0
         ddvap=0
         deef=0
         delef=0
         delf=0
         denei=0
         deni=0
         denvac=0
         denvae=0
         denvap=0
         depf=0
         deqh=0
         devef=0
         devf=0
         dgle=0
         dglp=0
         dgve=0
         dgvp=0
         dile=0
         dilp=0
         dive=0
         divp=0
         dmef=0
         dmpf=0
!         dpcef=0
         dq=0
         dqh=0
         dqt=0
         drc=0
         dstm=0
         dtpa=0
         dtpae=0
         dva=0
         enlf=0
         enva=0
         envf=0
         
      end if
      
      if (iz.eq.1) then
         endpos = kgmres + nbd + nstorepiv
         bigblock(1:endpos) = 0.d0
c         gmres = 0.d0
c         b = 0.d0
c         piv = 0.d0
      endif
      
      return
      end

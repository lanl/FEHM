      subroutine  plotc1( ifg,ii )
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
CD1 To write tracer concentration to the trc file for the current time
CD1 and solute number.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/plotc1.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:38   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:36   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:56   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2  
CD2    Rev 1.4   Thu Feb 01 14:13:18 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   08/07/95 11:37:32   awolf
CD2 Sets breakthrough output capability for Ptrack
CD2 
CD2    Rev 1.2   04/21/95 15:54:44   llt
CD2 added node information to *.trc output file for time history plots
CD2 
CD2    Rev 1.1   03/18/94 16:15:46   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:24   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Name       Type        Description
CD3 
CD3 ifg         I          Parameter indicating whether the
CD3                        information is to be written at this time
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                Description
CD3
CD3 tape number istrc   File number associated with trc file
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 istrc, nspeci, nsp, icns, days, m, nskw, npn, an, anv
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop index
CD5 mi           int         Node number
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following functions are carried out in this routine:
CD6 
CD6   The code checks to see if there is a trc file for this run, and
CD6   if there is:
CD6   
CD6     The code checks to see if information is to be written this
CD6     time step, or if this call is only for writing the node
CD6     information and the number of solutes at the top of the file.
CD6     
CD6       The code writes the number of solutes if required, or
CD6       
CD6       The code writes the time, current solute number, and
CD6       concentrations at each node number specified.
CD6     
CD6  The code thern returns to the calling routine.
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
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN plotc1
CPS 
CPS IF a trc file exists to write to for this run
CPS 
CPS   IF we are writing only the node information and the number of 
CPS   solutes at the top of the file
CPS     Write number of nodes being output
CPS       FOR each output node
CPS           write node number and coordinates
CPS       END FOR
CPS     Write number of solutes
CPS   ELSE we are writing concentrations
CPS     IF this is not a Henry's Law, vapor specified species
CPS       Write the time, number of the solute
CPS       Write the concentration at each specified node
CPS     ELSE we must calculate the vapor concentration from the liquid
CPS       Write the time, number of the solute
CPS       Write the concentration at each specified node
CPS     ENDIF
CPS   ENDIF
CPS 
CPS ENDIF
CPS END plotc1
CPS 
C**********************************************************************
c****--------------------------------------------------------------****c
c**** write to plot tape                                           ****c

      use combi
      use comchem
      use comci
      use comdi
      use comrxni
      use compart
      use comdti
      use comai
      use comsptr
      implicit none

      integer ifg,ii
      integer i, mi, iunit
      integer ic,im,iv,ix
      real*8, allocatable :: complex_temp(:)
      real*8 complex_conc, ptime
      character*80 form_string
      logical :: tracout_flag = .false.

      save tracout_flag, form_string

      if ((istrc.gt.0).and.(pout.ge.0))  then

         if (time_flag .eq. 1) then
            ptime = abs(days / 365.25d00)
         else if (time_flag .eq. 3) then
            ptime = abs(days * 86400.d00)
         else if (time_flag .eq. 4) then
            ptime = abs(days * 24.d00)
         else
            ptime = abs(days)
         end if

         if( ifg .eq. 0 )  then
c**** write number of plot nodes, node numbers and coordinates ****

            write (istrc, *)  m
            do i = 1, m
               mi = nskw(i)
               write(istrc, 6010) mi, cord(mi,1), cord(mi,2), cord(mi,3)
 6010          format(i8, 3e16.9)
            end do

            write(istrc ,   *)  nspeci,ncpnt,nimm,nvap,ncplx

            if (ishisc .ne. 0) then
               tracout_flag = .true.
 
               call tracout
               select case (form_flag)
               case (1)
c Tecplot
                  write (form_string, 300) m
               case (2)
c Surfer (csv)
                  write (form_string, 310) m
               case default
c Standard text
                  write (form_string, 300) m
               end select
 300     format ("(g16.9, ", i5, "(1x, g16.9))")
 310     format ("(g16.9, ", i5, '(", ", g16.9))')               
            end if
            call flush(istrc)
         else if (ifg.eq.1) then
            if (tracout_flag) then
c Shouldn't get here not applicable to particle tracking              
            else
               if ( icns(nsp) .gt. -1 ) then
                  write(istrc ,   *)  days, nsp
                  write(istrc ,   *)  ( an(nskw(i)+npn) , i=1,m )
               else
                  write(istrc ,   *)  days, nsp
                  write(istrc ,   *)  ( anv(nskw(i)+npn) , i=1,m )
               endif
            end if
         else if (ifg.eq.2) then
! Aqueous species
            ic=ii
            if (tracout_flag) then
               iunit = ishisc + nsp
               if ( icns(nsp) .gt. -1 ) then
                  write(iunit, form_string) ptime, 
     &                 (an(nskw(i)+npn), i = 1,m )
               else
                  write(iunit, form_string) ptime, 
     &                 (anv(nskw(i)+npn), i = 1,m )
               end if
               call flush(iunit)
            else
               if ( icns(nsp) .gt. -1 ) then
                  write(istrc ,   *)  days, nsp,'  ', cpntnam(ic)
                  write(istrc ,   *)  ( an(nskw(i)+npn) , i=1,m )
               else
                  write(istrc ,   *)  days, nsp,'  ', cpntnam(ic)
                  write(istrc ,   *)  ( anv(nskw(i)+npn) , i=1,m )
               endif
            end if
         else if (ifg.eq.3) then
! Immobile species
            im=ii
            if (tracout_flag) then
               iunit = ishisc + nsp
               write(iunit, form_string) ptime, (an(nskw(i)+npn),
     &              i = 1,m )
               call flush(iunit)
            else
               write(istrc ,   *)  days , nsp,'  ', immnam(im)
               write(istrc ,   *)(an(nskw(i)+npn), i = 1,m )
            end if
         else if (ifg.eq.4) then
! Vapor species
            iv=ii
            if (tracout_flag) then
               iunit = ishisc + nsp
               write(iunit, form_string) ptime, (an(nskw(i)+npn),
     &              i = 1,m )               
               call flush(iunit)
            else
               write(istrc ,   *)  days , nsp,'  ', vapnam(iv)
               write(istrc ,   *)(an(nskw(i)+npn), i = 1,m )
            end if
         else if (ifg.eq.5) then
! Free Ion
            ic=ii
            if (tracout_flag) then
               iunit = ishisc + ic + nspeci
               write(iunit, form_string) ptime, (cpntsv(ic,nskw(i)), 
     &              i = 1,m ) 
               call flush(iunit)
            else
               write(istrc, *) days, ic, ' Free Ion ', cpntnam(ic)
               write(istrc ,   *)(cpntsv(ic,nskw(i)), i = 1,m )
            end if       
         else if (ifg.eq.6) then
! Aqueous complex
            ix=ii
            if (.not. allocated(complex_temp)) then
               allocate (complex_temp(m))
            end if
            do i = 1,m
               call cplxcalc(nskw(i),ix,complex_conc)
               complex_temp(i)=complex_conc
            enddo
            if (tracout_flag) then
               iunit = ishisc + ix
               write(iunit, form_string) ptime, (complex_temp(i),i=1,m)
               call flush(iunit)
            else
               write(istrc, *)days, ix,' Complex ', cplxnam(ix)
               write(istrc, *)(complex_temp(i),i=1,m)
            end if
         endif
      endif
      if (.not. tracout_flag) call flush(istrc)
      return
      end subroutine plotc1

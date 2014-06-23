      subroutine  diskc(read_trac)
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
CD1 To read or write restart files for tracer variables.
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
CD2 $Log:   /pvcs.config/fehm90/src/diskc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:22   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:34   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Feb 15 10:03:44 1996   zvd
CD2 Added requirement and updated purpose.
CD2 
CD2    Rev 1.3   Mon Jan 29 15:01:46 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   03/15/95 17:04:44   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.1   03/18/94 16:15:36   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:02   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 None
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 Name                  Description
CD3 
CD3 File number iread     Input restart file
CD3 Tape number iout      Tape number for error messages
CD3 File number iatty     File number for error messages
CD3 File number isave     Output restart file
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
CD4 nsave, iread, iout, iatty, npn, npt, n an, anv, nspeci, an0, isave
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 Global Subprograms
CD4
CD4 None
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
CD5 isolflg      int         Flag used to determine if data are
CD5                             present in file
CD5 dumm         real*8      Dummy variable used in read
CD5 nspecl       int         Number of species written to output file
CD5 iq           int         Do loop index over each species
CD5 mi           int         Loop index over each node
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
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
CD9 2.7 Provide Restart Capability
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
CPS BEGIN diskc
CPS 
CPS IF this call is for reading solute data from file
CPS 
CPS   Read in dummy value to see if solute data exits
CPS   
CPS   IF solute data does not exist
CPS     Write error message
CPS   ELSE solute data exits
CPS     Read number of species written to output file
CPS     
CPS     FOR each species written to output file
CPS       Read the concentrations
CPS       Write information to check file
CPS     ENDFOR
CPS     
CPS     FOR each new species
CPS       Write information to check file
CPS     ENDFOR
CPS     
CPS   ENDIF
CPS 
CPS ELSE this call is for writing data to output file
CPS   IF writing this file is specified
CPS     Write number of species to file
CPS     
CPS     FOR each species
CPS     
CPS       IF this is a Henry's Law species with vapor concentrations...
CPS       ... written to the file
CPS         Write vapor concentrations
CPS       ELSE
CPS         Write liquid concentrations
CPS       ENDIF
CPS       
CPS     ENDFOR
CPS     
CPS   ENDIF
CPS ENDIF
CPS   
CPS END diskc
CPS 
C***********************************************************************

      use comai
      use combi
      use comchem, only : conc_read
      use comdi
      use comdti
      use comrxni
      implicit none

      integer isolflg,nspecl,iq,mi
      real*8 dumm
      character*4 tracword
      logical read_trac

      nspecl = 0
      if ( nsave .eq. 0 )  then
         isolflg = 0
         if (wdd1(5:8) .eq. 'trac') then
            if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
               read(iread, end = 10, err = 10) nspecl
            else
               read(iread ,   *, end=10, err = 10) nspecl
            end if
            isolflg = 1
         else if (read_trac) then
            if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
               read(iread, end = 10, err = 10) tracword
               read(iread, end = 10, err = 10) nspecl
            else
               read(iread ,   *, end=10, err = 10) tracword
               read(iread ,   *, end=10, err = 10) nspecl
            end if
            if (tracword .eq. 'trac') isolflg = 1
         end if
 10      continue
         if( isolflg .eq. 0) then
c     
c     either eof found or particle tracking info found (no
c     concentration data to read in
c     
            write(ierr, 6000)
            if (iout .ne. 0) write(iout, 6000)
            if (iptty .gt. 0 ) write(iptty, 6000)
 6000       format('Tracer data not found in restart file')
         else if (.not. read_trac) then
            write(ierr, 6001)
            if (iout .ne. 0) write(iout, 6001)
            if (iptty .gt. 0 ) write(iptty, 6001)
 6001       format('Tracer data found in restart file will not be used')
         else
            do iq=1,nspecl
               npn=npt(iq)
               if (bin_flag .eq. 1 .or. bin_flag .eq. 2) then
                  read (iread)  ( an(mi+npn) , mi=1,n )
               else
                  read (iread ,   *)  ( an(mi+npn) , mi=1,n )
               end if
               if (icns(iq) .eq. -2) then
                  anv(1+npn:n+npn) = an(1+npn:n+npn)
               else
                  anl(1+npn:n+npn) = an(1+npn:n+npn)
               end if
               if (ischk .ne. 0) then
                  write(ischk, *)
                  write(ischk, *) 'Concentrations for species ', iq,
     2                 ' read from restart file'
                  write(ischk, *)
               end if
               conc_read(iq) = .true.
            enddo
            do iq=nspecl+1,nspeci
               if (ischk .ne. 0) then
                  write(ischk, *)
                  write(ischk, *) 'Concentrations for species ', iq,
     2                 ' not present in restart file'
                  write(ischk, *) 'They are set to input values from', 
     2                 ' macro trac instead'
                  write(ischk, *)
               end if
            enddo
         endif
      else
         if ( isave .gt. 0 )  then
            if (bin_flag .eq. 1 .or. bin_flag .eq. 3) then
               if (header_flag .eq. 'new') write(isave) 'trac'
               write(isave) nspeci
               do iq=1,nspeci
                  npn    =  npt(iq)
                  if( icns(iq) .eq. -2 ) then
                     write(isave) (anv(mi+npn), mi = 1,n)
                  else
                     write(isave) (an(mi+npn), mi = 1,n)
                  end if
               enddo
            else
               if (header_flag .eq. 'new') write(isave, *) 'trac'
               write(isave, *) nspeci
               do iq=1,nspeci
                  npn    =  npt(iq)
                  if( icns(iq) .eq. -2 ) then
                     write(isave , 6100) (anv(mi+npn), mi = 1,n)
                  else
                     write(isave , 6100) (an(mi+npn), mi = 1,n)
                  end if
 6100             format(4g20.10)
               enddo
            endif
         endif
      endif
      
      return
      end
      

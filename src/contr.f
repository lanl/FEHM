      subroutine contr(inj)
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
!***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Write to contour plot tape.  
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-FEB-94    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/contr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:12 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Mon Jan 29 13:58:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.7   03/10/95 10:31:18   llt
CD2 changed to calculate velocity with 1 call - gaz
CD2 
CD2    Rev 1.6   12/01/94 08:45:32   llt
CD2 Changed max functions so would run on ibm.
CD2 
CD2    Rev 1.5   08/23/94 08:48:58   llt
CD2 Accidently lost revision 1.3. Added changes back in.
CD2 
CD2    Rev 1.4   08/23/94 08:24:12   llt
CD2 gaz changes
CD2 
CD2    Rev 1.3   06/20/94 10:55:08   zvd
CD2 Modified to remove contour output files if AVS output is used.
CD2
CD2    Rev 1.2   03/18/94 15:42:54   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/14/94 11:43:36   zvd
CD2 Modified so if contour files not present code won't attempt write.
CD2
CD2    Rev 1.0   01/20/94 10:22:30   pvcs
CD2 original version in process of being certified
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   inj             INT      I    Flag to indicate if this is first call to
CD3                                   routine
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   iscon                    O    Contour plot data file
CD3   iscon1                   O    Contour plot data file for dual or dpdp
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
CD4   altc            CHAR     faax   String used for denoting type of output
CD4                                     for contour plots
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   cpr             REAL*8   fdd    Rock specific heat at each node
CD4   icap            INT      fddi1  Capillary pressure model at each node
CD4   iccen           INT      faai   Parameter which indicates if the tracer
CD4                                     solution is enabled
CD4   ipqz            INT      faai   The number of nodes used for the contour
CD4                                     plot output
CD4   irlp            INT      fddi1  Relative permeability model at each node
CD4   iscon           INT      faai   Unit number for contour data file
CD4   iscon1          INT      faai   Unit number for dual porosity or dpdp
CD4                                     contour data file
CD4   istrs           INT      faai   Parameter indicating if the stress
CD4                                     solution is enabled
CD4   jdate           CHAR     faac1  Contains the date (mm/dd/yr)
CD4   jtime           CHAR     faac1  Contains the time (hr:mn:sc)
CD4   nelm            INT      fbb    Initially information about nodes in each
CD4                                     element, later nodal connectivity
CD4                                     information
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4   nmfil           CHAR ARY faax   I/O file names:
CD4                                   10 - Contour plot output file name
CD4                                   11 - Dual porosity or dpdp contour plot
CD4                                          output file name
CD4   nspeci          INT      fdd1i  Number of species for tracer solution
CD4   pcp             REAL*8   fdd    Capillary pressure at each node
CD4   phi             REAL*8   fdd    Pressure at each node
CD4   pnx             REAL*8   fdd    Permeability in the x-direction, liquid
CD4                                     velocity in the x-direction, vapor
CD4                                     velocity in the x-direction
CD4   pny             REAL*8   fdd    Permeability in the y-direction liquid
CD4                                     velocity in the y direction, vapor
CD4                                     velocity in the y-direction
CD4   pnz             REAL*8   fdd    Permeability in the z-direction, liquid
CD4                                     velocity in the z-direction, vapor
CD4                                     velocity in the z-direction
CD4   ps              REAL*8   fdd    Porosity at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   t               REAL*8   fdd    Temperature at each node
CD4   thx             REAL*8   fdd    Thermal conductivity x-direction
CD4   thy             REAL*8   fdd    Thermal conductivity y-direction
CD4   thz             REAL*8   fdd    Thermal conductivity z-direction
CD4   verno           CHAR     faac   Contains version number of FEHMN code
CD4   wdd             CHAR     faac   Problem title
CD4  
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   concen                   Control tracer simulation output
CD4   veloc                    Calculate fluid velocities
CD4   
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   Identifier      Type     Description
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
CD5   j               INT      Loop index
CD5   neq1            INT      Loop index for dual/dpdp
CD5   neq2            INT      Loop index for dual/dpdp
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
CPS BEGIN contr
CPS 
CPS   call subroutine veloc for phase velocities
CPS 
CPS   IF avs output is required
CPS      IF this is the first call
CPS         IF contour file exists
CPS            close contour file
CPS         END IF
CPS         IF dual/dpdp contour file exists
CPS            close dual/dpdp contour file
CPS         END IF
CPS      END IF
CPS      call avs_io to write output
CPS   ELSE IF ment output is required
CPS      [doesn't do anything]
CPS   ELSE IF ptrn output is required
CPS   
CPS      IF this is a liquid phase problem and the contour file exists
CPS         call contrj to write output
CPS      END IF
CPS      
CPS   ELSE IF fehm binary or free format output is wanted
CPS   
CPS      IF contour file exists
CPS         IF this is initial call
CPS            IF file is binary format
CPS               close and reopen contour file in unformatted mode
CPS               ****** unformatted writes ******
CPS               write program version, date, time, and problem title
CPS               IF tracer problem
CPS                  write 'trac' to file
CPS               ELSE
CPS                  write blank line
CPS               END IF
CPS
CPS               IF stress problem
CPS                  write 'strs' to file
CPS               ELSE
CPS                  write blank line
CPS               END IF
CPS               
CPS               write node coordinates and connectivity data
CPS               write permeabilities, conductivities, porosities, . . .
CPS               . . . specific heats, capillary pressures, relative . . .
CPS               . . . permeability model, and capillary pressure model
CPS               IF tracers are being used
CPS                  write number of species being used
CPS               END IF
CPS            ELSE
CPS               ****** free format writes ******
CPS               same as for unformatted writes (see above) except . . .
CPS               . . . program version, date, time, and problem title . . .
CPS               . . . have already been written to file
CPS            END IF
CPS         ELSE
CPS         
CPS            call veloc to compute velocities
CPS            
CPS            IF file is binary format
CPS               ****** unformatted writes ******
CPS               write simulation time and problem phase
CPS               IF all liquid
CPS                  write velocities, pressures, and temperatures
CPS               ELSE
CPS                  write velocities, pressures, and saturations
CPS               END IF
CPS            ELSE
CPS               ****** free format writes ******
CPS               same as unformatted writes (see above)
CPS            END IF
CPS            
CPS         END IF
CPS         
CPS      END IF
CPS      
CPS      IF dual/dpdp contour file exists
CPS         same as for contour file (see above), except for . . .
CPS         . . . dual posorsity / dual permeabilty nodes
CPS      END IF
CPS      
CPS      IF this isn't the initial call
CPS         call concen to write contour information for tracers
CPS      END IF
CPS   
CPS   END IF
CPS   
CPS END contr
CPS
C***********************************************************************

      use avsio
      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comflow, only : flag_heat_out
      use comhi
      use comii
      use compart
      use comwt
      use comxi
      use comsi
      use davidi

      implicit none
 
      real*8 px, py, pz
      real*8 anld,anvd
      real*8 phihd, sdum, hddum, phipcp
      real*8 dumconv, dumconv1, pdum, tdum, rolconv

      integer i,inj,j,neq1,neq2, dummyint,k,jp,kk,ll,length
      integer ialiq, iavap, ihead_sv
      
      if(icont.eq.0) return
c check for contour file output form
c call veloc calculations here
c both vapor and liquid velocities calculated with one call
c gaz 020122 call to vtk conversion      
      if(inj.eq.2) then
       call FEHM_tec_to_vtk(1)
       return
      endif  
c gaz 021523 testing flow rate vectors
c      if(inj.gt.0.and.compute_flow) call veloc
       if(inj.gt.0.and.compute_flow) then
        call flowrate_vectors(-1)
        call flowrate_vectors(1)
        call flowrate_vectors(-2)
       endif

c add call to heatloc to calculate heat fluxes for output
      if (inj.gt.0.and.flag_heat_out) call heatloc

c if head data requested call airctr(8,0)
c
      if(ichead.ne.0) then
c gaz 083119 don't turn head on exept to printout  
         ihead_sv = 1 
         if(ihead.eq.0) then
          ihead_sv = 0
          ihead = 1
         endif
         dumconv = crl(1,1)
         dumconv1 = crl(4,1)
         pdum = pres0+rol0*head0*(-grav)
         tdum = temp0        
         call water_density(tdum,pdum,rolconv)
         crl(1,1)=rolconv
         crl(4,1)=pres0
         call headctr(2,0,0.,0.)           
         crl(1,1)= dumconv
         crl(4,1)= dumconv1 
         if(ihead_sv.eq.0) then
          ihead = 0
         endif          
      else if(ihead.ne.0) then
         if (inj.ne.0) then
            rolconv = crl(1,1)
c
c convert from pressure to head
c
            call headctr(2,0,0.0,0.0)
	    if(ifree.ne.0) then
               do i = 1,n
                  if(rlxyf(i).lt.sattol+rlptol) head(i) = head_id
               enddo
	    endif
         else
	    head = pho
         endif
      else
         iohead=0
      endif
      if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'sur'.or. 
     &     altc(1:3) .eq. 'tec') then
c use avs form        
         if (inj .eq. 0) then
            inquire (file = nmfil(10), exist = ex)
            if (ex.and.iscon.ne.0) then
               close (iscon, status = 'delete')
            end if
            inquire (file = nmfil(11), exist = ex)
            if (ex.and.iscon1.ne.0) then
               close (iscon1, status = 'delete')
            end if
         end if       
         call avs_io(inj)
!      else if (altc(1:3) .eq. 'sur'.or. altc(1:3) .eq. 'tec') then
!         call surf_tec(inj)
      else if (altc .eq. 'ment') then
c     use mentat form
      else if (altc .eq. 'ptrn') then
c     use patran form
c     call only if liquid phase and contour plot file wanted
         if (inj .lt. 0 .and. iscon .gt. 0) call contrj(inj)
      else 
c use fehm binary or free format output
c**** 2-d, 3-d output ****
         if (iscon .gt. 0) then
            if (inj .eq. 0) then
               if(altc .eq. 'fehm') then
c for unformatted file close contour file and open again as an
c unformatted file
                  close (iscon,status='delete')
                  open(iscon, file = nmfil(10), status='unknown',
     *                 form = 'unformatted')
                  write(iscon) verno, jdate, jtime
                  write(iscon) wdd
                  ipqz = neq_primary
                  if (iccen .ne. 0) then
                     write(iscon)  'trac'
                  else
                     write(iscon)  '    '
                  end if
                  if (istrs .ne. 0)  then
                     write(iscon)  'strs'
                  else
                     write(iscon)  '    '
                  end if
                  write(iscon) neq_primary
                  write(iscon)  (corz(i, 1), corz(i, 2), corz(i,3),
     *                 i = 1, neq_primary)
                  write(iscon)  ns , nei
                  write(iscon) ((nelm((j-1)*ns+i), i =1, ns), 
     *                 j = 1, nei)
c**** selected input quantities ****
                  write(iscon)  (pnx(i)*1.0e-06,
     *                 pny(i)*1.0e-06, pnz(i)*1.0e-06, i=1,neq_primary)
                  if(ico2.ge.0) then
                     write(iscon )
     *                    ( thx(i) , thy(i) , thz(i) , i=1,neq_primary)
                  else
                     write(iscon )
     *                    ( 0.0d00 , 0.0d00 , 0.0d00 , i=1,neq_primary)
                  endif
                  if(ico2.ge.0 .and. irdof .ne. 13) then
                     write(iscon )
     *                    ( ps(i) , cpr(i) , pcp(i) , i=1,neq_primary )
                  else
                     write(iscon )
     *                    ( ps(i) , cpr(i) , 0.0d00 , i=1,neq_primary )
                  endif
                  write(iscon ) idof,igrav,grav
                  if ( iccen .ne. 0 )  then
                     write(iscon )  nspeci
                  end if
               else
                  ipqz   =  neq_primary
                  if (iccen .ne. 0)  then
                     write(iscon, *)  'trac'
                  else
                     write(iscon, *)  '    '
                  end if
                  if (istrs .ne. 0)  then
                     write(iscon, *)  'strs'
                  else
                     write(iscon, *)  '    '
                  end if
                  write(iscon, *)  neq_primary
                  write(iscon, 10) (corz(i, 1), corz(i, 2), corz(i, 3), 
     *                 i = 1, neq_primary)
 10               format(1x,3g20.8)
                  write(iscon, *)  ns , nei
                  write(iscon, *)  ((nelm((j-1)*ns+i), i = 1, ns), 
     *                 j = 1, nei)

c**** selected input quantities ****

                  write(iscon, 10)  (pnx(i)*1.0e-06, pny(i)*1.0e-06, 
     *                 pnz(i)*1.0e-06, i = 1, neq_primary)
                  if(ico2.ge.0) then
                     write(iscon,10)
     *                    ( thx(i) , thy(i) , thz(i) , i=1,neq_primary)
                  else
                     write(iscon,10)
     *                    ( 0.0d00 , 0.0d00 , 0.0d00 , i=1,neq_primary)
                  endif
                  if(ico2.ge.0 .and. irdof .ne. 13) then
                     write(iscon ,10)
     *                    (ps(i), cpr(i), pcp(i), i = 1, neq_primary)
                  else
                     write(iscon ,10)
     *                    (ps(i), cpr(i), 0.0d00, i = 1, neq_primary)
                  endif
                  write(iscon ,*) idof,igrav,grav
                  if (iccen .ne. 0)  then
                     write(iscon, *)  nspeci
                  end if
               end if
            else
               if (altc .eq. 'fehm') then
                  write(iscon )  days , inj
                  if (inj .gt. 0)  then
                     npn = n
                     write(iscon) (dil(i)/
     *                    max(rolf(i),0.000000001d+1),
     *                    rolf(i),  phi(i)-pcp(i), 
     *                    t(i), i = 1, neq_primary)
                  else
                     npn = n + n
                     write(iscon) (div(i)/
     *                    max(rovf(i),0.000000001d+1),
     *                    rovf(i),  phi(i), 
     *                    t(i), i = 1, neq_primary)
                  end if
               else
                  write(iscon, *)  days, inj
                  if (inj .gt. 0)  then
                     npn = n
                     write(iscon, 10) (dil(i)/
     *                    max(rolf(i),0.000000001d+1),
     *                    rolf(i),  phi(i)-pcp(i), 
     *                    t(i), i = 1, neq_primary)
                  else
                     npn    =  n+n
                     write(iscon, 10) (div(i)/
     *                    max(rovf(i),0.000000001d+1),
     *                    rovf(i),  phi(i), 
     *                    t(i), i = 1, neq_primary)
                  end if
               end if
            end if
         end if

c write out info for dp solution
         if((idualp .ne. 0 .or. idpdp .ne. 0) 
     *        .and. iscon1 .gt. 0) then

            neq1 = neq_primary + 1
            neq2 = neq_primary + neq_primary
            if (inj .eq. 0) then
               if(altc .eq. 'fehm') then
c for unformatted file close contour file and open again as an
c unformatted file
                  close (iscon1, status = 'delete')
                  open(iscon1, file = nmfil(11), status='unknown',
     *                 form = 'unformatted')
                  write(iscon1) verno, jdate, jtime
                  write(iscon1) wdd

                  if (iccen .ne. 0) then
                     write(iscon1)  'trac'
                  else
                     write(iscon1)  '    '
                  end if

                  if (istrs .ne. 0) then
                     write(iscon1)  'strs'
                  else
                     write(iscon1)  '    '
                  end if
                  write(iscon1) neq
                  write(iscon1) (corz(i,1), corz(i,2),
     *                 corz(i,3), i = 1, neq_primary)
                  write(iscon1) ns, nei
                  write(iscon1) ((nelm((j-1)*ns+i), i = 1, ns), 
     *                 j = 1, nei)

c**** selected input quantities ****
                  write(iscon1)  (pnx(i)*1.0e-06, pny(i)*1.0e-06, 
     *                 pnz(i)*1.0e-06, i = neq1, neq2)
                  if(ico2.ge.0) then
                     write(iscon1)  (thx(i), thy(i), thz(i), 
     *                    i = neq1, neq2)
                  else
                     write(iscon1)  (0.0d00, 0.0d00, 0.0d00, 
     *                    i = neq1, neq2)
                  endif
                  if(ico2.ge.0 .and. irdof .ne. 13) then
                     write(iscon1)  (ps(i), cpr(i), pcp(i), 
     &                    i = neq1, neq2)
                  else
                     write(iscon1)  (ps(i), cpr(i), 0.d0, 
     &                    i = neq1, neq2)
                  end if
                  write(iscon ) idof,igrav,grav
                  if (iccen .ne. 0) then
                     write(iscon1)  nspeci
                  end if
               else
                  if (iccen .ne. 0) then
                     write(iscon1, *) 'trac'
                  else
                     write(iscon1, *) '    '
                  end if
                  if (istrs .ne. 0) then
                     write(iscon1, *) 'strs'
                  else
                     write(iscon1, *) '    '
                  end if
                  write(iscon1, *) neq_primary
                  write(iscon1, 10) (corz(i,1), corz(i,2),
     *                 corz(i,3), i = 1, neq_primary)
                  write(iscon1, *) ns, nei
                  write(iscon1, *) ((nelm((j-1)*ns+i), i = 1, ns), 
     *                 j = 1, nei)

c**** selected input quantities ****
                  write(iscon1, 10) (pnx(i)*1.0e-06, pny(i)*1.0e-06, 
     *                 pnz(i)*1.0e-06, i = neq1, neq2 )
                  if(ico2.ge.0) then
                     write(iscon1,10)  (thx(i), thy(i), thz(i), 
     *                    i = neq1, neq2)
                  else
                     write(iscon1,10)  (0.0d00, 0.0d00, 0.0d00, 
     *                    i = neq1, neq2)
                  endif
                  if(ico2.ge.0 .and. irdof .ne. 13) then
                     write(iscon1, 10) (ps(i), cpr(i), pcp(i), 
     *                    i = neq1, neq2)
                  else
                     write(iscon1, 10) (ps(i), cpr(i), 0.d0, 
     *                    i = neq1, neq2)
                  end if
                  write(iscon ) idof,igrav,grav
                  if (iccen .ne. 0) then
                     write(iscon1, *) nspeci
                  end if
               end if
            else
               if(altc .eq. 'fehm') then
                  write(iscon1) days, inj
                  if (inj .gt. 0) then
                     npn    =  n
                     if (irdof .ne. 13) then
                        write(iscon1) (dil(i)/
     *                       max(rolf(i),0.000000001d+1),
     *                       rolf(i),  phi(i)-pcp(i), 
     *                       t(i), i = neq1, neq2)
                     else
                        write(iscon1) (dil(i)/
     *                       max(rolf(i),0.000000001d+1),
     *                       rolf(i),  phi(i), 
     *                       t(i), i = neq1, neq2)
                     end if
                  else
                     npn = n + n
                     write(iscon1) (div(i)/
     *                    max(rovf(i),0.000000001d+1),
     *                    rovf(i),  phi(i), 
     *                    t(i), i = neq1, neq2)
                  end if
               else
                  write(iscon1, *) days, inj
                  if (inj .gt. 0) then
                     npn = n
                     if (irdof .ne. 13) then
                        write(iscon1, *) (dil(i)/
     *                       max(rolf(i),0.000000001d+1),
     *                       rolf(i),  phi(i)-pcp(i), 
     *                       t(i), i = neq1, neq2)
                     else
                        write(iscon1, *) (dil(i)/
     *                       max(rolf(i),0.000000001d+1),
     *                       rolf(i),  phi(i), 
     *                       t(i), i = neq1, neq2)
                     end if
                  else
                     npn = n + n
                     write(iscon1, *) (div(i)/
     *                    max(rovf(i),0.000000001d+1),
     *                    rovf(i),  phi(i), 
     *                    t(i), i = neq1, neq2)
                  end if
               end if
            end if
         end if

         if (inj .ne. 0) then
            call  concen (5,0)
         end if
         if(ichead.ne.0) then
            ihead = 0
         endif
               
      endif
      end

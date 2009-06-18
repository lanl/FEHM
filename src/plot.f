      subroutine plot(igf, tmavg, pravg)
!***********************************************************************
!  Copyright, 1994, 2004,  The  Regents of the University of California.
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
CD1 Write to time history plot tape.   
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-FEB-94    Z. Dash        22      Add prolog / minor change.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/plot.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:38   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:12   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:56   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Feb 28 16:58:52 1996   robinson
CD2 Added keyword to the output to denote type of fluid flow solution
CD2 
CD2    Rev 1.5   Wed Feb 28 16:36:36 1996   robinson
CD2 Output now includes air par. press. and rel. humidity
CD2 
CD2    Rev 1.4   Thu Feb 01 14:11:58 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   06/03/94 15:36:18   gaz
CD2  
CD2 
CD2    Rev 1.2   03/18/94 15:43:12   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/14/94 11:44:46   zvd
CD2 Modified so if history plot file not present code won't attempt a write.
CD2
CD2    Rev 1.0   01/20/94 10:26:22   pvcs
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
CD3   igf             INT      I    Flag to indicate if this is first call to
CD3                                   routine
CD3   tmavg           REAL*8   ?    Not used                   
CD3   pravg           REAL*8   ?    Not used
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   ishis                    O    Plot history data file
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
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   iccen           INT      faai   Parameter which indicates if the tracer
CD4                                     solution is enabled
CD4   ishis           INT      faai   Unit number for history data file
CD4   istrs           INT      faai   Parameter indicating if the stress
CD4                                     solution is enabled
CD4   m               INT      faai   Total number of nodes used for output
CD4                                     information
CD4   nskw            INT      fddi   Contains nodes for print-out
CD4   pcp             REAL*8   fdd    Capillary pressure at each node
CD4   phi             REAL*8   fdd    Pressure at each node
CD4   qh              REAL*8   fdd    Energy source term at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   sk              REAL*8   fdd    Source strength of each node
CD4   t               REAL*8   fdd    Temperature at each node
CD4
CD4 Global Subprograms
CD4
CD4   None
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
CD5   i               INT      Loop index
CD5   mi              INT      Node number
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
CPS BEGIN plot
CPS
CPS   IF plot history file exists
CPS   
CPS      IF this is initial call to history data output routine
CPS         IF tracer problem
CPS            write 'trac' to file
CPS         ELSE
CPS            write blank line
CPS         END IF
CPS         
CPS         IF stress problem
CPS            write 'strs' to file
CPS         ELSE
CPS            write blank line
CPS         END IF
CPS         
CPS         write number of nodes being output
CPS         FOR each output node
CPS             write node number and coordinates
CPS         END FOR
CPS         
CPS         IF ngas problem
CPS           write output headings
CPS         ELSE
CPS           write output headings
CPS         END IF
CPS         
CPS      END IF
CPS      
CPS      write time
CPS      
CPS      IF ngas problem
CPS        FOR each output node
CPS            write node number and values of enthalpy, flow, . . .
CPS             . . . temperature, pressure, air pressure, . . .
CPS             . . . capillary pressure, saturation, and . . .
CPS             . . . relative humidity
CPS        END FOR
CPS      ELSE
CPS        FOR each output node
CPS            write node number and values of enthalpy, flow, . . .
CPS             . . . temperature, pressure, capillary pressure, . . .
CPS             . . . and saturation 
CPS        END FOR
CPS
CPS      END IF
CPS      
CPS   END IF
CPS   
CPS END plot
C***********************************************************************

      use combi
      use comdi
      use comfi
      use comii
      use comdti
      use comai
      use comco2
      use comwt, only : sattol, rlptol
      use davidi
      implicit none

      integer i,igf,mi 
      real*8 tmavg,pravg,relhum,pwatersat,psat,dpdummy
      real*8 qhdum, skdum, tdum, phidum, pcidum
      real*8 pcpdum, sdum, headdum
      real*8 dumconv, dumconv1, rolconv
      real*8 pdum

      if (ishis .gt. 0) then

         if (igf .eq. 0)  then

            if(ico2.gt.0) then
               write(ishis, '(a4)')  'ngas'
            elseif(ico2.lt.0) then
               write(ishis, '(a4)')  'airw'
            else
               write(ishis, '(a4)')  '    '
            end if

            if (iccen .ne. 0)  then
               write(ishis, '(a4)')  'trac'
            else
               write(ishis, '(a4)')  '    '
            end if

            if ( istrs .ne. 0 )  then
               write(ishis, '(a4)')  'strs'
            else
               write(ishis, '(a4)')  '    '
            end if

c**** write number of plot nodes, node numbers and coordinates ****

            write (ishis, *)  m
            do i = 1, m
               mi = nskw(i)
               if (mi .gt. neq*2) then
                  mi = mi - neq*2
               else if (mi .gt. neq) then
                  mi = mi - neq
               end if
               write(ishis, 6010) nskw(i), corz(mi,1), corz(mi,2),
     *              corz(mi,3)
 6010          format(i8, 3e16.9)
            end do
            if(ico2 .gt. 0) then
               write(ishis, 6016) 'headings', 'node', 
     *              'flow enthalpy(Mj/kg)', 'flow(kg/s)',
     *              'temperature(deg C)', 'total pressure(Mpa)',
     *              'air pressure(MPa)', 'capillary pressure(Mpa)',
     *              'saturation(kg/kg)', 'relative humidity'
            else
             if(ihead.ne.0.or.ichead.ne.0) then
               write(ishis, 6015) 'headings', 'node', 
     *              'flow enthalpy(Mj/kg)', 'flow(kg/s)',
     *              'temperature(deg C)', 'hydraulic head(m)',
     *              'total pressure(Mpa)', 'saturation(kg/kg)'
             else
               write(ishis, 6015) 'headings', 'node', 
     *              'flow enthalpy(Mj/kg)', 'flow(kg/s)',
     *              'temperature(deg C)', 'total pressure(Mpa)',
     *              'capillary pressure(Mpa)', 'saturation(kg/kg)'
             end if
            end if
 6015       format(a8,/,a4,1x,a20,1x,a10,1x,a18,1x,a19,/,a23,1x,a17)
 6016       format(a8,/,a4,1x,a20,1x,a10,1x,a18,1x,a19,/,a17,1x,
     *             a23,1x,a17, 1x, a17)

         end if

         write(ishis, *)  days
         if(ico2 .gt. 0) then
            do i = 1, m
               mi = nskw(i)
               pwatersat = psat(t(mi),dpdummy,0)
               relhum = max((phi(mi)-pci(mi))/pwatersat,1.d-20)
               qhdum=max(qh(mi),1.d-20)
               if(abs(sk(mi)).lt.1.d-20) then
                skdum=0.0
               else
                skdum=sk(mi)
               endif
               tdum=max(t(mi),1.d-20)
               phidum=max(phi(mi),1.d-20)
               pcidum=max(pci(mi),1.d-20)
               if(abs(pcp(mi)).lt.1.d-20) then
                pcpdum=0.0
               else
                pcpdum=pcp(mi)
               endif
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  sdum=max(s(mi),1.d-20)
               else
                  sdum = 1.0d0
               end if
               write(ishis, 6020)  mi, qhdum , skdum , tdum , phidum ,
     *              pcidum, pcpdum  ,sdum , relhum 
            end do
 6021          format(i8, 8(1x, g16.9))
         else
            do i = 1, m
               mi = nskw(i)
               qhdum=max(qh(mi),1.d-20)
               if(abs(sk(mi)).lt.1.d-20) then
                skdum=0.0
               else
                skdum=sk(mi)
               endif
               tdum=max(t(mi),1.d-20)
               phidum=max(phi(mi),1.d-20)
               if (irdof .ne. 13) then
                  if(abs(pcp(mi)).lt.1.d-20) then
                     pcpdum=0.0
                  else
                     pcpdum=pcp(mi)
                  endif
               else
                  pcpdum=0.0
               end if
               if (ifree .ne. 0) then
                  sdum=max(s(mi), rlptol)
                  if (sdum .le. rlptol) sdum = 0.d0
               else if (irdof .ne. 13) then
                  sdum = max (s(mi), 1.d-20)
               else
                  sdum = 1.0d0
               end if
               if(ihead.ne.0) then
                  call headctr(4,mi,phi(mi),headdum)
                  if(sdum .le. sattol) headdum=0.d0
                  headdum=max(headdum,0.0d00)
                  phidum = max(crl(4,1),phidum - crl(1,1)*head0*(-grav))
                  write(ishis, 6020)  mi, qhdum, skdum, tdum, 
     *                 headdum, phidum , sdum 
               else if(ichead.ne.0) then
                  ihead=1
                  dumconv = crl(1,1)
                  dumconv1 = crl(4,1)
                  pdum = pres0+rol0*head0*(-grav)
                  tdum = temp0
                  call water_density(tdum,pdum,rolconv)
                  crl(1,1)=rolconv
                  crl(4,1)=pdum
                  call headctr(4,mi   ,pho(mi   ),headdum)
                  crl(1,1)= dumconv
                  crl(4,1)= dumconv1
                  ihead=0
                  if(sdum.lt.1.0) headdum=0.0
                  headdum=max(headdum,0.0d00)
                  write(ishis, 6020)  mi, qhdum, skdum, tdum , 
     *                 headdum, phidum, sdum
               else
                  write(ishis, 6020)  mi, qhdum, skdum, tdum, 
     *                 phidum, pcpdum, sdum 
               endif
 6020          format(i8, 6(1x, g16.9))
            end do
         end if

         call flush (ishis)

      end if

      if (ice.ne.0) call plot_hyd(igf)

C RJP 04/10/07 added below
      if(icarb.eq.1) call plot_co2(igf)

      end

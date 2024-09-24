      subroutine  datchk
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Analyze input data and some generated quantities.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 23-NOV-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/datchk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:34   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:10   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:26 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Tue May 21 13:45:22 1996   gaz
CD2 changed loop parameter from neq to n near bottom
CD2 
CD2    Rev 1.8   Thu May 02 16:38:10 1996   gaz
CD2 set volume(i) array to sx1(i) for mdnodes
CD2 
CD2    Rev 1.7   Fri Apr 26 15:07:24 1996   gaz
CD2 set nodes to minimum volume if 0.0
CD2 
CD2    Rev 1.6   Mon Jan 29 14:19:22 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.5   Wed Jan 17 11:14:18 1996   zvd
CD2 Passed ierr value to check_sx from datchk
CD2 
CD2    Rev 1.4   09/08/95 11:49:00   zvd
CD2 Corrected write to iptty when unassigned.
CD2 
CD2    Rev 1.3   03/22/95 10:09:14   gaz
CD2 gaz: skipped over boundary cells(large volume) in prining max cell volumes
CD2 
CD2    Rev 1.2   01/28/95 13:54:16   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 15:53:22   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:50   pvcs
CD2 original version in process of being certified
CD2 
c 12/15/94
c added call checksum
c call reallocate 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   ischk                    O    File used for input data consistency check
CD3                                   output
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
CD4   cord            REAL*8   fbs    Coordinates of each node
CD4   cpr             REAL*8   fdd    Rock specific heat at each node
CD4   denr            REAL*8   fdd    Rock density at each node
CD4   icnl            INT      faai   Problem dimension
CD4   irdof           INT      david1 Reduced degree of freedom model used
CD4   ischk           INT      faai   Unit number for input data check
CD4   n               INT      faai   Total number of nodes
CD4   nmfil           CHAR ARY faax   I/O file names:
CD4                                   13 - Input check output file name
CD4   pci             REAL*8   co2    Gas pressure
CD4   phi             REAL*8   fdd    Pressure at each node
CD4   pnx             REAL*8   fdd    Permeability in the x-direction, liquid
CD4   pny             REAL*8   fdd    Permeability in the y-direction liquid
CD4   pnz             REAL*8   fdd    Permeability in the z-direction, liquid
CD4   ps              REAL*8   fdd    Porosity at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   sx1             REAL*8   fbc    Volume associated with each node
CD4   t               REAL*8   fdd    Temperature at each node
CD4   thx             REAL*8   fdd    Thermal conductivity x-direction
CD4   thy             REAL*8   fdd    Thermal conductivity y-direction
CD4   thz             REAL*8   fdd    Thermal conductivity z-direction
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4   
CD4   check_sx        N/A      Check volumes and finite element flow coefficients
CD4   min_max         N/A      Find the minimum and maximum parameter values 
CD4                              and their location
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
CD5   cprmmn          REAL*8   Minimum value of rock specific heat 
CD5   cprmmx          REAL*8   Maximum value of rock specific heat 
CD5   crxmax          REAL*8   Maximum x coordinate value 
CD5   crxmin          REAL*8   Minimum x coordinate value 
CD5   crymax          REAL*8   Maximum y coordinate value 
CD5   crymin          REAL*8   Minimum y coordinate value 
CD5   crzmax          REAL*8   Maximum z coordinate value 
CD5   crzmin          REAL*8   Minimum z coordinate value 
CD5   denmmn          REAL*8   Minimum value of rock density 
CD5   denmmx          REAL*8   Maximum value of rock density
CD5   dfctor          REAL*8   Factor for setting z coordinate to 0 in 2d
CD5                              problem
CD5   icprmn          INT      Location of minimum value of rock specific heat
CD5   icprmx          INT      Location of maximum value of rock specific heat 
CD5   icrxmn          INT      Location of minimum x coordinate
CD5   icrxmx          INT      Location of maximum x coordinate
CD5   icrymn          INT      Location of minimum y coordinate
CD5   icrymx          INT      Location of maximum y coordinate
CD5   icrzmn          INT      Location of minimum z coordinate
CD5   icrzmx          INT      Location of maximum z coordinate
CD5   idenmn          INT      Location of minimum value of rock density
CD5   idenmx          INT      Location of maximum value of rock density
CD5   ipermn          INT      Location of minimum value of permeability
CD5   ipermx          INT      Location of maximum value of permeability
CD5   ipormn          INT      Location of minimum value of porosity
CD5   ipormx          INT      Location of maximum value of porosity 
CD5   iprcmn          INT      Location of minimum value of gas pressure
CD5   iprcmx          INT      Location of maximum value of gas pressure
CD5   iprsmn          INT      Location of minimum value of pressure 
CD5   iprsmx          INT      Location of maximum value of pressure  
CD5   isatmn          INT      Location of minimum value of liquid saturation
CD5   isatmx          INT      Location of maximum value of liquid saturation
CD5   itemmn          INT      Location of minimum value of temperature 
CD5   itemmx          INT      Location of maximum value of temperature 
CD5   ithcmn          INT      Location of minimum value of thermal conductivity
CD5   ithcmx          INT      Location of maximum value of thermal conductivity
CD5   ivolmn          INT      Location of minimum value of volume
CD5   ivolmx          INT      Location of maximum value of volume
CD5   min_nr          INT      Minimum value of storage for fe coefficients
CD5   perimn          REAL*8   Minimum value of permeability at single node
CD5   permmn          REAL*8   Minimum value of permeability  
CD5   perimx          REAL*8   Maximum value of permeability at single node 
CD5   permmx          REAL*8   Maximum value of permeability  
CD5   pormax          REAL*8   Maximum value of porosity
CD5   pormin          REAL*8   Minimum value of porosity 
CD5   prcmax          REAL*8   Maximum value of gas pressure 
CD5   prcmin          REAL*8   Minimum value of gas pressure
CD5   prsmax          REAL*8   Maximum value of pressure  
CD5   prsmin          REAL*8   Minimum value of pressure  
CD5   satmax          REAL*8   Maximum value of liquid saturation
CD5   satmin          REAL*8   Minimum value of liquid saturation 
CD5   temmax          REAL*8   Maximum value of temperature
CD5   temmin          REAL*8   Minimum value of temperature
CD5   thcimn          REAL*8   Minimum value of thermal conductivity at single
CD5                              node
CD5   thcmmn          REAL*8   Minimum value of thermal conductivity 
CD5   thcimx          REAL*8   Maximum value of thermal conductivity at single
CD5                              node 
CD5   thcmmx          REAL*8   Maximum value of thermal conductivity  
CD5   volmax          REAL*8   Maximum value of volume
CD5   volmin          REAL*8   Minimum value of volume 
CD5   voltot          REAL*8   Total volume
CD5   zdiff           REAL*8   Distance between z coordinates
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
CD8 None
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9 2.8 Provide Multiple Realization Option
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
CPS BEGIN  datchk
CPS 
CPS   initialize minimum and maximum check parameter values
CPS   
CPS   FOR each node
CPS       find location and value of parameter minima and maxima 
CPS       [parameters being checked include: pressure, temperature, 
CPS        saturation, porosity, permeability, rock density, 
CPS        specific heat, thermal conductivity, coordinate values 
CPS        and volumes]
CPS   END FOR
CPS   
CPS   IF condition is satisfied [tracer, stress, noncondesible gas, dual 
CPS    porosity, ice, wellbore, heat or heat and mass, quadrature, 
CPS    upwind weight, porosity, reduced degree of freedom, finite 
CPS    volumes]
CPS      write miscellaneous problem information to check file
CPS   END IF
CPS   
CPS   write minimum/maximum values and locations to input check file
CPS   
CPS END  datchk
CPS
C***********************************************************************

      use davidi
      use comfi
      use comci, only : dstm
      use comdi
      use combi
      use comdti
      use comai
      use comxi
      use compart
      use comii
      implicit none

      real*8  cprmmn, cprmmx, crxmax, crxmin, crymax, crymin, crzmax,
     .	      crzmin, denmmn, denmmx, dfctor, perimn, perimx, permmn, 
     .	      permmx, pormax, pormin, prcmax, prcmin, prsmax, prsmin, 
     .	      satmax, satmin, temmax, temmin, thcimn, thcimx, thcmmn, 
     .          thcmmx, volmax, volmin, voltot, zdiff, vol_fac, den_vap,
     .          perm_fac
      
      integer i, icprmn, icprmx, icrxmn, icrxmx, icrymn, icrymx, icrzmn,
     .	      icrzmx, idenmn, idenmx, ipermn, ipermx, ipormn, ipormx,
     .	      iprcmn, iprcmx, iprsmn, iprsmx, isatmn, isatmx, itemmn,
     .	      itemmx, ithcmn, ithcmx, ivolmn, ivolmx, min_nr
      integer icsub, j, ivol_old
c gaz 022419        
      parameter (vol_fac = 1.d-8, perm_fac = 1.d3, ivol_old = 0)
      
      prsmin =  120.0
      prsmax = -120.0
      iprsmn =    1
      iprsmx =    1
      
      prcmin =  120.0
      prcmax = -120.0
      iprcmn =    1
      iprcmx =    1

      temmin =  370.0
      temmax = -300.0
      itemmn =    1
      itemmx =    1

      satmin =    2.0
      satmax =   -1.0
      isatmn =    1
      isatmx =    1

      pormin =    1.0
      pormax =    0.0
      ipormn =    1
      ipormx =    1

      permmn =    1.0d+20
      permmx =   -1.0d+20
      ipermn =    1
      ipermx =    1

      denmmn =    1.0d+20
      denmmx =   -1.0d+20
      idenmn =    1
      idenmx =    1

      cprmmn =    1.0d+20
      cprmmx =   -1.0d+20
      icprmn =    1
      icprmx =    1

      thcmmn =   1.d+20
      thcmmx =   -1.d+20
      ithcmn =    1
      ithcmx =    1

c**** find max and min of coordinate values and volumes ****

      crxmin =    1.0d+20
      crxmax =   -1.0d+20
      icrxmn =    1
      icrxmx =    1

      crymin =    1.0d+20
      crymax =   -1.0d+20
      icrymn =    1
      icrymx =    1

      crzmin =    1.0d+20
      crzmax =   -1.0d+20
      icrzmn =    1
      icrzmx =    1

      volmin =    1.0d+20
      volmax =   -1.0d+20
      ivolmn =    1
      ivolmx =    1
      voltot =    0.0

c**** search for maximum and minimum values ****

! If a check file is not requested skip checking (zvd 01/07/2005)
      if (ischk .ne. 0) then
            zdiff = 0.0d0
         do i = 1, n
            zdiff = max(zdiff,abs(cord(1, igrav) - cord(i, igrav)))
            call min_max (phi(i), i, prsmin, iprsmn, prsmax, iprsmx)
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               call min_max (pci(i), i, prcmin, iprcmn, prcmax, iprcmx)
            end if
            call min_max (t  (i), i, temmin, itemmn, temmax, itemmx)
c gaz 110120 
            if (irdof .ne. 13 .and. ifree .ne. 0.and.idoff.ne.-1) then
               call min_max (s  (i), i, satmin, isatmn, satmax, isatmx)
            else
               satmin = 1.0d0
               satmax = 1.0d0
            end if
            call min_max (ps (i), i, pormin, ipormn, pormax, ipormx)
            call min_max (denr(i), i, denmmn, idenmn, denmmx, idenmx)
            call min_max (cpr(i), i, cprmmn, icprmn, cprmmx, icprmx)

            if (icnl .eq. 0) then
               if (idoff .gt. 0) then
                  perimn = min(pnx(i), pny(i), pnz(i)) * 1.0d-06
                  perimx = max(pnx(i), pny(i), pnz(i)) * 1.0d-06 
               end if
               if(ico2.ge.0.or.ice.ne.0) then
                  thcimn = min(thx(i), thy(i), thz(i))
                  thcimx = max(thx(i), thy(i), thz(i))
               endif
            else  
               if (idoff .gt. 0) then
                  perimn = min(pnx(i), pny(i)) * 1.0d-06
                  perimx = max(pnx(i), pny(i)) * 1.0d-06
               end if
               if(ico2.ge.0) then
                  thcimn = min(thx(i), thy(i))
                  thcimx = max(thx(i), thy(i))
               endif
            end if

            call min_max (perimn, i, permmn, ipermn, permmx, ipermx)
            call min_max (perimx, i, permmn, ipermn, permmx, ipermx)
            call min_max (thcimn, i, thcmmn, ithcmn, thcmmx, ithcmx)
            call min_max (thcimx, i, thcmmn, ithcmn, thcmmx, ithcmx)
            if(i.le.neq) then
               call min_max (cord(i,1), i, crxmin, icrxmn, crxmax, 
     &              icrxmx)
               call min_max (cord(i,2), i, crymin, icrymn, crymax,  
     &              icrymx)
               call min_max (cord(i,3), i, crzmin, icrzmn, crzmax,  
     &              icrzmx)
            endif

            if(volume(i).ne.0.0) then
               call min_max (sx1(i), i, volmin, ivolmn, volmax, ivolmx)
               voltot = voltot + sx1(i)
            endif

         end do

c**** print out node numbers and other misc. info ****
         if (iout .ne. 0) write(iout, 6000)  nmfil(13)
         if (iptty .gt. 0)  write(iptty, 6000)  nmfil(13)
 6000    format(/, 1x, '**** analysis of input data on file ', a32, 
     .        ' ****', /)
         write(ischk, 6010)  neq
 6010    format(/, 1x, 'number of nodes in problem = ', i10, /)
      
         if (iccen .ne. 0)  write(ischk, 6020)
 6020    format(1x, 'tracer solution is enabled')
         
         if (istrs .ne. 0)  write(ischk, 6021)
 6021    format(1x, 'stress module is enabled')

         if (ico2 .ne. 0)  write(ischk, 6022)
 6022    format(1x, 'noncondesible gas module is enabled')

         if (idualp .ne. 0)  write(ischk, 6023)
 6023    format(1x, 'dual porosity solution is enabled')
         
         if (idpdp .ne. 0)  write(ischk, 6011)
 6011    format(1x, 'dual porosity/dual permeability solution is ',
     &        'enabled')

         if (ice .ne. 0)  write(ischk, 6024)
 6024    format(1x, 'ice module is enabled')

         if (iwelb .ne. 0)  write(ischk, 6025)
 6025    format(1x, 'wellbore module is enabled')

         if (idoff .le. 0)  then
            write(ischk, 6026)
 6026       format(1x, 'heat conduction only solution')
         else if (idof .eq. 2 .or. idof .eq. 4)  then 
            write(ischk, 6027)
 6027       format(1x, 'heat and mass transfer solution')
         end if

         if (intg .le. 0)  then
            write(ischk, 6028)
 6028       format(1x, 'lobatto (node point) quadrature used')
         else 
            write(ischk, 6029)
 6029       format(1x, 'gauss quadrature used')
         end if

         write(ischk, 6030)  1.0 - upwgt
 6030    format(1x, 'value of upwind weight = ', g12.5)

         if (iccen .ne. 0)  write(ischk, 6031)  1.0 - upwgta
 6031    format(1x, 'value of tracer upwind weight = ', g12.5)

         if (iporos .ne. 0)  write(ischk, 6032)
 6032    format(1x, 'porosity model is enabled')

         if (irdof .ne. 0)  write(ischk, 6033)
 6033    format(1x, 'reduced degree of freedom enabled')

         if (ivf .ne. 0)  write(ischk, 6034)
 6034    format(1x, 'finite volume calculations used')

c**** further analysis of storage ****

         if (istrs .eq. 0 .and. nq .gt. 1)  write(ischk, 6040)
 6040    format(1x, 'nq may be reduced to 1')

         if (ico2 .eq. 0 .and. n4 .gt. 1)  write(ischk, 6041)
 6041    format(1x, 'n4 may be reduced to 1')

         if (idualp .eq. 0 .and. n5 .gt. 1)  write(ischk, 6042)
 6042    format(1x, 'n5 may be reduced to 1')

         if (ice .eq. 0 .and. n6 .gt. 1)  write(ischk, 6043)
 6043    format(1x, 'n6 may be reduced to 1')

         if (iccen .eq. 0 .and. n7 .gt. 1)  write(ischk, 6044)
 6044    format(1x, 'n7 may be reduced to 1')

         if (iporos .eq. 0 .and. n8 .gt. 1)  write(ischk, 6046)
 6046    format(1x, 'n8 may be reduced to 1', /)

c**** check for possible error in dimension specification ****

         if (icnl .eq. 0)  then
            if (abs(zdiff) .lt. zero_t)  then
               write(ischk, 6050)
               if (iout .ne. 0) write(iout, 6050)
               if (iatty .gt. 0)  write(iatty, 6050)
 6050          format(1x, 'icnl (macro ctrl) check for 2-d problem')
            end if
         else
            if (abs(zdiff) .gt. zero_t .and. icnl .ne. 4)  then
               write(ischk, 6051)
               if (iout .ne. 0) write(iout, 6051)
               if (iatty .gt. 0)  write(iatty, 6051)
 6051          format(1x, 'icnl (macro ctrl) check for 3-d problem')
            end if
         end if

c**** zero out z coordinate in 2-d case ****

         if (icnl .eq. 0)  then
            dfctor = 1.0
         else
            dfctor = 0.0
         end if

C**** print reference values for air/water problem ****
        if(ico2.lt.0) then
 6053    format(1x, 'reference temperature, pressure for isothermal ',
     &        'EOS problems')      
         write(ischk, 6053)
c gaz 110819 changed crl(4,1) and crl(6,1) pref and tref         
         write(ischk, 6055) tref, pref
 6054    format(5x,"liquid : density, viscosity, compressibility ")
 6055    format(10x,3(1x, g15.7))
         write(ischk, 6054)
         write(ischk, 6055) crl(1,1), crl(2,1), crl(3,1) 
 6056    format(5x,"gas    : density, viscosity, compressibility ")
         write(ischk, 6056)
         den_vap = roc0*(273.0/(tref+273.0))*pref/0.101325
         write(ischk, 6055) den_vap, crl(5,1), den_vap/pref          
c**** print out max and min values ****
 6100    format(3x, g16.8, 3x, i8, 3(3x, g10.3))
        else
         write(ischk,*)
         write(ischk,*)'EOS data obtained from internal equations'
         write(ischk,*)
        endif

c**** pressure information ****
      
         if(iprsmx.le.neq) then
            icsub=0
         else if(iprsmx.ge.neq+1.and.iprsmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6101)
 6101    format(/, 1x, 'max total pressure, node number, ', 
     *     'x coordinate, y coordinate, z coordinate')
         write(ischk, 6100)  prsmax, iprsmx, cord(iprsmx-icsub, 1), 
     *        cord(iprsmx-icsub, 2), cord(iprsmx-icsub, 3)*dfctor

         if(iprsmn.le.neq) then
            icsub=0
         else if(iprsmn.ge.neq+1.and.iprsmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6103)
 6103    format(/, 1x, 'min total pressure, node number, ', 
     *        'x coordinate, y coordinate, z coordinate')
         write(ischk, 6100)  prsmin, iprsmn, cord(iprsmn-icsub, 1), 
     *        cord(iprsmn-icsub, 2), cord(iprsmn-icsub, 3)*dfctor

c**** temperature information ****

         if(itemmx.le.neq) then
            icsub=0
         else if(itemmx.ge.neq+1.and.itemmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6111)
 6111    format(/, 1x, 'max temperature, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  temmax, itemmx, cord(itemmx-icsub, 1), 
     *        cord(itemmx-icsub, 2), cord(itemmx-icsub, 3)*dfctor

         if(itemmn.le.neq) then
            icsub=0
         else if(itemmn.ge.neq+1.and.itemmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6113)
 6113    format(/, 1x, 'min temperature, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  temmin, itemmn, cord(itemmn-icsub, 1), 
     *        cord(itemmn-icsub, 2), cord(itemmn-icsub, 3)*dfctor

c**** saturation information ****

         if(isatmx.le.neq) then
            icsub=0
         else if(isatmx.ge.neq+1.and.isatmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6121)
 6121    format(/, 1x, 'max saturation, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  satmax, isatmx, cord(isatmx-icsub, 1), 
     *        cord(isatmx-icsub, 2), cord(isatmx-icsub, 3)*dfctor

         if(isatmn.le.neq) then
            icsub=0
         else if(isatmn.ge.neq+1.and.isatmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6123)
 6123    format(/, 1x, 'min saturation, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  satmin, isatmn, cord(isatmn-icsub, 1), 
     *        cord(isatmn-icsub, 2), cord(isatmn-icsub, 3)*dfctor

c**** air pressure information ****

         if (irdof .ne. 13 .or. ifree .ne. 0) then
            if (ico2 .ne. 0 ) then

               if(iprcmx.le.neq) then
                  icsub=0
               else if(iprcmx.ge.neq+1.and.iprcmx.le.neq+neq) then
                  icsub=neq
               else
                  icsub=neq+neq
               endif

               write(ischk, 6131)
 6131          format(/, 1x, 'max air pressure, node number, ', 
     *              'x coordinate, y coordinate, z coordinate')
               write(ischk, 6100) prcmax, iprcmx, cord(iprcmx-icsub,1), 
     *              cord(iprcmx-icsub, 2), cord(iprcmx-icsub, 3)*dfctor

               if(iprcmn.le.neq) then
                  icsub=0
               else if(iprcmn.ge.neq+1.and.iprcmn.le.neq+neq) then
                  icsub=neq
               else
                  icsub=neq+neq
               endif

               write(ischk, 6133)
 6133          format(/, 1x, 'min air pressure, node number, ', 
     *              'x coordinate, y coordinate, z coordinate')
               write(ischk, 6100) prcmin, iprcmn, cord(iprcmn-icsub,1), 
     *              cord(iprcmn-icsub, 2), cord(iprcmn-icsub, 3)*dfctor
            
            end if
         end if

c**** porosity information ****

         if(ipormx.le.neq) then
            icsub=0
         else if(ipormx.ge.neq+1.and.ipormx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6141)
 6141    format(/, 1x, 'max porosity, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  pormax, ipormx, cord(ipormx-icsub, 1), 
     *        cord(ipormx-icsub, 2), cord(ipormx-icsub, 3)*dfctor

         if(ipormn.le.neq) then
            icsub=0
         else if(ipormn.ge.neq+1.and.ipormn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6143)
 6143    format(/, 1x, 'min porosity, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  pormin, ipormn, cord(ipormn-icsub, 1), 
     *        cord(ipormn-icsub, 2), cord(ipormn-icsub, 3)*dfctor

c**** permeability information ****

         if(ipermx.le.neq) then
            icsub=0
         else if(ipermx.ge.neq+1.and.ipermx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6151)
 6151    format(/, 1x, 'max permeability, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  permmx, ipermx, cord(ipermx-icsub, 1), 
     *        cord(ipermx-icsub, 2), cord(ipermx-icsub, 3)*dfctor

         if(ipermn.le.neq) then
            icsub=0
         else if(ipermn.ge.neq+1.and.ipermn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6153)
 6153    format(/, 1x, 'min permeability, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  permmn, ipermn, cord(ipermn-icsub, 1), 
     *        cord(ipermn-icsub, 2), cord(ipermn-icsub, 3)*dfctor

c**** rock density information ****

         if(idenmx.le.neq) then
            icsub=0
         else if(idenmx.ge.neq+1.and.idenmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6161)
 6161    format(/, 1x, 'max rock density, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  denmmx, idenmx, cord(idenmx-icsub, 1), 
     *        cord(idenmx-icsub, 2), cord(idenmx-icsub, 3)*dfctor

         if(idenmn.le.neq) then
            icsub=0
         else if(idenmn.ge.neq+1.and.idenmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6163)
 6163    format(/, 1x, 'min rock density, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  denmmn, idenmn, cord(idenmn-icsub, 1), 
     *        cord(idenmn-icsub, 2), cord(idenmn-icsub, 3)*dfctor

c**** heat capacity information ****

         if(icprmx.le.neq) then
            icsub=0
         else if(icprmx.ge.neq+1.and.icprmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6171)
 6171    format(/, 1x, 'max heat capacity, node number, ', 
     *        'x coordinate, y coordinate, z coordinate')
         write(ischk, 6100)  cprmmx, icprmx, cord(icprmx-icsub, 1), 
     *        cord(icprmx-icsub, 2), cord(icprmx-icsub, 3)*dfctor

         if(icprmn.le.neq) then
            icsub=0
         else if(icprmn.ge.neq+1.and.icprmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6173)
 6173    format(/, 1x, 'min heat capacity, node number, ', 
     *        'x coordinate, y coordinate, z coordinate')
         write(ischk, 6100)  cprmmn, icprmn, cord(icprmn-icsub, 1), 
     *        cord(icprmn-icsub, 2), cord(icprmn-icsub, 3)*dfctor
    
c**** t conductivity information ****
     
         if(ithcmx.le.neq) then
            icsub=0
         else if(ithcmx.ge.neq+1.and.ithcmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6181)
 6181    format(/, 1x, 'max t conductivity, node number, ', 
     *        'x coordinate, y coordinate, z coordinate')
         write(ischk, 6100)  thcmmx, ithcmx, cord(ithcmx-icsub, 1), 
     *        cord(ithcmx-icsub, 2), cord(ithcmx-icsub, 3)*dfctor
      
         if(ithcmn.le.neq) then
            icsub=0
         else if(ithcmn.ge.neq+1.and.ithcmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6183)
 6183    format(/, 1x, 'min t conductivity, node number, ', 
     *        'x coordinate, y coordinate, z coordinate')
         write(ischk, 6100)  thcmmn, ithcmn, cord(ithcmn-icsub, 1), 
     *        cord(ithcmn-icsub, 2), cord(ithcmn-icsub, 3)*dfctor

c**** x coordinate information ****

         write(ischk, 6191)
 6191    format(/, 1x, 'max x coordinate, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  crxmax, icrxmx, cord(icrxmx, 1), 
     *        cord(icrxmx, 2), cord(icrxmx, 3)*dfctor

         write(ischk, 6193)
 6193    format(/, 1x, 'min x coordinate, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  crxmin, icrxmn, cord(icrxmn, 1), 
     *        cord(icrxmn, 2), cord(icrxmn, 3)*dfctor

c**** y coordinate information ****

         write(ischk, 6201)
 6201    format(/, 1x, 'max y coordinate, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  crymax, icrymx, cord(icrymx, 1), 
     *        cord(icrymx, 2), cord(icrymx, 3)*dfctor

         write(ischk, 6203)
 6203    format(/, 1x, 'min y coordinate, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  crymin, icrymn, cord(icrymn, 1), 
     *        cord(icrymn, 2), cord(icrymn, 3)*dfctor

c**** z coordinate information ****

         if (icnl .eq. 0) then
         
            write(ischk, 6211)
 6211       format(/, 1x, 'max z coordinate, node number, ', 
     *           'x coordinate, y coordinate, z coordinate')
            write(ischk, 6100)  crzmax, icrzmx, cord(icrzmx, 1), 
     *           cord(icrzmx, 2), cord(icrzmx, 3)*dfctor

            write(ischk, 6213)
 6213       format(/, 1x, 'min z coordinate, node number, ', 
     *           'x coordinate, y coordinate, z coordinate')
            write(ischk, 6100)  crzmin, icrzmn, cord(icrzmn, 1), 
     *           cord(icrzmn, 2), cord(icrzmn, 3)*dfctor

         end if

c**** cell volume information ****

         if(ivolmx.le.neq) then
            icsub=0
         else if(ivolmx.ge.neq+1.and.ivolmx.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6221)
 6221    format(/, 1x, 'max cell volume, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  volmax, ivolmx, cord(ivolmx-icsub, 1), 
     *        cord(ivolmx-icsub, 2), cord(ivolmx-icsub, 3)*dfctor

         if(ivolmn.le.neq) then
            icsub=0
         else if(ivolmn.ge.neq+1.and.ivolmn.le.neq+neq) then
            icsub=neq
         else
            icsub=neq+neq
         endif

         write(ischk, 6223)
 6223    format(/, 1x, 'min cell volume, node number, x coordinate, ', 
     *        'y coordinate, z coordinate')
         write(ischk, 6100)  volmin, ivolmn, cord(ivolmn-icsub, 1), 
     *        cord(ivolmn-icsub, 2), cord(ivolmn-icsub, 3)*dfctor

c**** print initial mass, steam mass, energy ****

         write(ischk, 6231)  am0, astmo, ame, voltot
 6231    format(/, 1x, 'initial mass = ', g10.3, 
     *        /, 1x, '       steam = ', g10.3, 
     *        /, 1x, '      energy = ', g10.3, 
     *        /, 1x, 'total volume = ', g10.3)
      end if
c
c if finite element coefficients are generated check them
c
      if(lda.le.0) then
         if(irun.eq.1.and.ianpe.eq.0) then
            call check_sx(neq,ianpe,nr,nelm,istrw,nelmdg,ischk,iptty,
     &           iout,ierr,min_nr,isox,isoy,isoz,intg,sx1,sx)
            if (ischk .ne. 0) then
               write(ischk,9000) min_nr,nr
               write(ischk,*) ' '
            end if
            if (iptty .ne. 0) then
               write(iptty,9000) min_nr,nr
               write(iptty,*) ' '
            end if
 9000       format(/,'storage for fe coefficients ',i10,
     &           ' allocated ',i10)
         else
            if (iptty .ne. 0) then
               write(iptty,*)'anisotrpic case: coefficients not checked'
               write(iptty,*) ' '
            endif
         endif
c     nbytes = lenreal * min_nr * 6
c     call reallocf(sx , nbytes)
c     
c set minimun volume for zero volume nodes
c this can happen for multiply defined nodes
c set the volume(i) array to non zero to distinguish it
c from reservoir nodes
c
         do i=1,n
c gaz 021919 added hi perm for flow pass through             
            if(sx1(i).le.0.0) then
               sx1(i)=vol_fac
               pnx(i) = perm_fac*pnx(i)
               pny(i) = perm_fac*pny(i)
               pnz(i) = perm_fac*pnz(i)
c gaz 022319 correct initial mass and energy
               ame = ame - deneh(i)*volume(i)
               am0 = am0 - denh(i)*volume(i)
               astmo = astmo - dstm (i)
               volume(i)=sx1(i)
            endif
         enddo
      else
c check  read-in connectivity 
         if(irun.eq.1.and.ianpe.eq.0) then
            call check_sx(neq,ianpe,nr,nelm,istrw,nelmdg,ischk,iptty,
     &           iout,ierr,min_nr,isox,isoy,isoz,intg,sx1,sx)
            if (ischk .ne. 0) then
               write(ischk,9000) min_nr,nr
               write(ischk,*) ' '
            end if
            if (iptty .ne. 0) then
               write(iptty,9000) min_nr,nr
               write(iptty,*) ' '
            end if
c 9000       format(/,'storage for fe coefficients ',i10,
c     &           ' allocated ',i10)
         else
            if (iptty .ne. 0) then
               write(iptty,*)'anisotrpic case: coefficients not checked'
               write(iptty,*) ' '
            endif
         endif 
c checking for negative volumes 
         write(ierr,9002)
 9002 format('ckecking for (and zeroing) negative volumes',/,
     &   t13,'node',t26,'zone',t37,'volume',t55,'X',t70,'Y',t85,'Z')
         j = 0
         do i=1,n
            if(sx1(i).le.0.0) then  
             j = j+1
             write(ierr,9001) i,izonef(i),sx1(i),cord(i,1),cord(i,2),
     &         cord(i,3)
             if(ivol_old.gt.0) then             
               sx1(i) = 0.
               volume(i) = 0.
               pnx(i) = 1.e-20
               pny(i) = 1.e-20
               pnz(i) = 1.e-20
c gaz 021919 added hi perm for flow pass through             
            else
               sx1(i)= vol_fac
               pnx(i) = perm_fac*pnx(i)
               pny(i) = perm_fac*pny(i)
               pnz(i) = perm_fac*pnz(i)
c gaz 022319 correct initial mass and energy
               ame = ame - deneh(i)*volume(i)
               am0 = am0 - denh(i)*volume(i)
               astmo = astmo - dstm (i)
               volume(i)=sx1(i)
            endif 
            endif
         enddo
         write(ierr,*)
         write(ierr,*) 'number of neg vol nodes = ',j
9001     format('neg vol',1x,i10,1x,i7,1x,1p,g14.3,1x,g14.3,
     &           1x,g14.3,1x,g14.3)
      end if

c     With the new particle tracking option (mptr),
c     we need to keep this open and let FORTRAN close
c     at the end of the run

      if(.not.ptrak) close  (ischk, status = 'keep')

      end

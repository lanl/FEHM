      subroutine input(cnum, simnum)
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
CD1 Control reading of input data file.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 22-DEC-93    Z. Dash        22      Add prolog/major cleanup.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/input.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.21   Wed Jun 12 15:30:20 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.20   Mon Jun 10 11:55:40 1996   hend
CD2 Moved input_msg to character*80
CD2 
CD2    Rev 1.19   Mon Jun 03 11:18:06 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.18   Fri May 31 10:39:26 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.17   Tue May 28 11:18:38 1996   gaz
CD2 diabled BOUS macro if icons=0,enabled otherwise
CD2
CD2    Rev 1.16   Wed Feb 07 10:18:32 1996   gaz
CD2 added macro bous for constant transmissibility fluid parameters
CD2 
CD2    Rev 1.15   Tue Jan 30 13:19:22 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.14   Wed Jan 17 14:13:14 1996   hend
CD2 Added use of parser for input lines to function on sgi
CD2 
CD2    Rev 1.13   12/13/95 10:29:04   robinson
CD2 Incorporated new zeolite hydration module
CD2 
CD2    Rev 1.12   12/13/95 08:43:22   gaz
CD2 deleted macro cap
CD2 
CD2    Rev 1.11   08/02/95 17:04:14   gaz
CD2 new call to mdnodes
CD2 
CD2    Rev 1.10   04/27/95 18:26:12   llt
CD2 added macro itup
CD2 
CD2    Rev 1.9   04/27/95 15:15:08   gaz
CD2 added -tmch option to macro ctrl 
CD2 
CD2    Rev 1.8   03/24/95 09:44:04   llt
CD2 added macro iupk
CD2 
CD2    Rev 1.7   03/10/95 10:57:30   llt
CD2 defined parser arguments in pseudocode
CD2 
CD2    Rev 1.6   02/21/95 16:41:26   llt
CD2 psuedcode for avs and corrected branch (changes by zvd)
CD2 
CD2    Rev 1.5   01/28/95 14:03:52   llt
CD2 modified for new particle tracking module
CD2 
CD2    Rev 1.4   01/03/95 13:24:38   llt
CD2 Added parser, so ibm could read characters in macro definitions.
CD2 
CD2    Rev 1.3   06/20/94 11:12:32   zvd
CD2 Added ierr unit number for error output.
CD2 
CD2    Rev 1.2   04/06/94 11:14:20   tam
CD2 for avsio: added line to skip comment lines starting with '#'
CD2 added read_avs_io from input or avs.in files
CD2 added error check for multiply defined formats
CD2 
CD2    Rev 1.1   03/18/94 16:03:46   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:25:16   pvcs
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
CD3   cnum            INT      I    Number of times zone has been called
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   inpt                     I    Main input data file
CD3   iout                     O    General output file
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
CD4   contim          REAL*8   faar   Interval (days) for contour plot output
CD4   delat           REAL*8   faar   Given tolerance for gas pressure
CD4   delpt           REAL*8   faar   Given tolerance for pressure
CD4   delst           REAL*8   faar   Given tolerance for saturation
CD4   deltt           REAL*8   faar   Given tolerance for temperature
CD4   depcng          REAL*8   faar   Depth at which temperature gradient
CD4                                     changes
CD4   emiss           REAL*8   faar   ?
CD4   g1              REAL*8   faar   Iteration accuracy control parameter
CD4   g2              REAL*8   faar   Iteration accuracy control parameter
CD4   g3              REAL*8   faar   Iteration accuracy control parameter
CD4   grad2           REAL*8   faar   Parameter in description of temperature
CD4                                     gradient
CD4   gradnt          REAL*8   faar   Parameter used in description of
CD4                                     temperature gradient
CD4   iadif           INT      faai   ?
CD4   iad_up          INT      faai   Upwind parameter
CD4   iback           INT      david1 LU factorization save parameter
CD4   icapp           INT      faai   Indicates capillary pressure model
CD4   ico2            INT      faai   Indicates if noncondensible gas solution
CD4                                     is enabled
CD4   icoupl          INT      david1 Number of SOR iterations
CD4   inpt            INT      faai   Unit number for input file
CD4   intg            INT      faai   Indicates integration type used
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   irdof           INT      david1 Reduced degree of freedom model used
CD4   islord          INT      david1 Parameter used in the reduced degree of
CD4                                     freedom model
CD4   ivfcal          INT      faai   Indicates if vfcal subroutine will be
CD4                                     called
CD4   ncntr           INT      faai   Contour plot interval
CD4   nrxns           INT      readrxn   Number of reactions
CD4   overf           REAL*8   faar   Over relaxation factor for sor equations
CD4   pein            REAL*8   faar   Initial pressure of problem
CD4   quad            REAL*8   faar   Parameter used in temperature gradient
CD4   rnmax           REAL*8   faar   Maximum run time allowed
CD4   rxn_flag        INT      henry  Flag denoting whether to call chemical
CD4                                     reaction routine
CD4   rxn_interval    INT      setup  Parameter for determining when to
CD4                                     perform full iteration of concentrations
CD4   sssol           CHAR     faac   Indicates if initial steady state solution
CD4                                     is needed
CD4   tin             REAL*8   faar   Parameter used in temperature gradient
CD4   tin0            REAL*8   faar   Initial problem temperature
CD4   tin1            REAL*8   faar   Parameter used in temperature gradient
CD4   tmch            REAL*8   faar   Machine tolerance
CD4   tort            REAL*8   faar   ?
CD4   wdd             CHAR     faac   Character input string
CD4   wdd1            CHAR     faac   Alternate character input string
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   null1            LOGICAL  Check for null lines or 0's in lines
CD4   zone                     Set zone information
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
CD5   cmsg            CHAR     Parsed string arguments of type character
CD5   idof            INT      Parameter indicating type of solution required
CD5   iieosd          INT      Flag denoting if simple thermodynamics are to
CD5                              be used
CD5   imsg            INT      Parsed string arguments of type integer
CD5   input_msg       CHAR     Input string to be parsed
CD5   kk              INT      Argument to use when calling user subroutine
CD5   macro           CHAR     Control statement designator
CD5   msg             INT      Type of agrument in parsed string
CD5   nwds            INT      Number of words in parsed string
CD5   xmsg            REAL*8   Parsed string arguments of type real
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
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
CPS BEGIN input
CPS 
CPS   reread input title
CPS   
CPS   REPEAT
CPS   
CPS     read input line
CPS     IF line is a comment
CPS        loop and read next input line
CPS     END IF
CPS     read input macro from input line just read
CPS
CPS     write macro identifier and input file unit to output file and tty 
CPS      if enabled
CPS     
CPS     EXIT IF macro read is stop
CPS       
CPS     ELSE IF macro read is adif
CPS        set iadif to 1 and read adif data
CPS     ELSE IF macro read is airwater
CPS        set ico2 to -1 and call airctr to read isothermal air-water 
CPS         transport input data
CPS     ELSE IF macro read is alti
CPS        [skip, alternate geometric data is read by incoord]
CPS     ELSE IF macro read is cap
CPS        set icapp to 1 and call cappr to read in capillary pressure 
CPS         information
CPS     ELSE IF macro read is coor
CPS        REPEAT
CPS        [skip, geometric data is read by incoord]
CPS        UNTIL end of coor data is found
CPS     ELSE IF macro read is cond
CPS        call incond to read thermal conductivity data
CPS     ELSE IF macro read is cont
CPS        read input line
CPS        IF contour file format is specified
CPS           reread line
CPS           call parse_string to get format type, time step interval, 
CPS            and time interval for contour plots
CPS           IF contour file format is avs
CPS              call read_avs_io to read and set avs io flags
CPS           END IF
CPS        ELSE
CPS           reread line to get time step interval, and time interval 
CPS            for contour plots
CPS        END IF
CPS        write contour format identifier to output file
CPS        IF time step interval or time interval are too small
CPS           set them to a large value
CPS        END IF
CPS     ELSE IF macro read is ctrl
CPS        call inctrl to read control data
CPS     ELSE IF macro read is dof
CPS        [doesn't do anything yet]
CPS     ELSE IF macro read is dpdp
CPS        set idpdp to 1 and call dpdp to read dual porosity/dual 
CPS         permeability data
CPS     ELSE IF macro read is dual
CPS        set idualp to 1 and call dual to read dual porosity
CPS     ELSE IF macro read is elem
CPS        REPEAT
CPS        [skip, geometric data is read by incoord]
CPS        UNTIL end of elem data is found
CPS     ELSE IF macro read is eos
CPS        set iieosd to 0 and call sther to read equation of state data
CPS     ELSE IF macro read is finv
CPS        set flags for finite volume calculations
CPS     ELSE IF macro read is flow
CPS        call inflow to read flow data
CPS     ELSE IF macro read is flo2
CPS        call inflo2 to read flow data for 3D planes
CPS     ELSE IF macro read is flo3
CPS        call inflo3 to read flow data for 3D planes
CPS     ELSE IF macro read is flx3
CPS        set nflx to 1 and call flxo to calculate intermode fluxes
CPS     ELSE IF macro read is hflx
CPS        call inhflx to read heat flux data
CPS     ELSE IF macro read is ice
CPS        set ice to 1 and call sice to read ice data
CPS     ELSE IF macro read is impf
CPS        set impf to 1 and read fractional variable change data
CPS     ELSE IF macro read is init
CPS        read initial value data for pressure and temperature
CPS     ELSE IF macro read is iter
CPS        read iteration parameters and set maximum run time
CPS     ELSE IF macro read is itup
CPS        read parameter for iad_up
CPS     ELSE IF macro read is iupk
CPS        set iupk to 1
CPS     ELSE IF macro read is ivfc
CPS        set ivfcal to 1 to enable vfcal routine
CPS     ELSE IF macro read is mdno
CPS        call mdnodes to read additional node connections
CPS     ELSE IF macro read is ngas or co2i
CPS        set ico2 to 3 and call co2ctr to read noncondensible gas data
CPS     ELSE IF macro read is node or nod2 or nod3
CPS        call innode to read nodal output values
CPS     ELSE IF macro read is num
CPS        [doesn't do anything yet]
CPS     ELSE IF macro read is perm
CPS        call inperm to read permeability data
CPS     ELSE IF macro read is ppor
CPS        set iporos to 1 and call porosi to read pressure and 
CPS         temperature dependent porosity and permeability data 
CPS     ELSE IF macro read is pres
CPS        call inpress to read pressure and temperature or saturation 
CPS        data
CPS     ELSE IF macro read is ptrk
CPS        call inptrk to read particle tracking input data
CPS     ELSE IF macro read is renm
CPS        call renum to read new node numbers
CPS     ELSE IF macro read is rflx
CPS        read and set emiss
CPS        call inhflx to read radiation flux data
CPS     ELSE IF macro read is rlp
CPS        call rlperm to read in relative permeability data
CPS     ELSE IF macro read is rock
CPS        call inrock to read rock propety data
CPS     ELSE IF macro read is rxn
CPS        call read_rxn to read reaction rate data
CPS        IF reactions are specified
CPS           set rxn_flag to 1
CPS        ELSE
CPS           set rxn_flag and rxn_interval to 0
CPS        ENDIF
CPS     ELSE IF macro read is sol
CPS        read solution and integration type
CPS        IF element integration type is <= 0
CPS           use Lobatto quadrature
CPS        ELSE IF 
CPS           use Gauss quadrature
CPS        END IF
CPS     ELSE IF macro read is solv
CPS        [doesn't do anything yet]
CPS     ELSE IF macro read is stea
CPS        set steady state solution to yes
CPS     ELSE IF macro read is strs
CPS        set istrs to 1 and call stress [does nothing in this version]
CPS     ELSE IF macro read is text
CPS        REPEAT
CPS          read input line
CPS          IF line is not null
CPS             write text to output file
CPS          END IF
CPS        UNTIL end of text (null line) is found
CPS        write macro identifier and input file unit to output file 
CPS         and tty if enabled
CPS     ELSE IF macro read is time
CPS        call intime to read time step control data
CPS     ELSE IF macro read is thic
CPS        call thickness to compute thickness information
CPS     ELSE IF macro read is trac
CPS        set iccen to 1 and call concen to read tracer concentration 
CPS        data
CPS     ELSE IF macro read is user
CPS        read user subroutine calling parameter
CPS        call user subroutine with calling parameter
CPS     ELSE IF macro read is vapl
CPS        set ivapl to 1
CPS     ELSE IF macro read is vcon
CPS        set ivcond to 1 and call vcon to read in nonlinear thermal 
CPS         conductivity information
CPS     ELSE IF macro read is dvel
CPS        set nflx to -1 and call flxo to calculate intermode 
CPS        velocities  
CPS     ELSE IF macro read is wlbr
CPS        set iwelb to 1 and call welbor to read wellbore style input
CPS        exit input
CPS     ELSE IF macro read is zone
CPS        call zone to read zone data
CPS     ELSE 
CPS        write input error message and terminate program 
CPS     END IF
CPS
CPS   UNTIL stop is found
CPS   
CPS   set total number of nodes
CPS   
CPS   IF dual and dpdp solutions are both enabled
CPS      disable dual solution and write message to data check file
CPS   END IF
CPS   
CPS   call dual, dpdp, and flxo to complete initialization of problem 
CPS   data
CPS
CPS END input
CPS
C***********************************************************************

      use comai
      use combi
      use comchem
      use comci
      use comco2
      use comdi
      use comdti
      use comevap, only : evaporation_flag
      use compart
      use comriv
      use comrxni
      use comsi, only : cnum_stress
      use comsplitts
      use comsptr
      use comwt
      use davidi
      implicit none

      integer i, izone, inode
      integer idum, j, jj
      real tmpli
      logical null1
      character*80 input_msg, dummy_line
      character*4 macro, macro1, chard
      integer cnum,iieosd,inptorig,kk,msg(4),nwds,imsg(4)
      real*8 xmsg(4), dummyreal, simnum
      character*32 cmsg(4)

      sssol = 'no  '
      altc = 'fehm'
      inptorig = inpt

      imsg = 0
      xmsg = 0.
      cmsg = ''

      read (inptorig, '(a80)') wdd

  100 continue

      read (inptorig, '(a80)') wdd1
      if (wdd1(1:1) .eq. '#') go to 100
      read (wdd1, '(a4)') macro
      if (iout .ne. 0) write(iout, 6000) macro, inptorig
      if (iptty .gt. 0) write(iptty, 6000) macro, inptorig
 6000 format(1x, '**** input title : ', a4, ' **** inpt = ', 
     .     i3, ' ****')

      if (macro .eq. 'stop') then
! zvd 06/02/03 If we have reached the end of the input, 
! set tracer porosity if necessary  
         if ((iccen .ne. 0) .or. sptrak) then
            if (tpor_flag) then
               do i = 1, n0
                  if (ps_trac(i) .eq. 10.) then
                     ps_trac(i) = ps(i)
                  end if
               end do
            else
               ps_trac(1:n0) = ps(1:n0)
            endif
         end if
c**** exit loop if stop is read ****
         go to 200
      else
         call start_macro(inptorig, inpt, macro)
      end if

      if (macro .eq. 'adif') then
c**** air-water vapor diffusion ****
        iadif = 1
        read(inpt,*) tort

      else if (macro(1:3) .eq. 'air') then
c**** isothermal air - water transport ****
         read(inpt,*) ico2
         ico2 = -2
         call airctr (0, 0)

      else if (macro(1:4) .eq. 'para') then
c**** parallel FEHM implimentation (isothermal only)***
         ipara = 1
         call paractr (0)

      else if (macro(1:4) .eq. 'wgtu') then
c**** areas, weights(user-defined) for boundary conditions ***
         call area_vol_weightctr (0)

      else if (macro(1:4) .eq. 'wtsi') then
c**** free surface calculation with head(saturated only)***
c**** water table simple
         call wtsictr (0)

      else if (macro(1:4) .eq. 'grad') then
c**** input af inital and boundary conditions with gradients ****
         igrad = 1
         call gradctr (0)

      else if (macro(1:4) .eq. 'conv') then
c**** input with mixed physics ****
         backspace inpt
         read(inpt,'(a80)') input_msg 
         read(input_msg,*,end= 899) macro, headconv
         go to 900 
 899     headconv = 0.0
 900     continue
         if(iconv.eq.0) iconv=1
         call convctr (0)

      else if (macro .eq. 'alti') then
c**** alternate element input, read in incoord subroutine ****

      else if (macro .eq. 'head') then
c**** input in terms of head, not pressures   ****
c bous macro enabled
c must use isothermal physics
c
         irdof = 13
         backspace inpt
         read(inpt,'(a80)') input_msg 
         read(input_msg,*,end= 999) macro, head0
         go to 1000
 999     head0 = 0.0
 1000    continue
         ihead=1
c         call headctr(0,0,0.0,0.0)

      else if (macro .eq. 'chea') then
c
c**** output in terms of head, not pressures   ****
c
         backspace inpt
         read(inpt,'(a80)') input_msg 
         read(input_msg,*,end= 995) macro, head0, temp0, pres0, sat_ich,
     &	    head_id
         go to 1005
 995     head0 = 0.0
         temp0 = 20.0
         pres0 = 0.1   
         sat_ich = 0.0
         head_id= 0.0
 1005    continue
         call water_density(temp0, pres0,rol0)
         ichead=1
         if(.not.allocated(head)) allocate(head(n0))

      else if (macro .eq. 'bous') then
c constant densities etc for flow terms
         read(inpt,*) icons
         if(icons.eq.0) then
            icons=10000
         else
            icons=1
         endif

      else if (macro .eq. 'nobr') then
c don't break connection between nodes with boundary conditions
         inobr = 1

      else if (macro .eq. 'cden') then
         if(nspeci.gt.0) then
            read(inpt,*) ispcden
            read(inpt,*) factcden
            if(ispcden.le.nspeci) then
               cden = .true.
            else
               if (iout .ne. 0) then
                  write(iout,*)'ispcden > nspeci, cden'
                  write(iout,*)'macro ignored'
               end if
               if(iptty.gt.0) then
                  write(iout,*)'ispcden > nspeci, cden'
                  write(iout,*)'macro ignored'
               end if
            end if
         else
            if (iout .ne. 0) then
               write(iout,*)'No solute transport, cden'
               write(iout,*)'macro ignored'
            end if
            if(iptty.gt.0) then
               write(iptty,*)'No solute transport, cden'
               write(iptty,*)'macro ignored'
            end if
         end if

      else if (macro .eq. 'cgdp') then
c**** rate-limited gdpm node identifcation (cgdpm)****
         igdpm_rate = 1
         call gdpm_corr(0)

      else if (macro .eq. 'coor') then
c**** node coordinate data, read in incoord subroutine ****
 110     continue
         read (inpt, '(a80)') wdd1
         if (null1(wdd1)) go to 115 
         go  to  110
 115     continue

      else if (macro .eq. 'cond') then
c**** thermal conductivity data ****
         call incond

      else if (macro .eq. 'conn') then
c**** thermal conductivity data ****
         connect_out = .true.

      else if (macro .eq. 'cont') then
         icont = 1
c**** contour plot information ****
         read(inpt, '(a80)') wdd1
         read(wdd1, '(a4)') altc
         backspace inpt
         if(altc .eq. 'ptrn' .or. altc .eq. 'ment' .or. altc .eq. 'fehm'
     &        .or. altc .eq. 'free' .or. altc(1:3) .eq. 'avs' .or. 
     &        altc(1:3) .eq. 'sur' .or. altc(1:3) .eq. 'tec') then
            read(inpt, '(a)') input_msg
            call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
            altc = cmsg(1)
            ncntr = imsg(2)
            if (msg(3).eq.1) then
               contim=imsg(3)
            else
               contim = xmsg(3)
            endif
            if(cmsg(4)(1:4).eq.'time') then
               contim = -abs(contim)
            endif
888         continue
            if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'sur'
     &           .or. altc(1:3) .eq. 'tec') then
               call read_avs_io(inpt)
            end if
         else
            read(inpt,*) ncntr, contim
            altc = 'fehm'
         endif
         if (iout .ne. 0) write(iout, 6010) altc  
 6010    format(1x, '**** contour format defined  : ', a4, ' ****')
         if (ncntr .le. 0) ncntr = 10000000
         if (abs(contim) .le. zero_t) then
            contim_rip = -abs(contim)
            contim = 1.0d+30
         else
            contim_rip = 1.0d+30
         end if

      else if (macro .eq. 'ctrl') then
c**** control variables ****
         call inctrl

      else if (macro .eq. 'dof ') then
c**** degrees of freedom (doesn't do anything yet) ****

      else if (macro .eq. 'dpdp') then
c**** dual porosity/dual permeability ****
         idpdp = 1
         call dpdp (0)

      else if (macro .eq. 'dual') then
c**** dual porosity ****
         idualp =  1
         call dual (0)

      else if (macro .eq. 'elem') then
c**** element node data, read in incoord subroutine  ****
 120     continue
         read (inpt, '(a80)') wdd1
         if (null1(wdd1)) go to 125 
         go  to  120
 125     continue

      else if (macro .eq. 'eos ') then
c**** forms simple water eos and/or changes the eos set number ****
         iieosd =  0
         call sther (iieosd)

      else if (macro .eq. 'evap') then
c**** set flag and input name of file with evaporation nodes ****
         evaporation_flag = .true.
         call evaporation(0)

      else if (macro .eq. 'fdm ') then
c finite difference input 
         jj = 0
 220     continue
         read (inpt, '(a80)') wdd1
         do i = 1,80
          if(wdd1(i:i).ne.' ') go to 220
         enddo
         jj = jj +1
         if(jj.le.2) go to 220

      else if (macro .eq. 'finv') then
c ****  finite volume calculations
c zvd 05-Jan-2007 -- this is now the default, do nothing but
c issue warning if FDM calculations are invoked
         if(ivf.eq.-1) then
            if (iout .ne. 0) write (iout,*) 
     &           '>>> Warning: finv not set for FDM calculations <<<' 
            if (iptty.ne.0) write (iptty,*)
     &           '>>> Warning: finv not set for FDM calculations <<<' 
         end if

      else if (macro .eq. 'nfin') then
c zvd 05-Jan-2007
c****  finite element calculations
! Default is finite volume calculations, so now I am turning them off
         ivf = 0
         mlz = 0

      else if (macro .eq. 'flow') then
c**** read in sources and sinks - original node ****
         call inflow

      else if (macro .eq. 'flwt') then
         call inflow_wt

      else if (macro .eq. 'flo2') then
c**** read in sources and sinks - planes in 3D ****
         call inflo2

      else if (macro .eq. 'flo3') then
c**** read in sources and sinks - seepage faces and drainage ****
         call inflo3

      else if (macro .eq. 'floa') then
c**** read in sources and sinks - seepage faces and drainage ****
         call infloa

      else if (macro .eq. 'flgh') then
c**** generalized head ****
         call inflogh(0)

      else if (macro .eq. 'ftsc') then
c**** flux correction for saturations over 1
         iflux_ts=1
         call cascade_sat(-1)

      else if (macro .eq. 'exrl') then
c**** explicit evaluation of relative permeability ****

         read(inpt,*) iexrlp

      else if (macro .eq. 'svar') then
c**** switch varables, limited now to pressure-enthalpy ****

         ivar=1
         call varctr(0)

      else if (macro .eq. 'boun') then
c**** read in sources and sinks - time dependent mode **** 
         iboun=1
         backspace (inpt)
         read(inpt,'(a80)') input_msg
         call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
         if (nwds .ge. 2) then
            if (msg(2) .eq. 3) then
               if (cmsg(2) .eq. 'nout' .or. cmsg(2) .eq. 'NOUT')
     &              boun_out = .false.
            end if
         end if
         if (mmodel .eq. 0) call flow_boundary_conditions(0)
         call flow_boundary_conditions(1)

      else if (macro .eq. 'flxn') then
c**** output source/sink node fluxes ****
         iflxn = 1

      else if (macro .eq. 'flxo') then
c**** calculate intermode fluxes ****
         ivelo = 1
         call flxo(0)

      else if (macro .eq. 'cflx') then
c**** calculate concentration flux ****
         backspace (inpt)
         read(inpt,'(a80)') input_msg
         call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
         if (nwds .gt. 1) then
            cflx_var = .false.
            do i = 2, nwds
               if (msg(i) .ne. 1) exit
               cflx_var(imsg(i)) = .true.
            end do
         else
            cflx_var = .true.
         end if 
         read(inpt,*) cflxz
         if(.not.allocated(icflxz)) allocate(icflxz(max(1,cflxz)))
         if(.not.allocated(izoncflxz)) allocate(izoncflxz(n0))   

         read(inpt,*)(icflxz(i),i=1,cflxz)

c     Loop over each zone for determining izoncflxz array

         izoncflxz = 0
         do izone = 1, cflxz
            do inode = 1, n0
               if(izonef(inode).eq.icflxz(izone)) then
                  izoncflxz(inode) = izone
               end if
            end do
         end do

      else if (macro .eq. 'flxz') then
c**** calculate intermode fluxes ****
         read(inpt,*) nflxz
         if(.not.allocated(iflxz)) allocate(iflxz(max(1,nflxz)))
         if(.not.allocated(izoneflxz)) allocate(izoneflxz(n0))   
         if(.not.allocated(totboun)) allocate(totboun(max(1,nflxz)))
         totboun = 0.0

         read(inpt,*)(iflxz(i),i=1,nflxz)

c     Loop over each zone for determining izoneflxz array

         izoneflxz = 0
         do izone = 1, nflxz
            do inode = 1, n0
               if(izonef(inode).eq.iflxz(izone)) then
                  izoneflxz(inode) = izone
               end if
            end do
         end do


c**** read in Generalized dual porosity information
      else  if (macro .eq. 'gdpm') then

         call ingdpm

      else  if (macro .eq. 'hflx') then
c**** read in heat source term(w) ****
         call inhflx(macro)

      else if (macro .eq. 'hist') then
c**** read in history output selections ****
         call inhist

      else  if (macro .eq. 'hyco') then
c**** read in hydraulic conductivity**
         if (idoff .eq. -1) then
            if (iout .ne. 0) write (iout, *) 
     &           'Hydraulic conductivity not read: ',
     .           'heat conduction only problem'
 123        read (inpt, '(a80)') wdd1
            if(null1(wdd1)) goto 321
            goto 123
         else
            call inhyco
         end if
 321     continue

      else if (macro .eq. 'ice ' .or. macro .eq. 'meth') then
c**** read in ice or hydrate information ****
         ice = 1
c also set to isothermal air-water to use some code structure
         ico2 = -2
c gaz 10-15-2001
         call icectr(0, 0)

c RJP 04/08/07 added CO2 keywords
      else if (macro .eq. 'carb') then
c read problem type
c iprtype = 1, water only problem (2 DOFs)
c iprtype = 2, co2 only problem (2 DOFs)
c iprtype = 3, CO2-water problem, no solubility (3 DOFs)
c iprtype = 4, co2-water problem, w/ soubility (4 DOFs)
c iprtype = 5, co2-water-air w/ solubility (5 DOFs)

         read(inpt,*) iprtype
         icarb = 1
         call icectrco2(0,0)

      else if (macro .eq. 'imex') then
c**** implicit-explicit solution ****
         read(inpt,*) isplitts, explicit_factor

      else if (macro .eq. 'impf') then
c**** calculation of fractional variable changes ****
         impf=1
         read(inpt, *) delpt, deltt, delst, delat
         if(delpt.le.0.0d00) delpt=1.0d00
         if(deltt.le.0.0d00) deltt=1.0d00
         if(delst.le.0.0d00) delst=1.0d00
         if(delat.le.0.0d00) delat=1.0d00

      else if (macro .eq. 'init') then
c**** initial values ****
         read (inpt, *) pein, tin0, tin, gradnt, depcng,
     *        tin1, grad2, quad

      else if (macro .eq. 'iter') then
c**** iteration parameters ****
         read (inpt, *) g1, g2, g3, tmch, overf
         read (inpt, *) irdof, islord, iback, icoupl, rnmax
         if (ihead .eq. 1 .and. irdof .ne. 13) then
            if (iout .ne. 0) write(iout,*) "Warning: irdof = 13 for ",
     &           "head problem"
            if (iptty .gt. 0) write(iptty,*) 
     .           "Warning: irdof = 13 for head problem"
            irdof = 13
         end if
         if (rnmax .lt. zero_t) rnmax = 1.0d+06
         rnmax = rnmax * 60.0
         if(tmch.lt.0.0) then
            fdum1=-1.0
            tmch=abs(tmch)
         else
            fdum1=0.0
         endif
         if(irdof.lt.0) then
            overf = 1.0
         endif

      else if (macro .eq. 'isot') then
c**** isotropic geometric coeficients
         isox=1
         isoy=1
         isoz=1
 
      else if (macro .eq. 'itfc') then
c**** define interface factors  ****
         interface_flag = 1
         call initfc
 
      else if (macro .eq. 'itup') then
c**** itup values ****
         read(inpt, *)iad_up
 
      else if (macro .eq. 'iupk') then
         iupk = 1

      else if (macro .eq. 'ivfc') then
c**** enable vfcal routine ****
         ivfcal =  1

      else if (macro .eq. 'mdno') then
c**** multiply defined nodes ****
         imdnode = 1 
         call md_nodes(0,0,0)                

      else if (macro .eq. 'mptr') then
c**** particle tracking ****
c     inmptr is called in part_track instead of here

         read(inpt,*)nspeci,maxlayers,max_particles
	   read(inpt,*)pout,prnt_rst
           read(inpt,'(a4)') macro1
           if(macro1.eq.'tcur') then
              read(inpt,'(a80)') dummy_line
              read(inpt,'(a80)') dummy_line
           else
              backspace(inpt)
           endif
           read(inpt,'(a4)')macro1
           if(macro1.eq.'zptr')then
              read(inpt,'(a4)') macro1
              read(inpt,'(a4)') macro1
           else
              backspace (inpt)
           end if
	   read(inpt,*)rseed
c*** water table rise modification
           read(inpt,'(a4)')macro
           if(macro.eq.'wtri')then
             read(inpt,*) water_table
           else
             water_table=-1.d+10
             backspace (inpt)
           end if
c*** water table rise modification
	   read(inpt,*)daycs,daycf,dayhf,dayhs

c     The rest of these reads are simply to get to the bottom
c     of mptr macro. They are read in for real in the call
c     of inmptr in part_track

           read(inpt,'(a4)')macro1
           if(macro1.eq.'file')then
              read(inpt,'(a4)') macro1
           else
              backspace (inpt)
           end if
c     Handle keyword afm if it is there
           read(inpt,'(a4)')macro1
           if(macro1(1:3).ne.'afm')then
              backspace (inpt)
           end if
c     Read through layer info
c     maxlayers+1 to make sure we reach the blank line
           do j = 1, maxlayers+1
              read(inpt,'(a80)') wdd1
              if(null1(wdd1)) goto 1320
           end do
 1320      continue
 1321      continue
c     Read through itrc info
           read(inpt,'(a80)') wdd1
           if(null1(wdd1)) goto 1322
           goto 1321
 1322      continue

c     loop over each species

           do j = 1, nspeci
c     ith, trak_type(ith), ...
              read(inpt,'(a80)') wdd1
c	Check for keyword size, read past it if encountered
              read(inpt,'(a80)') wdd1
              if(wdd1(1:4).eq.'size') then
c	Check for existence of keyword file
                 read(inpt,'(a80)') wdd1
                 if(wdd1(1:4).eq.'file') then
c	Read in file name, keep going
                    read(inpt,'(a80)') wdd1
                 else
                    backspace(inpt)
 4323               continue
c     read individual lines of size distribution input
                    read(inpt,'(a80)') wdd1
                    if(null1(wdd1)) goto 4324
                    goto 4323
 4324               continue
                 end if
              else
                 backspace(inpt)
              end if

CHari 01-Nov-06 Read past colloid diversity model
              read(inpt,'(a80)') wdd1
              if(wdd1(1:4).eq.'dive') then
c daughter colloid flag
                 read(inpt,*) idum
                 if (idum .eq. 0) then
c if not file, the tprpflag
                    read(inpt,'(a80)')wdd1
                    if (wdd1(1:4) .eq. 'file') then
c read past filename
                       read(inpt,'(a80)')wdd1
                    end if
c if file the tprpflag, else the CDF equation data
                    read(inpt,'(a80)')wdd1
c reversible/irreversible flag
                    read(inpt,'(a80)')wdd1
                 end if
              else
                 backspace(inpt)
              endif

              read(inpt,*) idum
              do jj = 1, idum
c     layeri, transflag, ...
                 read(inpt,'(a80)') wdd1
              end do
              read(inpt,*) idum
              do jj = 1, idum
c     tmpcnsk
                 read(inpt,'(a80)') wdd1
 1323            continue
c     pinmass, ...
                 read(inpt,'(a80)') wdd1
                 if(null1(wdd1)) goto 1324
                 goto 1323
 1324            continue
                 
              end do
           end do

      else if (macro .eq. 'ngas' .or. macro .eq. 'co2i') then
c**** noncondensible gas ****
         ico2 = 3
         call co2ctr (0)

      else if (macro .eq. 'node' .or. macro .eq. 'nod2'
     &        .or. macro .eq. 'nod3')  then
c**** node information ****
         call innode (macro)
         
      else if (macro .eq. 'perm') then
c**** permeabilities ****
         if (idoff .eq. -1) then
            if (iout .ne. 0) write (iout, *) 'Permeability not read: ',
     .           'heat conduction only problem'
 345        read (inpt, '(a80)') wdd1
            if(null1(wdd1)) goto 543
            goto 345
         else
            call inperm
         end if
 543     continue
      else if (macro .eq. 'anpe') then
c**** anisotropic permeabilities (off-diagonal terms) ****
         ianpe = 1
         call inanpe

      else if (macro .eq. 'fper') then
c**** permeabilities ****
         call scaled_perm

      else if (macro .eq. 'ppor') then
c**** pressure and temperature dependent porosity and permeability ****
         iporos = 1
         call porosi(0)

      else if (macro .eq. 'pres') then
c**** pressure information ****
         call inpres

      else if (macro .eq. 'ptrk') then
c**** particle tracking ****
         call inptrk(simnum)
         
      else if (macro .eq. 'renu') then
c**** renumber ****
         ireord = 1
         call renum(0)
      else if (macro .eq. 'rflo') then
c     Do nothing - flag to read in flux values already set in scanin
      else  if (macro .eq. 'rflx') then
c**** read in radiation source term(w) ****
         read(inpt,*) emiss
         emiss = emiss * 5.6697e-8
         call inhflx(macro)

      else if (macro .eq. 'rlp ') then
c**** read in relative permeability information ****
         call rlperm (0,0)

      else if (macro .eq. 'frlp') then
c**** read in relative permeability factors ****
c**** these are minimum relative permeabilities *****
         call infrlp        
      else if (macro .eq. 'rich') then
         jswitch = 1
         if (idpdp .eq. 1) joff = 4
         read (inpt,*) strd_rich,tol_phase 

c**** rock densities, etc ****
      else if (macro .eq. 'rock') then
c**** rock densities, etc ****
         call inrock

      else if(macro .eq. 'rxn ') then
c**** Multiple reaction data ****
         rxn_flag = 1
         call read_rxn

      else if (macro .eq. 'intg') then
c**** set integration type (node or gaussian quadrature)
         read (inpt, *) intg

      else if (macro .eq. 'hcon') then
c**** set solution to heat conduction ****
         idof = 1
         idoff = -1
      else if (macro .eq. 'sol ') then
c**** solution and integration type ****
         read (inpt, *) idoff, intg
         if (idoff .le. 0) then
            idoff = -1
            idof = 1
         else
            idof = 2
         end if

      else if (macro .eq. 'solv') then
c**** solv (doesn't do anything yet) ****

******sptr***groemer 5/1/98***********************
      else if (macro .eq. 'sptr') then
c**** streamline particle tracker
         call insptr (wdd)
******sptr***groemer 5/1/98************************* 

c gaz 11-27-2001
      else if (macro .eq. 'stea') then
c**** Steady state solution  ****
         isteady = -1
         call steady(0, 0., 0.)
         
      else if (macro .eq. 'subm') then
c**** print out sub model boundary conditions
         isubbc = 1
         call submodel_bc(0)
         call submodel_bc(1)

      else if (macro .eq. 'wflo') then
c**** print out sub model boundary conditions
         isubbc = 2
         call submodel_bc(0)
         call submodel_bc(1)       

      else if (macro .eq. 'exuz') then
c**** explicit solution of UZ problems ****
c         read(inpt,'(a80)') wdd1
c         read(wdd1,*,end= 881) sssol, time_ss, tims_trans
c         go to 882
c 881     continue  
c         time_ss = 0.0
c         tims_trans = 0.0 
c 882     continue      
c         isteady = -1

      else if (macro .eq. 'strs') then
c**** stress solution data ****
         istrs = 1
         cnum_stress = cnum
         call stressctr (0,0)
         cnum = cnum_stress

      else if (macro .eq. 'szna' .or. macro .eq. 'napl') then
c**** isothermal NAPL - water transport ****
         read(inpt,*) ico2
         ico2 = -3
         call airctr (0, 0)

! zvd 20-Jan-2004 Revamp for reading and writing restart files
      else if (macro .eq. 'rest') then
c**** restart file I/O parameters
         call inrestart

      else if (macro .eq. 'text') then
c**** read text from input file ****
 140     continue
        read (inpt, '(a80)') wdd1
        if (.not. null1(wdd1)) then 
           if (iout .ne. 0) write(iout, '(a80)') wdd1
           go  to  140
        end if
        if (iout .ne. 0) write(iout, 6000) macro, inpt
        if (iptty .gt. 0) write(iptty, 6000) macro, inpt

      else if (macro .eq. 'time') then
c**** time step information ****
         call intime

      else if (macro .eq. 'thic') then
c**** thickness information ****
         ithic = 1
         call thickness(0)                  

      else if (macro .eq. 'trac') then
c**** tracer data ****
         iccen = 1
         call concen (0,0,dummyreal)

      else if (macro .eq. 'user') then
c**** call user subroutine during input with argument kk ****
         read (inpt, *) kk
         call user(kk)

      else if (macro .eq. 'vapl') then
c**** enable vapor pressure lowering ****
        ivapl=1

      else if (macro .eq. 'vcon') then
c**** read in nonlinear thermal conductivity information ****
         ivcond = 1
         call vcon (0, 0)
         
      else if (macro .eq. 'vbou') then
c**** read in nonlinear thermal conductivity information ****
         ivcond = 1
         call vboun (0, 0)

      else if (macro .eq. 'dvel') then
c**** calculate intermode velocities ****
         ivelo = -1
         call flxo(0)

      else if (macro .eq. 'pest') then
c**** output for parameter estimation routine PEST ****
         ipest = 1
         call pest(0)
         
      else if (macro .eq. 'rive'.or. macro .eq. 'well') then
c RJP 12/13/06 modification to read wellbore/river specific data
c**** river and wellbore model input ****

         call river_ctr(0)
         call river_ctr(1)
         call river_ctr(5)
         read (inpt, '(a80)') wdd1

         if (.not. null1(wdd1)) backspace inpt

      else if (macro .eq. 'phys'.or. macro .eq. 'ndph') then
c GAZ 11/04/08 modification to read wellbore/river specific data
c**** wellbore physics or non darcy flow

         call wellphysicsctr(0, 0)

      else if (macro .eq. 'weli') then
c detailed well impedance calculations (peaceman-type well)
       iwelimd = 1
       call wellimped_ctr(0)
       
      else if (macro .eq. 'wlbr') then
c**** wellbore input ****
         iwelb = 1
         write(ierr,*) 'This option (welbor) not supported.'
         write(ierr,*) 'Stop in input'
         stop
c        call welbor (0)

c**** return to main program because all data is read ****
         go to 210

      else if (macro .eq. 'zeol') then
c**** Zeolite water balance input
         call inzeol(0)

      else if (macro .eq. 'zone') then
c**** set zone information ****
         cnum = cnum + 1
         call zone(cnum, inpt)

      else if (macro .eq. 'zonn') then
c**** set zone information ****
         cnum = cnum + 1
         call zone(cnum, inpt)

      else
c**** error occurred ****
         write(ierr, 6100) macro
         if (iout .ne. 0) write(iout, 6100) macro
         if (iptty .gt. 0) write(iptty, 6100) macro
 6100    format(/, 1x, '**** error in input deck  : ', a4,' ****', /)
         stop
      endif

      if (inpt .ne. inptorig) then
         call done_macro (inpt)
      end if

      goto  100

 200  continue
      inpt = inptorig

c**** set n = neq for standard problem, n = 3*neq for dual porosity ****
      n = neq

c**** if both dual and dpdp are enabled, disable dual ****
      if(idualp .ne. 0 .and. idpdp .ne. 0) then
         idualp = 0
         write (ierr, 6110)
         if (ischk .ne. 0) write (ischk, 6110)
 6110    format(/, 'Both dual and dpdp options were enabled,', /,
     .        'dual has been disabled, dpdp will be used') 
      endif

      call dual (3)
      call dpdp (6)
c**** check if gravity is present and modify head solution as necessary
      call headctr(0,0,0.0,0.0)
c**** check if air macro called if head macro called
      call headctr(1,0,0.0,0.0)

 210  continue
      call steady(-1,0.,0.)
      end








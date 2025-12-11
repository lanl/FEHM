      subroutine rdcon(iz)
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
CD1 To read the tracer information from the input file, initialize
CD1 tracer mass balance terms.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2                                     However, previous non-YMP
CD2                                     versions of FEHM exist, and
CD2                                     the current version differs
CD2                                     from these in minor ways.  
CD2
CD2 $Log:   /pvcs.config/fehm90/src/rdcon.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:40   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:56   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:22   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:38 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.23   Thu Jun 13 08:05:34 1996   hend
CD2 Added dummy argument for concadiff for consistency
CD2 
CD2    Rev 1.22   Wed Jun 12 15:21:14 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.21   Mon Jun 03 11:18:18 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.20   Fri May 31 11:01:18 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.19   Wed May 29 15:07:22 1996   hend
CD2 Added variable diffusion with water content
CD2 
CD2    Rev 1.18   Fri May 24 09:55:34 1996   hend
CD2 Updated trac for mdnodes
CD2 
CD2    Rev 1.17   Tue May 21 10:37:50 1996   hend
CD2 Updated for use with both ldsp and dspl,dspv,or dspb
CD2 
CD2    Rev 1.16   Wed May 08 13:34:56 1996   hend
CD2 Updated trac time step manipulations
CD2 
CD2    Rev 1.15   Thu Apr 25 13:34:58 1996   hend
CD2 Updated for use in long/trans dispersion option
CD2 
CD2    Rev 1.14   Thu Mar 21 13:18:54 1996   hend
CD2 Fixed to use indexing independent of istrw
CD2 
CD2    Rev 1.13   Mon Mar 04 16:12:06 1996   hend
CD2 Removed uneccesary calculations from coneq1 and added trac input option
CD2 
CD2    Rev 1.12   Fri Feb 02 10:22:08 1996   hend
CD2 Added to Requirements Traceability
CD2 
CD2    Rev 1.11   Thu Feb 01 15:31:16 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.10   Wed Jan 17 14:13:18 1996   hend
CD2 Added use of parser for input lines to function on sgi
CD2 
CD2    Rev 1.9   12/13/95 10:30:04   robinson
CD2 Alpha is changed to alpha**2, which is now used in coneq1
CD2 
CD2    Rev 1.8   08/18/95 10:36:12   llt
CD2 inpt already defined, removed for cray
CD2 
CD2    Rev 1.7   04/25/95 09:59:04   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.6   04/10/95 11:29:38   robinson
CD2 Code can now read keyword for solute user sub
CD2 
CD2    Rev 1.5   01/28/95 14:20:34   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.4   06/20/94 11:14:30   zvd
CD2 Added ierr unit number for error output.
CD2
CD2    Rev 1.3   04/04/94 16:05:56   robinson
CD2 Corrected input problem for solid species.
CD2 
CD2    Rev 1.2   03/18/94 16:15:48   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/16/94 09:20:46   zvd
CD2 Fixed bug, defined constant to set minimum dispersivity value.
CD2 
CD2    Rev 1.0   01/20/94 10:26:52   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier              Type     Use  Description
CD3
CD3 iz                      int       I   Parameter used to bypass the
CD3                                       entire routine if desired
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name           Description
CD3
CD3 input file     Contains all input data from FEHMN macros in ASCII
CD3                form
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 inpt, an0, awc, epc, upwgta, n0, an, daycs, daycf, dayhf, dayhs,
CD4 iaccmx, daycm, daycmm, daycmx, daymin, nspeci, nsp, qcout, qcin, qcrxn,
CD4 npn, npt, icns, tclx, tcly, tclz, tcvx, tcvy, tcvz, a1adfl, a2adfl,
CD4 betadfl, a1adfv, a2adfv, betadfv, wdd1, iadsfl, iadsfv, diffmfl,
CD4 diffmfv, narrays, pointer, itype, default, igroup, ireturn,
CD4 macroread, zero_t
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
CD4 exit    N/A      Terminates execution of the program after closing
CD4                  all files
CD4 null1    N/A      Checks for null line or zero as next input line
CD4 plotc1  N/A      Writes tracer information to the output tracer file
CD4 thermc  N/A      Computes the tracer storage, sorption, and
CD4                  reaction terms of the mass balance
CD4 initdata
CD4         N/A      Reads and sets input values defined at each node
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop parameter
CD5 ij           int         Do loop index parameter
CD5 jj           int         Do loop index parameter
CD5 iret         int         Value to determine if execution is to be
CD5                          terminated within this routine
CD5 ntpp         int         Number of spaces alloted in concentration
CD5                          arrays for each tracer
CD5 nsp          int         Do loop index parameter
CD5 n_node_sets  int         Number of sets of nodes
CD5 i_node_set   int         Current node set
CD5 mi           int         Current position in concentration array
CD5 mim          int         Current node number
CD5 ndummy       int         Parameter directing code to correct
CD5                              position in concentration array
CD5 npt_subst    int         Parameter directing code to correct
CD5                              position in concentration array
CD5 ispecies     int         Do loop index over all species
CD5 tcxd         real*8      Dispersion coefficient in the x-direction
CD5 tcyd         real*8      Dispersion coefficient in the y-direction
CD5 tczd         real*8      Dispersion coefficient in the z-direction
CD5                              liquid concentration
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6 First, the code determines if this data is to be read, and if not,
CD6 bypasses the section that reads the data.  Otherwise, the first
CD6 line of input is read.  If the initial concentration field is not
CD6 to be read from a file, then the concentration is initialized at
CD6 each position to the input value (which may be changed later if
CD6 the user chooses to input initial concentration node by node).
CD6 Next, the code checks for the reasonablenes of values of
CD6 implicitness factor and upweighting factors, chnaging the input if
CD6 necessary.  
CD6
CD6 Next, the second line is read, and the times input are changed if
CD6 they are inconsistent.  The third input line is then read, and the
CD6 time steps are set based on the input.  Then, the fourth line is
CD6 read (number of tracers).  The the code loops through for each
CD6 tracer and initializes parameter values.  The next two lines,
CD6 representing the sorption and reaction input parameters, are read
CD6 in.  Then, an IF-ELSEIF block determines which type of isotherm is
CD6 specified and automatrically resets sorption model parameters to
CD6 reflect the chosen sorption model.
CD6
CD6 Next, a REPEAT-UNTIL block reads in the dispersion coefficient
CD6 input.  The values are set to their minimum values if they are too
CD6 low.  Then, the code determines from the input whether it is by
CD6 zone or by node.  If by zone, the code loops through each node, and
CD6 sets the dispersion coefficient if the node is in the current zone.
CD6  If by node, the code loops through each of the specified nodes and
CD6 sets the dispersion coefficient values.  This block is exited when
CD6 a value of zero is encountered for the first entry in the line
CD6 specifying the node or zone to be handled is 0.
CD6
CD6 Next, a REPEAT-UNTIL block identical in structure to the previous
CD6 one is used to input initial concentration values that are input
CD6 by node or zone.  These values overwrite the initial
CD6 concentrations input above for the entire domain.  A third REPEAT-
CD6 UNTIL block then reads and sets the inlet concentration and
CD6 injection time information.  The outer loop for each tracer is
CD6 then exited after these operations are performed for each tracer.
CD6 This concludes the section of the code that reads the data.
CD6
CD6 Next, if this call to the routine is to initialize arrays, the
CD6 code checks to see if there is enough memory allocated, and if
CD6 there is not, ERROREXITS (note that the remainder of the routine
CD6 is bypassed if this is not the initialization step).  Otherwise,
CD6 the code loops through each tracer and calls thermc to set the
CD6 initial storage, source/sink, and reaction terms of the tracer
CD6 mass balance.  The code checks to see if this is a dual porosity
CD6 simulation, and if it is, calls thermc twice more in succession to
CD6 initialize the two groups of nodes within the matrix.  Then, the
CD6 code loops through each node and sets other mass balance arrays to
CD6 their initial values.  Finally, the code calls plotc1 to write the
CD6 initial information to the tracer output file.  Then, in the
CD6 ERRORSEGMENT section, the code checks to see if this point has
CD6 been reached because of an error and if this is true, calls exit
CD6 to terminate execution.  Otherwise, the code returns to the
CD6 calling routine.
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
CD9 2.3.4 Solute-transport equations
CD9 2.4.5 Adsorbing solutes
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN rdcon
CPS 
CPS IF the data are to be read in on this pass
CPS
CPS   Read in first line of input
CPS
CPS   IF the implicitness factor is higher than the highest...
CPS   ... acceptable value
CPS     Set implicitness factor to highest value
CPS   ENDIF
CPS
CPS   IF the implicitness factor is lower than the lowest acceptable...
CPS   ...  value
CPS     Set implicitness factor to lowest value
CPS   ENDIF
CPS
CPS
CPS   IF the upweighting factor is higher than the highest acceptable...
CPS   ...  value
CPS     Set upweighting factor to highest value
CPS   ENDIF
CPS
CPS
CPS   IF the upweighting factor is lower than the lowest acceptable...
CPS   ...  value
CPS     Set upweighting factor to lowest value
CPS   ENDIF
CPS
CPS
CPS   Read in the second line of input
CPS
CPS   IF the time when the heat and mass solution is enabled is less...
CPS   ...  than the time when it was to be disabled
CPS     Set the two time equal
CPS   ENDIF
CPS
CPS   Read in the third line of input
CPS   Set values for time steps based on the input
CPS   Set time step and implicitness factor
CPS   Read in the number of species and set pointers
CPS
CPS   FOR each tracer
CPS     Initialize parameter values for each tracer
CPS   ENDFOR
CPS
CPS   FOR  each tracer
CPS
CPS     Set tracer partition number
CPS     Read in phase designation for this tracer
CPS     IF the current tracer is solid
CPS       Initialize parameters
CPS     ELSE the tracer is vapor or liquid
CPS       FOR each region to be read
CPS         null1 - Check for null line or zero
CPS         IF there is another region to read
CPS           IF this is a liquid species
CPS             Read in sorption and dispersion data for this tracer
CPS           ELSEIF this is a vapor species
CPS             Read in sorption and dispersion data for this tracer
CPS           ELSEIF this is a Henry's Law  species
CPS             Read in sorption and dispersion data for this tracer
CPS           ENDIF
CPS           IF this is a vapor species
CPS             Handle special case of a vapor only species by...
CPS             ... setting parameter values
CPS           ENDIF
CPS           Set parameters to their limiting values when they are...
CPS           ... outside their ranges
CPS           IF this is a conservative liquid tracer 
CPS             Set sorption parameters to reflect a conservative tracer
CPS           ELSEIF this is a liquid linear isotherm
CPS             Set sorption parameters to reflect a linear isotherm
CPS           ELSEIF this is a liquid Freundlich isotherm
CPS             Set sorption parameters to reflect a Freundlich isotherm
CPS           ELSEIF this is a liquid Langmuir isotherm
CPS             Set sorption parameters to reflect a Langmuir isotherm
CPS           ENDIF
CPS           IF this is a conservative vapor tracer
CPS             Set sorption parameters to reflect a conservative tracer
CPS           ELSEIF this is a vapor linear isotherm
CPS             Set sorption parameters to reflect a linear isotherm
CPS           ELSEIF this is a vapor Freundlich isotherm
CPS             Set sorption parameters to reflect a Freundlich isotherm
CPS           ELSEIF this is a vapor Langmuir isotherm
CPS             Set sorption parameters to reflect a Langmuir isotherm
CPS           ENDIF
CPS         EXITIF there is no more data to read
CPS         ENDIF
CPS       ENDFOR each region to be read
CPS       
CPS       Set parameters for reading region numbers
CPS     
CPS       initdata - read region numbers and set at each node
CPS     
CPS       IF the solute is a Henry's Law species
CPS         Set the Henry's Law parameters
CPS       ENDIF
CPS       
CPS     ENDIF
CPS     
CPS     Set parameters for reading initial concentrations
CPS     initdata - read concentrations and set at each node
CPS       
CPS     FOR each node
CPS       Initialize concentration value for previous time step
CPS     ENDFOR
CPS       
CPS     Set parameters for inlet concentrations
CPS     initdata - read inlet concentrations and set at each node
CPS       
CPS   ENDFOR
CPS   
CPS   Set flag to denote that the trac macro has been called
CPS       
CPS   
CPS ENDIF the data are to be read in on this pass
CPS 
CPS IF this call to the routine is for initialization
CPS
CPS   Compute size of concentration array
CPS   
CPS   FOR each species
CPS     IF this species is a solid
CPS       FOR each unknown number
CPS         Set the region parameter to 1
CPS       ENDFOR
CPS     ENDIF
CPS   ENDFOR
CPS
CPS   IF there is not enough memory for the problem specified
CPS     Write error message
CPS     ERROREXIT
CPS   ENDIF
CPS   
CPS   FOR each tracer
CPS 
CPS     Initialize mass balance parameter values
CPS     
CPS     IF this solute is a Henry's Law species with vapor...
CPS     ... concentration specified
CPS     
CPS       IF this is a dpdp simulation
CPS         Set number of node sets to 2
CPS       ELSEIF it is dual porosity
CPS         Set number of node sets to 3
CPS       ELSE
CPS         Set number of node sets to 1
CPS       ENDIF
CPS   
CPS       FOR each node set
CPS         FOR each node
CPS           Compute liquid concentration corresponding to the input...
CPS           ... vapor concentration
CPS           IF there is no transport to the liquid phase
CPS             Write error message to the screen
CPS             Set flag to stop program
CPS             ERROREXIT
CPS           ENDIF
CPS         ENDFOR each node
CPS       ENDFOR each node set
CPS   
CPS     ENDIF
CPS     
CPS     thermc - set initial tracer storage and reaction terms in...
CPS     ...    mass balance
CPS
CPS     IF this is a dpdp simulation
CPS       thermc - set storage and reaction terms for the matrix nodes
CPS     ELSEIF is is a dual porosity simulation
CPS       thermc - set storage and reaction terms for the first set...
CPS       thermc - set storage and reaction terms for the second set...
CPS       ...      of matrix nodes
CPS     ENDIF
CPS
CPS     FOR each node
CPS       Set initial values of other mass balance arrays
CPS     ENDFOR
CPS     
CPS   ENDFOR each tracer
CPS   
CPS   plotc1 - write initial information to tracer output file
CPS ENDIF this call to the routine is for initialization
CPS
CPS ERRORSEGMENT
CPS   IF execution is to be terminated
CPS     Terminate program
CPS   ENDIF
CPS ENDSEGMENT
CPS END rdcon
C**********************************************************************

      use comflow
      use comcouple
      use comrxni
      use comfi
      use comdi
      use comci
      use comchem
      use combi
      use comdti
      use comai
      use comki
      use compart
      use davidi, only : irdof
	use trxnvars
      implicit none

      real*8 h_const, dvap_conc
      integer mi,mim,ndummy,iret,i,iz,jj,ij,ispecies,npt_subst
      integer n_node_sets,i_node_set, ja, jb, jc
      logical null1, readflag, alt_zone, found, mole_input
      character*5 user_macro
      character*4 zmacro
      character*3 nstring
c     SPC
      character*3 string3
c gaz 040420 
      character*4 back_trac
      integer iback_trac
c      
      character*80 input_msg
      character*100 zone_file
      integer msg(16)
      integer nwds
      real*8 xmsg(16)
      integer imsg(16)
      character*32 cmsg(16)
      integer kz,iz2,iz2p,i1,i2,ipivt,ii,sehindexl,sehindexv,neqp1
      integer izf, nin, num_nodes, nend, znum
      integer, allocatable :: node_list(:)
      real*8 tempx,tempy,tempz,templength
      real*8 cord1x,cord1y,cord2x,cord2y
      real*8 dispxavw,dispyavw,dispzavw,cord1z,cord2z,x
      integer icpnt,iimm,ivap
      integer idum, io_stat, open_file, mfile
      real*8 rdum1, rdum2, mvol, psdum, roldum, sdum, sctmp, frac
      real*8 ctmp, ctol, andum, antmp, anvtmp, mvolv
      real*8, allocatable :: rvol(:), lvol(:), vvol(:)
      integer, allocatable :: hflag(:)
      logical, allocatable :: seen(:)
	logical opend

      save mfile, mole_input
 	  do i=1,nspeci
	  if(species(i).eq.' ') write(species(i),'(i2)') i
	  end do
      ctol = 1.d-06
      iret = 0
      if(iz.eq.0) then
         if(trxn_flag .eq. 0) then
c     
c     input tracer data
c     
c     Hari  6/10/04  -------------
            read(inpt,'(a3)') string3
            if(string3(1:3) .eq. 'rip') then
               read(inpt,*) ripfehm
            else
               backspace inpt
            endif
c     Hari ----------------------------------
            read(inpt,'(a5)') user_macro
            if(user_macro(1:5) .eq. 'userc' ) then
               backspace inpt
               read(inpt,'(a80)') input_msg
               call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
               if (msg(2).eq.1) then
                  an0=imsg(2)
               else
                  an0=xmsg(2)
               endif
               if (msg(3).eq.1) then
                  awc=imsg(3)
               else
                  awc=xmsg(3)
               endif
               if (msg(4).eq.1) then
                  epc=imsg(4)
               else
                  epc=xmsg(4)
               endif
               if (msg(5).eq.1) then
                  upwgta=imsg(5)
               else
                  upwgta=xmsg(5)
               endif
c     zvd 04/26/07 Moves call to userc, after first line parameters are read
c     This way an optional file name for the userc data file can be read
c     and we only have to backspace once
               call userc(0, idum, rdum1, rdum2)
            else
               backspace inpt
               read(inpt,*) an0,awc,epc,upwgta
            end if

            if(awc.le.1.0) then
               awc=1.0
            end if

            if(awc.gt.1.0) then
               awc=0.5
            end if

            if(upwgta.gt.1.0) then
               upwgta=1.0
            end if

            if(upwgta.lt.0.5) then
               upwgta=0.5
            end if

            an = an0

            read(inpt,*) daycs,daycf,dayhf,dayhs
            if(dayhs.lt.dayhf) then
               dayhs=dayhf
            end if
c     Revised input to handle optional print out variable BAR 11-18-98
c     read(inpt,*) iaccmx,daycm,daycmm,daycmx
            read(inpt,'(a80)') input_msg
            call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)

            iaccmx = imsg(1)

            if(msg(2).eq.1) then
               daycm = imsg(2)
            else
               daycm = xmsg(2)
            end if
            
            if(msg(3).eq.1) then
               daycmm = imsg(3)
            else
               daycmm = xmsg(3)
            end if
            
            if(msg(4).eq.1) then
               daycmx = imsg(4)
            else
               daycmx = xmsg(4)
            end if

            if(nwds.gt.4) then
               nprttrc = imsg(5)
            else
               nprttrc = 1
            end if
            
            daycmx=min(daycmx,daymax)
            dtotc=daycmm * 86400.
            ayc=1.0-awc

c     read in optional tracer porosity
            read(inpt,'(a4)') input_msg
            if (input_msg(1:4).eq.'tpor') then
c     zvd 06/02/03 Set default porosity to bogus value, assign all unassigned
c     nodes the value set for rock macro after all data has been read 
c     ps_trac(1:n0)=ps(1:n0)
               narrays=1
               itype(1) = 8
               macro = "trac"
               igroup = 5
               default(1) = 10.
               call initdata2( inpt, ischk, n0, narrays, itype,
     2              default, macroread(5), macro, igroup, ireturn,
     3              r8_1=ps_trac(1:n0) )
            else
               backspace inpt
c     zvd 06/02/03 Moved to end of input check, assign tracer porosities
c     after all data has been read
c     ps_trac(1:n0) = ps(1:n0)
            endif
c     read in number of species and set pointers
            read(inpt,*) nspeci
            if(rxn_flag.eq.0)then
               ncpntprt=ncpnt
               do ii=1,ncpntprt
                  cpntprt(ii)=ii
                  cpntnam(ii)(1:16) = "Aqueous_Species_"
                  write (nstring, '(i3.3)') ii
                  cpntnam(ii)(17:19) = nstring
               enddo
               ncplxprt=ncplx
               do ii=101,ncplxprt+100
                  cplxprt(ii)=ii
                  cplxnam(ii)(1:16) = "Aqueous_Complex_"
                  write (nstring, '(i3.3)') ii
                  cplxnam(ii)(17:19) = nstring
               enddo
               nimmprt=nimm
               do ii=1,nimmprt
                  immprt(ii)=ii
                  immnam(ii)(1:17) = "Immobile_Species_"
                  write (nstring, '(i3.3)') ii
                  immnam(ii)(18:20) = nstring
               enddo
               nvapprt=nvap
               do ii=1,nvapprt
                  vapprt(ii)=ii
                  vapnam(ii)(1:14) = "Vapor_Species_"
                  write (nstring, '(i3.3)') ii
                  vapnam(ii)(15:17) = nstring
               enddo
            endif
            icpnt = 0
            iimm = 0
            ivap = 0

c     
c     partition memory for multiple tracers
c     
            ntpp=n7/nspeci
            do nsp=1,nspeci
               npt(nsp)=(nsp-1)*ntpp
c     zero out species mass balance parameters
               qcout(nsp)=0.0
               qcin(nsp)=0.0
               qcrxn(nsp)=0.0
            enddo

c     Check for keyword specifying to use longitudinal and transverse 
c     dispersion coefficients, seh
            read(inpt,'(a4)') input_msg
            if (input_msg(1:4).eq.'ldsp') then
               ldsp=1
            else
               ldsp=0
               backspace inpt
            endif

c     seh
C     Check for keyword specifying to use same dispersivity and diffusivity
C     for all liquid / vapor calculations
            read(inpt,'(a4)') input_msg
            dispsame=0
            hvliquid=0
            hvvapor=0
            if (input_msg(1:3).eq.'dsp') then
               dispsame = 1
               jj = 0
               if (ldsp .eq. 0) then
                  do
                     read(inpt,'(a80)') wdd1
                     if (null1(wdd1)) exit
                     backspace inpt
                     jj = jj + 1
                     call parse_string (wdd1, imsg, msg, xmsg, cmsg, 
     &                    nwds)

                     select case (input_msg(1:4))
                  case ('dspl')
                     if (msg(1) .eq. 1) then
                        read (inpt,*,err=2000,iostat=io_stat) 
     &                       mflagl(1,jj), sehdiff(jj),
     &                       tclx(1,jj), tcly(1,jj), tclz(1,jj)
                     else
c     Try reading using old format
                        read (inpt,*,err=2000,iostat=io_stat) 
     &                       sehdiff(jj),
     &                       tclx(1,jj), tcly(1,jj), tclz(1,jj)
                        mflagl(1,jj) = 0
                     end if
                  case('dspv')
                     if (msg(1) .eq. 1) then
                        read (inpt,*,err=2000,iostat=io_stat) 
     &                       mflagv(1,jj), sehdiffv(jj),
     &                       tcvx(1,jj), tcvy(1,jj), tcvz(1,jj)  
                     else
c     Try reading using old format
                        read (inpt,*,err=2000,iostat=io_stat) 
     &                       sehdiffv(jj),
     &                       tcvx(1,jj), tcvy(1,jj), tcvz(1,jj)  
                        mflagv(1,jj) = 0
                     end if
                  case('dspb')
                     if (msg(1) .eq. 1) then
                        read (inpt,*,err=2000,iostat=io_stat) 
     &                       mflagl(1,jj), sehdiff(jj),
     &                       tclx(1,jj), tcly(1,jj), tclz(1,jj),
     &                       mflagv(dispsame,jj), sehdiffv(jj),
     &                       tcvx(1,jj), tcvy(1,jj), tcvz(1,jj)
                     else
c     Try reading using old format
                        read (inpt,*,err=2000,iostat=io_stat) 
     &                       sehdiff(jj),
     &                       tclx(1,jj), tcly(1,jj), tclz(1,jj),
     &                       sehdiffv(jj),
     &                       tcvx(1,jj), tcvy(1,jj), tcvz(1,jj)
                        mflagl(1,jj) = 0
                        mflagv(1,jj) = 0
                     end if                 
                  end select
               end do
            else ! If ldsp != 0
               do
                  read(inpt,'(a80)') wdd1
                  if (null1(wdd1)) exit
                  backspace inpt
                  jj = jj + 1
                  call parse_string (wdd1, imsg, msg, xmsg, cmsg, 
     &                 nwds)
                  
                  select case (input_msg(1:4))
               case ('dspl')
                  if (msg(1) .eq. 1) then
                     read (inpt,*,err=2000,iostat=io_stat) 
     &                    mflagl(1,jj), sehdiff(jj),
     &                    tclx(1,jj), tcly(1,jj)
                  else
c     Try reading using old format
                     read (inpt,*,err=2000,iostat=io_stat) 
     &                    sehdiff(jj),
     &                    tclx(1,jj), tcly(1,jj)
                     mflagl(1,jj) = 0
                  end if                     
               case('dspv')
                  if (msg(1) .eq. 1) then
                     read (inpt,*,err=2000,iostat=io_stat) 
     &                    mflagv(1,jj), sehdiffv(jj),
     &                    tcvx(1,jj), tcvy(1,jj)
                  else
c     Try reading using old format
                     read (inpt,*,err=2000,iostat=io_stat) 
     &                    sehdiffv(jj),
     &                    tcvx(1,jj), tcvy(1,jj)
                     mflagv(1,jj) = 0
                  end if                     
               case('dspb')
                  if (msg(1) .eq. 1) then
                     read (inpt,*,err=2000,iostat=io_stat) 
     &                    mflagl(1,jj), sehdiff(jj),
     &                    tclx(1,jj), tcly(1,jj),
     &                    mflagv(dispsame,jj), sehdiffv(jj),
     &                    tcvx(1,jj), tcvy(1,jj)
                  else
c     Try reading using old format
                     read (inpt,*,err=2000,iostat=io_stat) 
     &                    sehdiff(jj),
     &                    tclx(1,jj), tcly(1,jj),
     &                    sehdiffv(jj),
     &                    tcvx(1,jj), tcvy(1,jj)
                     mflagl(1,jj) = 0
                     mflagv(1,jj) = 0
                  end if                 
               end select
            end do
         end if ! if ldsp == 0
         if ((sehdiff(jj).lt.0.) .or. (sehdiffv(jj).lt.0.)) then
            write (ierr, *) 'Old Conca model not supported'
            write (ierr,*)'See UM for new diffusion model options'
            goto 2000
         end if
         narrays=1
         itype(1)=4
         default(1)=1
         macro="trac"
         igroup=10
         call initdata2(inpt,ischk,n0,narrays,itype,default,
     &        macroread(5),macro,igroup,ireturn,
     &        i4_1=itrcdsp(1:n0)) 
      else ! If we are not using dispsame
         backspace inpt
      end if

c     
c     read in data for each species
c     
c     Define maximum saturation for vapor or Henry's species,
c     currently only a single value can be defined for all  
c     vapor/Henry's species so it only needs to be entered once
      strac_max = 0.99
      mole_input = .false.
      do nsp=1,nspeci
         npn=npt(nsp)
c     read(inpt,*) icns(nsp)
         read(inpt,'(a80)') input_msg
         call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
         icns(nsp) = imsg(1)

c     construct input necessary for chemod and read_rxn
         if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2)then
            icpnt = icpnt + 1
            pcpnt(icpnt)=nsp
            if (nwds .ge. 2) then
c     check to see if this is a component name
               if (rxn_flag .eq. 0 .and. msg(2) .eq. 3)  
     &              cpntnam(icpnt) = cmsg(2)
c     check to see if the maximum saturation is specified for Henry's
               if(abs(icns(nsp)).eq.2) then 
                  if (nwds .eq. 2 .and. msg(2) .eq. 2) then
                     strac_max = xmsg(2)
                  else if (nwds .eq. 3 .and. msg(3) .eq. 2) then
                     strac_max = xmsg(3)
                  end if
               end if
            end if
            if (cden_flag .eq. 2) then
               if (cpntnam(icpnt) .eq. cden_spnam) ispcden = nsp
            end if
         elseif(icns(nsp).eq.0)then
            iimm = iimm + 1
            pimm(iimm)=nsp
            if (nwds .ge. 2 .and. rxn_flag .eq. 0 .and. 
     &           msg(2) .eq. 3) immnam(iimm) = cmsg(2)
         else
            ivap = ivap + 1
            pvap(ivap)=nsp
            if (nwds .ge. 2) then
c     check to see if this is a component name
               if (rxn_flag .eq. 0 .and. msg(2) .eq. 3) 
     &              vapnam(ivap) = cmsg(2)
c     check to see if the maximum saturation is specified for vapor
               if (nwds .eq. 2 .and. msg(2) .eq. 2) then
                  strac_max = xmsg(2)
               else if (nwds .eq. 3 .and. msg(3) .eq. 2) then
                  strac_max = xmsg(3)
               end if
            end if
         endif

C     PS     IF the current tracer is solid
         if(icns(nsp).eq.0)then                
C     PS       Initialize some parameters
            tclx(nsp,1)=0
            tcly(nsp,1)=0
            tclz(nsp,1)=0
            a1adfl(nsp,1)=0
            a2adfl(nsp,1)=0
            betadfl(nsp,1)=1
            tcvx(nsp,1)=0
            tcvy(nsp,1)=0
            tcvz(nsp,1)=0
            a1adfv(nsp,1)=0
            a2adfv(nsp,1)=0
            betadfv(nsp,1)=1
C     PS     ELSE the tracer is a liquid
         else
            jj = 0
            do 
               read(inpt,'(a80)') wdd1
               if(null1(wdd1)) exit
               backspace inpt
               jj = jj + 1
               call parse_string (wdd1, imsg, msg, xmsg, cmsg, 
     &              nwds)

               select case (icns(nsp))
            case (1)
               hvliquid=1
               if (dispsame.eq.1) then
                  read(inpt,*) iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                 a2adfl(nsp,jj),betadfl(nsp,jj)
               else
                  if (ldsp.eq.0) then
                     if (msg(5) .eq. 1) then
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       mflagl(nsp,jj),diffmfl(nsp,jj),
     &                       tclx(nsp,jj),tcly(nsp,jj),
     &                       tclz(nsp,jj)
                     else
c     Try reading using old format
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       diffmfl(nsp,jj),tclx(nsp,jj),
     &                       tcly(nsp,jj),tclz(nsp,jj)
                        mflagl(nsp,jj) = 0
                     end if
                  else
                     if (msg(5) .eq. 1) then
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       mflagl(nsp,jj),diffmfl(nsp,jj),
     &                       tclx(nsp,jj),tcly(nsp,jj)
                     else
c     Try reading using old format
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       diffmfl(nsp,jj),
     &                       tclx(nsp,jj),tcly(nsp,jj)
                        mflagl(nsp,jj) = 0
                     end if
                  endif
               endif
            case (-1)
C     PS         ELSEIF this is a vapor species
               hvvapor=1
               if (dispsame.eq.1) then
                  read(inpt,*) iadsfv(nsp,jj),a1adfv(nsp,jj),
     2                 a2adfv(nsp,jj),betadfv(nsp,jj)
               else
                  if (ldsp.eq.0) then
                     if (msg(5) .eq. 1) then
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     &                       a2adfv(nsp,jj),betadfv(nsp,jj),
     &                       mflagv(nsp,jj),diffmfv(nsp,jj),
     &                       tcvx(nsp,jj),tcvy(nsp,jj),
     &                       tcvz(nsp,jj)
                     else
c     Try reading using old format
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     &                       a2adfv(nsp,jj),betadfv(nsp,jj),
     &                       diffmfv(nsp,jj),tcvx(nsp,jj),
     &                       tcvy(nsp,jj),tcvz(nsp,jj)
                        mflagv(nsp,jj) = 0
                     end if
                  else
                     if (msg(5) .eq. 1) then
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     &                       a2adfv(nsp,jj),betadfv(nsp,jj),
     &                       mflagv(nsp,jj),diffmfv(nsp,jj),
     &                       tcvx(nsp,jj),tcvy(nsp,jj)
                     else
c     Try reading using old format
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     &                       a2adfv(nsp,jj),betadfv(nsp,jj),
     &                       diffmfv(nsp,jj),
     &                       tcvx(nsp,jj),tcvy(nsp,jj)
                        mflagv(nsp,jj) = 0
                     end if
                  endif
               endif
C     PS   IF vapor species handle special case of a vapor only species
               a1adfl(nsp,jj) = a1adfv(nsp,jj)
               a2adfl(nsp,jj) = a2adfv(nsp,jj)
               betadfl(nsp,jj) = betadfv(nsp,jj)
               iadsfl(nsp,jj)= iadsfv(nsp,jj)
C     PS   ENDIF
            case (2, -2)
C     PS         ELSEIF this is a Henry's Law  species
               hvliquid=1
               hvvapor=1
               if (dispsame.eq.1) then
                  read(inpt,*) iadsfl(nsp,jj),a1adfl(nsp,jj),
     2                 a2adfl(nsp,jj),betadfl(nsp,jj),
     3                 iadsfv(nsp,jj),a1adfv(nsp,jj),
     5                 a2adfv(nsp,jj),betadfv(nsp,jj)
               else
                  if (ldsp.eq.0) then
                     if (msg(5) .eq. 1) then
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       mflagl(nsp,jj),diffmfl(nsp,jj),
     &                       tclx(nsp,jj),tcly(nsp,jj),
     &                       tclz(nsp,jj),iadsfv(nsp,jj),
     &                       a1adfv(nsp,jj),a2adfv(nsp,jj),
     &                       betadfv(nsp,jj),mflagv(nsp,jj),
     &                       diffmfv(nsp,jj),tcvx(nsp,jj),
     &                       tcvy(nsp,jj),tcvz(nsp,jj)
                     else
c     Try reading using old format
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       diffmfl(nsp,jj),tclx(nsp,jj),
     &                       tcly(nsp,jj),tclz(nsp,jj),
     &                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     &                       a2adfv(nsp,jj),betadfv(nsp,jj),
     &                       diffmfv(nsp,jj),tcvx(nsp,jj),
     &                       tcvy(nsp,jj),tcvz(nsp,jj)
                        mflagl(nsp,jj) = 0
                        mflagv(nsp,jj) = 0
                     end if
                  else
                     if (msg(5) .eq. 1) then
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       mflagl(nsp,jj),diffmfl(nsp,jj),
     &                       tclx(nsp,jj),tcly(nsp,jj),
     5                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     6                       a2adfv(nsp,jj),betadfv(nsp,jj),
     7                       mflagv(nsp,jj),diffmfv(nsp,jj),
     &                       tcvx(nsp,jj),tcvy(nsp,jj)
                     else
c     Try reading using old format
                        read(inpt,*,err=2000,iostat=io_stat)
     &                       iadsfl(nsp,jj),a1adfl(nsp,jj),
     &                       a2adfl(nsp,jj),betadfl(nsp,jj),
     &                       diffmfl(nsp,jj),
     &                       tclx(nsp,jj),tcly(nsp,jj),
     5                       iadsfv(nsp,jj),a1adfv(nsp,jj),
     6                       a2adfv(nsp,jj),betadfv(nsp,jj),
     7                       diffmfv(nsp,jj),
     &                       tcvx(nsp,jj),tcvy(nsp,jj)
                        mflagl(nsp,jj) = 0
                        mflagv(nsp,jj) = 0
                     end if
                  endif                                 
               endif

C     PS        ENDIF
            end select
            if ( (diffmfl(nsp,jj) .lt. 0.) .or. 
     &           (diffmfv(nsp,jj) .lt. 0.) ) then
               write (ierr, *) 'Old Conca model not supported'
               write (ierr, *) 'See UM for new diffusion ',
     &              'model options'
               goto 2000
            end if
c     
c**** check for minimum value of diffmf
c     
            if (dispsame.eq.0) then
               if (diffmfl(nsp,jj).ge.0.) diffmfl(nsp,jj)=
     &              max(diffmfl(nsp,jj),zero_t)
               tclx(nsp,jj)=  max( tclx(nsp,jj),zero_t )
               tcly(nsp,jj)=  max( tcly(nsp,jj),zero_t )
               tclz(nsp,jj)=  max( tclz(nsp,jj),zero_t )
               
c     Square of dispersivity is used in coneq1
               
               tclx(nsp,jj) = tclx(nsp,jj)*tclx(nsp,jj)
               tcly(nsp,jj) = tcly(nsp,jj)*tcly(nsp,jj)
               tclz(nsp,jj) = tclz(nsp,jj)*tclz(nsp,jj)
               
               if (diffmfv(nsp,jj).ge.0.) diffmfv(nsp,jj)=
     &              max(diffmfv(nsp,jj),zero_t)
               tcvx(nsp,jj)=  max( tcvx(nsp,jj),zero_t )
               tcvy(nsp,jj)=  max( tcvy(nsp,jj),zero_t )
               tcvz(nsp,jj)=  max( tcvz(nsp,jj),zero_t )
               
c     Square of dispersivity is used in coneq1
               
               tcvx(nsp,jj) = tcvx(nsp,jj)*tcvx(nsp,jj)
               tcvy(nsp,jj) = tcvy(nsp,jj)*tcvy(nsp,jj)
               tcvz(nsp,jj) = tcvz(nsp,jj)*tcvz(nsp,jj)
            endif                        

C     PS
C     PS     IF this is a liquid conservative tracer
c     
            if( iadsfl(nsp,jj) .eq. 0 ) then
c     
C     PS       Set sorption parameters to reflect a conservative tracer
c     
               a1adfl(nsp,jj) = 0.
               a2adfl(nsp,jj) = 0.
               betadfl(nsp,jj) = 1.
c     
C     PS     ELSEIF this is a liquid linear isotherm
c     
            else if( iadsfl(nsp,jj) .eq. 1 ) then
c     
C     PS       Set sorption parameters to reflect a linear isotherm
c     
               a2adfl(nsp,jj) = 0.
               betadfl(nsp,jj) = 1.
c     
C     PS     ELSEIF this is a liquid Freundlich isotherm
c     
            else if( iadsfl(nsp,jj) .eq. 2 ) then
c     
C     PS       Set sorption parameters to reflect a Freundlich isotherm
c     
               a2adfl(nsp,jj) = 0.
c     
C     PS     ELSEIF this is a liquid Langmuir isotherm
c     
            else if( iadsfl(nsp,jj) .eq. 4 ) then
c     
C     PS       Set sorption parameters to reflect a Langmuir isotherm
c     
               betadfl(nsp,jj) = 1.
c     
C     PS     ENDIF
c     
            end if
C     PS     IF this is a vapor conservative tracer 
c     
            if( iadsfv(nsp,jj) .eq. 0 ) then
c     
C     PS       Set sorption parameters to reflect a conservative tracer
c     
               a1adfv(nsp,jj) = 0.
               a2adfv(nsp,jj) = 0.
               betadfv(nsp,jj) = 1.
c     
C     PS     ELSEIF this is a vapor linear isotherm
c     
            else if( iadsfv(nsp,jj) .eq. 1 ) then
c     
C     PS       Set sorption parameters to reflect a linear isotherm
c     
               a2adfv(nsp,jj) = 0.
               betadfv(nsp,jj) = 1.
c     
C     PS     ELSEIF this is a vapor Freundlich isotherm
c     
            else if( iadsfv(nsp,jj) .eq. 2 ) then
c     
C     PS       Set sorption parameters to reflect a Freundlich isotherm
c     
               a2adfv(nsp,jj) = 0.
c     
C     PS     ELSEIF this is a vapor Langmuir isotherm
c     
            else if( iadsfv(nsp,jj) .eq. 4 ) then
c     
C     PS       Set sorption parameters to reflect a Langmuir isotherm
c     
               betadfv(nsp,jj) = 1.
c     
c     
            end if
C     PS     ENDIF
c     
c     
         end do

         narrays = 1
         itype(1) = 4
         default(1) = 1
         macro = "trac"
         igroup = 13

         call initdata2( inpt, ischk, n0, narrays,
     2        itype, default, macroread(5), macro, igroup, 
     3        ireturn, i4_1=itrc(1+(nsp-1)*n0:nsp*n0) )

C     PS     
C     PS       IF the solute is a Henry's Law species
C     PS         Set the Henry's Law parameters
C     PS       ENDIF
         if( abs(icns(nsp)) .eq. 2 ) then
            read(inpt, *)henry_model(nsp)
            if(henry_model(nsp).eq.1)then
               backspace inpt
               read(inpt,*)henry_model(nsp),a_henry(nsp),
     2              dh_henry(nsp)
            elseif(henry_model(nsp).eq.2)then
               backspace inpt
               read(inpt,*)henry_model(nsp),hawwa(nsp,1),
     2              hawwa(nsp,2),hawwa(nsp,3),hawwa(nsp,4),
     3              hawwa(nsp,5)
            elseif(henry_model(nsp).eq.3)then
               backspace inpt
               read(inpt,*)henry_model(nsp),a_henry(nsp),
     2              dh_henry(nsp)
            else 
               write(ierr,*)' ** Using Old Input '
               write(ierr,*)' Enter Temperature Dependency '
               write(ierr,*)' Model Number: 1 - Van Hoff '
               write(ierr,*)' 2 - awwa model, see manual'
               write(ierr,*)' for details **'
               write(iptty,*)' Using Old Input '
               write(iptty,*)' Enter Temperature Dependency '
               write(iptty,*)' Model Number: 1 - Van Hoff '
               write(iptty,*)' 2 - awwa model, see manual'
               write(iptty,*)' for details **'
               stop
            endif
         end if
      endif

      read (inpt, '(a80)') input_msg
      if (input_msg(1:5) .eq. 'moles') then
c     zvd 22-Jul-09, 12-Apr-10
c     Modify for option to read total moles as input and compute moles/kg
c     if keyword moles is found (only entered for first species)
c     for nodal volume or zone volume (need rolf and sx1 defined)
c     (read group 15 data and write to a temporary file)
c     Also, need to be able to read zone definitions on the fly to allow 
c     multiple sources to be accumulated at the same node for overlapping zones
         mole_input = .true.
      else
         backspace inpt
      end if
      if (mole_input) then
         if (nsp .eq. 1) then
            mfile = open_file ('moles_in.tmp', 'unknown')
         end if
         write (mfile, *) nsp
         do
            read (inpt, '(a80)') input_msg
            if (null1(input_msg)) then
               write (mfile, *)
               if (nsp .eq. nspeci) close (mfile)
               exit 
            end if
            write (mfile, '(a)') trim(input_msg)
         end do
      else
c     Input is in moles/kg fluid
c     read the original way
         narrays = 1
         itype(1) = 8
         default(1) = an0
         macro = "trac"
         igroup = 15
         call initdata2( inpt, ischk, n0, narrays,
     2        itype, default, macroread(5), macro, igroup, 
     3        ireturn, r8_1=an(1+(nsp-1)*n0:nsp*n0) )

         do ij = 1, n0
            anlo(ij+npn)=an(ij+npn)
         end do

      end if
c gaz 031520 added cnsk_background
      read(inpt,'(a4)') back_trac
      if(back_trac.eq.'back') then
c new code with extra background tracer inflow
      narrays = 4
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      itype(4) = 8
      default(1) = 0.
      default(2) = 0.
      default(3) = 0.
      default(4) = 0.
      macro = "trac"
      igroup = 16
      call initdata2( inpt, ischk, n0, narrays,
     2     itype, default, macroread(5), macro, igroup, 
     3     ireturn, r8_1=cnsk(1+(nsp-1)*n0:nsp*n0),
     4     r8_2=t1sk(1+(nsp-1)*n0:nsp*n0), 
     5     r8_3=t2sk(1+(nsp-1)*n0:nsp*n0),
     6     r8_4=cnsk_background(1+(nsp-1)*n0:nsp*n0))       
          
      else
c old code
      backspace inpt
      narrays = 3
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      default(1) = 0.
      default(2) = 0.
      default(3) = 0.
      macro = "trac"
      igroup = 16
      call initdata2( inpt, ischk, n0, narrays,
     2     itype, default, macroread(5), macro, igroup, 
     3     ireturn, r8_1=cnsk(1+(nsp-1)*n0:nsp*n0),
     4     r8_2=t1sk(1+(nsp-1)*n0:nsp*n0), 
     5     r8_3=t2sk(1+(nsp-1)*n0:nsp*n0))
      endif
      
      do mim = 1,n0
         mi=mim+(nsp-1)*n0

c     Add new option for solute remaining in system despite the
c     fluid leaving BAR 4-28-99

         if(cnsk(mi).lt.0..and.t2sk(mi).lt.0.) then
            write(ierr,*)'ERROR - solute accumulation option'
            write(ierr,*)'cannot be used with cnsk<0'
            stop
         end if
         
         if(cnsk(mi).lt.0) then
            pcnsk(mi)=-1.0
         elseif(t2sk(mi).lt.0.) then
            t2sk(mi) = abs(t2sk(mi))
            pcnsk(mi) = 1.
         end if
         
      enddo

C     PS ENDLOOP through each species
      end do

      if (cden_flag .ne. 0) then
         read (inpt, '(a4)', end = 9) macro
 9       if (macro .eq. 'cden') then
            call incden(inpt)
         else
            write (ierr, 10)
            if (iout .ne. 0) write (iout, 10)
            if (iptty .ne. 0) write (iptty, 10)
            cden = .false.
            backspace (inpt)
         end if
      end if
 10   format ('Molecular weight input not found for cden macro ,',
     &     ' cden will not be used')

      macroread(5) = .TRUE.

C     PS ENDIF the data are to be read in on this pass
      elseif(trxn_flag .eq. 1) then
         call rdtr
                  do npn=1,nspeci
                  do ij = 1, n0
                     anlo((npn-1)*n0+ij)=an((npn-1)*n0+ij)
                  end do
                  end do

	if (molein) mole_input = .true.
      else
         write(ierr, *) 'Internal error differentiating trac/trxn.'
         stop
      endif
	!call trxn_varcheck
      endif
C     PS IF this call to the routine is for initialization 
C     PS
      if( iz .lt. 0 ) then
         allocate (hflag(n7))
         hflag = 0
         if (mole_input) then
c     zvd 22-Jul-09
c     Modify for option to read total moles as input and compute moles/kg
c     for nodal volume or zone volume
            allocate (seen(n7))
            alt_zone = .false.
            mfile = open_file ('moles_in.tmp', 'old')
            if (.not. allocated(lvol)) 
     &           allocate(lvol(n0), rvol(n0), vvol(n0))
c     Compute rock, liquid, and vapor volume for each node
            lvol = 0.
            rvol = 0.
            vvol = 0.
            do i = 1, n0
               rvol(i) = sx1(i)*denr(i)*(1-ps_trac(i))
               if (irdof .ne. 13 .or. ifree.ne.0) then
                  sdum = s(i)
               else
                  sdum = 1.0
               end if
               lvol(i) = sx1(i)*sdum*rolf(i)*ps_trac(i)
               vvol(i) = sx1(i)*(1 - sdum)*rovf(i)*ps_trac(i) 
            end do

            do nsp = 1, nspeci
               if (.not. conc_read(nsp)) then
                  read (mfile, *) idum
                  if (idum .ne. nsp) then
c     We didn't save input for this species, this is an error
                     write (ierr, *) 
     &                    'Error:  Missing mole input data for species ', nsp, '.'
                     stop
                  end if
                  npn=npt(nsp)
                  do
                     read (mfile, '(a80)') input_msg
                     if (null1(input_msg)) exit
                     if (input_msg(1:) .eq. 'file') then
c     Read name of alternate zone file
                        read (mfile, '(a100)') zone_file
			if (izf .ne. 0) then
				close(izf)
				izf = 0
			endif
                        izf = open_file(zone_file, 'unknown')
                        read (mfile, '(a80)') input_msg
                        alt_zone = .true.
                     else if (input_msg(1:6) .eq. 'altend') then
c     Use currently defined zone list (default)
                        alt_zone = .false.
                     end if
                     read (input_msg, *) ja, jb, jc, andum
                     if ( ja .eq. 1 .and. jb .eq. 0 .and. jc .eq. 0 ) 
     &                    then
                        jb = n0
                        jc = 1
                     end if

c     Input is total moles
                     if (ja .lt. 0) then
c     Using zone definitions
                        mvol = 0.d0
                        mvolv = 0.d0
                        if (alt_zone) then
c     Assumes nnum style zone input
                           num_nodes = 0
                           rewind (izf)
                           read (izf, '(a4)', end = 77) zmacro
                           if (zmacro(1:3) .ne. 'zon') then
                              write (ierr, *) 
     &   'Error:  Zone file for mole input is not a zone file.'
				stop
                           end if
                           do
                              read (izf, *, end = 77) znum
                              read (izf, '(a4)') zmacro
                              if (zmacro .ne. 'nnum') then
                                 write (ierr, *) 
     &   'Error:  Zone file for mole input can only use method "nnum".'
                                 stop
                              end if
                              read (izf, *) num_nodes
                              if (num_nodes .gt. 0) then
                                 if (znum .eq. abs(ja)) then
                                    allocate(node_list(num_nodes))
                                    read (izf, *) (node_list(i), i = 1, 
     &                                   num_nodes)
                                    exit
                                 else
                                    read (izf, *) (rdum1, i = 1, 
     &                                   num_nodes)
                                    num_nodes = 0
                                 end if
                              end if
                           end do
 77                        nend = num_nodes
                        else
                           nend = n0
                        end if
c     Find the total volume for current zone to apportion moles
                        do nin = 1, nend
                           if (alt_zone) then
                              i = node_list(nin)
                              found = .true.
                           else
                              i = nin
                              if (izonef(nin) .eq. abs (ja)) then
                                 found = .true.
                              else
                                 found = .false.
                              end if
                           end if
                           if (found .and. ps_trac(i) .gt. 0.) then
c     Only include if porosity > 0
                              if (icns(nsp) .eq. 0) then
                                 mvol = mvol + rvol(i)
                              else if (icns(nsp) .eq. 1) then
                                 mvol = mvol + lvol(i)
                              else if (icns(nsp) .eq. -1) then
                                 mvol = mvol + vvol(i)
                              else
c     Henry's law species
                                 mvol = mvol + lvol(i)
                                 if(henry_model(nsp).eq.1) then
                                    h_const= a_henry(nsp)*
     2                                   exp(dh_henry(nsp)/
     3                                   (gas_const*0.001)*(1/298.16-1/
     4                                   (t(i)+temp_conv)))
                                 else if(henry_model(nsp).eq.2) then
                                    h_const = 10**(hawwa(nsp,1)+
     2                                   hawwa(nsp,2)*(t(i)+
     3                                   temp_conv)+hawwa(nsp,3)
     4                                   /(t(mim)+temp_conv)+hawwa
     5                                   (nsp,4)*dlog10(t(i)
     6                                   +temp_conv)+hawwa(nsp,5)
     7                                   /(t(i)+temp_conv)**2)
                                    h_const= (101325*rolf(i)*1e-3)/
     2                                   (h_const*1e6*mw_water)
                                 else if(henry_model(nsp).eq.3) then
                                    h_const= (phi(i) - pci(i)) * 
     2                                   dh_henry(nsp)
                                 endif
                                 dvap_conc = (mw_water*h_const) /
     2                                (phi(i)*avgmolwt(i))
                                 mvolv = mvolv + dvap_conc * vvol(i)
                              end if
                           end if
                        end do
                        found = .false.
c     Apportion mass to each node in zone
                        do nin = 1, nend
                           if (alt_zone) then
                              i = node_list(nin)
                              found = .true.
                           else
                              i = nin
                              if (izonef(nin) .eq. abs (ja)) then
                                 found = .true.
                              else
                                 found = .false.
                              end if
                           end if
                           if (found .and. ps_trac(i) .gt. 0.) then
c     Only include if porosity > 0
                              mi = i + npn
                              hflag(mi) = 1
                              dvap_conc = 0.
                              if ( abs(icns(nsp)) .eq. 2 ) then
c     Henry's law species
                                 if(henry_model(nsp).eq.1) then
                                    h_const= a_henry(nsp)*
     2                                   exp(dh_henry(nsp)/
     3                                   (gas_const*0.001)*(1/298.16-1/
     4                                   (t(i)+temp_conv)))
                                 else if(henry_model(nsp).eq.2) then
                                    h_const = 10**(hawwa(nsp,1)+
     2                                   hawwa(nsp,2)*(t(i)+
     3                                   temp_conv)+hawwa(nsp,3)
     4                                   /(t(mim)+temp_conv)+hawwa
     5                                   (nsp,4)*dlog10(t(i)
     6                                   +temp_conv)+hawwa(nsp,5)
     7                                   /(t(i)+temp_conv)**2)
                                    h_const= (101325*rolf(i)*1e-3)/
     2                                   (h_const*1e6*mw_water)
                                 else if(henry_model(nsp).eq.3) then
                                    h_const= (phi(i) - pci(i)) * 
     2                                   dh_henry(nsp)
                                 endif
                                 dvap_conc = (mw_water*h_const) /
     2                                (phi(i)*avgmolwt(i))
                              end if
                              if (.not. seen(mi)) then
                                 an(mi) = abs(andum) / (mvol + mvolv)
                                 seen(mi) = .true.
                              else
                                 an(mi) = an(mi) + abs(andum) / 
     &                                (mvol + mvolv)
                              end if
                           end if
                        end do
                        if (allocated(node_list)) deallocate (node_list)
                     else   
                        do i = ja, jb, jc
                           mi = i + npn
                           hflag(mi) = 1
                           if (ps_trac(i) .ne. 0) then
                              dvap_conc = 0.
                              if ( abs(icns(nsp)) .eq. 2 ) then
c     Henry's law species
                                 if(henry_model(nsp).eq.1) then
                                    h_const= a_henry(nsp)*
     2                                   exp(dh_henry(nsp)/
     3                                   (gas_const*0.001)*(1/298.16-1/
     4                                   (t(i)+temp_conv)))
                                 else if(henry_model(nsp).eq.2) then
                                    h_const = 10**(hawwa(nsp,1)+
     2                                   hawwa(nsp,2)*(t(i)+
     3                                   temp_conv)+hawwa(nsp,3)
     4                                   /(t(i)+temp_conv)+hawwa
     5                                   (nsp,4)*dlog10(t(i)
     6                                   +temp_conv)+hawwa(nsp,5)
     7                                   /(t(i)+temp_conv)**2)
                                    h_const= (101325*rolf(i)*1e-3)/
     2                                   (h_const*1e6*mw_water)
                                 else if(henry_model(nsp).eq.3) then
                                    h_const= (phi(i) - pci(i)) * 
     2                                   dh_henry(nsp)
                                 endif
                                 dvap_conc = (mw_water*h_const) /
     2                                (phi(i)*avgmolwt(i))
                              end if
                              if (.not. seen(mi)) then
                                 seen(mi) = .true.
                                 an(mi) = abs(andum) / (lvol(i) + 
     &                                vvol(i) * dvap_conc)
                              else
                                 an(mi) = an(mi) + abs(andum) / 
     &                                (lvol(i) + vvol(i) * dvap_conc)
                              end if
                           end if
                        end do
                     end if
                  end do

                  do i = 1, n0
c     Account for sorption (moles input is total)
                     mi = i + npn
                     if (seen(mi)) then
                        rdum1 = an(mi) * (lvol(i) +
     &                       vvol(i) * dvap_conc)
                        sctmp = 0.
                        if (iadsfl(nsp,itrc(mi)) .ne. 0 .or.
     &                       iadsfv(nsp,itrc(mi)) .ne. 0 ) then
                           ctmp = rdum1
                           antmp = an(mi)
                           do
                              if (icns(nsp) .gt. 0) then
                                 sctmp = (denr(i) * a1adfl(nsp,
     1                                itrc(mi)) * antmp**betadfl
     2                                (nsp,itrc(mi)) / (1.0 + 
     3                                a2adfl(nsp,itrc(mi)) * antmp
     4                                **betadfl(nsp,itrc(mi)) ))
                              else if (icns(nsp) .lt. 0) then
                                 if (icns(nsp) .eq. -2) then
                                    anvtmp = antmp * dvap_conc
                                 else
                                    anvtmp = antmp
                                 end if
                                 sctmp = (denr(i)*a1adfv(nsp,
     1                                itrc(mi)) * anvtmp**betadfv
     2                                (nsp,itrc(mi)) / (1.0 + 
     3                                a2adfv(nsp,itrc(mi)) *anvtmp
     4                                **betadfv(nsp,itrc(mi)) ))
                              end if
                              sctmp = sctmp * sx1(i)
                              rdum2 = ctmp + sctmp
                              if (abs(rdum1 - rdum2) .le. ctol) 
     &                             then
                                 an(mi) = antmp
                                 exit
                              end if
                              if (rdum2 .lt. rdum1) then
                                 ctmp = ctmp + (rdum1 - rdum2)
                                 antmp = ctmp / (lvol(i) + 
     &                                vvol(i) * dvap_conc)
                              else
                                 frac = rdum1/ rdum2
                                 antmp = antmp * frac
                                 ctmp = antmp * (lvol(i) 
     &                                + vvol(i) * dvap_conc)
                              end if 
                           end do
                        end if
                     end if
                  end do
                  do ij = 1, n0
                     anlo(ij+npn)=an(ij+npn)
                  end do
               end if
            end do
            close (mfile, status = 'delete')
		inquire(izf, opened=opend)
		if (opend) close(izf)
         end if
c     
c     seh
C     C Check for keyword specifying to use same dispersivity and diffusivity
C     C for all liquid / vapor calculations
         if (dispsame.eq.1) then
            if (hvliquid.eq.1) then
               do jj=1,numd
                  tclx(1,jj)=max(tclx(1,jj),zero_t)
                  tcly(1,jj)=max(tcly(1,jj),zero_t)
                  tclz(1,jj)=max(tclz(1,jj),zero_t)
c     Square of dispersivity is used in alpha calculation
                  tclx(1,jj)=tclx(1,jj)*tclx(1,jj)
                  tcly(1,jj)=tcly(1,jj)*tcly(1,jj)
                  tclz(1,jj)=tclz(1,jj)*tclz(1,jj)
               enddo
               neqp1=neq+1
               if (ldsp.eq.0) then
                  sehindexl=1
                  if (icnl.eq.0) then
                     do iz2=1,n0
                        if (iz2.gt.neq) then
                           iz2p=iz2-neq
                        else
                           iz2p=iz2
                        endif
                        i1=nelm(iz2p)+1
                        i2=nelm(iz2p+1)
                        ipivt=nelmdg(iz2p)
                        do jj=ipivt+1,i2
                           iw=istrw(jj-neqp1)
                           if (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
     &                          .ne.0.) then
                              kz=nelm(jj)
                              cord1x=cord(iz2p,1)
                              cord1y=cord(iz2p,2)
                              cord1z=cord(iz2p,3)
                              cord2x=cord(kz,1)
                              cord2y=cord(kz,2)
                              cord2z=cord(kz,3)
                              ii=itrcdsp(iz2)
                              dispxavw=tclx(1,ii)
                              dispyavw=tcly(1,ii)
                              dispzavw=tclz(1,ii)
                              tempx = (cord1x-cord2x)**2
                              tempy = (cord1y-cord2y)**2
                              tempz = (cord1z-cord2z)**2
                              templength = tempx+tempy+tempz
                              alphaconl(sehindexl) =
     +                             sqrt((dispxavw*dispyavw*dispzavw*
     +                             templength)/(dispyavw*dispzavw*tempx+
     3                             dispxavw*dispzavw*tempy +
     4                             dispxavw*dispyavw*tempz))
                              sehindexl=sehindexl+1
                           endif
                        enddo
                     enddo   
                  else 
                     do iz2=1,n0
                        if (iz2.gt.neq) then
                           iz2p=iz2-neq
                        else
                           iz2p=iz2
                        endif
                        i1=nelm(iz2p)+1
                        i2=nelm(iz2p+1)
                        ipivt=nelmdg(iz2p)
                        do jj=ipivt+1,i2
                           iw=istrw(jj-neqp1)
                           if (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
     &                          .ne.0.) then
                              kz=nelm(jj)
                              cord1x=cord(iz2p,1)
                              cord1y=cord(iz2p,2)
                              cord2x=cord(kz,1)
                              cord2y=cord(kz,2)
                              tempx = (cord1x-cord2x)**2
                              tempy = (cord1y-cord2y)**2
                              ii=itrcdsp(iz2)
                              alphaconl(sehindexl)=
     +                             sqrt((tclx(1,ii)*tcly(1,ii)*
     &                             (tempx+tempy))/(tclx(1,ii)*tempy
     &                             + tcly(1,ii)*tempx ))
                              sehindexl=sehindexl+1
                           endif
                        enddo
                     enddo
                  endif
               endif
            endif
            if (hvvapor.eq.1) then
               do jj=1,numd
                  tcvx(1,jj)=max(tcvx(1,jj),zero_t)
                  tcvy(1,jj)=max(tcvy(1,jj),zero_t)
                  tcvz(1,jj)=max(tcvz(1,jj),zero_t)
c     square of dispersivity used in alpha calculation
                  tcvx(1,jj)=tcvx(1,jj)*tcvx(1,jj)
                  tcvy(1,jj)=tcvy(1,jj)*tcvy(1,jj)
                  tcvz(1,jj)=tcvz(1,jj)*tcvz(1,jj)
               enddo
               neqp1=neq+1
               if (ldsp.eq.0) then
                  sehindexv=1
                  if (icnl.eq.0) then
                     do iz2=1,n0
                        if (iz2.gt.neq) then
                           iz2p=iz2-neq
                        else
                           iz2p=iz2
                        endif
                        i1=nelm(iz2p)+1
                        i2=nelm(iz2p+1)
                        ipivt=nelmdg(iz2p)
                        do jj=ipivt+1,i2
                           iw=istrw(jj-neqp1)
                           if (sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
     &                          .ne. 0.) then
                              kz=nelm(jj)
                              cord1x=cord(iz2p,1)
                              cord1y=cord(iz2p,2)
                              cord1z=cord(iz2p,3)
                              cord2x=cord(kz,1)
                              cord2y=cord(kz,2)
                              cord2z=cord(kz,3)
                              ii=itrcdsp(iz2)
                              dispxavw=tcvx(1,ii)
                              dispyavw=tcvy(1,ii)
                              dispzavw=tcvz(1,ii)
                              tempx = (cord1x-cord2x)**2
                              tempy = (cord1y-cord2y)**2
                              tempz = (cord1z-cord2z)**2
                              templength = tempx+tempy+tempz
                              alphaconv(sehindexv) =
     +                             sqrt((dispxavw*dispyavw*dispzavw*
     +                             templength)/(dispyavw*dispzavw*tempx+
     3                             dispxavw*dispzavw*tempy +
     4                             dispxavw*dispyavw*tempz))
                              sehindexv=sehindexv+1
                           endif
                        enddo
                     enddo
                  else
                     do iz2=1,n0
                        if (iz2.gt.neq) then
                           iz2p=iz2-neq
                        else
                           iz2p=iz2
                        endif
                        i1=nelm(iz2p)+1
                        i2=nelm(iz2p+1)
                        ipivt=nelmdg(iz2p)
                        do jj=ipivt+1,i2
                           iw=istrw(jj-neqp1)
                           if (sx(iw,1)+sx(iw,2)+sx(iw,3).ne.0.) then
                              kz=nelm(jj)
                              cord1x=cord(iz2p,1)
                              cord1y=cord(iz2p,2)
                              cord2x=cord(kz,1)
                              cord2y=cord(kz,2)
                              tempx = (cord1x-cord2x)**2
                              tempy = (cord1y-cord2y)**2
                              ii=itrcdsp(iz2)
                              alphaconv(sehindexv)=
     +                             sqrt((tcvx(1,ii)*tcvy(1,ii)*
     &                             (tempx+tempy))/(tcvx(1,ii)*tempy
     &                             + tcvy(1,ii)*tempx ))
                              sehindexv=sehindexv+1
                           endif
                        enddo
                     enddo
                  endif
               endif
            endif
         endif
         
         ntpp=n7/nspeci 
         do ispecies = 1,nspeci
            if( icns(ispecies) .eq. 0 ) then
               npt_subst =(ispecies-1)*ntpp
               do ij = 1,n
                  itrc(ij+npt_subst)=1
               enddo
            end if
         enddo

         if(ntpp.lt.n) then
            write (ierr, 6000)
            if (iout .ne. 0) write(iout  ,6000)
 1          if ( iptty .gt. 0 )  write(iptty ,6000)
            iret = -1
            goto 9000
         endif
 6000    format(/,1x,'**** memory too small for multiple tracers ****'
     *        ,/,1x,'****---------------------------------------****',
     *        /)

         do nsp=1,nspeci
            npn=npt(nsp)
            cm0(nsp)=0.
            dtotc=1.0
C     PS     
C     PS     IF this solute is a Henry's Law species with vapor...
C     PS     ... concentration specified
c     
            if( icns(nsp) .eq. -2  ) then
c     
C     PS     
C     PS       IF this is a dpdp simulation
c     
               if(idpdp.ne.0) then
c     
C     PS         Set number of node sets to 2
c     
                  n_node_sets = 2
               else if( idualp .ne. 0 ) then
                  n_node_sets = 3
c     
C     PS       ELSEIF it is dual porosity
c     
               else
c     
C     PS         Set number of node sets to 1
c     
                  n_node_sets = 1
c     
C     PS   ENDIF
c     
               end if
C     PS   
C     PS       FOR each node set
c     
               do i_node_set = 1, n_node_sets
                  ndummy = ( i_node_set - 1 ) * neq
c     
C     PS         FOR each node
c     
                  do i = 1, neq
c     
                     mi = i + ndummy + npt(nsp)
                     mim = i + ndummy
c     
C     PS           Compute liquid concentration corresponding to the input...
C     PS           ... vapor concentration
c     
                     if(henry_model(nsp).eq.1)then
                        h_const= a_henry(nsp)*
     2                       exp(dh_henry(nsp)/
     3                       (gas_const*0.001)*(1/298.16-1/
     4                       (t(mim)+temp_conv)))
                     elseif(henry_model(nsp).eq.2)then
                        h_const = 10**(hawwa(nsp,1)+
     2                       hawwa(nsp,2)*(t(mim)+temp_conv)+
     3                       hawwa(nsp,3)/(t(mim)+temp_conv)+
     4                       hawwa(nsp,4)*dlog10(t(mim)+temp_conv)+
     5                       hawwa(nsp,5)/(t(mim)+temp_conv)**2)
                        h_const= (101325*rolf(mim)*1e-3)/(h_const*1e6*
     2                       mw_water)
                     elseif(henry_model(nsp).eq.3)then
                        h_const= (phi(i) - pci(i)) * dh_henry(nsp) 
                     else
                        write(ierr,*)' ** Using Old Input '
                        write(ierr,*)' Enter Temperature Dependency '
                        write(ierr,*)' Model Number: 1 - Van Hoff '
                        write(ierr,*)' 2 - awwa model, see manual'
                        write(ierr,*)' for details **'
                        write(iptty,*)' Using Old Input '
                        write(iptty,*)' Enter Temperature Dependency '
                        write(iptty,*)' Model Number: 1 - Van Hoff '
                        write(iptty,*)' 2 - awwa model, see manual'
                        write(iptty,*)' for details **'
                        stop
                     endif
                     dvap_conc = ( mw_water * h_const ) /
     2                    ( phi(mim) * avgmolwt(mim) )
C     PS           IF there is no transport to the liquid phase
C     PS             Write error message to the screen
C     PS             Set flag to stop program
C     PS             ERROREXIT
C     PS           ENDIF
                     if( h_const .lt. 1.e-20 ) then
                        write(ierr, 999) nsp, mim
                        if (iptty .ne. 0) write(iptty, 999) nsp, mim
                        iret = -1
                        goto 9000
                     end if
 999                 format('Fatal error',
     .                    'You specified a Henrys Law species',
     .                    'with initial concentrations input for',
     .                    'the vapor phase (icns = -2), yet the',
     .                    'Henrys Constant is computed as 0',
     .                    'for species number ', i8,
     .                    'and',
     .                    'node number ', i8, '.',
     .                    'If you want to simulate a vapor-borne',
     .                    'species with no interphase transport,',
     .                    'then you must specify a gaseous',
     .                    'species (icns = -1).')
                     if (hflag(mi) .eq. 0) then                 
                        an(mi) = an(mi) / dvap_conc
                        anlo(mi) = an(mi)
                     end if
c     
C     PS         ENDFOR each node
c     
                  end do
C     PS       ENDFOR each node set
c     
               end do
C     PS   
C     PS     ENDIF
C     PS     
            end if
c     
            call thermc(0)
            if(idpdp.ne.0) then
               call thermc(neq)
            else if(idualp .ne. 0 ) then
               call thermc(neq)
               call thermc(neq+neq)
            endif

            do ij=1,n
               dench(ij+npn)=0.
               dench(ij+npn)=denci(ij+npn)*dtotc
               denci(ij+npn)=0.
               dencj(ij+npn)=0.0
               if (ps_trac(ij) .gt. 0) then
                  cm0(nsp)=cm0(nsp)+dench(ij+npn)*volume(ij)
               end if
            end do
         end do
c     zvd 05/16/2007
c     Move call to startup to ensure all data read prior to setup 
c     to accomodate new output options
c     call plotc1(0,0)
      end if

      if (allocated(lvol)) deallocate(lvol,rvol,vvol,hflag,seen)
      return

 9000 continue
      if( iret .ne. 0 ) then
         stop
      end if

 2000 if (io_stat .ne. 0) then
         if (dispsame .eq. 1) then
            write (ierr, *) 'Error reading group 9 in trac'
         else
            write (ierr, *) 'Error reading group 12 in trac'
         end if
         stop
      end if
      end

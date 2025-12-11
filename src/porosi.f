      subroutine  porosi(iz)
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
CD1 To calculate the pressure dependant porosity and permeability, and 
CD1 derivatives.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/porosi.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:40   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:58   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Mon Mar 31 08:40:52 1997   gaz
CD2 major changes for specific discharge and aquifer compressibility
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:14 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:18:16 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:58:22 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Fri Feb 16 10:56:58 1996   zvd
CD2 Updated purpose and added requirement.
CD2 
CD2    Rev 1.2   Thu Feb 01 14:18:00 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:11:32   gaz
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
CD3 Identifier              Type     Use  Description
CD3
CD3 iz                      int       I   Control parameter
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 Identifier    Use  Description
CD3
CD3 Number iout    O   Output from simulation
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 iporos, n0, amgang, pgangi, inpt, thexp, young, sigini, wdd1, psini,
CD4 phi, phini, dporp, dport, iad, iporf, ps, m, nskw, pho, pein, 
CD4 narrays, pointer, itype, default, igroup, ireturn,
CD4 macroread
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4 
CD4 Identifier   Type   Description
CD4 
CD4 welbor       N/A    Compute wellbore solution
CD4 initdata     N/A    Reads and set parameter values by zone or node
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
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
CD5 jji          int         Do loop index parameter
CD5 pso          real*8      Initial porosity at current node
CD5 delp         real*8      Change in pressure from previous time
CD5 
CD5 rlf          real*8      Relative permeability of liquid
CD5 drlpf        real*8      Derivative of relative permeability
CD5 drlef        real*8      Derivative of relative permeability
CD5 
CD5 closur       real*8      Porosity closure parameter
CD5 gt           real*8      Parameter used in porosity calculation
CD5 xrl          real*8      Relative permeability
CD5 dgtp         real*8      Parameter used in porosity calculation
CD5 dgtt         real*8      Parameter used in porosity calculation
CD5 drl          real*8      Parameter used in porosity calculation
CD5 dclosp       real*8      Derivative of closur with pressure
CD5 dclost       real*8      Derivative of closur with temperature
CD5 dpp          real*8      Derivative of porosity with pressure
CD5 md           int         Current node number for output
CD5 delpsb       real*8      Delta pressure output
CD5 permsb       real*8      Permeability
CD5 dsigt        real*8      Stress
CD5 tmpPor       real*8      porosity
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
CD8 Only liquid problems are implemented for pressure-dependent
CD8 porosity.
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.4.10 Stress-dependent properties
CD9 2.6    Provide Input/Output Data Files
CD9 3.0    INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN porosi
CPS 
CPS IF pressure dependent porosities are used
CPS 
CPS   IF this call is to read input
CPS   
CPS     FOR each node
CPS       Initialize porosity law parameters
CPS     ENDFOR
CPS     
CPS     Read in first line of input
CPS     
CPS     Set parameters for reading porosity parameters
CPS     
CPS     initdata - read porosity parameters and set at each node
CPS     
CPS     Set flag to denote that the ppor macro has been called
CPS       
CPS   ELSEIF this call is to compute porosities
CPS   
CPS     IF a linear compressibility model is used
CPS       FOR all node
CPS         IF the porosity is nonzero
CPS           Compute porosity and derivatives
CPS         ENDIF
CPS       ENDFOR
CPS     ELSEIF the nonlinear Gangi opening law is used
CPS     
CPS       FOR each node
CPS         IF this is the first iteration
CPS           Set flag values to indicate all nodes are open
CPS         ENDIF
CPS         IF porosity is unity or greater
CPS           Set porosity to unity and derivatives to 0
CPS         ELSEIF porosity parameter is greater than 0
CPS           Compute closure parameters
CPS           IF the closure is less than 0 or the node is already...
CPS           ... fully closed
CPS             Set parameters to zero out porosity and derivatives
CPS           ELSE
CPS             Compute parameters needed for porosity calculation
CPS           ENDIF
CPS           Compute porosity and derivatives
CPS         ENDIF
CPS       ENDFOR each node
CPS       
CPS     ELSEIF the well bore routines are used to compute porosity
CPS       welbor - call routine to compute porosity in wellbore module
CPS     ENDIF
CPS   
CPS   ELSEIF this is a call to write porosity information
CPS     IF this is not a wellbore simulation
CPS       FOR each output node
CPS         Compute parameters
CPS         Write to output file
CPS       ENDFOR
CPS     ELSE
CPS       welbor - write porosity information
CPS     ENDIF
CPS 
CPS   ELSEIF this is a call to initialize permeabilities and porosities
CPS   
CPS     IF porosities need to be initialized
CPS     
CPS       IF this is the first iteration
CPS         Set flag values to indicate all nodes are open
CPS       ENDIF
CPS       IF porosity is unity or greater
CPS         Set porosity to unity and derivatives to 0
CPS       ELSEIF porosity parameter is greater than 0
CPS         Compute closure parameters
CPS         IF the closure is less than 0 or the node is already...
CPS         ... fully closed
CPS           Set parameters to zero out porosity and derivatives
CPS         ELSE
CPS           Compute parameters needed for porosity calculation
CPS         ENDIF
CPS         Compute porosity and derivatives
CPS       ENDIF
CPS       
CPS     ENDIF
CPS   
CPS   ENDIF this is a call to initialize permeabilities and porosities
CPS 
CPS ENDIF pressure dependent porosities are used
CPS 
CPS END porosi
CPS 
C**********************************************************************
c
c calculate pressure dependant porosity and derivative

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comki
      use comxi
      use comwt
      use comnire
      use comchem
c
      implicit  none
c
      integer  iz, i, jji, md, i_file, i1, i2, ii, kb
      integer  yama1, yama2, open_file
c
      real*8  dpp ,dclost, pso, delp, closur, gt, xrl, dgtp, dgtt, drl
      real*8  dclosp, delpsb, permsb, dsigt, sy, alength, tmpPor
      real*8  den_w0, comp_w0,hmid,hmin,hmax,por_ll, pwv
      real*8  pnx_new, dum1
      real*8  days_last
      parameter(den_w0=997.808d00,comp_w0=5.687D-4,por_ll=1.d-8)
      save yama1, yama2, days_last

c
      if ( iporos .ne. 0 )  then
c check for steady-state
         if(isteady.ne.0.and.iz.ne.0) return
c
         if ( iz .eq. 0 )  then
c
****   read input data   ****
c
****   added Yamashita_Tenma-Gangi model   ****
c
            read (inpt,*)  iporos
c
            if ( iporos .eq. -4 )  then
c tenma_ww array now allocated in allocmem
c
               backspace  inpt
c
               read (inpt,      *)  iporos
c
               do   i_file=100,1,-1
                  if ( nmfil(5)(i_file:i_file) .eq. '.' )  then
                     yama_filename1( 1:i_file  ) =  nmfil(5)(1:i_file)
                     yama_filename2( 1:i_file  ) =  nmfil(5)(1:i_file)
                     yama_filename1(i_file+1:i_file+6) =  'gangi1'
                     yama_filename2(i_file+1:i_file+6) =  'gangi2'
                     go  to  10
                  end   if
               end   do
 10            continue
c
               yama1 = open_file(yama_filename1,'unknown')
               yama2 = open_file(yama_filename2,'unknown')

c            open (41, file=yama_filename1,
c     &            status='unknown', form='formatted' )
c            open (42, file=yama_filename2,
c     &            status='unknown', form='formatted' )
c
               if(m3.gt.0) then
                  write(yama2,3201)  0, 0.0, ( nskw3(i), i=1,m3 )
               endif
 3201          format(i10,g15.7,250i7)
c
            end   if
c
****   analyse input data   ****
c
            if ( ico2 .lt. 0 )  then
c
               if ( abs(iporos) .eq. 2 )  then
                  if (iout .ne. 0) then
                     write(iout,*)
                     write(iout, 100)
                     write(iout, 110)
                     write(iout,*)
                  end if
                  if ( iptty .ne. 0 )  then
                     write(iptty,*)
                     write(iptty, 100)
                     write(iptty, 110)
                     write(iptty,*)
                  end  if
               end  if
 100           format ('>>> Warning thermal term inoperative')
 110           format ('for gangi porosity model in isothermal mode')
c
            else if ( ico2 .eq. 0 )  then
c
               if ( iporos .eq. -1 )  then
                  if (iout .ne. 0) then
                     write(iout,*)
                     write(iout, 120)
                     write(iout, 130)
                     write(iout,*)
                  end if
                  if ( iptty .ne. 0 )  then
                     write(iptty,*)
                     write(iptty, 120) 
                     write(iptty, 130)
                     write(iptty,*)
                  end  if
               end  if
 120           format ('>>> Warning: Specific storage uses ')
 130           format ('properties at 20 C nonisothermal problem')
c
            else if ( ico2 .gt. 0 )  then

               if ( abs(iporos) .eq. 2 )  then
                  write(ierr,*)
                  write(ierr, 140)
                  write(ierr, 150)
                  write(ierr,*)
                  if (iout .ne. 0) then
                     write(iout,*)
                     write(iout, 140)
                     write(iout, 150)
                     write(iout,*)
                  end if
                  if ( iptty .ne. 0 )  then
                     write(iptty,*)
                     write(iptty, 140)
                     write(iptty, 150)
                     write(iptty,*)
                  end   if
                  stop
               end   if
 140           format ('>>> Gangi model not yet available')
 150           format (' for air-water-heat conditions : stopping')
c
               if ( iporos .eq. -1 )  then
                  write(ierr,*)
                  write(ierr, 160)
                  write(ierr, 170)
                  write(ierr,*)
                  if (iout .ne. 0) then
                     write(iout,*)
                     write(iout, 160)
                     write(iout, 170)
                     write(iout,*)
                  end if
                  if ( iptty .ne. 0 )  then
                     write(iptty,*)
                     write(iptty, 160)
                     write(iptty, 170)
                     write(iptty,*)
                  end  if
                  stop
               end   if
 160           format ('>>> Specific storage not available ')
 170           format ('for non isothermal conditions : stopping')
c
            end   if
c
****   GZ Gangi model   ****
c
            if ( abs(iporos) .eq. 2 )  then
c
               narrays    =  4
               itype  (1) =  8
               itype  (2) =  8
               itype  (3) =  8
               itype  (4) =  8
               default(1) =  1.0
               default(2) =  1.0e+30
               default(3) =  0.0   
               default(4) =  0.0   
               macro      =  "ppor"
               igroup     =  2
c
               call  initdata2  ( inpt, ischk, n0, narrays, itype, 
     *              default, macroread(4), macro, igroup, ireturn,
     1              r8_1=amgang(1:n0),
     2              r8_2=pgangi(1:n0),
     3              r8_3=strgan(1:n0),
     4              r8_4=alphae(1:n0) )
c
               macroread(4) =  .TRUE.
c
****   aquifer compressibility model & specific storege model   ****
c
            else if ( abs(iporos) .eq. 1 )  then
c
               narrays    =   1
               itype  (1) =   8
               default(1) =   0.0
               macro      =   "ppor"
               igroup     =   2
c
               call  initdata2  ( inpt, ischk, n0, narrays, itype, 
     *              default, macroread(4), macro, igroup, ireturn,
     1              r8_1=amgang(1:n0) )
c
               do i=1,n0
                  if(amgang(i).lt.0.0) then
                     amgang(i) = 1.0d00/10.0d00**(abs(amgang(i)))
                  endif
               enddo
c
               macroread(4) =  .TRUE.

c MODFLOW-like specific yield
            else if ( iporos .eq. -5 )  then
               narrays    =   1
               itype  (1) =   8
               default(1) =   0.0d0
               macro      =   "ppor"
               igroup     =   2
c
               call  initdata2  ( inpt, ischk, n0, narrays, itype, 
     *              default, macroread(4), macro, igroup, ireturn,
     *              r8_1=amgang(1:n0) )
c
               macroread(4) =  .TRUE.
c
****   Yamashita_Tenma-Gangi model   ****
c
            else if ( iporos .eq. -4 )  then
c
               narrays    =  5
               itype  (1) =  8
               itype  (2) =  8
               itype  (3) =  8
               itype  (4) =  8
               itype  (5) =  8
               default(1) =  0.5
               default(2) =  1.0e+20
               default(3) =  0.0
               default(4) =  16.6
               default(5) =  0.03
               macro      =  "ppor"
               igroup     =  2
c
               call  initdata2  ( inpt, ischk, n0, narrays, itype, 
     1              default, macroread(4), macro, igroup, ireturn,
     2              r8_1=amgang(1:n0),
     3              r8_2=pgangi(1:n0),
     4              r8_3=wgangi(1:n0),
     5              r8_4=sgangi(1:n0),
     6              r8_5=agangi(1:n0) )
c
               macroread(4) =  .TRUE.
c
            else if ( iporos .eq. 5 )  then
c               this option was added by ZL for temperature-dependent porosity

               narrays    =  4
               itype  (1) =  8
               itype  (2) =  8
               itype  (3) =  8
               itype  (4) =  8
               default(1) =  0.2
               default(2) =  10.0
               default(3) =  0.0
               default(4) =  0.0
               macro      =  "ppor"
               igroup     =  2

               call  initdata2  ( inpt, ischk, n0, narrays, itype,
     *              default, macroread(4), macro, igroup, ireturn,
     1              r8_1=porTemp1(1:n0),
     2              r8_2=porTemp2(1:n0),
     3              r8_3=porTemp3(1:n0),
     4              r8_4=porTemp4(1:n0) )

               macroread(4) =  .TRUE.

            else if ( iporos .eq. 6 ) then
                if( nspeci .le. 0 ) then
                if(iout.ne.0) then
                   write(iout,*) ""
                   write(iout,*) "Chemistry required for ppor model 6"
                   write(iout,*) ""
                endif
                if(iptty.ne.0) then
                   write(iptty,*) ""
                   write(iptty,*) "Chemistry required for ppor model 6"
                   write(iptty,*) ""
                endif
c               stop
            endif
c           DRH: updates based on ps_delta_rxn for salt
               narrays    =   1
               itype  (1) =   8
               default(1) =   0.0d0
               macro      =   "ppor"
               igroup     =   2
c
               call  initdata2  ( inpt, ischk, n0, narrays, itype,
     *              default, macroread(4), macro, igroup, ireturn,
     1              r8_1=amgang(1:n0) )
c
               macroread(4) =  .TRUE.

c
            else if ( iporos .eq. 7 ) then
c      DRH 1/2/2013: perm-poro relation taken from
c      FCRD-UFD-2013-000064, which references
c      Cinar et al. (2006) Transport in Porous Media 64,
c      pp. 199-228.
c      k = 4.866e-9*ps^4.637
c      valid up to porosity = 0.3
            if( nspeci .le. 0 ) then
                if(iout.ne.0) then
                   write(iout,*) ""
                   write(iout,*) "Chemistry required for ppor model 7"
                   write(iout,*) ""
                endif
                if(iptty.ne.0) then
                   write(iptty,*) ""
                   write(iptty,*) "Chemistry required for ppor model 7"
                   write(iptty,*) ""
                endif
c               stop
            endif
c           DRH: updates based on ps_delta_rxn for salt
               narrays    =   4
               itype  (1) =   8
               itype  (2) =   8
               itype  (3) =   8
               itype  (4) =   8
c              multiplier for exponential function
               default(1) =   0.d0 
c              power for exponential function
               default(2) =   0.d0
c              lower porosity limit for exponential function
               default(3) =   0.d0
c              upper porosity limit for exponential function
               default(4) =   0.d0
               macro      =   "ppor"
               igroup     =   2
c
               call  initdata2  ( inpt, ischk, n0, narrays, itype,
     *              default, macroread(4), macro, igroup, ireturn,
     1              r8_1=porTemp1(1:n0),
     2              r8_2=porTemp2(1:n0),
     3              r8_3=porTemp3(1:n0),
     4              r8_4=porTemp4(1:n0) )
c
               macroread(4) =  .TRUE.
            endif
c
         else if ( iz .eq. 1 )  then
c
****   Calculation handling division   ****
c
            if ( iporos .eq. -1 )  then
c
c specific storage input
c assume specific storage measured at ambient conditions
c T = 20 C, p= 0.1 Mpa
c density(den_w0) = 997.808 kg/m3
c water compressibility(comp_w0) = 5.687e-4 Mpa-1
c see parameter statement
c
               do jji=1,n
                  pso=psini(jji)
                  if (pso.gt.1.e-06) then
                     dpp=max(0.0d00,
     &                    1.d06*amgang(jji)/9.81/den_w0-pso*comp_w0)
                     delp=phi(jji)-phini(jji)
                     ps(jji)=pso+dpp*delp
                     if(ps(jji).gt.por_ll) then
                        dporp(jji)=dpp
                        dport(jji)=0.0
                     else
                        ps(jji) = por_ll
                        dporp(jji)=0.0
                        dport(jji)=0.0
                     endif
                  end if
               enddo
c
            else if ( iporos .eq. -5 )  then
c
c MODFLOW-like specific yield 
c need max and min pressures for each cell 
c delp is the max pressure change over the cell 
c
               do jji=1,n
                  pso=psini(jji)
                  if(amgang(jji).gt.0.0d0) then
                     if (pso.gt.1.e-06) then
                        dpp=max(0.0d00,
     &                       1.d06*amgang(jji)/9.81/den_w0-pso*comp_w0)
                        delp=phi(jji)-phini(jji)
                        ps(jji)=pso+dpp*delp
                        if(ps(jji).gt.por_ll) then
                           dporp(jji)=dpp
                           dport(jji)=0.0
                        else
                           ps(jji) = por_ll
                           dporp(jji)=0.0
                           dport(jji)=0.0
                        endif
                     end if
                  else               
                     if(pso.gt.0.0) then
                        sy = abs(amgang(jji))
                        delp=wgangi(jji)
                        dpp = (phi(jji)-phini(jji)+delp)/delp
                        ps(jji)=sy*dpp
                        dporp(jji)=sy/delp
                        dport(jji)=0.0
                     endif
                  end if
               enddo

            else if ( iporos .eq. 5 )  then
c
****   temperature dependent porosity
c
               do jji=1,n
                  pso = psini(jji)
                  if(pso > 1.e-6) then
                    tmpPor =porTemp1(jji)-porTemp2(jji)/(t(jji)+273.15)
                    dporp(jji) =  0.0
                    dport(jji) =  porTemp2(jji)
                    if(tmpPor < ps(jji) )  ps(jji) = tmpPor
                  endif
               enddo

            else if ( iporos .eq. 1 )  then
c
****   simple aquifer compressibility model   ****
c
               do   jji=1,n
c
                  pso    =  psini(jji)
c
                  if ( pso .gt. 1.e-06 )  then
                     dpp        =  amgang(jji)     
                     delp       =  phi(jji)-phini(jji)
                     ps(jji) =  pso+dpp*delp
                     dporp(jji) =  dpp
                     dport(jji) =  0.0
                  end   if
c
               end   do
c
            else if ( abs(iporos) .eq. 2 )  then
c
****   GZ gangi model   ****
c
               do   jji=1,n
c
                  if ( iad .eq. 0 )  then
                     iporf(jji) =  0
                  end   if
c
                  pso    =  psini(jji)
c
                  if ( ps(jji) .ge. 1.0 )  then
c
                     ps   (jji) =  1.0
                     rlf  (jji) =  1.0
                     dporp(jji) =  0.0
                     dport(jji) =  0.0
                     drlpf(jji) =  0.0
                     drlef(jji) =  0.0
c
                  else if ( pso .gt. 1.0e-06 )  then
c
                     dsigt  =  alphae(jji)*(tini(jji)-t(jji))
                     closur =  (-phi(jji)+phini(jji)+strgan(jji))-dsigt
c
                     if ( closur .lt. 0.0.or.iporf(jji) .ne. 0 )  then
                        iporf(jji) =  1
                        gt         =  1.0
                        xrl        =  1.0
                        dgtp       =  0.0
                        dgtt       =  0.0
                        drl        =  0.0
                     else
                        gt     = (1.0-(closur/pgangi(jji))**amgang(jji))
                        dclosp = -1.0
                        dclost =  alphae(jji)
                        dgtp   = -amgang(jji)*( closur/pgangi(jji) )**
     &                       ( amgang(jji)-1.0 )*dclosp/pgangi(jji)
                        dgtt   = -amgang(jji)*( closur/pgangi(jji) )**
     &                       ( amgang(jji)-1.0 )*dclost/pgangi(jji)
                        xrl    =  gt**3
                        drl    =  3.0*gt**2
                     end  if
c
                     ps   (jji) =  pso*gt
                     dporp(jji) =  pso*dgtp
                     dport(jji) =  pso*dgtt
                     rlf  (jji) =  xrl
                     drlpf(jji) =  drl*dgtp
                     drlef(jji) =  drl*dgtt
c
                  end   if
c
               end   do
c
            else if ( iporos .eq. 3 )  then
c
****   call to sub rock in welbore module ** non-use   ****
c
*           call  welbor  ( 11 )
c
            else if ( iporos .eq. -4 )  then
c
****   Yamashita_Tenma-Gangi model   ****
c
               do   jji=1,n
c
                  if ( wgangi(jji) .gt. 1.0e-04 )  then
c
                     if ( iad .eq. 0 )  then
                        iporf(jji) = 0
                     end   if
c
                     pso    =  psini(jji)
c
                     if ( ps(jji) .ge. 1.0 )  then
c
                        ps   (jji) =  1.0
                        rlf  (jji) =  1.0
                        dporp(jji) =  0.0
                        dport(jji) =  0.0
                        drlpf(jji) =  0.0
                        drlef(jji) =  0.0
c
                     else if ( pso .gt. 1.0e-06 ) then
c
                        dsigt  =  thexp*young*( tini(jji)-t(jji) )
                        closur =(-phi(jji)+phini(jji)+sgangi(jji))-dsigt
c
                        if (closur .lt. 0.0 .or. iporf(jji) .ne. 0) then
                           iporf(jji) =  1
                           gt         =  1.0
                           xrl        =  1.0
                           dgtp       =  0.0
                           dgtt       =  0.0
                           drl        =  0.0
                        else
                           gt = (1.0-(closur/pgangi(jji))**amgang(jji))
                           dclosp = -1.0
                           dclost =  thexp*young
                           dgtp   = -amgang(jji)*(closur/pgangi(jji))**
     &                          ( amgang(jji)-1.0 )*dclosp/pgangi(jji)
                           dgtt   = -amgang(jji)*(closur/pgangi(jji))**
     &                          ( amgang(jji)-1.0 )*dclost/pgangi(jji)
                           xrl    =  gt**3
                           drl    =  3.0*gt**2
                        end   if
c
                        ps   (jji) =  pso*gt
                        dporp(jji) =  pso*dgtp
                        dport(jji) =  pso*dgtt
                        rlf  (jji) =  xrl
                        drlpf(jji) =  drl*dgtp
                        drlef(jji) =  drl*dgtt
                     end   if
                  end   if
               end   do
c            end   if

            else if ( iporos .eq. 6 )  then
c
****   DRH: salt model 1: simple perm/poro relationship based   ****
****   on change in porosity ****
c

               do   jji=1,n
               		if( amgang(jji) .gt. 0.0d0 ) then
                  		pso = psini(jji)
                  		if(pso > 1.e-6) then
                  			ps(jji) =  ps(jji) + ps_delta_rxn(jji)
                  			pnx(jji) = pnxi(jji)*(ps(jji)/psini(jji))**2
               				pny(jji) = pnyi(jji)*(ps(jji)/psini(jji))**2
               				pnz(jji) = pnzi(jji)*(ps(jji)/psini(jji))**2
                  		endif
                  	endif
               end   do
c
           else if ( iporos .eq. 7 )  then
c
****   DRH: salt model 2 - Crushed salt:  ****
c      DRH 1/2/2013: perm-poro relation taken from
c      FCRD-UFD-2013-000064, which references
c      Cinar et al. (2006) Transport in Porous Media 64,
c      pp. 199-228.
c      valid up to porosity = 0.3
c      Scaled by 1.e6 to jive with rest of code.
c      Initialization below for (iz .eq. 3)
c      PHS 8/5/13  Adding check if negative porosity = 0.0
c                  Also removed zeroing out of ps_delta_rxn, done in csolve now.
c
c       PHS 1/20/14  adding check for cut time step days_last. subtract off previous 
c gaz 090113 added por_salt_min and averaging  
               psdelta = 0.0
               psvol = 0.0
               do   jji=1,n
                if( porTemp1(jji) .gt. 0. ) then
                
                    if(days.LT.days_last) then
                        ps_delta_rxn(jji) =
     &                      ps_delta_rxn(jji) - ps_delta_rxn_s(jji)
                    end if  
c                    
                 pso = ps(jji)
c                 
                 ps(jji) =  ps(jji) + ps_delta_rxn(jji)
                 if(ps(jji).LT.por_salt_min) ps(jji) = por_salt_min
c                 Calculate total change in porosity and volume of
c                 ppor 7 model nodes
c
                  pso = ps(jji)
c                  
                  psdelta = psdelta + (ps(jji) - psini(jji))*volume(jji)
                  psvol = psvol + volume(jji)
c                calculate permeabilitiy for new porosity
c
                   if(iaprf.eq.0) then
c
c use Cinar relationship
c
                    if(pso>=porTemp3(jji).and.pso<=porTemp4(jji)) then
                     pnx_new = porTemp1(jji)
     1                   *ps(jji)**porTemp2(jji)*1.e6
                    else if(pso > porTemp4(jji)) then
c                    Set to upper limit
                     pnx_new = porTemp1(jji)
     1                   *porTemp4(jji)**porTemp2(jji)*1.e6
                    else if(pso < porTemp3(jji)) then
c                    Set to lower limit
                     pnx_new = porTemp1(jji)
     1                   *porTemp3(jji)**porTemp2(jji)*1.e6
                    endif
                     pnx(jji) = pnx_new
                     pny(jji) = pnx(jji)
                     pnz(jji) = pnx(jji)   
                    else
c use Olivella relationship
c added another call to evaluate a single gridblock value
                    call perm_olivella(2,jji)
                   endif         
                endif
               end   do
 666    format(4G12.6)
               days_last = days
           endif
c
c
         else if ( iz .eq. 2 )  then
c additional printout for variable porosity
            if ( iporos .ne. 3 )  then
c
****   printout varible porosity information   ****
c
               if ( iporos .eq. -1 .or. iporos .eq. -5 )  then
                  if (iout .ne. 0 .and. ntty.eq.2) write(iout, 6119)
                  if (iatty .gt. 0) then
                     write(iatty,*)
                     write(iatty, 6119)
                  end if
c
                  do   i=1,m
c
                     md     =  nskw(i)
                     permsb =  pnx(md)*1.0e-6
                     delpsb =  pho(md)-phini(md)
c
                     if (iout .ne. 0 .and. ntty.eq.2) 
     &                    write(iout,6120) md, permsb, ps(md), delpsb
                     if (iatty .gt. 0) 
     &                    write(iatty,6120) md, permsb, ps(md), delpsb
c
                  end  do
 6119             format(1x,'  Node', 3x, 'permeability', 3x, 
     &                 'porosity', 7x, 'pressure change')
 6120             format(1x,i7,3(3x,g12.6))
               else if ( iporos .eq. 5 ) then
               else if ( iporos .eq. 1 ) then
               else if ( iporos .eq. -1 ) then

c   - - - - PHS 7/17/13 added output to screen and out file for model 7 SALT
               else if  ( iporos .eq. 7 ) then

                  if (iout .ne. 0) write(iout,6015)
                  if (iatty .ne. 0) write(iatty,6015)

                  if (iout .ne. 0) write(iout,6016)
                  if (iatty .ne. 0) write(iatty,6016)

                  do   i=1,m

                     md     =  nskw(i)
                     permsb =  pnx(md)*1.0e-6
                     pwv    =  phi(md) - pci(md)

                     if (iout .ne. 0)     
     &                  write(iout,6017) md,permsb,ps(md),thx(md)*1e6,
     &                     pwv,dvas(md),ps_delta_rxn_s(md)

                     if (iatty .ne. 0)
     &                  write(iatty,6017) md,permsb,ps(md),thx(md)*1e6,
     &                     pwv,dvas(md),ps_delta_rxn_s(md)

                  end  do
                 if(iout.ne.0) then 
                  write(iout,*) 'Total change in volume: ', psdelta,
     &                 ' m'
                  write(iout,*) 'Percent change in total volume: ',
     &                 psdelta/psvol*100, ' %'
                  write(iout,*) 'Total Volume involved in ppor  ',
     &               psvol
                 endif
                 if(iatty.ne.0) then 
                 write(iatty,*) 'Total change in volume: ', psdelta,
     &                 ' m'
                 write(iatty,*) 'Percent change in total volume: ',
     &                 psdelta/psvol*100, ' %'
                 write(iatty,*) 'Total Volume involved in ppor  ',
     &               psvol
                 endif



 6015        format(' - - - - - - - - - - - - - - - - - - - - - - - - ',
     &       '- - - - - - - - - - - - - - - - -')
 6016             format(1x,'  Node', 5x, 'perm (m2)', 6x,
     &                 'porosity', 7x, 'Kx W/(m K)', 5x,'Pwv (MPa)',
     &                6x, 'D*wv (m2/s)',2x,' ps_delta_rxn')
 6017             format(1x,i7,6(3x,g12.5))


              else if ( iporos .ne. -4.AND.iporos .ne.7 )  then
c
****   GZ gangi model   ****
c
                  if (iout .ne. 0) write(iout,6019)
c
                  do   i=1,m
c
                     md     =  nskw(i)
                     permsb =  pnx(md)*1.0e-6*rlf(md)
                     dsigt  =  thexp*young*( tini(md)-t(md) )
                     closur =  sigini-pho(md)-dsigt
                     delpsb =  pho(md)-pein
c
                     if (iout .ne. 0)
     &                    write(iout,6020) md,permsb,ps(md),dsigt,
     &                    closur,delpsb
c
                  end  do
 6019             format(1x,'  Node', 3x, 'permeability', 3x, 
     &                 'porosity', 7x, 't stress', 7x, 'closure s',
     &                 6x, 'dl pressure')
 6020             format(1x,i7,5(3x,g12.6))
c
               else if ( iporos .eq. -4 )  then
c
****   Yamashita_Tenma-Gangi model   ****
c
                  do   i=1,m3
c
                     md     =  nskw3(i)
                     dsigt  =  thexp*young*( tini(md)-t(md) )
                     closur =  ( -phi(md)+phini(md)+sgangi(md) )-dsigt
                     delpsb =  pho(md)-pein
c
                     if ( closur .lt. 0.0 .or. iporf(md) .ne. 0 )  then
                        yama_ww =  wgangi(md)*dexp(dlog(1.0/wgangi(md))*
     &                       ( phi(md)/( phini(md)+sgangi(md) ) ) )
                        yama_ww = min (0.008d0, yama_ww)
                     else
                        gt      =  1.0-( closur/pgangi(md) )**amgang(md)
                        yama_ww =  wgangi(md)*gt
                     end   if
c
                     if ( yama_ww .gt. 0.0 )  then 
                        yama_ff =  1.0+17.0*( agangi(md)/
     &                       ( 2.0*yama_ww ) )**1.5
                        permsb  =  yama_ww**2/( 12.0*yama_ff )
                     else
                        permsb  =  0.0
                     end   if
c
c                write(yama1,3100)  l, days, md, permsb, pho(md), yama_ww
c 3100           format(i10,g15.7,i10,3g15.7)
                     write(yama1,3100)  l, days, md, permsb, pho(md), 
     .                    yama_ww, tenma_ww(md)
 3100                format(i10,g15.7,i10,4g15.7)
c
                  end   do
c
                  if(m3.gt.0) then
                     write(yama2,3200)  l, days, (pho(nskw3(i)),i=1,m3)
                  endif
 3200             format(i10,g15.7,250f7.3)
c
               end  if
c
            else
c
*           call welbor(12)
c
            end  if    !  (iz = 2)
c
         else if ( iz .eq. 3 )  then
c
**** calculate initial permeabilities and porosities
c
            if ( iporos .eq. -2 )  then
c
****   GZ gangi model   ****
c
               do   jji=1,n
c
                  if ( iad .eq. 0 )  then
                     iporf(jji) =  0
                  end   if
c
                  pso    =  psini(jji)
c
                  if ( ps(jji) .ge. 1.0 )  then
c
                     ps (jji) =  1.0
                     rlf(jji) =  1.0
c
                  else if ( pso .gt. 1.0e-06 )  then
c
                     dsigt  =  alphae(jji)*( tini(jji)-t(jji) )
                     closur =  (-phi(jji)+phini(jji)+strgan(jji))-dsigt
c
                     if ( closur .lt. 0.0 .or. iporf(jji) .ne. 0 )  then
                        iporf(jji) =  1
                        gt         =  1.0
                        xrl        =  1.0
                     else
                        gt         =  ( 1.0-
     &                       ( closur/pgangi(jji) )**amgang(jji) )
                        xrl        =  gt**3
                     end   if
c
                     psini(jji) =  pso/gt
                     pnx  (jji) =  pnx(jji)/xrl
                     pny  (jji) =  pny(jji)/xrl
                     pnz  (jji) =  pnz(jji)/xrl
c
                  end   if
               end   do
c
            else if ( iporos .eq. 5 )  then
c              temperature-dependent porosity
c              added by ZL on 8/15/2008
               do jji=1,n
c                 (1) a linear function of temperature
                  pso = psini(jji)
                  if(pso > 1.0e-6) then
                    tmpPor=porTemp1(jji)-porTemp2(jji)/(t(jji)+273.15) 
                    if(tmpPor < ps(jji) )  ps(jji) = tmpPor
                  endif
c                  if ( ps(jji) .ge. 1.0 )  ps (jji) =  1.0
                  psini(jji) = ps(jji)
               enddo

****   Yamashita_Tenma-Gangi model   ****
c
            else if ( iporos .eq. -4 )  then
c
               if (iout .ne. 0) 
     &              write(iout,*)'=============  pnx ============='
c
               do   jji=1,n
c
                  if ( wgangi(jji) .gt. 1.0e-04 )  then
c
                     if ( iad .eq. 0 )  then
                        iporf(jji) =  0
                     end   if
c
                     pso     =  psini(jji)
                     yama_ww =  0.0
                     yama_ff =  0.0
                     gt      =  0.0
c
                     if ( ps(jji) .ge. 1.0 )  then
c
                        ps (jji) =  1.0
                        rlf(jji) =  1.0
                        yama_phi =  (  phi(jji)-phini(jji) )/sgangi(jji)
c
                        if ( yama_phi .gt. 1.0 )  then
                           closur  =  (-phi(jji)+phini(jji)+sgangi(jji))
     &                       -dsigt
                           yama_ww =  wgangi(jji)*
     &                          dexp( dlog( 1.0/wgangi(jji) )*yama_phi )
                           yama_ww = min (0.008d0, yama_ww)
                        else
                           dsigt   =  thexp*young*( tini(jji)-t(jji) )
                           closur  =  (-phi(jji)+phini(jji)+sgangi(jji))
     &                       -dsigt
                           gt      =  1.0-(closur/pgangi(jji))**
     &                          amgang(jji)
                           xrl     =  gt**3 
                           yama_ww =  wgangi(jji)*gt
                        end   if
c
                        if ( gt .gt. 1.0e-20 )  then
                           yama_ff  =  1.0+17.0*( agangi(jji)/
     &                          ( 2.0*yama_ww ) )**1.5
                           pnx(jji) =  yama_ww**2/(12.0*yama_ff)*1.0d+06
                        end   if
c
                     else if ( pso .gt. 1.0e-06 )  then
c     
                        dsigt    = thexp*young*( tini(jji)-t(jji) )
                        closur  =  (-phi(jji)+phini(jji)+sgangi(jji))
     &                       -dsigt
c
                        if (closur .lt. 0.0 .or. iporf(jji) .ne. 0) then
                           iporf(jji) =  1
                           gt         =  1.0
                           xrl        =  1.0
                           yama_ww    = wgangi(jji) *
     &                          dexp( dlog( 1.0/wgangi(jji) )*
     &                          ((-phi(jji)+phini(jji) )/sgangi(jji)))
                           yama_ww = min (0.008d0, yama_ww)
                        else
                           gt    =  1.0-( closur/pgangi(jji) )**
     &                          amgang(jji)
                           xrl   =  gt**3 
                           yama_ww  =  wgangi(jji)*gt
                           yama_ww = min (0.008d0, yama_ww)
                        end   if
c  
                        if ( gt .gt. 1.0e-20 )  then
                           if (yama_ww .ge. tenma_ww(jji)+1.0d-06) then
                              psini(jji)  = min( pso/gt,1.0d0 )
                           else 
                              psini(jji)  = 0.01
                           end if
                           yama_ff  = 1.0+17.0*( agangi(jji)/
     &                          ( 2.0*yama_ww ) )**1.5
                           pnx(jji) =  (yama_ww/0.008)*yama_ww**2/
     &                          ( 12.0*yama_ff )*1.0d+06
                           pny(jji) =  pnx(jji)
                           pnz(jji) =  pnx(jji)
                        end   if
c
                     end   if
c
c                     write(18,9999)  jji,  ps(jji), yama_ww, pnx(jji)
c Adding this write to yama1_file (?)
                     write(yama1,9999)  jji,  ps(jji), yama_ww, pnx(jji)
 9999                format(' pnx',i6,' ',3g9.3)
c
                  end   if
c
               end   do
c
c
            else if ( iporos .eq. 6 )  then
c              If model 6, set copy initial perms
               	pnxi = pnx
               	pnyi = pny
               	pnzi = pnz

c              If model 7, initialize perms
            else if ( iporos .eq. 7.and.iaprf.ne.0) then
c gaz 020216
c if model 7 and olivella prem model get initial perms
              call perm_olivella(1,0)
c
            else if ( iporos .eq. 7 ) then
             do   jji=1,n
              if( porTemp1(jji) .gt. 0. ) then
               pso = ps(jji)
               if(pso>=porTemp3(jji).and.pso<=porTemp4(jji)) then
                pnx(jji) = porTemp1(jji)*ps(jji)**porTemp2(jji)*1.e6
                pny(jji) = pnx(jji)
                pnz(jji) = pnx(jji)
               else if(pso > porTemp4(jji)) then
c                    Set to upper limit
c  gaz debug 082713 (line longer than 72-OK sometimes?)
c
                pnx(jji) = 
     &          porTemp1(jji)*porTemp4(jji)**porTemp2(jji)*1.e6
                pny(jji) = pnx(jji)
                pnz(jji) = pnx(jji)
               else if(pso < porTemp3(jji)) then
c                    Set to lower limit
                pnx(jji) = 
     &          porTemp1(jji)*porTemp3(jji)**porTemp2(jji)*1.e6
                pny(jji) = pnx(jji)
                pnz(jji) = pnx(jji)
               endif
              endif
             enddo
            endif
         else if(iz.eq.4) then
c 
c set inital pressures and temperatures (only)
c         
           phini = pho
           psini = ps
           tini = to        
         else  if ( iz .eq. 5 )  then

c 
c calculate intial quantities for model -5 if enabled
c
c must call iz= 4 first to set up phini and tini
c
            if ( iporos .eq. -5 )  then
c calculate max pressure change for specific yield model -5 
               pho = phini
               phi = phini
               do i = 1,neq
                  if(amgang(i).lt.0.0d0) then
                     i1=nelm(i)+1
                     i2=nelm(i+1)
                     hmid=cord(i,igrav)
                     hmin=0.
                     hmax=0.
                     do ii =i1,i2
                        kb=nelm(ii)
                        hmax=max(cord(kb,igrav)-hmid,hmax)
                        hmin=min(cord(kb,igrav)-hmid,hmin)
                     enddo	         
                     if(ivf.eq.-1) then
                        alength = max(hmax,abs(hmin))
                     else
                        alength = abs(hmax-hmin)/2.	            
                     endif 
                     wgangi(i) = alength*rho1grav
                  endif	                     
               enddo
            endif 
         endif
c
      end   if
c
      return
      end

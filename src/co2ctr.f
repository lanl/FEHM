      subroutine co2ctr(iflg)
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
CD1 To provide overall control for a nonisothermal air-water simulation.
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
CD2 $Log:   /pvcs.config/fehm90/src/co2ctr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:58   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:14   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.10   Wed Jun 12 15:20:58 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2
CD2    Rev 1.9   Mon Jun 03 11:17:44 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.8   Fri May 31 10:43:08 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.7   Thu Feb 15 09:25:12 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.6   Mon Jan 29 13:39:04 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.5   03/23/95 18:19:40   gaz
CD2 gaz added provision for ico2 to enable irdof
CD2 
CD2    Rev 1.4   03/10/95 10:28:40   llt
CD2 humidity BC added - gaz
CD2 
CD2    Rev 1.3   01/28/95 13:54:04   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.2   03/18/94 15:45:16   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   01/21/94 11:26:58   llt
CD2 ngas was not getting read - Zora corrected
CD2 
CD2    Rev 1.0   01/20/94 10:21:46   pvcs
CD2 original version in process of being certified
CD2 
c 1/9/95 gaz added specified source for ngas(used array dqpc)
c 1/9/95 gaz added humidity BC 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier    Type     Description
CD3 
CD3 iflg          int      Control flag for determining the purpose
CD3                           for calling this routine.
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 Note: the files below are identified by their tape numbers
CD3 
CD3 Name          Description
CD3 
CD3 inpt          Input file
CD3 iout          Output file
CD3 iptty         Alternate output file
CD3 iatty         Alternate output file
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
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 ico2, qtc, qtotc, amc, inpt, ischk, n, ippci, ipeskc, pci, ieos,
CD4 iout, iptty, pcp, dpcef, pho, to, tini, t, denpch, denpci,
CD4 dtot, volume, pcio, iac, acner, qc, dtotdm, sk, denpcj, idpdp, m,
CD4 idualp, nskw, sk, phi, bpc, difc, iatty, bp,narrays, pointer, itype,
CD4 default, igroup, ireturn,
CD4 macroread
CD4 
CD4 Global Subprograms
CD4
CD4 Name      Type     Description
CD4 
CD4 initdata  N/A      Reads input values and set values at nodes
CD4 psatl     real*8   Computes water vapor pressure
CD4 exit      N/A      Closes files, stops execution
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5 
CD5 tbnd         real*8      Lower limit for gas pressure
CD5 max_inputs   int         Maximum number of input arrays being read
CD5                              on a single line
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 ico2d        int         Input value of gas flag
CD5 i            int         Do loopp parameter over all nodes
CD5 pcid         real*8      Current value of gas pressure
CD5 istflag      int         Flag for determining if routine ended
CD5                              with error
CD5 pv           real*8      Water vapor pressure
CD5 dtsatp       real*8      Unused parameter returned from call to
CD5                              psatl
CD5 dpsats       real*8      Unused parameter returned from call to
CD5                              psatl
CD5 dpsatt       real*8      Unused parameter returned from call to
CD5                              psatl
CD5 tdumm        real*8      Temperature returned from call to psatl
CD5 mi           int         Do loop parameter over all nodes
CD5 it           int         Parameter used in referencing rlperm model
CD5 hum          real*8      Humidity
CD5 alp          real*8      Parameter in relative perm model
CD5 beta         real*8      Parameter in relative perm model
CD5 sr           real*8      Residual liquid saturation       
CD5 smax         real*8      Maximum liquid saturation       
CD5 qtcd         real*8      Parameter used in calculation
CD5 qtotci       real*8      Total mass injected in source term
CD5 dencht       real*8      Last time step value of mass accumulation
CD5                              at the current node
CD5 mlev1        int         Position in output array where first set
CD5                              of matrix nodes starts
CD5 mlev2        int         Position in output array where second set
CD5                              of matrix nodes starts
CD5 mdd          int         Current printout node
CD5 md           int         Current printout node
CD5 rqd          real*8      Source strength at the current node
CD5 qcd          real*8      Gas fraction in the source at current node
CD5 bpd          real*8      Term used in calculation
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 N/A
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
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations 
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
CPS BEGIN co2ctr
CPS 
CPS Set error flag to indicate no errors
CPS IF this is a two-phase simulation
CPS 
CPS   IF this call is for reading input
CPS     Read gas flag
CPS     
CPS     IF the input gas flag indicates no gas in simulation
CPS       Set flag
CPS     ELSE
CPS       Set flag
CPS     ENDIF
CPS     
CPS     Set input parameters for reading gas pressure
CPS     
CPS     initdata - read gas pressures and set values at nodes
CPS     Set input parameters for reading gas source strength
CPS     initdata - read gas source strengths and set values at nodes
CPS     
CPS     Set flag to denote that the ngas macro has been called
CPS     
CPS   ELSEIF this call is for checking data
CPS   
CPS     FOR each node
CPS       IF the gas pressure is less than the limiting value
CPS         IF this is not a two-phase simulation
CPS           Write error message
CPS           ERROREXIT
CPS         ENDIF
CPS         psatl - compute vapor pressure term
CPS         Compute new gas pressure
CPS         IF the gas pressure is less than the minimum value
CPS           Write error message
CPS           ERROREXIT
CPS         ENDIF
CPS         Set gas pressure in array
CPS       ELSEIF this is a two-phase simulation
CPS         IF the gas pressure is greater than the total pressure
CPS           Write error message
CPS           ERROREXIT
CPS         ELSEIF the gas pressure is less than the minimum value
CPS           Write error message
CPS           ERROREXIT
CPS         ENDIF
CPS       If the gas pressure is a specified boundary condition
CPS        Change the temperature if in two-phase region
CPS        Write a message saying what was done
CPS       ENDIF
CPS       If there is a gas flow term at a node with no flow specified
CPS       Identify the node as having sources
CPS       ENDIF
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is for computing storage terms at each time step
CPS   
CPS     FOR each node
CPS       Compute saturation for given humidity
CPS ????????? have to finish ???????????????????
CPS     ENDFOR
CPS   
CPS     FOR each node
CPS       Compute storage term
CPS       Set old gas pressure
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is for computing other storage terms
CPS   
CPS     FOR each node
CPS       Compute other storage terms
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is for computing mass balance parameter
CPS   
CPS     Compute mass balance parameter
CPS     
CPS   ELSEIF this call is for writing output
CPS   
CPS     IF writing is to be performed to first output file
CPS     
CPS       IF this is a dpdp solution
CPS         Set parameters used for writing output at certain nodes
CPS       ELSEIF this is a dual porosity simulation
CPS         Set parameters used for writing output at certain nodes
CPS       ELSE this is single porosity
CPS         Set parameters used for writing output at certain nodes
CPS       ENDIF
CPS       
CPS       FOR each output node
CPS         IF this node is is for matrix level 2
CPS           Write header to file
CPS         ENDIF
CPS         IF this node is is for matrix level 3
CPS           Write header to file
CPS         ENDIF
CPS         Compute mass balance terms
CPS         Write data to file
CPS       ENDFOR
CPS       
CPS       Write global mass and energy balance information
CPS     
CPS     ENDIF
CPS     IF writing is to be performed to second output file
CPS     
CPS       IF this is a dpdp solution
CPS         Set parameters used for writing output at certain nodes
CPS       ELSEIF this is a dual porosity simulation
CPS         Set parameters used for writing output at certain nodes
CPS       ELSE this is single porosity
CPS         Set parameters used for writing output at certain nodes
CPS       ENDIF
CPS       
CPS       FOR each output node
CPS         IF this node is is for matrix level 2
CPS           Write header to file
CPS         ENDIF
CPS         IF this node is is for matrix level 3
CPS           Write header to file
CPS         ENDIF
CPS         Compute mass balance terms
CPS         Write data to file
CPS       ENDFOR
CPS       
CPS       Write global mass and energy balance information
CPS     
CPS     ENDIF
CPS     
CPS   ENDIF
CPS 
CPS ENDIF this is a two-phase simulation
CPS 
CPS ERRORSEGMENT
CPS   IF an error occurred in this routine
CPS     exit - stop execution of the code
CPS   ENDIF
CPS ENDSEGMENT
CPS 
CPS END co2ctr
CPS 
C**********************************************************************
c
c subroutine to control noncondensible gas simulation
c
c code written by g zyvoloski feb 1980

      use comci
      use combi
      use davidi
      use comdi
      use comei
      use comii
      use comhi
      use comgi
      use comfi
      use comdti
      use comai
      use comki
      use comxi
      implicit none 

       real*8 bpd,tbnd,pcid,pv,dtsatp,dpsats,dpsatt,tdumm
       integer iflg,ico2d,i,istflag,mi,it
       real*8 hum,alp,beta,sr,smax,qtcd,dencht,tol_p
       real*8 hum_frac
       real(8) :: qtotci = 0.
       integer mlev1,mlev2,mdd,md
       real*8 rqd,qcd,psatl,qcmax       
       real*8 sat_bc, tol_sat_bc
       integer ngas_flag2
       logical ngas_flag
       character*80 dumstring
       integer msg(20)
       integer nwds
       real*8 xmsg(20)
       integer imsg(20)
       character*32 cmsg(20)
       logical :: old_input = .false.
       character*90 wdum
c gaz 102419 changed tol_p  from 1.e-6    
c applying tol_p to mass frac calc       
       parameter (tol_p = 1.d-8)
       parameter (tol_sat_bc = 1.d-10)
       real*8, allocatable :: aiped(:)      
c     set tbnd for pco2 change(see about line 580)
       parameter(tbnd = -1.0)
       save ngas_flag, ngas_flag2
c     
c     return if no noncondensible present
c     
      istflag = 0
      if(ico2.gt.0) then
c     iflg is used to tell if call is for initialization
         if(iflg.eq.0) then
c
c check macro line for more information
c
          backspace inpt
            ngas_flag2 = 0
            read(inpt,'(a90)') wdum 
            do i = 5,80
             if(wdum(i:i+5).eq.'normal') ngas_flag2 = 0
             if(wdum(i:i+6).eq.'reset T') ngas_flag2 = 1
             if(wdum(i:i+6).eq.'reset P') ngas_flag2 = 2
            enddo
            if(ngas_flag2.eq.1) then
             if(iout.ne.0) write(iout,12)
             if(iptty.ne.0) write(iptty,12)          
            elseif(ngas_flag2.eq.2) then 
             if(iout.ne.0) write(iout,13)
             if(iptty.ne.0) write(iptty,13)     
            endif
12     format('temperature reset to h2o vapor pressure at total', 
     &        'initial pressure')
13     format('total pressure reset to h2o vapor pressure at ', 
     &        'initial temperature if h2o vapor pressure is greater',
     &          'than initial total presssure')
            qtc=0.0
            qtotc=0.0
            amc=0.0
            ngas_flag = .false.
            read(inpt,*) ico2d
c     dont call co2ctr if ico2d=0
            if(ico2d.eq.0) then
               ico2=0
            else
               if(ico2d.eq.2) then
                  if(irdof.eq.0) irdof=2
                  if(icoupl.eq.0) icoupl=5
               else if(ico2d.eq.1) then
                  if(irdof.eq.0) irdof=1
                  if(icoupl.eq.0) icoupl=5
               endif
               
               ico2=3              
            endif
c
c allocate impedance array for noncondensible gas and humidity for output
c
            allocate (wellima(n0))
            allocate (humida(n0))
c     
c     read in initial ncon pressure
c     
            narrays = 1
            itype(1) = 8
            default(1) = 0.
            macro = "ngas"
            igroup = 2
            
            call initdata2(inpt,ischk, n0, narrays, itype,
     &           default, macroread(2), macro, igroup, ireturn,
     &           r8_1=pci(1:n0))
c pci () = -666 equal reset water vapor pressure
c gaz 070916 -1 < pci < 0. set initial relative humidity             
            do i = 1, n0
               if (pci(1) .eq. -666) then
                  ngas_flag = .true.
                  exit
               end if
            end do

c     Check to verify if new input is being used for group 3
            read (inpt, '(a80)') dumstring
            backspace inpt
            call parse_string2(dumstring,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .gt. 0 .and. nwds .lt. 5) then
               old_input = .true.
               write (ierr, 14)
               if (iout .ne. 0) write (iout, 14)
               if (iptty .ne. 0) write(iptty, 14)
            end if
 14         format ('**** Warning: Old style ngas input being used, ',
     &          'macro input should be updated ****')
 
            allocate (pflowa(n0))

            if (old_input) then
c     
c     read in specified pressure for  noncon gas
c     
               igroup = 3
c     
c     Other values are the same as above
c     
            
            call initdata2(inpt,ischk, n0, narrays, itype,
     &           default, macroread(2), macro, igroup, ireturn,
     &              r8_1=eskc(1:n0))

               do i = 1, n0
                  if (eskc(i) .le. 0.) then
c     Specified noncondnesible gas pressure
                     pflowa(i) = abs (eskc(i))
                     eskc(i) = 0.
                  else
c     Specified relative humidity
                     pflowa(i) = 0.0
                  end if
               end do
            else
c     
c     read in fixed humidity for noncon gas(positive value)
c     read in fixed saturation for noncon gas(negative value-saturation is fixed to absolute read-in value)
c
            default(1) = 1.d50
            default(2) = 1.d50
               narrays = 2
               itype(1) = 8
               itype(2) = 8 
               igroup = 3
c     
c     Other values are the same as above
c     
            
               call initdata2(inpt,ischk, n0, narrays, itype,
     &              default, macroread(2), macro, igroup, ireturn,
     &              r8_1=eskc(1:n0),r8_2=pflowa(1:n0))

               do i = 1,n0
                  if(eskc(i).eq.default(1)) then
                     eskc(i) = 0.0
                  endif
                  if(pflowa(i).eq.default(2)) then
                     pflowa(i) = 0.0
                  endif
               enddo

            end if
c     Check to verify if new input is being used for group 4
            if (.not. old_input) then
               read (inpt, '(a80)') dumstring
               backspace inpt
               call parse_string2(dumstring,imsg,msg,xmsg,cmsg,nwds)
               if (nwds .gt. 0 .and. nwds .lt. 5) then
                  old_input = .true.
                  write (ierr, 14)
                  if (iout .ne. 0) write (iout, 14)
                  if (iptty .ne. 0) write(iptty, 14)
               end if
            end if
c gaz debug  120814 
            allocate (aiped(n0))

            if (old_input) then
c     
c     read in source strength for noncon gas(kg/s)
c     
               default(1) = 1.d50
               narrays = 1
               igroup = 4
c     
c     Other values are the same as above
c     
            
               call initdata2(inpt,ischk, n0, narrays, itype,
     &              default, macroread(2), macro, igroup, ireturn,
     2              r8_1=dqpc(1:n0))

c
               do i=1,n0
                  if(dqpc(i).eq.default(1)) then
                     dqpc(i)=0.0
                  elseif(dqpc(i).eq.0.0) then
                     dqpc(i)=1.d-30
                  endif
               enddo
c     The impedance, by default is set to 0.
               qng = dqpc
               aiped = 0.

            else
c     
c     read in source strength for noncon gas(kg/s),air fraction, impendance
c     
               default(1) = 0.0
               default(2) = 0.0
               narrays = 2
               itype(1) = 8
               itype(2) = 8
c
               igroup = 4
c     
c     Other values are the same as above
c  
c         
               call initdata2(inpt,ischk, n0, narrays, itype,
     &              default, macroread(2), macro, igroup, ireturn,
     &              r8_1=qng(1:n0),r8_2=aiped(1:n0))

            end if

            macroread(2) = .TRUE.
c gaz 060820 new conditions
c if aiped ne.0.0, then imped associated   with pflowa          
            if(.not.allocated(xairfl)) allocate (xairfl(n0))   
             xairfl = 0.0
            do i=1,n0
               if(qng(i).eq.default(1)) then
                  qng (i)=0.0
               elseif(qng(i).eq.0.0) then
                  qng(i)=1.d-30
               else
                  xairfl(i) = 1.0  
               endif
               if(aiped(i).eq.default(2)) then
                  wellima(i)=0.0
                  xairfl(i) = 0.0
               elseif(aiped(i).eq.0.0) then
                  wellima(i)=0.0
                  xairfl(i) = 0.0
               else
                  xairfl(i) = qng(i)
                  qng(i) = 0.0
                  wellima(i) = aiped(i)*1.d06
               endif
            enddo
c
c gaz new allocate an if ngas enabled (031515)
c this is so salt vapor lowering can be consistent with psatl
c
c gaz 093017 (might need next 2 lines)
c           deallocate(an)  
c           allocate(an(max(n0,nspeci*n0)))            

           deallocate (aiped)
 
         elseif(iflg.eq.6) then
c     
c     check for temperature input and ncondensible pressure
c  
c gaz 110314
c     leave subroutine if iread.ne.0 (variables will be consistent)
c
         if(iread.ne.0) return   
            do i=1,n
               pcid=pci(i)
               if(pcid. eq. -999. .or. pcid .eq. -666.) then
                pcid=-to(i)
               endif            
               if(ieos(i).eq.3) then
c check for relative humidity IC
                  if(pcid.lt.0.0d0.and.abs(pcid).gt.1.0d0) then
                     write (ierr, 100)
                     if (iout .ne. 0) write(iout,100) 
                     if (iptty .ne. 0) write(iptty, 100) 
                     istflag = -1
                     goto 9000
 100                 format ('cannot input ngas temp in single phase')
c gaz 121218 correction to allow pci < 1.0 but > 0                      
                  else if(pcid.lt.0.0.and.abs(pcid).le.1.0d0
     &                .and.abs(pcid).ge.0.0d0)then
c rel humidity IC and calculate new pci
c 100 percent hmidity water vapor pressure
                   pv= psatl(to(i),pcp(i),dpcef(i),dpsatt,dpsats,0,
     &                        an(i))
                   humida(i) = abs(pcid)
                   pci(i) = pho(i)-abs(pcid)*pv
                  endif
c negative 0<= pci <= 1. indicates relative permeability   
c gaz 122018 need to allow small pci()                
                 else if(ieos(i).eq.2.and.abs(pcid).ge.1.0.
     %              and.pcid.lt.-tbnd) then
c set water vapor pressure                    
                  if(ngas_flag) then
                   pv = 0.0
                   pci(i)=pho(i)-pv
                  else
                    pv= psatl(-pcid,pcp(i),dpcef(i),dpsatt,dpsats,0,
     &                        an(i))
                    if (ngas_flag2.eq.1) then
                     tdumm=psatl(pho(i)-tol_p,pcp(i),dpcef(i),dtsatp,
     &               dpsats,1,an(i))
                     pv = pho(i)
                     pci(i)=tol_p
                     to(i) = tdumm
                    elseif (ngas_flag2.eq.2.and.pv.ge.pho(i)) then
                     pci(i) = tol_p
                     pho(i) = pv + tol_p
                     phi(i) = pho(i)
                     to(i) = -pcid
                    else
                     pci(i)=pho(i)-pv 
                     to(i) = -pcid              
                    end if
                  endif
                  if(pci(i).lt.tbnd) then
                     write(ierr, 110)
                     if (iout .ne. 0) write(iout, 110)
                     if (iptty .ne. 0) write(iptty, 110)
                     tdumm=psatl(pho(i),pcp(i),dpcef(i),dtsatp,dpsats,
     &                      1,an(i))
                     write(ierr, 120) tdumm
                     if (iout .ne. 0) write(iout, 120) tdumm
                     if (iptty .ne. 0) write(iptty, 120) tdumm
                     istflag = -1
                     goto 9000
 110                 format ('ngas pressure lt 0 at temp and total ',
     &                    'pressure given')
 120                 format ('max allowable temperature ', g15.4)
                  endif
c                  to(i)=-pcid
               else if(ieos(i).eq.2) then
c================================================================
c      PHS  9/1/2006   Modification so that code does not 
c           crash on restart.  Now assumes that total P is good
c           and resets air pressure to be total - water vapor
c================================================================
                  if(pcid.gt.pho(i)) then
c                     write(ierr, 130) i
c                     if (iout .ne. 0) write(iout, 130) i
c                     if (iptty .ne. 0) write(iptty, 130) i
	               pv= psatl(t(i),pcp(i),dpcef(i),dpsatt,dpsats,
     &                            0,an(i))
	               pci(i) = pho(i) - pv
c                     istflag = -1
c                     goto 9000
c 130                 format ('ngas pressure gt total pressure at i=',
c     &                    i8)
                  else if(pcid.lt.tbnd) then
                     write(ierr, 140) 
                     if (iout .ne. 0) write(iout, 140) 
                     if (iptty .ne. 0) write(iptty, 140)
                     istflag = -1
                     goto 9000
 140                 format ('ngas pressure lt 0.')
                  else if(ngas_flag) then 
                     pv= psatl(t(i),pcp(i),dpcef(i),dpsatt,dpsats,
     &                         0,an(i))
                     pho(i)= pv+ pcid   
                     phi(i) = pho(i)
                     to(i)=psatl(pv,pcp(i),dpcef(i),dtsatp,dpsats,
     &                         1,an(i))
                  else
                     pv=pho(i)-pcid
                     to(i)=psatl(pv,pcp(i),dpcef(i),dtsatp,dpsats,
     &                          1,an(i))
                  end if
               endif

               tini(i)= to(i)
               t(i)=to(i)
            enddo
         elseif(iflg.eq.1) then
c     
c     figure out humidity condition (VG cap pressure models only)
c     
            do i=1,n
               if(eskc(i).gt.0.0) then
                  it=irlp(i)
                  if(irlpt(it).eq.3.or.irlpt(it).eq.5) then 
                     hum=eskc(i)
                     alp=rp3f(irlp(it))
                     beta=rp4f(irlp(it))
                     sr=rp1f(irlp(it))
                     smax=rp2f(irlp(it))
                     call humidity(hum,alp,beta,sr,smax,to(i),s(i))
                  else if(irlpt(it).eq.4.or.irlpt(it).eq.6.or
     &                       .irlpt(it).eq.7) then
                     if(idpdp.eq.0) then
c     matrix values
                        hum=eskc(i)
                        alp=rp3f(irlp(it))
                        beta=rp4f(irlp(it))
                        sr=rp1f(irlp(it))
                        smax=rp2f(irlp(it))
                        call humidity(hum,alp,beta,sr,smax,to(i),s(i))
                     else
                        if(i.le.neq) then
c     use fracture values for 1-neq nodes
                           hum=eskc(i)
                           alp=rp10f(irlp(it))
                           beta=rp11f(irlp(it))
                           sr=rp8f(irlp(it))
                           smax=rp9f(irlp(it))
                           call humidity(hum,alp,beta,sr,smax,to(i),
     &                          s(i))
                        else
c     use matrix values for neq+1 to 2*neq nodes
                           hum=eskc(i)
                           alp=rp3f(irlp(it))
                           beta=rp4f(irlp(it))
                           sr=rp1f(irlp(it))
                           smax=rp2f(irlp(it))
                           call humidity(hum,alp,beta,sr,smax,to(i),
     &                          s(i))
                        endif
                     endif
                  else if(icap(i).eq.1) then
c WARNING : not yet implemented for icap
                  endif
                  if(s(i).gt.0.0.and.s(i).lt.1.0) then
                    if(eskc(i).gt.0.0) then
                     pflow(i)=-s(i)
                     ka(i)=-1
                     eskc(i)=0.0
                     wellim(i) = sx1(i)*1.e-6
                    endif
                  else
                     eskc(i)=0.0
                     if (iout .ne. 0) write(iout,*) 
     &                    'humidity not set at node = ',i
                     if (iptty .ne. 0) write(iptty,*) 
     &                    'humidity not set at node = ',i
                  endif
  
               else if(eskc(i).lt.0.0) then
c specified saturation 
                if(eskc(i).eq. -888) then
                 sat_bc = 888
                     if (iout .ne. 0) write(iout,*) 
     &                    'humidity BC disabled at node = ',i
                     if (iptty .ne. 0) write(iptty,*) 
     &                    'humidity BC disabled at node = ',i
                 pflow(i)=-sat_bc
                 sat_bc = 0.0
                 go to 112
                endif
                sat_bc = abs(eskc(i))  
                if(sat_bc.gt.1.0) then
                     if (iout .ne. 0) write(iout,*) 
     &                    'sat bc = ',sat_bc,
     &                    'sat_bc set to 1 at node = ',i
                     if (iptty .ne. 0) write(iptty,*) 
     &                    'sat bc = ',sat_bc,
     &                    'sat_bc set to 1 at node = ',i
                  sat_bc = 1.0
                  pflow(i)=-sat_bc
                endif
                if(sat_bc.lt.tol_sat_bc) then
                     if (iout .ne. 0) write(iout,*) 
     &                    'sat bc = ',sat_bc,
     &                    'sat_bc set to 0 at node = ',i
                     if (iptty .ne. 0) write(iptty,*) 
     &                    'sat bc = ',sat_bc,
     &                    'sat_bc set to 0 at node = ',i
                  sat_bc = tol_sat_bc
                  pflow(i)=-sat_bc
                endif
112              continue
                 ka(i)=-1
                 eskc(i)=0.0
c gaz debug 121315
c                 pflow(i) = -sat_bc  (commented out 010416)
c                 s(i) = sat_bc
c                 so(i) = sat_bc
c gaz debug 081715
c                wellim(i) = sx1(i)*1.e06    
                 wellim(i) = sx1(i)*1.e3
              endif             
            enddo
        elseif(iflg.eq.-1) then
c gaz 122314 added another iflag (-1) to split funcionality
c called from startup 
            do mi=1,n
               denpch(mi)=denpci(mi)*dtot
               amc=amc+denpch(mi)*volume(mi)
               pcio(mi)=pci(mi)
               denpci(mi)=0.
            enddo
         elseif(iflg.eq.2) then
            do i=1,n
               pci(i)=pcio(i)
            enddo
         elseif(iflg.eq.3) then
            acner=0.0
            qtotci = qtotc
            do i=1,n
               qtc=qtc+qc(i)*dtotdm
               if(qc(i).gt.0.0) then
                  qtotc=qtotc+qc(i)*dtotdm
               endif
               dencht=denpch(i)
               denpch(i)=denpci(i)*dtot+dencht
               acner=acner+denpch(i)*volume(i)
               denpcj(i)=denpci(i)
               denpci(i)=0.0
               pcio(i)=pci(i)
            enddo
         elseif(iflg.eq.4) then
            if(iac.eq.0) then
               difc = 0.0
               qtcd = qtotc-qtc
               qcmax = max(abs(qtotc), abs(qtcd))
               if(abs(qtcd).gt.0.) difc=(acner-amc+qtc)/qcmax
            endif
         elseif(iflg.eq.5) then
            if(ntty.eq.2)  then
c     
c     organize for dual and dpdp outputs
c     
               if(idpdp.ne.0) then
                  mlev1=m/2+1
                  mlev2=0 
               else if(idualp.ne.0) then
                  mlev1=m/3+1
                  mlev2=2*m/3+1
               else          
                  mlev1=0
                  mlev2=0
               endif

               if (iout .ne. 0) write(iout,803)
               if(iatty.gt.0) write(iatty,803)               
               do i=1,m
                  if(i.eq.mlev1) then
                     if (iout .ne. 0) write(iout,*) 
     &                    ' matrix level = ', 2 
                     if (iatty .ne. 0) write(iatty,*) 
     &                    ' matrix level = ', 2 
                  endif
                  if(i.eq.mlev2) then
                     if (iout .ne. 0) write(iout,*) 
     &                    ' matrix level = ', 3 
                     if (iatty .ne. 0) write(iatty,*) 
     &                    ' matrix level = ', 3 
                  endif
                  mdd=nskw(i)
                  md=mdd
c calculate the relative humidity fraction
                  if(ieos(md).ne.3) then
                   hum_frac = 1.
                  else
c gaz  debug 062516                    
c calculate the relative humidity here                      
                   hum_frac = humida(md)
                  endif
c gaz 110119 calculate air mass fraction
                  dencht = denpch(md)/(denpch(md)+denh(md)+tol_p)
                  rqd= sk(md)
                  qcd=0.0
                  bpd=bp(md+nrhs(3))
c                 if(rqd.ne.0) qcd=qc(md)/rqd
c gaz 5-3-2001 need to output air source/sink
                  qcd=qc(md)
                  if (iout .ne. 0) write(iout,804)
     &                 mdd,pci(md),pcp(md),phi(md)-pcp(md),qcd,
     &                 bpd, ieos(md), hum_frac, dencht
                  if (iatty .ne. 0) write(iatty,804)
     &                 mdd,pci(md),pcp(md),phi(md)-pcp(md),qcd,
     &                 bpd, ieos(md), hum_frac, dencht
               enddo
 803           format(/, 20x, 'Nodal Information (Gas)', /, 8x, 
     &              'Partial P', 3x, 'Capillary', 3x, 'Liquid', 6x, 
     &              'Gas source/sink', /, 3x, 'Node', 1x, 'Gas (MPa)', 
     &              3x, 'Pres (MPa)', 2x, 'Pres (MPa)', 5x, '(kg/s)', 
     &              7x, 'Residual','    State',
     &              x,'R humidity',3x,'air m frac')
 804           format(i7, 1x, g11.4, 1x, g11.4, 1x, g11.4, 1x,
     &              g11.4, 5x, g11.4, 1x, i5, 3x, f8.4, 2x, g11.4)
c     calculate global mass and energy flows
               if (iout .ne. 0) then
                  write(iout,689) qtc,difc
                  write(iout,691)
                  write(iout,699) qtotc-qtotci
                  write(iout,690)
                  write(iout,699) qtotc
               end if
               if(iatty.gt.0) then
                  write(iatty,689) qtc,difc
                  write(iatty,691)
                  write(iatty,699) qtotc-qtotci
                  write(iatty,690)
                  write(iatty,699) qtotc
               endif

 689           format(/,
     &              'net gas discharge ',g15.5,' balance error gas ',
     &               g15.5)
 691           format(/,'this time step discharges : gas')
 699           format(1x,1pe14.3,' kg ',1pe14.3,
     &              ' mj ',1pe14.3,' mw ')
 690           format(/,'cumulative discharges : gas')
            endif
         endif
      endif
 9000 continue
      if( istflag .ne. 0 ) then
         stop
      end if
      
      return
      end
      

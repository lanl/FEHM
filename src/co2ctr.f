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
CD1 To provide overall control for a nonisothermal air-water(AWH) simulation.
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
      use commass_AWH, only : imass_phase, ivar_mass
      implicit none 

       real*8 bpd,tbnd,pcid,pv,dtsatp,dpsats,dpsatt,tdumm
       integer iflg,ico2d,i,istflag,mi,it
       real*8 hum,alp,beta,sr,smax,qtcd,dencht,tol_p
       real*8 hum_frac
       real(8) :: qtotci = 0.
       integer mlev1,mlev2,mdd,md
       real*8 rqd,qcd,psatl,qcmax       
       real*8 sat_bc, tol_sat_bc
c  gaz 041721     
       real*8 dum_hum, pv_hum
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

c tam  avoid memory error when junk is passed to parse_string2
c      initialize parse_string2 parameters
       nwds = 0
       imsg = 0
       xmsg = 0.
       cmsg = ''
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
            qtotin = 0.0
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
c gaz 052322 total mass             
       if(ivar_mass.eq.0) then      
c     
c     read in initial noncondensible gas pressure
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
               else if (pci(i).gt.-1.0d0.and.pci(i).lt.0.0d0) then
c gaz added phase designation for humidity IC (gas,ieos(i) = 3 )   
                  ieos(i) = 3
               end if
            end do

c     Check to verify if new input is being used for group 3
            read (inpt, '(a80)') dumstring
            backspace inpt
            nwds=0
            call parse_string2(dumstring,imsg,msg,xmsg,cmsg,nwds)
            if (nwds .gt. 0 .and. nwds .lt. 5) then
               old_input = .true.
               write (ierr, 14)
               if (iout .ne. 0) write (iout, 14)
               if (iptty .ne. 0) write(iptty, 14)
            end if
 14         format ('**** Warning: Old style ngas input being used, ',
     &          'macro input should be updated ****')
         else
c gaz 052322    
c     read in initial noncondensible gas mass fraction
c  
         if(.not.allocated(zntotf)) allocate(zntotf(n0))  
            narrays = 1
            itype(1) = 8
            default(1) = 0.
            macro = "ngas"
            igroup = 2
            
            call initdata2(inpt,ischk, n0, narrays, itype,
     &           default, macroread(2), macro, igroup, ireturn,
     &           r8_1= zntotf(1:n0))
                 
            
          endif  
            
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
               nwds=0
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
c     read in source strength for noncon gas(kg/s)
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
c
            allocate (xairfl(n0))   
            do i=1,n0
               if(qng(i).eq.default(1)) then
                  qng (i)=0.0
               elseif(qng(i).eq.0.0) then
                  qng(i)=1.d-30
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
c gaz 052522 added call for phase determination with ngas mass frac variable
c gaz 022423  changed to ivar_mass
          if(ivar_mass.eq.0) then
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
          else 
c gaz 052422 initial phase determination with mass variable 
c T is always fixed 

            call mass_component_fractions(0,0,0,0,0,0)     
            call mass_component_fractions(1,1,neq,0,0,0)  

          endif
              elseif(iflg.eq.1) then
c     
c     figure out humidity condition 
c gaz 041721 removed VG based humidity for true relative humidity                      
c     
            do i=1,n
               if(eskc(i).gt.0.0) then
c gaz 041721 specified humidity new formulation
c uses formulation for specified (relative) humidity in flow_boundary_conditions
c this formulation uses T,P that are the initial conditions                   
                if (.not. allocated (huma)) allocate(huma(n0))
                if (.not. allocated (xnva)) allocate(xnva(n0))
                if (.not. allocated (entha)) allocate(entha(n0))
                if (.not. allocated (phuma)) allocate(phuma(n0))
                if (.not. allocated (thuma)) allocate(thuma(n0))
                if (.not. allocated (pchuma)) allocate(pchuma(n0))
                     iha = 1
                     ieos(i) = 3
                     phuma(i) = pho(i)
                     thuma(i) = to(i)
                     huma(i) = -eskc(i)
                     pv_hum = psatl(thuma(i),dum_hum,dum_hum,dum_hum,
     &                   dum_hum,0,dum_hum)                     
                     pchuma(i) = phuma(i)-pv_hum*eskc(i)
                     t(i) = thuma(i) 
                     phi(i) = phuma(i)
                     pci(i) = pchuma(i)
                     if(wellim(i).eq.0.0) wellim(i)=sx1(i)*1.d06
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
c gaz 122314 added another iflag (-1) to split functionality
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
               else
                  qtotin = qtotin+qc(i)*dtotdm
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
c gaz 021121 removed phase state (now in water output)                  
                  qcd=qc(md)
                  if (iout .ne. 0) write(iout,804)
     &                 mdd,pci(md),pcp(md),phi(md)-pcp(md),qcd,
     &                 bpd, hum_frac, dencht
                  if (iatty .ne. 0) write(iatty,804)
     &                 mdd,pci(md),pcp(md),phi(md)-pcp(md),qcd,
     &                 bpd, hum_frac, dencht
               enddo
 803           format(/, 20x, 'Nodal Information (Gas)', /, 8x, 
     &              'Partial P', 3x, 'Capillary', 3x, 'Liquid', 6x, 
     &              'Gas source/sink', /, 3x, 'Node', 1x, 'Gas (MPa)', 
     &              3x, 'Pres (MPa)', 2x, 'Pres (MPa)', 5x, '(kg/s)', 
     &              7x, '   Residual',
     &              x,'R humidity',2x,'ngas m frac')
 804           format(i7, 1x, g11.4, 1x, g11.4, 1x, g11.4, 1x,
     &              g11.4, 5x, g11.4, 3x, f8.4, 2x, g11.4)
c   air  calculate global mass and energy flows
c gaz 040621 modified for clarity and added information 
               if (iout .ne. 0) then
                  write(iout,785) 
                  write(iout,786) acner
                  write(iout,689) abs(qtotin)
                  write(iout,690) qtotc                  
                  write(iout,691) qtotc-qtotci
                 write(iout,699) difc
               end if
               if(iatty.gt.0) then
                  write(iatty,785)
                  write(iatty,786) acner
                  write(iatty,689) abs(qtotin)                  
                  write(iatty,690) qtotc
                  write(iatty,691) qtotc-qtotci                  
                  write(iatty,699) difc
               endif
785           format(/,20x,'Global Gas (Air Only) Balances')
786           format('Total air mass at this time (kg):',t35,
     &           1pe15.6)
               
 689           format('Total gas inflow (kg): ',t35,1pe15.6) 
 691           format('Gas discharge this time step (kg):',
     &              t35,1pe15.6) 
 699           format('Balance error gas: ',t35,   e15.6)
 690           format('Total gas discharge (kg):',t35,1pe15.6)
            endif
         endif 
      endif
 9000 continue
      if( istflag .ne. 0 ) then
         stop
      end if
      
      return
      end
      
  
      subroutine mass_frac_total(iflg,ij)
c switch variables is necessary
c gaz 121620  initial coding
c gaz 051722 initial installation in co2ctr.f      
       use comfi
       use comai
c       use comdi, only : s
       use comdi
       use com_prop_data
       implicit none
       real*8 tl,pl,pcl,sl
       real*8 xtol,alpha, dalpca,tsolfac, alpha_tol
       real*8 dalphat,dalphap,dtsolfacdt,dtsolfacdp
       real*8 tboil,dtsatp,dpsats,dtsolfacdpcl,dalphapcl
       real*8 psatl, dpsatt, pcl_partial, dpcl_partialt, dpcl_partialp
       real*8 z_total, z_massfrac,  xnv, rol, rov
       integer iflg,i,mi,ij,ieosd           
       parameter(xtol=1.d-16, alpha_tol = 1.d-9)   
       xnl_max =0.1             
        if(iflg.eq.0) then
c read input
c gaz 051722 starting with AWH
c read in mass fraction and convert to original variables 
         
       else if(iflg.eq.1) then  
c switch variables case 1  
         if(ieosd.eq.2) then
c possible change from 2 phase conditions  
c since calculating on an individual gridblock we can ignore porosity and volume

c xhua error#6631 A non-optional actual argument must be present when 
c                 invoking a procedure with an explicit interface
c        call mass_component_liq(ij,2)    this is old one
c        call mass_component_gas(ij,2)    this is old one
         call mass_component_liq(iflg,ij,2)
         call mass_component_gas(iflg,ij,2)

         sl = s(ij)
         sv =1.0-sl
         z_total = (xnl*rol*sl+xnv*rov*sv)
         z_massfrac = z_total/(rol*sl+rov*sv)
c check for liquid phase transition
         call mass_component_liq(iflg,ij,1)
c check for gas phase transition
         call mass_component_gas(iflg,ij,1)
c check for sc phase transition
         call mass_component_sc(iflg,ij,1)         
         else if(ieosd.eq.1) then         
c possible change from liquid conditions          
         else if(ieosd.eq.3) then
c possible change from gas conditions          
         else
c possible change from SC conditions          
         endif

       else if(iflg.eq.2) then
c gaz 041122 calculate mass of gridblock for each phase type
c called after variable NR update       
        do i = 1, neq
         pl = phi(i)
         tl = t(i)
         pcl = pci(i)
         
        enddo

       else
        
       endif
       
      return

      end
      subroutine mass_component_fractions(iflg,i1,i2,iphase,ndummy,
     &  ifail_phase_chk)
c gaz 051722 added       
c calculate non water liq component in binary mixture
      use comai
      use comci
      use comdi
      use comdti, only: n0
      use comfi
      use comii
      use com_prop_data
      implicit none
      integer iflg,mi,i,i1,i2,ieosd,maxit_phase,ndummy,ifail_phase_chk
      real*8 pcl,tl,pl,sl,air_mass0,air_mass1,pcl0
      real*8 water_mass1          
      real*8 rol,rov,ros,xnv,alpha
      real*8 pv, drocp
      real*8 sl_chng,pcl_in,residual,pcl1,roc1,tol_resid
      real*8 delpcl,dresid_dpcl,drovpc,drovt,pv1,rov1,pv_in
      real*8 pcl_orig,roc_orig,pv_h2o, rol_liq, rol_vap, por
      real*8 dresid_sl,delsl,var_dum,var_dum1,sl_best, sl_best32
      real*8 value(9), value_a(9)
      real*8 z_old, z, zl_ngas, roc_tmp, phi_tmp, dresid_masspc 
c gaz 072522 added local variables   
      real*8 pci_tmp, psatl, psatl_100, dpsatt_100, dpsats

      real*8 resid_mass, tol_mass_phase 
      integer  ij, iphase, istate, ieos_save, intv, intv_cnt
      integer inr0, max_inr0
      integer, allocatable :: node_intv(:)

c tam Error: Type mismatch between actual arguments  (CHARACTER(80)/REAL(8
c dumb, dumc changed from character*80 to real*8
      real*8 dum0, dum2, dum3, dumb, dumc
      character*80 dum1
      logical phase_nr(3), test_phase
      parameter (max_inr0 = 10, tol_mass_phase = 1.d-6, intv = 10)
c     parameter (test_phase = .true.)

cDEC$ FIXEDFORMLINESIZE:132      
      if(iflg.eq.0) then
c allocate memory            
c initialize and allocate memory   
       call fluid_props_control(0, 0, 0, fluid(1), 
     &       'all      ', '         ')
       call fluid_props_control(0, 0, 0, fluid(2), 
     &       'all      ', '         ')
      else if(iflg.eq.1.or.iflg.eq.-1) then
c calculate initial phase states   
c gaz 072522 check for gas phase conditions step 1
c check for liq phase conditions step 2
c if not then then 2-phase conditions
c if iflag = -1 then do not update variable (just check phase)
c
c
        test_phase = .false.
        if(l.ge.20) test_phase = .true.
c
        if(.not.allocated(ieos_prev)) then
         allocate(ieos_prev(n))
         allocate(node_intv(intv))
        endif
        if(.not.allocated(node_intv)) then
         allocate(node_intv(intv))
        endif
        ieos_prev = 0
        node_intv = 0
        n_phase_ch = 0
        do i = i1,i2
         mi = i + ndummy
         tl = t(mi)
         pl = phi(mi)
         phi_tmp = pl
         pci_tmp = pci(mi)
         zl_ngas = zntotf(mi)
         ieos_save = ieos(mi)
         phase_nr(1:3) = .false.
c check for gas only phase 100% humidity        
         ieos(mi) = 3  
c initial guess pci = 0.  
           if(t(mi).ge.tcrit_h2o) then
              psatl_100 = pcrit_h2o
           else
            psatl_100 = psatl(t(mi),pcp(mi),dpcef(mi),dpsatt_100,dpsats,
     &                   0,an(mi))
           endif
         pv = psatl_100
         pci(mi) = max(pl-pv,0.0d0)
         pl_last = pl
c removed fluid_prop_control (2h)             
c h2o_properties_new(iflg,iphase,var1,var2,var3,istate,var4,var5,var6)          
          call h2o_properties_new(4,3,pv,tl,dum1,istate,
     &                 dum2,value,dum3)
          rov	=  value(1) 

         do inr0 = 1, max_inr0
          pv = max(pl-pci(mi),0.0d0)
          pv = min(pv,psatl_100)
          call h2o_properties_new(4,3,pv,tl,dum1,istate,
     &                 dum2,value,dum3)
          rov	=  value(1) 
          drovpc = -value(3)
          call air_eos_sol_props(1,pci(mi),tl,istate,value_a)
          roc    = value_a(1)
          drocp  = 0.0
          drocpc = value_a(3)

          resid_mass = roc - zl_ngas*(roc + rov)
          dresid_masspc = drocpc - zl_ngas*(drocpc+drovpc)        
          pci(mi) = pci(mi) - resid_mass/dresid_masspc
          if(abs(resid_mass).le.tol_mass_phase) then   
           pl = pv+pci(mi) 
           if(iflg.eq.-1) pci(mi) = pci_tmp
           phase_nr(3) = .true.
           go to 100
          endif
         enddo
100      continue
c gaz pl check for gas phase     
         if(phi(mi).lt.pl.and.phase_nr(3)) then
c pci(mi) is the new partial pressure
          ieos(mi) = 3
          s(mi) = 0.0d0
          if(.not.test_phase) go to 102
         else
c reset to original partial pressure and phase state   
          ieos(mi) = ieos_save
          pci(mi) = pci_tmp 
         endif
c gaz check for liq only phase 
c solubility is henry's law (temperature dependent)type 
c pv is constant         
        ieos(mi) = 1
        pl = phi(mi)
        pci(mi) = max(pl-pv,0.0d0)
        do inr0 = 1, max_inr0
         call air_sol(tl,pl,pci(mi),xnl,dxnlp,dxnlpc,dxnlt)
         resid_mass = zl_ngas-xnl  
         dresid_masspc = -dxnlpc
         pci(mi) = pci(mi) - resid_mass/dresid_masspc
         if(abs(resid_mass).le.tol_mass_phase) then 
           pl = pv+pci(mi) 
           phase_nr(1) = .true.
           go to 101
         endif
        enddo
          ieos(mi) = ieos_save
          pci(mi) = pci_tmp 
          go to 104
101     continue       
         if(phi(mi).ge.pl.and.phase_nr(1)) then
          ieos(mi) = 1
          s(mi) = 1.0d0
          if(.not.test_phase) go to 102
         else
c reset to original partial pressure and phase state                
          ieos(mi) = ieos_save
          pci(mi) = pci_tmp 
         endif
c gaz at this point the phase state is 2
c still need to calculate pci(mi)
c water vapor in 2-phase is 100% humidity (pv = psatl_100)
c resid_mass = (xnl*rol*sl + xnv*ros*(1-sl) = zl_gas*(rol*sl+ros*(1-sl))
c
104      continue
          ieos(mi) = 2
           if(t(mi).ge.tcrit_h2o) then
              psatl_100 = pcrit_h2o
           else
            psatl_100 = psatl(t(mi),pcp(mi),dpcef(mi),dpsatt_100,dpsats,
     &                   0,an(mi))
           endif 
        pv = psatl_100
        pl = phi(mi)
        pci(mi) = max(pl-pv,0.0d0)
        sl = s(mi)
c gaz set s(mi) = 0.5 to get liq and vap props  
          s(mi) = 0.5
          call fluid_props_control(1, mi, mi,fluid(1), 
     &         'all      ', '         ') 
          call fluid_props_control(1, mi, mi,fluid(2),
     &         'all      ', '         ') 
          rol = den_h2o(mi,1)
          rov	= den_h2o(mi,4)          
          roc = den_ngas(mi,1) 
          xnl = xnl_ngas(mi,1)
          ros = roc+rov
          xnv = roc/(roc+rov) 
          s(mi) =sl 
        do inr0 = 1, max_inr0
c        
          resid_mass = xnl*rol*sl + xnv*ros*(1-sl) - 
     &                 zl_ngas*(rol*sl+ros*(1-sl)) 
          dresid_sl = xnl*rol - xnv*ros - zl_ngas*(rol-ros)
          sl = sl - resid_mass/dresid_sl
          if(abs(resid_mass).le.tol_mass_phase) then 
c don't update s() if initially ieos(id)=2           
           ieos(mi) = 2
           if(ieos_save.ne.2) then
            s(mi) = sl
           endif
           phase_nr(2) = .true.
           go to 102
          endif
          enddo   
c did not converge in max_inr iterations
        if(iptty.ne.0) write(iptty,201) mi,ieos(mi),phi(mi),t(mi),
     &    pci(mi),s(mi)  
        if(iout.ne.0) write(iout,201) mi,ieos(mi),phi(mi),t(mi),
     &    pci(mi),s(mi)  
        if(ierr.ne.0) write(ierr,201) mi,ieos(mi),phi(mi),t(mi),
     &    pci(mi),s(mi)    
201     format('>> did not converge in phase check',
     &   ' (sub mass_component_fractions()) <<',/,
     &   'node', t10,'phase state', t24,'pressure',t38,'temp',t52,
     &   'pc',t66,'saturation',/,
     &    i8,t10,i10,t20,g14.5,t34,g14.4,t48,g14.4,t62,g14.4)  
        if(iptty.ne.0) write(iptty,*) 'stopping or time step reduction'
        if(iout.ne.0) write(iout,*) 'stopping or time step reduction'
        if(ierr.ne.0) write(ierr,*) 'stopping or time step reduction'
        ifail_phase_chk = -1
        return
102     continue
        if(ieos_save.ne.ieos(mi)) then
          ieos_prev(mi) = ieos_save
          n_phase_ch = n_phase_ch +1
        endif
       enddo
        write(ierr,501) days,l,iad,n_phase_ch
        intv_cnt = 0
        do mi = i1,i2
        if(ieos_prev(mi).ne.0) then
         intv_cnt = intv_cnt + 1  
         if(mod(intv_cnt,intv).ne.0) then
          node_intv(intv_cnt) = mi
         elseif(mod(intv_cnt,intv).eq.0) then
          node_intv(intv_cnt) = mi     
          write(ierr,502)(node_intv(i),ieos_prev(node_intv(i)),
     &          ieos(node_intv(i)),i = 1,intv_cnt)
          intv_cnt = 0  
         endif
        endif
        enddo
        if(intv_cnt.ne.0) then
          write(ierr,502)(node_intv(i),ieos_prev(node_intv(i)),
     &          ieos(node_intv(i)),i = 1,intv_cnt)
          intv_cnt = 0 
        endif
501     format('time (days)',1x,g12.5,' tstep',1x,i6,
     &         ' iter',1x,i4,' phase chng',1x,i6)   
502     format(1x,10(1x,'node =',i7,',','(',i1,',',i1,')'))   
        if(iflg.eq.-1.and.phase_nr(1).eqv..false..and.
     &     phase_nr(2).eqv..false.
     &     .and.phase_nr(3).eqv..false.) then       
c something went wrong in. phase check    
         if(iptty.ne.0) write(iptty,*) 'phase check error - stopping'
         if(iout.ne.0) write(iout,*) 'phase check error - stopping'
         if(ierr.ne.0) write(ierr,*) 'phase check error - stopping'
         stop
        endif
      else if(iflg.eq.2) then
       do i = i1,i2
        ij = i +ndummy
       if(iphase.eq.1) then
c liquid phase solubility (lig mass frac xnl)
        tl = t(ij)
        pl = phi(ij)
        pcl = pci(ij)
        call air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)
       else if(iphase.eq.3) then
c gas mass fraction
        tl = t(ij)
        pl = phi(ij)
        pcl = pci(ij)
c air gas density
        call air_eos_sol_props(1,pcl,tl,istate,value_a)
        roc = value_a(1)
        droct = value_a(2)
        drocpc = value_a(3) 
c water vapor density       
        call h2o_properties_new(4,3,pl-pcl,tl,dum1,istate,
     &                 dumb,value,dumc)
        rov = value(1)
        xnv = roc/(roc+rov)
       else if(iphase.eq.2) then
c possible change from liq to two phase
        air_mass1=(denpci(mi)*dtot+denpch(mi))
        water_mass1 = (deni(ij)*dtot+denh(ij))
        z_old = air_mass1/(air_mass1+water_mass1)
c estimate saturation  
        call air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)
        call air_eos_sol_props(1,pcl,tl,istate,value_a)
        roc = value_a(1)
        droct = value_a(2)
        drocpc = value_a(3) 
        call water_eos_sol_props
     &        (1,ij,pl,pcl,tl,istate,value,value_a) 
        call air_sol(tl,pl,pcl,xnl,dxnlp,dxnlpc,dxnlt)  
        call h2o_properties_new(4,3,pl-pcl,tl,dum1,istate,
     &                 dumb,value,dumc)
        rov = value(1)
        xnv = roc/(roc+rov)
        sl = rov*(xnv-z)/(rov*(xnv-z)-rol*(xnl-z))
       endif
      enddo
      endif
      return 
      end
      subroutine mass_component_gas(iflg,ij,iphase)
c calculate non water gas component in binary mixture
      implicit none
      integer iflg, ij, iphase, inr, max_inr
      parameter (max_inr = 10)
      if(iflg.eq.1) then
c calculate initial phase state
       if(iphase.eq.3) then
        do inr = 1, max_inr
        enddo
       else
       endif
      else if(iflg.eq.2) then
      endif
      return 
      end 
      subroutine mass_component_liq(iflg,ij,iphase)
c calculate non water gas component in binary mixture
      integer iflg, ij, iphase
      real*8 xnl,rol
      
      return 
      end 
      subroutine mass_component_sc(iflg,ij,iphase)
c calculate non water gas component in binary mixture
      integer iflg, ij, iphase
      real*8 xnl,rol
      
      return 
      end       
      subroutine variable_switch_ngas(iflg)
c switch variables from(P,T,PC) to (P,T,XL or XV) (ieosd = 1 and 3)
c switch variables from(P,S,T) to (P,Z,T) (ieosd = 1 and 3)
      use comai
      use comci
      use combi
      use comdi
      use comei 
      use comfi
      use comgi
      use comii
      use davidi
      implicit none
      integer iflg,mi,i,ieosd,maxit_phase
      real*8 pcl,tl,pl,sl,air_mass0,air_mass1,pcl0,roc,drocpc
      real*8 water_mass1          
      real*8 droct,rol,rov,ros,xnl,xnv,alpha
      real*8 dxnlp,dxnlpc,dxnlt
      real*8 sl_chng,pcl_in,residual,pcl1,roc1,tol_resid
      real*8 delpcl,dresid_dpcl,drovpc,drovt,pv1,rov1,pv_in
      real*8 pcl_orig,roc_orig,pv_h2o, rol_liq, rol_vap, por
      real*8 dresid_sl,delsl,var_dum,var_dum1,sl_best, sl_best32
      real*8 value_a(9) 
c tam changed from 3 to 4 to match deriv in massfrac_derivatives()
      real*8 deriv(4,4)
      integer i1,i2,j,ij,istate,iphase,jmia,neqp1
      integer i3,kb,nr1,nr2,nr3
      real*8 ztol
      real*8 a11_base, a12_base, a13_base, a21_base, a22_base, a23_base
      real*8 a31_base, a32_base, a33_base
      real*8 a11_mod, a12_mod, a13_mod, a21_mod, a22_mod, a23_mod
      real*8 a31_mod, a32_mod, a33_mod
      parameter(ztol = 1.d-14)
      if(ivar_switch.le.0) return
      
      if(iflg.eq.0) then
c initialize z (xntot) here (after first call to thrmwc.f)   
       if(.not.allocated(zntotf)) then
       endif
       do i = 1, n
         if(denh(i).le.ztol) then
          zntotf(i) = 0.0
         else  
          zntotf(i) = dench(i)/(denh(i)+ztol)
         endif
       enddo
      else if(iflg.eq.1) then
c
c modify linear system variables from (P,T,PC) or (P,S,T) to (P,T,Z) or (P,Z,T) 
c called just after assembling NR equations
c      
      neqp1 = neq + 1

       do i = 1, n
        i1 = nelm(i)+1
        i2 = nelm(i+1)
        do j = i1,i2
        kb = nelm(j)
         jmia = j-neqp1
c         
         a11_mod = 0.0         
         a12_mod = 0.0     
         a13_mod = 0.0 
         a21_mod = 0.0     
	   a22_mod = 0.0     
         a23_mod = 0.0     
         a31_mod = 0.0     
	   a32_mod = 0.0     
         a33_mod = 0.0     
         
c chain rule
         if(ieos(kb).eq.1) then
         
c        deriv(1,3) = dpcdz      
         
          call massfrac_derivatives(1,ieos(kb),kb,deriv)
          a13_base =  a(jmia+nmat(3))
          a23_base =  a(jmia+nmat(6))
          a33_base =  a(jmia+nmat(9))	
          a13_mod  =  a13_base*deriv(1,3)  
          a23_mod  =  a23_base*deriv(1,3) 
          a33_mod  =  a33_base*deriv(1,3) 
          a(jmia+nmat(3))=a(jmia+nmat(3))+a13_mod
          a(jmia+nmat(6))=a(jmia+nmat(6))+a23_mod
          a(jmia+nmat(9))=a(jmia+nmat(9))+a33_mod
          
	 else if(ieos(kb).eq.2) then
	 
c	 deriv(2,1) = dsldp 
c        deriv(2,2) = dsldz
c        deriv(2,13) = dsldt

	  call massfrac_derivatives(1,ieos(kb),kb,deriv)

          a11_base =  a(jmia+nmat(1))
          a21_base =  a(jmia+nmat(4))
          a31_base =  a(jmia+nmat(7))
          a12_base =  a(jmia+nmat(2))
          a22_base =  a(jmia+nmat(5))
          a32_base =  a(jmia+nmat(8))
          a13_base =  a(jmia+nmat(3))
          a23_base =  a(jmia+nmat(6))
          a33_base =  a(jmia+nmat(9))
          a11_mod  =  a11_base*deriv(2,1)  
	  a21_mod  =  a21_base*deriv(2,1) 
          a31_mod  =  a31_base*deriv(2,1)          
          a12_mod  =  a12_base*deriv(2,2)  
	  a22_mod  =  a22_base*deriv(2,2) 
          a32_mod  =  a32_base*deriv(2,2) 
          a13_mod  =  a13_base*deriv(2,3)  
          a23_mod  =  a23_base*deriv(2,3) 
          a33_mod  =  a33_base*deriv(2,3) 
          a(jmia+nmat(1))=a(jmia+nmat(1))+a11_mod
	  a(jmia+nmat(2))= a12_mod
          a(jmia+nmat(3))=a(jmia+nmat(3))+a13_mod
          a(jmia+nmat(4))=a(jmia+nmat(4))+a21_mod
	  a(jmia+nmat(5))= a22_mod
          a(jmia+nmat(6))=a(jmia+nmat(6))+a23_mod 
          a(jmia+nmat(7))=a(jmia+nmat(7))+a31_mod
	  a(jmia+nmat(8))= a32_mod
          a(jmia+nmat(9))=a(jmia+nmat(9))+a33_mod
          
	 else if(ieos(kb).eq.3) then  
	 
c        deriv(3,3) = dpcdz  

	  call massfrac_derivatives(1,ieos(kb),kb,deriv)

          a13_base =  a(jmia+nmat(3))
          a23_base =  a(jmia+nmat(6))
          a33_base =  a(jmia+nmat(9))	
          a13_mod  =  a13_base*deriv(3,3)  
          a23_mod  =  a23_base*deriv(3,3) 
          a33_mod  =  a33_base*deriv(3,3) 
          a(jmia+nmat(3))=a(jmia+nmat(3))+a13_mod
          a(jmia+nmat(6))=a(jmia+nmat(6))+a23_mod
          a(jmia+nmat(9))=a(jmia+nmat(9))+a33_mod	  
	 endif

        enddo
       enddo
      else if(iflg.eq.2) then
c NR updates for  (P,T,Z) or (P,Z,T)

c
c     NR corrections for water and noncondensible
c
c     strd is passed through common
               nr1=nrhs(1)
               nr2=nrhs(2)
               nr3=nrhs(3)
c               
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  i3=i+nr3
                  ieosd=ieos(i)
                  if(ps(i).eq.0.0) then
                     t(i)=t(i)-bp(i2)*strd                    
                  elseif(ieosd.eq.1.or.ieosd.eq.4) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                     zntotf(i)=zntotf(i)-bp(i3)*strd
c calculate pci here  

                  elseif(ieosd.eq.2) then
                     phi(i)=phi(i)-bp(i1)*strd
                     zntotf(i)=zntotf(i)-bp(i2)*strd         
                     t(i)=t(i)-bp(i3)*strd                    
c calculate pci and s here 

                  elseif(ieosd.eq.3) then
                     phi(i)=phi(i)-bp(i1)*strd
                     t(i)=t(i)-bp(i2)*strd
                     zntotf(i)=zntotf(i)-bp(i3)*strd
c calculate pci and s here   

                  endif
c 
               enddo 

      else if(iflg.eq.3) then
      else if(iflg.eq.4) then
      endif
      return 
      end    

      subroutine massfrac_derivatives(iflg,ieosd,ij,deriv)
      
      use comci
      use comfi
      implicit none
      integer iflg, ieosd, ij
      real*8 deriv(4,4)
      real*8 xnv,xnl,rov,rol,dxnlp,dxnvp
      real*8 dsldz, denom, ddenomp, ddenomt 
      real*8 drolp, drovp, drolt, drovt, dsldp, dsldpc
      real*8 dsldt,dxnvt,dxnlt,z,sl
c zero out derivatives      
      deriv=0.0d0
c      
      if(iflg.eq.1) then
c change from PC to Z      
c     dpcdz = 1/dxnldpc
      deriv(1,3) = 1./dclcf(ij)
c      
      else if(iflg.eq.2) then
c change from S to Z
     
      rov = rovf(ij)
      rol = rolf(ij)
      xnl = cnlf(ij)
      xnv = cnvf(ij)
      z = zntotf(ij)
      sl = rov*(xnv-z)/(rov*(xnv-z)-rol*(xnl-z))
      dsldz = -rov/(rov*(xnv-z)-rol*(xnl-z)) -
     &    (rov*(xnv-z)/(rov*xnv-z)-rol*(xnl-z))*(-rov+rol)
      deriv(2,2) = dsldz
      denom = (rov*(xnv-z)-rol*(xnl-z))
      ddenomp = drovp*(xnv-z)-rov*dxnvp-drolp*(xnl-z)-rol*dxnlp
      dsldp = drovp*(xnv-z)+rov*dxnvp/denom - rov*(xnv-z)/denom**2*
     & ddenomp 
      deriv(2,1) = dsldp
      ddenomt = drovt*(xnv-z)-rov*dxnvt-drolt*(xnl-z)-rol*dxnlt     
      dsldt = drovt*(xnv-z)+rov*dxnvt/denom - rov*(xnv-z)/denom**2*
     & ddenomt 
      deriv(2,3) = dsldt
      dsldpc = 0.0
      
      else if(iflg.eq.3) then
c change from PC to Z  (XNV)    
c     dxnvp =   dcvf()
c     dxnvt =   dcvef()
c     dxnvpc =  dcvcf()
      deriv(3,3) =  dcvcf(ij)

      else
     
      endif
      return
      end
      
c *******************************************************************************      
      subroutine generate_balance_equations_ngas(iflg,i1,i2,ndummy)
c gaz 080622 equation assembly for mass-based phase change 
c form ratios for noncondensible mass fraction and energy fraction      
      use comai
      use combi,only : sx1,nelm,nelmdg
      use comci
c      use comdi,only : ps, denh, deneh, ieos, t, phi, s
      use comdi
      use comdti,only : n0
      use comei,only : a
      use comfi
      use comgi,only : bp
      use davidi,only : nmat,nrhs
      use com_exphase,only : i_ex_update, ieq_ex
      implicit none
      integer iflg,i,id,id1,i1,i2,nsizea1,nsizea,neqp1,ndummy   
      real*8 dtin, ztol, energy_norm, frac_ngas, sx1d
      real*8 delmax(3), tol_eq1, tol_eq2, tol_eq3
      real*8 fdum_phase_calc, fdum_phase_prev 
      real*8 tol_phase_calc, fdum_raw(3)
      real*8 dtot_orig, ts_fac
      real*8 strd_smpl
      real*8 phi_low,phi_high,t_low,t_high,pci_low,pci_high
      real*8 s_low,s_high
      logical test_091322
      integer ndex(3)
      integer inr, max_inr,ic_tol_count, inr_fail, ifail_phase_chk
      real*8 avg_inr
      integer i_ex_update_save, ieq_ex_save, ieos_save
      integer  ieos_in, num_phase
      real*8  weightavg, weight_21, dmpf_new, dmpf_old, dmpf_avg
      parameter (weight_21 = 1.0)
      parameter(ztol = 1.d-14, max_inr =10, tol_phase_calc = 1.d-4)
      parameter(ts_fac = 0.5d0)
c gaz 082622 debug line next
      id = iad+l+fdum+ieos(1)+t(1)+phi(1)+pci(1)+s(1)+so(1)
       test_091322 = .false.
       neqp1 = neq+1
       
       if(iflg.eq.0) then
c allocate memory  
        if(.not.allocated(zntotf)) allocate(zntotf(n0)) 
        if(.not.allocated(inode_chng)) allocate(inode_chng(n0)) 
        if(.not.allocated(mass_h2o)) then 
         allocate(mass_h2o(n0))
         allocate(mass_ngas(n0))
         allocate(energy_tot(n0))
         allocate(mass0_h2o(n0))
         allocate(mass0_ngas(n0))
         allocate(energy0_tot(n0))
         allocate(dum_calc_eos(n0))
         allocate(dtot_test(n0))
        endif
       else if(iflg.eq.-1) then
c deallocate memory
        deallocate(mass_h2o,mass_ngas,energy_tot)
        deallocate(mass0_h2o,mass0_ngas,energy0_tot)
       else if(iflg.eq.1) then
c reset Jacobian matrix
        nsizea1=nelm(neqp1)-neqp1
        nsizea=9*nsizea1
        do i=1,nsizea
         a(i)=0.0d00
        enddo
c gaz zero residuals
        bp = 0.0d0
       else if(iflg.eq.2) then
500     continue 
c reset Jacobian matrix
        nsizea1=nelm(neqp1)-neqp1
        nsizea=9*nsizea1
        do i=1,nsizea
         a(i)=0.0d00
        enddo        
        dtot_orig = dtot
        dtot_min = dtot
        call thrmwc(0)
c gaz 032523 moved ic_tol_count lower
c        ic_tol_count = 0   
        tol_eq1 = tol_phase_calc
        tol_eq2 = tol_phase_calc
        tol_eq3 = tol_phase_calc
c i1 = 1, i2 =neq (from varchk_AWH header)  
c gaz 032523 moved inr loop   
        do inr = 1, max_inr
c sub iteration loop
c     zero out arrays
        do i=1,neq
         bp(i+nrhs(1))=0.0d00
         bp(i+nrhs(2))=0.0d00
         bp(i+nrhs(3))=0.0d00
        enddo
        do id=i1,i2
c loop on all equations
        if(iad.eq.0) then    
         dtot_test(id) = dtot_orig
c save variables in case of timestep change
         phi_reset = phi(id)
         t_reset = t(id)
         pci_reset = pci(id)
         s_reset = s(id)
         ieos_reset = ieos(id)  
        endif
c     
c     form balance equations after NR update
c 
         dum_calc_eos(id)(1:4) = '0000'
         if(ps(id).gt.0.0) then
            call geneqc(id)
         else
            call geneqc(id)
            a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
            bp(id+nrhs(1))=0.0d00
            a(nelmdg(id)-neqp1+nmat(9))=sx1(id)
            bp(id+nrhs(3))=0.0d00
         endif
        enddo
c gaz 090222  dum_calc_eos(id) identifies non converged node 
        ic_tol_count = 0  
        do id =i1,i2
         fdum_raw(1) = abs(bp(id+nrhs(1)))
         if(fdum_raw(1).ge.tol_eq1) then
             dum_calc_eos(id)(1:4) = '1000'
             ic_tol_count = ic_tol_count + 1
         endif
         fdum_raw(2) = abs(bp(id+nrhs(2)))
         if(fdum_raw(2).ge.tol_eq1) then
             dum_calc_eos(id)(2:2) = '1'
             if(dum_calc_eos(id)(1:1).eq.'0') then
              ic_tol_count = ic_tol_count + 1
             endif
         endif
         fdum_raw(3) = abs(bp(id+nrhs(3)))
         if(fdum_raw(3).ge.tol_eq1) then
             dum_calc_eos(id)(3:3) = '1'
              if(dum_calc_eos(id)(1:2).eq.'00') then
              ic_tol_count = ic_tol_count + 1
             endif            
         endif
         continue
        enddo
c gaz 032923 if residuals are small, return
      if(ic_tol_count.le.0) then
       return
      endif
c solve diagonal equations
c need explicit update vvariables so thrmwc can be called with one node 
c not   doing explicit update, borrowing variables  
       write(ierr,200) l,days,dtot/86400.,iad,ic_tol_count
200    format(1x,'t step ',i4,1p,' days ',g13.5,' t step size',g13.6,       
     &   ' iad ',i3,' nonconverged ',i5)    
       i_ex_update_save = i_ex_update
       i_ex_update = 1
       ieq_ex_save = ieq_ex
        fdum_phase_calc = 0.0d0    
        weightavg = 1.0 
        idelp = 0
       inr_fail = 0
       inr_count =0
       node_awh = 0
       num_phase = 0
       inode_chng(1:neq) = 0
       do id1 = i1,i2
        id=id1+ndummy
c only loop on nonconverged nodes  
c        if(dum_calc_eos(id)(1:4).ne.'0000') then
count nonconverged nodes  (should total ic_tol_count)
        node_awh = node_awh +1
        sx1d = sx1(id)  
        strd_smpl =1.d0
        inode_chng(id) = 0
c gaz 082722 use explicit update for single node call to thrmwc         
         ieq_ex = id        
c remove accumulation terms from residual EQs        
        bpt0(1) =  bp(id+nrhs(1)) - (sx1d*deni(id)+sk(id))
        bpt0(2) =  bp(id+nrhs(2)) - (sx1d*denei(id)+qh(id))
        bpt0(3) =  bp(id+nrhs(3)) - (sx1d*denpci(id)+qc(id))
        mass_h2o(id) = (deni(id)*dtot + denh(id))*sx1d
        mass_ngas(id) = (denpci(id)*dtot + denpch(id))*sx1d
        energy_tot(id) = (denei(id)*dtot + deneh(id))*sx1d
        frac_ngas = mass_ngas(id)/(mass_ngas(id)+mass_h2o(id)+ztol)
        zntotf(id) = frac_ngas
        ieos_save = ieos(id) 
100      continue  
         inr_count = inr_count+1
         ieq_ex = id 
         call thrmwc(0)
         if(ps(id).gt.0.0) then
            call geneqc(id)
         else
            call geneqc(id)
            a(nelmdg(id)-neqp1+nmat(1))=sx1(id)
            bp(id+nrhs(1))=0.0d00
            a(nelmdg(id)-neqp1+nmat(9))=sx1(id)
            bp(id+nrhs(3))=0.0d00
         endif
         bpt(1) = bpt0(1) + sx1d*deni(id)+sk(id)
         bpt(2) = bpt0(2) + sx1d*denei(id)+qh(id)
         bpt(3) = bpt0(3) + sx1d*denpci(id)+qc(id)  
         fdum_phase_calc = abs(bpt(1))+abs(bpt(2))+abs(bpt(3))

c form derivatives with explicit flow terms (might include src/sinks here)     
         a33t(1,1) = sx1d*dmpf(id) + dq(id)
         a33t(1,2) = sx1d*dmef(id) + dqt(id)
         a33t(1,3) = sx1d*dmc(id) + dqpc(id)
         a33t(2,1) = sx1d*depf(id) + dqh(id)
         a33t(2,2) = sx1d*deef(id) + deqh(id)
         a33t(2,3) = sx1d*dec(id) + dcqh(id)
         a33t(3,1) = sx1d*dcp(id) + dqc(id)
         a33t(3,2) = sx1d*dce(id) + deqc(id)
         a33t(3,3) = sx1d*dcc(id) + dcqc(id)
         call inverse_matrix(3,a33t,a33ti,ndex,0)
         del(1) = a33ti(1,1)*bpt(1)+a33ti(1,2)*bpt(2)+a33ti(1,3)*bpt(3)
         del(2) = a33ti(2,1)*bpt(1)+a33ti(2,2)*bpt(2)+a33ti(2,3)*bpt(3)
         del(3) = a33ti(3,1)*bpt(1)+a33ti(3,2)*bpt(2)+a33ti(3,3)*bpt(3)
c 040223 variable update now in varchk_simple_awh(2,...)
         call varchk_simple_awh(2,id,id,0,0.d0,0.d0,strd_smpl)   
         mass_h2o(id) = (deni(id)*dtot + denh(id))*sx1d
         mass_ngas(id) = (denpci(id)*dtot + denpch(id))*sx1d
         energy_tot(id) = (denei(id)*dtot + deneh(id))*sx1d
         frac_ngas = mass_ngas(id)/(mass_ngas(id)+mass_h2o(id)+ztol)
         zntotf(id) = frac_ngas
         ifail_phase_chk = 0
c        
         strd_smpl = 1.0d0
         ieos_in = ieos(id)
         call varchk_simple_awh(1,id,id,0,0.d0,0.d0,strd_smpl)
         if(inode_chng(id).ne.0) then
           num_phase = num_phase +1
         endif
c      endif
1000  continue
      enddo
c gaz 032923 if no phase changes return
       if(num_phase.le.0) then
        return
       endif
      enddo
      write(ierr,310) l,days,iad, ic_tol_count, inr_fail
250   format('inr iter failed: id ',i5,'  bpt(1-3) ',1p,3(1x,g12.5),
     &  ' phi ',g12.5,' s ',g12.5,' t ',g12.5)
260   format('inr iter failed: id ',i5,' bpt(1-3) ',1p,3(1x,g12.5),
     &  g12.5,' t ',g12.5,' pci ',g12.5)      
300   format(1x,'node ',i6,' phase nr ',i3,' ieos_sv ',i2,
     &  ' phi ',' ieos ',i3,' resid ',1p,g12.5)
310   format(1x,'t step ',i6,' days ',1p, g13.5,' iad ',i3,
     & ' total unconverged nodes ',i6,' total failed inr iter ',i6)

c need explicit update to original setting     
      i_ex_update = i_ex_update_save
      ieq_ex =  ieq_ex_save       
c form ratios for mass and energy variables 
c calculate ngas mass fraction
c  like:    dencht = denpch(md)/(denpch(md)+denh(md)+tol_p)      
c calculate accum terms and use total balance eqs
       avg_inr = float(inr_count/max(node_awh,1))
       do id1=i1,i2
c     
c form fractions
c     
       id=id1+ndummy
        sx1d = sx1(id)
        mass_h2o(id) = (deni(id)*dtot + denh(id))*sx1d
        mass_ngas(id) = (denpci(id)*dtot + denpch(id))*sx1d
        energy_tot(id) = (denei(id)*dtot + deneh(id))*sx1d
        frac_ngas = mass_ngas(id)/(mass_ngas(id)+mass_h2o(id)+ztol)
        energy_norm= energy_tot(id)/(mass_ngas(id)+mass_h2o(id)+ztol)
        mass0_h2o(id) = bp(id+nrhs(1))*dtot
        mass0_ngas(id) = bp(id+nrhs(3))*dtot
        energy0_tot(id) = bp(id+nrhs(2))*dtot
        zntotf(id) = frac_ngas
      enddo
c
c       call mass_component_fractions(1,1,neq,0,0,0)  
       continue
      else if(iflg.eq.3) then
      
      endif
      return 
      end  

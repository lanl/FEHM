      subroutine airctr(iflg,ndummy)
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
CD1 To manage the isothermal air-water calculations.
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
CD2 $Log:   /pvcs.config/fehm90/src/airctr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:06   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.21   Mon Mar 31 08:28:24 1997   gaz
CD2 new iteration parameters
CD2 
CD2    Rev 1.20   Wed Jun 12 16:44:12 1996   zvd
CD2 Added missing comma to format statement
CD2 
CD2    Rev 1.19   Wed Jun 12 15:20:58 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.18   Mon Jun 03 11:17:42 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.17   Fri May 31 10:33:08 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.16   Wed May 08 13:51:52 1996   hend
CD2 Rearranged and Added Output
CD2 
CD2    Rev 1.15   Fri Apr 26 14:45:10 1996   gaz
CD2 minor changes to phase change criteria
CD2
CD2    Rev 1.14   Wed Feb 14 10:19:00 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.13   Wed Feb 07 10:07:44 1996   gaz
CD2 step length changes and phase change stradegy
CD2 
CD2    Rev 1.12   Mon Jan 29 13:10:44 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.11   Mon Jan 29 10:02:32 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.10   11/29/95 14:02:18   gaz
CD2 format change
CD2 
CD2    Rev 1.9   11/15/95 10:13:32   gaz
CD2 complimentary changes from air_rdof.f
CD2 
CD2    Rev 1.8   08/18/95 09:50:04   llt
CD2 irlp and irlpt were already defined, removed for cray
CD2 
CD2    Rev 1.7   06/02/95 10:16:58   llt
CD2 removed upwinding commons plus gaz changes
CD2 
CD2    Rev 1.6   05/01/95 15:16:14   gaz
CD2  phase change control modified (gaz)
CD2 
CD2    Rev 1.4   03/24/95 13:54:56   gaz
CD2 gaz added phi(I) lt crl(4,1) to phase change criteria
CD2 
CD2    Rev 1.3   03/23/95 18:14:52   gaz
CD2 gaz now determine and save phase state
CD2 also enable ico2 to determine irdof behavior
CD2 
CD2    Rev 1.1   03/18/94 15:40:54   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:12   pvcs
CD2 original version in process of being certified
CD2
c 17-mar-94
c got rid of b array access(done in wrtout.f
c 2/9/95 gaz initialized liquid pressure
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier              Type     Use  Description
CD3
CD3 iflg                    int       I   Parameter used to control
CD3                                           the execution of the
CD3                                           routine
CD3 ndummy                  int       I   Parameter used to obtain
CD3                                           correct node number for
CD3                                           dual porosity calculations
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name           Description
CD3
CD3 file number inpt
CD3                Contains all input data from FEHMN macros in ASCII
CD3                form
CD3 file number iout
CD3                Contains output data
CD3 file number iatty
CD3                Contains output data
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
CD4 ipsx, n0, ico2, qtc, qtotc, amc, inpt, pmin, pmax, tmin, tmax,
CD4 neq, dtot, ps, ieos, iieos, iwelb, s, phi, t, crl, dmpf, rolf, 
CD4 dmef, dq, ieos, s, dil, thx, thy, thz, nrhs, bp, ntty, iout, idualp,
CD4 m, idpdp, nskw, sk, qh, pcp, dte, dife, qtote, qtotei, n, qc, to,
CD4 tini, so, b, strd, irlp, irlpt
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
CD5 ico2d        int         Flag denoting the type of two-phase
CD5                              simulation
CD5 tref         real*8      Reference temperature for isothermal air-
CD5                              water calculation
CD5 pref         real*8      Reference pressure for isothermal air-
CD5                              water calculation
CD5 ndum         int         Temporary storage for variable neq
CD5 pssv         real*8      porosity
CD5 ieoss        int         Equation of state parameter
CD5 iieoss       int         Equation of state parameter
CD5 iwelbs       int         Denotes if wellbore solution is enabled
CD5 ssv          real*8      Saturation
CD5 phisv        real*8      Pressure
CD5 dmpfd        real*8      Derivative of mass accumulation term
CD5                              with respect to pressure
CD5 dmefd        real*8      Derivative of mass accumulation term
CD5                              with respect to energy variable
CD5 dqd          real*8      Derivative of mass source term with
CD5                              respect to pressure
CD5 i            int         Do loop parameter
CD5 mid          int         Do loop parameter
CD5 mi           int         Current node number
CD5 nr1          int         Flag for right hand side of equation
CD5 nr2          int         Flag for right hand side of equation
CD5 i1           int         Index for unknown number
CD5 i2           int         Index for unknown number
CD5 ilev         int         Number of sets of nodes
CD5 mlev         int         Number of nodes at which information is
CD5                              written
CD5 il           int         Do loop index over all node levels
CD5 md           int         Current node number being written
CD5 rqd          real*8      Term used in output calculation
CD5 qcd          real*8      Mass balance term
CD5 irlpsv       int         Temporary storage for irlp value
CD5 irlptsv      int         Temporary storage for irlpt value
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
CD9 2.5.1 Implement time-step mechanism
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
CPS BEGIN airctr
CPS 
CPS IF there is air present
CPS 
CPS   IF this call is for initialization
CPS   
CPS     Read model flag
CPS     IF isothermal air-water simulation is not specified
CPS       Set flag parameter
CPS     ELSE isothermal air-water simulation is not specified
CPS       Set flag parameter
CPS       Read reference temperature and pressure
CPS       Set minimum temperature and pressure values
CPS       Initialize parameter values
CPS       thermw - calculate density, compressibility, and viscosity...
CPS       ... of water
CPS       Set parameters for reference values
CPS       FOR each node
CPS         Set low values for thermal conductivities
CPS       ENDFOR
CPS       
CPS     ENDIF
CPS   
CPS   ELSEIF this call is to determine the phase state
CPS   
CPS     FOR each node
CPS       Determine flag for phase state
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is to update the solution
CPS   
CPS     FOR each node
CPS       Compute new pressure, saturation
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is to compute air water thermodynamics
CPS   
CPS     thrair - compute isothermal air water thermodynamic parameters
CPS     
CPS   ELSEIF this call is for output
CPS   
CPS     IF air output is being written
CPS     
CPS       IF this is a dual porosity simulation
CPS         Set parameter values for writing
CPS       ELSEIF this is a dpdp simulation
CPS         Set parameter values for writing
CPS       ELSE it is a single porosity simulation
CPS         Set parameter values for writing
CPS       ENDIF
CPS       
CPS       FOR each set of nodes
CPS         Write header denoting which set of nodes
CPS         FOR each node being written
CPS           Determine node number
CPS           Determine sink term written during output
CPS           Write information for this node
CPS         ENDFOR
CPS       ENDFOR
CPS       Write overall mass balance information
CPS       IF output is also going to a second file
CPS         Write overall mass balance information
CPS       ENDIF
CPS       
CPS     ENDIF
CPS     
CPS   ELSEIF this call is for initialization
CPS     
CPS     FOR all nodes
CPS       Set temperatures, source flow rates, and EOS flag
CPS       IF this is a compressed liquid
CPS         Set saturation to unity
CPS       ENDIF
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is to zero out enthalpy in equations for...
CPS   ... isotherms problems
CPS   
CPS     FOR each node
CPS       Zero out residual terms for enthalpy equations
CPS     ENDFOR
CPS     
CPS   ENDIF ends all choices for control parameter
CPS   
CPS ENDIF there is air present
CPS 
CPS END airctr
CPS 
C**********************************************************************
c
c subroutine to manage isothermal air calculations
c
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comxi
      use comwt
      use davidi
      use comflow 
      use comsplitts 
      implicit none

      integer iflg,ndummy,ico2d,ndum,nndum,ieoss,iieoss,iwelbs,i,mid
      integer mi,i1,i2,ilev,mlev,il,md,irlpsv,irlptsv,icapsv,nr1
      integer nr2,irdofsv 
      integer  ii,ij,im,inode,iwm,ja,k,kb
      real*8 tref,pref,pssv,ssv,phisv,dmpfd,dmefd,dqd,rqd,qcd
      real*8 inflow_thstime,inen_thstime,denht,deneht
      real*8 dels,delp,dfdum11,dfdum12,dfdum21,dfdum22
      real*8 dfdum11i,dfdum12i,dfdum21i,dfdum22i,detdf
      real*8 fdum01,fdum02,sx1d,phidum,phi_dif,phi_1,phi_2
      real*8 hmax, hmin, hmid
      real*8 cden_correction
      character*80 form_string
      real*8 pref_1_2,pref_2_1,s_1_2
c      parameter (pchng = 0.005,schng = 0.005)
c      gaz pchng and schng in comai 103009
      save tref,pref
      if (jswitch.ne.0) strd_iter = strd_rich

C     ? ico2d not passed in, not initialized
      ico2d = 0
c     return if no air present 
      if(ico2.lt.0) then
c     iflg is used to tell if call is for initialization
         if(iflg.eq.0) then
            qtc=0.0
            qtotc=0.0
            amc=0.0
c     
c     read in reference pressure and temperature
c     
            read(inpt,*) tref,pref
c     
c     set max and min values of p and t
c     
            pmin(1)=-1.e15
            pmax(1)=1.e15
            tmin(1)=-1.e05
            tmax(1)=1.e05

c     Read in NAPL properties for SZ NAPL case

            if(ico2.eq.-3) then
               read(inpt,*) dennapl, viscnapl
            end if
         elseif(iflg.eq.-1) then

c     
c     calculate density,compressibility,and viscosity of water
c     
            ndum=neq
            nndum = n
            irdofsv=irdof
            irdof=0
            neq=1
            n = 1
            dtot=1.0
            pssv=ps(1)
            ps(1)=1.0
            ieoss=ieos(1)
            iwelbs=iwelb
            iwelb=0
            ssv=s(1)
            irlpsv=irlp(1)
            irlptsv=irlpt(1)
            icapsv = icapp
            icapp = 0
            s(1)=1.0
            iieoss=iieos(1)
            ieos(1)=1
            iieos(1)=1
            phisv=phi(1)
            phi(1)=pref
            t(1)=tref
            irlp(1)=1
            irlpt(1)=1
            dmpfd=dmpf(1)
            dmefd=dmef(1)
            dqd=dq(1)
            call thermw(0)
            irdof=irdofsv
            neq=ndum
            n = nndum
            ieos(1)=ieoss
            iieos(1)=iieoss
            phi(1)=phisv
            ps(1)=pssv
            s(1)=ssv
            iwelb=iwelbs
            irlp(1)=irlpsv
            irlpt(1)=irlptsv
            icapp = icapsv
c     reference density is in crl(1,1)
c     reference liquid viscosity is in crl(2,1)
c     reference compressibility is in crl(3,1)
c     reference pressure in crl(4,1)
c     reference air viscosity is in crl(5,1)
c     subtract initial density calculation (added back in thrair.f)
            if(cden) then
             crl(1,1)=rolf(1)  - cden_correction(1)
            else
             crl(1,1)=rolf(1)
            endif
            crl(2,1)=1.0/(dil(1)/rolf(1))
            crl(3,1)=dmpf(1)/rolf(1)
            crl(4,1)=pref
            crl(5,1)=182.e-07
            crl(6,1)=tref
            if(ico2d.eq.-1) then
               crl(7,1)=1.0
            else
               crl(7,1)=0.0
            endif

            dmpf(1)=dmpfd
            dmef(1)=dmefd
            dq(1)=dqd
c     
c     make thernmal conductivities small
c     
c     GAZ 10-23-98 eliminate thermal conductivities for isothermal
c     do i=1,n0
c     thx(i)=1.e-30
c     thy(i)=1.e-30
c     thz(i)=1.e-30
c     enddo
c     
c     initialize 2-phase regions             
c     
            do i=1,n0
               ieos(i)=2
            enddo
         elseif(iflg.eq.-2) then
c     
c     for rich eq make "a" smaller
c     eliminate a_vxy
c     
            if(jswitch.ne.0) then
c               if (idpdp .eq. 0 .and. allocated(a_vxy)) 
c     &              deallocate(a_vxy)
               if (allocated(a_vxy)) deallocate(a_vxy)
               ndum = nelm(neq+1) - (neq+1)
               if(allocated(a)) deallocate(a)
               if (idpdp .eq. 0) then
                  allocate(a(2*ndum))
               else
                  allocate(a(8*ndum))
               end if
            end if
         elseif(iflg.eq.1) then
c     
c     determine phase state
c     
            if(ifree.ne.0) then
               
               do  mid=1,neq
                  mi=mid+ndummy
                  ieos(mi) = 2
               enddo
               continue
            else if(irdof.ne.13) then   
               if (iad.eq.0) strd = 1.0 
               pref = crl(4,1) 
               pref_1_2 =  pref - pchng
               pref_2_1 =  pref + pchng  
               s_1_2 = 1.0 - schng
	       if (jswitch.ne.0) then
                  do  mid=1,neq
                     mi=mid+ndummy
c     
                     if(phi(mi).lt.pref_1_2.and.ieos(mi)
     &                 .eq.1.and.days.ge.time_ieos(mi)) then
                        strd = strd_iter
                        s(mi) = s_1_2
                        phi(mi) = pref
                        ieos(mi)=2 
                        time_ieos(mi) = days + time_ch
                     else if(s(mi).lt.1.0-tol_phase.and.ieos(mi) 
     &                   .eq.1.and.days.ge.time_ieos(mi)) then 
                        strd = strd_iter
                        s(mi) = s_1_2
                        ieos(mi)=2 
                        time_ieos(mi) = days + time_ch
                     else if(s(mi).gt.tol_phase.and.ieos(mi)
     &                .eq.3.and.days.ge.time_ieos(mi)) then 
c   changed eq.2 to eq.3 gaz 103009     
                        strd = strd_iter
                        ieos(mi)=2                      
                     else if(s(mi).gt.1.0+tol_phase.and.ieos(mi)
     &                   .eq.2.and.days.ge.time_ieos(mi)) then     
                        strd = strd_iter
                        phi(mi) = pref_2_1
                        s(mi) = 1.
                        ieos(mi)=1
                        time_ieos(mi) = days + time_ch
                     else if(s(mi).le.-tol_phase.and.ieos(mi)
     &                   .eq.2.and.days.ge.time_ieos(mi)) then       
c                        s(mi)=0.0
                        strd = strd_iter
c     ieos(mi)=3
                     endif
c     
                  enddo
               else
c     
c     determine phase state
c     
                  if(irdof.ne.13) then
c     pnx used to be set in the following loops -- currently removed
                     if (iad.eq.0) strd = 1.0        
                     do  mid=1,neq
                        mi=mid+ndummy
c     
                        if(s(mi).lt.1.0.and.ieos(mi).eq.1) then
                           if(iad.eq.1) strd = strd_iter
                           strd = strd_iter
                           ieos(mi)=2
c     call phase_timcrl(days,day,dtot,dtotdm)
c     
                        else if(s(mi).gt.0.0.and.ieos(mi).eq.3) then
                           if(iad.eq.1) strd = strd_iter
                           strd = strd_iter
                           ieos(mi)=2
c     call phase_timcrl(days,day,dtot,dtotdm)
c     
                        else if(s(mi).gt.1.0.and.ieos(mi).eq.2) then
                           if(iad.eq.1) strd = strd_iter
                           ieos(mi)=1
c     call phase_timcrl(days,day,dtot,dtotdm)
c     
                        else if(s(mi).le.0.0.and.ieos(mi).eq.2) then
                           if(iad.eq.1) strd = strd_iter
                           s(mi)=0.0
                           strd = strd_iter
                           ieos(mi)=3
c     call phase_timcrl(days,day,dtot,dtotdm)
c     
                        endif
                     enddo
                     do  mid=1,neq
                        mi=mid+ndummy
                        
                        if(phi(mi).lt.1.d-2) then
                           if(iad.eq.1) strd = strd_iter
                           phi(mi)=1.d-2
                           if (so(mi).ge.1.0) then
                              s(mi)= 0.999
                              phi(mi)=5.d-2
                              strd = strd_iter
                           endif
                        endif
c     
                        if(s(mi).lt.0.0) then
                           if(iad.eq.1) strd = strd_iter
                           s(mi)=0.0
c     strd = strd_iter
                        endif
c     
                        if(s(mi).gt.1.0) then
                           if(iad.eq.1) strd = strd_iter
c     s(mi)=1.0
c     strd = strd_iter
                        endif
c     
                     enddo
                  endif
               endif
            endif
         elseif(iflg.eq.2) then
            if(abs(irdof).ne.14) then
c     
c     update solution
c     
               nr1=nrhs(1)
               nr2=nrhs(2)
c     strd is passed through common
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  phi(i)=phi(i)-bp(i1)*strd
c     gaz 041805
c     if (irdof .ne. 13 .or. ifree .ne. 0) 
                  if (irdof .ne. 13 ) 
     &                 s(i)=s(i)-bp(i2)*strd
               enddo
c     
c     call cascade redistribution if requested
c     
               if(iflux_ts.ne.0) then
                  call cascade_sat(1)
               endif
c     
            else
c     
c     update solution(with exact mass balance)
c     
               nr1=nrhs(1)
               nr2=nrhs(2)
c     strd is passed through common
               do i=1,neq
                  i1=i+nr1
                  i2=i+nr2
                  phi(i)=phi(i)-bp(i1)*strd
               enddo
c     
c     call wellrate(0,0) 
c     call thrair(0)
               call interblock_iso(0) 
               
            endif
         elseif(iflg.eq.3) then
            if(ico2.eq.-3) then
c     call isothermal SZ NAPL thermodynamics
               call thrsznapl(ndummy)
            else
c     call isothermal air water thermodynamics
               if(ifree.ne.0) then
                  call wtsictr(2)
                  call thrair(ndummy)
c     call wtsictr(4)
               else
                  call thrair(ndummy)
               endif
            end if
         elseif(iflg.eq.4) then
c     call equation generation and load a array
c     new:zero a_axy (anisotropy requires accumulation)
            a_axy=0.0d00
            if(jswitch.eq.1) then
               call gensl2_switch
            else 
               call gensl2
            endif
         elseif(iflg.eq.-4) then
c     solve explicit equations
            call gensl2_explicit       
c     
         elseif(iflg.eq.5 .and. irdof .ne. 13) then
            if(ntty.eq.2) then
c     
c     output for air
c     
               if (iout .ne. 0) write(iout,803)
               if (iatty .ne. 0) write(iatty,803)
c     
c     organize differing amounts of output for dpdp and dual solutions
c     
c     
               if(idualp.ne.0) then
                  ilev=3
                  mlev=m/3
               else if(idpdp.ne.0) then
                  ilev=2
                  mlev=m/2
               else
                  ilev=1
                  mlev=m
               endif
               
               do il=1,ilev
                  if(il.ne.1) then
                     if (iout .ne. 0) write(iout,702) il
                     if (iatty .gt. 0) write(iatty,702) il
 702                 format(2x,'Matrix Level = ',i1)
                  endif
                  do i=1,mlev
                     md=  nskw(i+(il-1)*mlev)
                     rqd= sk(md)
                     qcd=qh(md)
                     if(ihead.ne.0.and.ifree.ne.0) then
c     wtsi solution
                        if(s(md).le.0.0) then
                           phidum = crl(4,1) - phi_inc
                        else if(s(md).lt.1.0) then
                           phi_1 = head12(md,1)
                           phi_2 = head12(md,2)
                           phi_dif = phi_2-phi_1
                           phidum = crl(4,1) + s(md)*(phi_dif)
                        else
                           phi_1 = head12(md,1)
                           phi_2 = head12(md,2)
                           phidum = phi(md) - phi_inc + 
     &                          0.5*(phi_2-phi_1)
                        endif
                        if (iout .ne. 0) write(iout,804) 
     &                       md,phidum,0.0,phidum,qcd
                        if (iatty .ne. 0) 
     &                       write(iout,804) md,phidum,0.0,phidum,qcd
                     else if (ihead.ne.0 .or. 
     &                       (irdof .eq. 13 .and. abs(ifree) .ne. 1)) 
     &                       then
c     head solution
                        if (iout .ne. 0) write(iout,804) 
     &                       md,max(phi(md)-phi_inc,0.1d00),
     &                       0.d00,max(phi(md)-phi_inc,0.1d00),qcd
                        if (iatty .ne. 0) write(iatty,804)
     &                       md,max(phi(md)-phi_inc,0.1d00),
     &                       0.d00,max(phi(md)-phi_inc,0.1d00),qcd
                     else
c     two-phase solution
                        if (iout .ne. 0) write(iout,804) 
     &                       md,phi(md),pcp(md),
     &                       phi(md)-pcp(md),qcd
                        if (iatty .ne. 0)
     &                       write(iatty,804) md,phi(md),pcp(md),
     &                       phi(md)-pcp(md),qcd
                     endif
                  enddo
               enddo
 803           format(/,20x,'Nodal Information (Vapor)',/, 9x, 
     &              'Air(Vapor)',26x, 'source/sink', /, 3x,'Node',1x,
     *              '  P (MPa)',3x,'Cap P (MPa)',1x,'Liq P (MPa)'
     *              ,1x,'Air(vp) (kg/s)')
 804           format(i7,1x,g11.4,1x,g11.4,1x,g11.4,1x,g11.4)
c     calculate global mass and energy flows
               if (iout .ne. 0) then
                  write(iout,703) 
                  write(iout,704) qtotei
                  write(iout,705) qtote
                  write(iout,706) qte
                  write(iout,707) dife
               end if
c     calculate global mass and energy flows
               if(iatty.ne.0) then
                  write(iatty,703) 
                  write(iatty,704) qtotei
                  write(iatty,705) qtote
                  write(iatty,706) qte
                  write(iatty,707) dife
               endif
 703           format(/,20x,'Global Mass Balances (Vapor)')
 704           format(1x,'Vapor discharge this time step: ',e14.6,' kg')
 705           format(1x,'Total vapor discharge: ',9x,e14.6,' kg')
 706           format(/,1x,'Net kg vapor discharge (total out-total ',
     &              'in): ',e14.6)
 707           format(1x,'Conservation Error: ',25x,e14.6)
            endif
         elseif(iflg.eq.6) then
c     
c     store sk in qc
c     initialize t,to,tini,iieos
c     
            tref=crl(6,1)
            do i=1,n
               qc(i)=sk(i)
               if(to(i).eq.0.0d00) then
                  t(i)=tref
                  to(i)=tref
                  tini(i)=tref
               endif
               iieos(i)=1
               if(ieos(i).eq.1) then
                  ieos(i)=2
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     s(i)=1.0
                     so(i)=1.0
                  end if
               endif
               if(abs(irdof).eq.14) then
                  denj(i)=1.0-s(i)
               endif
            enddo
c     
c     check also for free surface calcs
c     
            if(ifree.ne.0) then
               call wtsictr(1)	
c               call wtsictr(12)
c     call wtsictr(9)
            endif
         elseif(iflg.eq.7) then
c     
c     this only applies if we have a moving water table, etc.
c     
            call headctr(3,0,0.0,0.0)

         elseif(iflg.eq.8) then
c     
c     convert from pressure to head
c     
            call headctr(2,0,0.0,0.0)
         elseif(iflg.eq.9) then
c     
c     convert from head to pressure boundary value
c     
            call headctr(6,0,0.0,0.0)

         elseif(iflg.eq.10) then
c     
c     convert from head to pressure boundary and initial 
c     conditions (based on pair = pref at max height)
c     correct for negative pressures
c     
            call head_2phase(0)           
            
         elseif(iflg.eq.11) then
c     
c     calculate the gridblock length in the gravity direction
c     
            if(.not.allocated(dzrg))then
               allocate(dzrg(neq))
            else
               deallocate(dzrg)
               allocate(dzrg(neq))               
            endif 
            do i = 1,neq_primary
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
c     distinguish between block and edge centered
               if(ivf.eq.-1) then
                  dzrg(i) = max(hmax,abs(hmin))
               else
                  dzrg(i) = abs(hmax-hmin)/2.	            
               endif      
            enddo
            
         elseif(iflg.eq.12) then
c     write wt outpt
c     do i=1,n_wt_cols
            iwm=0
            do ij=1,m
               i=wcol(nskw(ij)) 
               do k=1,iwm
                  if(i.eq.col_out(k)) goto 566
               end do
               iwm=iwm+1
               col_out(iwm)=i
               do im=n_col(i),1,-1
                  inode=col(i,im)
                  if(s(inode).lt.1.or.im.eq.1) then
                     wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                    cord(inode,igrav)
                     if(wt_elev.eq.0..and.im.ne.n_col(i)) then
                        inode=col(i,im+1)
                        wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                       cord(inode,igrav)
                     endif
                     if(wt_elev.eq.0..and.im.eq.n_col(i)) then
                        if (iptty .ne. 0) write(iptty, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                        if (iout .ne. 0) write(iout, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                     endif
                     goto 4009
                  endif
               end do
 4009          continue
               if (form_flag .eq. 2) then
                  write(ishiswt,4004) days, cord(inode,1),
     &                 cord(inode,2), cord(inode,3), wt_elev,
     &                 ps(inode), inode, nskw(ij)
               else
                  write(ishiswt,4005) days, cord(inode,1),
     &                 cord(inode,2), cord(inode,3), wt_elev,
     &                 ps(inode), inode, nskw(ij)
               end if
 566        end do
 4004       format(1x,6(g16.9,', '),i8,', ',i8)
 4005       format(1x,6(g16.9,1x),2(i8,1x))
 4006       format(1x,'all nodes are dry ', i8, 1x, 2(g16.9, 1x))
         else if (iflg.eq.14) then
c     write wt output for contours
            if (altc(1:4) .eq. 'avsx') then
               write (form_string, 4015) ' : ', ' : '
               write (isconwt, 4019)
            else if (altc(1:3) .eq. 'sur') then
               write (form_string, 4015) ', ', ', '
               write (isconwt, 4020)
            else if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'tec') 
     &              then
               write (form_string, 4015) ' ', ' '
               if (altc(1:3) .eq. 'tec') then
                  write (isconwt, 4018) days
               else
                  write (isconwt, '("04 1 1 1 1")')
                  write (isconwt, '(a)') 'X coordinate (m), (m)'
                  write (isconwt, '(a)') 'Y coordinate (m), (m)'
                  write (isconwt, '(a)') 'Z coordinate (m), (m)'
                  write (isconwt, '(a)') 
     &                 'Water table elevation (m), (m)'
               end if
            end if
            do i=1,n_wt_cols
               do im=n_col(i),1,-1
                  inode=col(i,im)
                  if(s(inode).lt.1.or.im.eq.1) then
                     wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                    cord(inode,igrav)
                     if(wt_elev.eq.0..and.im.ne.n_col(i)) then
                        inode=col(i,im+1)
                        wt_elev = (s(inode) - 0.5)*dzrg(inode) +
     &                       cord(inode,igrav)
                     endif
                     if(wt_elev.eq.0..and.im.eq.n_col(i)) then
                        if (iptty .ne. 0) write(iptty, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                        if (iout .ne. 0) write(iout, 4006) inode,
     &                       cord(inode,1), cord(inode,2)
                     endif
                     goto 4010
                  endif
               end do
 4010          continue
               write(isconwt,form_string) cord(inode,1), cord(inode,2),
     &              cord(inode,3), izonef(inode), wt_elev, 0.0
            end do

 4015       format("(1x, 3(g16.9, '", a, "'), i4, 2('", a, "', g16.9))")
 4020       format(1x, 'X (m), Y (m), Z (m), Zone, WT elev (m), ', 
     &           'WT elev2 (m)')
 4019       format(1x, 'X (m) : Y (m) : Z (m) : Zone : WT elev (m) : ',
     &           'WT elev2 (m)')           
 4018       format('variables = "X (m)" "Y (m)" "Z (m)" " Zone" ', 
     &           '"WT elev (m) "', '"WT elev2 (m)"', / 
     &           'zone t = "Simulation time ', g16.9, ' days"') 
         endif
      endif       

      return
      end
      subroutine phase_timcrl(days,day,dtot,dtotdm)
c     
c     adjusts timestep when phase change occurs
c     
      implicit none
      real*8 days,dtot,day,dtotdm
c     
      return
      if(dtot.gt.dtotdm) then
         days=days-day
         dtot=dtotdm  
         day = dtot/86400.0d00
         days=days+day
      endif
      return
      end

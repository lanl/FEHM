      subroutine sther(iieosd)
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
CD1 To set thermodynamic parameters when simple thermodynamics are
CD1 invoked.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 ?            G. Zyvoloski   N/A     Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/sther.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Jun 12 15:21:24 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.5   Mon Jun 03 11:18:40 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.4   Fri May 31 10:46:34 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.3   Fri Feb 16 11:35:42 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.2   Mon Feb 05 11:32:10 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:12:16   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:20   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 iieosd       int     I       Flag denoting if simple
CD3                                  thermodynamics are to be used
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 inpt                   I    File contains input data
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 inpt, ipsat, itsat,
CD4 ew1, ew2, ew3, ew4, ew5, ew6, ew7, ew8, ew9, ew10, ew11,
CD4 ev1, ev2, ev3, ev4, ev5, ev6, ev7, ev8, ev9, ev10, ev11,
CD4 tsa0, tspa1, tspa2, tspa3, tsb0, tspb1, tspb2, tspb3,
CD4 psa0, psta1, psta2, psta3, psb0, pstb1, pstb2, pstb3,
CD4 crl, crv, cvl, cvv, cel, cev, pmin, pmax, tmin, tmax, n0, iieos
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 None
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
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 tsat         real*8      Input saturation temperature
CD5 psat         real*8      Input saturation pressure
CD5 i            int         Do loop parameter
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
CD9 2.4.3 Equation-of-state models
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN sther
CPS 
CPS   IF thermodynamic parameters are read in
CPS     
CPS     Read in flags, parameter values
CPS     
CPS     IF input dictates we limit solution to all gas or all liquid
CPS       Set parameter to keep solution vapor if called for
CPS       Set parameter to keep solution liquid if called for
CPS       Set other thermodynamic parameters
CPS     ENDIF
CPS     
CPS   EXITIF we are now done
CPS   
CPS   ENDIF
CPS   
CPS   IF the EOS flag is negative
CPS     Set to positive value
CPS   ELSE
CPS     FOR all terms in thermodynamic expression
CPS       Initialize terms to zero
CPS     ENDFOR
CPS     Compute new thermodynamic terms
CPS     Set large limits for coefficient sets
CPS   ENDIF
CPS   
CPS ENDIF
CPS 
CPS FOR all nodes
CPS   Set EOS flag to value indicating simplified thermodynamics
CPS ENDFOR
CPS 
CPS END sther
CPS
C**********************************************************************

      use comii
      use comdi
      use comdti
      use comai
      use comxi, only: nmfil
      use comtable
      implicit none
C*****
C*****AF 11/15/10
C*****
c      include 'comtable.h'                ! phs 4/23/99 
C*****

      real*8 tsat,psat,dtsatp
      real*8 visl0,visl1,t0,t1,coef0,coef1
      integer i, iieosd, open_file
      character cdum*4, table_name*80
C*****
C*****AF 11/15/10
C*****
      integer ios                         ! phs 4/23/99
c      integer ipsat, itsat
C*****
      tableFLAG = 0.
C*****      

      if ( iieosd .eq. 0 )  then
c     read input when appropriate
c     only one simple thermo allowed 
         ipsat = 0
         read (inpt  ,'(a80)') wdd1
         if(wdd1(1:3).eq.'ice') then
c   simple ice thermo
c   no other input required
          read(wdd1,*) cdum,IceEndTemp,LiqEndTemp
          iieosd = 3  
          itsat = 11 
          do i=1,n0
           iieos(i)=3
          enddo 
          allocate(ieos_aux(n0))          
         else if(wdd1(1:4).eq.'boil') then
c   simple boiling model
c   no other input required       
          read(wdd1,*) cdum,LiqEndTemp, VapEndTemp
          iieosd = 3  
          itsat = 12
          do i=1,n0
           iieos(i)=3
          enddo  
c  for now start with ieos_aux = 1
          allocate(ieos_aux(n0))
          ieos_aux = 1   
          else if(wdd1(1:5).eq.'table') then  
c gaz 060720 check  for file name here (new input format)
c read water EOS table name (iwater_table set to 1 after table read and 'id'ed in scanin)
               if(num_eos_table.gt.1) then
                  do i = 1, num_eos_table-1
                   read(inpt,*)
                  enddo
               endif
               if(iwater_table.eq.1) then
                table_name = nmfil(31)
                 if(iout.ne.0) then
                  write(iout,81) 'water table EOS :', table_name
                 endif
                 if(iptty.ne.0) then
                  write(iptty,81) 'water table EOS:', table_name
                 endif
81               format('>>> ',a17,1x,a50) 
               endif  
c gaz 070720 check for air table               
               if(iair_table.eq.1) then
c gaz 070720 read table from input file   
                table_name = nmfil(32)
                 if(iout.ne.0) then
                  write(iout,81) 'air table EOS :', table_name
                 endif
                 if(iptty.ne.0) then
                  write(iptty,81) 'air table EOS:', table_name
                 endif                      
              endif
c gaz 081921 check for co2 table               
               if(ico2wh_table.eq.1) then
c gaz 081921 read table from input file   
                table_name = nmfil(33)
                 if(iout.ne.0) then
                  write(iout,81) 'co2 table EOS :', table_name
                 endif
                 if(iptty.ne.0) then
                  write(iptty,81) 'co2 table EOS:', table_name
                 endif                      
              endif
c gaz 110715
c
c set very high bounds for table
c  
          do i=1,n0
           iieos(i)=1
          enddo   
         else
c classic simple thermo         
          read(wdd1,*) iieosd,itsat
          if(itsat.eq.2) then
           read(wdd1,*) iieosd,itsat,tsat,psat,dtsatp
          endif
          read (inpt  ,   *)  ew1,ew2,ew3,ew4,ew5,ew6,ew7,ew8,ew9, 
     &        ew10,ew11
          read (inpt  ,   *)  ev1,ev2,ev3,ev4,ev5,ev6,ev7,ev8,ev9,
     &        ev10,ev11
         endif
         if(itsat.le.10) then
         if ( itsat .eq. 2)  then
c   linear saturation line 
            ipsat = 0
            tsa0=tsat-dtsatp*psat
            tspa1=dtsatp
            tspa2=0.0
            tspa3=0.0
            tspa4 = 0.0
            tsb0=1.0
            tspb1=0.0
            tspb2=0.0
            tspb3=0.0
            tspb4=0.0
            psa0=psat-tsat/dtsatp
            psta1=1./dtsatp
            psta2=0.0
            psta3=0.0
            psta4=0.0
            psb0=1.0
            pstb1=0.0
            pstb2=0.0
            pstb3=0.0

         else if ( itsat .ne. 0 )  then
c     set ipsat=1 always when itsat ne 0
            ipsat=1
c     if vapor phase(itsat<0) set tsat=-1000.(this will keep it vapor)
c     if liquid phase(itsat>0) set tsat=+1000.(this will keep it liquid)
            if ( itsat .lt. 0) tsat=-1000.
            if ( itsat .gt. 0) tsat=1000.
            tsa0=tsat
            tspa1=0.0
            tspa2=0.0
            tspa3=0.0
            tsb0=1.0
            tspb1=0.0
            tspb2=0.0
            tspb3=0.0
            psat=1000.
            psa0=psat
            psta1=0.0
            psta2=0.0
            psta3=0.0
            psb0=1.0
            pstb1=0.0
            pstb2=0.0
            pstb3=0.0
         end if        
         if ( iieosd .eq. 0) goto 9000
         
      if ( iieosd .lt. 0 )  then
         iieosd=-iieosd

C*****
C*****AF 11/15/10
C*****
C     else
C     
c     end if      ! if ( iieosd .lt. 0 )
c-------------------------------------------
c-------CASE 5 is the LOOKUP TABLE
c-------Uses comtable.h to store the
c-------values read in from lookup.in   
c-------tableFLAG = 1 means table is in use!
c--------------------------------------------LOOKUP TABLE

      else if (iieosd.EQ.5) then 

         tableFLAG = 1
         iieosd = 1
         ios    = 0
         n      = 0
         pmax(1)   =   0.
         tmax(1)   =   0.
         pmin(1)   = 100.
         tmin(1)   = 100.
         incp   = 0.0
         inct   = 0.0
         nump   = 0
         numt   = 0

c--------read in the table                     LOOKUP TABLE

         iolookup = open_file(lookup_file, 'old')
C     
         do while (ios.EQ.0)
            n = n + 1
            read(iolookup,101,iostat=ios) (PP(n,i), i=1,11)
            if(ios.EQ.0) then
               pmin(1) = min(PP(n,1),pmin(1))
               tmin(1) = min(PP(n,2),tmin(1))
               pmax(1) = max(PP(n,1),pmax(1))
               tmax(1) = max(PP(n,2),tmax(1))
            end if
         end do
         n=n-1
         close(iolookup)
c---------LOOKUP TABLE
         do i = 1,n
            if (PP(i,1).NE.PP(i+1,1))then
               nump = nump + 1
               if(incp.EQ.0) then
                  incp = ABS(PP(i,1) - PP(i+1,1))
               end if
            end if
C     
            if (PP(i,2).NE.PP(i+1,2)) then
               if(inct.EQ.0) then
                  inct = ABS(PP(i,2) - PP(i+1,2))
               end if
            end if
         end do
C     
         numt = n/nump
c     LOOKUP TABLE
 101     FORMAT(f7.3,1x,f9.4,1x,9(e13.6,1x))
c-------------------------------------------------------
c--------------------------------------------------DO LINEAR 
C*****
C  
         
      else
c     linear expressions for thermo functions
c     exceptions are liquid viscosity and vapor density 
c     ew1-liquid reference pressure
c     ew2-liquid reference temperature
c     ew3-liquid reference density
c     ew4-liquid d(den)/dp
c     ew5-liquid d(den)/dt
c     ew6-liquid reference enthalpy
c     ew7-liquid dh/dp
c     ew8-liquid dh/dt
c     ew9-liquid viscosity
c     ew10-liquid dv/dp
c     ew11-liquid dv/dt
c     
c     ev1-vapor reference pressure
c     ev2-vapor reference temperature
c     ev3-vapor reference density
c     ev4-vapor d(den)/dp
c     ev5-vapor d(den)/dt
c     ev6-vapor reference enthalpy
c     ev7-vapor dh/dp
c     ev8-vapor dh/dt
c     ev9-vapor viscosity
c     ev10-vapor dv/dp
c     ev11-vapor dv/dt
c     
         do i=1,10
c     first zero density terms
            crl(i,iieosd)=0.0
            crl(i+10,iieosd)=0.0
            crv(i,iieosd)=0.0
            crv(i+10,iieosd)=0.0
c     first zero viscosity terms
            cvl(i,iieosd)=0.0
            cvl(i+10,iieosd)=0.0
            cvv(i,iieosd)=0.0
            cvv(i+10,iieosd)=0.0
c     first zero enthalpy terms
            cel(i,iieosd)=0.0
            cel(i+10,iieosd)=0.0
            cev(i,iieosd)=0.0
            cev(i+10,iieosd)=0.0
         enddo

c set density
         crl(11,iieosd)=1.0
         crl( 1,iieosd)=ew3-ew4*ew1-ew5*ew2
         crl( 2,iieosd)=ew4
         crl( 5,iieosd)=ew5
c     reflect the perfect gas law
         crv(11,iieosd)=273.0*ev1
         crv( 2,iieosd)=ev3*(273.0+ev2)
         crv(15,iieosd)=ev1

c     set enthalpy
         cel(11,iieosd)=1.0
         cel( 1,iieosd)=ew6-ew7*ew1-ew8*ew2
         cel( 2,iieosd)=ew7
         cel( 5,iieosd)=ew8
         cev(11,iieosd)=1.0
         cev( 1,iieosd)=ev6-ev7*ev1-ev8*ev2
         cev( 2,iieosd)=ev7
         cev( 5,iieosd)=ev8

c     set viscosity
c        cvl(11,iieosd)=1.0
c        cvl( 1,iieosd)=ew9-ew10*ew1-ew11*ew2
c        cvl( 2,iieosd)=ew10
c        cvl( 5,iieosd)=ew11
c        assume derivative wrt t 1 degC
c        assume assume no pressure dependence of liquid viscosity
c        t0=ew2
c        t1=t0+1.0
c        visl0=ew9
c        visl1=visl0+ew11
c        coef1=-ew11/(visl1*t1-visl0*t0)
c        coef0=visl0*(1.0+coef1*t0)
c        cvl(11,iieosd)=1.0
c        cvl(15,iieosd)=coef1                         
c        cvl( 1,iieosd)=coef0  
      if(ew11.ne.0) then           
         t0=ew2
         t1=t0+1.0
         visl0=ew9
         visl1=visl0+ew11
         coef1=-(visl1*t1-visl0*t0)/ew11
         coef0=visl0*(coef1+t0)
         cvl(11,iieosd)=coef1
         cvl(15,iieosd)=1.0                           
         cvl( 1,iieosd)=coef0   
      else
         cvl(11,iieosd) = 1.
         visl0=ew9
         cvl( 1,iieosd)=visl0
      endif                      
         cvv(11,iieosd)=1.0
         cvv( 1,iieosd)=ev9-ev10*ev1-ev11*ev2
         cvv( 2,iieosd)=ev10
         cvv( 5,iieosd)=ev11
    
c     set maximum limits large
         pmin(iieosd)=-2.
         pmax(iieosd)=1000.0
         tmin(iieosd)=-200.
         tmax(iieosd)=1500.0

C     
C*****
C*****AF
C     
      end if                    ! if iieosd .eq. 5
C     
C     
c     change coefficient set to iieosd
      end if   
      do i=1,n0
         iieos(i)=iieosd
      enddo
 9000 continue
c      endif
      else if(itsat.eq.11) then
c     set maximum limits to include ice
         pmin(3)=-100.
         pmax(3)=1000.0
         tmin(3)=-30.
         tmax(3)=150.0  
c set boiling temp to very high(1000 C)         
            tsa0=1000.
            tspa1=0.0
            tspa2=0.0
            tspa3=0.0
            tsb0=1.0
            tspb1=0.0
            tspb2=0.0
            tspb3=0.0
            psat=1000.
            psa0=psat
            psta1=0.0
            psta2=0.0
            psta3=0.0
            psb0=1.0
            pstb1=0.0
            pstb2=0.0
            pstb3=0.0       
      else if(itsat.eq.12) then
c     set maximum limits to include water-steam
         pmin(3)=0.01
         pmax(3)=1000.0
         tmin(3)=-3.
         tmax(3)=500.0  
c set boiling temp to very high(1000 C)              
            tsa0=1000.
            tspa1=0.0
            tspa2=0.0
            tspa3=0.0
            tsb0=1.0
            tspb1=0.0
            tspb2=0.0
            tspb3=0.0
            psat=1000.
            psa0=psat
            psta1=0.0
            psta2=0.0
            psta3=0.0
            psb0=1.0
            pstb1=0.0
            pstb2=0.0
            pstb3=0.0     
      endif
c      endif
      continue
      return
      end
      subroutine eos_aux(itype,t,p,iprop,ipres_on,prop,dpropt,dpropp)
c
c   subroutine manages  simple eos for ice and water-vapor
c     
      integer itype,iprop,ipres_on
      real*8 t,p,prop,dpropt,dpropp
      if(itype.eq.11) then
c  ice-water  
       call WaterProps(t,p,iprop,ipres_on,prop,dpropt,dpropp)
      else if(itype.eq.12) then
c  water-watervapor      
      call WaterProps_Vap(t,p,iprop,ipres_on,prop,dpropt,dpropp)
      endif
      end

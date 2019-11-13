      subroutine vcon (iz,ndummy)
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
CD1 To calculate variable thermal conductivity.
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
CD2 $Log:   /pvcs.config/fehm90/src/vcon.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jun 12 15:21:24 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.4   Mon Jun 03 11:18:44 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.3   Fri May 31 10:51:30 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.2   Mon Feb 05 11:33:54 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:12:30   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:14   pvcs
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
CD3 iz           int     I       Control parameter
CD3 ndummy       int     I       Parameter to yield correct node
CD3                                 number for dual porosity nodes
CD3                                  
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
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 ivcond, inpt, wdd1, vc1f, vc2f, vc3f, ivcon, n, ivcn, n0, neq,
CD4 thx, thy, thz, t, s, narrays, pointer, itype, default, igroup,
CD4 ireturn, macroread
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4 
CD4 Name       Type      Description
CD4 
CD4 initdata   N/A       Reads and sets input values defined at each
CD4                          node
CD4 null1       LOGICAL   Checks for null line in input file
CD4 
CD4 None
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
CD5 i            int         Counter for number of lines of data read
CD5 ivcnd        int         Parameter used for reading
CD5 mid          int         Do loop index parameter
CD5 mi           int         Current node number
CD5 it           int         Current group number
CD5 itp          int         Current conductivity model number
CD5 vc1          real*8      Current thermal conductivity parameter 1
CD5 vc2          real*8      Current thermal conductivity parameter 2
CD5 vc3          real*8      Current thermal conductivity parameter 3
CD5 vc4          real*8      Current thermal conductivity parameter 4
CD5 vc5          real*8      Current thermal conductivity parameter 5
CD5 vc6          real*8      Current thermal conductivity parameter 6
CD5 vc7          real*8      Current thermal conductivity parameter 7
CD5 vc8          real*8      Current thermal conductivity parameter 8
CD5 vc12         real*8      vc1 - vc2
CD5 sqrsat       real*8      Parameter used in calculation
CD5 tmpPor       real*8      Parameter used in calculation
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
CD9 2.4.11 Variable thermal conductivity
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
CPS BEGIN vcon
CPS 
CPS IF variable thermal conductivities are required
CPS 
CPS   IF this call is for reading input
CPS   
CPS     LOOP to read first group of data
CPS       null1 - determine if there is more data to read
CPS       IF more data was found
CPS         Read next line of input
CPS       ELSE
CPS         Set flag appropriately for case of no data read in
CPS     EXITIF no more data is found
CPS     ENDLOOP
CPS     
CPS     
CPS     Set parameters for reading region numbers
CPS       
CPS     initdata - read region numbers and set at each node
CPS       
CPS     Set flag to denote that the vcon macro has been called
CPS       
CPS   ELSE this call is for computing thermal conductivities
CPS   
CPS     FOR each node
CPS     
CPS       IF linear thermal conductivity model is chosen
CPS         Compute thermal conductivities
CPS       ELSEIF square root model is chosen
CPS         Compute thermal conductivities
CPS       ENDIF
CPS     
CPS     ENDFOR
CPS   
CPS   ENDIF
CPS 
CPS ENDIF variable thermal conductivities are required
CPS 
CPS END vcon
CPS
C**********************************************************************

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comdti
      use comai
      use comki
      implicit none
      
      integer iz,ndummy,i,ivcnd,mid,mi,it,itp
      real*8 vc1,vc2,vc3,vc4, vc5, vc6, vc7, vc8, vc12,sqrsat, tmpPor
      real*8 tmpS, tmpT, thm_sat, thm_dry, thm_solid
      real*8 thx_tol
      parameter (thx_tol=1.0d-12)
      logical null1

c gaz 110318 used Peter Johnson code completely
c read in data
      if(ivcond.ne.0) then
         if(iz.eq.0) then
            i=0
 10         continue
            read(inpt,'(a80)') wdd1
            if(.not. null1(wdd1)) then
               backspace inpt
               read(inpt,*) ivcnd
               backspace inpt
               if(isalt.ne.0) then
                if(ivcnd.ne.4.and.ivcnd.ne.5) then
                 write(ierr,*) 'warning non salt vcon model entered '
                 if(iout.ne.0) write(iout,*) 
     &                          'warning non salt vcon model entered '
                endif
               endif
               i=i+1
               if(ivcnd .eq.  1) then
                  read(inpt,*) ivcon(i),vc1f(i)
               elseif(ivcnd .eq. 2 .or. ivcnd .eq.  3) then
                  read(inpt,*) ivcon(i),vc1f(i),vc2f(i),vc3f(i)
               else if(ivcnd .eq. 4) then
c                 for crushed salt
                  read(inpt,*) ivcon(i),vc1f(i),vc2f(i),vc3f(i),
     &                         vc4f(i),vc5f(i),vc6f(i), vc7f(i),vc8f(i)
               else if(ivcnd .eq. 5) then
c                 for unsaturated case(Olivella,2011)
                  read(inpt,*) ivcon(i),vc1f(i),vc2f(i),vc3f(i) 
               endif
            else
               if(i.eq.0) ivcond=0
               go to 20
            endif
            go to 10
c     read in nodal capillary type
 20         continue

            narrays = 1
            itype(1) = 4
            default(1) = 0
            macro = "vcon"
            igroup = 2
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(9), macro, igroup, ireturn,
     2           i4_1=ivcn(1:n0) )
            
            macroread(9) = .TRUE.
            
         else

c     load nodal thermal conductivities
            do mid=1,neq
               mi=mid+ndummy
               it=ivcn(mi)

               if(it.ne.0) then
               	itp=ivcon(it)
               	if(itp.eq.1) then
c     constant thermal conductivity (similar to cond macro)
                thx(mi) = vc1f(it)*1e-6
                thy(mi) = thx(mi)
                thz(mi) = thx(mi)

               	elseif(itp.eq.2) then

c     linear variation with temperature
c     vc1f=reference temperature,vc2f=reference conductivity
c     vc3f=d(cond.)/d(temp) at reference conditions
c     NOTE:no derivatives,used explicity
                  vc1=vc1f(it)
                  vc2=vc2f(it)
                  vc3=vc3f(it)
                  thx(mi)=(vc3*(t(mi)-vc1) +vc2)*1.e-6
                  thy(mi)=thx(mi)
                  thz(mi)=thx(mi)
               	else if(itp.eq.2) then

c     square root variation with saturation
c     vc1f=conductivity at stauration=1
c     vc2f=conductivity at saturation=0
c     NOTE:no derivatives,used explicity
                  vc1=vc1f(it)
                  vc2=vc2f(it)
                  vc12=vc1-vc2
                  sqrsat=sqrt(s(mi))
                  thx(mi)=(vc12*sqrsat +vc2)*1.e-6
                  thy(mi)=thx(mi)
                  thz(mi)=thx(mi)

               	else if(itp.eq.3) then

c     Thermal conductivity for intact salt (Munson et al, 1990)
c     lambda_salt (T) = lambda_300*(300 K/T)^gamma
c     where lambda_300 is material constant = 5.4 (W/m K)
c           gamma is material constant = 1.14 (at 300K)
c           T is temperature in Kelvin
c     If the reference temperature is 300K, in input file, vc1f= 26.85, vc2f = 5.4, vc3f =1.14
C     some other values should be used if the reference temperute is different from 300K
c     NOTE:no derivatives,used explicity
                  vc1=vc1f(it)+273.15
                  vc2=vc2f(it)
                  vc3=vc3f(it)
                  thx(mi)=vc2*( vc1/(t(mi)+273.15) )**vc3*1.e-6
                  thy(mi)=thx(mi)
                  thz(mi)=thx(mi)

               	else if(itp.eq.4) then

c     Thermal conductivity for crushed salt (Bechthold et al, 2004)
c     lambda_salt (T) = lambda_300*(300 K/T)^gamma
c     where gamma is material constant = 1.14 (at 300K)
c           T is temperature in Kelvin
c           lambda_300 (W/m K), depends on porosity (read in as vc2f):
c              lambda_300 (phi) = 1.08*(-270*phi^4 +370*phi^3-136*phi^2+1.5*phi + 5.0)
c     If the reference temperature is 300K, in input file, vc1f = 26.85, vc3f =1.14
C     some other values should be used if the reference temperature is different from 300K
c     NOTE:no derivatives,used explicity

                  vc1=vc1f(it) + 273.15
                  vc2=vc2f(it) 
                  vc3=vc3f(it) 
                  vc4=vc4f(it) 
                  vc5=vc5f(it) 
                  vc6=vc6f(it) 
                  vc7=vc7f(it) 
                  tmpPor=ps(mi)
c				  DRH 1/3/2013: Relationship is only valid up to 0.4
c				  conductivity held constant at porosities above 0.4
c                  if(tmpPor.gt.0.4) tmpPor = 0.4
c pjjohnson fix to DRH fix
c relationship begins increasing above about por=0.3873
c may be inducing error in test problems where heat flow is too good at high porosity
c putting a linear extrapolation, thermal cond goes to 0.028 W/mK (value for air at room temp)
c
c SCHU revise 08282019 
c  fix pjjohnson change so thermal conductivity for porosity above 0,3873 calculated correctly
                  if(tmpPor.gt.0.3873) then
                     tmpPor = vc2*(0.960308 - 0.930308*tmpPor)
c                      thx(mi)=(0.97 - (0.95333*tmpPor))
                  else
                     tmpPor = vc2*(vc3*tmpPor**4 +vc4*tmpPor**3 
     &                 +vc5*tmpPor**2+vc6*tmpPor + vc7)
                  endif

                  vc8=vc8f(it)
                  thx(mi)=(1e-6*tmpPor)*(vc1/(t(mi)+273.15))**vc8
                  if(thx(mi).le.thx_tol) thx(mi) = thx_tol
                  thy(mi)=thx(mi)
                  thz(mi)=thx(mi)
c                 endif
                  else if(itp.eq.5) then
                      
c end pjjohnson                      

c     Thermal conductivity for unsaturated rock salt (Olivella, 2011)
c      modified by gaz 1/30/2016
c     Input parameters for this option are:
c           thermal conductivity of gas (or air), vc1f
c           thermal conductivity of liquid (or water), vc2f
c           thermal conductivity of solid (or grain), vc3f
c           Olivella (2011) values:
c           air (vc1f) - 0.1 W/m-K
c           water (vc2f) - 0.6 W/m-K
c           grain (vc3f) - 5.734 W/m-K
c           T is temperature in C (GAZ changed to C)
c     NOTE:no derivatives,used explicity
                  vc1=vc1f(it) 
                  vc2=vc2f(it) 
                  vc3=vc3f(it) 
                  tmpPor=ps(mi)
                  tmpS=s(mi)
                  tmpT=t(mi)   ! tmpT in C

                  thm_solid = vc3 -1.83d-2*tmpT + 2.86d-5*tmpT**2 
     &                        - 1.51d-8*tmpT**3
                  thm_dry = thm_solid**(1.0-tmpPor) * vc1**tmpPor
                  thm_sat = thm_solid**(1.0-tmpPor) * vc2**tmpPor
                  thx(mi) =  thm_dry**(1.-tmpS) * thm_sat**tmpS * 1.d-6

                  if(thx(mi).le.thx_tol) thx(mi) = thx_tol
                  thy(mi)=thx(mi)
                  thz(mi)=thx(mi)

               	endif
               endif
            enddo

         end if
      end if
     

      return
      end




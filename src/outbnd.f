      subroutine outbnd
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
CD1 To check if thermodynamic variables are out of bounds.
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
CD2 $Log:   /pvcs.config/fehm90/src/outbnd.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Fri Feb 16 10:44:50 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.5   Thu Feb 01 14:04:38 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   11/27/95 17:01:04   gaz
CD2 format change
CD2 
CD2    Rev 1.3   11/27/95 15:49:18   gaz
CD2 added coordinates of out of bound node
CD2 
CD2    Rev 1.2   09/29/95 16:09:48   llt
CD2 mlz defined twice. removed for cray.
CD2 
CD2    Rev 1.1   03/18/94 15:45:18   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:06   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 None
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 None: file names are identified by their tape numbers
CD3 
CD3 Name         Description
CD3 
CD3 iout         Output file
CD3 iptty        Alternate output file
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
CD4 n, ps, phi, iieos, pmin, pmax, tmin, tmax, ieos, t, s, iout, l,
CD4 day, iptty,  idof
CD4 
CD4 Global Subprograms
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
CD5 smind        real*8      Minimum allowable saturation value
CD5 smaxd        real*8      Maximum allowable saturation value
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 mlz          int         Flag denoting if out of bounds node is
CD5                              found
CD5 i            int         Do loop index over all nodes
CD5 pl           real*8      Current fluid pressure
CD5 ieo          int         Current phase state flag
CD5 tl           real*8      Current temperature
CD5 sl           real*8      Current saturation
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
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
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
CPS BEGIN outbnd
CPS 
CPS FOR each node
CPS 
CPS   IF the porosity is nonzero and idof greater than 1
CPS   
CPS     Set parameter values
CPS     
CPS     IF this node is one phase
CPS       Set flag to denote out of bounds value when temperature or...
CPS       ... pressure is outside the acceptable range
CPS     ELSE this node is two pahse
CPS       Set flag to denote out of bounds value when saturation or...
CPS       ... pressure is outside the acceptable range
CPS     ENDIF
CPS     
CPS     IF a value was outside the acceptable range
CPS       Write header information for this node
CPS       Write error information
CPS       IF this information is also written to second output file
CPS         Write header information for this node
CPS         Write error information
CPS       ENDIF
CPS       ERROREXIT
CPS     ENDIF
CPS   
CPS   ENDIF
CPS 
CPS ENDFOR each node
CPS 
CPS ERRORSEGMENT
CPS ENDSEGMENT
CPS 
CPS END outbnd
CPS 
C**********************************************************************

      use combi
      use comdi
      use comfi
      use comii
      use comdti
      use comai
      use davidi
      implicit none

      real*8 smind,smaxd,smind_t,smaxd_t
      real*8 pmin_ice,pmax_ice
      real*8 tmin_ice,tmax_ice
      parameter(smind = -2.,smaxd = 5.,smind_t = -0.5,smaxd_t = 1.5)
      parameter(pmin_ice=-1.0,pmax_ice=200.0)
      parameter(tmin_ice=-10.0,tmax_ice=400.0)
      integer i, ii
      real*8 pl
      integer ieo
      real*8 tl
      real*8 sl
      mlz=0
      if(ico2.ge.0) then
      do i=1,n
      if(ps(i).ne.0.0.and.idof.gt.1) then
         if(i.le.neq) then
            ii=i
         else if(i.ge.neq+1.and.i.le.neq+neq) then
            ii=i-neq
         else
            ii=i-neq-neq
         endif
c
         pl=phi(i)
c
         ieo=iieos(i)
         if(ieo.gt.2) ieo=2
c     
            tl=t(i)
            sl=s(i)
            if(pl.lt.pmin(ieo)) mlz=1
            if(pl.gt.pmax(ieo)) mlz=1
            if(sl.lt.smind_t) mlz=1
            if(sl.gt.smaxd_t) mlz=1
            if(tl.lt.tmin(ieo)) mlz=1
            if(tl.gt.tmax(ieo)) mlz=1
            if(ico2.gt.0) then
               if(pci(i).lt.0.0d00) mlz=1
            endif
         if(mlz.eq.1) then
            if (iout .ne. 0) then
               write(iout,*)'time step= ',l,' time step size= ',day
               write(iout,9001) i,cord(ii,1),cord(ii,2),cord(ii,3)
               write(iout,*)'p= ',pl,' t= ',tl,' s= ',sl
               if(ico2.gt.0) write(iout,*)'pci(i)= ',pci(i)
            end if
            if(iptty.gt.0) then
               write(iptty,*)'time step= ',l,' time step size= ',
     &              day
               write(iptty,9001) i,cord(ii,1),cord(ii,2),cord(ii,3)
               write(iptty,*)'p= ',pl,' t= ',tl,' s= ',sl
               if(ico2.gt.0) write(iptty,*)'pci(i)= ',pci(i)
            endif
            goto 9000
         endif
      end if
      enddo
c gaz 10-20-2001
      else if(ice.eq.0) then
c isothermal air water
      do i=1,n

      if(ps(i).ne.0.0.and.idof.gt.1) then
         if(i.le.neq) then
            ii=i
         else if(i.ge.neq+1.and.i.le.neq+neq) then
            ii=i-neq
         else
            ii=i-neq-neq
         endif	
         pl=phi(i)
c
         ieo=iieos(i)
         if(ieo.gt.2) ieo=2
c     
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            sl=s(i)
         else
            sl = 1.0d0
         end if
         if(pl.lt.pmin(ieo)) mlz=1
         if(pl.gt.pmax(ieo)) mlz=1
         if(sl.lt.smind) mlz=1
         if(sl.gt.smaxd) mlz=1
c
         if(mlz.eq.1) then
            if (iout .ne. 0) then
               write(iout,*)'time step= ',l,' time step size= ',day
               write(iout,9001) i,cord(ii,1),cord(ii,2),cord(ii,3)
               write(iout,*)'p= ',pl,' s= ',sl
            end if
            if(iptty.gt.0) then
               write(iptty,*)'time step= ',l,' time step size= ',
     &              day
               write(iptty,9001) i,cord(ii,1),cord(ii,2),cord(ii,3)
               write(iptty,*)'p= ',pl,' s= ',sl
            endif
            goto 9000
         endif
      end if
      enddo
c gaz 10-20-2001
      else if(ice.ne.0) then
      do i=1,n

      if(ps(i).ne.0.0.and.idof.gt.1) then
         if(i.le.neq) then
          ii=i
         else if(i.ge.neq+1.and.i.le.neq+neq) then
          ii=i-neq
         else
          ii=i-neq-neq
         endif	
c
         pl=phi(i)
c
         ieo=iieos(i)
         if(ieo.gt.2) ieo=2
c     
            tl=t(i)
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               sl=s(i)
            else
               sl = 1.0d0
            end if
            if(pl.lt.pmin_ice ) mlz=1
            if(pl.gt.pmax_ice ) mlz=1
            if(sl.lt.smind_t) mlz=1
            if(sl.gt.smaxd_t) mlz=1
            if(tl.lt.tmin_ice ) mlz=1
            if(tl.gt.tmax_ice ) mlz=1
         if(mlz.eq.1) then
            if (iout .ne. 0) then
               write(iout,*)'time step= ',l,' time step size= ',day
               write(iout,9001) i,cord(ii,1),cord(ii,2),cord(ii,3)
               write(iout,*)'p= ',pl,' t= ',tl,' s= ',sl
            end if
            if(iptty.gt.0) then
               write(iptty,*)'time step= ',l,' time step size= ',
     &              day
               write(iptty,9001) i,cord(ii,1),cord(ii,2),cord(ii,3)
               write(iptty,*)'p= ',pl,' t= ',tl,' s= ',sl
            endif
            goto 9000
         endif
      end if
      enddo
      endif
 9000 continue
 9001 format('out of bounds:node ',i8
     $ ,' x= ',g12.4,' y= ',g12.4,' z= ',g12.4) 
      return
      end

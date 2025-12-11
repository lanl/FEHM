      subroutine crdpdp
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
CD1 To update solution for double porosity/double permeability solution.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/crdpdp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:18 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Wed Feb 07 10:36:32 1996   gaz
CD2 changes to NR step length, strd passed 
CD2 
CD2    Rev 1.3   Mon Jan 29 14:02:30 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   08/03/95 14:43:26   gaz
CD2 added strd=pnx(1+n) for newton step length (2 places
CD2 
CD2    Rev 1.1   03/18/94 15:47:10   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:38   pvcs
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
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier   Type        Description
CD4 
CD4 
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 Identifier   Type        Description
CD4 
CD4 ico2, neq, phi, bp, s, t, nrhs, pnx, n, ieos, pci
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type       Description
CD4 
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
CD5 i            int         Do loop index over all nodes
CD5 nr1          int         Index for solution arrays
CD5 nr2          int         Index for solution arrays
CD5 nr3          int         Index for solution arrays
CD5 i1           int         Index for solution arrays
CD5 i2           int         Index for solution arrays
CD5 i3           int         Index for solution arrays
CD5 ipneq        int         Number of the the matrix node
CD5                              corresponding to this fracture node
CD5 ieosd        int         Current equation of state flag
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
CD9 2.4.9 Double-porosity/double-permeability formulation
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
CPS BEGIN crdpdp
CPS 
CPS IF this is an isothermal air-water solution
CPS 
CPS   FOR each node
CPS     Compute updated pressure and saturation
CPS   ENDFOR
CPS   
CPS ELSEIF this is pure water solution
CPS 
CPS   FOR each node
CPS     IF this node is pure liquid
CPS       Compute updated pressure and temperature
CPS     ELSEIF this node is two-phase
CPS       Compute updated pressure and saturation
CPS     ELSEIF this node is pure gas
CPS       Compute updated pressure and temperature
CPS     ELSEIF a heat transfer only solution is being carried out
CPS       Compute updated temperature
CPS     ENDIF
CPS   ENDFOR
CPS 
CPS ELSE this is a water-noncondensible gas solution
CPS 
CPS   FOR each node
CPS     IF this node is pure liquid
CPS       Compute updated pressure, temperature, and gas pressure
CPS     ELSEIF this node is two-phase
CPS       Compute updated pressure, saturation, and temperature
CPS     ELSEIF this node is pure gas
CPS       Compute updated pressure, temperature, and gas pressure
CPS     ELSEIF a heat transfer only solution is being carried out
CPS       Compute updated temperature
CPS     ENDIF
CPS   ENDFOR
CPS   
CPS   Load array for correction of fracture variables
CPS 
CPS ENDIF
CPS 
CPS END crdpdp
CPS 
C**********************************************************************
c
c update solution for double porosity/double permeability solution
c smoothed rl perms so always take largest step
c
      use davidi
      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      implicit none

      integer i,nr1,nr2,nr3,i1,i2,i3,ipneq,ieosd

      if(ico2.lt.0) then
c
c     air / water isothermal corrections
c
         do i=1,neq
            phi(i+neq)=phi(i+neq)-bp(i+nrhs(3))
            s(i+neq)=s(i+neq)-bp(i+nrhs(4))
         enddo
      else if(ico2.eq.0) then
c     
c     pure water
c
         nr1=nrhs(3)
         nr2=nrhs(4)
         do i=1,neq
            i1=i+nr1
            i2=i+nr2
            ieosd=ieos(i+neq)
            if(ieosd.eq.1) then
               phi(i+neq)=phi(i+neq)-bp(i1)*strd
               t(i+neq)=t(i+neq)-bp(i2)*strd
            elseif(ieosd.eq.2) then
               phi(i+neq)=phi(i+neq)-bp(i1)*strd
               s(i+neq)=s(i+neq)-bp(i2)*strd
            elseif(ieosd.eq.3) then
               phi(i+neq)=phi(i+neq)-bp(i1)*strd
               t(i+neq)=t(i+neq)-bp(i2)*strd
            elseif(ieosd.eq.-1) then
               t(i+neq)=t(i+neq)-bp(i2)*strd
            endif
         enddo
      else
c     
c     water and noncondensible
c     
         nr1=nrhs(4)
         nr2=nrhs(5)
         nr3=nrhs(6)
         do i=1,neq
            ipneq=i+neq
            i1=i+nr1
            i2=i+nr2
            i3=i+nr3
            ieosd=ieos(ipneq)
            if(ieosd.eq.1) then
               phi(ipneq)=phi(ipneq)-bp(i1)*strd
               t(ipneq)=t(ipneq)-bp(i2)*strd
               pci(ipneq)=pci(ipneq)-bp(i3)*strd
            elseif(ieosd.eq.2) then
               phi(ipneq)=phi(ipneq)-bp(i1)*strd
               s(ipneq)=s(ipneq)-bp(i2)*strd
               t(ipneq)=t(ipneq)-bp(i3)*strd
            elseif(ieosd.eq.3) then
               phi(ipneq)=phi(ipneq)-bp(i1)*strd
               t(ipneq)=t(ipneq)-bp(i2)*strd
               pci(ipneq)=pci(ipneq)-bp(i3)*strd
            elseif(ieosd.eq.-1) then
               t(ipneq)=t(ipneq)-bp(i2)*strd
            endif
            pci(ipneq) = max(pci(ipneq),0.0d00)
         enddo
      endif

      return
      end

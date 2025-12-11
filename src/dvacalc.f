      subroutine dvacalc
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
CD1 To evaluate air water/vapor diffusion coefficients.
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
CD2 $Log:   /pvcs.config/fehm90/src/dvacalc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:00:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Thu Feb 15 10:00:26 1996   zvd
CD2 Updated requirement.
CD2 
CD2    Rev 1.4   Mon Jan 29 15:54:20 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   01/28/95 13:54:40   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.2   06/15/94 10:37:00   robinson
CD2 Corrected term for air-water diffusion and derivatives
CD2 
CD2    Rev 1.1   03/18/94 16:11:28   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:40   pvcs
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
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 iadif, n, phi, t, tort, ps, rovf, s, dgve, dgvc, dva, ieos, ddvap,
CD4 ddvae, ddvac
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
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
CD5 dva0         real*8      Diffusion coefficient at reference
CD5                             conditions
CD5 theta        real*8      Exponent used in calculation
CD5 p0           real*8      Reference pressure
CD5 t0           real*8      Reference temperature
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
CD5 rat          real*8      Parameter used in calculation
CD5 dratp        real*8      Derivative of rat with pressure
CD5 dratt        real*8      Derivative of rat with temperature
CD5 dva0d        real*8      Parameter used in calculation
CD5 dva0p        real*8      Parameter used in calculation
CD5 dva0c        real*8      Parameter used in calculation
CD5 dva0e        real*8      Parameter used in calculation
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
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9 
CD9 2.4.2 Properties of air and air/water vapor mixtures
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
CPS BEGIN dvacalc
CPS 
CPS IF air water vapor diffusion coefficients need to be calculated
CPS 
CPS   FOR each node
CPS   
CPS     Compute parameters needed later
CPS     IF this is not an isothermal air-water simulation
CPS       Compute diffusion coefficient and derivatives
CPS     ELSE
CPS       Compute diffusion coefficient and derivatives
CPS     ENDIF
CPS   
CPS   ENDFOR each node
CPS 
CPS ENDIF
CPS 
CPS END dvacalc
CPS 
C**********************************************************************
C ----- PHS  REVISIONS    3/10/2004    FIX   dva(i) 
C -------   Take out rovf from this part    
C --------   adding different theta = 1.81  and Do = 2.23e-5
c--------------------------------------------------------------------

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use davidi
      use comhi
      use comgi
      use comdti
      use comai
      use com_exphase
      implicit none

      real*8 dva0,theta,p0,t0,rat,dratp,dratt,dva0d,dva0p,dva0c,dva0e
      real*8 t_min, t_max, atort, dratc, temp, temp2, diffcoeff
      real*8 tort2, dvas_denom_min, dvas_denom, s_dva_term, dva_tiny
c      real*8,  allocatable :: dva_save(:)
      real*8 dpsatt,dpsats,pv,dum_dva, psatl
      real*8 dva_t
      integer i
c gaz 123021  for explicit update   
      integer loop_start, loop_end
      parameter(dva0=2.23e-5)
      parameter(theta=1.810)
      parameter(p0=0.1)
      parameter(t0=273.15)
      parameter(t_min=10.0,t_max=350.0)
      parameter(dvas_denom_min = 1.d-18, dva_tiny = 1.d-18)
c      save dva_save
      if(.not.allocated(dva_save)) allocate (dva_save(n))

      if(iadif.ne.0.and.tort.ge.0.0.and.tort.le.1.0) then
c gaz 123020 manage explicit update 
       if(i_ex_update.ne.0.and.ieq_ex.gt.0) then
         loop_start = ieq_ex
         loop_end = ieq_ex
       else
         loop_start = 1
         loop_end = n         
       endif    
       do i=loop_start,loop_end
c         do i=1,n

c              PHS attempting to put in dratc     3/18/04
c                  using phi = pci + pvap (water vapor pressure)

c
c  New stuff   PHS took out density from rat
c              now no derivatives wrt density!
c              dgvc dgvp dgve go away!
c gaz 081116
            if(gdkm_flag.eq.0) then
             dva_t = dva0
            else
             dva_t = dva0*gdkm_volume_fraction(i)   
            endif    
            if(t(i).ge.t_min.and.t(i).le.t_max) then
               rat=(p0/phi(i))*((t(i)+t0)/t0)**theta
               dratp=-rat/phi(i)
               dratt=(p0/phi(i))*theta*(((t(i)+t0)/t0)**(theta-1.0))/t0
c     dratc=-rat/phi(i)
               dratc = 0
            else if(t(i).lt.t_min) then
               rat=(p0/phi(i))*((t_min+t0)/t0)**theta
               dratp=-rat/phi(i)
               dratt=0.0
               dratc=-rat/phi(i)
            else if(t(i).gt.t_max) then
               rat=(p0/phi(i))*((t_max+t0)/t0)**theta
               dratp=-rat/phi(i)
               dratt=0.0
               dratc=-rat/phi(i)
            endif
c
c            dva0d=tort*ps(i)*(1.0-s(i))*dva0
            dva0d=tort*ps(i)*(1.0-s(i))*dva_t
c
c parts of derivatives  for ieos=2  TotPres, Temp, GasPres
c                                    P        E      C
c                        are primary variables
c

            dva0p=0.0
            dva0c=0.0
            dva0e=0.0
            dva(i)=dva0d*rat + dva_tiny
            if(ieos(i).ne.2) then
               ddvap(i)=dva0d*dratp+dva0p*rat
               ddvae(i)=dva0d*dratt+dva0e*rat
               ddvac(i)=dva0d*dratc+dva0c*rat
            else
               ddvac(i)=-tort*ps(i)*dva0*rat
               ddvae(i)=0.0
               ddvap(i)=dva0d*dratp
            endif
         enddo



c---------------------------------------------
c-- PHS  cleaned up the negative tort section.
c-       dva = atort * dva0 * rat
c-      NEEDS CHECKED
c---------------------------------------------

      else if(iadif.ne.0.and.tort.lt.0.0) then
         do i=1,n
c
c  New stuff
c
            if(t(i).ge.t_min.and.t(i).le.t_max) then
               rat=(p0/phi(i))*((t(i)+t0)/t0)**theta
               dratp=-rat/phi(i)
               dratt=(p0/phi(i))*theta*((t(i)+t0)/t0)**(theta-1.0)/t0
            else if(t(i).lt.t_min) then
               rat=(p0/phi(i))*((t_min+t0)/t0)**theta
               dratp=-rat/phi(i)
               dratt=0.0
            else if(t(i).gt.t_max) then
               rat=(p0/phi(i))*((t_max+t0)/t0)**theta
               dratp=-rat/phi(i)
               dratt=0.0
            endif
c
            atort=abs(tort)
            dva0d = atort*dva0

            dva(i)=dva0d*rat + dva_tiny
            ddvap(i)=dva0d*dratp
            ddvae(i)=dva0d*dratt
            ddvac(i)=0.0

         enddo
c-------------------------------------------
c   NEW SECTION TO TIE THE ADIF DIFFUSION TO TRAC 
c   SPECIES #1   vapor diffusion model parameters
c   these are the same as in concadiff but with extra 
c     (1-s)*porosity
c-------------------------------------------------------

      else if(iadif.ne.0.and.tort.GT.1.0.and.tort.NE.333) then
         tort2 = tort
         do i=1,n

            diffcoeff = diffmfv(1,1)
            temp=(1-s(i))*ps(i)

            select case (mflagv(1,1))
            case (0)            ! Constant diffusion
               dva(i) =diffcoeff * temp + dva_tiny
               ddvae(i)= -ps(i)*diffcoeff
            case (1)            ! Millington Quirk
               temp2 = temp**2.3333333
               dva(i) =diffcoeff*temp*temp2/(ps(i)**2) + dva_tiny
               ddvae(i)= -3.33*diffcoeff*temp2/ps(i)
            case (2)            ! alternate Millington-Quirk
               dva(i) =diffcoeff*temp*temp/(ps(i)**0.6666) + dva_tiny
               ddvae(i)= -2*diffcoeff*temp*(ps(i)**0.3333)
            case (3)            ! using old adif form 
               if(tort.GT.666.) tort2 = tort - 666.
               if(t(i).ge.t_min.and.t(i).le.t_max) then
                  rat=(p0/phi(i))*((t(i)+t0)/t0)**theta
                  dratp=-rat/phi(i)
                  dratt=(p0/phi(i))*theta*(((t(i)+t0)/t0)**
     &                 (theta-1.0))/t0
c     dratc=-rat/phi(i)
                  dratc = 0
               else if(t(i).lt.t_min) then
                  rat=(p0/phi(i))*((t_min+t0)/t0)**theta
                  dratp=-rat/phi(i)
                  dratt=0.0
                  dratc=-rat/phi(i)
               else if(t(i).gt.t_max) then
                  rat=(p0/phi(i))*((t_max+t0)/t0)**theta
                  dratp=-rat/phi(i)
                  dratt=0.0
                  dratc=-rat/phi(i)
               endif
c
               dva0d=tort2*ps(i)*(1.0-s(i))*dva0
c
c parts of derivatives  for ieos=2  TotPres, Temp, GasPres
c                                    P        E      C
c                        are primary variables
c

               dva0p=0.0
               dva0c=0.0
               dva0e=0.0
               dva(i)=dva0d*rat + dva_tiny
               if(ieos(i).ne.2) then
                  ddvap(i)=dva0d*dratp+dva0p*rat
                  ddvae(i)=dva0d*dratt+dva0e*rat
                  ddvac(i)=dva0d*dratc+dva0c*rat
               else
                  ddvac(i)=-tort2*ps(i)*dva0*rat
                  ddvae(i)=0.0
                  ddvap(i)=dva0d*dratp
               endif

            case default
               goto 1000

            end select

            ddvap(i)= 0.0
            ddvac(i)= 0.0
          
         end do
c-------------------------------------------
c   NEW SECTION  Water Vapor following MillingtonQuirk
c    to calculated tortuosity for the dva f(PT) preuss equation.
c   tort2 is the toruosity from MQ plugged into the original 
c   section of code from the top of this routine.
c-------------------------------------------------------

      else if(tort.EQ.333) then
       if(i_ex_update.ne.0.and.ieq_ex.gt.0) then
         loop_start = ieq_ex
         loop_end = ieq_ex
       else
         loop_start = 1
         loop_end = n          
       endif 
       do i = loop_start, loop_end
           temp=(1-s(i))*ps(i)
           if(ps(i).GT.0) then
            tort2 = temp**2.3333/(ps(i)**2)
           else
            tort2 = 0.0
           end if
c
c  New stuff   PHS took out density from rat
c              now no derivatives wrt density! 
c              dgvc dgvp dgve go away! 
c
           if(t(i).ge.t_min.and.t(i).le.t_max) then
              rat=(p0/phi(i))*((t(i)+t0)/t0)**theta
              dratp=-rat/phi(i)
              dratt=(p0/phi(i))*theta*(((t(i)+t0)/t0)**(theta-1.0))/t0
c     dratc=-rat/phi(i)
              dratc = 0
           else if(t(i).lt.t_min) then
              rat=(p0/phi(i))*((t_min+t0)/t0)**theta
              dratp=-rat/phi(i)
              dratt=0.0
c     dratc=-rat/phi(i)
              dratc = 0.0
           else if(t(i).gt.t_max) then
              rat=(p0/phi(i))*((t_max+t0)/t0)**theta
              dratp=-rat/phi(i)
              dratt=0.0
c              dratc=-rat/phi(i)
              dratc = 0.0
           endif
c
           dva0d=tort2*ps(i)*(1.0-s(i))*dva0
c
c parts of derivatives  for ieos=2  TotPres, Temp, GasPres
c                                    P        E      C
c                        are primary variables
c

           dva0p=0.0
           dva0c=0.0
           dva0e=0.0
         if(iad.eq.0) then
c explicit update
           dva(i)=dva0d*rat + dva_tiny
           dva_save(i) = dva(i)
c           dvas(i) = dva(i)/(temp)
c          if(iatty.NE.0) write(iatty,*) 'dvacalc dva=',dva(i),tort2,dva0

           if(ieos(i).ne.2) then
              ddvap(i)=dva0d*dratp+dva0p*rat
              ddvae(i)=dva0d*dratt+dva0e*rat
              ddvac(i)=dva0d*dratc+dva0c*rat
           else
c              if(isalt.ne.0) then
c               call saltctr(1,i,dpsatt,dpsats)
c              else
c               pv = psatl(t(i),0.0d0,dum_dva,dpsatt,dpsats,0)
c              endif
c using dpsatt = 1./dtsatp (dps/dts = 1/(dts/dps)
c gaz debug 012015 the ddvac and ddvae seem reversed but it works
c              ddvac(i)=-tort2*ps(i)*dva0*rat
cc              ddvae(i)=dva0d*(dratt)
c              ddvap(i)=dva0d*(dratp + dratt/dpsatt)
c              ddvap(i)=dva0d*(dratp )
           endif
              ddvae(i)=0.0
              ddvac(i)=0.0
              ddvap(i)=0.0
          else
          dva(i) = dva_save(i) 
          endif
        enddo      ! (i = 1,n)

      endif     !   (tort.EQ.333)

c - - - - - - - - 7/17/13  PHS  Load dva/(air content) into dvas - - - - - 
c gaz debug 090314
      if(irdof.ne.13) then
c gaz 123020 manage explicit update 
       if(i_ex_update.ne.0.and.ieq_ex.gt.0) then
         loop_start = ieq_ex
         loop_end = ieq_ex
       else
         loop_start = 1
         loop_end = n          
       endif    
       do i=loop_start,loop_end          
         if((s(i).LT.1).AND.(ps(i).GT.0)) then
          dvas(i) = dva(i)/(ps(i)*(1-s(i))+dvas_denom_min)
         else
          dvas(i) = 0.0
         end if 
       end do
      endif
      return

 1000 write(ierr,*) 'illegal vapor diffusion model'
      if (iptty .ne. 0) write(iptty,*) 'illegal vapor diffusion model'
      stop 

      end

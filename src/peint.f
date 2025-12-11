       subroutine peint
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine initializes the pressures and temperatures when 
CD1  a temperature gradient exists.  The pressures are initialized to
CD1  hydrostatic.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/peint.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:38   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Fri Feb 16 10:52:46 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.2   Thu Jan 11 08:56:00 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:43:12   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:26:16   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable.  See Special Comments
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  This is a general utility routine used to initialize pressure and 
CD4  temperature fields.
CD4  
C***********************************************************************
C*****
C*****AF 11/15/10 updated for lookup table
C*****

      use combi
      use comci
      use comdi
      use comfi
      use comii
      use comdti
      use comai
      use comtable
      implicit none
C*****
C*****AF 11/15/10 updated for lookup table
C*****
C      include 'comtable.h'        ! phs 4/24/99
C*****

      integer i, ij, iieosd
      integer ndep, ndep1
      real*8  cord0, cordd
      real*8  delx, depth, dla0, dlb0
      real*8  dlpa1, dlpa2, dlpa3, dlpb1, dlpb2, dlpb3
      real*8  dlpt2a, dlpt2b, dlp2ta, dlp2tb, dlpta, dlptb
      real*8  dlta1, dlta2, dlta3, dltb1, dltb2, dltb3
      real*8  pl, pl2, pl2tl, pl3, pltl, pltl2, press
      real*8  rol, rold, roln, rolpd, rolpn
      real*8  rolptd, rolptn, roltd, roltn
      real*8  tl, tl2, tl3
C*****
C*****AF 11/15/10 updated for lookup table
C*****
      real*8  zwp, zwt ! phs 4/25/99
      integer indexp, indext, point(4),izerrFLAG ! phs 4/25/99
C*****

      ndep   =  80
      ndep1  =  ndep+1
      depcng = abs( depcng )
      
c     find top elevation
      cord0=-10000000000000.0
      if ( icnl .eq. 0 )  then
         do i=1,neq
            cord0=max(cord0,cord(i,3))
         enddo
      else
         do i=1,neq
            cord0=max(cord0,cord(i,2))
         enddo
      endif
      cord0 =  abs( cord0 )

      do ij=1,neq
c**** find depth ****
         if ( icnl .eq. 0 )  then
c**** 3-dimension ****
            cordd  =  abs( cord(ij,3) )
         else
c**** 2-dimension ****
            cordd  =  abs( cord(ij,2) )
         endif
c**** heat conduction only ****
         press  =  pein
         iieosd =  iieos(ij)
c         tin    =  max( tin ,tmin(iieosd) )
c         tin1   =  max( tin1,tmin(iieosd) )

         if ( cordd .le. depcng )  tl =  tin +gradnt*cordd
         if ( cordd .gt. depcng )  tl =  tin1+grad2 *cordd+quad*cordd**2
c     loop to calc pressure
c**** liquid phase only ****
         delx   =  ( cordd-cord0 )/ndep
C     gaz 3-16-2002   depth  =  cord0+delx
         depth  =  cord0-delx
         press  =  pein

         do i=1,ndep1
            depth  =  depth+delx
c**** set temperture ****
            if ( depth .gt. depcng )  tl =  tin +gradnt*depth
            if ( depth .le. depcng )  tl =  tin1+grad2 *depth +
     &           quad*depth**2
            if ( i .ne. 1 )  then
c**** calculate density ****
               if ( abs( grav ) .gt. zero_t )  then
                  pl     =  press
c**** chose correct pressure region ****
C*****
C*****AF 11/15/10
C*****
c     
c-----------------------------------------------------
c-----------------use coefficients if tableFLAG.NE.1    phs 4/25/99
c     
                  if (tableFLAG.NE.1) then ! phs 4/25/99
C*****
C*****AF 11/15/10
C*****
c**** liquid density * numerator coefficients ****
                     dla0   =  crl( 1,iieosd)
                     dlpa1  =  crl( 2,iieosd)
                     dlpa2  =  crl( 3,iieosd)
                     dlpa3  =  crl( 4,iieosd)
                     dlta1  =  crl( 5,iieosd)
                     dlta2  =  crl( 6,iieosd)
                     dlta3  =  crl( 7,iieosd)
                     dlpta  =  crl( 8,iieosd)
                     dlp2ta =  crl( 9,iieosd)
                     dlpt2a =  crl(10,iieosd)
c**** liquid density * denomenator coefficients ****
                     dlb0   =  crl(11,iieosd)
                     dlpb1  =  crl(12,iieosd)
                     dlpb2  =  crl(13,iieosd)
                     dlpb3  =  crl(14,iieosd)
                     dltb1  =  crl(15,iieosd)
                     dltb2  =  crl(16,iieosd)
                     dltb3  =  crl(17,iieosd)
                     dlptb  =  crl(18,iieosd)
                     dlp2tb =  crl(19,iieosd)
                     dlpt2b =  crl(20,iieosd)

                     pl2    =  pl**2
                     pl3    =  pl**3
                     tl2    =  tl**2
                     tl3    =  tl**3
                     pltl   =  pl *tl
                     pl2tl  =  pl2*tl
                     pltl2  =  pl *tl2

c**** liquid density ****
                     rolpn  =  dla0+dlpa1*pl  +dlpa2 *pl2  +dlpa3 *pl3
                     roltn  =       dlta1*tl  +dlta2 *tl2  +dlta3 *tl3
                     rolptn =       dlpta*pltl+dlp2ta*pl2tl+dlpt2a*pltl2
                     roln   =  rolpn+roltn+rolptn
                     rolpd  =  dlb0+dlpb1*pl  +dlpb2 *pl2  +dlpb3 *pl3
                     roltd  =       dltb1*tl  +dltb2 *tl2  +dltb3 *tl3
                     rolptd =       dlptb*pltl+dlp2tb*pl2tl+dlpt2b*pltl2
                     rold   =  rolpd+roltd+rolptd
                     rol    =  roln/rold
C*****
C*****AF 11/15/10
C*****
c----------------------------------------
                  else          ! if tableFLAG.EQ.1          phs 4/25/99

                     izerrFLAG = 0.

                     indexp = pmin(1) + incp*dint((pl-pmin(1))/incp)
                     indext = tmin(1) + inct*dint((tl-tmin(1))/inct)
c---  find 4 points            LOOKUP

                     point(1) = 1 + (((indexp-pmin(1))/incp)*numt)
     x                    +  ((indext-tmin(1))/inct)
                     point(2) = point(1) + numt
                     point(3) = point(2) + 1
                     point(4) = point(1) + 1
c---  
                     if(PP(point(1),3).LT.0) izerrFLAG = 1.
                     if(PP(point(2),3).LT.0) izerrFLAG = 1.
                     if(PP(point(3),3).LT.0) izerrFLAG = 1.
                     if(PP(point(4),3).LT.0) izerrFLAG = 1.
                     if(izerrFLAG.EQ.1.) then
                        write(6,*)  'Out of bounds in peinit.f'
                        stop
                     endif


c---  compute weights for the P and T direction
                     zwp = (pl-indexp)/incp
                     zwt = (tl-indext)/inct
c---  find values as function of P+T                 LOOKUP

                     rol=(1-zwp)*(1-zwt)*PP(point(1),3) + (1-zwt)*zwp*
     &                    PP(point(2),3) + zwt*zwp*PP(point(3),3) +
     &                    (1-zwp)*zwt*PP(point(4),3)

                  end if        
c     tableFLAG or polynomial fit for water density.    phs 4/25/99
c--------------------------------------------------------------------------
C*****
C*****AF 11/15/10
C*****

                  press  =  press-rol*grav*delx
               endif
            endif
         enddo
c     end loop to calc press
         if ( abs( pho(ij) ) .lt. zero_t )  pho(ij) =  press
         if ( abs( to (ij) ) .lt. zero_t )  to (ij) =  tl
         t(ij) = to(ij)
         phi(ij) = pho(ij)
         if (idpdp.ne.0.or.idualp.ne.0) then
            pho(ij+neq)=pho(ij)
            to(ij+neq) =to(ij)
            t(ij+neq) = to(ij)
            phi(ij+neq) = pho(ij)           
         endif
         if(idualp.ne.0) then
            pho(ij+neq+neq)=pho(ij)
            to(ij+neq+neq) =to(ij)
            t(ij+neq+neq) = to(ij)
            phi(ij+neq+neq) = pho(ij)           
         endif
      enddo
      if ( iout .ne. 0 )  write(iout  ,6000)
      if ( iatty .gt. 0 )  write(iatty ,6000)
 6000 format(/,1x,'pressures and temperatures set by gradients')

      return
      end

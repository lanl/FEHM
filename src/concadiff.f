      real*8 function concadiff(flag,mflg,diffcoeff,poros,satr,vpr,tpr)
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
CD1 Compute diffusion coefficient (m^2/sec) based on volumetric water
CD1 content
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 $Log:   /pvcs.config/fehm90/src/concadiff.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:42   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:00:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:07:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:58:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:04 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Jun 27 13:51:18 1996   zvd
CD2 Corrected typo
CD2 
CD2    Rev 1.3   Thu Jun 27 13:34:32 1996   zvd
CD2 Added error write to ierr
CD2 
CD2    Rev 1.2   Wed Jun 26 12:31:04 1996   hend
CD2 Updated Prolog
CD2 
CD2    Rev 1.1   Fri Jun 07 11:25:16 1996   hend
CD2 Checked for diffusion coeff. droping under 0.
CD2 
CD2    Rev 1.0   Wed May 29 14:32:28 1996   hend
CD2 Added variable diffusion with water content
CD2 
C**********************************************************************
CD3
CD3 SPECIAL COMMENTS AND REFERENCES
CD3
CD3 Method Derived from Work of James L Conca, Pacific Northwest Lab.
CD3 "Diffusion Barrier Transport Properties of Unsaturated Paintbrush
CD3 Tuff Rubble Backfill", High Level Radioactive Waste Management 
CD3 Volume 1, pp. 394-401
CD3
CD3 Requirements from SDN: 10086-RD-2.20-00
CD3   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD3   FEHM Application Version 2.20
CD3
C**********************************************************************
CD4
CD4  REQUIREMENTS TRACEABILITY
CD4
CD4  2.3.4 Solute-transport equations
CD4
C**********************************************************************

!      use comdi
      use comai
!      use comdti
      implicit none
      
      real*8 diffcoeff, temp, temp2, poros, satr, vpr, tpr, satr2
      real*8 theta,p0,t0
      integer flag,mflg

      parameter(theta=1.810)
      parameter(p0=0.1)
      parameter(t0=273.15)

!  PHIL CHANGED TO include both MQ and Conca-Wright formulation
!  5/20/2003
!  mflg (input mflag) tells which diffusion model to use
!  Input from the read deck now give free water free air diffusion 
!  coefficient
C  Conca - from the fit given on the Neptune 2003 Diffusion write-up
c  Uses MQ for the air-phase since Conca is water only 
c  For MQ2 in the Vapor use Thetaresid Vapor = 666

      if (satr .gt. 1.0) then
         satr2 = 1.0
      else if (satr .lt. 0.) then
         satr2 = 0.
      else
         satr2 = satr
      end if

! flag = 1 Liquid
      if (flag.eq.1) then
         temp=satr2*poros
         select case (mflg)
         case (0)            ! Constant diffusion
            concadiff = diffcoeff
         case (1)            ! Millington Quirk
            temp2 = temp**2.3333333
            concadiff = diffcoeff*temp2/(poros**2)
         case (2)            ! Conca
            if(temp.EQ.0) temp=0.00001
            temp2 = -4.1 + 2.7*dlog10(temp) + 0.32*((dlog10(temp))**2)
            concadiff = (10**(temp2) / temp) * 1e-4
C Note 1e-4 term in line above is to convert to m2/s from cm2/s
         case default
            goto 1000
         end select
! flag = 2 Vapor
      else if (flag .eq. 2) then
         temp=(1-satr2)*poros
         select case (mflg)
         case (0)      ! Constant diffusion
            concadiff = diffcoeff
         case (1)      ! Millington Quirk
            temp2 = temp**2.3333333
            concadiff = diffcoeff*temp2/(poros**2)
         case (2)      ! alternate Millington-Quirk
            concadiff = diffcoeff*temp/(poros**0.6666) 
         case (3)      ! vapor diffusion same as water vapor
            concadiff = tort*2.23e-5*((p0/vpr)*((tpr+t0)/t0)**theta)
         case default
            goto 1000
         end select
      endif

      if (concadiff.lt.0.) concadiff=1e-30
      return

 1000 write(ierr, *) 'ERROR -- Illegal Flag to concadiff'
      write(ierr, *) 'Code Aborted in concadiff'
      if(iptty .ne. 0) then
         write(iptty, *) 'ERROR -- Illegal Flag to concadiff'
         write(iptty, *) 'Code Aborted in concadiff'
      end if
      stop

      end
